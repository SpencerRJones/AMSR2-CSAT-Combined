import numpy as np
import os
import xarray as xr
import glob
import matplotlib.pyplot as plt

#---------------------------------------------------------
#
#    Function: regrid_sfc_data
#
#    Regrids data from native ERA5 surface resolution
#        (0.28125 degrees) to profile data resolution
#        (0.5 degrees) using a footprint average.
#
#----------------------------------------------------------

def regrid_sfc_data(data, old_lats, old_lons, new_deg):

    new_lats = np.arange(-90.0, 90.0+new_deg, new_deg)
    new_lons = np.arange(0, 360, new_deg)

    new_data = np.zeros([new_lats.size,new_lons.size])
    new_data[:] = -9999.9

    factor = 1 / new_deg

    old_lats_rounded = np.round(old_lats * factor).astype(int)
    old_lons_rounded = np.round(old_lons * factor).astype(int)

    new_lats_xfactor = (new_lats * factor).astype(int)
    new_lons_xfactor = (new_lons * factor).astype(int)

    for i,ilat in enumerate(new_lats):

        target_lat      = new_lats_xfactor[i]
        newlat          = ilat
        lat_window      = np.where(old_lats_rounded == target_lat)
        beg_lat_indx    = lat_window[0][0]
        end_lat_indx    = lat_window[0][-1] + 1

        for j,ilon in enumerate(new_lons):

            target_lon   = new_lons_xfactor[j]
            newlon       = ilon
            lon_window   = np.where(old_lons_rounded == target_lon)
            beg_lon_indx = lon_window[0][0]
            end_lon_indx = lon_window[0][-1] + 1

            data_window = data[beg_lat_indx:end_lat_indx,beg_lon_indx:end_lon_indx]

            new_data[i,j] = np.mean(data_window)

    return new_data




def create_layer_means_log(level_var):

    nlevs = level_var.shape[-1]
    nlyrs = nlevs - 1

    layer_var = np.zeros([level_var.shape[0], nlyrs])

    for ilyr in np.arange(0,nlyrs):
        d1 = level_var[:,ilyr]
        d2 = level_var[:,ilyr+1]

        d1[np.where(d1 <= 0.)] = 1.0e-10 #avoid log(0)
        d2[np.where(d2 <= 0.)] = 1.0e-10

        where_same = np.where(d1==d2)[0] #Avoid divide by 0

        layer_var[:,ilyr] = (d2 - d1) / (np.log(d2) - np.log(d1))

        layer_var[where_same,ilyr] = d1[where_same]

    return layer_var



def create_layer_means(level_var):

    nlevs = level_var.shape[1]
    nlyrs = nlevs - 1

    layer_var = np.zeros([level_var.shape[0],nlyrs])

    for ilyr in np.arange(0,nlyrs):
        d1 = level_var[:,ilyr]
        d2 = level_var[:,ilyr+1]

        layer_var[:,ilyr] = np.mean([d1,d2], axis=0)

    return layer_var


def mass2wc(mass, T, p, q):

    p    = p*100.    #[hPa --> Pa]
    q    = q/1000.   #[g/kg --> kg/kg]

    R_d = 287.0

    vol = R_d * (1./p) * (1 + 0.61*q) * T

    wc = mass / vol

    return wc


#----------------------------------------------------------------

def interpolate_linear(old_hgt, new_hgt, data):

    new_data = np.zeros(new_hgt.shape, dtype='f')
    new_data[:] = -9999.9

    j = 0

    for i in np.arange(0,new_hgt.size):

        ihgt = new_hgt[i]

        if ihgt == 0.:
            new_data[i] = data[-1]
            continue

        j = np.where(ihgt > old_hgt)[0][0]

        lo_bound = old_hgt[j-1]
        up_bound = old_hgt[j]

        #print(lo_bound, ihgt, up_bound, i,j, ihgt, old_hgt[j])

        w1 = 1./np.abs(ihgt - up_bound)
        w2 = 1./np.abs(ihgt - lo_bound)
        new_data[i] = ((data[j] * w1) + (data[j-1] * w2)) / (w1 + w2)

    return new_data

#------------------------------------------------------------------------------


def interpolate_log(old_hgt, new_hgt, data):

#     old_hgt = np.flipud(old_hgt)
#     data    = np.flipud(data)

    data = np.log(data.copy())

    new_data = np.zeros(new_hgt.shape, dtype='f')
    new_data[:] = -9999.9

    j = 0

    for i in np.arange(0,new_hgt.size):

        ihgt = new_hgt[i]

        if ihgt == 0.:
            new_data[i] = data[-1]
            continue

        j = np.where(ihgt > old_hgt)[0][0]

        lo_bound = old_hgt[j-1]
        up_bound = old_hgt[j]

        #print(lo_bound, ihgt, up_bound, i,j, ihgt, old_hgt[j])

        w1 = 1./np.abs(ihgt - up_bound)
        w2 = 1./np.abs(ihgt - lo_bound)
        new_data[i] = ((data[j] * w1) + (data[j-1] * w2)) / (w1 + w2)

    return np.exp(new_data)


def q2RH(q, T, p):

    #------------------------------------------------------------
    #
    # e/e_s, where e_s(T) = 611.2exp((17.67T_C)/(T_C + 243.5))
    #
    # Petty, A First Course in Atm. Thermodynamics, p. 183, eq. 7.19
    #
    # Inputs: q[kg/kg], T[K], p[Pa]
    #
    # Output: RH[%]
    #
    #------------------------------------------------------------

    R_v = 461.5 #[J kg^-1 K^-1]
    R_d = 287.0

    eps = R_d/R_v

    T_C = T - 273.15

    e = (p * q) / (eps + ((1. - eps) * q))

    e_s = 611.2 * np.exp( (17.67 * T_C) / (T_C + 243.5) )

    RH = (e / e_s) * 100.

    return RH




#----------------------------------------------------
#
#   ERA5 NATIVE
#
#----------------------------------------------------


era5_dir = '/qdata2/archive/ERA5/'
era5wc_dir = '/pdata4/pbrown/ERA5/'
era5pt_dir = '/pdata4/pbrown/ERA5/'
out_dir  = '/edata2/spencer/thesis_data/era5_ocean_profiles/err_covar_amsr2_csat/'

beg_date = '2015-10-26'
end_date = '2015-10-31'

beg_date = np.datetime64(beg_date)
end_date = np.datetime64(end_date)

times = [0,6,12,18]
ntimes = len(times)


for idate in np.arange(beg_date,end_date+np.timedelta64(1, 'D')):

    era5_date = str(idate)
    era5_date = f'{era5_date[0:4]}{era5_date[5:7]}{era5_date[8:]}'

    #Get the files
    surf_file = glob.glob(f'{era5_dir}{era5_date[:-2]}/*{era5_date}*_surf.nc')[0]
    plev_file = glob.glob(f'{era5_dir}{era5_date[:-2]}/*{era5_date}*_plev.nc')[0]
    wc_file   = glob.glob(f'{era5wc_dir}{era5_date[:-2]}/*{era5_date}*_wc.nc')[0]
    pt_file   = glob.glob(f'{era5pt_dir}{era5_date[:-2]}/*{era5_date}*prcptype.nc')[0]

    print(surf_file, plev_file, wc_file, pt_file)

    print('Reading surface data...')
    with xr.open_dataset(surf_file) as f:
        sfclat = f.latitude.values
        sfclon = f.longitude.values
        sst_in = f.sst[times,:,:].values
        wsp_in = np.sqrt(f.u10**2 + f.v10**2).values
        t2m_in = f.t2m[times,:,:].values
        si_in  = f.siconc[times,:,:].values
        sp_in  = f.sp[times,:,:].values / 100.


    print('Reading profile data...')
    with xr.open_dataset(plev_file) as f:
        plvlat = f.latitude.values
        plvlon = f.longitude.values
        plev   = f.level.values

        z_in = f.z[times,:,:,:].values / 9.80665
        t_in = f.t[times,:,:,:].values
        q_in = f.q[times,:,:,:].values * 1000.
        q_cw_in = f.clwc[times,:,:,:].values * 1000.
        q_ci_in = f.ciwc[times,:,:,:].values * 1000.


    print('Reading wc data...')
    with xr.open_dataset(wc_file) as f:
        q_pw_in = f.crwc.values * 1000.
        q_pi_in = f.cswc.values * 1000.

#     print('Reading precip type data...')
#     with xr.open_dataset(pt_file) as f:
#         pt_in = f.ptype.values
#         pt_in = pt_in.round().astype('i')


    nslats = sfclat.size
    nslons = sfclon.size
    nplats = plvlat.size
    nplons = plvlon.size

    #Screen out sea ice and land:
    sst_in[np.where(si_in > 0.)] = np.nan


    #Regrid surface data:
    print('Regridding surface data...')
    sst_05 = np.zeros([ntimes,nplats,nplons])
    wsp_05 = np.zeros([ntimes,nplats,nplons])
    t2m_05 = np.zeros([ntimes,nplats,nplons])
    sp_05  = np.zeros([ntimes,nplats,nplons])

    for itime in np.arange(0,ntimes):
        sst_05[itime,:,:] = regrid_sfc_data(sst_in[itime,:,:], sfclat, sfclon, 0.5)
        wsp_05[itime,:,:] = regrid_sfc_data(wsp_in[itime,:,:], sfclat, sfclon, 0.5)
        t2m_05[itime,:,:] = regrid_sfc_data(t2m_in[itime,:,:], sfclat, sfclon, 0.5)
        sp_05[itime,:,:]  = regrid_sfc_data(sp_in[itime,:,:], sfclat, sfclon, 0.5)


    sst  = []
    wsp  = []
    t2m  = []
    sp   = []
    temp = []
    pres = []
    hum  = []
    hgt  = []
    clwc = []
    ciwc = []
    plwc = []
    piwc = []

    #Extract only good ocean data
    first = True
    print('Extracting good data...')
    for itime in np.arange(0,ntimes):

        good_lats = np.where(~np.isnan(sst_05[0,:,:]))[0]
        good_lons = np.where(~np.isnan(sst_05[0,:,:]))[1]

        if first:
            lat  = plvlat[good_lats]
            lon  = plvlon[good_lons]
            hgt  = z_in[itime,:,good_lats,good_lons]
            temp = t_in[itime,:,good_lats,good_lons]
            hum  = q_in[itime,:,good_lats,good_lons]
            q_cw = q_cw_in[itime,:,good_lats,good_lons]
            q_ci = q_ci_in[itime,:,good_lats,good_lons]
            q_pw = q_pw_in[itime,:,good_lats,good_lons]
            q_pi = q_pi_in[itime,:,good_lats,good_lons]
            sst  = sst_05[itime,good_lats,good_lons]
            wsp  = wsp_05[itime,good_lats,good_lons]
            t2m  = t2m_05[itime,good_lats,good_lons]
            sp   = sp_05[itime,good_lats,good_lons]

            first = False

        else:

            lat = np.append(lat, plvlat[good_lats])
            lon = np.append(lon, plvlon[good_lons])
            hgt = np.vstack((hgt, z_in[itime,:,good_lats,good_lons]))
            temp = np.vstack((temp, t_in[itime,:,good_lats,good_lons]))
            hum  = np.vstack((hum, q_in[itime,:,good_lats,good_lons]))
            q_cw = np.vstack((q_cw, q_cw_in[itime,:,good_lats,good_lons]))
            q_ci = np.vstack((q_ci, q_ci_in[itime,:,good_lats,good_lons]))
            q_pw = np.vstack((q_pw, q_pw_in[itime,:,good_lats,good_lons]))
            q_pi = np.vstack((q_pi, q_pi_in[itime,:,good_lats,good_lons]))
            sst  = np.append(sst, sst_05[itime,good_lats,good_lons])
            wsp  = np.append(wsp, wsp_05[itime,good_lats,good_lons])
            t2m  = np.append(t2m, t2m_05[itime,good_lats,good_lons])
            sp   = np.append(sp, sp_05[itime,good_lats,good_lons])

    hgt  = np.flip(hgt,  axis=1)
    temp = np.flip(temp, axis=1)
    hum  = np.flip(hum, axis=1)
    q_cw = np.flip(q_cw, axis=1)
    q_ci = np.flip(q_ci, axis=1)
    q_pw = np.flip(q_pw, axis=1)
    q_pi = np.flip(q_pi, axis=1)

    nprofs = sst.size
    nlevs  = plev.size

    pres = np.flipud(plev)


    #---Interpolate data:
    hgt_radar = np.arange(15000,-500,-500)
    nhgts     = hgt_radar.size
    temp_int  = np.zeros([nprofs, nhgts])
    pres_int  = np.zeros([nprofs, nhgts])
    hum_int   = np.zeros([nprofs, nhgts])
    q_cw_int  = np.zeros([nprofs, nhgts])
    q_ci_int  = np.zeros([nprofs, nhgts])
    q_pw_int  = np.zeros([nprofs, nhgts])
    q_pi_int  = np.zeros([nprofs, nhgts])


    print('Interpolating data...')
    percent = 0.1
    for iprof in np.arange(0, nprofs):
    #for iprof in np.arange(0, 10000):
        if iprof / nprofs >= percent:
            print(f'{percent*100.}%')
            percent += 0.1
        temp_int[iprof,:] = interpolate_linear(hgt[iprof], hgt_radar, temp[iprof])
        pres_int[iprof,:] = interpolate_log(hgt[iprof], hgt_radar, pres)
        hum_int[iprof,:]  = interpolate_log(hgt[iprof], hgt_radar, hum[iprof])
        q_cw_int[iprof,:] = interpolate_linear(hgt[iprof], hgt_radar, q_cw[iprof])
        q_ci_int[iprof,:] = interpolate_linear(hgt[iprof], hgt_radar, q_ci[iprof])
        q_pw_int[iprof,:] = interpolate_linear(hgt[iprof], hgt_radar, q_pw[iprof])
        q_pi_int[iprof,:] = interpolate_linear(hgt[iprof], hgt_radar, q_pi[iprof])

        #Assign surface values
        pres_int[iprof,-1] = sp[iprof]
        temp_int[iprof,-1] = t2m[iprof]


    #---Convert from cloud water mass mixing ratio to clwc [kg/kg --> g/m^3]
    clwc = mass2wc(q_cw_int, temp_int, pres_int, hum_int)
    ciwc = mass2wc(q_ci_int, temp_int, pres_int, hum_int)
    plwc = mass2wc(q_pw_int, temp_int, pres_int, hum_int)
    piwc = mass2wc(q_pi_int, temp_int, pres_int, hum_int)

    #---Create layer variables:
    print('Creating layer variables:')
    qlyr = create_layer_means_log(hum_int)
    clwclyr = create_layer_means(clwc)
    ciwclyr = create_layer_means(ciwc)
    plwclyr = create_layer_means(plwc)
    piwclyr = create_layer_means(piwc)
    tlyr = create_layer_means(temp_int)
    plyr = create_layer_means_log(pres_int)


    #---Create snow and rain dsd parameters:
    print('Creating DSD parameters:')
    N_0_rain = np.zeros([nprofs], dtype='f')
    N_0_snow = np.zeros([nprofs], dtype='f')
    mu_rain  = np.zeros([nprofs], dtype='f')
    rho_snow = np.zeros([nprofs], dtype='f')
    N_0_rain[:] = 1.1e+05
    N_0_snow[:] = 5100
    #mu is allowed to vary from 0 to 2.5, so assign random value between these
    mu_rain[:]  = np.random.choice(np.arange(0.,2.51,0.01), size=nprofs)
    #assign random snow particle density from 0.05 to 0.4 g/cm^3
    rho_snow[:] = np.random.choice(np.arange(50.,401., 1.), size=nprofs)


    #---Write output file
    print('Writing output file...')

    outfname = f'era5_ocean_profiles_{era5_date}_native.bin'
    outfile = f'{out_dir}{outfname}'

    temp_final = temp_int
    pres_final = pres_int
    qlyr_final = qlyr
    clwclyr_final = clwclyr
    plwclyr_final = plwclyr
    tiwclyr_final = ciwclyr + piwclyr #cloud ice and snow both represented by same PSD

    hgt_final = np.zeros([nprofs,nhgts])
    for iprof in np.arange(0,nprofs):
        hgt_final[iprof,:] = hgt_radar[:] / 1000.

    with open(outfile, 'wb') as f:
        np.array(nprofs).astype('i').tofile(f, sep='', format='unformatted')
        np.array(nhgts).astype('i').tofile(f, sep='', format='unformatted')
        lat.astype('f').tofile(f, sep='', format='unformatted')
        lon.astype('f').tofile(f, sep='', format='unformatted')
        hgt_final.T.astype('f').tofile(f, sep='', format='unformatted')
        temp_final.T.astype('f').tofile(f, sep='', format='unformatted')
        pres_final.T.astype('f').tofile(f, sep='', format='unformatted')
        qlyr_final.T.astype('f').tofile(f, sep='', format='unformatted')
        clwclyr_final.T.astype('f').tofile(f, sep='', format='unformatted')
        plwclyr_final.T.astype('f').tofile(f, sep='', format='unformatted')
        tiwclyr_final.T.astype('f').tofile(f, sep='', format='unformatted')
        sst.astype('f').tofile(f, sep='', format='unformatted')
        t2m.astype('f').tofile(f, sep='', format='unformatted')
        wsp.astype('f').tofile(f, sep='', format='unformatted')
        N_0_rain.astype('f').tofile(f, sep='', format='unformatted')
        mu_rain.astype('f').tofile(f, sep='', format='unformatted')
        N_0_snow.astype('f').tofile(f, sep='', format='unformatted')
        rho_snow.astype('f').tofile(f, sep='', format='unformatted')

    lat = []
    lon = []
    hgt = []
    temp = []
    hum = []
    pres = []
    qlyr = []
    clwclyr = []
    plwclyr = []
    tiwclyr = []
    sst = []
    t2m = []
    wsp = []
    N_0_rain = []
    mu_rain = []
    N_0_snow = []
    rho_snow = []

    print('Done!')
