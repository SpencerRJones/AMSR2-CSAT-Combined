#!/usr/local/anaconda/bin/python


'''

Preprocessor for AMSR2.

1) Combines and extracts CloudSat and AMSR2 coincident observations
2) Adds ERA5 ancillary data
3) Writes output file to be fed into CSU1DVAR combined retrieval

Spencer Jones, CSU Atmos., 11/2023


Flags:

    sfctype: Surface Type Flag
        0: Ocean
        2: Land
        3: Sea Ice

    qflg: Quality Flag
        0: Good
       -1: Too many bad CloudSat observations in kernel averaging window
        1: One or more Tbs missing
        2: One or more EIA missing
        3: Possible sunglint
        4: Possible sea ice (SST < 271)
        5: Non-ocean surface

    modiscf: MODIS Cloud Flag from 2B-GEOPROF
        0: Clear High Confidence
        1: Clear Low Confidence
        2: Cloudy Low Confidence
        3: Cloudy High Confidence

'''

import numpy as np
import xarray as xr
import glob
import os
import h5py as hdf5
import datetime
import sys


#-------------------------------------------------------------------------------------


averaging_length = 30   #number of cloudsat radar observations to along-track average

nlyrs = 30              #number of vertical layers

layer_thickness = 500.  #thickness of layers [m]

csat_geo_dir = '/edata2/spencer/cloudsat/2B-GEOPROF.P1_R05_nc4/'
csat_rp_dir  = '/edata2/spencer/cloudsat/2C-RAIN-PROFILE.P1_R05_nc4/'
csat_sp_dir  = '/edata2/spencer/cloudsat/2C-SNOW-PROFILE.P1_R05_nc4/'
era5_dir     = '/qdata2/archive/ERA5/'
reynolds_dir = '/edata2/spencer/Reynolds_SST/'
rss_dir      = '/edata2/spencer/RSS_AMSR2/'

#-------------------------------------------------------------------------------------

if len(sys.argv) != 3:
    raise ValueError(f'Number of args is wrong.')


infile  = sys.argv[1]
outfile = sys.argv[2]

print(f"Infile: {infile}")
print(f"Outfile: {outfile}")



#----------------------------------------------------------------
#
#   Function definitions
#
#----------------------------------------------------------------

def haversine(lat1, lon1, lat_arr, lon_arr):
    
    R = 6371000.0
    
    phi_1 = np.deg2rad(lat1)
    phi_2 = np.deg2rad(lat_arr)
    
    delta_phi = np.deg2rad(lat_arr - lat1)
    delta_lambda = np.deg2rad(lon_arr - lon1)

    a = np.sin(delta_phi / 2.0) ** 2 + np.cos(phi_1) * np.cos(phi_2) * np.sin(delta_lambda / 2.0) ** 2
    
    if np.any(a>1.0):
        a[np.where(a>1)] = 1.0
    
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    
    meters = R * c 
    km = meters / 1000.0
    
    return km


#-----------------------------------------------------------------

def open_l1c(infile):
    
    nchans = 10
    
    
    with xr.open_dataset(infile, group='S1') as f:
        nscans = f.dims['phony_dim_2']
        npixs  = f.dims['phony_dim_3']
        
        amsr2_Tbs    = np.zeros([nscans,npixs,nchans],dtype='f')
        amsr2_eia    = np.zeros([nscans,npixs,nchans],dtype='f')
        amsr2_snglnt = np.zeros([nscans,npixs,nchans], dtype='f')
        
        amsr2_Tbs[:,:,0:2] = f.Tc[:,:,:]
        amsr2_eia[:,:,0:2] = f.incidenceAngle[:,:,:]
        amsr2_snglnt[:,:,0:2] = f.sunGlintAngle[:,:,:]

        
    with xr.open_dataset(infile, group='S2') as f:
        amsr2_lats = np.zeros([nscans,npixs], dtype='f')
        amsr2_lons = np.zeros([nscans,npixs], dtype='f')
        amsr2_lats[:] = -9999.9
        amsr2_lons[:] = -9999.9
        
        amsr2_lats[:,:]       = f.Latitude[:,:]
        amsr2_lons[:,:]       = f.Longitude[:,:]
        
        amsr2_Tbs[:,:,2:4]    = f.Tc[:,:,:]
        amsr2_eia[:,:,2:4]    = f.incidenceAngle[:,:,:]
        amsr2_snglnt[:,:,2:4] = f.sunGlintAngle[:,:,:]
        
    with xr.open_dataset(infile, group='S3') as f:
        amsr2_Tbs[:,:,4:6]    = f.Tc[:,:,:]
        amsr2_eia[:,:,4:6]    = f.incidenceAngle[:,:,:]
        amsr2_snglnt[:,:,4:6] = f.sunGlintAngle[:,:,:]
            
    with xr.open_dataset(infile, group='S4') as f:
        amsr2_Tbs[:,:,6:8]    = f.Tc[:,:,:]
        amsr2_eia[:,:,6:8]    = f.incidenceAngle[:,:,:]
        amsr2_snglnt[:,:,6:8] = f.sunGlintAngle[:,:,:]
            
    with xr.open_dataset(infile, group='S5') as f:
        amsr2_Tbs[:,:,8:10]    = f.Tc[:,::2,:]
        amsr2_eia[:,:,8:10]    = f.incidenceAngle[:,::2,:]
        amsr2_snglnt[:,:,8:10] = f.sunGlintAngle[:,::2,:]
        
        
    with hdf5.File(infile, 'r') as f:
        scyear = f['S2/ScanTime/Year'][:]
        scmon  = f['S2/ScanTime/Month'][:]
        scday  = f['S2/ScanTime/DayOfMonth'][:]
        schr   = f['S2/ScanTime/Hour'][:]
        scmin  = f['S2/ScanTime/Minute'][:]
        scsec  = f['S2/ScanTime/Second'][:]
        
    sctime = np.zeros([nscans,6],dtype='i')
    
    sctime[:,0] = scyear
    sctime[:,1] = scmon
    sctime[:,2] = scday
    sctime[:,3] = schr
    sctime[:,4] = scmin
    sctime[:,5] = scsec
        
    return nscans, npixs, amsr2_lats, amsr2_lons, amsr2_Tbs, amsr2_eia, amsr2_snglnt, sctime

#---------------------------------------------------------------------------------------------

def get_corr_csat_files(beg_dttm, end_dttm, csat_dir):
    
    no_date = False
    
    #---Gets CloudSat granules near the AMSR2 granule
    
    corr_csat_files = []
    
    #Get next and previous day's CloudSat filenames
    currday   = beg_dttm.date()
    daybefore = (beg_dttm - datetime.timedelta(days=1)).date()
    dayafter  = (beg_dttm + datetime.timedelta(days=1)).date()
    
    cs_dte_bf = str(daybefore)[0:4]+str(daybefore)[5:7]+str(daybefore)[8:]
    cs_dte_af = str(dayafter)[0:4] +str(dayafter)[5:7]+str(dayafter)[8:]
    cs_dte    = str(currday)[0:4]  +str(currday)[5:7] + str(currday)[8:]
    
    cs_flist_dybf = glob.glob(csat_dir + cs_dte_bf[:-2] + '/' + cs_dte_bf + '/*.nc4')
    cs_flist_dyaf = glob.glob(csat_dir + cs_dte_af[:-2] + '/' + cs_dte_af + '/*.nc4')
    cs_flist_curr = glob.glob(csat_dir + cs_dte[:-2]    + '/' + cs_dte    + '/*.nc4')
    
    cs_flist_dybf.sort()
    cs_flist_dyaf.sort()
    cs_flist_curr.sort()
    
    if len(cs_flist_dybf) == 0:
        print('No cloudsat date for ', cs_dte_bf)
        no_date = True
    
    if len(cs_flist_dyaf) == 0:
        print('No cloudsat date for ', cs_dte_af)
        no_date = True
        
    if len(cs_flist_curr) == 0:
        print('No cloudsat date for ', cs_dte)
        no_date = True
        
#    if no_date == True:
#        corr_csat_files = []
#        return corr_csat_files, no_date
        
    
    cs_flist = np.append(cs_flist_dybf, cs_flist_curr)
    cs_flist = np.append(cs_flist, cs_flist_dyaf)
    
    #---Loop through files and get granule times near the AMSR2 file
    for ifile in cs_flist:
        
        fyr  = int(ifile.split('/')[-1].split('_')[0][0:4])
        fdoy = int(ifile.split('/')[-1].split('_')[0][4:7]) - 1
        
        begofyr = datetime.date(fyr,1,1).toordinal()
        
        fdate = datetime.date.fromordinal(begofyr+fdoy)
        
        fmon = fdate.month
        fdy = fdate.day
        
        ftime = ifile.split('/')[-1].split('_')[0][-6:]
        
        fhr  = int(ftime[0:2])
        fmin = int(ftime[2:4])
        fsec = int(ftime[4:])
        
        fdttm = datetime.datetime(fyr,fmon,fdy,fhr,fmin,fsec)
        
        
        if fdttm < beg_dttm - datetime.timedelta(hours=1):
            continue
            
        if fdttm > end_dttm + datetime.timedelta(hours=1):
            continue
            
        corr_csat_files = np.append(corr_csat_files, ifile)
        
    
    return corr_csat_files, no_date


#-----------------------------------------------------------------------------

def get_coinc_pixels(alats, alons, clats, clons, nrays1, nrays2):
    
    alats = alats.T.copy()
    alons = alons.T.copy()
    
    ntotrays = clats.size
    
    nscans = alats.shape[1]
    npixs  = alats.shape[0]
    
    if nrays2 == 0:  #only one file
        print(f'Warning: two CloudSat Files do not exist.')
        return
    
    #---Loop through CloudSat rays and get the first one corresponding with an AMSR2 pixel
    for iray in np.arange(9000, nrays1):
        
        dists = haversine(clats[iray], clons[iray], alats[:,:100], alons[:,:100])
        min_dist = dists.min()
        
        if min_dist < 5.:
            beg_indx = iray
            break
            
    for iray in np.arange(46000, ntotrays):
        
        dists = haversine(clats[iray], clons[iray], alats[:,-200:], alons[:,-200:])
        min_dist = dists.min()
        
        if min_dist > 10:
            end_indx = iray
            break
    
    
    size = end_indx - beg_indx
    
        
    #---Now, loop through the rays, beginning at the first one found previously:
    first = True
    
    ascns = np.zeros(size, dtype='i')
    apixs = np.zeros(size, dtype='i')
    ascns[:] = -9999
    apixs[:] = -9999
    
    
    
    for iray in np.arange(beg_indx, end_indx):
        
        clat = clats[iray]
        clon = clons[iray]
        
        ipix = iray - beg_indx
        
        
        if first: #If this is the first ray
            
            dists = haversine(clat, clon, alats[:,:20], alons[:,:20])
            min_dist = dists.min()
            where_min = np.where(dists == min_dist)
            
            if min_dist > 5.:
                raise ValueError(f'Minimum distance is {min_dist} for first scan.')
                
            
            first_ascn = where_min[1]
            first_apix = where_min[0]
            
            ascns[ipix] = first_ascn
            apixs[ipix] = first_apix
            
            first = False
        
        else:  
            
            #Check first 50 scans for closest pixels to this radar beam
            if ascns[ipix-1] < 50:
                dists = haversine(clat, clon, alats[40:-40,:60], alons[40:-40,:60])
                min_dist = dists.min()
                where_min = np.where(dists == min_dist)
                
                if min_dist > 10:
                    print(f'min_dist is {min_dist} for iray {iray}')
                    raise
                    
                if dists[where_min].size > 1: #Found more than one the same distance away? rare, but happens.
                    where_min = (where_min[0][0], where_min[1][0])
                
                
                ascns[ipix] = where_min[1]
                apixs[ipix] = where_min[0] + 40
                
                
            #After 50, keep checking a running window
            elif ascns[ipix-1] >= 50 and ascns[ipix-1] <= nscans - 50:
                
                last_scan = ascns[ipix-1]
                
                dists = haversine(clat, clon, 
                                       alats[40:-40,last_scan-10:last_scan+10],
                                       alons[40:-40,last_scan-10:last_scan+10])
                
                min_dist = dists.min()
                where_min = np.where(dists == min_dist)
                
                if min_dist > 10:
                    print(f'min_dist is {min_dist}')
                    raise
                    
                if dists[where_min].size > 1:
                    where_min = (where_min[0][0], where_min[1][0])
                    
                    
                ascns[ipix] = where_min[1] + last_scan - 10
                apixs[ipix] = where_min[0] + 40

                
                
            elif ascns[ipix-1] > nscans - 50:
                
                dists = haversine(clat, clon,
                                       alats[40:-40,-50:], alons[40:-40,-50:])
            
                min_dist = dists.min()
                where_min = np.where(dists == min_dist)
                
                if min_dist > 10:
                    print(f'min_dist is {min_dist}')
                    raise
                
                if dists[where_min].size > 1:
                    where_min = (where_min[0][0], where_min[1][0])
                
                ascns[ipix] = where_min[1] + nscans - 50
                apixs[ipix] = where_min[0] + 40
        
    
    
    return beg_indx, end_indx, ascns, apixs


#-------------------------------------------------------------------------------------

def add_ecmwf_surf(pix_sctime, pix_lats, pix_lons):
    
    npixs = pix_lats.size
    
    plons = pix_lons.copy()
    plats = pix_lats.copy()
    ptime = pix_sctime.copy()
    
    global u10
    global v10
    global t2m
    global d2m
    global si
    global sst
    global sp
    global tcwv
    global tciw
    global tclw
    global tp
    global cp
    
    u10  = np.zeros([npixs], dtype='f')
    v10  = np.zeros([npixs], dtype='f')
    t2m  = np.zeros([npixs], dtype='f')
    d2m  = np.zeros([npixs], dtype='f')
    si   = np.zeros([npixs], dtype='f')
    sst  = np.zeros([npixs], dtype='f')
    sp   = np.zeros([npixs], dtype='f')
    tcwv = np.zeros([npixs], dtype='f')
    tciw = np.zeros([npixs], dtype='f')
    tclw = np.zeros([npixs], dtype='f')
    tp   = np.zeros([npixs], dtype='f')
    cp   = np.zeros([npixs], dtype='f')

    u10[:]  = -9999.9
    v10[:]  = -9999.9
    t2m[:]  = -9999.9
    d2m[:]  = -9999.9
    si[:]   = -9999.9
    sst[:]  = -9999.9
    sp[:]   = -9999.9
    tcwv[:] = -9999.9
    tciw[:] = -9999.9
    tclw[:] = -9999.9
    tp[:]   = -9999.9
    cp[:]   = -9999.9
    
    #Get date:
    date = pix_sctime[0,:]
    yr = date[0]
    mon = date[1]
    dy = date[2]
    hr = date[3]
    mn = date[4]
    sec = date[5]
    
    date1 = datetime.datetime(yr, mon, dy, hr, mn, sec)
    date2 = date1 + datetime.timedelta(days=1)
    
    era5_date1 = str(date1.date())
    era5_date1 = f'{era5_date1[:4]}{era5_date1[5:7]}{era5_date1[8:]}'
    era5_date2 = str(date2.date())
    era5_date2 = f'{era5_date2[:4]}{era5_date2[5:7]}{era5_date2[8:]}'
    
    era5_flist = [glob.glob(f'{era5_dir}{era5_date1[:-2]}/*{era5_date1}*surf.nc')[0],
                  glob.glob(f'{era5_dir}{era5_date2[:-2]}/*{era5_date2}*surf.nc')[0]]
    
    with xr.open_dataset(era5_flist[0]) as f:
        elats = f.latitude.values
        elons = f.longitude.values
        etime = f.time.values
        
        nlats = elats.size
        nlons = elons.size
        ntimes = etime.size
        
        u10_in = np.zeros([ntimes+1, nlats, nlons], dtype='f')
        v10_in = np.zeros([ntimes+1, nlats, nlons], dtype='f')
        t2m_in = np.zeros([ntimes+1, nlats, nlons], dtype='f')
        d2m_in = np.zeros([ntimes+1, nlats, nlons], dtype='f')
        si_in  = np.zeros([ntimes+1, nlats, nlons], dtype='f')
        sst_in = np.zeros([ntimes+1, nlats, nlons], dtype='f')
        skt_in = np.zeros([ntimes+1, nlats, nlons], dtype='f')
        sp_in  = np.zeros([ntimes+1, nlats, nlons], dtype='f')
        tcwv_in = np.zeros([ntimes+1, nlats, nlons], dtype='f')
        tciw_in = np.zeros([ntimes+1, nlats, nlons], dtype='f')
        tclw_in = np.zeros([ntimes+1, nlats, nlons], dtype='f')
        tp_in = np.zeros([ntimes+1, nlats, nlons], dtype='f')
        cp_in = np.zeros([ntimes+1, nlats, nlons], dtype='f')
        
        u10_in[:] = -9999.9
        v10_in[:] = -9999.9
        t2m_in[:] = -9999.9
        d2m_in[:] = -9999.9
        si_in[:]  = -9999.9
        sst_in[:] = -9999.9
        skt_in[:] = -9999.9
        sp_in[:]  = -9999.9
        tcwv_in[:] = -9999.9
        tciw_in[:] = -9999.9
        tclw_in[:] = -9999.9
        tp_in[:] = -9999.9
        cp_in[:] = -9999.9
        
        u10_in[:ntimes,:,:]  = f.u10.values
        v10_in[:ntimes,:,:]  = f.v10.values
        t2m_in[:ntimes,:,:]  = f.t2m.values
        d2m_in[:ntimes,:,:]  = f.d2m.values
        si_in[:ntimes,:,:]   = f.siconc.values
        sst_in[:ntimes,:,:]  = f.sst.values
        skt_in[:ntimes,:,:]  = f.skt.values
        sp_in[:ntimes,:,:]   = f.sp.values
        tcwv_in[:ntimes,:,:] = f.tcwv.values
        tciw_in[:ntimes,:,:] = f.tciw.values
        tclw_in[:ntimes,:,:] = f.tclw.values
        tp_in[:ntimes,:,:]   = f.tp.values
        cp_in[:ntimes,:,:]   = f.cp.values
        
    with xr.open_dataset(era5_flist[1]) as f:
        etime = np.append(etime, f.time.values[0])
        u10_in[ntimes,:,:]  = f.u10.values[0,:,:]
        v10_in[ntimes,:,:]  = f.v10.values[0,:,:]
        t2m_in[ntimes,:,:]  = f.t2m.values[0,:,:]
        d2m_in[ntimes,:,:]  = f.d2m.values[0,:,:]
        si_in[ntimes,:,:]   = f.siconc.values[0,:,:]
        sst_in[ntimes,:,:]  = f.sst.values[0,:,:]
        skt_in[ntimes,:,:]  = f.skt.values[0,:,:]
        sp_in[ntimes,:,:]   = f.sp.values[0,:,:]
        tcwv_in[ntimes,:,:] = f.tcwv.values[0,:,:]
        tciw_in[ntimes,:,:] = f.tciw.values[0,:,:]
        tclw_in[ntimes,:,:] = f.tclw.values[0,:,:] 
        tp_in[ntimes,:,:]   = f.tp.values[0,:,:]
        cp_in[ntimes,:,:]   = f.cp.values[0,:,:]

    
    plons[np.where(plons < 0.)] = plons[np.where(plons < 0.)] + 360.
    
    for ipix in np.arange(0,npixs):
        ptime = datetime.datetime(pix_sctime[ipix,0],pix_sctime[ipix,1],pix_sctime[ipix,2],
                                  pix_sctime[ipix,3],pix_sctime[ipix,4],pix_sctime[ipix,5])
        
        ptime = np.datetime64(ptime)
        
        time_indx = np.abs(ptime - etime).argmin()
        lat_indx  = np.abs(plats[ipix] - elats).argmin()
        lon_indx = np.abs(plons[ipix] - elons).argmin()
        
    
        u10[ipix]  = u10_in[time_indx,lat_indx,lon_indx]
        v10[ipix]  = v10_in[time_indx,lat_indx,lon_indx]
        t2m[ipix]  = t2m_in[time_indx,lat_indx,lon_indx]
        d2m[ipix]  = d2m_in[time_indx,lat_indx,lon_indx]
        si[ipix]   = si_in[time_indx,lat_indx,lon_indx]
        sst[ipix]  = sst_in[time_indx,lat_indx,lon_indx]
        sp[ipix]   = sp_in[time_indx,lat_indx,lon_indx]
        tcwv[ipix] = tcwv_in[time_indx,lat_indx,lon_indx]
        tciw[ipix] = tciw_in[time_indx,lat_indx,lon_indx]
        tclw[ipix] = tclw_in[time_indx,lat_indx,lon_indx]
        tp[ipix]   = tp_in[time_indx,lat_indx,lon_indx]
        cp[ipix]   = cp_in[time_indx,lat_indx,lon_indx]
        
    
    return


#-----------------------------------------------------------------------------------------

def add_reynolds_sst(pix_sctime, pix_lats, pix_lons):
    
    npixs = pix_lats.size
    
    plons = pix_lons.copy()
    plats = pix_lats.copy()
    ptime = pix_sctime.copy()
    
    global rey_sst
    
    rey_sst = np.zeros([npixs], dtype='f')
    
    rey_sst[:] = -9999.9
    
    #Get date:
    date = pix_sctime[0,:]
    yr = date[0]
    mon = date[1]
    dy = date[2]
    hr = date[3]
    mn = date[4]
    sec = date[5]
    
    reyn_date = datetime.datetime(yr, mon, dy, hr, mn, sec)
    
    reyn_file = glob.glob(f'{reynolds_dir}sst.day.mean.{reyn_date.year}.nc')[0]
    
    plons[np.where(plons < 0.)] = plons[np.where(plons < 0.)] + 360.
    
    with xr.open_dataset(reyn_file) as f:
        time = f.time.values
        lat  = f.lat.values
        lon  = f.lon.values
        sst  = f.sst.values + 273.15
        
        
    for ipix in np.arange(0, npixs):
        
        pdttm = datetime.datetime(ptime[ipix,0], ptime[ipix,1], ptime[ipix,2], 
                                  ptime[ipix,3], ptime[ipix,4], ptime[ipix,5])
        pdttm = np.datetime64(pdttm)
        
        time_indx = np.abs(pdttm - time).argmin()
        lat_indx  = np.abs(plats[ipix] - lat).argmin()
        lon_indx  = np.abs(plons[ipix] - lon).argmin()
        
        rey_sst[ipix] = sst[time_indx, lat_indx, lon_indx]
            
    
    return


#-----------------------------------------------------------------------------------------


def add_rss_sst(pix_sctime, pix_lats, pix_lons):
    
    npixs = pix_lats.size
    
    plons = pix_lons.copy()
    plats = pix_lats.copy()
    ptime = pix_sctime.copy()
    
    global rss_sst
    
    rss_sst = np.zeros([npixs], dtype='f')
    
    rss_sst[:] = -9999.9
    
    #Get date:
    date = pix_sctime[0,:]
    yr = date[0]
    mon = date[1]
    dy = date[2]
    hr = date[3]
    mn = date[4]
    sec = date[5]
    
    rss_date1 = datetime.datetime(yr, mon, dy, hr, mn, sec)
    rss_date1 = np.datetime64(rss_date1)
    rss_date2 = rss_date1 + np.timedelta64(1, 'D')
    
    try:
        rss_file1 = glob.glob(f"{rss_dir}RSS_AMSR2_ocean_L3_daily_{str(rss_date1)[:10]}_v08.2.nc")[0]
        rss_file2 = glob.glob(f"{rss_dir}RSS_AMSR2_ocean_L3_daily_{str(rss_date2)[:10]}_v08.2.nc")[0]
    except:
        print(f'RSS SST missing for {rss_date1} or {rss_date2}.')
        return
    
    with xr.open_mfdataset([rss_file1, rss_file2], combine='nested', concat_dim='pass') as f:
        lat = f.lat.values
        lon = f.lon.values
        time = f.time.values
        sst  = f.SST.values + 273.15
    
    
    plons[np.where(plons < 0.)] = plons[np.where(plons < 0.)] + 360.
    
    isnat = time.view('i8') == np.datetime64('NaT').view('i8')
    time[np.where(isnat)] = np.datetime64('1970-01-01')
        
        
    for ipix in np.arange(0, npixs):
        
        pdttm = datetime.datetime(ptime[ipix,0], ptime[ipix,1], ptime[ipix,2], 
                                  ptime[ipix,3], ptime[ipix,4], ptime[ipix,5])
        pdttm = np.datetime64(pdttm)
        
        lat_indx = np.abs(plats[ipix] - lat).argmin()
        lon_indx = np.abs(plons[ipix] - lon).argmin()
        
        #Get pass index:
        times = time[:,lat_indx,lon_indx]
        time_diffs = np.abs(pdttm - times)
        time_indx  = time_diffs.argmin()
        if time_diffs[time_indx] > np.timedelta64(1, 'h'):
            #print('RSS pass time missing...')
            continue
        
        rss_sst[ipix] = sst[time_indx, lat_indx, lon_indx]
        
    
    return




#-----------------------------------------------------------------------------------------


def add_ecmwf_prof(pix_sctime, pix_lats, pix_lons):
    
    npixs = pix_lats.size
    
    plons = pix_lons.copy()
    plats = pix_lats.copy()
    ptime = pix_sctime.copy()
    
    nlevs  = 27
    
    global temp
    global hum
    global hgt
    global plev
    global clwc
    global ciwc
    
    temp = np.zeros([npixs,nlevs], dtype='f')
    hum  = np.zeros([npixs,nlevs], dtype='f')
    hgt  = np.zeros([npixs,nlevs], dtype='f')
    plev = np.zeros([nlevs], dtype='f')
    clwc = np.zeros([npixs,nlevs], dtype='f')
    ciwc = np.zeros([npixs,nlevs], dtype='f')
    
    temp[:] = -9999.9
    hum[:]  = -9999.9
    hgt[:]  = -9999.9
    clwc[:] = -9999.9
    ciwc[:] = -9999.9

    
    #Get date:
    date = pix_sctime[0,:]
    yr = date[0]
    mon = date[1]
    dy = date[2]
    hr = date[3]
    mn = date[4]
    sec = date[5]
    
    date1 = datetime.datetime(yr, mon, dy, hr, mn, sec)
    date2 = date1 + datetime.timedelta(days=1)
    
    era5_date1 = str(date1.date())
    era5_date1 = f'{era5_date1[:4]}{era5_date1[5:7]}{era5_date1[8:]}'
    era5_date2 = str(date2.date())
    era5_date2 = f'{era5_date2[:4]}{era5_date2[5:7]}{era5_date2[8:]}'
    
    era5_flist = [glob.glob(f'{era5_dir}{era5_date1[:-2]}/*{era5_date1}*plev.nc')[0],
                  glob.glob(f'{era5_dir}{era5_date2[:-2]}/*{era5_date2}*plev.nc')[0]]
    
    plons[np.where(plons < 0.)] = plons[np.where(plons < 0.)] + 360.
    
    with xr.open_dataset(era5_flist[0]) as f:
        elats = f.latitude.values
        elons = f.longitude.values
        etime = f.time.values
        eplev = f.level.values
        
        nlats  = elats.size
        nlons  = elons.size
        ntimes = etime.size
        nlevs  = eplev.size
        
        plev[:] = eplev
        
        etemp = np.zeros([ntimes+1,nlevs,nlats,nlons], dtype='f')
        ehum  = np.zeros([ntimes+1,nlevs,nlats,nlons], dtype='f')
        ehgt  = np.zeros([ntimes+1,nlevs,nlats,nlons], dtype='f')
        eclwc = np.zeros([ntimes+1,nlevs,nlats,nlons], dtype='f')
        eciwc = np.zeros([ntimes+1,nlevs,nlats,nlons], dtype='f')
        
        etemp[:] = -9999.9
        ehum[:]  = -9999.9
        ehgt[:]  = -9999.9
        eclwc[:] = -9999.9
        eciwc[:] = -9999.9
        
        etemp[:-1,:,:,:] = f.t.values
        ehum[:-1,:,:,:]  = f.q.values * 1000. 
        ehgt[:-1,:,:,:]  = f.z.values / 9.80665
        eclwc[:-1,:,:,:] = f.clwc.values * 1000.
        eciwc[:-1,:,:,:] = f.ciwc.values * 1000.




    with xr.open_dataset(era5_flist[1]) as f:
        etime = np.append(etime, f.time.values[0])
        etemp[-1,:,:,:] = f.t.values[0,:,:,:]
        ehum[-1,:,:,:]  = f.q.values[0,:,:,:] * 1000.
        ehgt[-1,:,:,:]  = f.z.values[0,:,:,:] / 9.80665
        eclwc[-1,:,:,:] = f.clwc.values[0,:,:,:] * 1000.
        eciwc[-1,:,:,:] = f.ciwc.values[0,:,:,:] * 1000.

        
    for ipix in np.arange(0,npixs):
        ptime = datetime.datetime(pix_sctime[ipix,0],pix_sctime[ipix,1],pix_sctime[ipix,2],
                                  pix_sctime[ipix,3],pix_sctime[ipix,4],pix_sctime[ipix,5])
        ptime = np.datetime64(ptime)
        
        time_indx = np.abs(ptime - etime).argmin()
        lat_indx  = np.abs(plats[ipix] - elats).argmin()
        lon_indx  = np.abs(plons[ipix] - elons).argmin()
            
        temp[ipix,:] = etemp[time_indx, :, lat_indx, lon_indx]
        hum[ipix,:]  = ehum[time_indx, :, lat_indx, lon_indx]
        hgt[ipix,:]  = ehgt[time_indx,:, lat_indx, lon_indx]
        clwc[ipix,:] = eclwc[time_indx, :, lat_indx, lon_indx]
        ciwc[ipix,:] = eciwc[time_indx, :, lat_indx, lon_indx]

    
    return


#---------------------------------------------------------------------------------------

#def interpolate_linear(hgt, hgt_radar, data):
    
    #-----------------------------------------------------
    #
    #   Inverse distance weighting interpolation
    #
    #-----------------------------------------------------
    
#    era_hgt = np.flipud(hgt)
#    data    = np.flipud(data)
    
#    j = 0
    
#    interp_data    = np.zeros(hgt_radar.size, dtype='f')
#    interp_data[:] = -9999.9
    
    #loop through radar gate heights
#    for i,ihgt in enumerate(hgt_radar):
        
#        cont = False
        
        #print('---------------------------')
        
        #Get boundaries for height levels
#        up_bound = era_hgt[j]
#        lo_bound = era_hgt[j+1]
        
        #print(up_bound, ihgt, lo_bound)
        
        #if current height is same as upper bound, grab data
#        if ihgt == up_bound:
#            interp_data[i] = data[j]
#            continue
            
        #Same for lower bound
#        elif ihgt == lo_bound:
#            interp_data[i] = data[j+1]
#            j = j + 1
#            continue
            
        #If current height is between two ECMWF heights (this is the case most of the time):
#        else:
            #Check if it is within the bounds
#            within_bounds = (ihgt < up_bound and ihgt > lo_bound)
            
            
            
            #If it is not, advance to next set of bounds...
#            while not within_bounds:

                #print('Bounds changed.')
                #print(up_bound,ihgt,lo_bound)
#                j = j + 1
                #Unless we are at the lowest bound...
#                if j+1 >= hgt.size and ihgt == 0.:
#                    interp_data[i] = data[-1]
#                    cont = True
#                    break
                #Recheck bounds now
#                up_bound = era_hgt[j]
#                lo_bound = era_hgt[j+1]
                #test if within bounds and return to loop
#                within_bounds = (ihgt < up_bound and ihgt > lo_bound)
                
#            if cont:
#                continue
                
            #Now we are within the bounds. Interpolate the data.
            #It is a weighted sum of the distance from the boundaries.
#            w1 = 1./np.abs(ihgt - up_bound)
#            w2 = 1./np.abs(ihgt - lo_bound)
#            interp_data[i] = ((data[j] * w1) + (data[j+1] * w2)) / (w1 + w2)
            
            #print(up_bound,ihgt,lo_bound)

#    return interp_data

def interpolate_linear(old_hgt, new_hgt, data):

    #-----------------------------------------------------
    #
    #   Inverse distance weighting interpolation
    #
    #-----------------------------------------------------

    old_hgt = np.flipud(old_hgt)
    data    = np.flipud(data)

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




#------------------------------------------------------------------------------------

def interpolate_log(old_hgt, new_hgt, data):
    
    old_hgt = np.flipud(old_hgt)
    data    = np.flipud(data)
    
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


#------------------------------------------------------------------------------


# Version 1.0 released by David Romps on September 12, 2017.
# Version 1.1 vectorized lcl.R, released on May 24, 2021.
# 
# When using this code, please cite:
# 
# @article{16lcl,
#   Title   = {Exact expression for the lifting condensation level},
#   Author  = {David M. Romps},
#   Journal = {Journal of the Atmospheric Sciences},
#   Year    = {2017},
#   Month   = dec,
#   Number  = {12},
#   Pages   = {3891--3900},
#   Volume  = {74}
# }
#
# This lcl function returns the height of the lifting condensation level
# (LCL) in meters.  The inputs are:
# - p in Pascals
# - T in Kelvins
# - Exactly one of rh, rhl, and rhs (dimensionless, from 0 to 1):
#    * The value of rh is interpreted to be the relative humidity with
#      respect to liquid water if T >= 273.15 K and with respect to ice if
#      T < 273.15 K. 
#    * The value of rhl is interpreted to be the relative humidity with
#      respect to liquid water
#    * The value of rhs is interpreted to be the relative humidity with
#      respect to ice
# - return_ldl is an optional logical flag.  If true, the lifting deposition
#   level (LDL) is returned instead of the LCL. 
# - return_min_lcl_ldl is an optional logical flag.  If true, the minimum of the
#   LCL and LDL is returned.

def LCL(p,T,rh=None,rhl=None,rhs=None,return_ldl=False,return_min_lcl_ldl=False):

    import math
    import scipy.special

    # Parameters
    Ttrip = 273.16     # K
    ptrip = 611.65     # Pa
    E0v   = 2.3740e6   # J/kg
    E0s   = 0.3337e6   # J/kg
    ggr   = 9.81       # m/s^2
    rgasa = 287.04     # J/kg/K 
    rgasv = 461        # J/kg/K 
    cva   = 719        # J/kg/K
    cvv   = 1418       # J/kg/K 
    cvl   = 4119       # J/kg/K 
    cvs   = 1861       # J/kg/K 
    cpa   = cva + rgasa
    cpv   = cvv + rgasv

    # The saturation vapor pressure over liquid water
    def pvstarl(T):
        return ptrip * (T/Ttrip)**((cpv-cvl)/rgasv) * \
            math.exp( (E0v - (cvv-cvl)*Ttrip) / rgasv * (1/Ttrip - 1/T) )

    # The saturation vapor pressure over solid ice
    def pvstars(T):
        return ptrip * (T/Ttrip)**((cpv-cvs)/rgasv) * \
            math.exp( (E0v + E0s - (cvv-cvs)*Ttrip) / rgasv * (1/Ttrip - 1/T) )

    # Calculate pv from rh, rhl, or rhs
    rh_counter = 0
    if rh  is not None:
        rh_counter = rh_counter + 1
    if rhl is not None:
        rh_counter = rh_counter + 1
    if rhs is not None:
        rh_counter = rh_counter + 1
    if rh_counter != 1:
        print(rh_counter)
        exit('Error in lcl: Exactly one of rh, rhl, and rhs must be specified')
    if rh is not None:
      # The variable rh is assumed to be 
      # with respect to liquid if T > Ttrip and 
      # with respect to solid if T < Ttrip
        if T > Ttrip:
            pv = rh * pvstarl(T)
        else:
            pv = rh * pvstars(T)
        rhl = pv / pvstarl(T)
        rhs = pv / pvstars(T)
    elif rhl is not None:
        pv = rhl * pvstarl(T)
        rhs = pv / pvstars(T)
        if T > Ttrip:
            rh = rhl
        else:
            rh = rhs
    elif rhs is not None:
        pv = rhs * pvstars(T)
        rhl = pv / pvstarl(T)
        if T > Ttrip:
            rh = rhl
        else:
            rh = rhs
    if pv > p:
        return NA

    # Calculate lcl_liquid and lcl_solid
    qv = rgasa*pv / (rgasv*p + (rgasa-rgasv)*pv)
    rgasm = (1-qv)*rgasa + qv*rgasv
    cpm = (1-qv)*cpa + qv*cpv
    if rh == 0:
        return cpm*T/ggr
    aL = -(cpv-cvl)/rgasv + cpm/rgasm
    bL = -(E0v-(cvv-cvl)*Ttrip)/(rgasv*T)
    cL = pv/pvstarl(T)*math.exp(-(E0v-(cvv-cvl)*Ttrip)/(rgasv*T))
    aS = -(cpv-cvs)/rgasv + cpm/rgasm
    bS = -(E0v+E0s-(cvv-cvs)*Ttrip)/(rgasv*T)
    cS = pv/pvstars(T)*math.exp(-(E0v+E0s-(cvv-cvs)*Ttrip)/(rgasv*T))
    lcl = cpm*T/ggr*( 1 - \
      bL/(aL*scipy.special.lambertw(bL/aL*cL**(1/aL),-1).real) )
    ldl = cpm*T/ggr*( 1 - \
      bS/(aS*scipy.special.lambertw(bS/aS*cS**(1/aS),-1).real) )

    # Return either lcl or ldl
    if return_ldl and return_min_lcl_ldl:
        exit('return_ldl and return_min_lcl_ldl cannot both be true')
    elif return_ldl:
        return ldl
    elif return_min_lcl_ldl:
        return min(lcl,ldl)
    else:
        return lcl
    
    
#---------------------------------------------------------------------------------

def RH(T, T_d):
    
    T_C = T - 273.15
    T_dC = T_d - 273.15
    
    e_s = 611.2*np.exp((17.67*T_C)/(T_C + 243.5))
    e   = 611.2*np.exp((17.67*T_dC)/(T_dC + 243.5))
    
    return e/e_s


#----------------------------------------------------------------------------------


def output(npixs_out, nlyrs_out, hgt_out, sctime_out, lats_out, 
           lons_out, Tbs_out, eia_out,
           refl_out, tcwv_out, tclw_out, tciw_out, sst_out, 
           rey_sst_out, rss_sst_out, t2m_out, d2m_out,
           wind_out, slp_out, tprof_out, pres_out, 
           hum_out, fl_out, fl_bin_out, cldbse_out, cldbse_bin_out, lcl_out,
           qflg_out, sfctype_out,
           modiscf_out, csatrr_out, csatsr_out):
    
    
    
    with open(outfile, 'wb') as f:
        npixs_out.tofile(f, sep='', format='unformatted')
        nlyrs_out.tofile(f, sep='', format='unformatted')
        hgt_out.tofile(f, sep='', format='unformatted')
        sctime_out.tofile(f, sep='', format='unformatted')
        lats_out.tofile(f, sep='', format='unformatted')
        lons_out.tofile(f, sep='', format='unformatted')
        Tbs_out.tofile(f, sep='', format='unformatted')
        eia_out.tofile(f, sep='', format='unformatted')
        refl_out.tofile(f, sep='', format='unformatted')
        tcwv_out.tofile(f, sep='', format='unformatted')
        tclw_out.tofile(f, sep='', format='unformatted')
        tciw_out.tofile(f, sep='', format='unformatted')
        sst_out.tofile(f, sep='', format='unformatted')
        rey_sst_out.tofile(f, sep='', format='unformatted')
        rss_sst_out.tofile(f, sep='', format='unformatted')
        t2m_out.tofile(f, sep='', format='unformatted')
        d2m_out.tofile(f, sep='', format='unformatted')
        wind_out.tofile(f, sep='', format='unformatted')
        slp_out.tofile(f, sep='', format='unformatted')
        tprof_out.tofile(f, sep='', format='unformatted')
        pres_out.tofile(f, sep='', format='unformatted')
        hum_out.tofile(f, sep='', format='unformatted')
        fl_out.tofile(f, sep='', format='unformatted')
        fl_bin_out.tofile(f, sep='', format='unformatted')
        cldbse_out.tofile(f, sep='', format='unformatted')
        cldbse_bin_out.tofile(f, sep='', format='unformatted')
        lcl_out.tofile(f, sep='', format='unformatted')
        qflg_out.tofile(f, sep='', format='unformatted')
        sfctype_out.tofile(f, sep='', format='unformatted')
        modiscf_out.tofile(f, sep='', format='unformatted')
        csatrr_out.tofile(f, sep='', format='unformatted')
        csatsr_out.tofile(f, sep='', format='unformatted')
    
    return


#======================================================================================




#---Create radar range gate height array
hgt_radar = np.arange(15000.,0.-layer_thickness,-layer_thickness)


#---Read in L1C file
print('Reading L1C...')
nscans, npixs, alats_in, alons_in, Tbs_in, eia_in, sunglint_in, sctime_in = open_l1c(infile)


#---Get Coincident CloudSat files
fdate = infile.split('/')[-1].split('.')[4].split('-')[0]
fdate = np.datetime64(f'{fdate[0:4]}-{fdate[4:6]}-{fdate[6:]}')

gran_beg = sctime_in[0,:]
gran_beg = datetime.datetime(gran_beg[0],gran_beg[1],gran_beg[2], 
                             gran_beg[3], gran_beg[4], gran_beg[5])
gran_end = sctime_in[-1,:]
gran_end = datetime.datetime(gran_end[0],gran_end[1],gran_end[2], 
                             gran_end[3], gran_end[4], gran_end[5])

print('Getting corresponding CloudSat data files...')
flist_csat_geo, nodate = get_corr_csat_files(gran_beg, gran_end, csat_geo_dir)
flist_csat_rp,  nodate = get_corr_csat_files(gran_beg, gran_end, csat_rp_dir)
flist_csat_sp,  nodate = get_corr_csat_files(gran_beg, gran_end, csat_sp_dir)

print('Corresponding CloudSat files:')
print(flist_csat_geo)
print(flist_csat_rp)
print(flist_csat_sp)

if len(flist_csat_geo) < 2 or len(flist_csat_rp) < 2:
    print('Error. Not enough Cloudsat files around.')
    exit()

#---Get coincident observations
#AMSR2 granules start at southernmost point in orbit, while
#CloudSat granules start at equator crossing. Therefore,
#we need two files
print('Finding coincident AMSR2 and CloudSat data...')
with xr.open_dataset(flist_csat_geo[0]) as f:
    nrays1 = f.dims['rays']
    ctime1 = f.profile_time.values

with xr.open_dataset(flist_csat_geo[1]) as f:
    nrays2 = f.dims['rays']
    ctime2 = f.profile_time.values

with xr.open_mfdataset(flist_csat_geo, combine='nested', concat_dim='rays') as f:
    clats_in = f.latitude.values
    clons_in = f.longitude.values
    
beg_indx, end_indx, coinc_scans, coinc_pixs = get_coinc_pixels(alats_in, alons_in, 
                                                                clats_in, clons_in, 
                                                                nrays1, nrays2)


size = coinc_scans.size
print(f'Total number of rays found: {size}')


coinc_rays = np.zeros([size], dtype='i')
first_chunk = np.arange(beg_indx, nrays1)
second_chunk = np.arange(0,coinc_rays.size-first_chunk.size)
coinc_rays[:first_chunk.size] = first_chunk
coinc_rays[first_chunk.size:] = second_chunk


#Grab coincident data:
amsr2_sctime = sctime_in[coinc_scans, :]
csat_lats    = clats_in[beg_indx:end_indx]
csat_lons    = clons_in[beg_indx:end_indx]
amsr2_lats   = alats_in[coinc_scans, coinc_pixs]
amsr2_lons   = alons_in[coinc_scans, coinc_pixs]
amsr2_Tbs    = Tbs_in[coinc_scans, coinc_pixs, :]
amsr2_eia    = eia_in[coinc_scans, coinc_pixs, :]
amsr2_snglnt = sunglint_in[coinc_scans, coinc_pixs, :]

with xr.open_mfdataset(flist_csat_geo, combine='nested', concat_dim='rays') as f:
    csat_refl = f.radar_reflectivity.values[beg_indx:end_indx] / 100.
    csat_hgt  = f.height.values[beg_indx:end_indx]
    csat_dqflg = f.data_quality.values[beg_indx:end_indx]
    csat_lsflg = f.Navigation_land_sea_flag[beg_indx:end_indx]
    mod_clflg = f.modis_cloud_flag[beg_indx:end_indx]

with xr.open_mfdataset(flist_csat_rp, combine='nested', concat_dim='rays') as f:
    csat_rr = f.rain_rate.values[beg_indx:end_indx]

with xr.open_mfdataset(flist_csat_sp, combine='nested', concat_dim='rays') as f:
    csat_snrate = f.snowfall_rate.values[beg_indx:end_indx]
    csat_snsfc  = f.snowfall_rate_sfc.values[beg_indx:end_indx]

print(f'Getting AMSR2 pixel locations...')
nrays = size

for iray in np.arange(0,nrays):
    if iray == 0:
        pix_indcs = [iray]
        j = 0
    else:
        lat_to_check = amsr2_lats[pix_indcs[j]]
        lon_to_check = amsr2_lons[pix_indcs[j]]
        cur_lat      = amsr2_lats[iray]
        cur_lon      = amsr2_lons[iray]
        
        if cur_lon != lon_to_check and cur_lat != lat_to_check:
            pix_indcs = np.append(pix_indcs, iray)
            j += 1
        else:
            continue

npixs = pix_indcs.size

pix_sctime = amsr2_sctime[pix_indcs,:]
pix_lats   = amsr2_lats[pix_indcs]
pix_lons   = amsr2_lons[pix_indcs]
pix_Tbs    = amsr2_Tbs[pix_indcs,:]
pix_eia    = amsr2_eia[pix_indcs,:]
pix_snglnt = amsr2_snglnt[pix_indcs,:]


#---Add ERA5 surface data to pixel locations:
print('Adding ERA5 surface data...')
add_ecmwf_surf(pix_sctime, pix_lats, pix_lons)


#---Add Reynolds SST:
print('Adding Reynolds SST...')
add_reynolds_sst(pix_sctime, pix_lats, pix_lons)

rey_sst[np.where(np.isnan(rey_sst))] = -9999.9

#---Add RSS SST:
print('Adding RSS SST...')
add_rss_sst(pix_sctime, pix_lats, pix_lons)

rss_sst[np.where(np.isnan(rss_sst))] = -9999.9


#---Add ERA5 profile data:
print('Adding ERA5 profile data...')
add_ecmwf_prof(pix_sctime, pix_lats, pix_lons)

#Get CloudSat beam locations that correspond to AMSR2 pixels
print('Getting center beams...')
ray_centers = []
for ipix in np.arange(0,npixs):
    plat = pix_lats[ipix]
    plon = pix_lons[ipix]
    
    dists    = haversine(plat,plon,csat_lats,csat_lons)
    min_dist = dists.argmin()
    
    if dists.min() > 6.:
        print(f'Warning: minimum distance found from pixel {ipix} was {dists.min()}.')
    
    corr_beam_indx = min_dist
    
    #Create array of indices for these center beams
    ray_centers = np.append(ray_centers, corr_beam_indx)
    
ray_centers = ray_centers.astype('i')    




qflg     = np.zeros([npixs], dtype='i')
sfc_type = np.zeros([npixs], dtype='i')
qflg[:]     = -99
sfc_type[:] = -99


    
    
#Loop through and kernel average reflectivities from center of pixel
print('Calculating mean reflectivity profiles...')
pix_refl = np.zeros([npixs,125])
pix_hgt  = np.zeros([npixs,125])
pix_rr   = np.zeros([npixs])
pix_mcf  = np.zeros([npixs,averaging_length+1])
pix_snsfc = np.zeros([npixs])
pix_rr[:] = -9999.9
pix_mcf[:] = -99
pix_snsfc[:] = -9999.9

for i,indx in enumerate(ray_centers):
    
    center = indx
    beg    = center - (averaging_length//2)
    end    = center + (averaging_length//2) + 1

    pix_hgt[i,:]  = csat_hgt[center]
    
    #Check pixel for missing sst (land)
    if np.isnan(sst[i]):
        pix_refl[i,:] = -9999.9
        continue
    #Check pixel for sea ice:
    if si[i] > 0.:
        pix_refl[i,:] = -9999.9
        continue
        
    
    refl_window = csat_refl[beg:end]
    rr_window   = csat_rr[beg:end]
    q_window    = csat_dqflg[beg:end]
    ls_window   = csat_lsflg[beg:end]
    snsfc_window = csat_snsfc[beg:end]
    
    #---Check for too many bad cloudsat profiles:
    good  = np.logical_and(ls_window==2, q_window==0)
    ngood = np.where(good)[0].size
    nbad  = q_window.size - ngood
    
    if nbad >= averaging_length/4:
        pix_refl[i,:] = -9999.9
        qflg[i] = -1  #Too many bad csat profiles
        continue
    
    #---Along-track average reflectivities:
    pix_refl[i,:] = 10.*np.log10(np.mean(10.**(refl_window/10.), axis=0))
    pix_rr[i]     = np.mean(rr_window)
    pix_snsfc[i]  = np.mean(snsfc_window)
    
    #---Save cloud flag:
    pix_mcf[i,:] = mod_clflg[beg:end]

#---Interpolate data into radar range gates:
print('Interpolating data into radar range gates...')
pix_pres = np.zeros([npixs, nlyrs+1], dtype='f')
pix_temp = np.zeros([npixs, nlyrs+1], dtype='f')
pix_hum  = np.zeros([npixs, nlyrs+1], dtype='f')
pix_clwc = np.zeros([npixs, nlyrs+1], dtype='f')
pix_ciwc = np.zeros([npixs, nlyrs+1], dtype='f')

for ipix in np.arange(0,npixs):
    pix_pres[ipix,:] = interpolate_log(hgt[ipix,:], hgt_radar, plev)
    pix_temp[ipix,:] = interpolate_linear(hgt[ipix,:], hgt_radar, temp[ipix,:])
    pix_hum[ipix,:]  = interpolate_log(hgt[ipix,:], hgt_radar, hum[ipix,:])
    pix_clwc[ipix,:] = interpolate_linear(hgt[ipix,:], hgt_radar, clwc[ipix,:])
    pix_ciwc[ipix,:] = interpolate_linear(hgt[ipix,:], hgt_radar, ciwc[ipix,:])
    
#---Set surface values:
pix_pres[:,-1] = sp/100.
pix_temp[:,-1] = t2m
    
    
#---Reduce resolution 
print(f'Reducing reflectivities into {nlyrs} layers...')
pix_refl_out = np.zeros([npixs, nlyrs], dtype='f')

for ipix in np.arange(0,npixs):
    
    refl_prof = pix_refl[ipix,:]
    hgt_prof  = pix_hgt[ipix,:]
    
    if np.all(refl_prof == -9999.9):
        pix_refl_out[ipix,:] = -9999.9
        continue
        
    for ilyr in np.arange(0,nlyrs):
        lyr_top = hgt_radar[ilyr]
        lyr_bot = hgt_radar[ilyr+1]
    
        within_window = np.where(np.logical_and(hgt_prof <= lyr_top, hgt_prof > lyr_bot))
        
        refl_in_lyr = refl_prof[within_window]
        
        lyr_mean_refl = 10. * np.log10(np.mean((10.**(refl_in_lyr/10.))))

        pix_refl_out[ipix,ilyr] = lyr_mean_refl


        
#---Do final quality check:
print('Creating quality and surface type flags...')


for ipix in np.arange(0,npixs):
    
    #Check for non-ocean surface:
    if np.isnan(sst[ipix]):
        sfc_type[ipix] = 2 #land
    elif si[ipix] > 0.:
        sfc_type[ipix] = 3 #sea ice contamination
    else:
        sfc_type[ipix] = 0 #Good ocean
        
    #Check for bad quality:
    if np.any(pix_Tbs[ipix]) < 0.:
        qflg[ipix] = 1 #Bad Tbs
        continue

    if np.any(pix_eia[ipix]) < 0.:
        qflg[ipix] = 2 #Bad eia
        continue

    if pix_snglnt[ipix,0] > 0. and pix_snglnt[ipix,0] < 20.:
        qflg[ipix] = 3 #Possible sunglint
        continue

    if sst[ipix] < 271.:
        qflg[ipix] = 4 #Possible sea ice
        continue

    if sfc_type[ipix] != 0:
        qflg[ipix] = 5 #Non-ocean surface
        continue
    
    if qflg[ipix] == -1:
        qflg[ipix] = -1 #Bad quality in csat data
        continue



    qflg[ipix] = 0 #Good quality



    #Make all negative rain rates missing:
    #if pix_rr[ipix] < 0:
    #    pix_rr[ipix] = -9999.9

    #Some negative rain rates are valid, so use -5 for the cutoff:
    if pix_rr[ipix] < 0. and pix_rr[ipix] > -5:
        pix_rr[ipix] = np.abs(pix_rr[ipix])



#---Calculate freezing level and get freezing level bin:
print('Calculating freezing level...')
frzlvl_bin = np.zeros([npixs], dtype='i')
frzlvl     = np.zeros([npixs], dtype='f')
frzlvl_bin[:] = -99
frzlvl[:]     = -9999.9


for ipix in np.arange(0,npixs):
    
    if sfc_type[ipix] != 0:
        continue
    
    frzlvl_bin[ipix] = np.where(pix_temp[ipix,:]<273.15)[0][-1] + 1
    if frzlvl_bin[ipix] == nlyrs + 1:
        frzlvl[ipix] = 0.
        continue
    temp_below = pix_temp[ipix,frzlvl_bin[ipix]-1]
    temp_above = pix_temp[ipix,frzlvl_bin[ipix]]
    
    if temp_below > 273.15:
        raise ValueError(f'temp_below {temp_below} is greater than 273.15')
    if temp_above < 273.15:
        raise ValueError(f'temp_above {temp_above} is less than 273.15')
    
    w1 = 1./np.abs(273.15 - temp_below)
    w2 = 1./np.abs(273.15 - temp_above)
    
    frzlvl[ipix] = (((hgt_radar[frzlvl_bin[ipix]-1] * w1) + 
                    (hgt_radar[frzlvl_bin[ipix]] * w2)) / (w1 + w2))


#---Calculate cloud base:
print('Calculating cloud base...')
cldbse     = np.zeros([npixs], dtype='f')
cldbse_bin = np.zeros([npixs], dtype='i')
lcl        = np.zeros([npixs], dtype='f')
cldbse[:]     = -9999.9
cldbse_bin[:] = -99
lcl[:]        = -9999.9

for ipix in np.arange(0,npixs):

    if sfc_type[ipix] != 0:
        continue
        
    if np.all(clwc[ipix,:] <= 0.001) and np.all(ciwc[ipix,:] <= 0.001):
        continue

    clw = clwc[ipix,:]
    ciw = ciwc[ipix,:]
    h   = hgt[ipix,:]
    
    cldbse[ipix]     = h[np.where(np.logical_or(clw >= 0.001, ciw >= 0.001))][0] / 1000. #km
    cldbse_bin[ipix] = np.abs(cldbse[ipix]*1000. - hgt_radar).argmin()
    lcl[ipix]        = LCL(p=sp[ipix],T=t2m[ipix],rh=RH(t2m[ipix],d2m[ipix]),
                           rhl=None,rhs=None,return_ldl=False,return_min_lcl_ldl=False) / 1000. #km
    





#---Write output data:
print('Writing output file...')

npixs_out = np.array(npixs, dtype='i')
nlyrs_out = np.array(nlyrs, dtype='i')
hgt_out  = hgt_radar.astype('f')
sctime_out = pix_sctime.T.astype('i')
lats_out = pix_lats.astype('f')
lons_out = pix_lons.astype('f')
Tbs_out  = pix_Tbs.T.astype('f')
eia_out  = pix_eia.T.astype('f')
snglnt_out = pix_snglnt.T.astype('f')
refl_out = pix_refl_out.T.astype('f')
tcwv_out = tcwv.astype('f')
tclw_out = tclw.astype('f')
tciw_out = tciw.astype('f')
sst_out = sst.astype('f')
rey_sst_out = rey_sst.astype('f')
rss_sst_out = rss_sst.astype('f')
t2m_out = t2m.astype('f')
d2m_out = d2m.astype('f')
wind_out = np.sqrt(u10**2 + v10**2).astype('f')
slp_out = (sp / 100.).astype('f') #hPa
tprof_out = pix_temp.T.astype('f')
pres_out = pix_pres.T.astype('f')
hum_out = pix_hum.T.astype('f')
fl_out = frzlvl.astype('f')
fl_bin_out = frzlvl_bin.astype('i')
cldbse_out = cldbse.astype('f') ###
cldbse_bin_out = cldbse_bin.astype('i') ###
lcl_out = lcl.astype('f') ###
qflg_out = qflg.astype('i')
sfctype_out = sfc_type.astype('i')
modiscf_out = pix_mcf.T.astype('i')
csatrr_out = pix_rr.astype('f')
csatsr_out = pix_snsfc.astype('f')


output(npixs_out, nlyrs_out, hgt_out, sctime_out, lats_out, 
           lons_out, Tbs_out, eia_out,
           refl_out, tcwv_out, tclw_out, tciw_out, sst_out, 
           rey_sst_out, rss_sst_out, t2m_out, d2m_out,
           wind_out, slp_out, tprof_out, pres_out, 
           hum_out, fl_out, fl_bin_out, cldbse_out, cldbse_bin_out, lcl_out, 
           qflg_out, sfctype_out, 
           modiscf_out, csatrr_out, csatsr_out)


print('Done!')


#====================================================================================


exit()
