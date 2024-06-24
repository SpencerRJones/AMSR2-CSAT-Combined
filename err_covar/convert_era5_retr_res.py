#----------------------------------------------------
#
#   RETRIEVAL RESOLUTION
#
#----------------------------------------------------

import numpy as np
import glob

datadir = '/edata2/spencer/thesis_data/era5_ocean_profiles/err_covar_amsr2_csat/'

beg_date = '2015-10-26'
end_date = '2015-10-31'

beg_date = np.datetime64(beg_date)
end_date = np.datetime64(end_date)

flist = []

for idate in np.arange(beg_date, end_date+np.timedelta64(1 ,"D"), np.timedelta64(1, "D")):
    datestr = f"{str(idate).replace('-', '')}"
    ifile = glob.glob(f'{datadir}*_{datestr}_native.bin')
    flist = np.append(flist, ifile)

#flist = glob.glob(f'{datadir}*_native.bin')
#flist.sort()


for ifile in flist:
    print(ifile)
    
    #---Read file:
    with open(ifile, 'rb') as f:
        nprofs = np.fromfile(f, sep='', count=1, dtype='i')[0]
        nhgts  = np.fromfile(f, sep='', count=1, dtype='i')[0]
        nlyrs  = nhgts - 1
        lat    = np.fromfile(f, sep='', count=nprofs, dtype='f')
        lon    = np.fromfile(f, sep='', count=nprofs, dtype='f')
        height = np.fromfile(f, sep='', count=nprofs*nhgts, dtype='f').reshape(nhgts,nprofs).T
        temp   = np.fromfile(f, sep='', count=nprofs*nhgts, dtype='f').reshape(nhgts,nprofs).T
        pres   = np.fromfile(f, sep='', count=nprofs*nhgts, dtype='f').reshape(nhgts,nprofs).T
        qlyr   = np.fromfile(f, sep='', count=nprofs*nlyrs, dtype='f').reshape(nlyrs,nprofs).T
        clwc   = np.fromfile(f, sep='', count=nprofs*nlyrs, dtype='f').reshape(nlyrs,nprofs).T
        plwc   = np.fromfile(f, sep='', count=nprofs*nlyrs, dtype='f').reshape(nlyrs,nprofs).T
        tiwc   = np.fromfile(f, sep='', count=nprofs*nlyrs, dtype='f').reshape(nlyrs,nprofs).T
        sst    = np.fromfile(f, sep='', count=nprofs, dtype='f')
        t2m    = np.fromfile(f, sep='', count=nprofs, dtype='f')
        wsp    = np.fromfile(f, sep='', count=nprofs, dtype='f')
        N_0_rain = np.fromfile(f, sep='', count=nprofs, dtype='f')
        mu_rain  = np.fromfile(f, sep='', count=nprofs, dtype='f')
        N_0_snow = np.fromfile(f, sep='', count=nprofs, dtype='f')
        rho_snow = np.fromfile(f, sep='', count=nprofs, dtype='f')
        
        
    #---Add noise to temperature profile to simulate departure from real:
    #sigma = 2 K
    print('Adding noise to temperature profiles...')
    temp_r   = temp.copy()
    eps_temp = np.random.normal(scale=np.sqrt(2.0), size=temp.shape)
    temp_r   = temp + eps_temp
    
    
    #---Reconstruct humidity from EOF1:
    print('Reconstructing humidity profiles...')
    sst_bin                        = sst.round().astype('i') - 271
    sst_bin[np.where(sst_bin < 0)] = 0
    
    eof_file = '/edata2/spencer/OE/amsr2_radar/OE_amsr2_noradar/binary/eof_mr.03.v3_30lyrs.bin'

    nbins = 33
    npc   = 1

    with open(eof_file, 'rb') as f:
        mprof = np.fromfile(f, sep='', count=nlyrs*nbins, dtype='f').reshape(nlyrs,nbins).T
        eofs  = np.fromfile(f, sep='', count=nbins*nlyrs*6, dtype='f').reshape(6,nlyrs,nbins).T

    sst_bin                         = sst.astype('i') - 271
    sst_bin[np.where(sst_bin < 0)]  = 0
    sst_bin[np.where(sst_bin > 32)] = 32

    X_raw = qlyr.copy()

    X    = np.zeros([nlyrs])
    V    = np.zeros([nlyrs,npc])
    X_re = np.zeros(np.shape(X_raw))

    eof_coeffs = np.arange(-25.0, 25.0, 0.1)

    for iprof in np.arange(0, nprofs):
        
        if iprof%10000 == 0:
            print(iprof)
        
        q = np.flipud(X_raw[iprof]*0.001)
        
        dp = np.flipud(np.diff(pres[iprof]*100.))
        
        mean_profile = mprof[sst_bin[iprof]]
    
        V[:,0] = eofs[sst_bin[iprof],:,0]
    
        X[:] = X_raw[iprof,:] - mean_profile
    
        tpw = (1./(9.80665*1000.)) * np.sum(q * dp)
        
        mp         = np.tile(mean_profile, (eof_coeffs.size,1))
        zs         = np.multiply(V, eof_coeffs).T
        reconstr   = np.add(mp, zs)
        
        reconstr[np.where(reconstr < 0.)] = 0.
        
        tpw_r = (1./(9.80665*1000.)) * np.sum((np.flip(reconstr, axis=1)*0.001 * dp), axis=1)
        
        best = np.abs(tpw - tpw_r).argmin()
        
        Z = np.array([eof_coeffs[best]])
        
        X_hat = np.matmul(Z, V.T)
        x_re = X_hat + mean_profile
        x_re[np.where(x_re < 0.)] = 0.
        X_re[iprof,:] = x_re
        
        
    qlyr_r = X_re.copy()
    qlyr_r[np.where(qlyr_r < 0.)] = 0.
    
    
    #---Cloudwater is integrated and redistributed evenly from the cloud base to the freezing level:
    print('Redistributing cloud water...')
    clwc_r  = np.zeros(clwc.shape)
    dhgt    = np.zeros([nlyrs])
    dhgt[:] = 500.

    for iprof in np.arange(0, nprofs):
        cloudwater = clwc[iprof].copy()
        t          = temp[iprof]

        if np.all(cloudwater < 1.0e-06):
            continue
            
        if np.any(np.isnan(cloudwater)):
            sst[iprof] = np.nan
            continue


        clw = np.sum(cloudwater*dhgt) #[g m^-2]


        frzl         = np.where(t < 273.15)[0][-1]
        cloudtop_bin = frzl - 1
        cloudbot_bin = np.where(cloudwater > 0.)[0][-1]

        if cloudtop_bin > cloudbot_bin:
            cloudbot_bin = cloudtop_bin

        ncloud_bins = (cloudbot_bin - cloudtop_bin) + 1

        clw_lyr = clw / ncloud_bins / 500.

        cloudwater[:] = 0.
        cloudwater[cloudtop_bin:cloudbot_bin+1] = clw_lyr
        
        clwc_r[iprof,:] = cloudwater
        
    #---Rain water is remodeled using noise added to lwc and fixed parameters:
    print('Remodeling rain water...')
    #gaussian noise added to plwc to emulate effect of averaging and interpolation
    plwc_r     = plwc.copy()
    perc_noise = np.random.normal(scale=0.03, size=[nprofs,nlyrs])
    plwc_r     = plwc_r + (perc_noise * plwc_r)
    
    #N_0 is fixed in retrieval, therefore we assume that we get N_0 wrong randomly
    N_0_rain_r = np.random.uniform(low=1.1e+04, high=1.1e+06, size=nprofs)
    
    #---Snow water is remodeled by adding noise to iwc and randomly varying fixed parameters:
    print('Remodeling snow water...')
    tiwc_r     = tiwc.copy()
    perc_noise = np.random.normal(scale=0.03, size=[nprofs,nlyrs])
    tiwc_r     = tiwc_r + (perc_noise * tiwc_r)
    
    #Allow N_0_snow to vary by one order of magnitude
    N_0_snow_r = np.random.uniform(low=5.1e+02, high=5.1e+04, size=nprofs)
    
    #Allow snow density to vary randomly
    rho_snow_r = np.random.uniform(low=50., high=400., size=nprofs)
    
    print('Writing output file...')
    
    outfile = f'{ifile[:-10]}retr.bin'
    
    print(outfile)
    
    pres_r = pres.copy()
    
    with open(outfile, 'wb') as f:
        np.array(nprofs).astype('i').tofile(f, sep='', format='unformatted')
        np.array(nhgts).astype('i').tofile(f, sep='', format='unformatted')
        lat.astype('f').tofile(f, sep='', format='unformatted')
        lon.astype('f').tofile(f, sep='', format='unformatted')
        height.T.astype('f').tofile(f, sep='', format='unformatted')
        temp_r.T.astype('f').tofile(f, sep='', format='unformatted')
        pres_r.T.astype('f').tofile(f, sep='', format='unformatted')
        qlyr_r.T.astype('f').tofile(f, sep='', format='unformatted')
        clwc_r.T.astype('f').tofile(f, sep='', format='unformatted')
        plwc_r.T.astype('f').tofile(f, sep='', format='unformatted')
        tiwc_r.T.astype('f').tofile(f, sep='', format='unformatted')
        sst.astype('f').tofile(f, sep='', format='unformatted')
        t2m.astype('f').tofile(f, sep='', format='unformatted')
        wsp.astype('f').tofile(f, sep='', format='unformatted')
        N_0_rain_r.astype('f').tofile(f, sep='', format='unformatted')
        mu_rain.astype('f').tofile(f, sep='', format='unformatted')
        N_0_snow_r.astype('f').tofile(f, sep='', format='unformatted')
        rho_snow_r.astype('f').tofile(f, sep='', format='unformatted')
        
    print('Done!')
    print('-------------------------------------------------')
