      MODULE optest

      USE define_csu1dvar
      USE dsd
      USE subr_csu1dvar

      implicit  none

      real :: wspd, wdir, slp, salinity, tavg
      real :: plev(nz+1)
      real :: tlev(nz+1)
      real :: hlev(nz+1)
      real :: rmix(nz)
      real :: dp(nz)
      real :: frzl
      integer :: frzl_bin

      real, allocatable :: tbsim(:)
      REAL, ALLOCATABLE :: sim_obs(:)
      real, allocatable :: oe_tbs(:)
      REAL, ALLOCATABLE :: oe_refl(:)
      REAL, ALLOCATABLE :: oe_obs(:)
      INTEGER :: ct_bin, cb_bin
      real :: prcp_rate(nz)
      real :: rn_rate(nz)
      real :: sn_rate(nz)

      integer :: sy_flag

!----------------------------------------------------------------------
!     non-precipitating 1DVAR retrieval algorithm for microwave imagers
!      
!     D Duncan, CSU, Dec 2016 
!     W Berg, CSU,   Apr 2017 : Replaced CRTM with MonoRTM lookup table
!     S Jones, CSU,  Dec 2023 : Adapted to minimize total observation vector
!                               CPR and AMSR2 simultaneously
!----------------------------------------------------------------------

      contains

      subroutine retrieval_nr
       
      real         :: oe_lat, oe_lon, oe_sst, oe_wspd
      real         :: worst
      integer*2    :: lodex,ladex, eclo,ecla,pma,pmo !index for lon,lat
      real         :: eclod,eclad ! diff from grid midpoint and lat/lon geo
      integer      :: lev,z
      logical      :: flag

      ! Arrays for Apriori Parameters and associated errors

      real         :: xa(nvar),  noretrieval 
      real         :: xmin(nretvar), xmax(nretvar)
      
      ! Retrieval/Error diagnostics

      real         :: Amatrix(nvar), sigma(nvar)

      ! Counters, etc.

      integer      :: niter, a, b, c, d, e, mc, i, ichan
      integer      :: tot_iter=0,noret_count=0,notrun_count=0 ! counters for oe
      real         :: dmr(nz), tpwerr_lyr(nz,3), tpw_posterr
      real         :: extra_ssd, extra_mr(nz), extra_eo(nz,6)
      real         :: offs(maxchans),fssd
      real(8)         :: psy(nobs,nobs), psy_i(nobs,nobs)    ! New: Sy can change with pixel
      LOGICAL :: sy_diag = .FALSE.   !!!  !Whether or not we force Sy to be a diagonal matrix



      allocate(oe_tbs(nch))
      ALLOCATE(oe_refl(nz))
      ALLOCATE(oe_obs(nch+nz))
      allocate(tbsim(nch))
      ALLOCATE(sim_obs(nobs))
      


        IF (sy_diag) THEN
                DO a = 1, nobs
                        DO b = 1, nobs
                                IF (a .NE. b) sy(a,b) = 0.
                        ENDDO
                ENDDO
        ENDIF



  mc = maxchans


        
      ! initialize oe_output fields, tb diff to miss_flt (not run!)

!      sa(:,:)            = 0.0
      iter_out(:)      =   0
      incr_out(:,:)      = miss_flt
      oe_output(:,:)   = miss_flt
      tbdif_out(:,:)   = miss_flt
      tbsim_out(:,:)   = miss_flt
      rmix_out(:,:)    = miss_flt
      clwc_out(:,:)    = miss_flt
      plwc_out(:,:)    = miss_flt
      tiwc_out(:,:)    = miss_flt
      saln_out(:)      = miss_flt
      rain_rate(:,:)   = miss_flt
      snow_rate(:,:)   = miss_flt
      sfcprcp_out(:)   = miss_flt
      retr_flag_out(:) = -99

      cltop_bin(:) = -99
      clbot_bin(:) = -99


      !--- Begin loop over pixels/scans ---!


        do oepix = pxstart, npix
          write(6,*) '---  pix = ',oepix


          sy_flag = -1

          sfc_type(oepix) = sfc_in(oepix)

          !!! Set pix to retrieve here if write_test == .True.
          !if (oepix .lt. 1588) cycle
          !!!

          !write(*,*) 'lat: ', lat_in(oepix), 'lon: ', lon_in(oepix)
          !write(*,*) 'ScanTime: ', dtime_in(oepix,:)

          IF (lat_in(oepix) .GT. -40. .AND. lat_in(oepix) .LT. 40.) then
                  write(*,*) 'not high lat'
                  cycle
          ENDIF

       
        
        !Check for ocean surface type, cycle if not:


            !All qc flagging done in preprocessor now

            IF (sfc_type(oepix) .NE. 0) then
                    write(*,*) 'not ocean'
                    retr_flag_out(oepix) = -9
                    cycle
            ENDIF
            IF (qflg_in(oepix) .NE. 0) then
                    write(*,*) 'qflag is ', qflg_in(oepix)
                    retr_flag_out(oepix) = -8
                    cycle
            ENDIF
 

            oe_refl(:) = refl_in(oepix,:)


            !Check for clear profiles
            IF (ALL(oe_refl .LE. -26.)) THEN
                    WRITE(*,*) 'Clear'
                    retr_flag_out(oepix) = 1
                    sfcprcp_out(oepix)  = 0.
                    sy_flag = 1 !Clear sy

            ELSE
                    !Check for cloudy but not precipitating profiles
                    IF (prcp_in(oepix) .LT. -10. .OR. snrt_sfc_in(oepix) .LT. 0.) THEN
                            WRITE(*,*) 'Bad CSAT precip rate: ', prcp_in(oepix), snrt_sfc_in(oepix)
                            retr_flag_out(oepix) = -99
                            sfcprcp_out(oepix)   = miss_flt
                            CYCLE
                    ENDIF
                    IF (ABS(prcp_in(oepix)) + snrt_sfc_in(oepix) .LT. 0.01) THEN
                            WRITE(*,*) 'CSAT shows no precip ', prcp_in(oepix), snrt_sfc_in(oepix)
                            retr_flag_out(oepix) = 2
                            sfcprcp_out(oepix) = 0.
                            sy_flag = 2 !Cloudy sy
                    ELSEIF (ABS(prcp_in(oepix)) + snrt_sfc_in(oepix) .GT. 0.01) THEN
                            WRITE(*,*) 'CSAT rain rate: ', prcp_in(oepix), ' CSAT snow rate: ', snrt_sfc_in(oepix)
                            retr_flag_out(oepix) = 3
                            sfcprcp_out(oepix) = -99.
                            sy_flag = 3 !Prcp sy
                    ELSE !something is wrong
                            WRITE(*,*) 'Something wrong with precip rate.', prcp_in(oepix), snrt_sfc_in(oepix)
                            retr_flag_out(oepix) = -99
                            sfcprcp_out(oepix) = miss_flt
                            CYCLE
                            STOP
                    ENDIF

            ENDIF

            IF (sy_flag .LT. 0) THEN
                    WRITE(*,*) 'Sy flag not assigned!'
                    STOP
            ENDIF

            
            !Set reflectivities below noise threshold to -30.
            DO i = 1, nz
                IF (oe_refl(i) .LE. -26.) &
                        oe_refl(i) = -30.
            ENDDO



            !Get cloud top and bottom bin:
            ct_bin = -99
            cb_bin = -99
            DO i = 1, nz
                IF (oe_refl(i) .GT. -26.) THEN
                        cb_bin = i
                        IF (ct_bin .GT. 0) CYCLE
                        ct_bin = i
                ENDIF
            ENDDO
                
    

            !Try to get rid of surface echo by setting cloud base to ECMWF cloud base
            IF (cldbse_bin_in(oepix) .GT. 0. .AND. cldbse_bin_in(oepix) .LT. cb_bin) &
                    cb_bin = cldbse_bin_in(oepix)
            IF (cb_bin .LE. ct_bin) ct_bin = cb_bin

            !Set reflectivities below cloud base to noise floor:

            oe_refl(cb_bin:nz) = -30.

            !Check again to see if this got rid of all clouds (sfc echo)
            IF (ALL(oe_refl .LE. -26.)) THEN
                    WRITE(*,*) 'Clear'
                    retr_flag_out(oepix) = 1
                    sfcprcp_out(oepix)   = 0.
                    sy_flag = 1
            ENDIF

            !Now, find cloud base bin again
            DO i = 1, nz - 1
                IF (oe_refl(i) .GT. -26.) THEN
                        cb_bin = i + 1
                ENDIF
            ENDDO

            !write(*,*) 'ct and cb: ', ct_bin, cb_bin
            !write(*,*) oe_refl

            !Read in appropriate sy matrix for this pixel
            CALL read_sy(sy_flag, radar_stat)


            !Initialize sy for this pixel
            psy(:,:) = -9999.9
            DO a = 1, nobs
                DO b = 1, nobs
                    psy(a,b) = sy(a,b)
                    IF (sy_flag .EQ. 3) THEN !If a prcp retrieval set uncertainties below cloud base really high
                        IF (a .GE. nch+cb_bin .OR. b .GE. nch+cb_bin) THEN
                                psy(a,b) = 0.
                                IF (a .EQ. b) psy(a,b) = 10000.
                        ENDIF
                    ENDIF
                ENDDO
            ENDDO


        cltop_bin(oepix) = ct_bin
        clbot_bin(oepix) = cb_bin

        ! define lat,sst,tbs for this pixel before calling OE
        oe_lat = lat_in(oepix) 
        oe_lon = lon_in(oepix)


        if (oe_lon .lt. 0) oe_lon = lon_in(oepix) + 360.0


        !Use ERA5 sst and wind speed as first guess
        !oe_sst     = sst_in(oepix)
        oe_sst     = rey_sst_in(oepix)
        oe_wspd    = wspd_in(oepix)

        IF (oe_sst .LE. 271.5) THEN !Sea ice possible
                WRITE(*,*) 'Sea ice possible, SST is: ', oe_sst
                retr_flag_out = -9
                CYCLE
        ENDIF

        


          ssdex = floor((oe_sst - 271.0) / 1.0) + 1 ! rounds down. 271.9 yields 1, etc.
          if (ssdex .lt. 1) ssdex = 1               ! first sst bin if <271K...
          if (ssdex .gt. nbins) ssdex = nbins       ! max sst bin

        !---Apply Tb offsets:

        oe_tbs(:) = tbobs_in(oepix,:) + toffsets(2:11)
        oe_tbs(1) = tbobs_in(oepix,1) + 2.  !Boost 10V
        oe_tbs(5) = tbobs_in(oepix,5) + 0.75
        oe_tbs(6) = tbobs_in(oepix,6) - 5.5 !Reduce 23H
        oe_tbs(9) = tbobs_in(oepix,9) - 1.2
        oe_tbs(10) = tbobs_in(oepix,10) - 1.2
        
          ! assign pixel values from analysis (in pp file)

          slp     = slp_in(oepix)
          wspd    = wspd_in(oepix)
          plev(:) = plev_in(oepix,:)
          tlev(:) = tlev_in(oepix,:)
          hlev(:) = hlev_in(:)
          frzl   = frzlvl_in(oepix)
          frzl_bin = frzlvl_bin_in(oepix)

          !Check for slp < lowest pressure value
          if (slp .le. plev(nz)) then
            plev(nz+1) = plev(nz) + 1.0
          else
            plev(nz+1) = slp
          endif

          do i = 1,nz ! get layer averages 
            dp(i) = plev(i+1) - plev(i)
            tavg  = (tlev(i) + tlev(i+1)) / 2.0
            tavg_out(oepix,i) = tavg
          enddo



          !Get mean mixing ratio for layer
          DO i = 1, nz
             !rmix(i) = EXP((LOG(rmix_in(oepix,i)) + LOG(rmix_in(oepix,i+1))) / 2.0)

             IF (rmix_in(oepix,i) .EQ. rmix_in(oepix,i+1)) THEN
                     rmix(i) = rmix_in(oepix,i)
                     CYCLE
             ENDIF

             rmix(i) = (rmix_in(oepix,i) - rmix_in(oepix,i+1)) / LOG(rmix_in(oepix,i) / rmix_in(oepix,i+1))
             
             IF (rmix(i) .LT. 0.) rmix(i) = 0.
          ENDDO



          if (minval(tavg_out(oepix,:)) .lt. 190) then
            print*,'Tprof too low!',oepix,oelin,tavg_out(oepix,:)
            stop
          endif



        salinity              = 35.0  !fix salinity as climatology


          incr_out(oepix,1) = tpw_in(oepix)
          incr_out(oepix,2) = wspd
          incr_out(oepix,3) = sst_in(oepix)
          incr_out(oepix,4) = sum(1.0/(1000.0*9.81) * rmix(:)*(dp(:))*100.0)



          ! Apriori error variances 
          ! Values of 0.00001 mean the parameters are fixed in the retrieval, ie not retrieved.
          ! Must not set any values to 0.0 since internal matrix 
          ! inversion will fail in the retrieval algorithm.

          !Latest version: noretrieval is not allowed
          noretrieval = 0.00001 


        !call inverse(DBLE(sy),sy_i,nobs)
        call inverse(psy,psy_i,nobs)



          ! Apriori state vector
          ! NOTE! -- Prior values for eof coefficients can matter a lot!
          !       -- Have to be non-zero. Big values will bias the answer, and 
          !       -- values too near zero cause more avg iterations needed!

          !Newest variable order: rwc 1-30, iwc 31-60, log10(cwp), mu, rho_snow, eof1, SST, wind

          xa(1:nz+nz) = 0. !re-initialize cloud water and ice content


          IF (sy_flag .EQ. 2) THEN !Cloudy retrieval
                !Only put liquid and ice water content in where we have reflectivities
                DO i = 1, frzl_bin - 1
                        IF (oe_refl(i) .GT. -26.) xa(nz+i) = 0.001 !ice water content (g/m^3) 
                ENDDO
                DO i = frzl_bin, nz
                        IF (oe_refl(i) .GT. -26.) xa(i) = 0.001 !precip water content (g/m^3)
                ENDDO
                IF (frzl_bin .GE. nz) xa(nz) = 0 !If freezing level is at or below surface, no precip water

                !If we are not using radar information, assume distribution of ice clouds (6-7km)
                !and no precip water for cloud-only scenes
                IF (radar_stat .EQ. 1) THEN
                        xa(nz+1:nz+nz)  = 0.
                        xa(nz+17:nz+18) = 0.001
                        xa(1:nz)        = 0.
                ENDIF

           ELSEIF (sy_flag .EQ. 3) THEN !prcp retrieval
                !Only put liquid and ice water content in where we have reflectivities
                DO i = 1, frzl_bin - 1
                        IF (oe_refl(i) .GT. -26.) xa(nz+i) = 0.001 !ice water content (g/m^3)
                ENDDO
                DO i = frzl_bin, nz
                        IF (oe_refl(i) .GT. -26.) xa(i) = 0.001 !precip water content (g/m^3)
                ENDDO
                !Assume distribution of precip below cloud base:
                IF (frzl_bin .GE. nz) THEN
                  xa(nz) = 0.  !If freezing level is at or below surface, no precip water.
                  xa(nz+nz) = 0.001 !Instead, put ice there
                ENDIF
                xa(cb_bin:nz) = xa(cb_bin-1)
                xa(cb_bin+nz:nz+nz) = xa(cb_bin+nz-1)
                !Sometimes the freezing level and cloud base are in the same bin,
                ! so put rain water in below cloud base
                IF (ALL(xa(cb_bin:nz) .EQ. 0.) .AND. ALL(xa(nz+cb_bin:nz+nz) .EQ. 0)) THEN
                  xa(cb_bin:nz) = 0.001
                ENDIF
                !Assume snow melts below freezing level:
                IF (frzl_bin .LT. nz) THEN
                        DO i = frzl_bin, nz
                                IF (xa(nz+i) .GT. 0.) THEN
                                        xa(nz+i) = 0.
                                        xa(i)    = 0.001
                                ENDIF
                        ENDDO
                ENDIF
           ENDIF  !retrieval type


           IF (radar_stat .EQ. 1) THEN !Radiometer-only retrieval: assume ice cloud layer
                xa(nz+1:nz+nz)  = 0.
                xa(nz+17:nz+18) = 0.001
                xa(1:nz)        = 0.
           ENDIF


          xa(nvar-5) = 1.         !log10(cloud water path [g/m^3])
          xa(nvar-4) = 1.5         !DSD shape parameter mu
          xa(nvar-3) = 0.1        !Snow particle density [g/cm^3]
          xa(nvar-2) = 0.02        !EOF1
          xa(nvar-1) = oe_sst      !SST [K]
          xa(nvar)   = oe_wspd     !Surface wind speed [m/s]



          !mrmp(:)    = rmix(:)         ! ERA5 water vapor profile is apriori
          mrmp(:)    = mprof(ssdex,:)
          peofs(:,:) = eofs(ssdex,:,:) ! initialize
	     
          !oe_obs is new total observation vector y
          oe_obs(1:nch)        = oe_tbs
          oe_obs(nch+1:nch+nz) = oe_refl



           call opt_est(nz, noretrieval, oe_obs, xa, sa, xr, xmax, xmin, &
                        psy, psy_i, flag, Amatrix, chisq, sigma, niter)

          
          !Calculate precip rates on final state:
          IF (flag) THEN !successful retrieval
                IF (sy_flag .EQ. 3) THEN !prcp retrieval
                        rn_rate = -9999.9
                        sn_rate = -9999.9
                        CALL precip_rate(xr, n0_lut, 5100.,rn_rate, sn_rate)

                        rain_rate(oepix,:) = rn_rate(:)
                        snow_rate(oepix,:) = sn_rate(:)
                        sfcprcp_out(oepix) = rn_rate(nz) + sn_rate(nz)

                ELSE !clear or cloudy retrieval
                        rain_rate(oepix,:) = 0.
                        snow_rate(oepix,:) = 0.
                        sfcprcp_out(oepix) = 0.
                ENDIF
          ELSE !unsuccessful retrieval
                rain_rate(oepix,:) = miss_flt
                snow_rate(oepix,:) = miss_flt
                sfcprcp_out(oepix) = miss_flt
          ENDIF



          ! posterior error on lwp is stddev(log10(lwp)). convert back to
          ! units of mm!
          if (flag) then ! lwp posterior error approximation
            sigma(nvar-5) = (10**(xr(nvar-5)+sigma(nvar-5))-10**(xr(nvar-5)-sigma(nvar-5)))/2.0 
            sigma(nvar-5) = sigma(nvar-5)*1000.0 ! convert into g/m^2


            write(*,*) 'Retrieved precip rate: ', sfcprcp_out(oepix)

            ! tpw posterior error:
            tpwerr_lyr(:,:) = 0.0
            do i = 1, nz
              do c = 1, npc
                dmr(i) = sigma(c)*peofs(i,c)
                tpwerr_lyr(i,c) = 1.0/(10.0*9.81) * dmr(i)*dp(i)
              enddo
            enddo
            if (npc==3) tpw_posterr = sqrt( (sum(tpwerr_lyr(:,1)))**2 + &
                                      (sum(tpwerr_lyr(:,2)))**2 + &
                                      (sum(tpwerr_lyr(:,3)))**2)
            if (npc==2) tpw_posterr = sqrt( (sum(tpwerr_lyr(:,1)))**2 + &
                                      (sum(tpwerr_lyr(:,2)))**2 )
          endif


          if (chisq .ge. 0.0) chisq = chisq / nobs ! normalize ChiSq by nchannels

          !Final output of retrieved state
          if (flag .and. chisq .lt. chisq_out_thresh) then

            post_out(oepix,1:nvar)  = sigma(:)    ! posterior errors
            post_out(oepix,nvar+1) = tpw_posterr  ! separate!
            run_count = run_count+1
            tot_iter = tot_iter + niter                 ! sum iterations for successful retrievals
            iter_out(oepix) = niter
            last_oeout(:) = xr(:)                       ! save for fg in next pixel
            loeop = oepix                               ! save last converged pixel number
            oe_output(oepix,1:nz) = xr(1:nz)        !Layers of precip liquid water content
            oe_output(oepix,nz+1:nz+nz) = xr(nz+1:nz+nz) !Layers of ice 
            oe_output(oepix,nvar-5) = xr(nvar-5)    !Cloud water path
            oe_output(oepix,nvar-4) = xr(nvar-4)    ! DSD mu
            oe_output(oepix,nvar-3) = xr(nvar-3)    !snow density
            oe_output(oepix,nvar-2) = xr(nvar-2)    !wv EOF1
            oe_output(oepix,nvar-1) = xr(nvar-1)    !SST
            oe_output(oepix,nvar)   = xr(nvar)      !Wspd
            !Additional output parameters
            oe_output(oepix,nvar+1) = chisq         !Chi squared (y-F(x))^T S_y^{-1} (y-F(x))
            oe_output(oepix,nvar+2) = final_tpw     !TPW from integrating wv profile


          endif


          ! Don't write out values for high chi squared values:

          if (flag .and. chisq .ge. chisq_out_thresh) then
            !if (outdebug) write(*,*),'=== ChiSqr too big!'
            oe_output(oepix,:) = miss_flt
            post_out(oepix,:) = sigma(:)
            iter_out(oepix) = niter
            noret_count = noret_count+1
            retr_flag_out(oepix) = -retr_flag_out(oepix)
            cycle
          endif
          if (flag .eq. 0 .or. final_tpw .gt. 75) then 
            !if (outdebug) write(6,*) '=Max iterations or tpw > 75 ',oepix,oelin
            oe_output(oepix,:) = miss_flt
            post_out(oepix,:)  = miss_flt  ! posterior errors
            noret_count = noret_count+1
            iter_out(oepix) = niter
            retr_flag_out(oepix) = -retr_flag_out(oepix)
            cycle
          endif


          ! for successful retrievals only, save Tbs

          tbdif_out(oepix,:) = oe_tbs(:) - tbsim(:) 
          tbsim_out(oepix,:) = tbsim(:)

          !write(6,'("pix   = ",i5)') oepix
          !write(6,'("scan  = ",i5)') oelin
          !write(6,'("slp   = ",f9.2)') slp
          !write(6,'("tpw   = ",f9.2)') incr_out(oepix,oelin,1)
          !write(6,'("wsp   = ",f9.2)') incr_out(oepix,oelin,2)
          !write(6,'("wdir  = ",f9.2)') wdir
          !write(6,'("tsfc  = ",f9.2)') incr_out(oepix,oelin,3)
          !write(6,'("tpw   = ",f9.2)') incr_out(oepix,oelin,4)
          !write(6,'("saln  = ",f9.2)') salinity
          !write(6,'("rmix  = ",20(f7.2))') rmix(:)
          !write(6,'("dp    = ",20(f7.2))') dp(:)
          !write(6,'("plev  = ",20(f7.2))') plev(:)
          !write(6,'("tlev  = ",20(f7.2))') tlev(:)
          !write(6,'("hlev  = ",20(f7.2))') hlev(:)
          !write(6,'("tbobs = ",20(f7.2))') oe_tbs(:)
          !write(6,'("tbsim = ",20(f7.2))') tbsim(:) 
          !write(6,'("chisq = ",f9.2)') chisq
        enddo  !pixels

      avg_iter = real(tot_iter) / real(run_count) ! iterations per retrieval

      return
      end subroutine retrieval_nr



!---------------------------------------------------------------------------

!-------------------------------------------------------------------
!
!       NO SCATTER
!
!-------------------------------------------------------------------


      subroutine noscatter(X,sim_out)

      !-------------------------------------------------------------
      ! ROUTINES CALLED:
      ! compute_tb:  calls monortm-lut version of Eddington RT model
      !-------------------------------------------------------------

      real    :: tavg(nz)
      real    :: tsfc, wind, lwp
      REAL    :: mu
      real    :: icthick, icbase, iwp, tiwc(nz), clwc(nz)
      REAL    :: plwc(nz)
      REAL    :: n0_rain, mu_rain, n0_snow, rho_snow
      REAL    :: cwp
      real    :: X(nvar),tbout(nch)
      REAL    :: sim_out(nobs)
      real    :: tbs(nch), tpw_lyr(nz)
      real    :: z_e(nz), z_a(nz), atten(nact_freqs)
      real    :: rmix(nz)
      integer :: cldbindex, cldtindex, icbindex,ictindex
      logical :: ice
      real    :: eof_coef(3) ! eof coefficients
      real    :: rh(nz), reh, mixr
      integer :: i,j,k,n,m,z


      !new variable order: liquid precip water, ice water, cloud water path, mu, rho_sn,  EOF1-3, sst, wind
      plwc(:)         = X(1:nz)
      tiwc(:)         = X(nz+1:nz*2)
      cwp             = 10**X(nvar-5)
      mu_rain         = X(nvar-4)
      rho_snow        = X(nvar-3)
      eof_coef(1:npc) = X(nvar-2)
      tsfc            = X(nvar-1)
      wind            = X(nvar)

      clwc(:)         = 0.  !Reinitialize cloud water

        IF (radar_stat .EQ. 1) THEN

             !If not using radar information, assume distribution of clouds
             
             !Cloud water is 1-2km
             clwc(27:28) = cwp / 2. / 500.

             !Unless freezing level is at the surface
             IF (frzl_bin .GE. nz) clwc(nz) = cwp / 500.


             !Cloud ice is 6-7 km. Done in xa.
             

         ELSE      
             !Distribute cloud liquid water between freezing level and cloud bottom
   
                IF (frzl_bin .LE. cb_bin) THEN
                        clwc(frzl_bin-1:cb_bin) = cwp / (cb_bin - frzl_bin + 1) / 500.
                ELSE
                        clwc(:) = 0.
                ENDIF

                IF (frzl_bin .GE. nz) clwc(nz) = cwp / 500.

        ENDIF

        IF (sy_flag .eq. 3) THEN
                !Assume distribution of precip below cloud base:
                plwc(cb_bin:nz) = plwc(cb_bin-1) !copy plwc from cloud base to surface
                tiwc(cb_bin:nz) = tiwc(cb_bin-1) !copy swc from cloud base to surface
                IF (ALL(plwc(cb_bin:nz) .EQ. 0.) .AND. ALL(tiwc(cb_bin:nz) .EQ. 0.)) THEN
                        plwc(cb_bin:nz) = tiwc(cb_bin-1)
                ENDIF
                IF (tiwc(frzl_bin) .GT. 0.) THEN
                        plwc(frzl_bin:nz) = tiwc(frzl_bin) + plwc(frzl_bin)
                        tiwc(frzl_bin:nz) = 0.
                ENDIF
        ENDIF



      tavg(:) = tavg_out(oepix,:) ! calculated above!


      ! start new water vapor profile code here:
      ! the methodology: 
      ! 1. read in mean WV Mixing Ratio profile for given sst bin.
      ! 2. upon iterating, allow first X EOFs to vary in
      !    magnitude to modify MR profile
      ! 3. output final profile in MR space, plus tpw value
      
      do i = 1, nz
        rmix(i) = mrmp(i) + &
                  eof_coef(1) * peofs(i,1)! + & ! changed to peofs (BCMZ)
                  !eof_coef(2) * peofs(i,2) + & 
                  !eof_coef(3) * peofs(i,3) 
        if (rmix(i) .lt. 0) rmix(i) = 0.0 ! can't have negative mass!
      enddo


      rmix_out(oepix,:) = rmix(:)
      plwc_out(oepix,:) = plwc(:)
      clwc_out(oepix,:) = clwc(:)
      tiwc_out(oepix,:) = tiwc(:)


      final_tpw=sum(1.0/(1000.0*9.81) * rmix(:)*(dp(:))*100.0)
      
      !write(*,*),'tpw: ',final_tpw
  

!----- pass lyr temperature, wv mixing ratio,
!       lwp/iwp per layer, sst, wind speed

 !     write(*,*) '-----------------------------------'
 !     write(6,'("plev   = ",20(f9.2))') plev
 !     write(6,'("tlev   = ",20(f9.2))') tlev
 !     write(6,'("hlev   = ",20(f9.2))') hlev
 !     write(6,'("rmix = ",20(f9.2))') rmix
 !     write(6,'("clwc   = ",20(f9.4))') clwc
 !     write(6,'("tiwc   = ",20(f9.4))') tiwc
 !     write(6,'("tsfc   = ",20(f9.2))') tsfc
 !     write(6,'("saln   = ",20(f0.2))') salinity
 !     write(6,'("wind   = ",20(f9.2))') wind
 !     write(6,'("eia    = ",20(f9.2))') eia_in(oepix,oelin,:)
 !     write(*,*) '------------------------------------'




n0_rain = n0_lut
n0_snow = 5100 ![mm^-1 m^-3]

IF (radar_stat .EQ. 1) rho_snow = 0.1 !set snow particle density to fixed value if no radar data



       call forward_model(hlev/1000., plev, tlev, rmix, clwc, plwc, tiwc, &
                          tsfc, salinity, wind, n0_rain, mu_rain, n0_snow, rho_snow, &
                          tbs, z_e, z_a, atten)

      tbout(:) = tbs(:)
      tbsim(:) = tbs(:)
     
      !Tbout is now sim_out and tbsim is now sim_obs
      sim_out(1:nch)      = tbs(:)
      sim_out(nch+1:nobs) = z_a(:)
      sim_obs(1:nch)      = tbs(:)
      sim_obs(nch+1:nobs) = z_a(:)

      !write(6,'("tbobs = ",15(f8.2))') tbobs_in(oepix,:)
      !write(6,'("tbsim = ",15(f8.2))') tbsim(:)
      !write(6,*)
      !write(6, '("Z_e = ", 30(f10.2))') z_e(:)
      !write(6, '("Z_in = ",15(f10.2))') oe_refl(:) 
      !write(6, '("Z_a = ", 30(f10.2))') z_a(:)
      !write(6, '("PIA = ", f8.2)')     atten
      !write(6,*) sim_out
      
      return

      end subroutine noscatter 




!-------------------------------------------------------------------------------



      subroutine opt_est(nz, noretrieval, y, xa2, sa2,x2, xmax, &
                         xmin, psy, psy_i, flag, Amatrix2, chisq, &
			             sigma2, niter)

      implicit none

      integer :: nz   ! number of layers (nz) in atmosphere.  
      logical :: flag ! flag indicates whether convergence was reached 

      ! Variables used in optimal estimation framework .

      integer :: n,niter ! defined above too!!
      integer :: a,b,c,d,i,count,index
      integer :: checkpol,end_flag
      real    :: noretrieval
      real    :: x(nretvar),x2(nvar)
      real    :: xnext(nretvar),xprime2(nvar)
      real    :: xa(nretvar),xa2(nvar)
      real    :: xmax(nretvar),xmin(nretvar)
      real(8)    :: sa(nretvar,nretvar)
      real    :: sa2(nvar,nvar)             !input Sa matrix
      real(8)    :: sa_i(nretvar,nretvar)
      !real    :: psy(nch,nch)		     ! New: pixel Sy
      !real    :: psy_i(nch,nch)		     ! Sy inverse
      REAL(8)    :: psy(nobs,nobs)
      REAL(8)    :: psy_i(nobs,nobs)
      !real    :: K(nch,nretvar)
      REAL(8)    :: K(nobs,nretvar)
      !real    :: K_t(nretvar,nch)
      REAL(8)    :: K_t(nretvar,nobs)
      real(8)    :: sx(nretvar,nretvar)         ! posteriori error matrix
      real(8)    :: sx_i(nretvar,nretvar)
      real    :: sigma2(nvar)
      real    :: AP(nretvar,nretvar)
      real    :: Amatrix2(nvar)
      !real    :: Fprime(nch),Fdbprime(nch)
      !real    :: y(nch),dF(nch),dx(nretvar)
      !real    :: F(nch),Fout(nch)
      !real    :: Foutprime(nch)              ! changed!
      REAL    :: Fprime(nobs), Fdbprime(nobs)
      REAL    :: y(nobs), dF(nobs), dx(nretvar)
      REAL    :: F(nobs), Fout(nobs)
      REAL    :: Foutprime(nobs)
      real    :: sum1(nretvar,1)             ! xa-x
      real    :: sum2(nobs,1) !sum2(nch,1)  ! y-F
      real    :: sum3(nretvar,1)             ! prod4+prod3
      real(8)    :: prod1(nretvar,nobs) !prod1(nretvar,nch)  ! K_t*sy_i
      real(8)    :: prod2(nretvar,nretvar) ! K_t*sy_i*K
      real    :: prod3(nretvar)              ! sa_i*(sum1)
      real    :: prod4(nretvar)              ! prod1*sum2
      real    :: prod5(nretvar)              ! sx*sum3
      real    :: xdiff(nretvar,1)            ! xnext-x
      real    :: xdiff_t(1,nretvar)          ! transposed
      real    :: prod6(1,nretvar)            ! xdiff_t*sx_i
      real    :: prod7(1,1)                  ! prod6*xdiff ! CHANGED!!
      real    :: sum4(nobs,1) !sum4(nch,1) ! F-y
      real    :: sum4_t(1,nobs) !sum4_t(1,nch) ! (F-y) transpose
      real    :: sum5(nretvar,1)             ! x-xa
      real    :: sum5_t(1,nretvar)           ! (x-xa) transpose
      real    :: prod8(1,nobs) !prod8(1,nch) ! sum4_t*sy_i
      real    :: prod9(1,1)                  ! prod8*sum2
      real    :: prod10(1,nretvar)           ! sum5_t*sa_i
      real    :: prod11(1,1)                 ! prod10*sum5 ! CHANGED!!!
      real    :: chisqtot(1,1)               ! chi squared (prod9 + prod11) ! CHANGED!!
      real    :: chisq                       ! cost function (apriori + measurement fit)
      real    :: sum6(nretvar,nretvar)       ! (IM-AP)
      real    :: prod12(nretvar,nobs) !prod12(nretvar,nch) ! Sx*K_t
      real    :: prod13(nretvar,nobs) !prod13(nretvar,nch) ! Sx*K_t*Sy_i
      real    :: prod14(nretvar)             ! prod13*sum2 - contribution from obs
      real    :: prod15(nretvar)             ! sum6*sum1 - contribution from apriori
      real    :: prod16(nretvar)             ! prod14+prod16 = xnext-x

      !real    :: Fold(nch)                   ! old fwd modeled Tbs
      !real    :: diffF(nch,1),diffF_t(1,nch) ! F-F_old, transpose
      !real    :: Sdy(nch,nch),Sdy_i(nch,nch) ! (Rodgers Eq5.27) 
      !real    :: Ksa(nch,nretvar)            ! K*Sa
      !real    :: Ksakt(nch,nch)              ! K*Sa*K_t
      !real    :: Sdysum(nch,nch)             ! (K*Sa*K_t+Sy)
      !real    :: Sdysum_i(nch,nch)           ! (K*Sa*K_t+Sy) *_i
      !real    :: almost(nch,nch)             ! Sy*(Ksakt+Sy)_i (most of 5.27)
      !real    :: almost3(1,nch)              ! diffF_t*Sdy_t
      real    :: Fold(nobs)
      real    :: diffF(nobs,1),diffF_t(1,nobs) ! F-F_old, transpose
      real(8)    :: Sdy(nobs,nobs),Sdy_i(nobs,nobs) ! (Rodgers Eq5.27) 
      real(8)    :: Ksa(nobs,nretvar)            ! K*Sa
      real(8)    :: Ksakt(nobs,nobs)              ! K*Sa*K_t
      real(8)    :: Sdysum(nobs,nobs)             ! (K*Sa*K_t+Sy)
      real(8)    :: Sdysum_i(nobs,nobs)           ! (K*Sa*K_t+Sy) *_i
      real(8)    :: almost(nobs,nobs)             ! Sy*(Ksakt+Sy)_i (most of 5.27)
      real(8)    :: almost3(1,nobs)              ! diffF_t*Sdy_t
      real*8  :: disq(1,1)                   ! di squared (Eq5.33), used for conv test

      real    :: psdiff                      ! pix/scan distance from last convergence reached 
      integer :: good(nretvar)
      
      ! Added by Rick S. March 2018
      real    :: mincost, chisqmincost
      real    :: Fmincost(nobs)
      real    :: x2mincost(nvar)
      integer :: iter_mincost

        logical :: write_test = .FALSE.  !Whether or not to write out iteration information to test_OE_pix.bin
        character(len=100) :: write_test_file

      mincost = 999.9

      

      !-----Declare Apriori errors 

        

      !---Reorder variables, get rid of noretrievals
      sa(:,:) = 0.0
      b = 1
      do a=1,nvar ! less complex if not using covariances...
        if (sa2(a,a) .ne. noretrieval) then 
          good(b) = a
          b = b + 1
        endif
      enddo
      do a=1,nretvar
       do c=1,nretvar
         sa(a,c) = sa2(good(a),good(c))    !New sa is (nretrvars x nretrvars)
       enddo
      enddo

        write_test_file = 'test_oepix.bin'

        
      IF (write_test) THEN
              CALL SYSTEM('rm '//TRIM(write_test_file))
              OPEN(unit=500, file=write_test_file, access='stream', status='unknown')
              WRITE(500) oepix
              WRITE(500) prcp_in(oepix)
              WRITE(500) snrt_sfc_in(oepix)
              WRITE(500) nvar
              WRITE(500) nretvar
              WRITE(500) nobs
              WRITE(500) sa
              WRITE(500) psy
              CLOSE(500)
      ENDIF


        b=1
        do a = 1, nvar
          x2(a) = xa2(a)
          if (sa2(a,a) .ne. noretrieval) then
            x(b)  = xa2(a)  ! simply use a priori values
            xa(b) = xa2(a)
            b=b+1
          endif
        enddo


      IF (write_test) THEN
              OPEN(unit=500, file=write_test_file, access='stream', status='old', position='append')
              WRITE(500) xa
              WRITE(500) mrmp
              WRITE(500) peofs
              WRITE(500) y
              WRITE(500) frzl
              WRITE(500) frzl_bin
              WRITE(500) cltop_bin(oepix)
              WRITE(500) clbot_bin(oepix)
              CLOSE(500)
      ENDIF

        !--- Newtonian iteration

      end_flag=0
      

!      write(*,*),'Y: ',y

      do n = 1, nn


      !write(*,*) 'Iteration: ', n

      IF (write_test) THEN
              OPEN(unit=500, file=write_test_file, access='stream', status='old', position='append')
              WRITE(500) n
              CLOSE(500)
      ENDIF


      IF (write_test) THEN
              OPEN(unit=500, file=write_test_file, access='stream', status='old', position='append')
              WRITE(500) x2
              CLOSE(500)
      ENDIF


        !---Call forward model on first guess
       ! write(*,*),'n iter',n
       ! write(*,*),'X: ',x2(:)
       ! write(*,*) 'Obs: ', y
        !write(*,*) 
        call noscatter(x2,Fout)
       ! write(*,*),'tpw: ',final_tpw
       ! write(*,*),'Fout: ',Fout

      IF (write_test) THEN
              OPEN(unit=500, file=write_test_file, access='stream', status='old', position='append')
              WRITE(500) Fout
              CLOSE(500)
      ENDIF

        F(:) = Fout(:)



        !--- calculate Jacobian K

        !Try absolute perturbations instead of relative:
        
        !initialize xprime
        DO a = 1, nvar
           xprime2(a) = x2(a)
        ENDDO

        a = 1
        DO d = 1, nz
                IF (xa2(d) .EQ. 0.) THEN !no reflectivity, so nothing in this layer 
                        K(:,a) = 0.
                        a = a + 1
                        CYCLE
                ENDIF
                xprime2(d) = x2(d) + 0.02*x2(d)  !2% perturbation
                if (xprime2(d) .le. 0.) xprime2(d) = 1.0E-06  !Don't let rain water go negative
                dx(a) = xprime2(d) - x2(d)
                CALL noscatter(xprime2, Foutprime)
                Fprime(:) = Foutprime(:)
                dF(:) = Fprime(:) - F(:)
                IF (dx(a) .EQ. 0.) THEN
                        K(:,a) = 0.
                ELSE
                        K(:,a) = dF(:) / dx(a)
                ENDIF
                a = a + 1
                xprime2(d) = x2(d)
        ENDDO
        DO d = nz+1, nz+nz
                IF (xa2(d) .EQ. 0.) THEN
                        K(:,a) = 0.
                        a = a + 1
                        CYCLE
                ENDIF
                xprime2(d) = x2(d) + 0.02*x2(d)
                if (xprime2(d) .le. 0.) xprime2(d) = 1.0E-06 !Don't let ice go negative
                dx(a) = xprime2(d) - x2(d)
                CALL noscatter(xprime2, Foutprime)
                Fprime(:) = Foutprime(:)
                dF(:) = Fprime(:) - F(:)
                IF (dx(a) .EQ. 0.) THEN
                        K(:,a) = 0.
                ELSE
                        K(:,a) = dF(:) / dx(a)
                ENDIF
                a = a + 1
                xprime2(d) = x2(d)
        ENDDO

        !CWP
        xprime2(nvar-5) = x2(nvar-5) + 0.02*x2(nvar-5)
        dx(a) = xprime2(nvar-5) - x2(nvar-5)
        CALL noscatter(xprime2, Foutprime)
        Fprime(:) = Foutprime(:)
        dF(:) = Fprime(:) - F(:)
        IF (dx(a) .EQ. 0.) THEN
                K(:,a) = 0.
        ELSE
                K(:,a) = dF(:) / dx(a)
        ENDIF
        a = a + 1
        xprime2(nvar-5) = x2(nvar-5)

        !DSD mu
        xprime2(nvar-4) = x2(nvar-4) + 0.05
        if (xprime2(nvar-4) .lt. 0.)  xprime2(nvar-4) = 0.    !Set limits for K calculation
        if (xprime2(nvar-4) .gt. 2.5) xprime2(nvar-4) = 2.5
        dx(a) = xprime2(nvar-4) - x2(nvar-4)
        CALL noscatter(xprime2, Foutprime)
        Fprime(:) = Foutprime(:)
        dF(:) = Fprime(:) - F(:)
        IF (dx(a) .EQ. 0.) THEN
                K(:,a) = 0.
        ELSE
                K(:,a) = dF(:) / dx(a)
        ENDIF
        a = a + 1
        xprime2(nvar-4) = x2(nvar-4)

        !Snow density
        xprime2(nvar-3) = x2(nvar-3) + 0.02*x2(nvar-3)
        if (xprime2(nvar-3) .lt. 0.05) xprime2(nvar-3) = 0.05 !Set limits for K calculation
        if (xprime2(nvar-3) .gt. 0.4)  xprime2(nvar-3) = 0.4
        dx(a) = xprime2(nvar-3) - x2(nvar-3)
        CALL noscatter(xprime2, Foutprime)
        Fprime(:) = Foutprime(:)
        dF(:) = Fprime(:) - F(:)
        IF (dx(a) .EQ. 0.) THEN
                K(:,a) = 0.
        ELSE
                K(:,a) = dF(:) / dx(a)
        ENDIF
        a = a + 1
        xprime2(nvar-3) = x2(nvar-3)

        !EOF1
        xprime2(nvar-2) = x2(nvar-2) + 0.02
        dx(a) = xprime2(nvar-2) - x2(nvar-2)
        CALL noscatter(xprime2, Foutprime)
        Fprime(:) = Foutprime(:)
        dF(:) = Fprime(:) - F(:)
        IF (dx(a) .EQ. 0.) THEN
                K(:,a) = 0.
        ELSE
                K(:,a) = dF(:) / dx(a)
        ENDIF
        a = a + 1
        xprime2(nvar-2) = x2(nvar-2)

        !SST
        xprime2(nvar-1) = x2(nvar-1) + 0.25
        dx(a) = xprime2(nvar-1) - x2(nvar-1)
        CALL noscatter(xprime2, Foutprime)
        Fprime(:) = Foutprime(:)
        dF(:) = Fprime(:) - F(:)
        IF (dx(a) .EQ. 0.) THEN
                K(:,a) = 0.
        ELSE
                K(:,a) = dF(:) / dx(a)
        ENDIF
        a = a + 1
        xprime2(nvar-1) = x2(nvar-1)

        !wspd
        xprime2(nvar) = x2(nvar) + 0.1
        dx(a) = xprime2(nvar) - x2(nvar)
        CALL noscatter(xprime2, Foutprime)
        Fprime(:) = Foutprime(:)
        dF(:) = Fprime(:) - F(:)
        IF (dx(a) .EQ. 0.) THEN
                K(:,a) = 0.
        ELSE
                K(:,a) = dF(:) / dx(a)
        ENDIF
        xprime2(nvar) = x2(nvar)


      IF (write_test) THEN
              OPEN(unit=500, file=write_test_file, access='stream', status='old', position='append')
              WRITE(500) K
              CLOSE(500)
      ENDIF

        !--- evaluate newtonian matrices

        K_t(:,:) = transpose(K(:,:))
        call inverse(sa,sa_i,nretvar)
        prod1 = matmul(K_t,psy_i)
        prod2 = matmul(prod1,K)
        sx_i = sa_i + prod2
        call inverse(sx_i,sx,nretvar)
        

      IF (write_test) THEN
              OPEN(unit=500, file=write_test_file, access='stream', status='old', position='append')
              WRITE(500) sx_i
              WRITE(500) sx
              CLOSE(500)
      ENDIF



        !--- perform newtonian step

        sum1(:,1) = xa(:) - x(:)
        sum2(:,1) = y(:) - F(:)
        prod3     = matmul(sa_i,sum1(:,1))
        prod4     = matmul(prod1,sum2(:,1))
        sum3(:,1) = prod4(:) + prod3(:)
        prod5     = matmul(sx(:,:),sum3(:,1))
        xnext     = x + prod5

      

            !--- Check limits on xnext


        !The below indices have to be changed if nretvar is changed

        !Rain
        DO a = 1, nz
                IF (xa(a) .GT. 0. .AND. xnext(a) .LE. 0.) xnext(a) = 1.0E-06
        ENDDO
        !ice
        DO a = nz+1, nz+nz
                IF (xa(a) .GT. 0. .AND. xnext(a) .LE. 0.) xnext(a) = 1.0E-06
        ENDDO

        IF (sy_flag .EQ. 3) THEN
                xnext(cb_bin:nz)       = xnext(cb_bin-1)    !Copy rwc from layer above cloud base to surface
                xnext(cb_bin+nz:nz+nz) = xnext(cb_bin+nz-1) !Copy swc to sfc
                xnext(frzl_bin+nz:nz+nz) = 0.               !Set ice below freezing level to 0.
                IF (frzl_bin .GE. nz) xnext(nz+nz) = xnext(nz+nz-1) !If the frzl is at the ground, put snow there
                IF (ALL(xnext(cb_bin:nz) .EQ. 0.) .AND. ALL(xnext(cb_bin+nz:nz+nz) .EQ. 0.)) &
                        xnext(cb_bin:nz) = xnext(cb_bin+nz-1) !Don't let both rwc and swc to be 0 in precip scene
                IF (frzl_bin .LT. nz) THEN !Check for snow below freezing level and melt it
                        IF (xnext(frzl_bin+nz) .GT. 0.) THEN
                                xnext(frzl_bin:nz) = xnext(frzl_bin+nz) + xnext(frzl_bin) !transfer snow water to rain water content
                                xnext(frzl_bin+nz:nz+nz) = 0.                             !set swc to 0
                        ENDIF
                ENDIF
        ENDIF

        !CWP
        IF (xnext(nretvar-5) .LT. x2(nvar-5) - 1.) xnext(nretvar-5) = x2(nvar-5) - 1. 
        IF (xnext(nretvar-5) .GT. x2(nvar-5) + 1.) xnext(nretvar-5) = x2(nvar-5) + 1.

        !DSD mu
        IF (xnext(nretvar-4) .LT. 0.) xnext(nretvar-4) = 0.
        IF (xnext(nretvar-4) .GT. 2.5) xnext(nretvar-4) = 2.5
        !Snow particle density
        IF (xnext(nretvar-3) .LT. 0.05) xnext(nretvar-3) = 0.05
        IF (xnext(nretvar-3) .GT. 0.4) xnext(nretvar-3) = 0.4
        IF (xnext(nretvar-3) .GT. x2(nvar-3) + 0.1) xnext(nretvar-3) = x2(nvar-3) + 0.1
        IF (xnext(nretvar-3) .LT. x2(nvar-3) - 0.1) xnext(nretvar-3) = x2(nvar-3) - 0.1
        !sst
        IF (xnext(nretvar-1) .LT. 271.5) xnext(nretvar-1) = 271.5
        IF (xnext(nretvar-1) .LT. xa(nretvar-1)-2.) xnext(nretvar-1) = xa(nretvar-1) - 2. !!!
        IF (xnext(nretvar-1) .GT. xa(nretvar-1)+2.) xnext(nretvar-1) = xa(nretvar-1) + 2. !!!
        !wind
        IF (xnext(nretvar) .LT. 0.5) xnext(nretvar) = 0.5


      IF (write_test) THEN
              OPEN(unit=500, file=write_test_file, access='stream', status='old', position='append')
              WRITE(500) xnext
              CLOSE(500)
      ENDIF 

        !--- error diagnostics

        AP            = matmul(sx,prod2)
        sum4(:,1)     = F(:) - y(:)
        sum4_t        = transpose(sum4)
        sum5(:,1)     = x(:) - xa(:)
        sum5_t        = transpose(sum5)
        prod8(:,:)    = matmul(sum4_t(:,:),psy_i(:,:))
        prod9(:,:)    = matmul(prod8(:,:),sum4(:,:))
        prod10(:,:)   = matmul(sum5_t(:,:),sa_i(:,:))
        prod11(:,:)   = matmul(prod10(:,:),sum5(:,:))
        chisqtot(1,1) = prod9(1,1) + prod11(1,1)
	    chisq = prod9(1,1)

      IF (write_test) THEN
              OPEN(unit=500, file=write_test_file, access='stream', status='old', position='append')
              WRITE(500) chisqtot
              IF (radar_stat .EQ. 0) THEN     !Combined
                      WRITE(500) chisq / nobs
              ELSEIF (radar_stat .EQ. 1) THEN !Radiometer only
                      WRITE(500) chisq / nch
              ELSEIF (radar_stat .EQ. 2) THEN !Radar only
                      WRITE(500) chisq / nz
              ELSE
                      WRITE(500) -99.
              ENDIF
              CLOSE(500)
      ENDIF


	    ! New: is this the lowest chisq yet?
	    !if (chisqtot(1,1) .lt. mincost) then
            if (chisq .lt. mincost) then
	      !mincost = chisqtot(1,1)
              mincost = chisq
	      chisqmincost = chisq
	      Fmincost = F(:)
	      x2mincost = x2(:)
              iter_mincost = n
          ! Calculate Amatrix, posterior error standard deviation
	      b = 1
	      do a = 1,nvar 
            if (sa2(a,a).NE.noretrieval) then	  
              Amatrix2(a)=AP(b,b) 
	          sigma2(a)  =sqrt(sx(b,b)) ! -- posterior error standard deviations
	          b = b + 1
	        else
              Amatrix2(a) = 0.0
	          sigma2(a)   = miss_flt 
	        endif
          enddo
	    endif  	  
        
            


        !--- Evaluate closeness of Xnext and X

        xdiff(:,1)   = xnext(:) - x(:)
        xdiff_t(:,:) = transpose(xdiff(:,:))
        prod6(1,:)   = matmul(xdiff_t(1,:),sx_i(:,:))
        prod7(:,:)   = matmul(prod6(:,:),xdiff(:,:))
	
!	    write(*,*),'xnext: ',xnext

        b = 1
        do a = 1, nvar
          if (sa2(a,a) .ne. noretrieval) then
            x(b)  = xnext(b)
            x2(a) = xnext(b)
            b = b + 1
          endif
        enddo
        
        
        !--- As final retrieved state is found, TBs need to be recomputed
        call noscatter(x2,Fout)

        IF (write_test) THEN
                OPEN(unit=500, file=write_test_file, access='stream', status='old', position='append')
                WRITE(500) Fout
                rn_rate(:) = -9999.9
                sn_rate(:) = -9999.9
                call precip_rate(x2, n0_lut, 5100., rn_rate, sn_rate)
                WRITE(500) rn_rate
                WRITE(500) sn_rate
                CLOSE(500)
        ENDIF

        !write(6,'(i2,") tbout: ",13(f8.2))') n,Fout(:)
        F(:) = Fout(:)
        if (n .gt. 1) then
          diffF(:,1) = F(:) - Fold(:)
          diffF_t = transpose(diffF)
          Ksa = matmul(K,sa)
          Ksakt = matmul(Ksa,K_t)
          Sdysum = Ksakt + psy
          call inverse(Sdysum,Sdysum_i,nobs)
          almost = matmul(psy,Sdysum_i)
          Sdy = matmul(almost,psy)
          call inverse(Sdy,Sdy_i,nobs)
          almost3(1,:) = matmul(diffF_t(1,:),Sdy_i(:,:)) ! try this????
          disq = matmul(almost3, diffF)




          IF (write_test) THEN
              OPEN(unit=500, file=write_test_file, access='stream', status='old', position='append')
              WRITE(500) disq
              WRITE(500) iter_mincost
              CLOSE(500)
          ENDIF

!          write(*,*) prod7

        !TEST 1
          ! Convergence criteria
          if (prod7(1,1) .lt. nretvar/100.) then ! Eq 5.29 in Rodgers 2000 
                write(*,*) 'Convergence reached.'
                if (write_test) then
                        write(*,*) 'Convergence reached.'
                        stop 888
                endif
                  end_flag=1
          endif

        !TEST 2
        ! New: check if niter > 5 and chisq not improving (if so then use lowest available)
        !if (n.gt.5 .and. chisqtot(1,1).gt.mincost .and. chisqmincost .lt. chisq_out_thresh) then
        IF (sy_flag .EQ. 1 .OR. sy_flag .EQ. 2) THEN
          if (n .gt. 5 .and. chisq .gt. mincost .and. chisqmincost .lt. chisq_out_thresh) then  
              end_flag=2		! Use older value with lowest cost function
          endif
        ELSEIF (sy_flag .EQ. 3) THEN
          if (n .gt. 15 .and. chisq .gt. mincost .and. chisqmincost .lt. chisq_out_thresh) then
              end_flag=2                ! Use older value with lowest cost function
          endif
        ENDIF

          if (disq(1,1) .lt. 0.0) then ! should be impossible?
            print*,'convergence error -- negative value not allowed!'
            
            !if (write_test) stop 888
            end_flag=2
          endif
        endif  !if (n .gt. 1)

        if (end_flag == 1) then
          flag=1
          niter = n
	      !write(*,*),'Final Tbsim: ',tbsim
	      !write(*,*),'Final x: ',x2
	      !write(*,*),'Final chisq: ',chisqtot(1,1)
          return
        endif
	    if (end_flag .eq. 2) then
	      flag=1
	      niter = n
	      !tbsim = Fmincost
              sim_obs = Fmincost
	      x2 = x2mincost
	      chisq = chisqmincost
              write(*,*) 'Jumped out after ', n, ' iterations.'
              write(*,*) 'chisqmincost: ', chisq / nobs
              if (write_test) then
                write(*,*) 'Jumped out after ', n, ' iterations.' 
                write(*,*) 'chisqmincost: ', chisq / nobs
                stop 888
              endif
	      return
	    endif
	 
        Fold(:) = Fout(:)
      enddo  ! end of do loop starting way up above


      write(*,*) 'Solution not found in ', nn, ' steps'
      IF (write_test) THEN
        STOP 888
      ENDIF

      ! No solution found upon reaching max iterations?  Set the following to -99.99
      do a = 1,nvar
        Amatrix2(a) = -99.99
        sigma2(a)   = -99.99
      enddo
      chisq      = -99.99
      flag       = 0
      niter      = n
      return

      end subroutine opt_est

      subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
        implicit none 

        integer :: n
        real*8  :: a(n,n), c(n,n)
        real*8  :: aa(n,n),cc(n,n)
        real*8  :: L(n,n), U(n,n), b(n), d(n), x(n)
        real*8  :: coeff
        integer :: i, j, k

        ! step 0: initialization for matrices L and U and b
        ! Fortran 90/95 aloows such operations on matrices

        i=0.0
        U=0.0
        b=0.0

        aa(:,:) = dble(a(:,:)) ! do double precision!

        ! step 1: forward elimination

        do k=1, n-1
          do i=k+1,n
            coeff=aa(i,k)/aa(k,k)
            L(i,k) = coeff
            do j=k+1,n
              aa(i,j) = aa(i,j)-coeff*aa(k,j)
            enddo
          enddo
        enddo

        ! Step 2: prepare L and U matrices 
        ! L matrix is a matrix of the elimination coefficient
        ! + the diagonal elements are 1.0

        do i=1,n
          L(i,i) = 1.0
        enddo

        ! U matrix is the upper triangular part of A

        do j=1,n
          do i=1,j
            U(i,j) = aa(i,j)
          enddo
        enddo

        ! Step 3: compute columns of the inverse matrix C

        do k=1,n
          b(k)=1.0
          d(1) = b(1)

          ! Step 3a: Solve Ld=b using the forward substitution

          do i=2,n
            d(i)=b(i)
            do j=1,i-1
              d(i) = d(i) - L(i,j)*d(j)
            enddo
          enddo

          ! Step 3b: Solve Ux=d using the back substitution

          x(n)=d(n)/U(n,n)
          do i = n-1,1,-1
            x(i) = d(i)
            do j=n,i+1,-1
              x(i)=x(i)-U(i,j)*x(j)
            enddo
            x(i) = x(i)/u(i,i)
          enddo

          ! Step 3c: fill the solutions x(n) into column k of C

          do i=1,n
            cc(i,k) = x(i)
          enddo
          b(k) = 0.0
        enddo
        !c(:,:) = real(cc(:,:))
        c(:,:) = cc(:,:)
      end subroutine inverse

!---------------------------------------------------------------------------------
! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
function inv(A) result(Ainv)
  real(8), dimension(:,:), intent(in) :: A
  real(8), dimension(size(A,1),size(A,2)) :: Ainv

  real(8), dimension(size(A,1)) :: work  ! work array for LAPACK
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer :: n, info

  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  n = size(A,1)

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
end function inv

!--------------------------------------------------------------------------------
      end module optest
