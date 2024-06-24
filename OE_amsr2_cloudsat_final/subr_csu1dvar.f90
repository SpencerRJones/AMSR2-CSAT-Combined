MODULE subr_csu1dvar

      !---------------------------------------------------------------
      ! contains subroutines for reading in pre-processor files,
      ! reading and processing ancillary data
      ! added geos5_pp read routine and updated to f90 syntax 12/14/16
      !
      ! Modified to include only necessary components: Spencer Jones 05/2023
      !
      ! Modified to read AMSR2-CloudSat preprocessor output files: Spencer Jones 04/2024
      !---------------------------------------------------------------

USE define_csu1dvar
USE NETCDF


IMPLICIT NONE

CONTAINS

!-----------------------------------------------------------------------


SUBROUTINE read_pp(input_file)

CHARACTER(LEN=256) :: input_file
INTEGER :: iost, ipix


write(*,*) input_file

OPEN(UNIT=101, FILE=input_file, STATUS='OLD', &
     ACCESS='STREAM', FORM='UNFORMATTED', IOSTAT=iost)


IF (iost .NE. 0) THEN
        WRITE(*,*) 'Error opening pp file.'
        STOP
ENDIF


READ(101) npix

WRITE(*,*) 'Number of pixels in:', npix

READ(101) nlyrs

IF (nlyrs .NE. nz) THEN
        WRITE(*,*) 'Error. Number of layers in pp file is not correct.'
        STOP
ENDIF

ALLOCATE(hlev_in(nlevs))
ALLOCATE(dtime_in(npix,6))
ALLOCATE(lat_in(npix))
ALLOCATE(lon_in(npix))
ALLOCATE(tbobs_in(npix,nchans))
ALLOCATE(eia_in(npix,nchans))
ALLOCATE(refl_in(npix,nlyrs))
ALLOCATE(tpw_in(npix))
ALLOCATE(clwp_in(npix))
ALLOCATE(ciwp_in(npix))
ALLOCATE(sst_in(npix))
ALLOCATE(rey_sst_in(npix))
ALLOCATE(rss_sst_in(npix))
ALLOCATE(t2m_in(npix))
ALLOCATE(d2m_in(npix))
ALLOCATE(wspd_in(npix))
ALLOCATE(slp_in(npix))
ALLOCATE(tlev_in(npix,nlevs))
ALLOCATE(plev_in(npix,nlevs))
ALLOCATE(rmix_in(npix,nlevs))
ALLOCATE(frzlvl_in(npix))
ALLOCATE(frzlvl_bin_in(npix))
ALLOCATE(cldbse_in(npix))
ALLOCATE(cldbse_bin_in(npix))
ALLOCATE(lcl_in(npix))
ALLOCATE(qflg_in(npix))
ALLOCATE(sfc_in(npix))
ALLOCATE(modiscf_in(npix,31))
ALLOCATE(prcp_in(npix))
ALLOCATE(snrt_sfc_in(npix))

READ(101) hlev_in,       &  !Height bins
          dtime_in,      &  !Scan Time
          lat_in,        &  !Latitude
          lon_in,        &  !Longitude
          tbobs_in,      &  !Observed AMSR2 Tbs
          eia_in,        &  !AMSR2 EIA
          refl_in,       &  !Observed CloudSat reflectivities
          tpw_in,        &  !ERA5 TPW
          clwp_in,       &  !ERA5 cloud water path
          ciwp_in,       &  !ERA5 cloud ice path
          sst_in,        &  !ERA5 SST
          rey_sst_in,    &  !Reynolds SST (if available)
          rss_sst_in,    &  !RSS SST      (if available)
          t2m_in,        &  !ERA5 2m temperature
          d2m_in,        &  !ERA5 2m dewpoint
          wspd_in,       &  !ERA5 2m wind speed
          slp_in,        &  !ERA5 sea level pressure
          tlev_in,       &  !ERA5 temperature profile
          plev_in,       &  !ERA5 pressure profile
          rmix_in,       &  !ERA5 humidity profile (specific humidity q)
          frzlvl_in,     &  !ERA5 freezing level
          frzlvl_bin_in, &  !Height bin of freezing level
          cldbse_in,     &  !Cloud base from lowest ERA5 cloud layer
          cldbse_bin_in, &  !Height bin of cloud base
          lcl_in,        &  !Calculated lifted condensation level
          qflg_in,       &  !Preprocessor quality flag
          sfc_in,        &  !Surface type
          modiscf_in,    &  !MODIS Cloud flag from 2B-GEOPROF
          prcp_in,       &  !2C-RAIN-PROFILE surface precip rate (30 km average)
          snrt_sfc_in       !2C-SNOW-PROFILE surface precip rate (30 km average)

CLOSE(101)



RETURN

END SUBROUTINE read_pp

!------------------------------------------------------------------------------





!--------------------------------------------------------------------------------------------------

      subroutine prep_oe

      implicit none

      real :: freq,temp,pres,rho
      integer :: bin,era_lun1,era_lun2,era_lun3,rsslun,ios,syl,i,j,k
      character(len=4) :: yr
      character(len=2) :: mos(12)=(/ '1', '2', '3', '4', '5', '6', &
                                     '7', '8', '9','10','11','12' /)
      character(len=2) :: das(31)=(/ '1', '2', '3', '4', '5', '6', &
                                     '7', '8', '9','10','11','12', &
                                    '13','14','15','16','17','18', &
                                    '19','20','21','22','23','24', &
                                    '25','26','27','28','29','30','31'/)
      character(len=2) :: mo, da
      character(len=1) :: spc
      real :: sim_toff(nbins,nang,nch)
      CHARACTER(LEN=100) :: offset_file




      ALLOCATE(toffsets(11))
      ALLOCATE(s_toffsets(nbins,11))
      ALLOCATE(s_sy(11,11,nbins))
      ALLOCATE(sy_dummy(11,11))
      allocate(sy(nobs,nobs))
      allocate(sy_i(nobs,nobs))
      ALLOCATE(sa(nvar,nvar)) 

      write(spc,'(I1)') npc

       

      !Read in Duncan's Tb offsets. Modified in optest
        offset_file = 'binary/AM2.SY.Fas.3.11.v2.1.bin'

        OPEN(UNIT=550, FILE=TRIM(offset_file), ACCESS='STREAM', STATUS='OLD', &
             FORM='UNFORMATTED')
        READ(550) sy_dummy
        READ(550) toffsets
        READ(550) s_sy
        READ(550) s_toffsets
        CLOSE(550)




     
      !--- allocate and read in EOF LUT data

      allocate(era_wm(nlon,nlat))
      allocate(era_ws(nlon,nlat))
      allocate(mprof(nbins,nz))
      allocate(eofs(nbins,nz,6)) ! npc
      allocate(peofs(nz,6))      ! npc
      



        !Read in interpolated wv EOF arrays

      erafile = trim(bin_dir) // 'eof_mr.03.v3_30lyrs.bin' ! v3 has consistency across bins
      call gprof_lun(era_lun1)
      open(unit=era_lun1,file=erafile,access='stream',status='old',form='unformatted',iostat = ios)
      if (ios .ne. 0) write(*,*)' unable to open EOF LUT file', erafile
        read(era_lun1) mprof
        read(era_lun1) eofs
      close(era_lun1)




!---Read in Sy:

        !CALL read_sy



!---Read in Sa:

        CALL read_sa


        CALL alloc

RETURN


      end subroutine prep_oe


!--------------------------------------------------------------------------------------------------


      subroutine alloc

      ! --- Allocate variables used in read procedure from pp (no longer
      !     separate in 'read_l1*' subroutine.

!      allocate(dtime_in(nscans,6))
!      allocate(sclat_in(nscans))
!      allocate(sclon_in(nscans))
!      allocate(scalt_in(nscans))
!      allocate(lat_in(npix,nscans))
!      allocate(lon_in(npix,nscans))
!      allocate(land_in(npix,nscans,4))
!      allocate(sfc_in(npix,nscans))      
!      allocate(tpw_in(npix,nscans))
!      allocate(tsfc_in(npix,nscans))
!      allocate(wspd_in(npix,nscans))
!      allocate(slp_in(npix,nscans))
!      allocate(wdir_in(npix,nscans))
!      allocate(eia_in(npix,nscans,maxchans))
!      allocate(alpha_in(npix,nscans))
!      allocate(azim_in(npix,nscans))
!      allocate(tbobs_in(npix,nscans,maxchans))
!      allocate(tlev_in(npix,nscans,nz+1))
!      !allocate(hlev_in(npix,nscans,nz+1))
!      ALLOCATE(hlev_in(nz+1))
!      allocate(rmix_in(npix,nscans,nz)) ! different from rmix_out
!      allocate(sglint_in(npix,nscans))
!
!      allocate(tai93time(nscans))
!      allocate(sfc_type(npix,nscans))
!      allocate(lo_flag(npix,nscans))
!      allocate(sst(npix,nscans))
!      allocate(tavg_out(npix,nscans,nz))
!      allocate(rmix_out(npix,nscans,nz))
!      allocate(clwc_out(npix,nscans,nz))
!      allocate(ciwc_out(npix,nscans,nz))
!      allocate(tbdif_out(npix,nscans,nch))
!      allocate(tbsim_out(npix,nscans,nch))   ! all
!      allocate(iter_out(npix,nscans))
!      allocate(saln_out(npix,nscans))
!      allocate(post_out(npix,nscans,nvar+1)) 
!      allocate(incr_out(npix,nscans,4))   ! save pp TPW/Wind/LWP if desired
!      allocate(oe_output(npix,nscans,9))
!      allocate(screen(npix,nscans))
!


ALLOCATE(sfc_type(npix))
ALLOCATE(tavg_out(npix,nz))
ALLOCATE(saln_out(npix))
ALLOCATE(incr_out(npix,4))
ALLOCATE(rmix_out(npix,nz))
ALLOCATE(clwc_out(npix,nz))
ALLOCATE(plwc_out(npix,nz))
ALLOCATE(tiwc_out(npix,nz))
!ALLOCATE(z_e(npix,nz,nact_freqs))
!ALLOCATE(z_a(npix,nz,nact_freqs))
!ALLOCATE(atten(npix))
ALLOCATE(oe_output(npix,nvar+2))
ALLOCATE(post_out(npix,nvar+1))
ALLOCATE(iter_out(npix))
ALLOCATE(tbdif_out(npix,nch))
ALLOCATE(tbsim_out(npix,nch))
ALLOCATE(cltop_bin(npix))
ALLOCATE(clbot_bin(npix))
ALLOCATE(sfcprcp_out(npix))
ALLOCATE(rain_rate(npix,nz))
ALLOCATE(snow_rate(npix,nz))
ALLOCATE(retr_flag_out(npix))


        return
      end subroutine alloc


      subroutine dealloc

      ! --- Deallocate variables used in output procedure, now ncdf

      deallocate(dtime_in)
      deallocate(sclat_in)
      deallocate(sclon_in)
      deallocate(scalt_in)
      deallocate(lat_in)
      deallocate(lon_in)
      deallocate(eia_in)
      deallocate(alpha_in)
      deallocate(azim_in)
      deallocate(tbobs_in)
      deallocate(tlev_in)
      deallocate(hlev_in)
      deallocate(slp_in)
      deallocate(wdir_in)
      deallocate(land_in)
      deallocate(sfc_in) 
      deallocate(tpw_in)
      deallocate(tsfc_in)
      deallocate(wspd_in)
      deallocate(rmix_in)

      deallocate(sst)
      deallocate(oe_output)
      deallocate(screen)
      deallocate(tai93time)
      deallocate(tbdif_out)
      deallocate(tbsim_out)
      deallocate(post_out)
      deallocate(rmix_out)
      deallocate(clwc_out)
      deallocate(tiwc_out)
      deallocate(saln_out)
      deallocate(sfc_type)
      deallocate(lo_flag)
      deallocate(iter_out)
      deallocate(incr_out)
      deallocate(era_wm)
      deallocate(era_ws)
      deallocate(eofs)
      deallocate(mprof,peofs)
      deallocate(sy)
      deallocate(sy_i)
!      deallocate(s_sy)
!      deallocate(toffsets)
!      deallocate(s_toffsets)

      end subroutine dealloc

      subroutine gprof_lun(ilun)

      !--- This routine gets an open logical unit number starting with 100

      integer  :: ilun
      logical  :: llun
       
      do ilun = 100,201
        inquire(unit=ilun, opened=llun)
        if (.not. llun) exit
      enddo
      return
      end subroutine gprof_lun


!------------------------------------------------------------------------------------

!--------------------------------------------------
!       Reads matrix Sy from ascii file
!
!       -Spencer Jones, CSU Atmos, 04/2023
!--------------------------------------------------




SUBROUTINE read_sy(syflg,rdrstat)

INTEGER :: syflg, rdrstat
CHARACTER(LEN=8)   :: sy_el_char
CHARACTER(LEN=(100)) :: line
INTEGER :: iost, iobs, jobs
REAL(8)    :: sy_el


IF (rdrstat .EQ. 0) THEN
        IF (syflg .EQ. 1) sy_file = 'sy/sy_clear.mtrx'
        IF (syflg .EQ. 2) sy_file = 'sy/sy_cloudy.mtrx'
        IF (syflg .EQ. 3) sy_file = 'sy/sy_prcp.mtrx'
ELSEIF (rdrstat .EQ. 1) THEN
        IF (syflg .EQ. 1) sy_file = 'sy/sy_clear.mtrx'
        IF (syflg .EQ. 2) sy_file = 'sy/sy_cloudy_noradar.mtrx'
        IF (syflg .EQ. 3) sy_file = 'sy/sy_prcp_noradar.mtrx'
ELSEIF (rdrstat .EQ. 2) THEN
        IF (syflg .EQ. 1) sy_file = 'sy/sy_clear.mtrx'
        IF (syflg .EQ. 2) sy_file = 'sy/sy_cloudy_radaronly.mtrx'
        IF (syflg .EQ. 3) sy_file = 'sy/sy_prcp_radaronly.mtrx'
ELSE
        WRITE(*,*) 'Invalid value for radar_stat.'
        STOP
ENDIF


sy(:,:) = -9999.9


OPEN(UNIT=600, FILE=sy_file, STATUS='OLD', ACCESS='STREAM', &
     FORM='UNFORMATTED', IOSTAT=iost)
        IF (iost .NE. 0) THEN
                WRITE(*,*) 'Unable to open Sy file.'
                STOP
        ENDIF


READ(600) sy


!IF (debug_mode .EQ. .TRUE.) THEN
!        WRITE(*,*) ''
!        WRITE(*,*) 'Sy matrix: '
!        WRITE(*,'(10F10.2)') sy
!        WRITE(*,*) ''
!ENDIF



!Error check
DO iobs = 1, nobs
        DO jobs = 1, nobs
                IF (sy(iobs,jobs) .EQ. -9999.9) THEN
                        WRITE(*,*) 'Missing Sy element.'
                        STOP
                ENDIF
        ENDDO
ENDDO

CLOSE(600)

RETURN

END SUBROUTINE read_sy


!-----------------------------------------------------------------------------------

SUBROUTINE read_sa

!---------------------------------------------------
!
!       Similar to read_sy.
!
!       Spencer Jones, CSU Atmos., 04/2023
!
!---------------------------------------------------

CHARACTER(LEN=8)   :: sa_el_char
CHARACTER(LEN=(100)) :: line
INTEGER :: iost, ivar, jvar
REAL    :: sa_el


sa_file = 'sa/sa.mtrx'

sa = -9999.9


OPEN(UNIT=700, FILE=sa_file, STATUS='OLD', ACCESS='STREAM', &
     FORM='UNFORMATTED', IOSTAT=iost)
        IF (iost .NE. 0) THEN
                WRITE(*,*) 'Unable to open Sa file.'
                STOP
        ENDIF


READ(700) sa


IF (debug_mode .EQ. .TRUE.) THEN
        WRITE(*,*) ''
        WRITE(*,*) 'Sa matrix:'
        DO ivar = 1, nvar
                WRITE(*,'(20F7.2)') sa(ivar,ivar)
        ENDDO
        WRITE(*,*) ''
ENDIF

CLOSE(700)


RETURN


END SUBROUTINE read_sa

!----------------------------------------------------------------------------------


!--------------------------------------------
!
!    Ice lookup tables by Spencer Jones, 04/2024
!
!--------------------------------------------

subroutine read_ice_lut

        character(len=5)  :: freqs(6) = (/ '10.65', '18.70', '23.80', '36.50', '89.00', '94.00' /)
        character(len=100) lut_file
        integer :: ifrq
        real, allocatable :: icelut(:,:,:,:)
        integer :: ichan

        nfreq_icelut = 6

        !First time through, just get sizes for allocation
        do ifrq = 1, 1
                

                lut_file = 'LUT/ice/ice_LUT_' // trim(freqs(ifrq)) // '.bin'

                open(unit=998, file=lut_file, status='old', access='stream', &
                     form='unformatted')
                
                read(998) ntemp_icelut
                read(998) niwc_icelut
                read(998) nrho_icelut

                allocate(temp_icelut(ntemp_icelut))
                allocate(iwc_icelut(niwc_icelut))
                allocate(rho_icelut(nrho_icelut))


                close(998)

        enddo

        allocate(icelut(ntemp_icelut, niwc_icelut, nrho_icelut, 4))
        allocate(ksca_icelut(nfreq_icelut, ntemp_icelut, niwc_icelut, nrho_icelut))
        allocate(asca_icelut(nfreq_icelut, ntemp_icelut, niwc_icelut, nrho_icelut))
        allocate(gsca_icelut(nfreq_icelut, ntemp_icelut, niwc_icelut, nrho_icelut))
        allocate(pbck_icelut(nfreq_icelut, ntemp_icelut, niwc_icelut, nrho_icelut))


        do ifrq = 1, nfreq_icelut

                
                lut_file = 'LUT/ice/ice_LUT_' // trim(freqs(ifrq)) // '.bin'

                open(unit=998, file=lut_file, status='old', access='stream', &
                     form='unformatted')


                read(998) ntemp_icelut
                read(998) niwc_icelut
                read(998) nrho_icelut
                read(998) temp_icelut
                read(998) iwc_icelut
                read(998) rho_icelut


                read(998) icelut


                ksca_icelut(ifrq,:,:,:) = icelut(:,:,:,1)
                asca_icelut(ifrq,:,:,:) = icelut(:,:,:,2)
                gsca_icelut(ifrq,:,:,:) = icelut(:,:,:,3)
                pbck_icelut(ifrq,:,:,:) = icelut(:,:,:,4)

                close(998)


        enddo

        stepsize_T = temp_icelut(2) - temp_icelut(1)
        stepsize_iwc = iwc_icelut(2) - iwc_icelut(1)
        stepsize_rho = rho_icelut(2) - rho_icelut(1)

        return

end subroutine read_ice_lut


!-----------------------------------------------------------------------------------      


SUBROUTINE output_nc(outfile)

        INTEGER, PARAMETER :: noutpix = 175   !Don't want to write out entire swath
        INTEGER, PARAMETER :: start_buff = 21            !Change based on swath width (how many pixels to cut off at swath beginning?)

!-----   Modified to make variables names more intuitive and file specs
!----           up-to-date with normal practices. -Spencer Jones 04/2023


CHARACTER(LEN=256), INTENT(IN) :: outfile


      integer, parameter :: CMPRS = 1 !level of compression!
      integer            :: ctime(8)
      integer            :: i, ipix, iscan, ich, ichan,cchan, itime, ieia
      integer            :: time_varid, sclat_varid, sclon_varid, scalt_varid, scorient_varid
      integer            :: lat_varid, lon_varid, press_varid, eia_varid, tbsim_varid, tbobs_varid
      integer            :: tbdif_varid, sal_varid, slp_varid, tai93_varid
      integer            :: tpw_varid, wsp_varid, lwp_varid, tlwp_varid, wdir_varid, sst_varid, reysst_varid
      integer            :: chi_varid, chan_varid, iwp_varid !, eof_varid
      integer            :: sun_varid, land_varid, qual_varid, niter_varid
      INTEGER            :: ppqflg_varid, hgt_varid, mu_varid, rhosn_varid, eof1_varid, eof2_varid, eof3_varid
      integer            :: wprof_varid, tprof_varid, lprof_varid, rprof_varid, iprof_varid, post_varid
      integer            :: rr_varid, sr_varid, sfcp_varid, retrflg_varid
      integer            :: ncid, pix_dims(1), time_dims(2), prof_dims(2)!, eof_dims(3),
      integer            :: eia_dims(2), tb_dims(2), chan_dims(2), post_dims(2), lev_dims(2), hgt_dims(1)
      integer            :: pix_dimid, time_dimid, scan_dimid, str_dimid, chan_dimid, eia_dimid
      integer            :: npc_dimid, var_dimid, post_dimid
      integer            :: nlev_dimid, nz_dimid
      integer            :: pix_dim, time_dim, str_dim, scan_dim, chan_dim, eia_dim
      integer            :: npc_dim, nz_dim, nlev_dim, post_dim, var_dim
      character(len=40)  :: sdate, sprior
      character(len=40)  :: sfreq
      character(len=10)  :: sgran,siter,smiss,sthresh,sthresh1, em_version
      character(len=6)   :: schan(maxchans)


      integer,   allocatable :: scantime(:,:)
      integer*1, allocatable :: scorientout(:)
      integer*1, allocatable :: qflag(:)
      integer*1, allocatable :: landout(:)
      integer*1, allocatable :: iterout(:)
      real, allocatable      :: sunout(:)

      real, allocatable      :: latout(:)
      real, allocatable      :: lonout(:)
      real, allocatable      :: tpwout(:)
      real, allocatable      :: wspout(:)
      real, allocatable      :: lwpout(:)
      real, allocatable      :: iwpout(:)
      real, allocatable      :: sstout(:)
      real, allocatable      :: chiout(:)
      REAL, ALLOCATABLE      :: muout(:)
      REAL, ALLOCATABLE      :: sndensout(:)
      REAL, ALLOCATABLE      :: eof1out(:)
      REAL, ALLOCATABLE      :: eof2out(:)
      REAL, ALLOCATABLE      :: eof3out(:)

      REAL, ALLOCATABLE      :: hgt_out(:)

      real, allocatable      :: reyout(:,:)
      real, allocatable      :: salout(:,:)
      real, allocatable      :: slpout(:)
      integer*2, allocatable :: dirout(:,:)
      real, allocatable      :: post(:,:)
      real, allocatable      :: eia(:,:)
      real, allocatable      :: tbsim(:,:)
      real, allocatable      :: tbobs(:,:)
      real, allocatable      :: tbdif(:,:)
      real, allocatable      :: eofout(:,:)
      real, allocatable      :: wvp_prof(:,:)
      real, allocatable      :: tmp_prof(:,:)
      real, allocatable      :: lwp_prof(:,:)
      REAL, ALLOCATABLE      :: plw_prof(:,:)
      !real, allocatable      :: iwp_prof(:,:)
      REAL, ALLOCATABLE      :: ice_prof(:,:)
      REAL, ALLOCATABLE      :: rr_out(:,:)
      REAL, ALLOCATABLE      :: sr_out(:,:)
      REAL, ALLOCATABLE      :: prcp_out(:)

        INTEGER(KIND=1), ALLOCATABLE :: pp_qflg(:)

      !---Allocate arrays---

      allocate (scantime(npix,6))
      allocate (qflag(npix))
      allocate (landout(npix))
      allocate (iterout(npix))
      allocate (sunout(npix))

      allocate (latout(npix))
      allocate (lonout(npix))
      allocate (tpwout(npix))
      allocate (wspout(npix))
      !allocate (dirout(npix))
      allocate (lwpout(npix))
      allocate (iwpout(npix))
      allocate (sstout(npix))
      allocate (chiout(npix))
      ALLOCATE (muout(npix))
      ALLOCATE (sndensout(npix))
      ALLOCATE (eof1out(npix))
      ALLOCATE (eof2out(npix))
      ALLOCATE (eof3out(npix))

      !allocate (salout(npix))
      allocate (slpout(npix))
      allocate (post(npix,4))
      allocate (tbsim(npix,nch))
      allocate (tbobs(npix,nch))
      allocate (tbdif(npix,nch))
      allocate (eia(npix,nch))
      allocate (wvp_prof(npix,nz))
      allocate (tmp_prof(npix,nz))
      allocate (lwp_prof(npix,nz))
      ALLOCATE (plw_prof(npix,nz))
      !allocate (iwp_prof(npix,nz))
      ALLOCATE (ice_prof(npix,nz))
      ALLOCATE (pp_qflg(npix))
      !ALLOCATE (hgt_out(nz+1))
      ALLOCATE (rr_out(npix,nz))
      ALLOCATE (sr_out(npix,nz))
      ALLOCATE (prcp_out(npix))


            !---Copy data into output arrays---


       ! hgt_out = hlev_in(1,:)

      do ipix = 1, npix

        do itime=1,6
          scantime(ipix,itime) = dtime_in(ipix,itime)
        enddo


                latout(ipix)      = lat_in(ipix)
                lonout(ipix)      = lon_in(ipix)
                !tpwout(ipix)      = ftrunc(oe_output(ipix,1))
                if (oe_output(ipix,nvar-5) .eq. miss_flt) then
                        lwpout(ipix) = miss_flt
                else
                        lwpout(ipix)      = 10.**(oe_output(ipix,nvar-5))
                endif
                wspout(ipix)      = ftrunc(oe_output(ipix,nvar), 1000.)
                chiout(ipix)      = ftrunc(oe_output(ipix,nvar+1), 1000.)
                sstout(ipix)      = ftrunc(oe_output(ipix,nvar-1), 100.)
                tpwout(ipix)      = ftrunc(oe_output(ipix,nvar+2), 100.)
                pp_qflg(ipix)     = qflg_in(ipix)
                muout(ipix)       = ftrunc(oe_output(ipix,nvar-4), 100.)
                sndensout(ipix)   = ftrunc(oe_output(ipix,nvar-3), 100.)
                eof1out(ipix)     = ftrunc(oe_output(ipix,nvar-2), 100.)
                eof2out(ipix)     = -9999.9 !ftrunc(oe_output(ipix,nvar-3), 100.)
                eof3out(ipix)     = -9999.9 !ftrunc(oe_output(ipix,nvar-2), 100.)
                prcp_out(ipix)    = ftrunc(sfcprcp_out(ipix), 100.)

                !write(*,*) ipix, prcp_out(ipix)

                DO i = 1, nz
                        wvp_prof(ipix,i) = ftrunc(rmix_out(ipix,i), 100.)
                        tmp_prof(ipix,i) = ftrunc(tavg_out(ipix,i), 100.)
                        lwp_prof(ipix,i) = ftrunc(clwc_out(ipix,i), 1000.)
                        plw_prof(ipix,i) = ftrunc(plwc_out(ipix,i), 1000.)
                        ice_prof(ipix,i) = ftrunc(tiwc_out(ipix,i), 1000.)
                        rr_out(ipix,i)   = ftrunc(rain_rate(ipix,i), 1000.)
                        sr_out(ipix,i)   = ftrunc(snow_rate(ipix,i), 1000.)
                        
                ENDDO

                slpout(ipix) = ftrunc(slp_in(ipix), 100.)
                iterout(ipix) = iter_out(ipix)
                
                DO ichan = 1, nch
                        tbobs(ipix,ichan) = ftrunc(tbobs_in(ipix,ichan), 100.)
                        tbsim(ipix,ichan) = ftrunc(tbsim_out(ipix,ichan), 100.)
                        tbdif(ipix,ichan) = ftrunc(tbdif_out(ipix,ichan), 100.)
                        eia(ipix,ichan)   = ftrunc(eia_in(ipix,ichan), 100.)
                ENDDO
               
!          latout(iscan,ipix)     = lat_in(ipix+start_buff,iscan)
!          lonout(iscan,ipix)     = lon_in(ipix+start_buff,iscan)
!          tpwout(iscan,ipix)     = ftrunc(oe_output(ipix+start_buff,iscan,1), 100.0)
!          lwpout(iscan,ipix)     = ftrunc(oe_output(ipix+start_buff,iscan,2),1000.0)
!          wspout(iscan,ipix)     = ftrunc(oe_output(ipix+start_buff,iscan,3), 100.0)
!          dirout(iscan,ipix)     = int(wdir_in(ipix+start_buff,iscan))
!          iwpout(iscan,ipix)     = ftrunc(oe_output(ipix+start_buff,iscan,4),1000.0)
!          chiout(iscan,ipix)     = ftrunc(oe_output(ipix+start_buff,iscan,5),1000.0)
!          sstout(iscan,ipix)     = ftrunc(oe_output(ipix+start_buff,iscan,6), 100.0)
!          post(iscan,ipix,1)     = ftrunc(post_out(ipix+start_buff,iscan,8), 100.0)
!          post(iscan,ipix,2)     = ftrunc(post_out(ipix+start_buff,iscan,4), 100.0)
!          post(iscan,ipix,3)     = ftrunc(post_out(ipix+start_buff,iscan,5), 100.0)
!          post(iscan,ipix,4)     = ftrunc(post_out(ipix+start_buff,iscan,7), 100.0)
!          do i=1,nz
!            wvp_prof(iscan,ipix,i) = ftrunc(rmix_out(ipix+start_buff,iscan,i), 100.0)
!            tmp_prof(iscan,ipix,i) = ftrunc(tavg_out(ipix+start_buff,iscan,i), 100.0)
!            lwp_prof(iscan,ipix,i) = ftrunc(clwc_out(ipix+start_buff,iscan,i), 1000.0)
!            iwp_prof(iscan,ipix,i) = ftrunc(ciwc_out(ipix+start_buff,iscan,i), 1000.0)
!          enddo

!          salout(iscan,ipix)       = ftrunc(saln_out(ipix+start_buff,iscan),100.0)
!          slpout(iscan,ipix)       = ftrunc(slp_in(ipix+start_buff,iscan),100.0)
!          landout(iscan,ipix)      = lo_flag(ipix+start_buff,iscan)
!          iterout(iscan,ipix)      = iter_out(ipix+start_buff,iscan)
!          cchan=1
!          do ichan=1,maxchans !nch
!            if (chan_avail(ichan) .eq. 1) then
!              tbobs(iscan,ipix,cchan) = ftrunc(tbobs_in(ipix+start_buff,iscan,ichan), 100.0)
!              tbsim(iscan,ipix,cchan) = ftrunc(tbsim_out(ipix+start_buff,iscan,cchan), 100.0)
!              tbdif(iscan,ipix,cchan) = ftrunc(tbdif_out(ipix+start_buff,iscan,cchan), 100.0)
!              eia(iscan,ipix,cchan)   = ftrunc(eia_in(ipix+start_buff,iscan,ichan), 100.0)
!              cchan=cchan+1
!            endif
!          enddo
      enddo



!      ich = 1
!      do ichan=1,maxchans
 !       if (chan_avail(ichan) .eq. 1) then
 !         schan(ich) = freqs(ichan)
 !         ich = ich + 1
 !       endif
 !     enddo



      !---open output file---

  call check(nf90_create(outfile,NF90_NETCDF4,ncid))

      !---specify diminsions---

      !pix_dim   = noutpix
      pix_dim   = npix
      !scan_dim  = nscans
      chan_dim  = nch
      time_dim  = 6
      str_dim   = 7
      eia_dim   = nch
      var_dim   = nvar
      post_dim  = 4 !nvar+1
      nz_dim    = nz
      nlev_dim  = nz+1

      call check(nf90_def_dim(ncid, "npix",    pix_dim,  pix_dimid))
      !call check(nf90_def_dim(ncid, "nscan",  scan_dim, scan_dimid))
      call check(nf90_def_dim(ncid, "nchan",  chan_dim, chan_dimid))
      call check(nf90_def_dim(ncid, "ntime",  time_dim, time_dimid))
      call check(nf90_def_dim(ncid, "strlen",  str_dim,  str_dimid))
      call check(nf90_def_dim(ncid, "neia",    eia_dim,  eia_dimid))
      call check(nf90_def_dim(ncid, "npost",  post_dim, post_dimid))
      call check(nf90_def_dim(ncid, "nvar",    var_dim,  var_dimid))
      call check(nf90_def_dim(ncid, "nlayer",   nz_dim,   nz_dimid))
      call check(nf90_def_dim(ncid, "nlevel", nlev_dim, nlev_dimid))
      time_dims = (/ pix_dimid, time_dimid /)
      pix_dims  = (/  pix_dimid /)
      prof_dims = (/  pix_dimid, nz_dimid /)
      lev_dims  = (/  pix_dimid, nlev_dimid /)
      post_dims = (/  pix_dimid, post_dimid /)
      eia_dims  = (/  pix_dimid, eia_dimid /)
      tb_dims   = (/  pix_dimid, chan_dimid /)
      chan_dims = (/ str_dimid,  chan_dimid /)
      hgt_dims  = (/ nlev_dimid /)

     !---write global metadata---

      write(sgran,'(i6.6)') granule
      write(siter,'(f5.2)') avg_iter
      write(smiss,'(f7.1)') miss_flt
      write(sthresh,'(f4.1)') chisq_out_thresh
      write(sthresh1,'(f4.1)') chis1
      call date_and_time(values=ctime)

      write(sdate,'(i4.4,"-",i2.2,"-",i2.2,"T",i2.2,":",i2.2,":",i2.2)') &
            ctime(1),ctime(2),ctime(3),ctime(5),ctime(6),ctime(7)
      sprior=trim("GEOS5")
      em_version = "FASTEM6"

      call check(nf90_put_att(ncid, NF90_GLOBAL, "Satellite",         trim(satellite)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "Granule",               trim(sgran)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "Version",        trim(code_version)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "RT_Model",        trim(rtm_version)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "Emis_Model",       trim(em_version)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "Preprocessor_Version",trim(ppversion)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "L1_File",           trim(orig_file)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "SY_File",             trim(sy_file)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "Missing_Value",         trim(smiss)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "ChiSquare_Threshold", trim(sthresh)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "Source_of_Prior",      trim(sprior)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "Average_Iterations",    trim(siter)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "Creation_Date",         trim(sdate)))


      !---declare variables---
      ! -- modified for ver2.3 to CF-compliant, D Duncan, 9/28/16
      call check(nf90_def_var(ncid, "ScTime",           NF90_SHORT,   time_dims,     time_varid))
      call check(nf90_def_var(ncid, "SCLatitude",       NF90_FLOAT,  pix_dimid,    sclat_varid))
      call check(nf90_def_var(ncid, "SCLongitude",      NF90_FLOAT,  pix_dimid,    sclon_varid))
      call check(nf90_def_var(ncid, "SCAltitude",       NF90_FLOAT,  pix_dimid,    scalt_varid))
      call check(nf90_def_var(ncid, "ChannelName",      NF90_CHAR,    chan_dims,     chan_varid))
      CALL check(NF90_DEF_VAR(ncid, "Height",           NF90_FLOAT,    hgt_dims,   hgt_varid))

      call check(nf90_def_var(ncid, "Latitude",          NF90_FLOAT,    pix_dims,      lat_varid))
      call check(nf90_def_var(ncid, "Longitude",         NF90_FLOAT,    pix_dims,      lon_varid))
      call check(nf90_def_var(ncid, "WaterVaporMass", NF90_FLOAT, pix_dims, tpw_varid))
      call check(nf90_def_var(ncid, "WindSpeed",        NF90_FLOAT,    pix_dims,      wsp_varid))
      call check(nf90_def_var(ncid, "CLW", NF90_FLOAT, pix_dims, lwp_varid))
      call check(nf90_def_var(ncid, "IWP",            NF90_FLOAT,    pix_dims,      iwp_varid))
      call check(nf90_def_var(ncid, "ChiSquared",       NF90_FLOAT,    pix_dims,      chi_varid))
      call check(nf90_def_var(ncid, "SeaSurfaceTemperature", NF90_FLOAT,    pix_dims,      sst_varid))
      call check(nf90_def_var(ncid, "NumberOfIterations",        NF90_BYTE,      pix_dims,    niter_varid))
      call check(nf90_def_var(ncid, "LandAreaFraction",NF90_BYTE,     pix_dims,     land_varid))
      call check(nf90_def_var(ncid, "PosteriorError",   NF90_FLOAT,   post_dims,     post_varid))
      call check(nf90_def_var(ncid, "TbDifference",     NF90_FLOAT,     tb_dims,    tbdif_varid))
      
      call check(nf90_def_var(ncid, "AirPressureProfile", NF90_SHORT,  lev_dims,    press_varid))
      call check(nf90_def_var(ncid, "EIA",          NF90_FLOAT,    eia_dims,      eia_varid))
      call check(nf90_def_var(ncid, "WindDirection",         NF90_SHORT,    pix_dims,     wdir_varid))
      call check(nf90_def_var(ncid, "WaterVaporProfile",     NF90_FLOAT,   prof_dims,    wprof_varid))
      call check(nf90_def_var(ncid, "TemperatureProfile",     NF90_FLOAT,   prof_dims,    tprof_varid))
      call check(nf90_def_var(ncid, "CloudWaterProfile",     NF90_FLOAT,   prof_dims,    lprof_varid))
      CALL check(NF90_DEF_VAR(ncid, "RainWaterProfile", NF90_FLOAT, prof_dims, rprof_varid))
      call check(nf90_def_var(ncid, "IceWaterProfile",     NF90_FLOAT,   prof_dims,    iprof_varid))
      call check(nf90_def_var(ncid, "RainRate", NF90_FLOAT, prof_dims, rr_varid))
      call check(nf90_def_var(ncid, "SnowRate", NF90_FLOAT, prof_dims, sr_varid))
      call check(nf90_def_var(ncid, "SurfacePrecip", NF90_FLOAT, pix_dims, sfcp_varid))
      call check(nf90_def_var(ncid, "TbsObs",        NF90_FLOAT,     tb_dims,    tbobs_varid))
      call check(nf90_def_var(ncid, "TbsSim",        NF90_FLOAT,     tb_dims,    tbsim_varid))
      call check(nf90_def_var(ncid, "SeaSurfaceSalinity",        NF90_FLOAT,    pix_dims,      sal_varid))
      call check(nf90_def_var(ncid, "SeaLevelPressure",          NF90_FLOAT,    pix_dims,      slp_varid))
      CALL check(NF90_DEF_VAR(NCID, "PPQuality", NF90_BYTE, pix_dims, ppqflg_varid))

      CALL check(NF90_DEF_VAR(ncid, "Mu_liq", NF90_FLOAT, pix_dims, mu_varid))
      CALL check(NF90_DEF_VAR(ncid, "IceParticleDensity", NF90_FLOAT, pix_dims, rhosn_varid))
      CALL check(NF90_DEF_VAR(ncid, "EOF1", NF90_FLOAT, pix_dims, eof1_varid))
      CALL check(NF90_DEF_VAR(ncid, "EOF2", NF90_FLOAT, pix_dims, eof2_varid))
      CALL check(NF90_DEF_VAR(ncid, "EOF3", NF90_FLOAT, pix_dims, eof3_varid))

      CALL check(NF90_DEF_VAR(ncid, 'RetrievalFlag', NF90_BYTE, pix_dims, retrflg_varid))

      !---compress variables---

      if (CMPRS .gt. 0) then
        call check(nf90_def_var_deflate(ncid,     time_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,    sclat_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,    sclon_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,    scalt_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,      lat_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,      lon_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,      tpw_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,      wsp_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,      lwp_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,      iwp_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,      chi_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,      sst_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,    niter_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,     land_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,     post_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,    tbdif_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,  press_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,    eia_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,   wdir_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,  wprof_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,  tprof_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,  lprof_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,  iprof_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,  tbobs_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,  tbsim_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,    sal_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,    slp_varid, 0, 1, CMPRS))
      endif

      call check(nf90_enddef(ncid))

      !---specify units---

      call check(nf90_put_att(ncid,   sclat_varid, "units", "degree_north"))
      call check(nf90_put_att(ncid,   sclon_varid, "units",  "degree_east"))
      call check(nf90_put_att(ncid,   scalt_varid, "units",           "km"))
      call check(nf90_put_att(ncid,     lat_varid, "units", "degree_north"))
      call check(nf90_put_att(ncid,     lon_varid, "units",  "degree_east"))
      call check(nf90_put_att(ncid,     tpw_varid, "units",       "kg m-2"))
      call check(nf90_put_att(ncid,     wsp_varid, "units",        "m s-1"))
      call check(nf90_put_att(ncid,     lwp_varid, "units",        "g m-2"))
      call check(nf90_put_att(ncid,     iwp_varid, "units",        "g m-2"))
      call check(nf90_put_att(ncid,     sst_varid, "units",            "K"))
      call check(nf90_put_att(ncid,    land_varid, "units",            "1"))
      call check(nf90_put_att(ncid,   tbdif_varid, "units",            "K"))
      call check(nf90_put_att(ncid, press_varid, "units",           "Pa"))
      call check(nf90_put_att(ncid,   eia_varid, "units",      "degrees"))
      call check(nf90_put_att(ncid,  wdir_varid, "units",      "degrees"))
      call check(nf90_put_att(ncid, wprof_varid, "units",         "g/kg"))
      call check(nf90_put_att(ncid, tprof_varid, "units",            "K"))
      call check(nf90_put_att(ncid, lprof_varid, "units", "kg/m^2/layer"))
      CALL check(NF90_PUT_ATT(ncid, rprof_varid, "units",           "g m^-3"))
      call check(nf90_put_att(ncid, iprof_varid, "units",         "g m^-3"))
      call check(nf90_put_att(ncid, rr_varid,    "units",         "mm hr^-1"))
      call check(nf90_put_att(ncid, sr_varid,    "units",         "mm hr^-1"))
      call check(nf90_put_att(ncid, sfcp_varid,  "units",         "mm hr^-1"))
      call check(nf90_put_att(ncid, tbobs_varid, "units",            "K"))
      call check(nf90_put_att(ncid, tbsim_varid, "units",            "K"))
      call check(nf90_put_att(ncid, sal_varid,   "units",          "psu"))
      call check(nf90_put_att(ncid, slp_varid,   "units",           "mb"))
      CALL check(NF90_PUT_ATT(ncid, hgt_varid,   "units",           "m"))
      CALL check(NF90_PUT_ATT(ncid, mu_varid,    "units",           ""))
      CALL check(NF90_PUT_ATT(ncid, rhosn_varid, "units",           "g cm^-3"))
      CALL check(NF90_PUT_ATT(ncid, eof1_varid,    "units",           ""))
      CALL check(NF90_PUT_ATT(ncid, eof2_varid,    "units",           ""))
      CALL check(NF90_PUT_ATT(ncid, eof3_varid,    "units",           ""))


      
      !---specify long name---

      call check(nf90_put_att(ncid,     time_varid, "long_name", "Scan time (year, month, day, hour, minute, second)"))
      call check(nf90_put_att(ncid,     chan_varid, "long_name", "Channel specification"))
      call check(nf90_put_att(ncid,    sclat_varid, "long_name", "Spacecraft latitude"))
      call check(nf90_put_att(ncid,    sclon_varid, "long_name", "Spacecraft longitude"))
      call check(nf90_put_att(ncid,    scalt_varid, "long_name", "Spacecraft altitude"))
      call check(nf90_put_att(ncid,      lat_varid, "long_name", "Latitude"))
      call check(nf90_put_att(ncid,      lon_varid, "long_name", "Longitude"))
      call check(nf90_put_att(ncid,      tpw_varid, "long_name", "Total precipitable water"))
      call check(nf90_put_att(ncid,      wsp_varid, "long_name", "Wind speed"))
      call check(nf90_put_att(ncid,      lwp_varid, "long_name", "Residual cloud Liquid water path"))
      call check(nf90_put_att(ncid,      iwp_varid, "long_name", "Ice water path"))
      call check(nf90_put_att(ncid,      chi_varid, "long_name", "Chi-squared"))
      call check(nf90_put_att(ncid,      sst_varid, "long_name", "Sea surface temperature"))
      call check(nf90_put_att(ncid,    niter_varid, "long_name", "Number of iterations"))
      call check(nf90_put_att(ncid,     land_varid, "long_name", "Percentage of land coverage"))
      call check(nf90_put_att(ncid,    post_varid, "long_name", "Posteriori errors (sigma)"))
      call check(nf90_put_att(ncid,    tbdif_varid, "long_name", "Observed - Simulated Tb"))

      call check(nf90_put_att(ncid,  press_varid, "long_name", "Level pressure"))
      call check(nf90_put_att(ncid,    eia_varid, "long_name", "Earth incidence angle"))
      call check(nf90_put_att(ncid,   wdir_varid, "long_name", "Wind direction"))
      call check(nf90_put_att(ncid,  wprof_varid, "long_name", "Water vapor mixing ratio"))
      call check(nf90_put_att(ncid,  tprof_varid, "long_name", "Mean layer temperature"))
      call check(nf90_put_att(ncid,  lprof_varid, "long_name", "Cloud liquid water content"))
      CALL check(NF90_PUT_ATT(ncid,  rprof_varid, "long_name", "Rain Water Content"))
      call check(nf90_put_att(ncid,  iprof_varid, "long_name", "Cloud ice water content"))
      CALL check(nf90_put_att(ncid,  rr_varid,    "long_name", "Rain rate"))
      CALL check(nf90_put_att(ncid,  sr_varid,    "long_name", "Liquid equivalent snowfall rate"))
      CALL check(nf90_put_att(ncid,  sfcp_varid,  "long_name", "Surface precipitation rate"))
      call check(nf90_put_att(ncid,  tbobs_varid, "long_name", "Observed brightness temperatures"))
      call check(nf90_put_att(ncid,  tbsim_varid, "long_name", "Simulated brightness temperatures"))
      call check(nf90_put_att(ncid,    sal_varid, "long_name", "Sea Surface Salinity"))
      call check(nf90_put_att(ncid,    slp_varid, "long_name", "Sea Level Pressure"))
      CALL check(NF90_PUT_ATT(NCID,  ppqflg_varid, "long_name", "Quality flag from preprocessor"))
      CALL check(NF90_PUT_ATT(ncid,  hgt_varid,   "long_name", "Height above sea surface"))
      CALL check(NF90_PUT_ATT(ncid,  mu_varid,    "long_name", "Gamma DSD shape parameter mu"))
      CALL check(NF90_PUT_ATT(ncid,  rhosn_varid,   "long_name", "Ice Particle Density"))
      CALL check(NF90_PUT_ATT(ncid,  eof1_varid,   "long_name", "Water Vapor Profile EOF1 coefficient"))
      CALL check(NF90_PUT_ATT(ncid,  eof2_varid,   "long_name", "Water Vapor Profile EOF2 coefficient"))
      CALL check(NF90_PUT_ATT(ncid,  eof3_varid,   "long_name", "Water Vapor Profile EOF3 coefficient"))
      CALL check(NF90_PUT_ATT(ncid, retrflg_varid, 'long_name', 'Retrieval flag: -99=no retrieval, -9=not ocean surface, -8=bad &
      quality, 1=clear retrieval, 2=cloudy retrieval, 3=precip retrieval, -1=failed clear, -2=failed cloudy, -3=failed precip'))


        !---write variables---

        
      call check(nf90_put_var(ncid,     time_varid,     scantime))
      call check(nf90_put_var(ncid,     chan_varid, schan(1:nch)))
      call check(nf90_put_var(ncid,    sclat_varid,     sclat_in))
      call check(nf90_put_var(ncid,    sclon_varid,     sclon_in))
      call check(nf90_put_var(ncid,    scalt_varid,     scalt_in))
      call check(nf90_put_var(ncid,      lat_varid,       latout))
      call check(nf90_put_var(ncid,      lon_varid,       lonout))
      call check(nf90_put_var(ncid,      tpw_varid,       tpwout))
      call check(nf90_put_var(ncid,      wsp_varid,       wspout))
      call check(nf90_put_var(ncid,      lwp_varid,       lwpout))
      call check(nf90_put_var(ncid,      iwp_varid,       iwpout))
      call check(nf90_put_var(ncid,      chi_varid,       chiout))
      call check(nf90_put_var(ncid,      sst_varid,       sstout))
      call check(nf90_put_var(ncid,    niter_varid,      iterout))
      call check(nf90_put_var(ncid,     land_varid,      landout))
      call check(nf90_put_var(ncid,     post_varid,         post))
      call check(nf90_put_var(ncid,    tbdif_varid,        tbdif))
      CALL check(nf90_put_var(ncid,   ppqflg_varid,      pp_qflg))
      CALL check(NF90_PUT_VAR(ncid,   hgt_varid,         hlev_in))
      
      call check(nf90_put_var(ncid,  press_varid,      plev_in))
      call check(nf90_put_var(ncid,    eia_varid,          eia))
      call check(nf90_put_var(ncid,   wdir_varid,       dirout))
      call check(nf90_put_var(ncid,  wprof_varid,     wvp_prof))
      call check(nf90_put_var(ncid,  tprof_varid,     tmp_prof))
      call check(nf90_put_var(ncid,  lprof_varid,     lwp_prof))
      CALL check(NF90_PUT_VAR(ncid,  rprof_varid,     plw_prof))
      call check(nf90_put_var(ncid,  iprof_varid,     ice_prof))
      call check(nf90_put_var(ncid,  rr_varid,        rr_out))
      call check(nf90_put_var(ncid,  sr_varid,        sr_out))
      call check(nf90_put_var(ncid,  sfcp_varid,      prcp_out))
      call check(nf90_put_var(ncid,  tbobs_varid,        tbobs))
      call check(nf90_put_var(ncid,  tbsim_varid,        tbsim))
      call check(nf90_put_var(ncid,    sal_varid,       salout))
      call check(nf90_put_var(ncid,    slp_varid,       slpout))
      CALL check(NF90_PUT_VAR(ncid,   mu_varid,         muout))
      CALL check(NF90_PUT_VAR(ncid,   rhosn_varid,      sndensout))
      CALL check(NF90_PUT_VAR(ncid,   eof1_varid,         eof1out))
      CALL check(NF90_PUT_VAR(ncid,   eof2_varid,         eof2out))
      CALL check(NF90_PUT_VAR(ncid,   eof3_varid,         eof3out))
      CALL check(NF90_PUT_VAR(ncid,  retrflg_varid,     retr_flag_out))

      !---close file---

      call check(nf90_close(ncid))


      !--- Deallocate arrays ---!
      !--- currently lots of variables are not deallocated ---!

      deallocate (scantime)
      deallocate (qflag)
      deallocate (landout)
      deallocate (sunout)
      deallocate (latout)
      deallocate (lonout)
      deallocate (tpwout)
      deallocate (wspout)
      !deallocate (dirout)
      deallocate (lwpout)
      deallocate (iwpout)
      deallocate (sstout)
      deallocate (chiout)
      !deallocate (salout)
      deallocate (slpout)
      deallocate (eia)
      deallocate (post)
      deallocate (tbobs)
      deallocate (tbsim)
      deallocate (tbdif)
      deallocate (wvp_prof)
      deallocate (tmp_prof)
      deallocate (lwp_prof)
      !deallocate (iwp_prof)
      
      deallocate(clbot_bin) 
      deallocate(cltop_bin)

      DEALLOCATE( hlev_in, dtime_in, lat_in, lon_in, tbobs_in, eia_in, &
                 refl_in, tpw_in, clwp_in, ciwp_in, sst_in, rey_sst_in, rss_sst_in, t2m_in, d2m_in,  wspd_in, &
                 slp_in, tlev_in, plev_in, rmix_in, frzlvl_in, frzlvl_bin_in, cldbse_in, cldbse_bin_in, &
                 lcl_in, qflg_in, sfc_in, modiscf_in, prcp_in, &
                 toffsets, s_toffsets, s_sy, sy_dummy, sy, sy_i, sa, &
                 era_wm, era_ws, mprof, eofs, peofs, &
                 sfc_type, tavg_out, saln_out, incr_out, rmix_out, clwc_out, &
                 plwc_out,  &
                 tiwc_out, oe_output, post_out, iter_out, tbdif_out, tbsim_out, &
                 ksca_icelut, asca_icelut, gsca_icelut, pbck_icelut, &
                 rain_rate, snow_rate, sfcprcp_out, retr_flag_out)


RETURN

 END SUBROUTINE output_nc


      subroutine check(status)

      integer, intent ( in) :: status

      if (status /= nf90_noerr) then
         print *, trim(nf90_strerror(status))
         stop "Stopped"
      end if

      end subroutine check

      function ftrunc(rval, rscale)

      real, intent (in) :: rval
      real, intent (in) :: rscale
      real              :: ftrunc

      if (rval .lt. 0.0) then
        ftrunc = real( int( rval * rscale - 0.5)) / rscale
      else
        ftrunc = real( int( rval * rscale + 0.5)) / rscale
      endif

      end function ftrunc






!-----------------------------------------------------------------------------------

END MODULE subr_csu1dvar
