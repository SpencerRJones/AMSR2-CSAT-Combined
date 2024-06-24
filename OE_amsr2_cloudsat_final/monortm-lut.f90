      MODULE MonoRTM

      USE define_csu1dvar
      implicit  none

      contains



!------------------------------------------------------------------


SUBROUTINE read_monortm_lut

!----------------------------------------------------
!
!   Mono-RTM gaseous absorption lookup tables
!      -Original tables by Wes Berg
!      -Modified and read routine by Spencer Jones
!      -Notes: These tables are very complicated and
!       are unintuitive to work with. The table form
!       is as follows: (record by record)
!
!          1. nfreqs     - Integer: number of unique frequencies for sensor
!          2. freq_lut   - Real:    Array of unique freqencies from lowest to highest
!          3. nchans     - Integer: number of radiometer channels (e.g. 13 - GMI)
!          4. ifreq_lut  - Integer: Array (size=nchans) of frequency indexers.
!                            -e.g. ifreq_lut = (/1,1,2,2/) indexes for 10.65, 10.65, 18.7, 18.7
!          5. ipol_lut   - Integer: Array (size=nchans) of polarity indexers.
!                            -e.g. ipol_lut = (/0,1,0,1/) indexes V, H, V, H
!          6. npres      - Integer: number of pressure levels
!          7. pres_lut   - Real:    Array (size=npres) of pressure levels in lut
!          8. ntemp      - Integer: number of temperature values
!          9. temp_lut   - Real:    Array (size=ntemp) of unique temperatures in lut
!         10. nrmix      - Integer: number of mixing ratio values
!         11. rmix_lut   - Real:    Array (size=nrmix) of mixing ratios in lut
!         12. kabs_lut   - Real:    Array (shape=(ntemp, npres, nrmix, nfreqs)) of k_ext values
!
!   Spencer Jones, CSU
!
!-----------------------------------------------------


        CHARACTER(LEN=100) :: lut_file_pas, lut_file_act
        INTEGER :: isens, nact_sensors, npas_sensors, freq_end, nf
        CHARACTER(LEN=8) :: sensor
        CHARACTER(LEN=1) :: sens_type
        INTEGER :: nfreqs


        npres_lut = 101
        ntemp_lut = 141
        nrmix_lut = 301

        nfreqs = 6 !5 unique AMSR2 frequencies and 1 CPR radar frequency

        ALLOCATE(pres_lut(npres_lut), temp_lut(ntemp_lut), rmix_lut(nrmix_lut))
        ALLOCATE(freq_lut(nfreqs), ifreq_lut(nchans+nact_freqs), ipol_lut(nchans+nact_freqs))

        lut_file_pas = 'LUT/AMSR2/MonoRTM_v5.3-AMSR2.tbl'
        lut_file_act = 'LUT/CPR/MonoRTM_v5.3-CPR.tbl'


        OPEN(UNIT=10, FILE=lut_file_pas, ACCESS='STREAM')
        READ(10) nfreq_lut_pas
        ALLOCATE(freq_lut_pas(nfreq_lut_pas))
        READ(10) freq_lut_pas
        
        READ(10) nchan_lut_pas
        ALLOCATE(ifreq_lut_pas(nchan_lut_pas))
        ALLOCATE(ipol_lut_pas(nchan_lut_pas))
        READ(10) ifreq_lut_pas
        READ(10) ipol_lut_pas

        READ(10) npres_lut
        READ(10) pres_lut
        READ(10) ntemp_lut
        READ(10) temp_lut
        READ(10) nrmix_lut
        READ(10) rmix_lut
        
        ALLOCATE(kabs_lut_pas(npres_lut, ntemp_lut, nrmix_lut, nfreq_lut_pas))
        READ(10) kabs_lut_pas

        CLOSE(10)



        OPEN(UNIT=10, FILE=lut_file_act, ACCESS='STREAM')
        READ(10) nfreq_lut_act
        ALLOCATE(freq_lut_act(nfreq_lut_act))
        READ(10) freq_lut_act

        READ(10) npres_lut
        READ(10) pres_lut
        READ(10) ntemp_lut
        READ(10) temp_lut
        READ(10) nrmix_lut
        READ(10) rmix_lut

        ALLOCATE(kabs_lut_act(npres_lut, ntemp_lut, nrmix_lut, nfreq_lut_act))
        READ(10) kabs_lut_act

        CLOSE(10)

        
        freq_lut(1:nfreq_lut_pas)        = freq_lut_pas
        freq_lut(nfreq_lut_pas+1:nfreqs) = freq_lut_act

        ifreq_lut(1:2)  = 1
        ifreq_lut(3:4)  = 2
        ifreq_lut(5:6)  = 3
        ifreq_lut(7:8)  = 4
        ifreq_lut(9:10) = 5
        ifreq_lut(11)   = 6 
 
        ipol_lut(1:nchans) = ipol_lut_pas
        ipol_lut(nchans+1:nchans+nact_freqs) = 0 
        
        !write(*,*) nfreq_lut_pas
        !write(*,*) freq_lut
        !write(*,*) ifreq_lut
        !write(*,*) ipol_lut
        
        ALLOCATE(kabs_lut(npres_lut, ntemp_lut, nrmix_lut, nfreqs))

        kabs_lut(:,:,:,1:nfreq_lut_pas)        = kabs_lut_pas
        kabs_lut(:,:,:,nfreq_lut_pas+1:nfreqs) = kabs_lut_act



RETURN


END SUBROUTINE read_monortm_lut


!--------------------------------------------------------------------

      subroutine monortm_lut(ichan, pavg, tavg, ravg,  kext)

!     Look-up table version of MonoRTM

!---Split active and passive k_ext tables -SJ


      integer :: ichan
      real    :: pavg
      real    :: tavg
      real    :: ravg
      real    :: kext
      CHARACTER :: stype

      integer :: freq_index
      integer :: i,j,k,n
      integer :: p1,p2
      integer :: t1,t2
      integer :: r1,r2
      real    :: pw1,pw2
      real    :: tw1,tw2
      real    :: rw1,rw2

      p1 = 1
      do n=1,npres_lut
        if (pres_lut(n) .gt. pavg) then
          p2 = n
          exit
        else
          p1 = n
          p2 = n
        endif
      enddo
      if (p1 .eq. p2) then
        pw1 = 0.5
        pw2 = 0.5
      else
        pw1 = (pres_lut(p2) - pavg) / (pres_lut(p2) - pres_lut(p1))
        pw2 = (pavg - pres_lut(p1)) / (pres_lut(p2) - pres_lut(p1))
      endif


      t1 = 1
      do n=1,ntemp_lut
        if (temp_lut(n) .gt. tavg) then
          t2 = n
          exit
        else
          t1 = n
          t2 = n
        endif
      enddo
      if (t1 .eq. t2) then
        tw1 = 0.5
        tw2 = 0.5
      else
        tw1 = (temp_lut(t2) - tavg) / (temp_lut(t2) - temp_lut(t1))
        tw2 = (tavg - temp_lut(t1)) / (temp_lut(t2) - temp_lut(t1))
      endif

      r1 = 1
      do n=1,nrmix_lut
        if (rmix_lut(n) .gt. ravg) then
          r2 = n
          exit
        else
          r1 = n
          r2 = n
        endif
      enddo
      if (r1 .eq. r2) then
        rw1 = 0.5
        rw2 = 0.5
      else
        rw1 = (rmix_lut(r2) - ravg) / (rmix_lut(r2) - rmix_lut(r1))
        rw2 = (ravg - rmix_lut(r1)) / (rmix_lut(r2) - rmix_lut(r1))
      endif



      freq_index = ifreq_lut(ichan)
      kext = (rw1 * tw1 * pw1 * kabs_lut(p1,t1,r1,freq_index)) + &
             (rw1 * tw1 * pw2 * kabs_lut(p2,t1,r1,freq_index)) + &
             (rw1 * tw2 * pw1 * kabs_lut(p1,t2,r1,freq_index)) + &
             (rw1 * tw2 * pw2 * kabs_lut(p2,t2,r1,freq_index)) + &
             (rw2 * tw1 * pw1 * kabs_lut(p1,t1,r2,freq_index)) + &
             (rw2 * tw1 * pw2 * kabs_lut(p2,t1,r2,freq_index)) + &
             (rw2 * tw2 * pw1 * kabs_lut(p1,t2,r2,freq_index)) + &
             (rw2 * tw2 * pw2 * kabs_lut(p2,t2,r2,freq_index))




      return
      end subroutine monortm_lut

!------------------------------------------------------------------------

      
      END MODULE MonoRTM
