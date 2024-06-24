      subroutine forward_model(hgt, plevel, tlevel, mix_ratio, cloud_water, rain_water, ice_water, tsfc, saln, wspd, &
                               N_0_plw, mu_plw, N_0_ice, rho_ice, &
                               tbout, z_eff, z_app, pia)

      USE define_csu1dvar
      USE MonoRTM
      USE FASTEM
      USE dsd
      USE mie

!     Compute brightness temperatures for absorbing atmosphere.
!     Subroutines called:
!
!                    **  RADTRAN I/O SPECIFICATIONS  **


!Version 2.1

!-----------------------------------------------------
!
!   Computes brightness temperatures for absorbing
!       and scattering atmosphere. Scattering in 
!       radtran added by Rick Schulte.
!   Modified by Spencer Jones to include gaseous
!       absorption and scattering by ice and 
!       liquid particles according to PSDs for
!       AMSR2 channels and CPR operating frequency.
!
!---------------------------------------------------

      implicit none

      integer, parameter :: diffuse = 1
      integer, parameter :: offsets = 0
      integer, parameter :: verbose = 0

        REAL, PARAMETER :: pi = 3.141592654

      real    :: hgt(nz+1)
      real    :: plevel(nz+1)
      real    :: tlevel(nz+1)
      real    :: mix_ratio(nz)
      real    :: cloud_water(nz)
      REAL    :: rain_water(nz)
      real    :: ice_water(nz)
      real    :: tsfc
      real    :: saln
      real    :: wspd
      REAL    :: N_0_plw
      REAL    :: mu_plw
      REAL    :: N_0_ice
      REAL    :: rho_ice
      real    :: wdir

      real    :: tbout(nch)
      REAL    :: z
      REAL    :: lamb_plw

      real    :: rhgt(0:nz)
      real    :: rtmp(0:nz)
      real    :: rprs(0:nz)
      real    :: rmix(nz)
      !real    :: clwp(nz)
      !real    :: clwc(nz)
      !real    :: ciwp(nz)
      !real    :: ciwc(nz)
      REAL    :: rclwc(nz)
      REAL    :: rplwc(nz)
      REAL    :: rtiwc(nz)
      !REAL    :: rciwc(nz)
      real    :: pavg
      real    :: tavg
      real    :: dhgt
      real    :: view_angle
      real    :: azim_angle
      real    :: emis(nch)
      real    :: refl(nch)
      real    :: salin
      real    :: freq
      integer :: pol
      real    :: femis(4), frefl(4)
      real    :: temis, ebar
      real    :: kext(nz), salb(nz), asym(nz)
      real    :: atm_ext, clw_ext
      real    :: kext_clw, salb_clw, asym_clw, pback_clw
      REAL    :: kext_plw, salb_plw, asym_plw, pback_plw
      REAL    :: kext_ice, salb_ice, asym_ice, pback_ice
      real    :: tb, tbdown
      real    :: omega(2)
      real    :: e0(2)
      real    :: ewind(2)
      real    :: edirstokes(4)
      real    :: trans, optdepth
      integer :: ilyr, ichan, n2, ifrq
      integer :: chan_index(10) = (/1,2,3,4,5,6,7,8,9,10/)
      ! These variables added for handling mixed polarization of cross track scanners (Rick, Feb 2018)
      real    :: delta, alpha
      real    :: vemis, hemis
      real    :: vrefl, hrefl
      ! These variables added for addition of cloud ice (Rick, Feb 2018)
      !real    :: kext_ice, salb_ice, asym_ice


        REAL :: z_eff(nz,nact_freqs)
        REAL :: z_app(nz,nact_freqs)
        REAL :: pia(nact_freqs)

        INTEGER :: nskipped 

      REAL :: k_dial
      REAL :: tau, tau_lyr
      REAL :: backsca(nz)
      REAL :: lamb
      REAL :: transmittance


      ! Loop over each channel
      azim_angle = 0.0
      salin      = saln

!        write(*,*) 'Into compute_tb:'


      if (first) then
        write(6,'("hgt        = ",30(f7.2))') hgt(1:nz+1)
        write(6,'("plevel     = ",30(f7.2))') plevel(1:nz+1)
        write(6,'("tlevel     = ",30(f7.2))') tlevel(1:nz+1)
        write(6,'("rmix       = ",30(f7.2))') mix_ratio(1:nz)
        write(6,'("clwc       = ",30(f7.4))') cloud_water(1:nz)
        !write(6,'("ciwc       = ",30(f7.4))') cloud_ice(1:nz)
        write(6,'("plwc       - ",30(f7.4))') rain_water(1:nz)
        write(6,'("tiwc       - ",30(f7.4))') ice_water(1:nz)
        write(6,'("frzlvl     - ",f7.4)') frzlvl_in(oepix)
        write(6,'("frzl_bin   - ",i7)') frzlvl_bin_in(oepix)
        write(6,'("cb_bin     - ",i7)') clbot_bin(oepix)
        write(6,'("nchan      = ",i7)') nch
        write(6,'("nlyr       = ",i7)') nz
        write(6,'("view_angle = ",30(f7.2))') eia_in(oepix,chan_index(1:nch))
        write(6,'("azim_angle = ",f7.2)') azim_angle
        write(6,'("salin      = ",f7.2)') salin
        write(6,'("tsfc       = ",f7.2)') tsfc
        write(6,'("wspd       = ",f7.2)') wspd
        write(6,'("wdir       = ",f7.2)') wdir
        write(6,'("freq       = ",30(f7.2))') freq_lut
        write(6,'("ifreq      = ",30(i7))') ifreq_lut
        write(6,'("ipol       = ",30(i7))') ipol_lut
        write(6,*)
        first = .FALSE.
      endif

      

      ! reverse hgt coordinate to go from surface to TOA

      rhgt(0) = hgt(nz+1)
      rtmp(0) = tlevel(nz+1)
      rprs(0) = plevel(nz+1)

      do ilyr = 1, nz
        n2 = nz - ilyr + 1
        rhgt(n2) = hgt(ilyr)
        rtmp(n2) = tlevel(ilyr)
        rprs(n2) = plevel(ilyr)
        rmix(n2) = mix_ratio(ilyr)
        rclwc(n2) = cloud_water(ilyr)
        rplwc(n2) = rain_water(ilyr)
        rtiwc(n2) = ice_water(ilyr)
      enddo


        dhgt = 0.5 ![km]

      ! loop though channels to compute tb

kext = 0.
salb = 0.
asym = 0.


DO ichan = 1, nch

        view_angle = eia_in(oepix, chan_index(ichan))
        freq       = freq_lut(ifreq_lut(ichan))
        pol        = ipol_lut(ichan) + 1

        optdepth = 0.


        DO ilyr = 1, nz

        
                tavg = 0.5 * (rtmp(ilyr) + rtmp(ilyr - 1))
                pavg = (rprs(ilyr-1) - rprs(ilyr)) / LOG(rprs(ilyr-1)/rprs(ilyr))


                !Gas absorption
                CALL monortm_lut(ichan, pavg, tavg, rmix(ilyr), atm_ext)
                
                !k_ext for cloud water
                CALL mie_clw(freq, tavg, rclwc(ilyr), kext_clw, salb_clw, asym_clw, pback_clw)
                
                kext_clw = kext_clw * 1000. ![m^-1] ---> [km^-1]


                !write(*,*) kext_clw, salb_clw, asym_clw, pback_clw

                !k_ext for liquid hydrometeors
                IF (rplwc(ilyr) .GT. 1.0E-06) THEN


                        CALL gamma_dsd(rplwc(ilyr), N_0_plw, mu_plw, lamb_plw)
                        
                       ! write(*,*) N_0_plw, mu_plw, lamb_plw


                        CALL mie_rain(freq, tavg, N_0_plw, lamb_plw, mu_plw, &
                                      kext_plw, salb_plw, asym_plw, pback_plw)
                        !write(*,*) kext_plw, salb_plw, asym_plw, pback_plw

                       

                ELSE
                        kext_plw = 0.
                        salb_plw = 0.
                        asym_plw = 0.
                        pback_plw = 0.

                ENDIF

                !k_ext for ice
                IF (rtiwc(ilyr) .GT. 1.0E-06) THEN

                        !Full routine
                        !CALL mie_snow(freq, tavg, rtiwc(ilyr), N_0_ice, rho_ice, &
                        !              kext_ice, salb_ice, asym_ice, pback_ice)
                        !write(*,*) kext_ice, salb_ice, asym_ice, pback_ice

                        !Lookup-table version
                        CALL mie_ice_lut(ichan, tavg, rtiwc(ilyr), N_0_ice, rho_ice, &
                                         kext_ice, salb_ice, asym_ice, pback_ice)

                        !write(*,*) kext_ice, salb_ice, asym_ice, pback_ice
                        !write(*,*) '--------'
                ELSE
                        kext_ice = 0.
                        salb_ice = 0.
                        asym_ice = 0.
                        pback_ice = 0.

                ENDIF


                !Total k_ext for layer is sum of k_ext from gases, cloud water, and hydrometeors
                kext(ilyr) = atm_ext + kext_clw + kext_plw + kext_ice

                !Total single-scatter albedo is a weighted sum of liquid and ice particle salb
                salb(ilyr) = ((salb_clw*kext_clw) + (salb_plw*kext_plw) + (salb_ice*kext_ice)) / kext(ilyr)

                !Total asymmetry parameter is a weighted sum as well
                IF (salb(ilyr) .GT. 0) THEN
                        asym(ilyr) = ((asym_clw*salb_clw*kext_clw) +  &
                                      (asym_plw*salb_plw*kext_plw) +  &
                                      (asym_ice*salb_ice*kext_ice)) / &
                                      (salb(ilyr) * kext(ilyr))
                ELSE
                        asym(ilyr) = 0.
                ENDIF

                !Finally, optical depth is a running total of k_ext for all layers
                optdepth = optdepth + (kext(ilyr) * dhgt)

        ENDDO !layers

        trans = exp(-optdepth)

        CALL compute_fastem(6, offsets, diffuse, freq, view_angle, tsfc, salin, &
                            wspd, azim_angle, trans, femis, frefl, verbose=verbose)

                !Select correct emissivity based on polarization

        IF (pol .EQ. 1) THEN     !Vertical polarization
                emis(ichan) = femis(1)
                refl(ichan) = frefl(1)
        ELSEIF (pol .EQ. 2) THEN !Horizontal polarization
                emis(ichan) = femis(2)
                refl(ichan) = frefl(2)
        ELSE !Should not have any other value for AMSR2
                WRITE(*,*) 'Polarization value out of bounds.'
                STOP
        ENDIF


        ebar = emis(ichan)
        
        !WRITE(*,*) 'Into radtran:'
        !WRITE(*,*) 'nz:          ', nz
        !WRITE(*,*) 'view_angle:  ', view_angle
        !WRITE(*,*) 'tsfc:        ', tsfc
        !WRITE(*,*) 'rhgt:        ', rhgt
        !WRITE(*,*) 'rtmp:        ', rtmp
        !WRITE(*,*) 'kext:        ', kext
        !WRITE(*,*) 'salb:        ', salb
        !WRITE(*,*) 'emis(ichan): ', emis(ichan)
        !WRITE(*,*) 'ebar:        ', ebar
        !WRITE(*,*) 'refl(ichan): ', refl(ichan)


        CALL radtran(nz,view_angle,tsfc,rhgt,rtmp,kext,salb,asym,emis(ichan),ebar,refl(ichan),tb,tbdown)

        tbout(ichan) = tb

ENDDO !channels


!write(*,*) 'Tbs: ', tbout


!---Radar frequencies:
!
! Reflectivites must be calculated from TOA to sfc

z_app = -9999.9
z_eff = -9999.9
pia   = -9999.9

k_dial = 0.75

DO ifrq = 1, nact_freqs

        freq = 94.0

        lamb = 300. / freq

        tau = 0.

        DO ilyr = 1, nlyrs

                
                tavg = 0.5 * (tlevel(ilyr) + tlevel(ilyr+1))
                pavg = (plevel(ilyr) - plevel(ilyr+1)) / LOG(plevel(ilyr)/plevel(ilyr+1))


                CALL monortm_lut(nch+ifrq, pavg, tavg, mix_ratio(ilyr), atm_ext)

                !Get cloud water backscatter
                CALL mie_clw(freq, tavg, cloud_water(ilyr), kext_clw, salb_clw, asym_clw, pback_clw)

                kext_clw = kext_clw * 1000.


                !Get liquid hydrometeor backscatter
                IF (rain_water(ilyr) .GT. 1.0E-06) THEN

                        CALL gamma_dsd(rain_water(ilyr), N_0_plw, mu_plw, lamb_plw)

                        CALL mie_rain(freq, tavg, N_0_plw, lamb_plw, mu_plw, &
                                      kext_plw, salb_plw, asym_plw, pback_plw)


                ELSE
                        kext_plw = 0.
                        salb_plw = 0.
                        asym_plw = 0.
                        pback_plw = 0.

                ENDIF


                !Get ice backscatter
                IF (ice_water(ilyr) .GT. 1.0E-06) THEN


                        !write(*,*) freq, tavg, ice_water(ilyr), rho_ice
                        !CALL mie_snow(freq, tavg, ice_water(ilyr), N_0_ice, rho_ice, &
                        !              kext_ice, salb_ice, asym_ice, pback_ice)

                        !write(*,*) kext_ice, salb_ice, asym_ice, pback_ice

                        CALL mie_ice_lut(nch+ifrq, tavg, ice_water(ilyr), N_0_ice, rho_ice, &
                                         kext_ice, salb_ice, asym_ice, pback_ice)

                        
                        !write(*,*) kext_ice, salb_ice, asym_ice, pback_ice
                        !write(*,*) '---------------'




                ELSE

                        kext_ice = 0.
                        salb_ice = 0.
                        asym_ice = 0.
                        pback_ice = 0.

                ENDIF


                kext(ilyr) = atm_ext + kext_clw + kext_plw + kext_ice

                tau_lyr         = kext(ilyr) * dhgt
                tau             = tau + tau_lyr
                !opt_depth(ilyr) = tau

                backsca(ilyr) = (pback_clw * salb_clw * (kext_clw/1000.)) + &
                                (pback_plw * salb_plw * (kext_plw/1000.)) + &
                                (pback_ice * salb_ice * (kext_ice/1000.)) !Backscatter is calculated in m for conversion below


                IF (backsca(ilyr) .GT. 0.0) THEN
                        
                        z_eff(ilyr,ifrq) = (((lamb*0.001)**4) / (pi**5 * K_dial)) &
                                              * backsca(ilyr) * 1.0E+18      !Units have to be mm^6 m^-3


                        transmittance = EXP(-2.*tau)


                        z_app(ilyr,ifrq) = Z_eff(ilyr,ifrq) * transmittance

                        z_eff(ilyr,ifrq) = 10. * LOG10(Z_eff(ilyr,ifrq)) !convert to dBZe
                        z_app(ilyr,ifrq) = 10. * LOG10(Z_app(ilyr,ifrq)) !convert to dBZ

                        IF (Z_eff(ilyr,ifrq) .LT. -30.) Z_eff(ilyr,ifrq) = -30. 
                        IF (Z_app(ilyr,ifrq) .LT. -30.) Z_app(ilyr,ifrq) = -30. !Set bounds on attenuated reflectivities

                ELSE
                        nskipped = nskipped + 1
                        z_eff(ilyr,ifrq) = -30.
                        z_app(ilyr,ifrq) = -30.

                ENDIF


        ENDDO !layers

        pia(ifrq) = 10.*LOG10(2.*tau)        


!write(*,*) 'Z_eff: ', Z_eff
!write(*,*) 'Z_app: ', Z_app
!stop 434

ENDDO !radar frequencies



RETURN


      end subroutine forward_model
