      module mie
      
      USE define_csu1dvar
      !USE cloudice_lut

      implicit none
      
      contains
      
      subroutine mie_clw(freqy, temp, lwc, ksca, asca, gsca, pbck)   

c     Compute the extinction, absorption, asymmetry parameter and
c     backscatter for a given water content of cloud water in [g/m^3]. 
c     Code assumes that cloud water drops are mono-disperse with an 
c     effective radius reff supplied in the parameter file.

c     Input:
c     freqy		frequency of radiation [GHz]
c     temp		temperature of particles [K]
c     lwc		water content of cloud water distribution [g/m**3]

c     Output:
c     ksca		extinction coefficient [1/km]
c     asca		single-scatter albedo []
c     gsca		asymmetry factor []
c     pbck		backscatter phase function/(4*pi) []
c

      implicit none

      real    freqy, temp, lwc, reff_cloudwater
      real    ksca, asca, gsca, pbck

      real    pi, wavel, x
      real    densliq, density
      real    rad, dropmass, vol_lw, vol_drop
      real    num
      real(8)    qext, qsca, asym, qbsca
      real    bext, bsca, bsym, bq11
      
      real     epsreal, epsimag
      complex  ewat!, cref
      complex(8)  cref

c      
c     Assign some useful constants
      data       pi /3.14159265/
      data       densliq /1.0e+3/                 ! kg/m^3
      !data       reff_cloudwater  / 0.0122 /      ! [12.2 um expressed in mm]
      data       reff_cloudwater / 0.01 /  ![mm]

      wavel = 300./freqy  !wavelength in mm for size parameter calculation below

c
c     Begin by checking if hydrometeors of this species are present.
c     If not, set scattering parameters to zero and return.
c
      if (lwc .lt. 0.00001) then
        ksca=0.0
        asca=0.0
        gsca=0.0
        pbck=0.0
        return
      endif
c
c     Initialize the scattering parameters
c
      bext=0.0
      bsca=0.0
      bsym=0.0
      bq11=0.0
c
c     Compute scattering properties
c




      rad = reff_cloudwater
      density = densliq
!      dropmass = 4./3.*pi*density*rad*rad*rad*1.0E-09
!      num = lwc/dropmass  !Is this a bug? Should it be lwc*0.001?

c     Compute number of scatterers from monodisperse assumption

      vol_lw   = (lwc*0.001) / density
      vol_drop = (4./3.)*pi*(rad*0.001)**3

      num = vol_lw / vol_drop




c
c     Get complex refractive index of liquid water
c
      call watoptic(freqy,temp,0.0,epsreal,epsimag)
      ewat = cmplx(epsreal,epsimag)
      cref = csqrt(ewat)

c
c     call Mie program
c
      x = 2.*pi*rad/wavel   !rad in mm and wavel in mm


      call mie_sphere(dble(x),cref,qsca,qext,asym,qbsca)

      bext=num*qext*pi*rad*rad*1.e-6 
      bsca=num*qsca*pi*rad*rad*1.e-6
!      bsym=num*qsca*asym*pi*rad*1.e-6
      bsym = num * qsca * asym * pi * rad * rad * 1.E-06
      bq11=num*qbsca*pi*rad*rad*1.e-6

c
c     check for distribution with very small extinction;
c     set parameters to zero to avoid numerical problems
c
!      if( bext .gt. 1.e-10) then
        ksca=bext
        asca=bsca/bext
        gsca=bsym/bsca
        pbck=bq11/bsca
!      else
!              write(*,*) 'set to 0.', bext
!        ksca=0.0
!        asca=0.0
!        gsca=0.0
!        pbck=0.0
!      end if

      return
      end subroutine mie_clw
      

!-------------------------------------------------------------------------------

      
      subroutine mie_rain(freqy,temp,n0,lam,mu,ksca,asca,gsca,pbck)
c      
c     Compute the extinction, absorption, asymmetry parameter and
c     backscatter for a given rain DSD with modified gamma parameters
c     n0, lam, and mu. Depth of layer (vres) defined in definition file

c     Input:
c     freqy		frequency of radiation [GHz]
c     temp		temperature of particles [K]
c     n0		Gamma DSD scale parameter [mm^(-1-mu) m^(-3)]
c     lam		Gamma DSD slope parameter [mm^-1]
c     mu		Gamma DSD shape parameter []

c     Output:
c     ksca		extinction coefficient [1/km]
c     asca		single-scatter albedo []
c     gsca		asymmetry factor []
c     pbck		backscatter phase function/(4*pi) []
c
      implicit none
c
      ! Input variables      
      real      freqy, temp, n0, lam, mu, plwc
      ! Output variables
      real      ksca, asca, gsca, pbck
      ! Specified parameters
      real      densliq
      ! Internal variables
      integer   i
      real      pi, wavel, diam, diam_increment, x
      real      density, num
      real(8)      qext, qsca, asym, qbsca
      real      bext, bsca, bsym, bq11      
      real      epsreal, epsimag      
      complex   ewat!, cref
      complex(8) cref
c
      data pi /3.14159265/
      data densliq /1.0e+3/        ! kg/m^3
      diam_increment = 0.05        ! mm
c
      wavel = 300./freqy	!mm
      
      ! Initialize scattering parameters
      bext = 0
      bsca = 0
      bsym = 0
      bq11 = 0

      
      ! Loop over particle sizes
      do i = 12, 200
        diam = 0.05 + diam_increment*float(i)
	x = pi*diam/wavel
	num = n0 * diam**mu * exp(-lam*diam) 
	
	! Complex refractive index of liquid water
	call watoptic(freqy,temp,0.0,epsreal,epsimag)
	ewat = cmplx(epsreal,epsimag)
	cref = csqrt(ewat)
	
	! call Mie program
	call mie_sphere(dble(x),cref,qsca,qext,asym,qbsca)
	
	! Integrate over particle size distribution
	bext=bext+num*qext*pi*0.25*diam*diam*diam_increment*1.e-3	! 1e-6 m^2/mm^2 * 1e3 m/km = 1e-3; final units [km^-1]
	bsca=bsca+num*qsca*pi*0.25*diam*diam*diam_increment*1.e-3
        bsym=bsym+num*qsca*asym*pi*0.25*diam*diam*diam_increment*1.e-3
        bq11=bq11+num*qbsca*pi*0.25*diam*diam*diam_increment*1.e-3
      end do
      
c     check for distribution with very small extinction;
c     set parameters to zero to avoid numerical problems
      if( bext .gt. 1.e-10) then
        ksca=bext
        asca=bsca/bext
        gsca=bsym/bsca
        pbck=bq11/bsca
      else
        ksca=0.
        asca=0.
        gsca=0.
        pbck=0.
      end if      
      
      return
      end subroutine mie_rain
      
!------------------------------------------------------------------


       subroutine mie_snow(freqy, temp, swc, n0_sn, rho_sn,
     >                     ksca, asca, gsca, pbck)
c      
c     Compute the extinction, absorption, asymmetry parameter and
c     backscatter for a given water content of snow in [g/m^3], and
c     a particle size distribution n(D) with intercept n0s.

c     Input:
c     freqy             frequency of radiation [GHz]
c     temp              temperature of particles [K]
c     iwc               ice water content of snow distribution [g/m**3]
c     rho_sn            snow particle density [g/cm^3]  !new -SJ

c     Output:
c     ksca              extinction coefficient [1/km]
c     asca              single-scatter albedo []
c     gsca              asymmetry factor []
c     pbck              backscatter phase function/(4*pi) []
c

      implicit none

      integer  i
      real     freqy, temp, swc
      REAL     n0_sn, rho_sn
      real     ksca, asca, gsca, pbck
      real     wavel, pi, densice, diam, diam_increment, density, fincl
      real     n0_snow, density_snow, lam, num, x
      real(8)     qext, qsca, asym, qbsca
      real     bext, bsca, bsym, bq11
      real     epsreal, epsimag
      complex  eice, eair, emg!, cref
      complex(8) cref



      data     pi /3.14159265/
      data     densice /0.917e+3/               ! [kg/m^3]  
      !data     n0_snow / 1.6e+7 /               ! [1/m^4]
      !data     density_snow / 0.1E+3 /          ! [kg/m^3]  
      data     diam_increment / 0.10 /          ! mm  

      n0_snow      = n0_sn * 1000.   ![mm^-1 m^-3] ---> [m^-4]
      density_snow = rho_sn * 1000.  ![g cm^-3]    ---> [kg m^-3]
!!!
        !rho_sn = 0.35
        !swc    = 10**(-0.1)
        !temp   = 262.2

!!!

c      
c     Assign some useful constants
      wavel = 300./freqy
c
c     Begin by checking if hydrometeors of this species are present.
c     If not, set scattering parameters to zero and return.
c
      !if(swc .lt. 0.0025) then
      if (swc .lt. 1.0e-06) then
        ksca=0.
        asca=0.
        gsca=0.
        pbck=0.
        return
      endif


c
c     If hydrometeors are present, initialize the scattering parameters
c
      bext=0.
      bsca=0.
      bsym=0.
      bq11=0.



c
c     Loop over particle sizes:

c     increments of particle diameter are 0.10 mm; the particle
c     size distribution is expressed as a particle number density,
c     num, per diameter increment; the original psd is
c     n(D) = n0s * exp(-lam * D) where n is the number density
c     per diameter increment, n0s is the distribution intercept,
c     lam is the slope of the distribution of ln(n(D)), and D is
c     the particle diameter. It is assumed here that n0_snow 
c     and swc are prescribed, and lam is therefore constrained 
c     to yield the prescribed water content:
c     lwc = integral {n(D) * pi * (D**3) * density(D) * dD/6}
c     therefore:
c     lam=(n0s*pi*density_snow/lwc)**(0.25)

c


        call iceoptic(freqy,temp,epsreal,epsimag)
        eice = cmplx(epsreal,epsimag)
        eair = cmplx(1.0006,0.0)

        fincl = 1. - (density_snow/densice)
        call mg_ellips(fincl, eice, eair, emg)
        cref = csqrt(emg)


      do i=0,200
        diam = 0.050 + diam_increment*float(i)
        x =  pi*diam/wavel
        lam = (N0_snow*pi*density_snow/(swc*(1.e-3)))**(0.25)
        num = N0_snow*exp(-lam*diam*(1.e-3))

c       Below code moved to outside loop for speed
c
c       complex refractive index of snow
c
c        call iceoptic(freqy,temp,epsreal,epsimag)
c        eice = cmplx(epsreal,epsimag)
c        eair = cmplx(1.0006,0.0)


c       calculate dielectric constant of snow as an ice matrix  
c       with air inclusions, using Maxwell-Garnett mixing for  
c       2-component media w. elliptical inclusions 
c        fincl = 1. - (density_snow/densice)
c        call mg_ellips(fincl, eice, eair, emg)
c        cref = csqrt(emg)
c
c       call Mie program
c



        call mie_sphere(dble(x),cref,qsca,qext,asym,qbsca)

        qext = qext
        qsca = qsca
        asym = asym
        qbsca = qbsca



c
c       integrate over particle size distribution;

        bext=bext+num*qext*pi*0.25*diam*diam*diam_increment*1.e-6 !units = [km^-1]
        bsca=bsca+num*qsca*pi*0.25*diam*diam*diam_increment*1.e-6
        bsym=bsym+num*qsca*asym*pi*0.25*diam*diam*diam_increment*1.e-6
        bq11=bq11+num*qbsca*pi*0.25*diam*diam*diam_increment*1.e-6

      end do

c
c     check for distribution with very small extinction;
c     set parameters to zero to avoid numerical problems
      if( bext .gt. 1.e-10) then

        ksca=bext
        asca=bsca/bext
        gsca=bsym/bsca
        pbck=bq11/bsca

      else
              !write(*,*) 'set to 0.'

        ksca=0.
        asca=0.
        gsca=0.
        pbck=0.

      end if

      return
      end subroutine mie_snow


!----------------------------------------------------------------

c
c     Lookup table version of mie_snow
c     Spencer Jones CSU, Atmos.,  03/2024
c


      subroutine mie_ice_lut(ifreq, temp, swc, n0_sn, rho_sn,
     >                       ksca, asca, gsca, pbck)
        implicit none

      integer  ifreq
      real     temp, swc
      REAL     n0_sn, rho_sn
      real     ksca, asca, gsca, pbck
      real     ksca1, ksca2, asca1, asca2 
      real     gsca1, gsca2, pbck1, pbck2
      real     T1, T2, dT, iwc1, iwc2, diwc, rho1, rho2, drho
      integer  freq_indx(1), temp_indx, iwc_indx, rho_indx

        !rho_sn = 0.35
        !swc    = 10**(-0.1)
        !temp   = 262.2


c      Check limits of LUT:
      IF (swc .GT. 10.) swc = 9.99





c      Get indices:


      
      freq_indx = ifreq_lut(ifreq)
      temp_indx = int((dble(temp) - temp_icelut(1)) /
     >            stepsize_T) + 1
      iwc_indx  = int((log10(dble(swc)) - iwc_icelut(1)) / 
     >            stepsize_iwc) + 1
      rho_indx  = int((dble(rho_sn) - rho_icelut(1)) / 
     >            stepsize_rho) + 1

      !if (temp_indx .gt. size(temp_icelut)) then
      !       write(*,*) 'reached temp limit in icelut'
      !      stop
      !endif 
      !if (iwc_indx .gt. size(iwc_icelut)) then
      !              write(*,*) swc
      !       write(*,*) 'reached iwc limit in icelut'
      !      stop
      !endif
      !      if (rho_indx .gt. size(rho_icelut)) then
      !       write(*,*) 'reached rho limit in icelut'
      !      stop
      !endif

      if (temp_indx .lt. 1) temp_indx = 1
      if (temp_indx .ge. size(temp_icelut))
     >       temp_indx = size(temp_icelut) - 1
      if (iwc_indx  .lt. 1) iwc_indx = 1
      if (iwc_indx  .ge. size(iwc_icelut))
     >       iwc_indx = size(iwc_icelut) - 1
      if (rho_indx  .lt. 1) rho_indx = 1
      if (rho_indx  .ge. size(rho_icelut)) 
     >       rho_indx = size(rho_icelut) - 1

      if (temp .lt. temp_icelut(temp_indx)) temp_indx = temp_indx + 1
      if (temp .gt. temp_icelut(temp_indx+1)) temp_indx = temp_indx - 1
      if (swc  .lt. 10**(iwc_icelut(iwc_indx))) iwc_indx = iwc_indx - 1
      if (swc  .gt. 10**(iwc_icelut(iwc_indx+1))) iwc_indx = iwc_indx+1
      if (rho_sn .lt. rho_icelut(rho_indx)) rho_indx = rho_indx + 1
      if (rho_sn .gt. rho_icelut(rho_indx+1)) rho_indx = rho_indx - 1
      


      if (temp .lt. temp_icelut(temp_indx) .or. 
     >    temp .gt. temp_icelut(temp_indx+1)) then
             write(*,*) 'indexing problem:'
             write(*,*) 'temp: ', temp
             write(*,*) 'axis value above: ', temp_icelut(temp_indx+1)
             write(*,*) 'axis value below: ', temp_icelut(temp_indx) 
             stop
      endif

      if (swc .lt. 10**(iwc_icelut(iwc_indx)) .or.
     >    swc .gt. 10**(iwc_icelut(iwc_indx+1))) then
           write(*,*) 'indexing problem:'
           write(*,*) 'swc: ', swc
           write(*,*) 'axis value above: ', 10**(iwc_icelut(iwc_indx+1))
           write(*,*) 'axis value below: ', 10**(iwc_icelut(iwc_indx))
           stop
      endif

      if (rho_sn .lt. rho_icelut(rho_indx) .or.
     >    rho_sn .gt. rho_icelut(rho_indx+1)) then
             write(*,*) 'indexing problem:'
             write(*,*) 'rho: ', rho_sn
             write(*,*) 'axis value above: ', rho_icelut(rho_indx+1)
             write(*,*) 'axis value below: ', rho_icelut(rho_indx)
             stop
      endif



!      write(*,*) temp_indx, size(temp_icelut), 
!     >           temp, temp_icelut(temp_indx),temp_icelut(temp_indx + 1)
!      write(*,*) iwc_indx,  size(iwc_icelut),
!     >           log10(swc),iwc_icelut(iwc_indx),iwc_icelut(iwc_indx+1)
!      write(*,*) rho_indx, size(rho_icelut),
!     >           rho_sn, rho_icelut(rho_indx),rho_icelut(rho_indx+1)




        call interpolate_trilinear(temp, log10(swc), rho_sn, 
     >              temp_icelut(temp_indx), temp_icelut(temp_indx+1), 
     >              iwc_icelut(iwc_indx), iwc_icelut(iwc_indx+1),
     >              rho_icelut(rho_indx), rho_icelut(rho_indx+1),
     >              ksca_icelut(freq_indx,
     >                          temp_indx:temp_indx+1,
     >                          iwc_indx:iwc_indx+1,
     >                          rho_indx:rho_indx+1),
     >                          ksca)

        call interpolate_trilinear(temp, log10(swc), rho_sn,
     >              temp_icelut(temp_indx), temp_icelut(temp_indx+1),
     >              iwc_icelut(iwc_indx), iwc_icelut(iwc_indx+1),
     >              rho_icelut(rho_indx), rho_icelut(rho_indx+1),
     >              asca_icelut(freq_indx,
     >                          temp_indx:temp_indx+1,
     >                          iwc_indx:iwc_indx+1,
     >                          rho_indx:rho_indx+1),
     >                          asca)


        call interpolate_trilinear(temp, log10(swc), rho_sn,
     >              temp_icelut(temp_indx), temp_icelut(temp_indx+1),
     >              iwc_icelut(iwc_indx), iwc_icelut(iwc_indx+1),
     >              rho_icelut(rho_indx), rho_icelut(rho_indx+1),
     >              gsca_icelut(freq_indx,
     >                          temp_indx:temp_indx+1,
     >                          iwc_indx:iwc_indx+1,
     >                          rho_indx:rho_indx+1),
     >                          gsca)


         call interpolate_trilinear(temp, log10(swc), rho_sn,
     >              temp_icelut(temp_indx), temp_icelut(temp_indx+1),
     >              iwc_icelut(iwc_indx), iwc_icelut(iwc_indx+1),
     >              rho_icelut(rho_indx), rho_icelut(rho_indx+1),
     >              pbck_icelut(freq_indx,
     >                          temp_indx:temp_indx+1,
     >                          iwc_indx:iwc_indx+1,
     >                          rho_indx:rho_indx+1),
     >                          pbck)



       
   
   
         return

      end subroutine mie_ice_lut

!-----------------------------------------------------------------      

C
C     Trilinear interpolation routine
C
C     Interpolates to an arbitrary value in linear space
C        for any 3-dimensional lookup table. Takes point
C        in 3-D space and 2x2x2 cube containing point.
C
C     Inputs:
C        x:       location on x-axis of point of interest
C        y:       location on y-axis
C        z:       location on z-axis
C        x0:      lower x value from lookup table
C        y0:      lower y value
C        z0:      lower z value
C        x1:      upper x value from lookup table
C        y1:      upper y value
C        z1:      upper z value
C        p_table: 2x2x2 cube containing point of interest.
C             corners of cube p_table are at:
C                (x0,y0,z0), (x1,y0,z0), (x0,y1,z0), (x1,y1,z0),
C                (x0,y0,z1), (x1,y0,z1), (x0,y1,z1), (x1,y1,z1)
C
C    Output:
C        p:        value of point p
C
C
C    Spencer Jones, CSU Atmos., --original coder
C


      SUBROUTINE interpolate_trilinear(x, y, z, 
     >             x0, x1, y0, y1, z0, z1, p_table, p)

        REAL :: x, y, z
        REAL :: x0, x1, y0, y1, z0, z1
        REAL :: dx, dy, dz
        REAL :: c0, c1, c2, c3, c4, c5, c6, c7
        REAL :: p000, p100, p010, p110, p001, p101, p011, p111
        REAL :: p_table(2,2,2)
        REAL :: p

        dx = (x - x0) / (x1 - x0)
        dy = (y - y0) / (y1 - y0)
        dz = (z - z0) / (z1 - z0)

        p000 = p_table(1,1,1)
        p100 = p_table(2,1,1)
        p010 = p_table(1,2,1)
        p110 = p_table(2,2,1)
        p001 = p_table(1,1,2)
        p101 = p_table(2,1,2)
        p011 = p_table(1,2,2)
        p111 = p_table(2,2,2)
        
        c0 = p000
        c1 = p100 - p000
        c2 = p010 - p000
        c3 = p001 - p000
        c4 = p110 - p010 - p100 + p000
        c5 = p011 - p001 - p010 + p000
        c6 = p101 - p001 - p100 + p000
        c7 = p111 - p011 - p101 - p110 + p100 + p001 + p010 - p000

        p = c0 + (c1*dx) + (c2*dy) + (c3*dz) +
     >      (c4*dx*dy) + (c5*dy*dz) + (c6*dz*dx) + (c7*dx*dy*dz)


        RETURN
        END SUBROUTINE interpolate_trilinear

!----------------------------------------------------------------

      SUBROUTINE MIE_SPHERE (X, MIN, QSCAT, QEXTI, ASYM, QBSCAT)
**
      implicit    none
      SAVE
**
**    Mie Routine P. Bauer 
**
      integer    limitx
      PARAMETER (LIMITX = 1500)

**
      REAL(8)        X
      REAL(8)        MR, MI, N1, N2
      REAL(8)        QSCAT, QEXTI, QABSO, ASYM, QBSCAT
**
      REAL(8)        RFAC1, RFAC2
      REAL(8)        RHELP1(2), RHELP2(2)
**
      COMPLEX(8)     M, MX, MIN
      COMPLEX(8)     CHELP1, CHELP2, CFAC1, CFAC2, CBSCAT
**
      COMPLEX(8)     DN(0:LIMITX), WN(-1:LIMITX)
      COMPLEX(8)     AN(LIMITX), BN(LIMITX)
**
      INTEGER     NEND
      INTEGER     I100, I101
**
      EQUIVALENCE (CHELP1, RHELP1 (1))
      EQUIVALENCE (CHELP2, RHELP2 (1))
**
************************************************************************
**
      M      = CONJG (MIN)
      CHELP1 = M
      MR     =        RHELP1 (1)
      MI     = -1.0 * RHELP1 (2)
**      
      MX   = M  * X
      N1   = MR * X
      N2   = MI * X
**
      IF (X .LE. 20000.0) NEND = X + 4.00 * X ** (1.0 / 3.0) + 2.0
      IF (X .LE.  4200.0) NEND = X + 4.05 * X ** (1.0 / 3.0) + 2.0
      IF (X .LE.     8.0) NEND = X + 4.00 * X ** (1.0 / 3.0) + 1.0
      IF (NEND .LE.    5) NEND = 5
      IF (NEND .GT. LIMITX) NEND = LIMITX
**
      RFAC1      = SIN  (N1) * SIN  (N1) + SINH (N2) * SINH (N2)
      RHELP1 (1) = SIN  (N1) * COS  (N1) / RFAC1
      RHELP1 (2) = SINH (N2) * COSH (N2) / RFAC1
**
      DN (0) = CHELP1
**
      RHELP1 (1) =             COS (X)
      RHELP1 (2) = -1.0 E+00 * SIN (X)
      RHELP2 (1) =             SIN (X)
      RHELP2 (2) =             COS (X)
**
      WN (-1) = CHELP1
      WN ( 0) = CHELP2
**
      QEXTI  = 0.0
      QSCAT  = 0.0
      QBSCAT = 0.0
      QABSO  = 0.0
      ASYM   = 0.0 
      CBSCAT = CMPLX (0.0,0.0)
**
      DO 100 I100 = 1, NEND
         DN (I100) = -1.0 * I100 / MX
     +             +  1.0 / (I100 / MX - DN (I100 - 1))
         WN (I100) = WN (I100 - 1) * (2.0 * I100 - 1.0) / X
     +             - WN (I100 - 2)
**
         CFAC1 = DN (I100) / M + I100 / X
         CFAC2 = M * DN (I100) + I100 / X
**
         CHELP1 = WN (I100)
         CHELP2 = WN (I100 - 1)
**
         AN (I100) = (CFAC1 * RHELP1 (1) - RHELP2 (1))
     +             / (CFAC1 * CHELP1     - CHELP2    )
         BN (I100) = (CFAC2 * RHELP1 (1) - RHELP2 (1))
     +             / (CFAC2 * CHELP1     - CHELP2    )
**
         CHELP1 = AN (I100)
         CHELP2 = BN (I100)
**
         RFAC1 = RHELP1 (1) + RHELP2 (1)
         RFAC2 = ABS (AN (I100)) * ABS (AN (I100))
     +         + ABS (BN (I100)) * ABS (BN (I100))
**
         QEXTI  = QEXTI  + (2.0 * I100 + 1.0) * RFAC1
         QSCAT  = QSCAT  + (2.0 * I100 + 1.0) * RFAC2
         CBSCAT = CBSCAT + (2.0 * I100 + 1.0) * (-1.0) ** I100
     +          * (AN (I100) - BN (I100))
**
         IF (I100 .EQ. 1) GO TO 100
**
         CHELP1 = AN (I100 - 1) * CONJG (AN (I100))
     +          + BN (I100 - 1) * CONJG (BN (I100))
         CHELP2 = AN (I100 - 1) * CONJG (BN (I100 - 1))
**
         I101 = I100 - 1
         RFAC1  = I101 * (I101 + 2) / (I101 + 1.0)
         RFAC2  = (2.0 * I101 + 1.0) / (I101 * (I101 + 1.0))
**
         ASYM = ASYM + RFAC1 * RHELP1 (1) + RFAC2 * RHELP2 (1)
100   CONTINUE
**
      QEXTI  = QEXTI * 2.0 / (X * X)
      QSCAT  = QSCAT * 2.0 / (X * X)
      ASYM   = ASYM  * 4.0 / (X * X * QSCAT)
      QBSCAT = ABS (CBSCAT) * ABS (CBSCAT) / (X * X)
      IF (QSCAT .GT. QEXTI) QSCAT = QEXTI
**
      RETURN
      END SUBROUTINE MIE_SPHERE
      
      end module mie
       
