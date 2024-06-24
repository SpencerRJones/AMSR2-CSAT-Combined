


MODULE dsd

USE define_csu1dvar


IMPLICIT NONE


!======================================================================

CONTAINS

!--------------------------------------------------------

!--------------------------------------------------
!
!    Gamma lut consists of an array of
!    values for mu, lambda, and the
!    total water content for each value pair.
!    More information can be found in
!    lwc_gamma_LUT.
!
!    Spencer Jones, CSU Atmos., 10/2023
!
!--------------------------------------------------

SUBROUTINE read_gamma_lut

CHARACTER(LEN=100) :: gamma_lut_file = 'LUT/gamma/lwc_gamma_LUT_N0110000.tbl'

INTEGER :: nN0, nmu, nlamb

OPEN(UNIT=99, FILE=gamma_lut_file, STATUS='OLD', ACCESS='STREAM', &
        FORM='UNFORMATTED')

READ(99) N0_lut
READ(99) nmu
ALLOCATE(mu_lut(nmu))
READ(99) mu_lut(:)
READ(99) nlamb
ALLOCATE(lamb_lut(nlamb))
READ(99) lamb_lut(:)
ALLOCATE(gamma_lwc_lut(nmu,nlamb))
READ(99) gamma_lwc_lut(:,:)


CLOSE(99)


RETURN

END SUBROUTINE read_gamma_lut


!--------------------------------------------------------

SUBROUTINE gamma_dsd(lwc, N_0, mu, lamb)

!Inputs
REAL :: lwc ![g m^-3]       
REAL :: N_0 ![number mm^-1 m^-3]
REAL :: mu  ![]

!Output
REAL :: lamb ![mm^-1]

REAL    :: mu_min
INTEGER :: mu_indx(1)
INTEGER :: lwc_indx(1)
REAL    :: dmu
INTEGER :: N_0_indx(1)
REAL    :: log_n0


mu_indx  = MINLOC(ABS(mu - mu_lut(:)))

lwc_indx = MINLOC(ABS(lwc - gamma_lwc_lut(mu_indx(1),:)))


lamb = lamb_lut(lwc_indx(1))


RETURN

END SUBROUTINE gamma_dsd


!----------------------------------------------------------

!-------------------------------------------------
!
!    Calculates precip rates from RWC and SWC
!
!    -Spencer Jones, CSU Atmos., 04/2024
!
!-------------------------------------------------



SUBROUTINE precip_rate(x, n0_r, n0_s, rain_rate, snow_rate)

REAL    :: x(nvar)
REAL    :: rwc
REAL    :: swc
REAL    :: rain_rate(nz)
REAL    :: snow_rate(nz)
REAL    :: n0_r, n0_s
REAL    :: mu, rho_s
REAL    :: lamb_r, lamb_s
INTEGER :: ilyr, id
REAL    :: diam
REAL    :: diam_increment
REAL    :: num
REAL    :: vol
REAL    :: mass_s
REAL    :: v
REAL    :: c_d
REAL    :: rho_l
REAL    :: rho_a
REAL    :: rr
REAL    :: sr
REAL, PARAMETER :: g = 9.80665
REAL, PARAMETER :: pi = 3.141592654

diam_increment = 0.05 !mm
c_d = 0.5
rho_l = 1000.
rho_a = 1.225 !Air density is constant. Not too bad of an assumption for precip hopefully?



!---Rain rate:

DO ilyr = 1, nz
        

        rwc = x(ilyr)
        if (rwc .eq. 0.) then
                rain_rate(ilyr) = 0.
                cycle
        endif
        mu  = x(nvar-4)
        CALL gamma_dsd(rwc, n0_r, mu, lamb_r)

        rr = 0.

        !Integrate over the DSD:
        DO id = 12, 200
                diam = 0.05 + (diam_increment * float(id))                         !mm
                num  = n0_r * (diam**mu) * exp(-lamb_r * diam)                     !mm^-1 m^-3
                vol  = (pi/6.) * (diam**3) * 1.0E-09                               !m^3   
                !Calculate terminal velocity of raindrops:
                !From Villermaux and Eloi 2011.          
                v    = SQRT((4./3.) * (1./c_d) * (rho_l/rho_a) * g * (diam*0.001)) !m s^-1
                rr   = rr + (num * vol * v * diam_increment * 3.6E+06)             !mm hr^-1
        ENDDO

        rain_rate(ilyr) = rr

ENDDO

!---Snow rate:

DO ilyr = 1, nz

        swc = x(ilyr+nz)
        IF (swc .eq. 0.) THEN
                snow_rate(ilyr) = 0.
                CYCLE
        ENDIF
        rho_s = x(nvar-4) * 1000.  !kg m^-3
        
        lamb_s = (n0_s*1000.*pi*rho_s/(swc*(1.e-3)))**(0.25) !m^-1

        sr = 0.

        DO id = 0, 200
                diam = 0.05 + (diam_increment * float(id))                        !mm
                num  = n0_s * exp(-lamb_s * diam * 0.001)                         !mm^-1 m^-3
                mass_s = rho_s  * (pi/6.) * ((diam*0.001)**3)                     !kg 
                !Terminal velocity of snow:
                !From GPM-DPR algorithm
                v    = 8.8 * SQRT(0.1*diam*(rho_s - rho_a)*0.001)                 !m s^-1
                sr   = sr + (num * (mass_s/rho_l) * v * diam_increment * 3.6E+06) !mm hr^-1
        ENDDO

        snow_rate(ilyr) = sr
         

ENDDO


RETURN

END SUBROUTINE precip_rate

!----------------------------------------------------------


!=========================================================================

END MODULE dsd
