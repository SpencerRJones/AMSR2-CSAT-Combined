!--------------------------------------------------------------------
!
!       CSU1DVAR.f90
!
!       Type: Main Program
!
!       -Originally written by Dave Duncan, 2016
!       -Adapted for TEMPEST-D by Rick Schulte
!       -Code recompiled, simplified, and re-adapted for
!               GMI by Spencer Jones 04/2023
!
!       -Combined AMSR2-CloudSat retrieval finalized 06/20/2024
!               -Spencer Jones, Department of Atmospheric Science, 
!                               Colorado State University
!
!---------------------------------------------------------------------



PROGRAM csu1dvar

USE define_csu1dvar
USE subr_csu1dvar
USE MonoRTM
USE optest
USE dsd

IMPLICIT NONE

!----------------------------------------------------------------



INTEGER            :: nargs
REAL               :: cputime
INTEGER            :: cpu1, cpu2
CHARACTER(LEN=256) :: input_file, output_file


!---Check command line arguments:

nargs = IARGC()

IF (nargs .EQ. 2) THEN
        CALL GETARG(1, input_file)
        CALL GETARG(2, output_file)
ELSE
        WRITE(*,*) 'Wrong number of command line args!'
        STOP
ENDIF




!---Print startup information:


WRITE(*,*) '|----------------------------------------------------------|'
WRITE(*,*) '|                                                          |'
WRITE(*,*) '|                                                          |'
WRITE(*,*) '|   CSU 1DVAR Retrieval                                    |'
WRITE(*,*) '|                                                          |'
WRITE(*,*) '|   Optimal Estimation Algorithm                           |'
WRITE(*,*) '|                                                          |'
WRITE(*,*) '|----------------------------------------------------------|'
WRITE(*,*)
WRITE(*,*) 'Input file:  ', input_file
WRITE(*,*) 'Output file: ', output_file



!---Read in MonoRTM LUT. Subroutine contained in monortm-lut.f90.
!---gamma_lut and ice_lut read routines in dsd.f90

CALL read_monortm_lut
CALL read_gamma_lut
CALL read_ice_lut


!---Read input file:

CALL read_pp(input_file)



!---Prep OE; in subr_csu1dvar.f90. Creates variables for OE code.

CALL prep_oe


!---Begin OE retrieval. OE code found in optest.f90
WRITE(*,*) 'Starting retrieval...'
CALL retrieval_nr
WRITE(*,*) 'Retrieval finished.'


CALL CPU_TIME(cputime)
cputime = cputime / 60.0
cpu1 = INT(cputime)
cpu2 = NINT((cputime - cpu1) * 60)

IF (debug_mode .EQ. .TRUE.) THEN
        WRITE(*,*) 'CPU time (m:s) = ',cpu1,':',cpu2
        WRITE(*,*) ''
ENDIF



!---Write output file.
WRITE(*,*) 'Finished. Creating output file.'

CALL output_nc(output_file)



END PROGRAM csu1dvar
