      SUBROUTINE MG_ELLIPS (FINCL, EMATRIX, EINCL, EMG)

      IMPLICIT NONE
      SAVE

** Maxwell-Garnett formula for effective permittivity of 2-component media 
** (elliptical inclusions) P. Bauer 1996
**
** FINCL     volume fraction of inclusions
** EMATRIX   permittivity of matrix
** EINCL     permittivity of inclusions
**
** EMG       effective permittivity

      REAL    FINCL
      COMPLEX EMATRIX, EINCL, EMG, GAMMA, Q

      Q     = (EINCL / (EINCL - EMATRIX)) 
     +      * CLOG (EINCL / EMATRIX) - 1.0
      GAMMA = 2.0 * EMATRIX * Q / (EINCL - EMATRIX)

      EMG = ((1.0 - FINCL) * EMATRIX + FINCL * GAMMA * EINCL)
     +    / (1.0 - FINCL + FINCL * GAMMA) 

      RETURN
      END
