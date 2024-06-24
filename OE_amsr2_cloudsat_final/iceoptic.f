      SUBROUTINE ICEOPTIC (freqy, temp, epsreal, epsimag)

**    Hufford (1991), see Brussard and Watson (1995), p.297
**
      implicit none

C     Input & output variables; eps the dielectric constant, epsilon
      real     freqy, temp
      real     epsreal, epsimag
      
C     internal variables 
      real     t_ice, theta, A , B    
**
      epsreal = 3.15

      if (temp .gt. 273.16) then
        t_ice = 273.16
      else
        t_ice = temp
      endif

      theta  = 300.0 / t_ice
      A      = 1.0 E-04 * (50.4 + 62.0 * (theta - 1.0))
     +         * EXP (-22.1 * (theta - 1.0))
      B      = 1.0 E-04 * (0.633 / theta - 0.131)
     +         + (7.36 E-04 * theta / (theta - 0.9927))
     +         * (7.36 E-04 * theta / (theta - 0.9927))
      epsimag = A / freqy + B * freqy
**
      return
      end
**
