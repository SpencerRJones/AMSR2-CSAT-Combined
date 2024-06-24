      subroutine watoptic(freqy, temp, salinity, epsreal, epsimag)
**
      implicit   none

C     Input & output variables; eps the dielectric constant, epsilon  
      real     freqy, temp, salinity
      real     epsreal, epsimag

C     internal variables      
      real     freqhz, ctemp, pi
      real     omega, epsstat, trelax, fac1
      real     epshigh
      
      data  epshigh  / 4.90 /
      data  pi  / 3.141592654  /
**
      freqhz = freqy*1.0E+09
      ctemp  = temp - 273.16
      omega  = 2.*pi*freqhz
**
      epsstat = (87.134 E+00 - 1.949 E-01 * ctemp
     +        - 1.276 E-02 * ctemp * ctemp
     +        + 2.491 E-04 * ctemp * ctemp * ctemp)
     +        * (1.0 E+00 + 1.613 E-05 * salinity * ctemp
     +        - 3.656 E-03 * salinity
     +        + 3.210 E-05 * salinity * salinity
     +        - 4.232 E-07 * salinity * salinity * salinity)
**
**
      trelax = (1.768 E-11 - 6.086 E-13 * ctemp
     +       + 1.104 E-14 * ctemp * ctemp
     +       - 8.111 E-17 * ctemp * ctemp * ctemp)
     +       * (1.0  E+00 + 2.282 E-05 * salinity * ctemp
     +       - 7.638 E-04  * salinity
     +       - 7.760 E-06  * salinity * salinity
     +       + 1.105 E-08  * salinity * salinity * salinity)
**
      fac1    = 1.0 + omega*omega*trelax*trelax
      epsreal = epshigh + (epsstat - epshigh) / fac1
      epsimag = ((epsstat - epshigh) * omega * trelax) / fac1
  
      return
      end
