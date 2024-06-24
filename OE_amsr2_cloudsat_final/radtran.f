      SUBROUTINE RADTRAN(nlyr, view_angle, Tsfc, lyrhgt, lyrtemp, kext,
     >                   salb, asym, emis, ebar, refl, Tb, Tbdown)

C
C     Based on Kummerow PHD THESIS PLANE PARALLEL ATMOS, BUT
C     A) ASSUMING NO SCATTERING  
C     B) DOWNWELLING RADIANCE RECALCULATED AND CHANGED
C        BY GREG ELSAESSER (MAR 2005)
C
C     New: Modified by Rick Schulte to once again include scattering

      implicit  none
C      include  'parameters.inc'

      integer   nlyr
      real      view_angle, Tsfc, emis, refl, Tb, Tbdown
      real      lyrhgt(0:nlyr), lyrtemp(0:nlyr), kext(nlyr)    
      real      B0(nlyr), B1(nlyr)
      real      Z(0:nlyr), Iout(0:nlyr)
      real      umu

      real    XA, XB, YA, YB, dz, TERM1, TERM2, TERM3, TERM4, TERM5
      real    L(nlyr), H(nlyr)
      integer j
      
      ! New variables
      real      salb(nlyr), asym(nlyr)
      real      ebar
      real      Iin(nlyr+1,100)
      real  W(2*nlyr,2*nlyr), BB(2*nlyr), DP(nlyr), DM(nlyr)
      real  MU, NU
      real  FISOT, RCOND
      real  XNU, XC, XD, YC, DDZ, XIUP
      integer I, NANG, NN
      
      data FISOT / 2.7 /
C

      umu = cos(view_angle*3.14159/180.0)
      !write(6,'("nlyr = ",i2)') nlyr
      !write(6,'("view_angle  = ",f7.4)') view_angle
      !write(6,'("tsfc = ",f7.2)') tsfc
      !write(6,'("hgt  = ",20(f7.2))') lyrhgt
      !write(6,'("temp = ",20(f7.2))') lyrtemp
      !write(6,'("kext = ",20(f7.4))') kext
      !write(6,'("salb = ",20(f7.4))') salb
      !write(6,'("asym = ",20(f7.4))') asym
      !write(6,'("emis = ",f7.4)') emis
      !write(6,'("refl = ",f7.4)') refl
      W(:,:) = 0.
      BB(:) = 0.
      RCOND = 0.
      Z(0) = lyrhgt(0)
      do j = 1, nlyr
        Z(j)  = lyrhgt(j)
        B0(j) = lyrtemp(j-1)
        B1(j) = (lyrtemp(j) - lyrtemp(j-1))/((lyrhgt(j) - lyrhgt(j-1) + 
     $           0.00001))
	    L(j) = sqrt(3.*kext(j)*kext(j)*(1.-salb(j))*(1.-salb(j)*asym(j)))
	    H(j) = 1.5 * kext(j) * (1.-salb(j)*asym(j))
      end do

C     FILL IN THE NON-ZERO MATRIX ELEMENTS
      W(1,1)   = ((EBAR - 2.)*L(1)/H(1)) + EBAR
      W(1,2)   = ((2. - EBAR)*L(1)/H(1)) + EBAR
      do I = 2,2*(NLYR-1),2
        W(I,I-1)   =  (1. - L(I/2)/H(I/2))*EXP(+L(I/2)*(Z(I/2)-
     $                Z(I/2-1)))
        W(I,I  )   =  (1. + L(I/2)/H(I/2))*EXP(-L(I/2)*(Z(I/2)-
     $                Z(I/2-1)))
        W(I,I+1)   = -(1. - L(I/2+1)/H(I/2+1))
        W(I,I+2)   = -(1. + L(I/2+1)/H(I/2+1))

        W(I+1,I-1) =  (1. + L(I/2)/H(I/2))*EXP(+L(I/2)*(Z(I/2)-
     $                Z(I/2-1)))
        W(I+1,I)   =  (1. - L(I/2)/H(I/2))*EXP(-L(I/2)*(Z(I/2)-
     $                Z(I/2-1)))
        W(I+1,I+1) = -(1. + L(I/2+1)/H(I/2+1))
        W(I+1,I+2) = -(1. - L(I/2+1)/H(I/2+1))
      end do
      W(2*NLYR,2*NLYR-1) =  (1. + L(NLYR)/H(NLYR))*EXP(+L(NLYR)*
     $                      (Z(NLYR)-Z(NLYR-1)))
      W(2*NLYR,2*NLYR)   =  (1. - L(NLYR)/H(NLYR))*EXP(-L(NLYR)*
     $                      (Z(NLYR)-Z(NLYR-1)))

C     FILL IN THE ROW OF CONSTANTS IN THE LINEAR EQUATIONS
      BB(1)    = EBAR*Tsfc - EBAR*B0(1) - 
     $           (EBAR - 2.)*B1(1)/H(1)
      do I = 2,2*(NLYR-1),2
        BB(I)   =  + B1(I/2)/H(I/2) - B1(I/2+1)/H(I/2+1)
        BB(I+1) =  - B1(I/2)/H(I/2) + B1(I/2+1)/H(I/2+1)
      end do
      BB(2*NLYR)  =  FISOT - B0(NLYR) - B1(NLYR)*(Z(NLYR) - 
     $               Z(NLYR-1) + 1/H(NLYR))
C
C     MATRIX INVERSION IN DONE IN SUBROUTINE LINPAK
      CALL LINPAK(NLYR, W, BB, RCOND)
C
      do I = 1,NLYR
        DP(I) = BB(2*I-1)
        DM(I) = BB(2*I)
      end do
C     AFTER D'S ARE KNOWN, CALCULATE SURFACE RADIANCE
      MU =  umu
      NU = -MU
C      
C     FOR THE FOLLOWING CALCULATIONS, REFER TO APPENDIX B OF THESIS
C     *************************************************************
C     CALCULATE THE DOWNWELLING FLUX AT ANGLE MU ONLY (Assuming "lambert" = false)
      NN = 22            ! THIS IS A DUMMY INDEX FOR I_IN
      IIN(NLYR+1,NN) = FISOT
C     LOOP THROUGH THE REMAINING LAYERS
      do J = NLYR,1,-1
C       CALCULATE RADIANCE FROM TOP OF LAYER "J"
        XA = B0(J) - 1.5*SALB(J)*ASYM(J)*NU*B1(J)/H(J)
        XB = B1(J)
        XC = SALB(J)*DP(J)*(1. - 1.5*ASYM(J)*NU*L(J)/H(J))               
        XD = SALB(J)*DM(J)*(1. + 1.5*ASYM(J)*NU*L(J)/H(J))
        YA = KEXT(J)/NU
        YB = YA + L(J)
        YC = YA - L(J)
        DDZ = Z(J) - Z(J-1)
	      
        TERM1 = IIN(J+1,NN)*EXP(YA*DDZ)
        TERM2 = XA*(1. - EXP(YA*DDZ))
        IF (ABS(YA*DDZ) .GT. 1.E-05) THEN
          TERM3 = XB/YA*(EXP(YA*DDZ) - YA*DDZ - 1.)
          ! previous to 7/2007 fix was:
          ! TERM3 = XB/YA*(EXP(YA*DZ)*(1. - YA*DZ) - 1.)
        ELSE
          TERM3 = 0.
          ! previous to 7/2007 was: TERM3=-XB*YA*DZ*DZ
        ENDIF
        IF (ABS(YB*DDZ) .GT. 1.E-05 ) THEN
          TERM4 = XC*YA/YB*(1. - EXP(YB*DDZ))
        ELSE
          TERM4 = -XC*YA*DDZ
        ENDIF
        IF (ABS(YC*DDZ) .GT. 1.E-05 ) THEN
          TERM5 = XD*YA/YC*(1. - EXP(YC*DDZ))
        ELSE
          TERM5 = -XD*YA*DDZ
        ENDIF
        IIN(J,NN) = TERM1 + TERM2 + TERM3 + TERM4 + TERM5
      end do
     
      XIUP = IIN(1,22)
      IOUT(0) = EMIS*Tsfc + refl*XIUP
      
!      write(6,'("TB_down = ",f8.2,", emis = ",f8.2,", refl = ",f8.2)')
!     >          Iin(1,NN),EMIS*Tsfc,Refl*Iin(1,NN)
     
      do J = 1,NLYR
C       CALCULATE THE UPWELLING RADIANCES AT THE TOP OF EACH LAYER J
        XA = B0(J) - 1.5*SALB(J)*ASYM(J)*MU*B1(J)/H(J)
        XB = B1(J)
        XC = SALB(J)*DP(J)*(1. - 1.5*ASYM(J)*MU*L(J)/H(J))               
        XD = SALB(J)*DM(J)*(1. + 1.5*ASYM(J)*MU*L(J)/H(J))
        YA = KEXT(J)/MU
        YB = YA + L(J)
        YC = YA - L(J)
        DDZ = Z(J) - Z(J-1)
	    	    
        TERM1 = IOUT(J-1)*EXP(-YA*DDZ)
        TERM2 = XA*(1. - EXP(-YA*DDZ))
        IF ( ABS(YA*DDZ) .GT. 1.E-05 ) THEN
          TERM3 = XB/YA*(EXP(-YA*DDZ) + YA*DDZ - 1.)
        ELSE
          TERM3 = 0.
        ENDIF
        IF ( ABS(YB*DDZ) .GT. 1.E-05 ) THEN
          TERM4 = XC*YA/YB*(EXP( (YB-YA)*DDZ ) - EXP(-YA*DDZ) )
        ELSE            
          TERM4 = XC*YA*DDZ*EXP(-YA*DDZ)
        ENDIF
        IF (ABS(YC*DDZ) .GT. 1.E-05 ) THEN
          TERM5 = XD*YA/YC*EXP(-YA*DDZ)*(EXP(YC*DDZ) - 1.)
        ELSE
          TERM5 = XD*YA*DDZ*EXP(-YA*DDZ)
        ENDIF
        IOUT(J) = TERM1 + TERM2 + TERM3 + TERM4 + TERM5
      end do
C
      TB = IOUT(NLYR)

      RETURN
      END
