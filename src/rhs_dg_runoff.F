C***********************************************************************
C     
C     SUBROUTINE RHS_DG_RUNOFF()
C     
C     This subroutine does:
C
C     1. Computes the area integrals for the DG-SAKE and
C     adds them into the RHS.
C     2. Add lateral inflows into 1D KWE.
C     
C     Written by Younghun Kang
C     
C***********************************************************************
      SUBROUTINE RHS_DG_RUNOFF(IT)
      
C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE NodalAttributes, ONLY : TAU, IFLINBF, IFHYBF, HBREAK, FTHETA,
     &     FGAMMA,LoadManningsN,ManningsN,CF

      USE sizes, ONLY: myproc,layers

#ifdef CMPI
      USE MESSENGER_ELEM 
#endif

      IMPLICIT NONE
      
C.....Declare local variables
      INTEGER IT
      INTEGER L,k,i,ll, J, II
      REAL(SZ) DPSIDX(3), DPSIDY(3)
      REAL(SZ) AREA, IMASS, TKX, TKY, Xpart, Ypart,H_0,C_1
      REAL(SZ) PHI_AREA_KI,MN_IN, MassAction1st,MassAction2nd,fx,fy
      REAL(SZ) LZ_XX, LZ_XY, LZ_YX, LZ_YY, rate, s_mass, s_sed,b_0
      REAL(SZ) DEPTH, F1_NL, FU_NL, FV_NL, FG_NL, FH_NL, FW_NL
      REAL(SZ) HUU, HVV, HUV, GH2,MZ_X(layers),MZ_Y(layers), fgauss, sig
      REAL(SZ) DEPTH_C, FH_NL_C, UX_C, UY_C, UMAG_C, DTDPH,SFACQUAD
      Real(SZ) discharge_modelX_IN,discharge_modelY_IN
      Real(SZ) DH_X,DH_Y,phi_tot,C_0,HZ_X,HZ_Y,TZ_X,TZ_Y

      Real(SZ),allocatable :: XBCbt(:),YBCbt(:)

      INTEGER GED, LED_IN, LED_EX, GP_IN, GP_EX, GED_EX
      REAL(SZ) HWD_IN, HWD_EX, HQX_IN, HQY_IN, HQX_EX, HQY_EX
      REAL(SZ) AWD_IN, AWD_PLUS, AWD_MINUS, AQ_IN, AQ_MINUS, AQ_PLUS
      REAL(SZ) AWD_H_IN, AWD_H_EX, P_PLUS, P_MINUS
      REAL(SZ) QL1
      REAL(SZ) WCED, CB, CSS
      REAL(SZ) mx,my
      
      ! INTEGER N1,N2
      
      Allocate ( XBCbt(MNE),YBCbt(MNE) )
      
      
      
      DTDPH = 1.D0/DTDP
         
      DO L = 1, NE
         
C.......Adjust the p values for constants
         
         pa = PDG_EL(L)

#ifdef P0
         if (pa.eq.0) then
            pa = 1
         endif
#endif
         
C.......If element is dry then skip calculations
         
         ! IF (WDFLG(L).EQ.0) then
            ! GOTO 1000
         ! endif
         
C.......Retrieve the global node numbers for the element

         N1 = NM(L,1)
         N2 = NM(L,2)
         N3 = NM(L,3)

C.......Compute avaraged values
C.......These will be used later when bottom friction is computed

         ! DEPTH_C = HB(1,L,1) + ZE(1,L,IRK)
         ! FH_NL_C = 1.D0/(NLEQ*DEPTH_C + LEQ)
         ! UX_C = QX(1,L,IRK)*FH_NL_C
         ! UY_C = QY(1,L,IRK)*FH_NL_C
         ! UMAG_C = SQRT(UX_C*UX_C + UY_C*UY_C)
         
C.......Compute derivatives of Lagrange basis functions at nodes
         
         IF ((NWS.NE.0).OR.(NTIP.NE.0)) THEN
            DPSIDX(1) = DRPSI(1)*DRDX(L) + DSPSI(1)*DSDX(L)
            DPSIDX(2) = DRPSI(2)*DRDX(L) + DSPSI(2)*DSDX(L)
            DPSIDX(3) = DRPSI(3)*DRDX(L) + DSPSI(3)*DSDX(L)
            DPSIDY(1) = DRPSI(1)*DRDY(L) + DSPSI(1)*DSDY(L)
            DPSIDY(2) = DRPSI(2)*DRDY(L) + DSPSI(2)*DSDY(L)
            DPSIDY(3) = DRPSI(3)*DRDY(L) + DSPSI(3)*DSDY(L)
         ENDIF

C.......Compute ZE, QX, QY, and HB at each area Gauss quadrature point
         IF (RAINFALL.GT.0) THEN
            CALL PRECIPITATION(IT*DTDP,L,SOURCE_R)
         ENDIF
         
         DO I = 1,NAGP(pa)
            
            HWD_IN = HWD(1,L,IRK)

            DO K = 2,DOFS(L)
               
               HWD_IN = HWD_IN + HWD(K,L,IRK)*PHI_AREA(K,I,pa)

            ENDDO
            
            IF (HWD_IN .lt. 0.D0) THEN
               HWD_IN = 0.D0
               ! DO K = 1,DOFS(L)
                  ! HWD(K,L,IRK) = 0.D0
               ! ENDDO
            ENDIF

            ! manning_n = 0.015
            HQX_IN = S0X(L)*(HWD_IN**(5.D0/3.D0)) / (manningKWE2(L)*sqrt(S0mag(L)))
            HQY_IN = S0Y(L)*(HWD_IN**(5.D0/3.D0)) / (manningKWE2(L)*sqrt(S0mag(L)))

            DO K = 1,DOFS(L)
               RHS_HWD(K,L,IRK) = RHS_HWD(K,L,IRK) + XFAC(K,I,L,pa)*HQX_IN
     &              + YFAC(K,I,L,pa)*HQY_IN  + SRFAC(K,I,L,pa)*SOURCE_R
            ENDDO
         ENDDO  
16111    format(4(8F,1X))
      ENDDO
      
C***********************************************************************
C
C Compute 1D KWE
C
C***********************************************************************

C.....Compute lateral inflow for 1D KWE
      
      DO I = 1,NCEDS
         GED = NCEDN(I)
         
         if (flipNEDEL(GED).eq.0) then
            EL_IN = NEDEL(1,GED)
            EL_EX = NEDEL(2,GED)
            
            LED_IN = NEDSD(1,GED)
            LED_EX = NEDSD(2,GED)
            
            IF (EL_EX.EQ.0) then
               EL_EX = EL_IN
               LED_EX = LED_IN
            ENDIF
            
            EL = EL_EX
              
            if (DOFS(EL_EX).LT.DOFS(EL_IN)) then
               EL = EL_IN
            endif
            
            pa = PDG_EL(EL)
              
            NX = COSNX(GED)
            NY = SINNX(GED)
            
            DO J = 1,NEGP(pa)
               if (NCEDNO(3,I).eq.0) then
                  GP_IN = J
                  GP_EX = NEGP(pa) - J + 1
               else
                  GP_IN = NEGP(pa) - J + 1
                  GP_EX = J
               endif
                  
               HWD_IN = HWD(1,EL_IN,IRK)
               HWD_EX = HWD(1,EL_EX,IRK)
               
               DO K = 2,DOFS(EL_IN)
                  HWD_IN = HWD_IN + HWD(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN,pa)
                  HWD_EX = HWD_EX + HWD(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX,pa)
               ENDDO
               
               IF (HWD_IN.lt.0.D0) THEN
                  HWD_IN = 0.D0
               ENDIF
               
               IF (HWD_EX.lt.0.D0) THEN
                  HWD_EX = 0.D0
               ENDIF
               
               HQX_IN = S0X(EL_IN)*(HWD_IN**(5.D0/3.D0)) / (manningKWE2(EL_IN)*sqrt(S0mag(EL_IN)))
               HQY_IN = S0Y(EL_IN)*(HWD_IN**(5.D0/3.D0)) / (manningKWE2(EL_IN)*sqrt(S0mag(EL_IN)))
               HQX_EX = S0X(EL_EX)*(HWD_EX**(5.D0/3.D0)) / (manningKWE2(EL_EX)*sqrt(S0mag(EL_EX)))
               HQY_EX = S0Y(EL_EX)*(HWD_EX**(5.D0/3.D0)) / (manningKWE2(EL_EX)*sqrt(S0mag(EL_EX)))
               
               IF (EL_EX.EQ.EL_IN .AND. LED_EX.EQ.LED_IN) THEN
                  HQX_EX = 0.d0
                  HQY_EX = 0.d0
               ENDIF
               ! IF (NEDCEDN(GED).eq.0) THEN
                  ! QL1 = 0.D0
               ! ENDIF
               IF (ET(GED).eq.-1) THEN
                  ! QL1 = (HQX_IN*NX + HQY_IN*NY)
                  QL1 = 0.D0
               ELSEIF (ET(GED).eq.-2) THEN
                  QL1 = (HQX_IN*NX + HQY_IN*NY) - (HQX_EX*NX + HQY_EX*NY)
               ELSE
                  QL1 = 0.D0
               ENDIF
               
            
               ! QL1 = (HQX_IN*NX + HQY_IN*NY) - (HQX_EX*NX + HQY_EX*NY)
               QL(J,I) = QL1
               ! write(9999,*),J,I,QL(J,I)
               
               ! if (isnan(QL(j,i))) then
                  ! print*, HWD_IN, HWD_EX
                  ! stop
               ! endif
               
            ENDDO
         
         
         elseif (flipNEDEL(GED).eq.1) then
            EL_IN = NEDEL(2,GED)
            EL_EX = NEDEL(1,GED)
            
            LED_IN = NEDSD(2,GED)
            LED_EX = NEDSD(1,GED)
            
            IF (EL_EX.EQ.0) then
               EL_EX = EL_IN
               LED_EX = LED_IN
            ENDIF
            
            EL = EL_EX
              
            if (DOFS(EL_EX).LT.DOFS(EL_IN)) then
               EL = EL_IN
            endif
            
            pa = PDG_EL(EL)
              
            NX = -COSNX(GED)
            NY = -SINNX(GED)
            
            DO J = 1,NEGP(pa)
               if (NCEDNO(3,I).eq.1) then
                  GP_IN = J
                  GP_EX = NEGP(pa) - J + 1
               else
                  GP_IN = NEGP(pa) - J + 1
                  GP_EX = J
               endif
                  
               HWD_IN = HWD(1,EL_IN,IRK)
               HWD_EX = HWD(1,EL_EX,IRK)
               
               DO K = 2,DOFS(EL_IN)
                  HWD_IN = HWD_IN + HWD(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN,pa)
                  HWD_EX = HWD_EX + HWD(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX,pa)
               ENDDO
               
               IF (HWD_IN.lt.0.D0) THEN
                  HWD_IN = 0.D0
               ENDIF
               
               IF (HWD_EX.lt.0.D0) THEN
                  HWD_EX = 0.D0
               ENDIF
               
               HQX_IN = S0X(EL_IN)*(HWD_IN**(5.D0/3.D0)) / (manningKWE2(EL_IN)*sqrt(S0mag(EL_IN)))
               HQY_IN = S0Y(EL_IN)*(HWD_IN**(5.D0/3.D0)) / (manningKWE2(EL_IN)*sqrt(S0mag(EL_IN)))
               HQX_EX = S0X(EL_EX)*(HWD_EX**(5.D0/3.D0)) / (manningKWE2(EL_EX)*sqrt(S0mag(EL_EX)))
               HQY_EX = S0Y(EL_EX)*(HWD_EX**(5.D0/3.D0)) / (manningKWE2(EL_EX)*sqrt(S0mag(EL_EX)))
            
               ! IF (NEDCEDN(GED).eq.0) THEN
                  ! QL1 = 0.D0
               ! ENDIF
               
               IF (EL_EX.EQ.EL_IN .AND. LED_EX.EQ.LED_IN) THEN
                  HQX_EX = 0.d0
                  HQY_EX = 0.d0
               ENDIF
               
               IF (ET(GED).eq.-1) THEN
                  ! QL1 = (HQX_IN*NX + HQY_IN*NY)
                  QL1 = 0.D0
               ELSEIF (ET(GED).eq.-2) THEN
                  QL1 = (HQX_IN*NX + HQY_IN*NY) - (HQX_EX*NX + HQY_EX*NY)
               ELSE
                  QL1 = 0.D0
               ENDIF
               
               ! QL1 = (HQX_IN*NX + HQY_IN*NY) - (HQX_EX*NX + HQY_EX*NY)
               QL(J,I) = QL1
               
            ENDDO
         endif
      ENDDO
      
      ! ! Add artifical lateral inflow (to test Iwagaki experiment with lateral inflow condition)
      ! if (IT*DTDP .le. 20) then
         ! DO J = 1,NEGP(pa)
            ! DO I = 1,NCEDS
               ! GED = NCEDN(I)
               ! if (abs(CEDSLOPE(I) - 0.02) .lt. 1e-3) then
                  ! QL(J,I) = 0.1080
               ! elseif (abs(CEDSLOPE(I) - 0.015) .lt. 1e-3) then
                  ! QL(J,I) = 0.0638
               ! elseif (abs(CEDSLOPE(I) - 0.01) .lt. 1e-3) then
                  ! QL(J,I) = 0.0800
               ! endif
               ! QL(J,I) = QL(J,I)*0.196*1e-2
            ! ENDDO
         ! enddo
      ! else
         ! ! DO J = 1,NEGP(pa)
            ! ! DO I = 1,NCEDS
              ! ! QL(J,I) = 0.d0
            ! ! ENDDO
         ! ! enddo
      ! endif
      ! print*,QL
      

C.....Compute RHS for 1D KWE

      DO I = 1,NCEDS
         GED = NCEDN(I)
         
         EL_IN = NEDEL(1,GED)
         EL_EX = NEDEL(2,GED)
         
         LED_IN = NEDSD(1,GED)
         LED_EX = NEDSD(2,GED)
         
         ! EL = EL_EX
         EL = EL_IN
           
         ! if (DOFS(EL_EX).LT.DOFS(EL_IN)) then
            ! EL = EL_IN
         ! endif
         
         pa = PDG_EL(EL)
         
         N2 = NCEDNO(2,I)
         
         DO J = 1,NEGP(pa)
            GP_IN = J
            AWD_IN = AWD(1,GED,IRK)
            DO K = 2,pa+1
               AWD_IN = AWD_IN + AWD(K,GED,IRK)*PHI1D(K,GP_IN,pa)
            ENDDO
            
            if (AWD_IN .lt. 1e-16) then
               AWD_IN = 0.D0
            endif
            
            WCED = (XEGP(J,pa) + 1.D0)/2.D0
            CB = (1.D0 - WCED)*CBWIDTH(1,I) + WCED*CBWIDTH(2,I)
            CSS = (1.D0 - WCED)*CSSLOPE(1,I) + WCED*CSSLOPE(2,I)
            
            if (abs(CSS).ge.1e8) then
               AWD_H_IN = AWD_IN/CB
               ! print*,CB
            elseif (CSS.gt.0) then
               AWD_H_IN = CSS/2.D0*(-CB + sqrt(CB**2.D0 + 4.D0/CSS*AWD_IN))
            else
               AWD_H_IN = CSS/2.D0*(-CB - sqrt(CB**2.D0 + 4.D0/CSS*AWD_IN))
            endif
            
            ! Compute wetted perimeter
            P_MINUS = CB + 2.D0*AWD_H_IN*sqrt(1.D0+1.D0/(CSS**2.D0))            
               
            if (AWD_IN.eq.0.D0) then
               AQ_IN = 0.D0
            else
               ! manning = 0.15, alpha = 1/n, beta = 2/3
               AQ_IN = (1.D0/manningKWE1(I))*sqrt(CEDSLOPE(I))*((AWD_IN/P_MINUS)**(2.D0/3.D0)) * AWD_IN
            endif
            
            DO K = 1,ph+1
               RHS_AWD(K,GED,IRK) = RHS_AWD(K,GED,IRK) + (2/XLEN(GED))*M1D_INV(K,pa)
     &                             *DXDPHI1D(K,GP_IN,pa)*WEGP(GP_IN,pa)*AQ_IN
     &                             + M1D_INV(K,pa)*PHI1D(K,GP_IN,pa)*WEGP(GP_IN,pa)*QL(GP_IN,I)
            ENDDO
            
         ENDDO
         
C.....Compute Numerical Flux for 1D KWE
         
         AWD_MINUS = AWD(1,GED,IRK)
         ! AWD_PLUS = AWD(1,II,IRK)
         
         DO K = 2,ph+1
            AWD_MINUS = AWD_MINUS + AWD(K,GED,IRK)*PHI1D_MINUS(K,pa)
            ! AWD_PLUS = AWD_PLUS + AWD(K,II,IRK)*PHI1D_PLUS(K,pa)
         ENDDO
         
         if (AWD_MINUS .lt. 1e-16) then
            AWD_MINUS = 0.D0
         endif
         
         
         if (abs(CSSLOPE(2,I)).ge.1e8) then
            AWD_H_IN = AWD_MINUS/CBWIDTH(2,I)
         elseif (CSSLOPE(2,I).gt.0) then
            AWD_H_IN = CSSLOPE(2,I)/2.D0*(-CBWIDTH(2,I)
     &               + sqrt(CBWIDTH(2,I)**2.D0 + 4.D0/CSSLOPE(2,I)*AWD_MINUS))
         else
            AWD_H_IN = CSSLOPE(2,I)/2.D0*(-CBWIDTH(2,I)
     &               - sqrt(CBWIDTH(2,I)**2.D0 + 4.D0/CSSLOPE(2,I)*AWD_MINUS))
         endif
     
         P_MINUS = CBWIDTH(2,I) + 2.D0*AWD_H_IN*sqrt(1.D0+1.D0/(CSSLOPE(2,I)**2.D0))
         
         ! manning = 0.15, alpha = 1/n, beta = 2/3
         if (AWD_H_IN.eq.0.D0) then
            AQ_MINUS = 0.D0
         else
            AQ_MINUS = (1.D0/manningKWE1(I))*sqrt(CEDSLOPE(I))
     &                         *((AWD_MINUS/P_MINUS)**(2.D0/3.D0)) * AWD_MINUS
         endif
         
         
         ! if (isnan(AQ_MINUS)) then
            ! print*,'flux',AWD_MINUS,AWD_H_IN,P_MINUS
            ! stop
         ! endif
C         AWD_H_EX = (-HYDR_INFO(1,N2) + 
C     &                 sqrt(HYDR_INFO(1,N2)**2.D0 + 4.D0*HYDR_INFO(2,N2)*AWD_PLUS))
C     &                    / (2.D0*HYDR_INFO(2,N2))
C     
C         P_PLUS = HYDR_INFO(1,N2) + 2.D0*AWD_H_EX*sqrt(1+HYDR_INFO(2,N2)**2.D0)
C         ! manning = 0.15, alpha = 1/n, beta = 2/3
C         AQ_PLUS = (1.D0/manningKWE1(II))*sqrt(CEDSLOPE(II))
C     &                      *((AWD_PLUS/P_PLUS)**(2.D0/3.D0)) * AWD_PLUS
         ! WRITE(*,*) AQ_MINUS
         DO K = 1,ph+1
            RHS_AWD(K,GED,IRK) = RHS_AWD(K,GED,IRK) - (2/XLEN(GED))*M1D_INV(K,pa)*PHI1D_MINUS(K,pa)*AQ_MINUS
            ! WRITE(*,*) K,RHS_AWD(K,GED,IRK)
         ENDDO
         
         ! II = NCEDED(2,I)
         
         if (NCEDEX(1,I).gt.0) then
            do J = 1,NCEDEX(1,I)
               GED_EX = NCEDEX(J+1,I)
               DO K = 1,ph+1
                  RHS_AWD(K,GED_EX,IRK) = RHS_AWD(K,GED_EX,IRK) + (2/XLEN(GED_EX))*M1D_INV(K,pa)*PHI1D_PLUS(K,pa)*AQ_MINUS
               ENDDO
            enddo
         endif
         
         ! if (I.eq.3) then
            ! ! do J = 1,NCEDEX(1,I)
               ! ! GED_EX = NCEDEX(J+1,I)
               ! ! print*,GED_EX
               ! DO K = 1,ph+1
                  ! RHS_AWD(K,GED,IRK) = RHS_AWD(K,GED,IRK) + (2/XLEN(GED))*M1D_INV(K,pa)*PHI1D_PLUS(K,pa)*AQ_MINUS
               ! ENDDO
            ! ! enddo
         ! endif
         
         ! IF (II .gt. 0) THEN
            ! DO K = 1,ph+1
               ! RHS_AWD(K,II,IRK) = RHS_AWD(K,II,IRK) + (2/XLEN(II))*M1D_INV(K,pa)*PHI1D_PLUS(K,pa)*AQ_MINUS
            ! ENDDO
         ! ENDIF
      ENDDO
      
      RETURN
      END SUBROUTINE
