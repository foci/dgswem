C***********************************************************************
C     
C     SUBROUTINE CHANNEL_EDGE_RUNOFF( )
C     
C     This subroutine does the following:
C     
C     
C***********************************************************************

      SUBROUTINE CHANNEL_EDGE_RUNOFF(IT)

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE NodalAttributes, ONLY :  ESLM
      use sizes, only:  sz, myproc, layers

#ifdef CMPI
      USE MESSENGER_ELEM 
#endif

      IMPLICIT NONE

C.....Declare local variables

      INTEGER I,J,K,L
      INTEGER LED_IN, LED_EX, GED, GP_IN, GP_EX,ll,IT
      REAL(SZ), PARAMETER :: ZERO = 1.D-12
      REAL(SZ) TX, TY, W_IN, W_EX

      REAL(SZ) ARK, BRK
      REAL(SZ) MAX_BOA          ! Maximum of beta_il/alpha_il for all l
      REAL(SZ) NLEQG_TMP, G_TMP
      REAL(SZ) F_HAT_O, G_HAT_O, H_HAT_O,i_hat_o,j_hat_o
      REAL(SZ) G_HAT_IN, H_HAT_IN
      REAL(SZ) G_HAT_EX, H_HAT_EX
      REAL(SZ) K_HAT_O
      REAL(SZ) EDFAC_IN, EDFAC_EX

      INTEGER  GED_EX
      REAL(SZ) HWD_IN, HWD_EX, HQX_IN, HQY_IN, HQX_EX, HQY_EX, HQ_HAT_IN, HQ_HAT_EX
      REAL(SZ) AWD_IN, AWD_PLUS, AWD_MINUS, AQ_IN, AQ_MINUS, AQ_PLUS
      REAL(SZ) AWD_H_IN, AWD_H_EX, P_PLUS, P_MINUS
      REAL(SZ) ZE_1D_IN, ZE_1D_PLUS, ZE_1D_MINUS
      REAL(SZ) ZE_1D_H_IN, ZE_1D_H_EX
      REAL(SZ) QX_1D_IN, QX_1D_PLUS, QX_1D_MINUS
      REAL(SZ) HB_1D_IN, HB_1D_PLUS, HB_1D_MINUS
      REAL(SZ) QL1
      REAL(SZ) WCED, CB, CSS
      REAL(SZ) mx,my
      
      INTEGER DOFS1

      REAL(SZ) DEPTH, F1_NL, FU_NL, FV_NL, FG_NL, FH_NL, FW_NL
      REAL(SZ) HUU, HVV, HUV, GH2,MZ_X(layers),MZ_Y(layers), fgauss, sig
      REAL(SZ) SFACQUAD
      

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
            
C.....Set the components for the tangential vector to the edge

            EDFAC_IN = XLEN(GED)/AREAS(EL_IN)
            EDFAC_EX = XLEN(GED)/AREAS(EL_EX)
            
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
               
C.....Compute the numerical flux (upwind scheme)
                            
               HQ_HAT_IN = HQX_IN*NX + HQY_IN*NY
               HQ_HAT_EX = HQX_EX*NX + HQY_EX*NY
               
               DO K = 1,DOFS(EL)
                  
                  W_IN = EDFAC_IN*EDGEQ(K,GP_IN,LED_IN,pa)
                  W_EX = EDFAC_EX*EDGEQ(K,GP_EX,LED_EX,pa)
                  
C                  W_IN = 2.0d0*XLEN(GED)/AREAS(EL_IN)*M_INV(K,pa)*
C     &                     PHI_EDGE(K,GP_IN,LED_IN,pa)*WEGP(GP_IN,pa)
C                  W_EX = 2.0d0*XLEN(GED)/AREAS(EL_EX)*M_INV(K,pa)*
C     &                     PHI_EDGE(K,GP_EX,LED_EX,pa)*WEGP(GP_EX,pa)
     
                  IF (ET(GED).eq.-1) THEN
                     RHS_HWD(K,EL_IN,IRK) = RHS_HWD(K,EL_IN,IRK) - W_IN*HQ_HAT_IN
                     ! RHS_HWD(K,EL_EX,IRK) = RHS_HWD(K,EL_EX,IRK) + W_EX*HQ_HAT_IN
                     IF (NEDCEDN(GED) .gt. 0) THEN
                        W_EX = EDFAC_EX*EDGEQ(K,GP_EX,LED_EX,pa)
                        ! RHS_HWD(K,EL_EX,IRK) = RHS_HWD(K,EL_EX,IRK) + W_EX*HQ_HAT_EX
                        ! stop
                     else
                        RHS_HWD(K,EL_EX,IRK) = RHS_HWD(K,EL_EX,IRK) + W_EX*HQ_HAT_IN
                     ENDIF
                     
                  ELSEIF (ET(GED).eq.0) THEN
                     ! No flow
                  ELSEIF (ET(GED).eq.-2) THEN
                     RHS_HWD(K,EL_IN,IRK) = RHS_HWD(K,EL_IN,IRK) - W_IN*HQ_HAT_IN
                     RHS_HWD(K,EL_EX,IRK) = RHS_HWD(K,EL_EX,IRK) + W_EX*HQ_HAT_EX
                  ENDIF

               ENDDO
               
               ! IF (NEDCEDN(GED).eq.0) THEN
                  ! QL1 = 0.D0
               ! ENDIF
               IF (ET(GED).eq.-2) THEN
                  QL1 = HQ_HAT_IN - HQ_HAT_EX
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
            
C.....Set the components for the tangential vector to the edge

            EDFAC_IN = XLEN(GED)/AREAS(EL_IN)
            EDFAC_EX = XLEN(GED)/AREAS(EL_EX)
            
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
               
C.....Compute the numerical flux (upwind scheme)
                            
               HQ_HAT_IN = HQX_IN*NX + HQY_IN*NY
               HQ_HAT_EX = HQX_EX*NX + HQY_EX*NY
               
               DO K = 1,DOFS(EL)
                  
                  W_IN = EDFAC_IN*EDGEQ(K,GP_IN,LED_IN,pa)
                  W_EX = EDFAC_EX*EDGEQ(K,GP_EX,LED_EX,pa)
                  
C                  W_IN = 2.0d0*XLEN(GED)/AREAS(EL_IN)*M_INV(K,pa)*
C     &                     PHI_EDGE(K,GP_IN,LED_IN,pa)*WEGP(GP_IN,pa)
C                  W_EX = 2.0d0*XLEN(GED)/AREAS(EL_EX)*M_INV(K,pa)*
C     &                     PHI_EDGE(K,GP_EX,LED_EX,pa)*WEGP(GP_EX,pa)
     
                  IF (ET(GED).eq.-1) THEN
                     RHS_HWD(K,EL_IN,IRK) = RHS_HWD(K,EL_IN,IRK) - W_IN*HQ_HAT_IN
                     ! RHS_HWD(K,EL_EX,IRK) = RHS_HWD(K,EL_EX,IRK) + W_EX*HQ_HAT_IN
                     IF (NEDCEDN(GED) .gt. 0) THEN
                        W_EX = EDFAC_EX*EDGEQ(K,GP_EX,LED_EX,pa)
                        ! RHS_HWD(K,EL_EX,IRK) = RHS_HWD(K,EL_EX,IRK) + W_EX*HQ_HAT_EX
                        ! stop
                     else
                        RHS_HWD(K,EL_EX,IRK) = RHS_HWD(K,EL_EX,IRK) + W_EX*HQ_HAT_IN
                     ENDIF
                     
                  ELSEIF (ET(GED).eq.0) THEN
                     ! No flow
                  ELSEIF (ET(GED).eq.-2) THEN
                     RHS_HWD(K,EL_IN,IRK) = RHS_HWD(K,EL_IN,IRK) - W_IN*HQ_HAT_IN
                     RHS_HWD(K,EL_EX,IRK) = RHS_HWD(K,EL_EX,IRK) + W_EX*HQ_HAT_EX
                  ENDIF

               ENDDO
               
               IF (ET(GED).eq.-2) THEN
                  QL1 = HQ_HAT_IN - HQ_HAT_EX
               ELSE
                  QL1 = 0.D0
               ENDIF
               
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
      
#ifdef CMPI
      CALL UPDATER_LATERALINFLOW(QL,HWD,IRK+1,1)
#endif

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
         
         IF (DOFS(EL).EQ.1) THEN
            DOFS1 = 1
         ELSE
            DOFS1 = pa+1
         ENDIF

         N2 = NCEDNO(2,I)
         
         DO J = 1,NEGP(pa)
            GP_IN = J
            AWD_IN = AWD(1,GED,IRK)
            DO K = 2,DOFS1
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
            if (abs(CSS).ge.1e8) then
            P_MINUS = CB + 2.D0*AWD_H_IN
            else
            P_MINUS = CB + 2.D0*AWD_H_IN*sqrt(1.D0+1.D0/(CSS**2.D0))
            endif
            
               
            if (AWD_IN.eq.0.D0) then
               AQ_IN = 0.D0
            else
               ! manning = 0.15, alpha = 1/n, beta = 2/3
               AQ_IN = (1.D0/manningKWE1(I))*sqrt(CEDSLOPE(I))*((AWD_IN/P_MINUS)**(2.D0/3.D0)) * AWD_IN
            endif
            
            DO K = 1,DOFS1
               RHS_AWD(K,GED,IRK) = RHS_AWD(K,GED,IRK) + (2.d0/XLEN(GED))*M1D_INV(K,pa)
     &                             *DXDPHI1D(K,GP_IN,pa)*WEGP(GP_IN,pa)*AQ_IN
     &                             + M1D_INV(K,pa)*PHI1D(K,GP_IN,pa)*WEGP(GP_IN,pa)*QL(GP_IN,I)
            ENDDO
            
         ENDDO
         
C.....Compute Numerical Flux for 1D KWE
         
         AWD_MINUS = AWD(1,GED,IRK)
         ! AWD_PLUS = AWD(1,II,IRK)
         
         DO K = 2,DOFS1
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
         
         if (abs(CSSLOPE(2,I)).ge.1e8) then
         P_MINUS = CBWIDTH(2,I) + 2.D0*AWD_H_IN
         else
         P_MINUS = CBWIDTH(2,I) + 2.D0*AWD_H_IN*sqrt(1.D0+1.D0/(CSSLOPE(2,I)**2.D0))
         endif
         
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
         DO K = 1,DOFS1
            RHS_AWD(K,GED,IRK) = RHS_AWD(K,GED,IRK) - (2.d0/XLEN(GED))*M1D_INV(K,pa)*PHI1D_MINUS(K,pa)*AQ_MINUS
            ! WRITE(*,*) K,RHS_AWD(K,GED,IRK)
         ENDDO
         
         ! II = NCEDED(2,I)
         
         if (NCEDEX(1,I).gt.0) then
            do J = 1,NCEDEX(1,I)
               GED_EX = NCEDEX(J+1,I)
               DO K = 1,DOFS1
                  RHS_AWD(K,GED_EX,IRK) = RHS_AWD(K,GED_EX,IRK)
     &             + (2.d0/XLEN(GED_EX))*M1D_INV(K,pa)*PHI1D_PLUS(K,pa)*AQ_MINUS/NCEDEX(1,I)
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
      
      END SUBROUTINE
