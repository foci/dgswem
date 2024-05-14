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
      
      INTEGER DOFS1
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
            ENDIF

            HQX_IN = S0X(L)*(HWD_IN**(5.D0/3.D0)) / (manningKWE2(L)*sqrt(S0mag(L)))
            HQY_IN = S0Y(L)*(HWD_IN**(5.D0/3.D0)) / (manningKWE2(L)*sqrt(S0mag(L)))

            DO K = 1,DOFS(L)
               RHS_HWD(K,L,IRK) = RHS_HWD(K,L,IRK) + XFAC(K,I,L,pa)*HQX_IN
     &              + YFAC(K,I,L,pa)*HQY_IN  + SRFAC(K,I,L,pa)*SOURCE_R
            ENDDO
         ENDDO  
16111    format(4(8F,1X))
      ENDDO
      
      RETURN
      END SUBROUTINE
