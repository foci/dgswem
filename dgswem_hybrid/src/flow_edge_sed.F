C***********************************************************************
C
C     SUBROUTINE FLOW_EDGE_SED( )
C
C     This subroutine does the following:
C
C       1.  Calculates the values of the necessary variables at the edge
C           gauss points for NON-ZERO FLUX edges
C       2.  Calls the appropriate subroutine to compute the flux at
C           these points.
C       3.  Calls the appropriate subroutine to compute the boundary
C           integrals.
C
C     Written by Ethan Kubatko (03-05-2005)
C
C***********************************************************************

      SUBROUTINE FLOW_EDGE_SED(L)

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE SED

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LED, GED,i,k

C.....Retrieve the global and local edge number

      GED = NFEDN(L)
      LED = NEDSD(1,GED)

C.....Retrieve the elements which share the edge

      EL_IN = NEDEL(1,GED)

C.....Retrieve the nodes of the given element

      N1 = NEDNO(1,GED)
      N2 = NEDNO(2,GED)

C.....Retrieve the components of the normal vector to the edge

      NX = COSNX(GED)
      NY = SINNX(GED)

C.....Compute ZE, QX, QY, and HB at each edge Gauss quadrature point

      DO 111 I = 1,NEGP

        ZE_IN = 0.D0
        QX_IN = 0.D0
        QY_IN = 0.D0

        HB_IN = 0.D0
        HB_EX = 0.D0

C.....Compute the solution at the interior state

        DO K = 1,DOF

          ZE_IN = ZE_IN + ZE(K,EL_IN,1  )*PHI_EDGE(K,I,LED)
          QX_IN = QX_IN + QX(K,EL_IN,1  )*PHI_EDGE(K,I,LED)
          QY_IN = QY_IN + QY(K,EL_IN,1  )*PHI_EDGE(K,I,LED)

          HB_IN = HB_IN + HB(K,EL_IN,IRK)*PHI_EDGE(K,I,LED)
          HB_EX = HB_EX + HB0(K,EL_IN,1  )*PHI_EDGE(K,I,LED)
          
        ENDDO

C.....If sediment tranpsort due to waves is on compute the wave period,
C.....height, and angle at the given interior gauss point

        IF (SEDFLAG.EQ.2) THEN

          WAVE_RATIO = (XEGP(I) - (-1.D0))/2.D0

          T_WAVE = WAVE_T(N1) + WAVE_RATIO*(WAVE_T(N2) - WAVE_T(N1))
          A_WAVE = WAVE_A(N1) + WAVE_RATIO*(WAVE_A(N2) - WAVE_A(N1))
          H_WAVE = RAMPDG*(WAVE_H(N1)
     &                           + WAVE_RATIO*(WAVE_H(N2) - WAVE_H(N1)))

        ENDIF

C.....Set exterior variables equal to the interior variables

        ZE_EX = ZE_IN
        QX_EX = QX_IN
        QY_EX = QY_IN

        HB_EX = HB_IN

C.....Compute the total height of the water column

        HT_IN = ZE_IN*IFNLFA + HB_IN
        HT_EX = ZE_EX*IFNLFA + HB_EX

C.....Compute the velocities in the x and y directions

        U_IN = QX_IN/HT_IN
        U_EX = QX_EX/HT_EX

        V_IN = QY_IN/HT_IN
        V_EX = QY_EX/HT_EX
        
C.....Compute the magnitudes of the velocities

        UMAG_IN = SQRT(U_IN**(2.D0) + V_IN**(2.D0))
        UMAG_EX = SQRT(U_EX**(2.D0) + V_EX**(2.D0))

C.....Compute the Roe averaged velocities

        U_ROE = (U_IN*SQRT(HT_IN) + U_EX*SQRT(HT_EX))/
     &                                     (SQRT(HT_IN) + SQRT(HT_EX))
        V_ROE = (V_IN*SQRT(HT_IN) + V_EX*SQRT(HT_EX))/
     &                                     (SQRT(HT_IN) + SQRT(HT_EX))

C.....If the roe velocity vector and normal > 0 use the interior element
C.....values to compute the sediment transport

        IF (((U_ROE*NX + V_ROE*NY).GT.0).AND.(UMAG_IN.NE.0)) THEN

          CALL LUND_FORMULA(HT_IN, U_IN, V_IN, HB_IN, QSX, QSY)

C.....If the roe velocity vector and normal < 0 use the interior element
C.....values to compute the sediment transport

        ELSEIF (((U_ROE*NX + V_ROE*NY).LT.0).AND.(UMAG_EX.NE.0)) THEN

          CALL LUND_FORMULA(HT_EX, U_EX, V_EX, HB_EX, QSX, QSY)

C.....Else skip the sediment transport and edge integral calculations

        ELSE
          GOTO 111
        ENDIF

C.....Compute the sediment flux normal to the edge

        QS_HAT = QSX*NX + QSY*NY

C.....Compute the edge integral

        DO K = 1,DOF
          CALL EDGE_INT_SED(EL_IN,LED,GED,I,QS_HAT)
        ENDDO

111   CONTINUE

      RETURN
      END SUBROUTINE

