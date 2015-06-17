!***********************************************************************
!
!     SUBROUTINE FLOW_EDGE_SED( )
!
!     This subroutine does the following:
!
!       1.  Calculates the values of the necessary variables at the edge
!           gauss points for NON-ZERO FLUX edges
!       2.  Calls the appropriate subroutine to compute the flux at
!           these points.
!       3.  Calls the appropriate subroutine to compute the boundary
!           integrals.
!
!     Written by Ethan Kubatko (03-05-2005)
!
!***********************************************************************

      SUBROUTINE FLOW_EDGE_SED(L)

!.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE SED

      IMPLICIT NONE

!.....Declare local variables

      INTEGER L, LED, GED,i,k

!.....Retrieve the global and local edge number

      GED = NFEDN(L)
      LED = NEDSD(1,GED)

!.....Retrieve the elements which share the edge

      EL_IN = NEDEL(1,GED)

!.....Retrieve the nodes of the given element

      N1 = NEDNO(1,GED)
      N2 = NEDNO(2,GED)

!.....Retrieve the components of the normal vector to the edge

      NX = COSNX(GED)
      NY = SINNX(GED)

!.....Compute ZE, QX, QY, and HB at each edge Gauss quadrature point

      DO 111 I = 1,NEGP

        ZE_IN = 0.D0
        QX_IN = 0.D0
        QY_IN = 0.D0

        HB_IN = 0.D0
        HB_EX = 0.D0

!.....Compute the solution at the interior state

        DO K = 1,DG%DOF

          ZE_IN = ZE_IN + ZE(K,EL_IN,1  )*PHI_EDGE(K,I,LED)
          QX_IN = QX_IN + QX(K,EL_IN,1  )*PHI_EDGE(K,I,LED)
          QY_IN = QY_IN + QY(K,EL_IN,1  )*PHI_EDGE(K,I,LED)

          HB_IN = HB_IN + HB(K,EL_IN,IRK)*PHI_EDGE(K,I,LED)
          HB_EX = HB_EX + HB0(K,EL_IN,1  )*PHI_EDGE(K,I,LED)
          
        ENDDO

!.....If sediment tranpsort due to waves is on compute the wave period,
!.....height, and angle at the given interior gauss point

        IF (SEDFLAG.EQ.2) THEN

          WAVE_RATIO = (XEGP(I) - (-1.D0))/2.D0

          T_WAVE = WAVE_T(N1) + WAVE_RATIO*(WAVE_T(N2) - WAVE_T(N1))
          A_WAVE = WAVE_A(N1) + WAVE_RATIO*(WAVE_A(N2) - WAVE_A(N1))
          H_WAVE = RAMPDG*(WAVE_H(N1)
     &                           + WAVE_RATIO*(WAVE_H(N2) - WAVE_H(N1)))

        ENDIF

!.....Set exterior variables equal to the interior variables

        ZE_EX = ZE_IN
        QX_EX = QX_IN
        QY_EX = QY_IN

        HB_EX = HB_IN

!.....Compute the total height of the water column

        HT_IN = ZE_IN*IFNLFA + HB_IN
        HT_EX = ZE_EX*IFNLFA + HB_EX

!.....Compute the velocities in the x and y directions

        U_IN = QX_IN/HT_IN
        U_EX = QX_EX/HT_EX

        V_IN = QY_IN/HT_IN
        V_EX = QY_EX/HT_EX
        
!.....Compute the magnitudes of the velocities

        UMAG_IN = SQRT(U_IN**(2.D0) + V_IN**(2.D0))
        UMAG_EX = SQRT(U_EX**(2.D0) + V_EX**(2.D0))

!.....Compute the Roe averaged velocities

        U_ROE = (U_IN*SQRT(HT_IN) + U_EX*SQRT(HT_EX))/
     &                                     (SQRT(HT_IN) + SQRT(HT_EX))
        V_ROE = (V_IN*SQRT(HT_IN) + V_EX*SQRT(HT_EX))/
     &                                     (SQRT(HT_IN) + SQRT(HT_EX))

!.....If the roe velocity vector and normal > 0 use the interior element
!.....values to compute the sediment transport

        IF (((U_ROE*NX + V_ROE*NY).GT.0).AND.(UMAG_IN.NE.0)) THEN

          CALL LUND_FORMULA(HT_IN, U_IN, V_IN, HB_IN, QSX, QSY)

!.....If the roe velocity vector and normal < 0 use the interior element
!.....values to compute the sediment transport

        ELSEIF (((U_ROE*NX + V_ROE*NY).LT.0).AND.(UMAG_EX.NE.0)) THEN

          CALL LUND_FORMULA(HT_EX, U_EX, V_EX, HB_EX, QSX, QSY)

!.....Else skip the sediment transport and edge integral calculations

        ELSE
          GOTO 111
        ENDIF

!.....Compute the sediment flux normal to the edge

        QS_HAT = QSX*NX + QSY*NY

!.....Compute the edge integral

        DO K = 1,DG%DOF
          CALL EDGE_INT_SED(EL_IN,LED,GED,I,QS_HAT)
        ENDDO

111   CONTINUE

      RETURN
      END SUBROUTINE

