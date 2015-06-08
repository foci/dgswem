C***********************************************************************
C
C     SUBROUTINE OCEAN_EDGE_HYDRO( )
C
C     This subroutine does the following:
C
C       1.  Calculates the values of the necessary variables at the edge
C           gauss points for ELEVATION SPECIFIED edges
C       2.  Calls the appropriate subroutine to compute the flux at
C           these points.
C       3.  Calls the appropriate subroutine to compute the boundary
C           integrals.
C
C     Written by Shintaro Bunya based on ocean_edge_hydro.F (10-11-2005) 
C
C***********************************************************************

      SUBROUTINE OCEAN_EDGE_HYDRO_STRICT( )

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE SED

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LED, GED, LN(3)
      REAL(SZ) ZP(3),ZAVG

      DO 1000 L=1, needs
      
C.....Retrieve the global and local edge number

      GED = NEEDN(L)
      LED = NEDSD(1,GED)

C.....Retrieve the elements which share the edge

      EL_IN = NEDEL(1,GED)
        
C.....Retrieve the components of the normal vector to the edge
        
      NX = COSNX(GED)
      NY = SINNX(GED)
      
      N1 = NEDNO(1,GED)
      N2 = NEDNO(2,GED)

C      print *,'X=',Y(N1),Y(N2),N1,N2

      LN(1) = MOD(LED+0,3) + 1
      LN(2) = MOD(LED+1,3) + 1
      LN(3) = MOD(LED+2,3) + 1

      ZAVG = ZE(1,EL_IN,IRK+1)
      
      ZP(:) = 0.D0

      DO I = 1,2

C.....Compute the specified open ocean elevation
        
        DO JJ=1,NBFR
        
          IF (PER(JJ).EQ.0.D0) THEN
            NCYC = 0.D0
          ELSE
            NCYC = INT(TIMEDG/PER(JJ))
          ENDIF

          ARGJ = AMIG(JJ)*(TIMEDG - NCYC*PER(JJ)) + FACE(JJ)
          RFF = FF(JJ)*RAMPDG
          
          EFA_GP = EFA_DG(JJ,L,I)
          EMO_GP = EMO_DG(JJ,L,I)

          ARG = ARGJ - EFA_GP
          
          ZP(LN(I)) = ZP(LN(I)) + EMO_GP*RFF*COS(ARG)

        ENDDO

      ENDDO

C      print*, 'ZP=',ZP(1),TIMEDG,ARG,EMO_GP,RFF,COS(0.810733248)
C      print*,'NM=',NM(EL_IN,LN1),N1,NM(EL_IN,LN2),N2

CC      ZP(LN(3)) = 3.D0*ZAVG - ZP(LN(1)) - ZP(LN(2))
C      ZP(LN(3)) = 0.5D0*(ZP(LN(1))+ZP(LN(2)))

C.....Compute elevations at node LN(3)

      DO KK = 1,DOF
        ZP(LN(3))  = ZP(LN(3)) + PHI_CORNER(KK,LN(3))*ZE(KK,EL_IN,IRK+1)
      ENDDO

      ZP(LN(3)) = 0.5D0*(ZP(LN(1))+ZP(LN(2)))

#if 0
      ZE(1,EL_IN,IRK+1) =
     &     1.D0/3.D0*(ZP(1)+ZP(2)+ZP(3))
      ZE(2,EL_IN,IRK+1) = 
     &     -1.D0/6.D0*(ZP(1)+ZP(2))+1.D0/3.D0*ZP(3)
      ZE(3,EL_IN,IRK+1) = 
     &     -0.5D0*ZP(1)+0.5D0*ZP(2)
#endif
      ZE(2:DOF,EL_IN,IRK+1) = 0.D0

#if 1
      QX(2:DOF,EL_IN,IRK+1) = 0.D0
      QY(2:DOF,EL_IN,IRK+1) = 0.D0
#endif

1000  CONTINUE
      RETURN
      END SUBROUTINE
