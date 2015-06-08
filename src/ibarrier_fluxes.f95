C***********************************************************************
C
C     SUBROUTINE IBARRIER_FLUXES( )
C
C     This subroutine computes the inward/outward normal flow over
C     specified internal barrier boundary (permeable or not) nodes.
C     Note that the averaged nodal values of the surface elevation are
C     used.
C
C     (taken from timestep.f)
C
C     Modifications, comments, and general clean-up
C     by Ethan Kubatko (06-02-2005)
C
C     Modifications for wet/dry
C     by Ethan Kubatko (09-02-2005)
C
C***********************************************************************

      SUBROUTINE IBARRIER_FLUXES()
      
C.....Use appropriate modules

      USE GLOBAL
      USE DG

      IMPLICIT NONE
      
      INTEGER i,k
      REAL(SZ) BRAMP

      BRAMP = 1.0
      
C.....Zero out the node code

      NIBNODECODE(:) = 0
        
C.....Loop over the flow boundary nodes
        
      DO 1034 I = 1,NVEL
        
C.....If the node is part of an internal barrier compute the flux based
C.....on weir formulas
        
        IF ((LBCODEI(I).EQ.4).OR.(LBCODEI(I).EQ.24).OR.
     &      (LBCODEI(I).EQ.5).OR.(LBCODEI(I).EQ.25)) THEN
     
C.....Retrieve the global node number for the flow boundary node

          NNBB1 = NBV(I)
            
C.....Retrieve the global node number on opposite side of the barrier

          NNBB2 = IBCONN(I)
            
C.....Compute the water level on each side of the barrier
            
          IF (IBSTART.EQ.0) THEN
            RBARWL1AVG(I) = ETA2(NNBB1) - BARINHT(I)
            RBARWL2AVG(I) = ETA2(NNBB2) - BARINHT(I)
            IBSTART = 1
          ELSE
            RBARWL1AVG(I) = (ETA2(NNBB1) - BARINHT(I) + BARAVGWT
     &                                    *RBARWL1AVG(I))/(1 + BARAVGWT)
            RBARWL2AVG(I) = (ETA2(NNBB2) - BARINHT(I) + BARAVGWT
     &                                    *RBARWL2AVG(I))/(1 + BARAVGWT)
          ENDIF
          RBARWL1 = RBARWL1AVG(I)
          RBARWL2 = RBARWL2AVG(I)
          RBARWL1F = 2.0D0*RBARWL1/3.0D0
          RBARWL2F = 2.0D0*RBARWL2/3.0D0
          QIB(NNBB1) = 0.0D0
            
C.....Case 0:  The side with a higher Water level is dry 
C.....I want to use WDFLG(EL) here in the future. sb

          IF (RBARWL1.GT.RBARWL2) THEN
            IF((ETA2(NNBB1)+DP(NNBB1)).LT.H0) THEN
              QIB(NNBB1) = 0.0D0
              GOTO 1034
            ENDIF
          ELSE 
            IF ((ETA2(NNBB2)+DP(NNBB2)).LT.H0) THEN
              QIB(NNBB1) = 0.0D0
              GOTO 1034
            ENDIF
          ENDIF

C.....Case 1:  Water level on both sides of barrier are below barrier

Csb-debug
#if 0
          IF ((RBARWL1.LT.0.0).AND.(RBARWL2.LT.0.0)) THEN
#else
          IF (.TRUE.) THEN
#endif
            QIB(NNBB1) = 0.0D0
            GOTO 1034
          ENDIF

C.....Case 2:  Water level equal on both sides of the barrier

          IF(RBARWL1.EQ.RBARWL2) THEN
            QIB(NNBB1) = 0.0D0
            GOTO 1034
          ENDIF
              
C.....Cases 3 and 4:  Overtopping occurring and water level HIGHER on
C.....this side of barrier
              
          IF ((RBARWL1.GT.RBARWL2).AND.(RBARWL1.GT.BARMIN)) THEN
              
C.....Case 3:  Outward subcritical flow

            IF (RBARWL2.GT.RBARWL1F) THEN
              QIB(NNBB1) = -BRAMP*RBARWL2*BARINCFSB(I)*
     &                                 (2.D0*G*(RBARWL1-RBARWL2))**0.5D0
              NIBNODECODE(NNBB1) = 1

C.....Case 4:  Outward supercritical flow

            ELSE
              QIB(NNBB1) = -BRAMP*BARINCFSP(I)*RBARWL1F*(RBARWL1F*G)
     &                                                           **0.5D0
              NIBNODECODE(NNBB1) = 1
            ENDIF

            GOTO 1034

          ENDIF
              
C.....Cases 5 and 6:  Overtopping occurring and water level LOWER on
C.....this side of the barrier
              
          IF ((RBARWL2.GT.RBARWL1).AND.(RBARWL2.GT.BARMIN)) THEN

C.....Case 5:  Inward subcritical

            IF (RBARWL1.GT.RBARWL2F) THEN
              QIB(NNBB1) = BRAMP*RBARWL1*BARINCFSB(I)*
     &                                (2.0D0*G*(RBARWL2-RBARWL1))**0.5D0
              NIBNODECODE(NNBB1) = 1

C.....Case 6:  Inward supercritical flow

            ELSE
              QIB(NNBB1) = BRAMP*BARINCFSP(I)*RBARWL2F*(RBARWL2F*G)
     &                                                           **0.5D0
              NIBNODECODE(NNBB1) = 1
            ENDIF
                
            GOTO 1034
                
          ENDIF
        ENDIF
1034  CONTINUE
      
      RETURN
      END SUBROUTINE
