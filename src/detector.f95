C***********************************************************************
C
C     SUBROUTINE DETECTOR()
C
C     Under Construction!!!
C
C***********************************************************************

      SUBROUTINE DETECTOR(L,IT)
      
C.....Use appropriate modules

      USE GLOBAL
      USE DG

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L,GED,LED_IN,LED_EX,II,IT,GP_IN(NEGP(ph)),GP_EX(NEGP(ph)),I,K,JJ
      REAL(SZ) DEN
      
C.....Set the inflow integrals and length to zero

      INFLOW_ZE = 0.D0
      
      INFLOW_LEN = 0.D0
      
C.....Loop over the three edges of the element
      
      DO II=1,3
      
C.....Retrieve the global edge number

        GED = NELED(II,L)

         PA = PDG_EL(GED)

#ifdef P0
         if (pa.eq.0) then
            pa = 1
         endif
#endif        

C.....Retrieve the edge length

        H_LEN = XLEN(GED)

C.....Retrieve the components of the normal vector to the edge

        EL_IN = NEDEL(1,GED)
        
        IF (EL_IN.EQ.L) THEN
          EL_EX  = NEDEL(2,GED)
          LED_IN = NEDSD(1,GED)
          LED_EX = NEDSD(2,GED)
          NX = COSNX(GED)
          NY = SINNX(GED)
          DO I=1,NEGP(pa)
            GP_IN(I) = I
            GP_EX(I) = NEGP(pa)-I+1
          ENDDO
        ELSE
          EL_EX  = NEDEL(1,GED)
          LED_IN = NEDSD(2,GED)
          LED_EX = NEDSD(1,GED)
          NX = -COSNX(GED)
          NY = -SINNX(GED)
          DO I=1,NEGP(pa)
            GP_IN(I) = NEGP(pa)-I+1
            GP_EX(I) = I
          ENDDO
        ENDIF

C.....Compute the solution at the midpoint of the edge

        QX_IN = 0.D0
        QY_IN = 0.D0
      
        DO K = 1,DOFS(L)

          QX_IN = QX_IN + QX(K,L,IRK+1)*PHI_MID(K,LED_IN,pa)
          QY_IN = QY_IN + QY(K,L,IRK+1)*PHI_MID(K,LED_IN,pa)
          
        ENDDO
         
C.....If the edge is an inflow boundary compute the edge integral

        IF ((QX_IN*NX + QY_IN*NY).LT.0) THEN
        
C.....Sum of the length of the outflow boundaries for the given element

          INFLOW_LEN = INFLOW_LEN + XLEN(GED)

          DO I = 1,NEGP(pa)

            ZE_IN = 0.D0
            ZE_EX = 0.D0

            DO K = 1,DOF
              ZE_IN = ZE_IN + ZE(K,L,IRK+1)*PHI_EDGE(K,GP_IN(I),LED_IN,pa)
            ENDDO

C.....If the edge is an open ocean edge compute the specified elevation

!            IF (OCEAN_EDGE_FLAG(GED)) THEN
            
              N1 = NEDNO(1,GED)
              N2 = NEDNO(2,GED)
            
              DO JJ=1,NBFR

                IF (PER(JJ).EQ.0.D0) THEN
                  NCYC = 0.D0
                ELSE
                  NCYC = INT(TIMEDG/PER(JJ))
                ENDIF

                ARGJ = AMIG(JJ)*(TIMEDG - NCYC*PER(JJ)) + FACE(JJ)
                RFF = FF(JJ)*RAMPDG
                EFA_GP = 0.5D0*(EFA(JJ,N1) + EFA(JJ,N2))
     &                 + 0.5D0*(EFA(JJ,N2) - EFA(JJ,N1))*XEGP(I,pa)
                EMO_GP = 0.5D0*(EMO(JJ,N1) + EMO(JJ,N2))
     &                 + 0.5D0*(EMO(JJ,N2) - EMO(JJ,N1))*XEGP(I,pa)
                ARG = ARGJ - EFA_GP
                ZE_EX = ZE_EX + EMO_GP*RFF*COS(ARG)

              ENDDO

C.....If the edge is a land edge compute the exterior value
              
!            ELSEIF (LAND_EDGE_FLAG(GED)) THEN
            
              ZE_EX = ZE_IN

!            ELSE
              DO K = 1,DOFS(L)
                ZE_EX = ZE_EX
     &                   + ZE(K,EL_EX,IRK+1)*PHI_EDGE(K,GP_EX(I),LED_EX,pa)
              ENDDO
!            ENDIF
            
C.....Sum the edge integral

            INFLOW_ZE = INFLOW_ZE +
     &                       XLEN(GED)/2.D0*(ZE_IN-ZE_EX)*WEGP(GP_IN(I),pa)
            
C.....Store the minimum edge length of the given element

            H_LEN = MIN(H_LEN,XLEN(GED))
            
          ENDDO
          
        ENDIF
        
      ENDDO
      
C.....If the given element has inflow boundaries estimate the max norm

      IF (INFLOW_ZE.NE.0) THEN
      
        ZE_NORM = 0.D0
     
        DO I=1,NAGP(ph)
        
          ZE_IN = 0.D0

          DO K = 1,DOFS(L)

            ZE_IN = ZE_IN + ZE(K,L,IRK)*PHI_AREA(K,I,pa)

          ENDDO
          
          ZE_NORM = MAX(ZE_NORM,ZE_IN)
            
        ENDDO
          
      ENDIF
      
C.....Compute the discontinuity detector

      DEN = H_LEN**((pa+1)/2.D0)*INFLOW_LEN*ZE_NORM
      IF (DEN.NE.0) THEN
        ZE_DECT= ABS(INFLOW_ZE)/DEN
      ELSE
        ZE_DECT = 0.D0
      ENDIF
      
      IF (ZE_DECT.GT.1) THEN
        !WRITE(89,*) L,ZE_DECT,TIME,IT
        DO K=2,DOFS(L)
          ZE(K,L,IRK+1) = 0
          QX(K,L,IRK+1) = 0
          QY(K,L,IRK+1) = 0
        ENDDO
      ENDIF
      
      IF (ZE_DECT.GT.1) THEN
        WRITE(87,*) 'LIMITING ON ELEMENT #',L
      ENDIF
      
      RETURN
      END SUBROUTINE

      
      
