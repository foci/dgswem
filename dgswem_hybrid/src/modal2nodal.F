      SUBROUTINE Modal2Nodal
#ifdef SWAN
Casey 110422: Create new data structure with discontinuous information for SWAN.
      USE Couple2Adcirc, ONLY: PASS2SWAN
#endif
      USE DG,     ONLY: DOF,
     &                  HB,
     &                  ph,
     &                  PHI_CORNER,
     &                  QX,
     &                  QY,
     &                  WDFLG,
     &                  ZE
      USE GLOBAL, ONLY: AREA_SUM,
     &                  AREAS,
     &                  CEN_SUM,
     &                  DP,
     &                  EL_COUNT,
     &                  ELETAB,
     &                  ETA2,
     &                  HB_DG,
     &                  LEQ,
     &                  MNP,
     &                  N1,
     &                  NBOR_EL,
     &                  NLEQ,
     &                  NO_NBORS,
     &                  NM,
     &                  QX_DG,
     &                  QY_DG,
     &                  SEDFLAG,
     &                  UU2,
     &                  VV2,
     &                  ZE_DG
      USE SIZES,  ONLY: SZ
      IMPLICIT NONE
      INTEGER  :: I
      INTEGER  :: J
      INTEGER  :: K
      INTEGER  :: KK
      REAL(SZ) :: ANGLE_SUM
      REAL(SZ) :: AREA
      REAL(SZ) :: DEPTH
      REAL(SZ) :: FH_NL
#ifdef SWAN
Casey 110111: Create new data structure with discontinuous information for SWAN.
      IF(.NOT.ALLOCATED(PASS2SWAN)) ALLOCATE(PASS2SWAN(1:MNP))
#endif
      DO I = 1,MNP
        NO_NBORS = EL_COUNT(I)
#ifdef SWAN
Casey 110111: Create new data structure with discontinuous information for SWAN.
        PASS2SWAN(I)%NO_NBORS = NO_NBORS
        IF(.NOT.ALLOCATED(PASS2SWAN(I)%NBOR_EL))THEN
           ALLOCATE(PASS2SWAN(I)%NBOR_EL(1:NO_NBORS))
           PASS2SWAN(I)%NBOR_EL = 0
        ENDIF
        IF(.NOT.ALLOCATED(PASS2SWAN(I)%ETA1))THEN
           ALLOCATE(PASS2SWAN(I)%ETA1(1:NO_NBORS))
           PASS2SWAN(I)%ETA1 = 0.
        ENDIF
        IF(.NOT.ALLOCATED(PASS2SWAN(I)%ETA2))THEN
           ALLOCATE(PASS2SWAN(I)%ETA2(1:NO_NBORS))
           PASS2SWAN(I)%ETA2 = 0.
        ENDIF
        IF(.NOT.ALLOCATED(PASS2SWAN(I)%UU2))THEN
           ALLOCATE(PASS2SWAN(I)%UU2(1:NO_NBORS))
           PASS2SWAN(I)%UU2 = 0.
        ENDIF
        IF(.NOT.ALLOCATED(PASS2SWAN(I)%VV2))THEN
           ALLOCATE(PASS2SWAN(I)%VV2(1:NO_NBORS))
           PASS2SWAN(I)%VV2 = 0.
        ENDIF
        IF(.NOT.ALLOCATED(PASS2SWAN(I)%VELMAG))THEN
           ALLOCATE(PASS2SWAN(I)%VELMAG(1:NO_NBORS))
           PASS2SWAN(I)%VELMAG = 0.
        ENDIF
#endif
        AREA_SUM = 0
        ANGLE_SUM = 0
        CEN_SUM = 0
        DO 333 J = 1,NO_NBORS
          NBOR_EL = ELETAB(I,1+J)
#ifdef SWAN
Casey 110111: Create new data structure with discontinuous information for SWAN.
          PASS2SWAN(I)%NBOR_EL(J) = NBOR_EL
#endif
          IF(WDFLG(NBOR_EL).EQ.0) CYCLE
          DO K = 1,3
            N1 = NM(NBOR_EL,K)
            IF (N1.EQ.I) THEN
              ZE_DG(J) = ZE(1,NBOR_EL,1)
              QX_DG(J) = QX(1,NBOR_EL,1)
              QY_DG(J) = QY(1,NBOR_EL,1)
              HB_DG(J) = HB(1,NBOR_EL,1)
              DO KK = 2,DOF
                ZE_DG(J) = ZE_DG(J) + PHI_CORNER(KK,K,ph)*ZE(KK,NBOR_EL,1)
                QX_DG(J) = QX_DG(J) + PHI_CORNER(KK,K,ph)*QX(KK,NBOR_EL,1)
                QY_DG(J) = QY_DG(J) + PHI_CORNER(KK,K,ph)*QY(KK,NBOR_EL,1)
                HB_DG(J) = HB_DG(J) + PHI_CORNER(KK,K,ph)*HB(KK,NBOR_EL,1)
              ENDDO
              AREA = 0.5*AREAS(NBOR_EL)
              AREA_SUM = AREA_SUM + AREA
              GOTO 333
            ENDIF
          ENDDO
333     CONTINUE
        ETA2(I) = 0.D0
        UU2(I)  = 0.D0
        VV2(I)  = 0.D0
        DO J = 1,NO_NBORS
          NBOR_EL = ELETAB(I,1+J)
          IF(WDFLG(NBOR_EL).EQ.0) CYCLE
          AREA = 0.5*AREAS(NBOR_EL)/AREA_SUM
          DEPTH = ZE_DG(J) + HB_DG(J)
          FH_NL = 1.D0/(NLEQ*DEPTH + LEQ)
          ETA2(I) = ETA2(I) + AREA*ZE_DG(J)
          UU2(I)  = UU2(I)  + AREA*QX_DG(J)*FH_NL
          VV2(I)  = VV2(I)  + AREA*QY_DG(J)*FH_NL
#ifdef SWAN
Casey 110111: Create new data structure with discontinuous information for SWAN.
!         PASS2SWAN(I)%ETA1(J) = PASS2SWAN(I)%ETA2(J)
          PASS2SWAN(I)%ETA2(J) = REAL( ZE_DG(J) )
          PASS2SWAN(I)%UU2(J)  = REAL( QX_DG(J)*FH_NL )
          PASS2SWAN(I)%VV2(J)  = REAL( QY_DG(J)*FH_NL )
          PASS2SWAN(I)%VELMAG(J)  = REAL( SQRT( PASS2SWAN(I)%UU2(J) *
     &             PASS2SWAN(I)%UU2(J) + PASS2SWAN(I)%VV2(J) *
     &             PASS2SWAN(I)%VV2(J) ) )
#endif
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE
