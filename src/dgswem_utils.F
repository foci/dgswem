      SUBROUTINE NEIGHB(NE,NP,NM,NNEIGH,NEIGH,NEIMIN,NEIMAX,
     &     X,Y,NSCREEN,NNEIGH_ELEM,NEIGH_ELEM)
      USE SIZES
C     
      INTEGER NP,NE,NEIMIN,NEIMAX,NSCREEN,N,NN,EN1,EN2,EN3,I,J
      INTEGER :: NEIGH(MNP,MNEI),NNEIGH(MNP),NNEIGH_ELEM(MNP)
      INTEGER NM(MNE,3),NEIGH_ELEM(MNP,MNEI)
      real(sz) X(MNP),Y(MNP),DELX,DELY,DIST
      real(sz) ANGLELOW,ANGLEMORE,RAD2DEG
      real(sz), ALLOCATABLE :: ANGLE(:)
      INTEGER,ALLOCATABLE :: NEITEM(:)
C     
      ALLOCATE ( ANGLE(MNEI) )
      ALLOCATE ( NEITEM(MNP) )
C     
      RAD2DEG=45.0/ATAN(1.0)
C     
      DO N=1,NP
         NNEIGH(N) = 0
         NNEIGH_ELEM(N) = 0
         DO NN=1,MNEI
            NEIGH(N,NN) = 0
            NEIGH_ELEM(N,NN) = 0
         ENDDO
      ENDDO
C     
      DO 10 N=1,NE
         EN1 = NM(N,1)
         EN2 = NM(N,2)
         EN3 = NM(N,3)
         NNEIGH_ELEM(EN1) = NNEIGH_ELEM(EN1)+1
         NEIGH_ELEM(EN1,NNEIGH_ELEM(EN1)) = N
         NNEIGH_ELEM(EN2) = NNEIGH_ELEM(EN2)+1
         NEIGH_ELEM(EN2,NNEIGH_ELEM(EN2)) = N
         NNEIGH_ELEM(EN3) = NNEIGH_ELEM(EN3)+1
         NEIGH_ELEM(EN3,NNEIGH_ELEM(EN3)) = N
         DO 20 J=1,NNEIGH(EN1)
 20         IF(EN2.EQ.NEIGH(EN1,J)) GOTO 25
            NNEIGH(EN1)=NNEIGH(EN1)+1
            NNEIGH(EN2)=NNEIGH(EN2)+1
         IF((NNEIGH(EN1).GT.MNEI-1).OR.(NNEIGH(EN2).GT.MNEI-1)) GOTO 999
            NEIGH(EN1,NNEIGH(EN1))=EN2
            NEIGH(EN2,NNEIGH(EN2))=EN1
 25         DO 30 J=1,NNEIGH(EN1)
 30            IF(EN3.EQ.NEIGH(EN1,J)) GOTO 35
               NNEIGH(EN1)=NNEIGH(EN1)+1
               NNEIGH(EN3)=NNEIGH(EN3)+1
         IF((NNEIGH(EN1).GT.MNEI-1).OR.(NNEIGH(EN3).GT.MNEI-1)) GOTO 999
               NEIGH(EN1,NNEIGH(EN1))=EN3
               NEIGH(EN3,NNEIGH(EN3))=EN1
 35            DO 50 J=1,NNEIGH(EN2)
 50               IF(EN3.EQ.NEIGH(EN2,J)) GOTO 10
                  NNEIGH(EN2)=NNEIGH(EN2)+1
                  NNEIGH(EN3)=NNEIGH(EN3)+1
         IF((NNEIGH(EN2).GT.MNEI-1).OR.(NNEIGH(EN3).GT.MNEI-1)) GOTO 999
                  NEIGH(EN2,NNEIGH(EN2))=EN3
                  NEIGH(EN3,NNEIGH(EN3))=EN2
 10            CONTINUE
C     
C     INSERT NODE ITSELF IN PLACE #1 and SORT other NEIGHBORS by increasing cw angle from East
C     
               DO I=1,NP
                  DO J=1,NNEIGH(I)
                     NEITEM(J)=NEIGH(I,J)
                     DELX=X(NEITEM(J))-X(I)
                     DELY=Y(NEITEM(J))-Y(I)
                     DIST=SQRT(DELX*DELX+DELY*DELY)
                     IF(DIST.EQ.0.0) GOTO 998
                     IF(DELY.NE.0.0) THEN
                        ANGLE(J)=RAD2DEG*ACOS(DELX/DIST)
                        IF(DELY.GT.0.0) ANGLE(J)=360.0-ANGLE(J)
                     ENDIF
                     IF(DELY.EQ.0.0) THEN
                        IF(DELX.GT.0.0) ANGLE(J)=0.0
                        IF(DELX.LT.0.0) ANGLE(J)=180.0
                     ENDIF
                  END DO
                  ANGLEMORE=-1.0
                  DO JJ=1,NNEIGH(I)
                     ANGLELOW=400.0
                     DO J=1,NNEIGH(I)
                 IF((ANGLE(J).LT.ANGLELOW).AND.(ANGLE(J).GT.ANGLEMORE))
     &                       THEN
                           ANGLELOW=ANGLE(J)
                           JLOW=J
                        ENDIF
                     END DO
                     NEIGH(I,JJ+1)=NEITEM(JLOW)
                     ANGLEMORE=ANGLELOW
                  END DO
                  NEIGH(I,1)=I
                  NNEIGH(I)=NNEIGH(I)+1
               ENDDO

C     
C     DETERMINE THE MAXIMUM AND MINIMUM NUMBER OF NEIGHBORS
C     
               NEIMAX = 0
               NEIMIN = 1000
               DO 60 N=1,NP
                  IF(NNEIGH(N).LT.NEIMIN) NEIMIN=NNEIGH(N)
                  IF(NNEIGH(N).GT.NEIMAX) NEIMAX=NNEIGH(N)
 60            CONTINUE
C     
               RETURN

C     TERMINATE PROGRAM IF MAXIMUM NUMBER OF NEIGHBORS SET TOO SMALL

 999           CONTINUE
               IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,99311)
               WRITE(16,99311)
99311      FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &           //,1X,'THE DIMENSIONING PARAMETER MNEI IS TOO SMALL'
     &           /,1X,'USER MUST RE-DIMENSION PROGRAM',
     &        //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
               STOP

 998           CONTINUE
             IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,99312) I,NEITEM(J)
               WRITE(16,99312) I,NEITEM(J)
99312      FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &              //,1X,'NODES ',I7,' AND ',I7,
     &         ' HAVE THE SAME COORDINATES'
     &       //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
               STOP
               END         
      


C******************************************************************************
C                                                                             *
C    Transform from lon,lat (lamda,phi) coordinates into CPP coordinates.     *
C    Lon,Lat must be in radians.                                              *
C                                                                             *
C******************************************************************************

      SUBROUTINE CPP(X,Y,RLAMBDA,PHI,RLAMBDA0,PHI0)
      use sizes, only : sz
      real(sz) X,Y,RLAMBDA,PHI,RLAMBDA0,PHI0,R
      R=6378206.4d0
      X=R*(RLAMBDA-RLAMBDA0)*COS(PHI0)
      Y=PHI*R
      RETURN
      END


C******************************************************************************
C                                                                             *
C    Transform from CPP coordinates to lon,lat (lamda,phi) coordinates        *
C    Lon,Lat is in radians.                                                   *
C                                                                             *
C******************************************************************************

      SUBROUTINE INVCP(XXCP,YYCP,RLAMBDA,PHI,RLAMBDA0,PHI0)
      use sizes, only : sz
      real(sz) XXCP,YYCP,RLAMBDA,PHI,RLAMBDA0,PHI0,R
      R=6378206.4d0
      RLAMBDA=RLAMBDA0+XXCP/(R*COS(PHI0))
      PHI=YYCP/R
      RETURN
      END
