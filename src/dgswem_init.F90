SUBROUTINE DGSWEM_INIT(s,dg_here,global_here,nodalattr_here)
  USE SIZES
  USE GLOBAL
  USE DG
  USE NodalAttributes
  IMPLICIT NONE
     
  type (sizes_type) :: s
  type (dg_type) :: dg_here
  type (global_type) :: global_here
  type (nodalattr_type) :: nodalattr_here

  REAL(4) CPU_TIME,CPU_SEC(2)
  REAL(4) TARRAY(2)
  INTEGER TESTFLAG,OUTITER,i,ModetoNode,time_here,ie
  character*80 tecfile, tecfile_max  

!  print*, 'init fileunits'

  call init_fileunits(s)

  
!  print*, 'make_dirname'

  CALL MAKE_DIRNAME(s)       ! Establish Working Directory Name
!  print*, 'read_input'
  CALL READ_INPUT(s,dg_here,global_here,nodalattr_here)         ! Establish sizes by reading fort.14 and fort.15

  
  !...  
  !...  ******************* START PROGRAM SETUP SECTION ********************
  !...  
  !     
      IF (global_here%IHOT.EQ.0) THEN
!         print*, 'coldstart'
         CALL COLDSTART(s,global_here)
      ELSE
#ifdef HOTSTART
         CALL HOTSTART()
#else
         print*, "Hotstart not supported."
         stop
#endif
      ENDIF


!...  
!...  DETERMINE THE NUMBER OF ACTIVE ELEMENTS (global_here%MJU) and total number of 
!...  elements (global_here%NODELE) ATTACHED TO EACH NODE
!...  
      DO I=1,global_here%NP
         global_here%MJU(I)=0
         global_here%NODELE(I)=0
         global_here%NODECODE(I)=global_here%NNODECODE(I)
      END DO

      DO IE=1,global_here%NE
         global_here%IE=IE
         global_here%NM1=global_here%NM(global_here%IE,1)
         global_here%NM2=global_here%NM(global_here%IE,2)
         global_here%NM3=global_here%NM(global_here%IE,3)
         global_here%NCELE=global_here%NODECODE(global_here%NM1)*global_here%NODECODE(global_here%NM2)*global_here%NODECODE(global_here%NM3)
         global_here%MJU(global_here%NM1)=global_here%MJU(global_here%NM1)+global_here%NCELE
         global_here%MJU(global_here%NM2)=global_here%MJU(global_here%NM2)+global_here%NCELE
         global_here%MJU(global_here%NM3)=global_here%MJU(global_here%NM3)+global_here%NCELE
         global_here%NODELE(global_here%NM1)=global_here%NODELE(global_here%NM1)+1
         global_here%NODELE(global_here%NM2)=global_here%NODELE(global_here%NM2)+1
         global_here%NODELE(global_here%NM3)=global_here%NODELE(global_here%NM3)+1
      END DO

      DO I=1,global_here%NP
         IF(global_here%MJU(I).EQ.0) global_here%MJU(I)=1
      END DO
      
      
!...  
!...  ************* SET FLAGS AND COEFFICIENTS USED IN TIME STEPPING ***********
!...  

!...  NONLINEAR FLAGS

      IF(nodalattr_here%NOLIBF.EQ.0) THEN
         nodalattr_here%IFNLBF=0
         nodalattr_here%IFLINBF=1
         nodalattr_here%IFHYBF=0
      ENDIF
      IF(nodalattr_here%NOLIBF.EQ.1) THEN
         nodalattr_here%IFNLBF=1
         nodalattr_here%IFLINBF=0
         nodalattr_here%IFHYBF=0
      ENDIF
      IF(nodalattr_here%NOLIBF.EQ.2) THEN
         nodalattr_here%IFNLBF=0
         nodalattr_here%IFLINBF=0
         nodalattr_here%IFHYBF=1
      ENDIF
      IF(global_here%NOLIFA.EQ.0) THEN
         global_here%IFNLFA=0
      ELSE
         global_here%IFNLFA=1
      ENDIF
      IF(global_here%NOLICA.EQ.0) THEN
         global_here%IFNLCT=0
         global_here%NLEQ = 0.D0
         global_here%LEQ = 1.D0
      ELSE
         global_here%IFNLCT=1
         global_here%NLEQ = 1.D0
         global_here%LEQ = 0.D0
      ENDIF
      IF(global_here%NOLICAT.EQ.0) THEN
         global_here%IFNLCAT=0
         global_here%NLEQ = 0.D0
         global_here%LEQ = 1.D0
      ELSE
         global_here%IFNLCAT=1
         global_here%NLEQ = 1.D0
         global_here%LEQ = 0.D0
      ENDIF
      global_here%NLEQG = global_here%NLEQ*global_here%G
      dg_here%FG_L = global_here%LEQ*global_here%G

      global_here%IFWIND=1
      IF(global_here%IM.EQ.1) global_here%IFWIND=0
 
!...  CONSTANT COEFFICIENTS
!jj   w - version m10
!jj   w      global_here%TT0L=((1.0+0.5*global_here%DT*nodalattr_here%TAU0)/global_here%DT)/global_here%DT
      global_here%GA00=global_here%G*global_here%A00
!jj   w - version m10
!jj   w      global_here%TT0R=((0.5*nodalattr_here%TAU0*global_here%DT-1.0)/global_here%DT)/global_here%DT
      global_here%GC00=global_here%G*global_here%C00
      global_here%TADVODT=global_here%IFNLCAT/global_here%DT
      global_here%GB00A00=global_here%G*(global_here%B00+global_here%A00)
      global_here%GFAO2=global_here%G*global_here%IFNLFA/2.D0
      global_here%GO3=global_here%G/3.D0
      global_here%DTO2=global_here%DT/2.D0
      global_here%DT2=global_here%DT*2.D0
      global_here%GDTO2=global_here%G*global_here%DT/2.D0
      global_here%SADVDTO3=global_here%IFNLCT*global_here%DT/3.D0

      
!      print*, 'prep_dg'
      CALL PREP_DG(s,dg_here,global_here,nodalattr_here)
!      print*, 'write_results'
      CALL WRITE_RESULTS(s,dg_here,global_here,0,.FALSE.)

!.....Write heading to unit 16
      WRITE(16,1112)
      WRITE(16,17931)
      IF (global_here%NSCREEN.EQ.1) WRITE(6,1112)
      IF (global_here%NSCREEN.EQ.1) WRITE(6,17931)

!     sb...Write initial conditions
!      print*, 'write_dg_ic'
      CALL WRITE_DG_IC(dg_here)

 1112 FORMAT(/,1X,79('_'))
17931 FORMAT(//,1X,'LIMITED RUNTIME INFORMATION SECTION ',//)


    end SUBROUTINE DGSWEM_INIT

      SUBROUTINE NEIGHB(s,NE,NP,NM,NNEIGH,NEIGH,NEIMIN,NEIMAX,X,Y,NSCREEN,NNEIGH_ELEM,NEIGH_ELEM)
      USE SIZES
      type (sizes_type) :: s
!     
      INTEGER NP,NE,NEIMIN,NEIMAX,NSCREEN,N,NN,EN1,EN2,EN3,I,J,JJ,JLOW
      INTEGER :: NEIGH(S%MNP,S%MNEI),NNEIGH(S%MNP),NNEIGH_ELEM(S%MNP)
      INTEGER NM(S%MNE,3),NEIGH_ELEM(S%MNP,S%MNEI)
      REAL(8) X(S%MNP),Y(S%MNP),DELX,DELY,DIST
      REAL(8) ANGLELOW,ANGLEMORE,RAD2DEG
      REAL(8), ALLOCATABLE :: ANGLE(:)
      INTEGER,ALLOCATABLE :: NEITEM(:)

      integer :: myproc_here
      myproc_here = 0
!     
      ALLOCATE ( ANGLE(S%MNEI) )
      ALLOCATE ( NEITEM(S%MNP) )
!     
      RAD2DEG=45.0d0/ATAN(1.0d0)
!     

      print*, "NP=", np, "s%MNP=", s%MNP
      DO N=1,NP
         NNEIGH(N) = 0
         NNEIGH_ELEM(N) = 0
         DO NN=1,S%MNEI
            NEIGH(N,NN) = 0
            NEIGH_ELEM(N,NN) = 0
         ENDDO
      ENDDO
!     
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
            IF((NNEIGH(EN1).GT.S%MNEI-1).OR.(NNEIGH(EN2).GT.S%MNEI-1)) GOTO 999
            NEIGH(EN1,NNEIGH(EN1))=EN2
            NEIGH(EN2,NNEIGH(EN2))=EN1
 25         DO 30 J=1,NNEIGH(EN1)
 30            IF(EN3.EQ.NEIGH(EN1,J)) GOTO 35
               NNEIGH(EN1)=NNEIGH(EN1)+1
               NNEIGH(EN3)=NNEIGH(EN3)+1
               IF((NNEIGH(EN1).GT.S%MNEI-1).OR.(NNEIGH(EN3).GT.S%MNEI-1)) GOTO 999
               NEIGH(EN1,NNEIGH(EN1))=EN3
               NEIGH(EN3,NNEIGH(EN3))=EN1
 35            DO 50 J=1,NNEIGH(EN2)
 50               IF(EN3.EQ.NEIGH(EN2,J)) GOTO 10
                  NNEIGH(EN2)=NNEIGH(EN2)+1
                  NNEIGH(EN3)=NNEIGH(EN3)+1
                  IF((NNEIGH(EN2).GT.S%MNEI-1).OR.(NNEIGH(EN3).GT.S%MNEI-1)) GOTO 999
                  NEIGH(EN2,NNEIGH(EN2))=EN3
                  NEIGH(EN3,NNEIGH(EN3))=EN2
 10            CONTINUE
!     
!     INSERT NODE ITSELF IN PLACE #1 and SORT other NEIGHBORS by increasing cw angle from East
!     
               DO I=1,NP
                  DO J=1,NNEIGH(I)
                     NEITEM(J)=NEIGH(I,J)
                     DELX=X(NEITEM(J))-X(I)
                     DELY=Y(NEITEM(J))-Y(I)
                     DIST=SQRT(DELX*DELX+DELY*DELY)
                     IF(DIST.EQ.0.0d0) GOTO 998
                     IF(DELY.NE.0.0d0) THEN
                        ANGLE(J)=RAD2DEG*ACOS(DELX/DIST)
                        IF(DELY.GT.0.0) ANGLE(J)=360.0d0-ANGLE(J)
                     ENDIF
                     IF(DELY.EQ.0.0d0) THEN
                        IF(DELX.GT.0.0d0) ANGLE(J)=0.0d0
                        IF(DELX.LT.0.d0) ANGLE(J)=180.0d0
                     ENDIF
                  END DO
                  ANGLEMORE=-1.d0
                  DO JJ=1,NNEIGH(I)
                     ANGLELOW=400.d0
                     DO J=1,NNEIGH(I)
                        IF((ANGLE(J).LT.ANGLELOW).AND.(ANGLE(J).GT.ANGLEMORE))&
                            THEN
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

!     
!     DETERMINE THE MAXIMUM AND MINIMUM NUMBER OF NEIGHBORS
!     
               NEIMAX = 0
               NEIMIN = 1000
               DO 60 N=1,NP
                  IF(NNEIGH(N).LT.NEIMIN) NEIMIN=NNEIGH(N)
                  IF(NNEIGH(N).GT.NEIMAX) NEIMAX=NNEIGH(N)
 60            CONTINUE
!     
               RETURN

!     TERMINATE PROGRAM IF MAXIMUM NUMBER OF NEIGHBORS SET TOO SMALL

 999           CONTINUE
               IF(NSCREEN.EQ.1.AND.MYPROC_HERE.EQ.0) WRITE(6,99311)
               WRITE(16,99311)
99311          FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',              //,1X,'THE DIMENSIONING PARAMETER MNEI IS TOO SMALL'              /,1X,'USER MUST RE-DIMENSION PROGRAM',              //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
               STOP

 998           CONTINUE
               IF(NSCREEN.EQ.1.AND.MYPROC_HERE.EQ.0) WRITE(6,99312) I,NEITEM(J)
               WRITE(16,99312) I,NEITEM(J)
99312          FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',              //,1X,'NODES ',I7,' AND ',I7,         ' HAVE THE SAME COORDINATES'              //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)

               STOP
               END         
      


!******************************************************************************
!                                                                             *
!    Transform from lon,lat (lamda,phi) coordinates into CPP coordinates.     *
!    Lon,Lat must be in radians.                                              *
!                                                                             *
!******************************************************************************

      SUBROUTINE CPP(X,Y,RLAMBDA,PHI,RLAMBDA0,PHI0)
      REAL*8 X,Y,RLAMBDA,PHI,RLAMBDA0,PHI0,R
      R=6378206.4d0
      X=R*(RLAMBDA-RLAMBDA0)*COS(PHI0)
      Y=PHI*R
      RETURN
      END


!******************************************************************************
!                                                                             *
!    Transform from CPP coordinates to lon,lat (lamda,phi) coordinates        *
!    Lon,Lat is in radians.                                                   *
!                                                                             *
!******************************************************************************

      SUBROUTINE INVCP(XXCP,YYCP,RLAMBDA,PHI,RLAMBDA0,PHI0)
      REAL*8 XXCP,YYCP,RLAMBDA,PHI,RLAMBDA0,PHI0,R
      R=6378206.4d0
      RLAMBDA=RLAMBDA0+XXCP/(R*COS(PHI0))
      PHI=YYCP/R
      RETURN
      END
