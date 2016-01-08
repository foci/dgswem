C******************************************************************************
C     DGSWEM RELEASE VERSION 1 11/2013   
C     
C******************************************************************************

      PROGRAM DGSWEM
C     
C     sb-disabled for Linux&Intel compiler
#ifdef VF
      USE DFPORT
#endif
      USE GLOBAL
      USE HARM
      USE DG
      USE NodalAttributes, ONLY :
     &     NoLiBF, NWP, Tau0, HBreak, FTheta, FGamma, Tau, CF, IFNLBF,
     &     InitNAModule, ReadNodalAttr, InitNodalAttr, ESLM, ESLC,
     &     IFLINBF, IFHYBF
#ifdef SWAN
Casey 110422: Create a new data structure with discontinuous information.
Casey 121126: DEBUG.
!     USE Couple2Adcirc, ONLY: PASS2SWAN
Casey 101118: We need the following information to couple to unstructured SWAN.
      USE Couple2Swan, ONLY: CouplingInterval,
     &                       PADCSWAN_FINAL,
     &                       PADCSWAN_INIT,
     &                       PADCSWAN_RUN
#endif
#ifdef CMPI
      USE MESSENGER_ELEM
C--   
      IMPLICIT NONE
      
      INTEGER TESTFLAG,OUTITER,istop,i,ModetoNode
C     sb-PDG1
      REAL(4) CPU_TIME,CPU_SEC(2)
      REAL(4) TARRAY(2)
c     nd
      character*80 tecfile,tecfile_max

      CALL MESSAGE_INIT()       ! Init MPI and get MPI-rank of this cpu

                                ! <ezpp-begin>
      IF(MYPROC.EQ.0)THEN
         write(*,*) 'mnproc ',mnproc
      ENDIF
      CALL MAKE_DIRNAME()       ! Establish Working Directory Name
      CALL READ_INPUT()         ! Establish sizes by reading fort.14 and fort.15
#else
      IMPLICIT NONE
      

      REAL(4) CPU_TIME,CPU_SEC(2)
      REAL(4) TARRAY(2)
      INTEGER TESTFLAG,OUTITER,i,ModetoNode
      character*80 tecfile, tecfile_max

      CALL MAKE_DIRNAME()       ! Establish Working Directory Name
      CALL READ_INPUT()         ! Establish sizes by reading fort.14 and fort.15
#endif         

      
      
C...  
C...  ******************* START PROGRAM SETUP SECTION ****************************
C...  
C     
      IF (IHOT.EQ.0) THEN
         CALL COLDSTART()
      ELSE
         CALL HOTSTART()
      ENDIF


C...  
C...  DETERMINE THE NUMBER OF ACTIVE ELEMENTS (MJU) and total number of 
C...  elements (NODELE) ATTACHED TO EACH NODE
C...  
      DO I=1,NP
         MJU(I)=0
         NODELE(I)=0
         NODECODE(I)=NNODECODE(I)
      END DO

      DO IE=1,NE
         NM1=NM(IE,1)
         NM2=NM(IE,2)
         NM3=NM(IE,3)
         NCELE=NODECODE(NM1)*NODECODE(NM2)*NODECODE(NM3)
         MJU(NM1)=MJU(NM1)+NCELE
         MJU(NM2)=MJU(NM2)+NCELE
         MJU(NM3)=MJU(NM3)+NCELE
         NODELE(NM1)=NODELE(NM1)+1
         NODELE(NM2)=NODELE(NM2)+1
         NODELE(NM3)=NODELE(NM3)+1
      END DO

      DO I=1,NP
         IF(MJU(I).EQ.0) MJU(I)=1
      END DO
      
      
C...  
C...  ************* SET FLAGS AND COEFFICIENTS USED IN TIME STEPPING ***********
C...  

C...  NONLINEAR FLAGS

      IF(NOLIBF.EQ.0) THEN
         IFNLBF=0
         IFLINBF=1
         IFHYBF=0
      ENDIF
      IF(NOLIBF.EQ.1) THEN
         IFNLBF=1
         IFLINBF=0
         IFHYBF=0
      ENDIF
      IF(NOLIBF.EQ.2) THEN
         IFNLBF=0
         IFLINBF=0
         IFHYBF=1
      ENDIF
      IF(NOLIFA.EQ.0) THEN
         IFNLFA=0
      ELSE
         IFNLFA=1
      ENDIF
      IF(NOLICA.EQ.0) THEN
         IFNLCT=0
         NLEQ = 0.D0
         LEQ = 1.D0
      ELSE
         IFNLCT=1
         NLEQ = 1.D0
         LEQ = 0.D0
      ENDIF
      IF(NOLICAT.EQ.0) THEN
         IFNLCAT=0
         NLEQ = 0.D0
         LEQ = 1.D0
      ELSE
         IFNLCAT=1
         NLEQ = 1.D0
         LEQ = 0.D0
      ENDIF
      NLEQG = NLEQ*G
      FG_L = LEQ*G

      IFWIND=1
      IF(IM.EQ.1) IFWIND=0
 
C...  CONSTANT COEFFICIENTS
cjj   w - version m10
cjj   w      TT0L=((1.0+0.5*DT*TAU0)/DT)/DT
      GA00=G*A00
cjj   w - version m10
cjj   w      TT0R=((0.5*TAU0*DT-1.0)/DT)/DT
      GC00=G*C00
      TADVODT=IFNLCAT/DT
      GB00A00=G*(B00+A00)
      GFAO2=G*IFNLFA/2.D0
      GO3=G/3.D0
      DTO2=DT/2.D0
      DT2=DT*2.D0
      GDTO2=G*DT/2.D0
      SADVDTO3=IFNLCT*DT/3.D0

      
C*************************DG SWEM*******************************

c     write(200+myproc,*) 'call prep_dg'
c     write(200+myproc,*) 'back from prep_dg'

      CALL PREP_DG()
      CALL WRITE_RESULTS(0,.FALSE.)

                               !cnd...for tecplot output
!Casey 120813: Begin the OUT_TEC conditional.
#ifdef OUT_TEC

      ModetoNode = 0

#ifdef CMPI
      !CALL SYSTEM('mkdir'//'/'//'tecplot_output{00..'//'4'//'}')
      tecfile = DIRNAME//'/'//'..'//'/'//'tecplot_output'//
     &     '/'//'tecplot.out_pe00000'
      tecfile_max = DIRNAME//'/'//'..'//'/'//'tecplot_max'//
     &     '/'//'tecplot.max_pe00000'

                                !call iwrite(tecfile,14,19,myproc)

      write(tecfile(45:49), "(i4.4)") myproc
      write(tecfile_max(42:46), "(i4.4)") myproc
#else
      tecfile = "tecplot.out"
      tecfile_max = "tecplot.max"
#endif

      open(777, file = tecfile)
      open(778, file = tecfile_max)
      write(777,*)  'TITLE = "dgswem output"'

      if (ModetoNode.eq.0) then
         write(777,*) 
     $        'VARIABLES = "x", "y", "b", "ze", "H", "u", "v", "|v|","|w|","p","iota","iota2","sum","diff"'
         write(777,*) 'ZONE ZONETYPE=FETRIANGLE ',
     $        'NODES=', np, 
     $        ' ELEMENTS=', ne, 
!     $     ' DATAPACKING=POINT ','SOLUTIONTIME=',time_a
     $        ' DATAPACKING=BLOCK ',' VARLOCATION=([3,4,5,6,7,8,10,11,12,13,14]=CELLCENTERED) ',
     $        'SOLUTIONTIME=',time_a 

c$$$  write(777,7777) x(i), y(i),  
c$$$  $        dp(i), eta2(i), eta2(i)+dp(i),uu2(i),vv2(i),
c$$$  $        sqrt(uu2(i)**2+vv2(i)**2),sqrt(wsx2(i)**2+wsy2(i)**2),
c$$$  $        myproc,tracer(i),tracer2(i),abs(tracer(i)+tracer2(i)),
c$$$  $        abs(tracer(i)-tracer2(i))
c$$$  enddo

         do i=1,np
            write(777,7777)  x(i)
         enddo
         do i=1,np
            write(777,7777)  y(i)
         enddo
         do i=1,ne
            write(777,7777)  dpe(i)
         enddo
         do i=1,ne 
            write(777,7777)  eta2(i)
         enddo
         do i=1,ne 
            write(777,7777)  eta2(i)+dpe(i)
         enddo
         do i=1,ne 
            write(777,7777)  uu2(i)
         enddo
         do i=1,ne 
            write(777,7777)  vv2(i)
         enddo
         do i=1,ne 
            write(777,7777)  sqrt(uu2(i)**2+vv2(i)**2)
         enddo
         do i=1,np 
            write(777,7777)  sqrt(wsx2(i)**2+wsy2(i)**2)
         enddo
         do i=1,ne 
            write(777,7777)  DBLE(pdg_el(i))
         enddo
         do i=1,ne 
            write(777,7777)  tracer(i)
         enddo
         do i=1,ne 
            write(777,7777)  tracer2(i)
         enddo
         do i=1,ne 
            write(777,7777)  abs(tracer(i)+tracer2(i))
         enddo
         do i=1,ne 
            write(777,7777)  abs(tracer(i)-tracer2(i))    
         enddo
      else
         write(777,*) 
     $'VARIABLES = "x", "y", "b", "ze", "H", "u", "v", "|v|","|w|","PE"'
         write(777,*) 'ZONE ZONETYPE=FETRIANGLE ',
     $        'NODES=', np, 
     $        ' ELEMENTS=', ne, 
     $        ' DATAPACKING=POINT ','SOLUTIONTIME=',time_a
         do i=1,np
            if (ics.eq.2) then
               write(777,7777) slam(i)/deg2rad, sfea(i)/deg2rad,  
     $              dp(i), eta2(i), eta2(i)+dp(i),uu2(i),vv2(i),
     $              sqrt(uu2(i)**2+vv2(i)**2),sqrt(wsx2(i)**2+wsy2(i)**2),
     $              myproc
            else
               write(777,7777) x(i), y(i),
     $              dp(i), eta2(i), eta2(i)+dp(i),uu2(i),vv2(i),
     $              sqrt(uu2(i)**2+vv2(i)**2),sqrt(wsx2(i)**2+wsy2(i)**2),
     $              myproc
            endif
         enddo

      endif

 7777    format(9f20.8,i10)
         do i=1,ne
            write(777,"(3i12)") nm(i,1), nm(i,2), nm(i,3)
         enddo

!Casey 120813: End the OUT_TEC conditional.
#endif

#ifdef SWAN
Casey 101118: Allow SWAN to initialize stuff before the start
C             of the time step loop.  This subroutine is inside
C             the 'couple2swan.F' src file.
      CALL PADCSWAN_INIT
#endif

C.....Write heading to unit 16
      WRITE(16,1112)
      WRITE(16,17931)
      IF (NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,1112)
      IF (NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,17931)

C     sb...Write initial conditions
      CALL WRITE_DG_IC()
#ifdef CMPI
!     istop=1
!     if (istop.eq.1) then
!     call message_fini()
!     stop
!     endif
#endif
C.....Begin time stepping
      DO 200 ITIME_A = ITHS+1,NT
c$$$         if (mod(itime_a,1000).eq.1) then
c$$$            if (myproc.eq.0) write(*,*) 'timestep ',itime_a
c$$$c     write(200+myproc,*) 'timestep ',itime_a,myproc
c$$$         endif
         CALL DG_TIMESTEP(ITIME_A)
#ifdef SWAN
Casey 090302: If it is time, then call the following subroutine
C             to then call the SWAN time-stepping subroutine.
         IF(MOD(ITIME_A,CouplingInterval).EQ.0)THEN
           CALL PADCSWAN_RUN(ITIME_A)
Casey 121126: DEBUG.
!          IF(ALLOCATED(PASS2SWAN))THEN
!             DO I=1,NP
!                DO J=1,PASS2SWAN(I)%NO_NBORS
!                   PASS2SWAN(I)%ETA1(J) = PASS2SWAN(I)%ETA2(J)
!                ENDDO
!             ENDDO
!           ENDIF
         ENDIF
#endif
 200  CONTINUE

 1112 FORMAT(/,1X,79('_'))
17931 FORMAT(//,1X,'LIMITED RUNTIME INFORMATION SECTION ',//)
      

C...  
C...  ****************** TIME STEPPING LOOP ENDS HERE ********************
C...  
C...  
C...  IF IHARIND=1 SOLVE THE HARMONIC ANALYSIS PROBLEM AND WRITE OUTPUT
C...  
      IF ((IHARIND.EQ.1).AND.(ITIME_A.GT.ITHAS)) THEN

C...  LINES COMPUTE MEANS AND VARIANCES
C...  FOR CHECKING THE HARMONIC ANALYSIS RESULTS.
C...  ACCUMULATE VARIANCE AND MEAN OF RECORD AT NODES.
C     
         if (CHARMV) then
            IF (FMV.NE.0.) THEN
               DO I=1,NP
                  ELAV(I)   = ELAV(I)/NTSTEPS
                  XVELAV(I) = XVELAV(I)/NTSTEPS
                  YVELAV(I) = YVELAV(I)/NTSTEPS
                  ELVA(I)   = ELVA(I)/NTSTEPS   - ELAV(I)*ELAV(I)
                  XVELVA(I) = XVELVA(I)/NTSTEPS - XVELAV(I)*XVELAV(I)
                  YVELVA(I) = YVELVA(I)/NTSTEPS - YVELAV(I)*YVELAV(I)
               END DO
               TIMEBEG=ITMV*DTDP + (STATIM-REFTIM)*86400.D0
               OPEN(55,FILE=DIRNAME//'/'//'fort.55')
               WRITE(55,*) NP
            ENDIF
         endif                  !  charmv
C...  
C......Fill out and decompose the LHS harmonic analaysis matrix
C...  
         CALL FULSOL(0)
C...  
C......Solve the harmonic analysis problem and write the output
C...  
         IF(NHAGE.EQ.1) CALL LSQSOLEG(NP,DIRNAME,LNAME,ELAV,ELVA)
C     
         IF(NHAGV.EQ.1) CALL LSQSOLVG(NP,DIRNAME,LNAME,
     &        XVELAV,YVELAV,XVELVA,YVELVA)
C     
         IF(NHASE.EQ.1) CALL LSQSOLES(NSTAE,DIRNAME,LNAME)
C     
         IF(NHASV.EQ.1) CALL LSQSOLVS(NSTAV,DIRNAME,LNAME)
C     
      ENDIF
      
#ifdef SWAN
Casey 101118: Let SWAN clean up stuff.
      CALL PADCSWAN_FINAL
#endif

#ifdef CMPI
      CALL MESSAGE_FINI()
#endif
C     
C     
#ifdef VF
C     CALL ETIME(TARRAY)
C     CPU_TIME = TARRAY(1) + TARRAY(2)
C     PRINT*,'CPU_TIME = ',CPU_TIME
#endif

      CALL ETIME(TARRAY)
      CPU_TIME = TARRAY(1) + TARRAY(2)
      PRINT*,'CPU_TIME = ',CPU_TIME
      STOP
      END

C******************************************************************************
C     *
C     Subroutine to generate a neighbor table from a connectivity table.     *
c     *
c     NOTE:the node itself is listed as neighbor #1                          *
c     NOTE:all other neighbors are sorted and placed in cw order from east   *
c     *
c     R.L.       4/26/95                                    *
C******************************************************************************
C     *
C     -  PARAMETERS WHICH MUST BE SET TO CONTROL THE DIMENSIONING OF ARRAYS   *
C     ARE AS FOLLOWS:                                                   *
C     *
C     MNP = MAXIMUM NUMBER OF NODAL POINTS                               *
C     MNE = MAXIMUM NUMBER OF ELEMENTS                                   *
C     MNEI= 1+MAXIMUM NUMBER OF NODES CONNECTED TO ANY ONE NODE IN THE   *
C     FINITE ELEMENT GRID                                       *
C     *
C******************************************************************************
C     *
C     VARIABLE DEFINITIONS:                                                 *
C     NE - NUMBER OF ELEMENTS                                               *
C     NP - NUMBER OF NODES                                                  *
C     NM(MNE,3) - NODE NUMBERS ASSOCIATED WITH EACH ELEMENT                 *
C     NNEIGH(MNP) NUMBER OF NEIGHBORS FOR EACH NODE                         *
C     NNEIGH_ELEM(MNP) NUMBER OF NEIGHBORS BY ELEMENT                       *
C     NEIGH(MNP,NEIMAX) 2D ARRAY OF NEIGHBORS FOR EACH NODE                 *
C     NEIMIN - 1+MINIMUM NUMBER OF NEIGHBORS FOR ANY NODE                   *
C     NEIMAX - 1+MAXIMUM NUMBER OF NEIGHBORS FOR ANY NODE                   *
C     *
C******************************************************************************
C     
      SUBROUTINE NEIGHB(NE,NP,NM,NNEIGH,NEIGH,NEIMIN,NEIMAX,
     &     X,Y,NSCREEN,NNEIGH_ELEM,NEIGH_ELEM)
      USE SIZES
C     
      INTEGER NP,NE,NEIMIN,NEIMAX,NSCREEN,N,NN,EN1,EN2,EN3,I,J
      INTEGER :: NEIGH(MNP,MNEI),NNEIGH(MNP),NNEIGH_ELEM(MNP)
      INTEGER NM(MNE,3),NEIGH_ELEM(MNP,MNEI)
      REAL(8) X(MNP),Y(MNP),DELX,DELY,DIST
      REAL(8) ANGLELOW,ANGLEMORE,RAD2DEG
      REAL(8), ALLOCATABLE :: ANGLE(:)
      INTEGER,ALLOCATABLE :: NEITEM(:)
C     
      ALLOCATE ( ANGLE(MNEI) )
      ALLOCATE ( NEITEM(MNP) )
C     
      RAD2DEG=45.0d0/ATAN(1.0d0)
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
99311          FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &              //,1X,'THE DIMENSIONING PARAMETER MNEI IS TOO SMALL'
     &              /,1X,'USER MUST RE-DIMENSION PROGRAM',
     &              //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
               STOP

 998           CONTINUE
               IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,99312) I,NEITEM(J)
               WRITE(16,99312) I,NEITEM(J)
99312          FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &              //,1X,'NODES ',I7,' AND ',I7,
     &         ' HAVE THE SAME COORDINATES'
     &              //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)

               STOP
               END         
      


C******************************************************************************
C                                                                             *
C    Transform from lon,lat (lamda,phi) coordinates into CPP coordinates.     *
C    Lon,Lat must be in radians.                                              *
C                                                                             *
C******************************************************************************

      SUBROUTINE CPP(X,Y,RLAMBDA,PHI,RLAMBDA0,PHI0)
      REAL*8 X,Y,RLAMBDA,PHI,RLAMBDA0,PHI0,R
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
      REAL*8 XXCP,YYCP,RLAMBDA,PHI,RLAMBDA0,PHI0,R
      R=6378206.4d0
      RLAMBDA=RLAMBDA0+XXCP/(R*COS(PHI0))
      PHI=YYCP/R
      RETURN
      END
