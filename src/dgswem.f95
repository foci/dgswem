!******************************************************************************
!     DGSWEM RELEASE VERSION 1 11/2013   
!     
!******************************************************************************

      PROGRAM DGSWEM
!     
!     sb-disabled for Linux&Intel compiler
#ifdef VF
      USE DFPORT
#endif
      USE GLOBAL
#ifdef HARM
      USE HARM
#endif
      USE DG
      USE NodalAttributes, ONLY :&
          NoLiBF, NWP, Tau0, HBreak, FTheta, FGamma, Tau, CF, IFNLBF,&
          InitNAModule, ReadNodalAttr, InitNodalAttr, ESLM, ESLC,&
          IFLINBF, IFHYBF
#ifdef CMPI
      USE MESSENGER_ELEM
!--   
      IMPLICIT NONE

      type (sizes_type) :: s
      type (dg_type) :: dg_here
      type (global_type) :: global_here
      
      INTEGER TESTFLAG,OUTITER,istop,i,ModetoNode
!     sb-PDG1
      REAL(4) CPU_TIME,CPU_SEC(2)
      REAL(4) TARRAY(2)
!     nd
      character*80 tecfile,tecfile_max

      CALL MESSAGE_INIT()       ! Init MPI and get MPI-rank of this cpu

                                ! <ezpp-begin>
      IF(MYPROC.EQ.0)THEN
         write(*,*) 'mnproc ',mnproc
      ENDIF
      CALL MAKE_DIRNAME(s)       ! Establish Working Directory Name
      CALL READ_INPUT(s,dg_here)         ! Establish sizes by reading fort.14 and fort.15
#else
      IMPLICIT NONE
     
      type (sizes_type) :: s
      type (dg_type) :: dg_here

      REAL(4) CPU_TIME,CPU_SEC(2)
      REAL(4) TARRAY(2)
      INTEGER TESTFLAG,OUTITER,i,ModetoNode
      character*80 tecfile, tecfile_max

      CALL MAKE_DIRNAME(s)       ! Establish Working Directory Name
      CALL READ_INPUT(s,dg_here)         ! Establish sizes by reading fort.14 and fort.15
#endif         

      
      
!...  
!...  ******************* START PROGRAM SETUP SECTION ****************************
!...  
!     
      IF (global_here%IHOT.EQ.0) THEN
         CALL COLDSTART(s)
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

      DO global_here%IE=1,global_here%NE
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
!jj   w      global_here%TT0L=((1.0+0.5*global_here%DT*TAU0)/global_here%DT)/global_here%DT
      global_here%GA00=global_here%G*global_here%A00
!jj   w - version m10
!jj   w      global_here%TT0R=((0.5*TAU0*global_here%DT-1.0)/global_here%DT)/global_here%DT
      global_here%GC00=global_here%G*global_here%C00
      global_here%TADVODT=global_here%IFNLCAT/global_here%DT
      global_here%GB00A00=global_here%G*(global_here%B00+global_here%A00)
      global_here%GFAO2=global_here%G*global_here%IFNLFA/2.D0
      global_here%GO3=global_here%G/3.D0
      global_here%DTO2=global_here%DT/2.D0
      global_here%DT2=global_here%DT*2.D0
      global_here%GDTO2=global_here%G*global_here%DT/2.D0
      global_here%SADVDTO3=global_here%IFNLCT*global_here%DT/3.D0

      
!*************************DG SWEM*******************************

!     write(200+myproc,*) 'call prep_dg'
!     write(200+myproc,*) 'back from prep_dg'

      CALL PREP_DG(s,dg_here,global_here)
      CALL WRITE_RESULTS(s,dg_here,global_here,0,.FALSE.)

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

                                !call dg%iwrite(tecfile,14,19,myproc)

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
     $        'VARIABLES = "global_here%x", "global_here%y", "b", "dg_here%ze", "H", "u", "v", "|v|","|w|","p","dg_here%iota","dg_here%iota2","sum","diff"'
         write(777,*) 'ZONE ZONETYPE=FETRIANGLE ',
     $        'NODES=', global_here%np, 
     $        ' ELEMENTS=', global_here%ne, 
!     $     ' DATAPACKING=POINT ','SOLUTIONTIME=',global_here%time_a
     $        ' DATAPACKING=BLOCK ',' VARLOCATION=([3,4,5,6,7,8,10,11,12,13,14]=CELLCENTERED) ',
     $        'SOLUTIONTIME=',global_here%time_a 

!$$$  write(777,7777) global_here%x(i), global_here%y(i),  
!$$$  $        global_here%dp(i), global_here%eta2(i), global_here%eta2(i)+global_here%dp(i),global_here%uu2(i),global_here%vv2(i),
!$$$  $        sqrt(global_here%uu2(i)**2+global_here%vv2(i)**2),sqrt(global_here%wsx2(i)**2+global_here%wsy2(i)**2),
!$$$  $        myproc,global_here%tracer(i),global_here%tracer2(i),abs(global_here%tracer(i)+global_here%tracer2(i)),
!$$$  $        abs(global_here%tracer(i)-global_here%tracer2(i))
!$$$  enddo

         do i=1,global_here%np
            write(777,7777)  global_here%x(i)
         enddo
         do i=1,global_here%np
            write(777,7777)  global_here%y(i)
         enddo
         do i=1,global_here%ne
            write(777,7777)  global_here%dpe(i)
         enddo
         do i=1,global_here%ne 
            write(777,7777)  global_here%eta2(i)
         enddo
         do i=1,global_here%ne 
            write(777,7777)  global_here%eta2(i)+global_here%dpe(i)
         enddo
         do i=1,global_here%ne 
            write(777,7777)  global_here%uu2(i)
         enddo
         do i=1,global_here%ne 
            write(777,7777)  global_here%vv2(i)
         enddo
         do i=1,global_here%ne 
            write(777,7777)  sqrt(global_here%uu2(i)**2+global_here%vv2(i)**2)
         enddo
         do i=1,global_here%np 
            write(777,7777)  sqrt(global_here%wsx2(i)**2+global_here%wsy2(i)**2)
         enddo
         do i=1,global_here%ne 
            write(777,7777)  DBLE(global_here%pdg_el(i))
         enddo
         do i=1,global_here%ne 
            write(777,7777)  global_here%tracer(i)
         enddo
         do i=1,global_here%ne 
            write(777,7777)  global_here%tracer2(i)
         enddo
         do i=1,global_here%ne 
            write(777,7777)  abs(global_here%tracer(i)+global_here%tracer2(i))
         enddo
         do i=1,global_here%ne 
            write(777,7777)  abs(global_here%tracer(i)-global_here%tracer2(i))    
         enddo
      else
         write(777,*) 
     $'VARIABLES = "global_here%x", "global_here%y", "b", "dg_here%ze", "H", "u", "v", "|v|","|w|","PE"'
         write(777,*) 'ZONE ZONETYPE=FETRIANGLE ',
     $        'NODES=', global_here%np, 
     $        ' ELEMENTS=', global_here%ne, 
     $        ' DATAPACKING=POINT ','SOLUTIONTIME=',global_here%time_a
         do i=1,global_here%np
            if (global_here%ics.eq.2) then
               write(777,7777) global_here%slam(i)/deg2rad, global_here%sfea(i)/deg2rad,  
     $              global_here%dp(i), global_here%eta2(i), global_here%eta2(i)+global_here%dp(i),global_here%uu2(i),global_here%vv2(i),
     $              sqrt(global_here%uu2(i)**2+global_here%vv2(i)**2),sqrt(global_here%wsx2(i)**2+global_here%wsy2(i)**2),
     $              myproc
            else
               write(777,7777) global_here%x(i), global_here%y(i),
     $              global_here%dp(i), global_here%eta2(i), global_here%eta2(i)+global_here%dp(i),global_here%uu2(i),global_here%vv2(i),
     $              sqrt(global_here%uu2(i)**2+global_here%vv2(i)**2),sqrt(global_here%wsx2(i)**2+global_here%wsy2(i)**2),
     $              myproc
            endif
         enddo

      endif

 7777    format(9f20.8,i10)
         do i=1,global_here%ne
            write(777,"(3i12)") global_here%nm(i,1), global_here%nm(i,2), global_here%nm(i,3)
         enddo

!Casey 120813: End the OUT_TEC conditional.
#endif

#ifdef SWAN
!asey 101118: Allow SWAN to initialize stuff before the start
!             of the time step loop.  This subroutine is inside
!             the 'couple2swan.dg%F' src file.
      CALL PADCSWAN_INIT
#endif

!.....Write heading to unit 16
      WRITE(16,1112)
      WRITE(16,17931)
#ifdef CMPI
      IF (global_here%NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,1112)
      IF (global_here%NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,17931)
#else      
      IF (global_here%NSCREEN.EQ.1) WRITE(6,1112)
      IF (global_here%NSCREEN.EQ.1) WRITE(6,17931)
#endif

!     sb...Write initial conditions
      CALL WRITE_DG_IC(dg_here)
#ifdef CMPI
!     istop=1
!     if (istop.eq.1) then
!     call message_fini()
!     stop
!     endif
#endif
!.....Begin time stepping
      DO 200 global_here%ITIME_A = global_here%ITHS+1,global_here%NT
!$$$         if (mod(global_here%itime_a,1000).eq.1) then
!$$$            if (myproc.eq.0) write(*,*) 'timestep ',global_here%itime_a
!$$$c     write(200+myproc,*) 'timestep ',global_here%itime_a,myproc
!$$$         endif
         CALL DG_TIMESTEP(s,dg_here,global_here%ITIME_A)
#ifdef SWAN
!asey 090302: If it is time, then call the following subroutine
!             to then call the SWAN time-stepping subroutine.
         IF(MOD(global_here%ITIME_A,CouplingInterval).EQ.0)THEN
           CALL PADCSWAN_RUN(global_here%ITIME_A)
!asey 121126: DEBUG.
!          IF(ALLOCATED(PASS2SWAN))THEN
!             DO I=1,global_here%NP
!                DO J=1,PASS2SWAN(I)%global_here%NO_NBORS
!                   PASS2SWAN(I)%global_here%ETA1(J) = PASS2SWAN(I)%global_here%ETA2(J)
!                ENDDO
!             ENDDO
!           ENDIF
         ENDIF
#endif
 200  CONTINUE

 1112 FORMAT(/,1X,79('_'))
17931 FORMAT(//,1X,'LIMITED RUNTIME INFORMATION SECTION ',//)
      

!...  
!...  ****************** TIME STEPPING LOOP ENDS HERE ********************
!...  
!...  
!...  IF global_here%IHARIND=1 SOLVE THE HARMONIC ANALYSIS PROBLEM AND WRITE OUTPUT
!...  

#ifdef HARM
      IF ((global_here%IHARIND.EQ.1).AND.(global_here%ITIME_A.GT.global_here%ITHAS)) THEN

!...  LINES COMPUTE MEANS AND VARIANCES
!...  FOR CHECKING THE HARMONIC ANALYSIS RESULTS.
!...  ACCUMULATE VARIANCE AND MEAN OF RECORD AT NODES.
!     
         if (CHARMV) then
            IF (global_here%FMV.global_here%NE.0.) THEN
               DO I=1,global_here%NP
                  global_here%ELAV(I)   = global_here%ELAV(I)/global_here%NTSTEPS
                  global_here%XVELAV(I) = global_here%XVELAV(I)/global_here%NTSTEPS
                  global_here%YVELAV(I) = global_here%YVELAV(I)/global_here%NTSTEPS
                  global_here%ELVA(I)   = global_here%ELVA(I)/global_here%NTSTEPS   - global_here%ELAV(I)*global_here%ELAV(I)
                  global_here%XVELVA(I) = global_here%XVELVA(I)/global_here%NTSTEPS - global_here%XVELAV(I)*global_here%XVELAV(I)
                  global_here%YVELVA(I) = global_here%YVELVA(I)/global_here%NTSTEPS - global_here%YVELAV(I)*global_here%YVELAV(I)
               END DO
               global_here%TIMEBEG=global_here%ITMV*global_here%DTDP + (global_here%STATIM-global_here%REFTIM)*86400.D0
               OPEN(55,FILE=DIRNAME//'/'//'fort.55')
               WRITE(55,*) global_here%NP
            ENDIF
         endif                  !  charmv
!...  
!......Fill out and decompose the LHS harmonic analaysis matrix
!...  
         CALL FULSOL(0)
!...  
!......Solve the harmonic analysis problem and write the output
!...  
         IF(global_here%NHAGE.EQ.1) CALL LSQSOLEG(global_here%NP,DIRNAME,LNAME,global_here%ELAV,global_here%ELVA)
!     
         IF(global_here%NHAGV.EQ.1) CALL LSQSOLVG(global_here%NP,DIRNAME,LNAME,global_here%XVELAV,global_here%YVELAV,global_here%XVELVA,global_here%YVELVA)
!     
         IF(global_here%NHASE.EQ.1) CALL LSQSOLES(global_here%NSTAE,DIRNAME,LNAME)
!     
         IF(global_here%NHASV.EQ.1) CALL LSQSOLVS(global_here%NSTAV,DIRNAME,LNAME)
!     
      ENDIF
#endif      
#ifdef SWAN
!asey 101118: Let SWAN clean up stuff.
      CALL PADCSWAN_FINAL
#endif

#ifdef CMPI
      CALL MESSAGE_FINI()
#endif
!     
!     
#ifdef VF
!     CALL ETIME(TARRAY)
!     CPU_TIME = TARRAY(1) + TARRAY(2)
!     PRINT*,'CPU_TIME = ',CPU_TIME
#endif

      CALL ETIME(TARRAY,CPU_TIME)
!      CPU_TIME = TARRAY(1) + TARRAY(2)
      PRINT*,'CPU_TIME = ',CPU_TIME
      STOP
   END  ! I think this is the end of the program? -zdb

!******************************************************************************
!     *
!     Subroutine to generate a neighbor table from a connectivity table.     *
!     *
!     NOTE:the node itself is listed as neighbor #1                          *
!     NOTE:all other neighbors are sorted and placed in cw order from east   *
!     *
!     R.L.       4/26/95                                    *
!******************************************************************************
!     *
!     -  PARAMETERS WHICH MUST BE SET TO CONTROL THE DIMENSIONING OF ARRAYS   *
!     ARE AS FOLLOWS:                                                   *
!     *
!     MNP = MAXIMUM NUMBER OF NODAL POINTS                               *
!     MNE = MAXIMUM NUMBER OF ELEMENTS                                   *
!     MNEI= 1+MAXIMUM NUMBER OF NODES CONNECTED TO ANY ONE NODE IN THE   *
!     FINITE ELEMENT GRID                                       *
!     *
!******************************************************************************
!     *
!     VARIABLE DEFINITIONS:                                                 *
!     NE - NUMBER OF ELEMENTS                                               *
!     NP - NUMBER OF NODES                                                  *
!     NM(MNE,3) - NODE NUMBERS ASSOCIATED WITH EACH ELEMENT                 *
!     NNEIGH(MNP) NUMBER OF NEIGHBORS FOR EACH NODE                         *
!     NNEIGH_ELEM(MNP) NUMBER OF NEIGHBORS BY ELEMENT                       *
!     NEIGH(MNP,NEIMAX) 2D ARRAY OF NEIGHBORS FOR EACH NODE                 *
!     NEIMIN - 1+MINIMUM NUMBER OF NEIGHBORS FOR ANY NODE                   *
!     NEIMAX - 1+MAXIMUM NUMBER OF NEIGHBORS FOR ANY NODE                   *
!     *
!******************************************************************************
!     
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
