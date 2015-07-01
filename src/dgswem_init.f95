SUBROUTINE DGSWEM_INIT(s,dg_here,global_here,nodalattr_here,NT)
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
  INTEGER TESTFLAG,OUTITER,i,ModetoNode,time_here,ie,NT
  character*80 tecfile, tecfile_max
  
  CALL MAKE_DIRNAME(s)       ! Establish Working Directory Name
  CALL READ_INPUT(s,dg_here,global_here,nodalattr_here)         ! Establish sizes by reading fort.14 and fort.15

  
  
  !...  
  !...  ******************* START PROGRAM SETUP SECTION ********************
  !...  
  !     
      IF (global_here%IHOT.EQ.0) THEN
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

      
      CALL PREP_DG(s,dg_here,global_here,nodalattr_here)
      CALL WRITE_RESULTS(s,dg_here,global_here,0,.FALSE.)

!.....Write heading to unit 16
      WRITE(16,1112)
      WRITE(16,17931)
      IF (global_here%NSCREEN.EQ.1) WRITE(6,1112)
      IF (global_here%NSCREEN.EQ.1) WRITE(6,17931)

!     sb...Write initial conditions
      CALL WRITE_DG_IC(dg_here)

      NT = global_here%NT

 1112 FORMAT(/,1X,79('_'))
17931 FORMAT(//,1X,'LIMITED RUNTIME INFORMATION SECTION ',//)


    end SUBROUTINE DGSWEM_INIT
