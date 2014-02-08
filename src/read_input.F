C***********************************************************************
C     
C     SUBROUTINE READ_INPUT()
C     
C     Surprisingly enough this subroutine reads in the input
C     
C     Written by (too) many people
C     
C***********************************************************************

      SUBROUTINE READ_INPUT()

C.....Use appropriate modules

      USE GLOBAL
      USE HARM
      USE WIND
      USE DG
      USE NodalAttributes, ONLY :
     &     NoLiBF, NWP, Tau0, HBreak, FTheta, FGamma, Tau, CF,
     &     InitNAModule, ReadNodalAttr, InitNodalAttr, ESLM, ESLC,
     &     TAU0VAR, STARTDRY, FRIC, EVM, EVC
#ifdef CMPI
      USE MESSENGER_NEW
#endif
#ifdef SWAN
Casey 101118: Enable hot-start file generation by SWAN.
      USE Couple2Swan, ONLY: SwanHotStartUnit
#endif

      IMPLICIT NONE
      
C.....Declare local variables

      INTEGER NIBP, IBN1, IK, NDISC, NBBN, NVEL2, II, i,j,jj,k
      CHARACTER(256) LINE
      
C     sb-PDG1 added
#ifdef CMPI
      INTEGER IDUM80
      CHARACTER CDUM80
#endif
C--   
C     ek...Zero out all the variables in the Nodal Attributes Module
C     ek...Added from version 46

      CALL InitNAModule()

C     sb-PDG1 added
#ifdef CMPI
C     
C     When running in parallel check to make sure that the number of
C     processors (MNPROC - obtained from the job control script via a call
C     to the messenger module) is the same as that used in ADCPREP and written
C     in the fort.80 file
C     
      OPEN(80,FILE='fort.80')
      READ(80,'(A)') CDUM80     !Skip RUNDES
      READ(80,'(A)') CDUM80     !Skip RUNID
      READ(80,'(A)') CDUM80     !Skip AGRID
      READ(80,*) IDUM80         !Skip NELG & NNODG
      READ(80,*) IDUM80         !Read in NPROC
      CLOSE(80)
      IF(IDUM80.NE.MNPROC) THEN
         IF(MYPROC.EQ.0) THEN
            WRITE(*,'(A)') '*** ERROR IN PARALLEL SETUP!'
            WRITE(*,'(2A,I4,A)') '*** Number of CPUS for submitted job ',
     &           '(NCPU = ',MNPROC,') is not equal to the'
            WRITE(*,'(2A,I4,A)') '*** number of CPUS specified during',
     &           ' ADCPREP (see fort.80: NCPU = ',IDUM80,').'
            WRITE(*,'(A)') '*** dgswem will now quit!'
         ENDIF
         CALL MESSAGE_FINI()
         STOP
      ENDIF
#endif
C--   

      ScreenUnit = 6

C.....Initialize all runtime option logicals to false

      C2DDI  = .FALSE.
      C3D    = .FALSE.
      vertexslope = .False.
      C3DDSS = .FALSE.
      C3DVS  = .FALSE.
      CLUMP  = .FALSE.
      CTIP   = .FALSE.
      CHARMV = .FALSE.

C.....Open statement for unit 14, 15, and 25 (fort.dg) input files

      OPEN(14,FILE=DIRNAME//'/'//'fort.14')
      OPEN(15,FILE=DIRNAME//'/'//'fort.15')
      OPEN(17,FILE=DIRNAME//'/'//'fort.17')
      OPEN(25,FILE=DIRNAME//'/'//'fort.dg')

C.....Open statement for unit 16 output file
      
      OPEN(16,FILE=DIRNAME//'/'//'fort.16')

C.....General purpose format statements

 1112 FORMAT(/,1X,79('_'))
 9972 FORMAT(////,1X,'!!!!!!!!!! INPUT ERROR !!!!!!!!!',/)
 9973 FORMAT(/,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
 9974 FORMAT(/,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!!!',//)

C.....Read in the fort.dg file

      READ(25,*) DGSWE
      READ(25,*) padapt,pflag
      READ(25,*) gflag,diorism
      READ(25,*) pl,ph,px
      READ(25,*) slimit
      READ(25,*) plimit
      READ(25,*) pflag2con1,pflag2con2,lebesgueP 
      READ(25,*) FLUXTYPE
      READ(25,*) RK_STAGE, RK_ORDER
      READ(25,*) Lts_Flag, Mbar, MaxNdom, dt_min0, dtcheck, cflfac ! CJT
      READ(25,*) LTSmodflag
C.....LTS
      if (LTSmodflag.eq.1) then
         READ(25,*) nmodtimes
         do k=1,nmodtimes
            READ(25,*) tmod1(k),tmod2(k),TsGroupmod(k)
            tmod1(k)=tmod1(k)*86400.
            tmod2(k)=tmod2(k)*86400.
         enddo
      endif

      READ(25,*) DG_TO_CG
      READ(25,*) MODAL_IC
      READ(25,*) DGHOT, DGHOTSPOOL
      READ(25,"(A256)") LINE
      READ(LINE,*) SLOPEFLAG
      IF(SLOPEFLAG.EQ.2) THEN
         READ(LINE,*) SLOPEFLAG, SL2_M, SL2_NYU
      ENDIF
      IF(SLOPEFLAG.EQ.3) THEN
         READ(LINE,*) SLOPEFLAG, SL2_M, SL2_NYU, SL3_MD
      ENDIF
      IF(SLOPEFLAG.EQ.4) THEN
         READ(LINE,*) SLOPEFLAG,slope_weight
         vertexslope = .True.
      ENDIF
      IF(SLOPEFLAG.EQ.5) THEN
         READ(LINE,*) SLOPEFLAG
         vertexslope = .True.
      ENDIF
      IF(SLOPEFLAG.EQ.6) THEN
         READ(LINE,*) SLOPEFLAG,slope_weight
         vertexslope = .True.
      ENDIF
      IF(SLOPEFLAG.EQ.7) THEN
         READ(LINE,*) SLOPEFLAG,slope_weight
         vertexslope = .True.
      ENDIF
      IF(SLOPEFLAG.EQ.8) THEN
         READ(LINE,*) SLOPEFLAG,slope_weight
         vertexslope = .True.
      ENDIF
      IF(SLOPEFLAG.EQ.9) THEN
         READ(LINE,*) SLOPEFLAG,slope_weight
         vertexslope = .True.
      ENDIF
      IF(SLOPEFLAG.EQ.10) THEN
         READ(LINE,*) SLOPEFLAG,slope_weight
         vertexslope = .True.
      ENDIF
      READ(25,*) SEDFLAG,porosity,SEVDM,layers
      READ(25,*) reaction_rate
      READ(25,*) MNES
      READ(25,*) artdif,kappa,s0,uniform_dif,tune_by_hand
      READ(25,'(a)') sed_equationX
      READ(25,'(a)') sed_equationY

      RHOWAT0 = 1000.D0
      IF(FLUXTYPE.NE.1.AND.FLUXTYPE.NE.2.AND.FLUXTYPE.NE.3.AND.FLUXTYPE.NE.4) THEN
         PRINT *, 'SPECIFIED FLUXTYPE (=', FLUXTYPE,') IS NOT ALLOWED.'
         PRINT *, 'EXECUTION WILL BE TERMINATED.'
         STOP 'SPECIFIED FLUXTYPE IS NOT ALLOWED.'
      ENDIF

C.....Print out header for output including version number and copyright

      WRITE(16,1112)
      WRITE(16,1112)
      WRITE(16,1114)
      WRITE(16,1112)
      IF (MYPROC.EQ.0) THEN
         WRITE(6,1112)
         WRITE(6,1114)
         WRITE(6,1112)
      ENDIF

 1114 FORMAT(//,19X,'dgswem.11.13 ',
     &     //,5X,'A disconintuous Galerkin shallow water equation model' 
     &     /,10X,'for coastal seas and estuarine research',
     &     ///,7X,'-  Based off of the ADCIRC model created by',
     &     //,10X,'R.A. Luettich, Jr., University of North Carolina',
     &     /,10X,'J.J. Westerink, University of Notre Dame',
     &     //,7X,'-  The ADCIRC source code has been copyrighted by',
     &     /,10X,'R.A. Luettich, Jr. & J.J. Westerink, 1994-2001',
     &     /,10X,'No part of the adcirc base code may be reproduced',
     &     /,10X,'or redistributed without the written permission of',
     &     /,10X,'the above authors.',
     &     ///,5X,'The DG version of the code was written largely ab initio,',   
     &     /,5X,'though some data structures from the CG version were used.',
     &     //,7X,'Developed by:',
     &     /,10X,'Ethan J. Kubatko, The Ohio State University (EJK)',
     &     /,10X,'Clint N. Dawson, UT ICES (cnd)',
     &     /,10X,'Shintaro Bunya, The University of Tokyo (sb)',
     &     /,10X,'Craig Michoski, UT ICES (cem)',
     &     /,10X,'Christopher Mirabito, MIT',
     &     /,10X,'Damrongsak Wirasaet, University of Notre Dame',
     &     /,10X,'Casey Dietrich, North Carolina State University',//)


C.....Write out header information describing how code has been set up

      WRITE(16,1210)
 1210 FORMAT(//,1X,'THE SOURCE CODE HAS BEEN CONFIGURED ',
     &     'BY THE PREPROCESSOR AS FOLLOWS:',/)

#ifdef C3DDSS
      WRITE(16,*) '      - 3D DSS MODEL OPTION'
#endif 

#ifdef C3DVS
      WRITE(16,*) '      - 3D VS MODEL OPTION'
#else
      WRITE(16,*) '      - 2D DEPTH INTEGRATED MODEL OPTION'
#endif 

#ifdef CMACHSUN
      WRITE(16,*) '      - CODE SETUP TO RUN ON SUN 4 OR SPARC ',
     &     'COMPUTERS'
#endif

#ifdef REAL4  
      WRITE(16,*) '      - CODE SETUP TO RUN WITH 4 byte REALS'
#else
      WRITE(16,*) '      - CODE SETUP TO RUN WITH 8 byte REALS'
#endif

#ifdef CVEC
      WRITE(16,*) '      - CODE OPTIMIZED FOR A VECTOR COMPUTER'
#endif

#ifdef CSCA
      WRITE(16,*) '      - CODE OPTIMIZED FOR A SCALAR COMPUTER'
#endif

      WRITE(16,*) '      - NONVECTORIZABLE PARTS OF CODE OPTIMIZED FOR',
     &     ' MEMORY'
      WRITE(16,*) '      - CODE WILL USE JCG ITERATIVE GWCE SOLVER'
      WRITE(16,1112)

C.....Input from unit 15 and output to unit 16 rundescription and run ID

      READ(15,'(A32)') RUNDES
      READ(15,'(A24)') RUNID
      WRITE(16,1) RUNDES
 1    FORMAT(//,1X,'RUN DESCRIPTION : ',A32)
      WRITE(16,209) RUNID
 209  FORMAT(/,1X,'RUN IDENTIFICATION : ',A24)

C.....Read and process NFOVER - nonfatal error override otion

      READ(15,*) NFOVER
      WRITE(16,1112)
      WRITE(16,1250)
 1250 FORMAT(//,1X,'GENERAL RUN INFORMATION',/)
      IF(NFOVER.EQ.1) THEN
         WRITE(16,1951) NFOVER
 1951    FORMAT(5X,'NFOVER = ',I2,
     &        /,9X,'IF NON-FATAL ERRORS ARE DETECTED, THEY WILL BE ',
     &        'CORRECTED AND EXECUTION CONTINUED')
      ELSE
         WRITE(16,1952) NFOVER
 1952    FORMAT(/,5X,'NFOVER = ',I3,
     &        /,9X,'NON-FATAL ERRORS WILL STOP EXECUTION ',/)
      ENDIF

C.....Read and process NABOUT - abbreviated unit 16 output option

      READ(15,*) NABOUT
      IF (NABOUT.EQ.1) THEN
         WRITE(16,3501) NABOUT
 3501    FORMAT(5X,'NABOUT = ',I2,
     &        /,9X,'ABREVIATED OUTPUT WILL BE PROVIDED TO UNIT 16',
     &        /,9X,'UNIT 14, 21, 22 INPUT DATA WILL NOT BE ECHO PRINTED',/)
      ELSE
         WRITE(16,3502) NABOUT
 3502    FORMAT(/,5X,'NABOUT = ',I3,
     &        /,9X,'DETAILED OUTPUT WILL BE PROVIDED TO UNIT 16',
     &        /,9X,'UNIT 14, 15, 21, 22 INPUT DATA WILL BE ECHO PRINTED',/)
      ENDIF

C.....Read and process NSCREEN - screen ouput option

      READ(15,*) NSCREEN
      NSCREEN_INC = NSCREEN
      IF (NSCREEN.NE.0) NSCREEN = 1
      IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) THEN
         WRITE(16,3561) NSCREEN
 3561    FORMAT(5X,'NSCREEN = ',I2,
     &        /,9X,'SCREEN OUTPUT WILL BE PROVIDED TO UNIT 6',/)
      ELSE
         WRITE(16,3562) NSCREEN
 3562    FORMAT(/,5X,'NSCREEN = ',I3,
     &        /,9X,'SCREEN OUTPUT WILL NOT BE PROVIDED TO UNIT 6',/)
      ENDIF
      
C.....Read and process IHOT - hot start option

      READ(15,*) IHOT
      IF ((IHOT.NE.0).AND.(IHOT.NE.67).AND.(IHOT.NE.68)) THEN
         IF ((NSCREEN.EQ.1).AND.(MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'IHOT =',IHOT
            WRITE(6,9732)
            WRITE(6,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'IHOT =',IHOT
         WRITE(16,9732)
         WRITE(16,9973)       
 9732    FORMAT(/,1X,'Your selection of IHOT (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
         STOP
      ENDIF
      IF (IHOT.NE.0) THEN
         WRITE(16,9733) IHOT
 9733    FORMAT(/,5X,'dgswem will be hot started using information ',
     &        'on UNIT ',I2)
      ELSE
         WRITE(16,9734)
 9734    FORMAT(/,5X,'dgswem will be cold started')
      ENDIF
#ifdef SWAN
Casey 100205: Enable hot-start file generation by SWAN.
      SwanHotStartUnit = IHOT
#endif     


C.....Read and process ICS - cartesian/spherical coordinate option

      READ(15,*) ICS
      IF ((ICS.NE.1).AND.(ICS.NE.2)) THEN
         IF ((NSCREEN.EQ.1).AND.(MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'ICS =',ICS
            WRITE(6,9735)
            WRITE(6,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'ICS =',ICS
         WRITE(16,9735)
         WRITE(16,9973)
 9735    FORMAT(/,1X,'Your selection of ICS (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
         STOP
      ENDIF
      IF (ICS.EQ.1) THEN
         WRITE(16,9736) ICS
 9736    FORMAT(/,5X,'ICS = ',I2,
     &        /,9X,'Governing equations are in Cartesian coordinates')
      ELSE
         WRITE(16,9737) ICS
 9737    FORMAT(/,5X,'ICS = ',I2,
     &        /,9X,'Governing equations are in Spherical coordinates',
     &        /,9X,'mapped using a CPP projection')
      ENDIF

C.....Read and process IM - 2D/3D model option

      READ(15,*) IM
      IF (IM.EQ.0) THEN
         C2DDI = .TRUE.
      ELSEIF (IM.EQ.1) THEN
         C3D  = .TRUE.
         C3DVS  = .TRUE.
      ELSEIF (IM.EQ.2) THEN
         C3D  = .TRUE.
         C3DDSS = .TRUE.
      ELSEIF (IM.EQ.10) THEN
         C2DDI = .TRUE.
      ELSE
         IF ((NSCREEN.EQ.1).AND.(MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'IM =',IM
            WRITE(6,9721)
            WRITE(6,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'IM =',IM
         WRITE(16,9721)
         WRITE(16,9973)
 9721    FORMAT(/,1X,'Your selection of IM (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
         STOP
      ENDIF

C.....Read and process NOLIBF - nonlinear bottom friction option

      READ(15,*) NOLIBF
      IF ((NOLIBF.LT.0).OR.(NOLIBF.GT.2)) THEN
         IF ((NSCREEN.EQ.1).AND.(MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'NOLIBF =',NOLIBF
            WRITE(6,9722)
            WRITE(6,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NOLIBF =',NOLIBF
         WRITE(16,9722)
         WRITE(16,9973)
 9722    FORMAT(/,1X,'Your selection of NOLIBF (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
         STOP
      ENDIF
      WRITE(16,9845) NOLIBF
 9845 FORMAT(/,5X,'NOLIBF = ',I3)
      IF (NOLIBF.EQ.0) WRITE(16,2050)
 2050 FORMAT(9X,'THE MODEL WILL USE LINEAR BOTTOM FRICTION')
      IF (NOLIBF.EQ.1) WRITE(16,2051)
 2051 FORMAT(9X,'THE MODEL WILL USE STANDARD QUADRATIC BOTTOM FRICTION')
      IF (NOLIBF.EQ.2) WRITE(16,2052)
 2052 FORMAT(9X,'THE MODEL WILL USE STANDARD QUADRATIC BOTTOM FRICTION',
     &     'IN DEEP WATER ',
     &     /,9X,'AND A FRICTION FACTOR THAT INCREASES AS THE DEPTH ',
     &     'DECREASES IN SHALLOW WATER')

C.....Read and process NOLIFA - nonlinear finite amplitude option

      READ(15,*) NOLIFA
      IF ((NOLIFA.LT.0).OR.(NOLIFA.GT.3)) THEN
         IF ((NSCREEN.EQ.1).AND.(MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'NOLIFA =',NOLIFA
            WRITE(6,9723)
            WRITE(6,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NOLIFA =',NOLIFA
         WRITE(16,9723)
         WRITE(16,9973)
 9723    FORMAT(/,1X,'Your selection of NOLIFA (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
         STOP
      ENDIF
      WRITE(16,9846) NOLIFA
 9846 FORMAT(/,5X,'NOLIFA = ',I3)
      IF (NOLIFA.EQ.0) WRITE(16,2053)
 2053 FORMAT(9X,'THE MODEL WILL NOT USE FINITE AMPLITUDE TERMS OR ',
     &     'WETTING AND DRYING')
      IF (NOLIFA.EQ.1) WRITE(16,2054)
 2054 FORMAT(9X,'THE MODEL WILL USE FINITE AMPLITUDE TERMS BUT NO ',
     &     'WETTING AND DRYING')
      IF (NOLIFA.EQ.2) WRITE(16,2049)
 2049 FORMAT(9X,'THE MODEL WILL USE FINITE AMPLITUDE TERMS AND ',
     &     'WETTING AND DRYING')
      IF (NOLIFA.EQ.3) WRITE(16,2048)
 2048 FORMAT(9X,'THE MODEL WILL USE FINITE AMPLITUDE TERMS AND ',
     &     'WETTING AND DRYING',/,10X,
     &     'AND INCLUDES THE ABILITY TO INITIALIZE ',
     &     'NODES WITH DEPTHS GREATER THAN H0 AS DRY')         
      NSTARTDRY = 0
      IF (NOLIFA.EQ.3) THEN
         NOLIFA = 2
         NSTARTDRY = 1
      ENDIF

C.....Read and process NOLICA - advective term spatial gradinet

      READ(15,*) NOLICA
      IF ((NOLICA.LT.0).OR.(NOLICA.GT.1)) THEN
         IF ((NSCREEN.EQ.1).AND.(MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'NOLICA =',NOLICA
            WRITE(6,9724)
            WRITE(6,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NOLICA =',NOLICA
         WRITE(16,9724)
         WRITE(16,9973)
 9724    FORMAT(/,1X,'Your selection of NOLICA (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
         STOP
      ENDIF
      WRITE(16,9847) NOLICA
 9847 FORMAT(/,5X,'NOLICA = ',I3)
      IF (NOLICA.EQ.0) WRITE(16,2055)
 2055 FORMAT(9X,'THE MODEL WILL NOT USE SPATIAL DERIVATIVE ',
     &     'COMPONENTS OF THE ADVECTIVE TERMS')
      IF (NOLICA.EQ.1) WRITE(16,2056)
 2056 FORMAT(9X,'THE MODEL WILL USE SPATIAL DERIVATIVE ',
     &     'COMPONENTS OF THE ADVECTIVE TERMS')

C.....Read and process NOLICAT - GWCE advective term time derivative

      READ(15,*) NOLICAT
      IF ((NOLICAT.LT.0).OR.(NOLICAT.GT.1)) THEN
         IF ((NSCREEN.EQ.1).AND.(MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'NOLICAT =',NOLICAT
            WRITE(6,9725)
            WRITE(6,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NOLICAT =',NOLICAT
         WRITE(16,9725)
         WRITE(16,9973)
 9725    FORMAT(/,1X,'Your selection of NOLICAT (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
         STOP
      ENDIF
      IF ((NOLIFA.GE.1).AND.(NOLICAT.EQ.0)) THEN
         IF ((NSCREEN.EQ.1).AND.(MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'NOLICAT =',NOLICAT
            WRITE(6,9726)
            IF (NFOVER.EQ.1) THEN
               WRITE(6,9974)
            ELSE
               WRITE(6,9973)
            ENDIF
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NOLICAT =',NOLICAT
         WRITE(16,9726)
         WRITE(16,9974)
 9726    FORMAT(/,1X,'Your selection of NOLICAT (a UNIT 15 input ',
     &        'parameter) is inconsistent with your ',
     &        /,1X,'selection of NOLIFA and may lead to mass ',
     &        'balance problems')
         IF (NFOVER.EQ.1) THEN
            WRITE(6,9974)
         ELSE
            WRITE(6,9973)
            STOP
         ENDIF
      ENDIF
      IF ((NOLIFA.EQ.0).AND.(NOLICAT.EQ.1)) THEN
         IF ((NSCREEN.EQ.1).AND.(MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'NOLICAT =',NOLICAT
            WRITE(6,9726)
            IF (NFOVER.EQ.1) THEN
               WRITE(6,9974)
            ELSE
               WRITE(6,9973)
            ENDIF
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NOLICAT =',NOLICAT
         WRITE(16,9726)
         WRITE(16,9974)
         IF (NFOVER.EQ.1) THEN
            WRITE(6,9974)
         ELSE
            WRITE(6,9973)
            STOP
         ENDIF
      ENDIF
      IF (NOLICA.NE.NOLICAT) THEN
         IF ((NSCREEN.EQ.1).AND.(MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'NOLICAT =',NOLICAT
            WRITE(6,9727)
            IF (NFOVER.EQ.1) THEN
               WRITE(6,9974)
            ELSE
               WRITE(6,9973)
            ENDIF
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NOLICAT =',NOLICAT
         WRITE(16,9727)
         WRITE(16,9974)
 9727    FORMAT(/,1X,'Your selection of NOLICAT (a UNIT 15 input ',
     &        'parameter) is inconsistent with your ',
     &        /,1X,'selection of NOLICA and may lead to mass ',
     &        'balance problems')
         IF (NFOVER.EQ.1) THEN
            WRITE(6,9974)
         ELSE
            WRITE(6,9973)
            STOP
         ENDIF
      ENDIF
      WRITE(16,9848) NOLICAT
 9848 FORMAT(/,5X,'NOLICAT = ',I3)
      IF (NOLICAT.EQ.0) WRITE(16,2057)
 2057 FORMAT(9X,'THE MODEL WILL NOT USE TIME DERIVATIVE COMPONENTS ',
     &     /,9X,'OF THE ADVECTIVE TERMS IN THE GWCE')
      IF (NOLICAT.EQ.1) WRITE(16,2058)
 2058 FORMAT(9X,'THE MODEL WILL USE TIME DERIVATIVE COMPONENTS ',
     &     /,9X,'OF THE ADVECTIVE TERMS IN THE GWCE')

C.....Read and process NWP - spatially varying bottom friction

      READ(15,*) NWP
      CALL ReadNodalAttr(NSCREEN, ScreenUnit, MYPROC, NABOUT) ! Ek added call to nodalatt

C.....Read and process NCOR - spatially varying Coriolis parameter

      READ(15,*) NCOR
      IF ((NCOR.NE.0).AND.(NCOR.NE.1)) THEN
         IF ((NSCREEN.EQ.1).AND.(MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'NCOR =',NCOR
            WRITE(6,9729)
            WRITE(6,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NCOR =',NCOR
         WRITE(16,9729)
         WRITE(16,9973)
 9729    FORMAT(/,1X,'Your selection of NCOR (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
         STOP
      ENDIF
      IF ((ICS.EQ.1).AND.(NCOR.EQ.1)) THEN
         IF ((NSCREEN.EQ.1).AND.(MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'NCOR =',NCOR
            WRITE(6,9730)
            WRITE(6,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NCOR =',NCOR
         WRITE(16,9730)
         WRITE(16,9973)
 9730    FORMAT(/,1X,'Your selection of NCOR (a UNIT 15 input ',
     &        'parameter) is inconsistent with your ',
     &        /,1X,'selection of coordinate systems.  Spatially ',
     &        'variable Coriolis should be used only with ',
     &        /,1X,'Spherical coordinates')
         STOP
      ENDIF
      IF (NCOR.EQ.0) THEN
         WRITE(16,233) NCOR
 233     FORMAT(/,5X,'NCOR = ',I2,
     &        /,9X,'A CONSTANT VALUE OF THE CORIOLIS PARAMETER WILL BE ',
     &        /,9X,'USED THROUGHOUT THE DOMAIN')
      ELSE
         WRITE(16,234) NCOR
 234     FORMAT(/,5X,'NCOR = ',I2,
     &        /,9X,'SPATIALLY VARYING CORIOLIS VALUES WILL BE COMPUTED ',
     &        'FROM INPUT LATITUDES')
      ENDIF

C.....Read and process NTIP - tidal potential forcing

      READ(15,*) NTIP
      IF ((NTIP.LT.0).OR.(NTIP.GT.2)) THEN
         IF ((NSCREEN.EQ.1).AND.(MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'NTIP =',NTIP
            WRITE(6,9710)
            WRITE(6,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NTIP =',NTIP
         WRITE(16,9710)
         WRITE(16,9973)
 9710    FORMAT(/,1X,'Your selection of NTIP (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
         STOP
      ENDIF

      IF ((ICS.EQ.1).AND.(NTIP.GE.1)) THEN
         IF ((NSCREEN.EQ.1).AND.(MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'NTIP =',NTIP
            WRITE(6,9711)
            WRITE(6,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NTIP =',NTIP
         WRITE(16,9711)
         WRITE(16,9973)
 9711    FORMAT(/,1X,'Your selection of NTIP (a UNIT 15 input ',
     &        'parameter) is inconsistent with your ',
     &        /,1X,'selection of coordinate systems.  Tidal',
     &        'potential forcing should be used only with ',
     &        /,1X,'Spherical coordinates')
         STOP
      ENDIF
      IF (NTIP.NE.0) CTIP = .TRUE.
      IF (NTIP.EQ.0) THEN
         WRITE(16,235) NTIP
 235    FORMAT(/,5X,'NTIP = ',I2,
     &    /,9X,'TIDAL POTENTIAL FORCING IS NOT USED IN THE COMPUTATION')
      ENDIF
      IF (NTIP.GE.1) THEN
        WRITE(16,236) NTIP
 236     FORMAT(/,5X,'NTIP = ',I2,
     &        /,9X,'TIDAL POTENTIAL FORCING IS USED IN THE COMPUTATION ',
     &        'BASED ON INPUT LONGITUDES/LATITUDES')
      ENDIF
      IF (NTIP.EQ.2) THEN
         WRITE(16,239)
 239     FORMAT(9X,'SELF ATTRACTION/LOAD TIDE FORCING IS ALSO USED ',
     &        'IN THE COMPUTATION')
      ENDIF

C.....Read and process NWS - wind and pressure forcing & wave rad stress

      READ(15,*) NWS
      IF ( (NWS.NE.0)  .AND.(NWS.NE.1)       .AND.(ABS(NWS).NE.2)  .AND.
     &     (NWS.NE.3)  .AND.(ABS(NWS).NE.4)  .AND.(ABS(NWS).NE.5)  .AND.
     &     (NWS.NE.6)  .AND.(NWS.NE.10)      .AND.(NWS.NE.11)      .AND.
     &     (ABS(NWS).NE.12 ).AND.
     &     (NWS.NE.100).AND.(NWS.NE.101)     .AND.(ABS(NWS).NE.102).AND.
     &     (NWS.NE.103).AND.(ABS(NWS).NE.104).AND.(ABS(NWS).NE.105).AND.
     &     (NWS.NE.106).AND.(NWS.NE.110)     .AND.(NWS.NE.111)     .AND.
     &     (NWS.NE.200).AND.(NWS.NE.201)     .AND.(ABS(NWS).NE.202).AND.
     &     (NWS.NE.203).AND.(ABS(NWS).NE.204).AND.(ABS(NWS).NE.205).AND.
     &     (NWS.NE.8).AND.
     &     (NWS.NE.206).AND.(NWS.NE.210)     .AND.(NWS.NE.211) .AND.
Casey 101118: Added the following cases for coupling to unstructured SWAN.
     &   (ABS(NWS).NE.300).AND.(ABS(NWS).NE.303).AND.(ABS(NWS).NE.304).AND.
     &   (ABS(NWS).NE.305).AND.(ABS(NWS).NE.306).AND.(ABS(NWS).NE.308).AND.
     &   (ABS(NWS).NE.309).AND.(ABS(NWS).NE.310).AND.(ABS(NWS).NE.311).AND.
     &   (ABS(NWS).NE.312)) THEN
         IF ((NSCREEN.EQ.1).AND.(MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'NWS =',NWS
            WRITE(6,9712)
            WRITE(6,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NWS =',NWS
         WRITE(16,9712)
         WRITE(16,9973)
 9712    FORMAT(/,1X,'Your selection of NWS (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
         STOP
      ENDIF

C.....Set wave radiation stress flag and adjust NWS accordingly
Casey 101118: Had to make some changes in this section.
      NRS = 0
      FRW = 0
      IF (ABS(NWS/100).EQ.1) THEN
        NRS=1
        NWS = (ABS(NWS) - 100)*(NWS/ABS(NWS))
      ENDIF
      IF (ABS(NWS/100).EQ.2) THEN
        NRS = 1
        FRW = 1
        NWS = (ABS(NWS) - 200)*(NWS/ABS(NWS))
      ENDIF
#ifdef SWAN
Casey 101118: Added the option for coupling directly to SWAN.
      IF(ABS(NWS/100).EQ.3) THEN
        NRS = 3
        NWS = (ABS(NWS) - 300)*(NWS/ABS(NWS))
      ENDIF
#endif

      IF (NWS.EQ.0) THEN
         WRITE(16,237) NWS
 237     FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'WIND STRESS OR SURFACE PRESSURE ARE NOT USED TO FORCE',
     &         'THE COMPUTATION')
      ENDIF
      IF (NWS.EQ.1) THEN
        WRITE(16,238) NWS
 238    FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'VALUES ARE READ AT ADCIRC GRID NODES FROM UNIT 22',
     &    /,9X,' EVERY MODEL TIME STEP')
      ENDIF
      IF (NWS.EQ.2) THEN
        WRITE(16,2381) NWS
 2381   FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'VALUES ARE READ AT ADCIRC GRID NODES FROM UNIT 22',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT TIME=STATIM.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.')
      ENDIF

      IF (NWS.EQ.-2) THEN
        WRITE(16,2380) NWS
 2380   FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'VALUES ARE READ AT ADCIRC GRID NODES FROM UNIT 22',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT THE TIME OF THE HOT START.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &       /,9X,'WITH THE MODEL TIME STEP.')
      ENDIF

      IF (NWS.EQ.3) THEN
         WRITE(16,2382) NWS
 2382    FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'WIND STRESS ONLY IS USED TO FORCE THE COMPUTATION.',
     &    /,9X,'WIND SPEEDS AND DIRECTIONS ARE READ FROM A FLEET ',
     &    /,9X,'NUMERIC FORMAT FILE AT UNIT 22 AND INTERPOLATED TO',
     &    /,9X,'THE ADCIRC GRID. ',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.',
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &        'DRAG LAW.')
      ENDIF
      IF (NWS.EQ.4) THEN
         WRITE(16,2383) NWS
 2383    FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ AT SELECTED',
     &    /,9X,'ADCIRC GRID NODES FROM A PBL FILE AT UNIT 22.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT TIME=STATIM.',
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &        'DRAG LAW.')
      ENDIF
      IF (NWS.EQ.-4) THEN
         WRITE(16,2388) NWS
 2388    FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ AT SELECTED',
     &    /,9X,'ADCIRC GRID NODES FROM A PBL FILE AT UNIT 22.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT THE TIME OF THE HOT START.',
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &        'DRAG LAW.')
      ENDIF
      IF (NWS.EQ.5) THEN
         WRITE(16,2384) NWS
 2384    FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ AT ADCIRC ',
     &    /,9X,'GRID NODES FROM UNIT 22',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT TIME=STATIM.',
     &     /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &     'DRAG LAW.')
      ENDIF
      IF (NWS.EQ.-5) THEN
         WRITE(16,2389) NWS
 2389    FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ AT ADCIRC ',
     &    /,9X,'GRID NODES FROM UNIT 22',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT THE TIME OF THE HOT START.',
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &        'DRAG LAW.')
      ENDIF
      IF (NWS.EQ.6) THEN
         WRITE(16,2385) NWS
 2385    FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ FROM A ',
     &    /,9X,'REGULARLY SPACED GRID FROM UNIT 22',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP AND IN SPACE TO BRING THE ',
     &    /,9X,'MET DATA FROM A REGULAR GRID TO THE ADCIRC GRID.'
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &        'DRAG LAW.')
      ENDIF
      IF (NWS.EQ.10) THEN
         WRITE(16,2386) NWS
 2386    FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ EVERY N',
     &    /,9X,' HOURS FROM A DIFFERENT FILE AT UNITS 200, 200+N,',
     &    ' 200+2N, ETC.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP AND IN SPACE TO BRING THE ',
     &    /,9X,'MET DATA FROM A GAUSSIAN GRID TO THE ADCRIC GRID.',
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &        'DRAG LAW.')
      ENDIF
      IF (NWS.EQ.11) THEN
         WRITE(16,2387) NWS
 2387    FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ EVERY 3 ',
     &    /,9X,'HOURS FROM ETA-29 FILES AT UNITS 200, 201, 202, ETC.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP AND IN SPACE TO BRING THE ',
     &    /,9X,'WIND DATA FROM THE 29 KM E GRID TO THE ADCRIC GRID.',
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &        'DRAG LAW.')
      ENDIF
      IF (NRS.EQ.0) THEN
         WRITE(16,2390) NRS
 2390    FORMAT(/,5X,'NRS = ',I2,
     &        /,9X,'WAVE RADIATION STRESS IS NOT USED TO FORCE THE ',
     &        'COMPUTATION')
      ENDIF
      
C.....ek added for NWS=-12,12

      IF(NWS.EQ.12) THEN
         WRITE(16,12384) NWS
12384    FORMAT(/,5X,'NWS = ',I2,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ FROM RAW ',
     &    /,9X,'OWI DATA FILES (UNIT 221-224).',
     &    /,9X,'META DATA IS READ FROM UNIT 220.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT TIME=STATIM.',
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &        'DRAG LAW.')
      ENDIF
      IF(NWS.EQ.-12) THEN
         WRITE(16,12389) NWS
12389    FORMAT(/,5X,'NWS = ',I3,
     &    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',
     &    /,9X,' THE COMPUTATION',
     &    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ FROM RAW ',
     &    /,9X,'OWI DATA FILES (UNIT 221-224).',
     &    /,9X,'META DATA IS READ FROM UNIT 220.',
     &    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',
     &    /,9X,'WITH THE MODEL TIME STEP.',
     &    /,9X,'THE UNIT 22 FILE BEGINS AT THE TIME OF THE HOT START.',
     &    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',
     &        'DRAG LAW.')
      ENDIF
      IF (NRS.EQ.1) THEN
         WRITE(16,2391) NRS
 2391    FORMAT(/,5X,'NRS = ',I2,
     &    /,9X,'WAVE RADIATION STRESS IS USED TO FORCE THE COMPUTATION',
     &    /,9X,'STRESSES ARE READ AT SELECTED ADCIRC GRID NODES FROM A',
     &    /,9X,'PBL TYPE FILE AT UNIT 23.  INTERPOLATION IN TIME IS ',
     &    /,9X,'DONE TO SYNC THE STRESS DATA WITH THE MODEL TIME STEP.',
     &    /,9X,'FOR A COLD START, THE UNIT 23 FILE BEGINS AT THE TIME ',
     &    /,9X,'OF THE COLD START.  FOR A HOT START, THE UNIT 23 FILE ',
     &    /,9X,'BEGINS AT THE TIME OF THE HOT START.')
      ENDIF
#ifdef SWAN
Casey 101118: Added the following lines.
      IF(NRS.EQ.3) THEN
         WRITE(16,2393) NRS
 2393    FORMAT(/,5X,'NRS = ',I2,
     &    /,9X,'WAVES WILL BE COUPLED TO SWAN!')
      ENDIF
#endif

C.....Read and process NRAMP - whether a ramp function will be used

      READ(15,*) NRAMP
      IF ((NRAMP.NE.0).AND.(NRAMP.GT.7)) THEN
         IF ((NSCREEN.EQ.1).AND.(MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'NRAMP =',NRAMP
            WRITE(6,9713)
            WRITE(6,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NRAMP =',NRAMP
         WRITE(16,9713)
         WRITE(16,9973)
 9713    FORMAT(/,1X,'Your selection of NRAMP (a UNIT 15 input ',
     &        'parameter) is not an allowable value')
         STOP
      ENDIF
      IF (NRAMP.EQ.0) THEN
         WRITE(16,240) NRAMP
 240     FORMAT(/,5X,'NRAMP = ',I2,
     &   /,9X,'NO RAMP FUNCTION IS USED IN THE COMPUTATION')
      ELSE
         WRITE(16,241) NRAMP
 241     FORMAT(/,5X,'NRAMP = ',I2,
     &   /,9X,'A HYPERBOLIC TANGENT RAMP IS APPLIED TO THE FORCING ',
     &        'FUNCTIONS')
      ENDIF

C.....Read and process G - gravity
C...  
      READ(15,*) G
      G2ROOT = SQRT(G/2.0d0)
      IF ((ICS.EQ.2).AND.(abs(G-9.81).gt.0.01)) THEN
         IF ((NSCREEN.EQ.1).AND.(MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'G =',G
            WRITE(6,9714)
            WRITE(6,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'G =',G
         WRITE(16,9714)
         WRITE(16,9973)
 9714    FORMAT(/,1X,'Your specification of the gravitational ',
     &        'constant, G, (a UNIT 15 input) is not ',
     &        /,1X,'consistant with the use of spherical coordinates.',
     &        '  G must be in units of m/s^2')
         STOP
      ENDIF
      WRITE(16,5) G
    5 FORMAT(///,5X,'GRAVITATIONAL CONSTANT G =',F10.5,/)

C.....Read and process TAU0 - weighting coefficient in the GWCE

      READ(15,*) TAU0
      IF (TAU0.LT.0) THEN
         WRITE(16,6)
 6       FORMAT(/,5X,'WEIGHTING COEFFICIENT FOR THE GENERALIZED',
     &        ' WAVE CONTINUITY EQUATION :',
     &        /,5x,'THIS VALUE WILL BE  SELECTED BASED ON NODAL DEPTH',
     &        ' ONCE DEPTHS HAVE BEEN PROCESSED',    
     &        /,5X,' DEPTH > 10       -> TAU0 = 0.005  ',    
     &        /,5X,' 10 >/ DEPTH        -> TAU0 = 0.020 ',  
     &        /,5X,' STARTDRY VALUE = -77777  -> TAU0 = 0.020 ',  
     &        /,5X,' STARTDRY VALUE = -88888  -> TAU0 = 0.020 ',/)  
      ELSE
         WRITE(16,7) TAU0
 7       FORMAT(/,5X,'WEIGHTING COEFFICIENT FOR THE GENERALIZED',
     &        ' WAVE CONTINUITY EQUATION :',
     &        /,5X, 'TAU0 = ',E14.8,2X,'1/sec',/)
      ENDIF

C.....Input from unit 15 and output to unit 16 time integration info

      WRITE(16,1112)
      WRITE(16,245)
 245  FORMAT(//,1X,'TIME INTEGRATION INFORMATION',//)

C.....Read and process DT - model time step

      READ(15,*) DTDP
      DT = DTDP
      WRITE(16,9) DTDP
    9 FORMAT(5X,'TIME STEP =',F12.6,5X,'SECONDS',/)

C.....Read and process STATIM - simulation starting time

      READ(15,*) STATIM
      WRITE(16,1113) STATIM
 1113 FORMAT(5X,'STARTING TIME FOR SIMULATION = ',F14.6,' DAYS',/)

C.....Read anf process REFTIM - harmonic reference time

      READ(15,*) REFTIM
      WRITE(16,1115) REFTIM
 1115 FORMAT(5X,'Harmonic REFERENCE TIME = ',F14.6,' DAYS',/)

C.....Read in and process additional timing information for wind.
Casey 101118: Changed NRS.EQ.1 to NRS.GE.1 throughout this section.
      IF ((NWS.EQ.0).AND.(NRS.GE.1)) READ(15,*) RSTIMINC
      IF ((NWS.EQ.1).AND.(NRS.GE.1)) READ(15,*) RSTIMINC
      IF (ABS(NWS).EQ.2) THEN
         IF (NRS.EQ.0) READ(15,*) WTIMINC
         IF (NRS.GE.1) READ(15,*) WTIMINC,RSTIMINC
      ENDIF
      IF (NWS.EQ.3) THEN
         READ(15,*) IREFYR,IREFMO,IREFDAY,IREFHR,IREFMIN,REFSEC
         WRITE(16,1116) IREFMO,IREFDAY,IREFYR,IREFHR,IREFMIN,REFSEC
 1116    FORMAT(5X,'WIND REFERENCE TIME FOR SIMULATION = ',
     &        I2,'/',I2,'/',I2,'  ',I2,':',I2,':',f7.4,/)
         CALL TIMECONV(IREFYR,IREFMO,IREFDAY,IREFHR,IREFMIN,REFSEC,
     &        WREFTIM, MyProc, NScreen, ScreenUnit )
         IF (NRS.EQ.0) READ(15,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC,
     &        WLONINC,WTIMINC
         IF (NRS.GE.1) READ(15,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC,
     &        WLONINC,WTIMINC,RSTIMINC
      ENDIF
      IF (ABS(NWS).EQ.4) THEN
         IF (NRS.EQ.0) READ(15,*) WTIMINC
         IF (NRS.GE.1) READ(15,*) WTIMINC,RSTIMINC
      ENDIF

      IF (ABS(NWS).EQ.5) THEN
         IF (NRS.EQ.0) READ(15,*) WTIMINC
         IF (NRS.GE.1) READ(15,*) WTIMINC,RSTIMINC
      ENDIF
      IF (NWS.EQ.6) THEN
         IF (NRS.EQ.0) READ(15,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC,
     &        WLONINC,WTIMINC
         IF (NRS.GE.1) READ(15,*) NWLAT,NWLON,WLATMAX,WLONMIN,WLATINC,
     &        WLONINC,WTIMINC,RSTIMINC
      ENDIF
      IF(ABS(NWS).EQ.8) THEN
         IF(NRS.EQ.0) THEN
            READ(15,*) IREFYR,IREFMO,IREFDAY,IREFHR,StormNumber,BLAdj
         ELSEIF (NRS.GE.1) THEN
            READ(15,*) IREFYR,IREFMO,IREFDAY,IREFHR,StormNumber,BLAdj,
     &           RSTIMINC
            WRITE(16,6111) IREFMO,IREFDAY,IREFYR,IREFHR
 6111       FORMAT(5X,'WIND REFERENCE TIME FOR SIMULATION = ',
     &           I2,'/',I2,'/',I2,'  ',I2,'H',/)
         ENDIF
         CALL TIMECONV(IREFYR,IREFMO,IREFDAY,IREFHR,0,0.0d0,
     &        WindRefTime, MyProc, NScreen, ScreenUnit)
      ENDIF

      IF (NWS.EQ.10) THEN
         NWLAT=190
         NWLON=384
         IF (NRS.EQ.0) READ(15,*) WTIMINC
         IF (NRS.GE.1) READ(15,*) WTIMINC,RSTIMINC
      ENDIF
      IF (NWS.EQ.11) THEN
         NWLAT=271
         NWLON=181
         WTIMINC=10800.
         IF (NRS.GE.1) READ(15,*) RSTIMINC
      ENDIF
      IF (NWS.EQ.11) THEN
         NWLAT=271
         NWLON=181
         WTIMINC=10800.
         IF (NRS.GE.1) READ(15,*) RSTIMINC
      ENDIF
C.....ek added NWS=12 (OWI format) from version 46

      IF(ABS(NWS).EQ.12) THEN
         IF(NRS.EQ.0) READ(15,*) WTIMINC
         IF(NRS.GE.1) READ(15,*) WTIMINC,RSTIMINC ! sb46.28sb03
      ENDIF

      IF (NWS.NE.0) WRITE(16,1117) WTIMINC
 1117 FORMAT(5X,'WIND TIME INCREMENT (SEC) = ',F10.2,/)
      IF (NRS.NE.0) WRITE(16,1118) RSTIMINC
 1118 FORMAT(5X,'RADIATION STRESS TIME INCREMENT (SEC) = ',F10.2,/)

C.....Read and process RNDAY - simulation duration i days

      READ(15,*) RNDAY
      WRITE(16,10) RNDAY
 10   FORMAT(5X,'TOTAL LENGTH OF NUMERICAL SIMULATION =',F12.4,
     &     5X,'DAYS',/)

C.....Read and process effective length of hyperbolic tangent ramp

C...  
C...  READ AND PROCESS EFFECTIVE LENGTH OF THE HYPERBOLIC TANGENT RAMP(S)
C...  IN DAYS
C...  
C     jgf46.08 Add fine-grained ramp functions.
C     jgf46.21 Add FluxSettlingTime for IBTYPE=52 to accomodate
C     MS river during Katrina, split ramps for flux b.c.s into internal
C     and external.
      FluxSettlingTime = 0.0d0
      DRamp = 1.0d0
      SELECT CASE(NRamp)
C     ---------
      CASE(0,1)                 ! Either no ramp, or same ramp for all forcings
C     ---------
         READ(15,*) DRamp
         DRampIntFlux = DRamp
         DRampExtFlux = DRamp
         DRampElev    = DRamp
         DRampTip     = DRamp
         DRampMete    = DRamp
         DRampWRad    = DRamp
C     -------
      CASE(2)                   ! Ramp for external flux boundary conditions.
C     -------
         READ(15,*) DRamp,DRampExtFlux,FluxSettlingTime
         DRampIntFlux = DRamp
         DRampElev    = DRamp
         DRampTip     = DRamp
         DRampMete    = DRamp
         DRampWRad    = DRamp
C     -------
      CASE(3)                   ! Ramp for internal flux boundary conditions.
C     -------
         READ(15,*) DRamp,DRampExtFlux,FluxSettlingTime,DRampIntFlux
         DRampElev = DRamp
         DRampTip  = DRamp
         DRampMete = DRamp
         DRampWRad = DRamp
C     -------
      CASE(4)                   ! Ramp for surface elevation specified boundary conditions.
C     -------
         READ(15,*) DRamp,DRampExtFlux,FluxSettlingTime,DRampIntFlux,
     &        DRampElev
         DRampTip  = DRamp
         DRampMete = DRamp
         DRampWRad = DRamp
C     -------
      CASE(5)                   ! Ramp for tidal potential
C     -------
         READ(15,*) DRamp,DRampExtFlux,FluxSettlingTime,DRampIntFlux,
     &        DRampElev,DRampTip
         DRampMete = DRamp
         DRampWRad = DRamp
C     -------
      CASE(6)                   ! Ramp for wind and atmospheric pressure
C     -------
         READ(15,*) DRamp,DRampExtFlux,FluxSettlingTime,DRampIntFlux,
     &        DRampElev,DRampTip,DRampMete
         DRampWRad = DRamp
C     -------
      CASE(7)                   ! Ramp for wave radiation stress
C     -------
         READ(15,*) DRamp,DRampExtFlux,FluxSettlingTime,DRampIntFlux,
     &        DRampElev,DRampTip,DRampMete,DRampWRad
C     ------------
      CASE DEFAULT              ! fall-through
C     ------------
         IF(NSCREEN.NE.0.AND.MYPROC.EQ.0) THEN
            WRITE(ScreenUnit,9972)
            WRITE(ScreenUnit,*) 'NRAMP =',NRAMP
            WRITE(ScreenUnit,9713)
            WRITE(ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'NRAMP =',NRAMP
         WRITE(16,9713)
         WRITE(16,9973)
#ifdef CMPI
         call MESSAGE_FINI()
#endif
         STOP
      END SELECT

C.....Compute the total number of timesteps the ramping function is used
 
      IRAMPING = INT(DRAMP*(86400.D0/DTDP) + 0.5D0)

C.....Read GWCE time weighting factors

      READ(15,*) A00,B00,C00
      WRITE(16,14)
 14   FORMAT(//,5X,'TIME WEIGHTING FACTORS IN THE WAVE EQUATION :'/)
      WRITE(16,15) A00,B00,C00
 15   FORMAT(9X,'AT TIME LEVEL K+1 : ',F8.5,
     &     /,9X,'AT TIME LEVEL K   : ',F8.5,
     &     /,9X,'AT TIME LEVEL K-1 : ',F8.5,/)

C.....Read minimum depth or wet/dry parameters from unit 15

      IF (NOLIFA.NE.2) THEN
         READ(15,*) H0
         WRITE(16,16) H0
 16      FORMAT(//,5X,'THE BATHYMETRIC DEPTH AT ALL NODES WILL BE ',
     &        'INCREASED TO H0= ',F12.4,' IF NECESSARY'/)
      ENDIF
      IF (NOLIFA.EQ.2) THEN
         READ(15,*) H0,NODEDRYMIN,NODEWETMIN,VELMIN
         WRITE(16,17) H0,NODEWETMIN,VELMIN,NODEDRYMIN
 17      FORMAT(//,5X,'DRYING WILL OCCUR WHEN THE WATER DEPTH < H0',
     &        /,5X,'H0 = ',E16.8,
     &        /,5X,'AND NODEREP > NODEWETMIN = ',I6,' TIME STEPS',
     &        /,5X,'NODEREP = NUMBER OF TIME STEPS SINCE A NODE ',
     &        'CHANGED STATE (EITHER WETTED OR DRIED)',
     &        //,5X,'WETTING WILL OCCUR WHEN THERE IS A FAVORABLE ',
     &        'PRESSURE GRADIENT THAT',
     &        /,5X,'WOULD DRIVE A STEADY VELOCITY TOWARDS A DRY NODE',
     &        /,5X,'THAT IS GREATER THAN VELMIN = ',F10.5,
     &        /,5X,'AND NODEREP > NODEDRYMIN = ',I6,' TIME STEPS',/)
      ENDIF

C.....Read grid information from units 14 and 15

      READ(14,'(A24)') AGRID
      READ(14,*) NE,NP
      MNP = NP
      MNE = NE
      READ(15,*) SLAM0,SFEA0
      SLAM0 = SLAM0*DEG2RAD
      SFEA0 = SFEA0*DEG2RAD
      WRITE(16,1112)
      WRITE(16,246)
 246  FORMAT(//,1X,'GRID INFORMATION',//)

C.....Allocate arrays dimensioned by MNP and MNE

      CALL ALLOC_MAIN1()

C.....If ICS = 1 input nodal coordinates and bathymetry from unit 14
C.....If either NTIP or NCOR = 1 compute the inverse CPP projection

      IF (ICS.EQ.1) THEN
         DO I = 1,NP
            READ(14,*) JKI,X(JKI),Y(JKI),DP(JKI)
            IF(JKI.NE.I) THEN
               IF ((NSCREEN.EQ.1).AND.(MYPROC.EQ.0)) WRITE(6,99801)
               WRITE(16,99801)
99801          FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',
     &              'INPUT ERROR  !!!!!!!!!',
     &              //,1X,'YOUR NODE NUMBERING IS NOT SEQUENTIAL ',
     &              'CHECK YOUR UNIT 14 INPUT FILE CAREFULLY',//)
            ENDIF
C     IF ((NTIP.GE.1).OR.(NCOR.EQ.1))
            CALL INVCP(X(JKI),Y(JKI),SLAM(JKI),SFEA(JKI),SLAM0,SFEA0)
C     ENDIF
         ENDDO
      ENDIF
      
C.....If ICS = 2 input nodal coordinates and bathymetry from unit 14
C.....and compute CPP projection

      IF (ICS.EQ.2) THEN
         DO I = 1,NP
            READ(14,*) JKI,SLAM(JKI),SFEA(JKI),DP(JKI)
            IF (JKI.NE.I) THEN
               IF ((NSCREEN.EQ.1).AND.(MYPROC.EQ.0)) WRITE(6,99801)
               WRITE(16,99801)
            ENDIF
            SLAM(JKI) = DEG2RAD*SLAM(JKI)
            SFEA(JKI) = DEG2RAD*SFEA(JKI)
            CALL CPP(X(JKI),Y(JKI),SLAM(JKI),SFEA(JKI),SLAM0,SFEA0)
         ENDDO
      ENDIF

C.....If ICS = 1 set the SFAC equal to unity

      IF(ICS.EQ.1) THEN
         DO I=1,NP
            SFAC(I)=1.0d0
         ENDDO
      ENDIF
      
C.....If ICS = 2 compute SFAC to adjust equations to CPP coordinates

      IF (ICS.EQ.2) THEN
         DO I = 1,NP
            SFAC(I) = COS(SFEA0)/COS(SFEA(I))
         ENDDO
      ENDIF

C.....If wetting and drying will not be used make sure all bathymetric
C.....depths are > or = to H0

      IF ((NOLIFA.EQ.0).OR.(NOLIFA.EQ.1)) THEN
         DO I = 1,NP
            IF (DP(I).LT.H0) DP(I)=H0
         ENDDO
      ENDIF

C.....Read the global connectivity table from unit 14, compute element
C.....areas, check that sufficient accuracy is provided by the code to
C.....handle the input grid, check to make sure correct convention has
C.....been used for inputting the connectivity table

      DO I = 1, NP
         NNEIGH(I) = 0
      ENDDO

      DO I = 1,NE
         READ(14,*) JKI,NHY,NM(JKI,1),NM(JKI,2),NM(JKI,3)
         NNEIGH(NM(JKI,1)) = NNEIGH(NM(JKI,1)) + 1
         NNEIGH(NM(JKI,2)) = NNEIGH(NM(JKI,2)) + 1
         NNEIGH(NM(JKI,3)) = NNEIGH(NM(JKI,3)) + 1
         IF (JKI.NE.I) THEN
            IF ((NSCREEN.EQ.1).AND.(MYPROC.EQ.0)) WRITE(6,99802)
            WRITE(16,99802)
99802       FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',
     &           'INPUT ERROR  !!!!!!!!!',
     &           //,1X,'YOUR ELEMENT NUMBERING IS NOT SEQUENTIAL ',
     &           /,1X,'CHECK YOUR UNIT 14 INPUT FILE CAREFULLY',//)
         ENDIF
         X1 = X(NM(JKI,1))
         X2 = X(NM(JKI,2))
         X3 = X(NM(JKI,3))
         Y1 = Y(NM(JKI,1))
         Y2 = Y(NM(JKI,2))
         Y3 = Y(NM(JKI,3))
         AVGXY = (ABS(X1)+ABS(X2)+ABS(X3)+ABS(Y1)+ABS(Y2)+ABS(Y3))/6.D0
         DIF1R = AVGXY/(((X2-X1)**2+(Y2-Y1)**2)**0.5d0)
         DIF2R = AVGXY/(((X3-X2)**2+(Y3-Y2)**2)**0.5d0)
         DIF3R = AVGXY/(((X3-X1)**2+(Y3-Y1)**2)**0.5d0)
         DIF1R = LOG10(DIF1R)
         DIF2R = LOG10(DIF2R)
         DIF3R = LOG10(DIF3R)
         IF((DIF1R.GT.NPREC).OR.(DIF2R.GT.NPREC).OR.(DIF3R.GT.NPREC))THEN
            IF ((NSCREEN.EQ.1).AND.(MYPROC.EQ.0)) WRITE(6,9898) JKI
            WRITE(16,9898) JKI
 9898       FORMAT(////,1X,'!!!!!!!!!!  WARNING  !!!!!!!!!',
     &           //,1X,'IF THE GRID COORDINATES HAVE 32 BITS ',
     &           '(APPROX 7 DIGITS) OF PRECISION',
     &           /,1X,'A ROBUST MODEL SOLUTION CAN NOT BE GUARANTEED',
     &           'AT ELEMENT NO. ',I10,
     &           //,1X,'MORE PRECISION MUST BE USED IN THE GRID',//)
         ENDIF
         
C.....NOTE: This is 2 times the actual element area (why?)

         AREAS(JKI) = (X1 - X3)*(Y2 - Y3) + (X3 - X2)*(Y1 - Y3)
         IF (AREAS(JKI).LT.0.0) THEN
            IF((NSCREEN.EQ.1).AND.(MYPROC.EQ.0)) WRITE(6,9899) JKI
            WRITE(16,9899) JKI
 9899       FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &           //,1X,'THE CONNECTIVITY FOR ELEMENT ',I6,
     &           '  HAS BEEN INCORRECTLY SPECIFIED ',
     &           /,1X,'CHECK INPUT AND ENSURE THAT COUNTERCLOCKWISE',
     &           ' CONVENTION HAS BEEN USED ',
     &           //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
            STOP
         ENDIF
      ENDDO

C.....If baroclinic 2D run, read in initial density field

C     IF (IM.EQ.100) THEN
C     OPEN(11,FILE=DIRNAME//'/'//'fort.11')
C     READ(11,*)
C     READ(11,*)
C     READ(11,*) NP2
C     IF (NP2.NE.NP) THEN
C     IF ((NSCREEN.EQ.1).AND.(MYPROC.EQ.0)) WRITE(6,9943)
C     WRITE(16,9943)
C99   43     FORMAT(////,' !!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
C     &              //,' THE NUMBER OF NODES (NP2) IN THE BAROCLINIC',
C     &                 ' INITIAL CONDITION FILE (UNIT 11) ',
C     &               /,' MUST EQUAL THE NUMBER OF NODES (NP) IN ',
C     &                 'THE ADCIRC GRID FILE (UNIT 14)'
C     &              //,' !!!!! EXECUTION WILL NOW BE TERMINATED !!!!!')
C     STOP
C     ENDIF
C     
C     DO I = 1,NP
C     READ(11,*) JKI,DASIGT(JKI),DATEMP(JKI),DASAL(JKI)
C     ENDDO
C     CLOSE(11)
C     ENDIF

! CJT ***************** CALCULATE LOCAL CFL TIME-STEPS *************************

#ifdef CMPI
      CALL ALLOC_MULTI_RATE(); 
      OPEN(366,FILE=DIRNAME//'/'//'EL.TS')
      READ(366,*) Ndom,(Msteps(j),j=1,Ndom),(DtDom(j),j=1,Ndom)
      do j=1,ne
         read(366,*) i,TsGroup(j)
      enddo
      Nelem(:)=0
      ElemNum(:,:)=0
      do i=1,Ndom
         do j=1,ne
            if (TsGroup(j).eq.i) then
               Nelem(i)=Nelem(i)+1
               ElemNum(i,Nelem(i))=j
            endif
         enddo
      enddo
      if (Lts_flag==0) then
         do i=1,Ndom
            dtDom(i)=max(dtDom(1),dt_min0)
            dtDom(i)=min(dtDom(i),dtdp)
            Msteps(i)=1
         enddo
      endif
#else
      CALL ALLOC_MULTI_RATE(); 
      CALL cfl_init(NE,NP,ics,X,Y,dp,NM,SLAM0,SFEA0,DtLoc,slopeflag,myproc)
      write(*,*) 'back from computing cfl '
      dt_max = MAXVAL(Dtloc)
      dt_min = MINVAL(Dtloc)
      dt_min = max(dt_min,dt_min0)
      Mint = dt_max/dt_min;
      WRITE(*,*) "Before Scale ---------------------------------------"
      WRITE(*,*) "Largest stable time-step over grid from CFL: ",dt_max
      WRITE(*,*) "Smallest stable time-step over grid from CFL: ",dt_min
      WRITE(*,*) "Largest to smallest ratio: ",Mint


      ! Scale Local time-steps to user input
      !IF (dtmax > DTDP) THEN
      DO i=1,NE
!        DtLoc(i) = min((DtLoc(i)/dt_max)*DT,dtloc(i))
        DtLoc(i) = min(DTDP,dtloc(i))
        WRITE(69,*) DtLoc(i)
      ENDDO
      dt_max = MAXVAL(Dtloc);
      dt_min = MINVAL(Dtloc);
      dt_min = max(dt_min,dt_min0)
      Mint = dt_max/dt_min;
      WRITE(*,*) "After Scale ---------------------------------------"
      WRITE(*,*) "Largest stable time-step over grid from CFL: ",dt_max
      WRITE(*,*) "Smallest stable time-step over grid from CFL: ",dt_min
      WRITE(*,*) "Largest to smallest ratio: ",Mint,Mbar

!  CND
      if (lts_flag /= 1) then
         Mint=1
      endif
!  
      if (Mint==1) then
         Ndom=1
         DtDom(1)=dt_min
         Msteps(1)=1
      else
         Ndom=MaxNdom
         do i=1,Ndom-1
            Msteps(i)=Mint**(1./(Ndom-1))
            Msteps(i)=min(Msteps(i),Mbar)
         enddo
         Msteps(Ndom)=1
         write(*,*) (Msteps(i),i=1,Ndom)
         if (Msteps(1)<2) then
            Ndom=1
            Dtdom(1)=dt_min
         else
            DtDom(1) = dt_min
            do i=2,Ndom
               DtDom(i)=Msteps(i-1)*dtdom(i-1)
            enddo
         endif
      endif
!!!!
      Nelem(:) = 0
      TsGroup(:)=1
      if (nDom > 1) then
         DO i=1,NE
            DO j=2,nDom
               IF ( (DtLoc(i)>=DtDom(j-1)) .and. (DtLoc(i)<DtDom(j)) ) THEN
                  DtLoc(i) = DtDom(j-1)
                  TsGroup(i) = j-1
                  Nelem(j-1) = Nelem(j-1) + 1
                  ElemNum(j-1,Nelem(j-1)) = i
               Endif
            enddo
            if (Dtloc(i) >= DtDom(Ndom)) then
               dtloc(i)=dtdom(Ndom)
               TsGroup(i)=Ndom
               Nelem(Ndom) = Nelem(Ndom) + 1
               ElemNum(Ndom,Nelem(Ndom)) = i
            endif
         ENDDO
      else
         DO i=1,NE
            DtLoc(i) = DtDom(1)
            TsGroup(i) = 1
            Nelem(1) = Nelem(1) + 1
            ElemNum(1,Nelem(1)) = i
         enddo
      endif
#endif        

      OPEN(69,FILE=DIRNAME//'/'//'elementTSGroups.out')
      DO i=1,NE
        WRITE(69,'(2i10,f20.10)')
     $        i,TsGroup(i),DtDom(TsGroup(i))
      ENDDO

      WRITE(69,*) "Number of Domains: ",Ndom,(Nelem(i),i=1,Ndom)
      CLOSE(69)

      IF (SUM(Nelem) .ne. NE) THEN
	WRITE(*,*) "ERROR: Local time-step binning error (Adcirc.f)"
	WRITE(*,*) "Total Number of Elements Binned: ",SUM(Nelem),"  NE: ",NE
	STOP
      ENDIF
 
 
!      ENDIF
!  Reset DTDP to maximum time step
      DTDP=DtDom(Ndom)
      DT=DTDP

C.....Compute total number of time steps NT

      NT = INT(RNDAY*(86400.D0/DTDP) + 0.5D0)
      WRITE(16,1920) NT
1920  FORMAT(5X,'NUMBER OF TIME STEPS  =',I8,/)

C.....Process startdry info from unit 12 (if NOLIFA = 3 -> NSTARTDRY = 1
C.....STARTDRY now set in nodal attributes

      IF ((NSTARTDRY.EQ.1).AND.(NWP.EQ.0)) THEN
         
         ALLOCATE(STARTDRY(NP))

C.....Open unit 12 file

         OPEN(12,FILE=DIRNAME//'/'//'fort.12')

C.....Read startdry info from unit 12

         READ(12,'(A24)') AGRID2
         READ(12,*) NE2,NP2

C.....Check that NE2 and NP2 mathe with grid file

         IF ((NE2.NE.NE).OR.(NP2.NE.NP)) THEN
            IF ((NSCREEN.EQ.1).AND.(MYPROC.EQ.0)) WRITE(6,9900)
            WRITE(16,9900)
 9900       FORMAT(////,1X,'!!!!!!!!!!  FATAL ERROR  !!!!!!!!!',
     &           //,1X,'THE PARAMETER NE2 AND NP2 MUST MATCH NE AND NP ',
     &           /,1X,'USER MUST CHECK FORT.12 INPUT FILE ',
     &           //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
            STOP
         ENDIF

C.....Read in startdry code values

         DO I = 1,NP
            READ(12,*) JKI,DUM1,DUM2,STARTDRY(JKI)
            IF (MODAL_IC.NE.3) THEN
               IF (STARTDRY(JKI).EQ.-88888) THEN
                  STARTDRY(JKI) = 1
               ELSE
                  STARTDRY(JKI) = 0
               ENDIF
            ENDIF
            IF (JKI.NE.I) THEN
               IF ((NSCREEN.EQ.1).AND.(MYPROC.EQ.0)) WRITE(6,99805)
               WRITE(16,99805)
99805          FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',
     &              'INPUT ERROR  !!!!!!!!!',
     &              //,1X,'YOUR NODE NUMBERING IS NOT SEQUENTIAL ',
     &              'CHECK YOUR UNIT 12 INPUT FILE CAREFULLY',//)
            ENDIF
         ENDDO

C.....Close unite 12 file

         CLOSE(12)

C.....If not using startup elevation file

      ENDIF


C.....Reset tau0var values based on input values of startdry and
C.....automatic selection of tau0 on a processor

C     DO I = 1,NP
C     IF (TAU0.LT.0.D0) THEN
C     IF (DP(I).LE.10.D0) TAU0VAR(I) = 0.02D0
C     IF (DP(I).GT.10.D0) TAU0VAR(I) = 0.005D0
C     IF (STARTDRY(I).EQ.-77777) TAU0VAR(I) = 0.02D0
C     IF (STARTDRY(I).EQ.-88888) TAU0VAR(I) = 0.02D0
C     WRITE(16,248) MYPROC,I,TAU0VAR(I)
C     248       FORMAT(/,' myproc = ',I6,' node = ',I8,
C     &             ' tau0 set to ',F12.6,/)
C     ELSE
C     C          TAU0VAR(I) = TAU0
C     ENDIF
C     ENDDO

C.....Output to unit 16 grid information including AGRID, NE, NP, H0,
C.....and nodal coordinates and bathymetry

      WRITE(16,2039) AGRID
 2039 FORMAT(/,5X,'GRID IDENTIFICATION : ',A24,/)
      IF(NSTARTDRY.EQ.1) WRITE(16,2038) AGRID2
 2038 FORMAT(5X,'STARTDRY FILE IDENTIFICATION : ',A24,/)
      WRITE(16,3) NP
    3 FORMAT(5X,'TOTAL NUMBER OF NODES =',I6,/)
      WRITE(16,4) NE
    4 FORMAT(5X,'TOTAL NUMBER OF ELEMENTS =',I6,/)
      IF(ICS.EQ.2) WRITE(16,13) SLAM0*RAD2DEG,SFEA0*RAD2DEG
 13   FORMAT(5X,'LONGITUDE ABOUT WHICH CPP PROJECTION IS CENTERED',
     &     '  SLAM0 = ',F9.4,' DEGREES',
     &     /,5X,'LATITUDE  ABOUT WHICH CPP PROJECTION IS CENTERED',
     &     '  SFEA0 = ',F9.4,' DEGREES',/)
      IF (NSTARTDRY.EQ.0) THEN
         IF (NABOUT.NE.1) THEN
            WRITE(16,24)
 24         FORMAT(/,1X,'NODAL COORDINATES AND BATHYMETRY :')
            IF (ICS.EQ.1) THEN
               IF ((NTIP.EQ.0).AND.(NCOR.EQ.0)) THEN
                  WRITE(16,25)
 25               FORMAT(/,10X,'NODE NO.',10X,'X',20X,'Y',15X,'DP',/)
                  DO I = 1,NP
                     WRITE (16,2008) I,X(I),Y(I),DP(I)
 2008                FORMAT(5X,I6,2(2X,F20.2),2X,F12.2)
                  ENDDO
               ELSE
                  WRITE(16,9195)
 9195             FORMAT(/,1X,'   NODE ',7X,'X',14X,'Y',9X,
     &                 'LAMBDA(DEG)',6X,'FEA(DEG)',9X,'DP',/)
                  DO I = 1,NP
                     WRITE (16,9197) I,X(I),Y(I),SLAM(I)*RAD2DEG,
     &                    SFEA(I)*RAD2DEG,DP(I)
 9197                FORMAT(1X,I6,2(1X,F14.1),1X,2(1X,E15.7),1X,F8.2)
                  ENDDO
               ENDIF
            ELSE
               WRITE(16,9225)
 9225          FORMAT(/,1X,'   NODE ',2X,'LAMBDA(DEG)',5X,'FEA(DEG)',11X,
     &              'XCP',14X,'YCP',11X,'DP',/)
               DO I = 1,NP
                  WRITE (16,9228) I,SLAM(I)*RAD2DEG,SFEA(I)*RAD2DEG,
     &                 X(I),Y(I),DP(I)
 9228             FORMAT(1X,I6,2(1X,F14.8),2(1X,F15.1),1X,F10.2)
               ENDDO
            ENDIF
         ELSE
            WRITE(16,3511)
 3511       FORMAT(/,5X,'NODAL COORDINATES AND BATHYMETRY',
     &           ' INFORMATION IS AVAILABLE IN THE',
     &           /,6X,'UNIT 14 INPUT FILE')
         ENDIF
      ELSE
         IF (NABOUT.NE.1) THEN
            WRITE(16,24)
            IF (ICS.EQ.1) THEN
               IF ((NTIP.EQ.0).AND.(NCOR.EQ.0)) THEN
                  WRITE(16,3527)
 3527             FORMAT(/,10X,'NODE NO.',10X,'X',20X,'Y',15X,'DP',
     &                 5X,'STARTDRY',/)
                  DO I = 1,NP
                     IF (STARTDRY(I).EQ.-88888.D0) THEN
                        WRITE (16,3529) I,X(I),Y(I),DP(I),STARTDRY(I)
 3529                   FORMAT(5X,I6,2(2X,F20.2),2X,F12.2,2X,F12.0)
                     ELSE
                        WRITE (16,2008) I,X(I),Y(I),DP(I)
                     ENDIF
                  ENDDO
               ELSE
                  WRITE(16,3530)
 3530             FORMAT(/,1X,'   NODE ',7X,'X',14X,'Y',9X,
     &                 'LAMBDA(DEG)',6X,'FEA(DEG)',9X,'DP',
     &                 5X,'STARTDRY',/)
                  DO I = 1,NP
                     IF (STARTDRY(I).EQ.-88888.D0) THEN
                        WRITE (16,3531) I,X(I),Y(I),SLAM(I)*RAD2DEG,
     &                       SFEA(I)*RAD2DEG,DP(I),STARTDRY(I)
 3531                   FORMAT(1X,I6,2(1X,F14.1),1X,2(1X,E15.7),1X,F8.2,
     &                       1X,F10.0)
                     ELSE
                        WRITE (16,9197) I,X(I),Y(I),SLAM(I)*RAD2DEG,
     &                       SFEA(I)*RAD2DEG,DP(I)
                     ENDIF
                  ENDDO
               ENDIF
            ELSE
               WRITE(16,3535)
 3535          FORMAT(/,1X,'   NODE ',2X,'LAMBDA(DEG)',5X,'FEA(DEG)',11X,
     &              'XCP',14X,'YCP',11X,'DP',
     &              5X,'STARTDRY',/)
               DO I = 1,NP
                  IF (STARTDRY(I).EQ.-88888.D0) THEN
                     WRITE (16,3537) I,SLAM(I)*RAD2DEG,SFEA(I)*RAD2DEG,
     &                    X(I),Y(I),DP(I),STARTDRY(I)
 3537                FORMAT(1X,I6,2(1X,F14.8),2(1X,F15.1),1X,F10.2,2X,F10.0)
                  ELSE
                     WRITE (16,9228) I,SLAM(I)*RAD2DEG,SFEA(I)*RAD2DEG,
     &                    X(I),Y(I),DP(I)
                  ENDIF
               ENDDO
            ENDIF
         ELSE
            WRITE(16,3540)
 3540       FORMAT(/,5X,'NODAL COORDINATES AND BATHYMETRY',
     &           ' INFORMATION IS AVAILABLE IN THE',
     &           /,6X,'UNIT 14 AND 12 INPUT FILES')
         ENDIF
      ENDIF

C.....Output to unit 16 the global connectivity table (node numbers for
C.....elements

      IF (NABOUT.NE.1) THEN
         WRITE(16,26)
 26      FORMAT(//,5X,'GLOBAL NODE NUMBERS FOR EACH ELEMENT :')
         WRITE(16,27)
 27      FORMAT(/,9X,'ELEMENT',8X,'N1',9X,'N2',10X,'N3',/)
         DO I = 1,NE
            WRITE(16,2009) I,NM(I,1),NM(I,2),NM(I,3)
 2009       FORMAT(8X,4(I7,4X))
         ENDDO
      ELSE
         WRITE(16,3512)
 3512    FORMAT(/,5X,'THE GLOBAL CONNECTIVITY TABLE',
     &        ' INFORMATION IS AVAILABLE IN THE',
     &        /,6X,'UNIT 14 INPUT FILE')
      ENDIF
      
C.....Read information concerning bottom friction coefficient
C.....If NWP = 1, input nodal friction coefficients from unit 21
C.....If NWP = 2, set nodal friction coefficients equal to Cf
C.....If NWP = 3

C...  
C...  READ INFORMATION CONCERNING BOTTOM FRICTION COEFFICIENT
C...  IF NWP=1, INPUT NODAL FRICTION COEFFICIENTS FROM UNIT 21
C...  IF NWP=0, SET NODAL FRICTION COEFFICIENTS EQUAL TO CF
C...  IF NWP=2, READ ADDITIONAL FRICTIONAL PARAMETERS FOR BRIDGE PILINGS
C...  
      WRITE(16,1112)
      WRITE(16,2045)
 2045 FORMAT(//,' BOTTOM FRICTION INFORMATION',//)

      HBREAK=1.
      FTHETA=1.
      FGAMMA=1.
      IF(NOLIBF.EQ.0) READ(15,*) TAU
      CF=TAU
      IF(NOLIBF.EQ.1) READ(15,*) CF
      IF(NOLIBF.EQ.2) READ(15,*) CF,HBREAK,FTHETA,FGAMMA

      IF (NWP.EQ.0) THEN
         ALLOCATE(FRIC(NP))
         DO I=1,NP
            FRIC(I)=CF
         END DO
         IF(NOLIBF.EQ.2) THEN
            WRITE(16,101) CF,HBREAK,FTHETA,FGAMMA
 101        FORMAT(5X,'HYBRID FRICTION RELATIONSHIP PARAMTERS, CFMIN =',
     &           F12.8,'  HBREAK = ',F8.2,
     &           /,5X,'FTHETA = ',F8.2,'  FGAMMA = ',F10.4,//)
         ENDIF
         IF(NOLIBF.EQ.1) THEN
            WRITE(16,8) CF
 8          FORMAT(5X,'NONLINEAR FRICTION COEFFICIENT CF =',F12.8,/)
         ENDIF
         IF(NOLIBF.EQ.0) THEN
            WRITE(16,106) TAU
 106        FORMAT(5X,'LINEAR BOTTOM FRICTION TAU =',F12.8,5X,'1/sec'/)
            IF(TAU.NE.TAU0) THEN !CHECK TAU VALUE AGAINST TAU0
               IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9951)
               WRITE(16,9951)
 9951          FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',
     &              'INPUT ERROR  !!!!!!!!!',
     &              //,1X,'TYPICALLY YOUR INPUT VALUE FOR ',
     &              'TAU0 SHOULD BE SET EQUAL TO TAU')
            ENDIF
         ENDIF
      ENDIF


C     IF(NWP.EQ.1) THEN
C     OPEN(21,FILE=DIRNAME//'/'//'fort.21')
C     READ(21,'(A20)') AFRIC
C     DO I=1,NP
C     READ(21,*) NHG,FRIC(NHG)
C     IF(NHG.NE.I) THEN
C     IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,99803)
C     WRITE(16,99803)
C99   803       FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ',
C     &                     'INPUT ERROR  !!!!!!!!!',
C     &        //,1X,'YOUR NODAL FRICTION NUMBERING IS NOT SEQUENTIAL ',
C     &        /,1X,'CHECK YOUR UNIT 21 INPUT FILE CAREFULLY',
C     &        //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
C     STOP
C     ENDIF
C     END DO
C     WRITE(16,3601) AFRIC
C     3601   FORMAT(/,5X,'FRICTION FILE IDENTIFICATN : ',A20,/)
C     IF(NABOUT.NE.1) THEN
C     WRITE(16,2080)
C     2080     FORMAT(/,10X,'NODE',5X,'BOTTOM FRICTION FRIC',5X,/)
C     DO I=1,NP
C     WRITE(16,2087) I,FRIC(I)
C     2087       FORMAT(7X,I6,6X,E15.10)
C     END DO
C     ELSE
C     WRITE(16,3504)
C     3504     FORMAT(/,5X,'NODAL BOTTOM FRICTION VALUES ARE AVAILABLE',
C     &           /,6X,' IN UNIT 21 INPUT FILE')
C     ENDIF
C     ENDIF

C     IF(NWP.EQ.2) THEN
C     CALL ALLOC_MAIN13()   !allocate bridge piling arrays
C     DO I=1,NP
C     NBNNUM(I)=0
C     BK(I)=0.d0
C     BALPHA(I)=0.d0
C     BDELX(I)=1.d0
C     ENDDO
C     OPEN(21,FILE=DIRNAME//'/'//'fort.21')
C     READ(21,'(A20)') AFRIC
C     READ(21,*) NBPNODES
C     DO I=1,NBPNODES
C     READ(21,*) NBNNUM(I),BK(I),BALPHA(I),BDELX(I),POAN
C     BDELX(I)=4.d0*BDELX(I)/POAN
C     ENDDO
C     WRITE(16,3602) AFRIC
C     3602   FORMAT(/,5X,'BRIDGE PIER FRICTION FILE IDENTIFICATN : ',A20,/)
C     IF(NABOUT.NE.1) THEN
C     WRITE(16,2081)
C     2081     FORMAT(/,10X,'NODE',3X,'PIER SHAPE FACTOR',3X,
C     &                 'CONSTRICTION FRACTION',3X,'EFFECTIVE DELX'/)
C     DO I=1,NBPNODES
C     WRITE(16,2082) NBNNUM(I),BK(I),BALPHA(I),BDELX(I)
C     2082       FORMAT(5X,I8,10X,F7.3,12X,F7.3,13X,F9.3)
C     END DO
C     ELSE
C     WRITE(16,2083)
C     2083     FORMAT(/,5X,'BRIDGE PILING FRICTION VALUES ARE AVAILABLE',
C     &           /,6X,' IN UNIT 21 INPUT FILE')
C     ENDIF
C     ENDIF

C...  
C...  READ IN AND WRITE OUT EDDY VISCOSITY/DIFFUSIVITY COEFFICIENTS
C...  


      IF (IM.EQ.10) THEN
         READ(15,*) ESLM,ESLC
         DO I=1,NP
            EVM(I)=ESLM
            EVC(I)=ESLC
         END DO
         WRITE(16,111) ESLM,ESLC
 111     FORMAT(5X,'EVM, EDDY VISCOSITY COEFFICIENT =',E15.8,/,
     &        5X,'EVC, EDDY DIFFUSIVITY COEFFICIENT =',E15.8,//)
      ELSE
         READ(15,*) ESLM
         IF (NWP.EQ.0) THEN
            ALLOCATE(EVM(NP))
            DO I=1,NP
               EVM(I)=ESLM
            END DO
            WRITE(16,11) ESLM
 11         FORMAT(5X,'EVM, EDDY VISCOSITY COEFFICIENT =',E15.8,//)
         ENDIF
      ENDIF


C...  02/19/2007 s.b.
      EVMSUM = 0.D0
      IF (NWP.EQ.0) THEN
         DO I=1,NP
            EVMSUM = EVMSUM + ABS(EVM(I))
         ENDDO
      ENDIF
      
C     
C     ek: Initialize nodal attributes, now that grid has been read
C     in from unit 14 file.

      IF (NWP.GT.0)
     &     CALL InitNodalAttr(DP, NP, G, NScreen, ScreenUnit,MyProc,NAbOut)
      


C...  
C...  READ CORIOLIS INFORMATION AND COMPUTE THE CORIOLIS VECTOR
C...  OUTPUT RESULTING CORIOLIS INFORMATION
C...  
      WRITE(16,1112)
      WRITE(16,2090)
 2090 FORMAT(//,1X,'CORIOLIS INFORMATION ',//)

      READ(15,*) CORI
      IF(NCOR.EQ.0) THEN
         DO I=1,NP
            CORIF(I)=CORI
         END DO
      ENDIF
      IF(NCOR.EQ.1) THEN
         DO I=1,NP
            CORIF(I)=2.0d0*7.29212d-5*SIN(SFEA(I))
         END DO
      ENDIF

      IF(NCOR.EQ.0) THEN
         WRITE(16,12) CORI
 12      FORMAT(5X,'CONSTANT CORIOLIS COEFFICIENT =',E15.8,5X,'1/SEC',/)
      ENDIF
      IF(NCOR.EQ.1) THEN
         WRITE(16,3604)
 3604    FORMAT(/,5X,'LATITUDES ARE USED TO COMPUTE VARIABLE CORIOLIS',
     &        /,7X,'AND ARE BASED ON INPUT NODAL COORDINATES',/)
         IF(NABOUT.NE.1) THEN
            WRITE(16,2092)
 2092       FORMAT(/,10X,' NODE ',5X,'NODAL CORIOLIS CORIF',/)
            DO I=1,NP
               WRITE(16,2096) I,CORIF(I)
 2096          FORMAT(7X,I6,10X,E15.9)
            END DO
         ENDIF
      ENDIF

C...  
C...  READ AND PROCESS INFORMATION ABOUT THE TIDAL POTENTIAL CONSTITUENTS
C...  

      READ(15,*) NTIF
      mntif = ntif
      if (ntif .eq. 0) mntif = 1

C.... allocate tidal potential arrays

      call alloc_main4a()

C.... READ TIDAL POTENTIAL AMPLITUDE, FREQUENCIES, NODAL FACTORS,
C.... EQUILIBRIUM ARGUMENTS AND ALPHANUMERIC LABEL
C.... 
      DO I=1,NTIF
         READ(15,'(A5)')  TIPOTAG(I)
         READ(15,*)  TPK(I),AMIGT(I),ETRF(I),FFT(I),FACET(I)
         IF(AMIGT(I).EQ.0.) THEN
            PERT(I)=0.
         ELSE
            PERT(I)=2.D0*PI/AMIGT(I)
         ENDIF
      END DO

C...  LINES TO USE EARTH LOAD/SELF-ATTRACTION PART OF TIDAL POTENTIAL FORCING

      CALL ALLOC_MAIN4b()
      IF(NTIP.EQ.2) THEN
         OPEN(24,FILE='fort.24')
         DO I=1,NTIF
            READ(24,9930)
 9930       FORMAT(///)
            DO J=1,NP
               READ(24,*) JJ,SALTAMP(I,JJ),SALTPHA(I,JJ)
               SALTPHA(I,JJ)=SALTPHA(I,JJ)*DEG2RAD
            END DO
         END DO
      ELSE
         DO I=1,NTIF
            DO J=1,NP
               SALTAMP(I,J)=0.d0
               SALTPHA(I,J)=0.d0
            END DO
         END DO
         CLOSE(24)
      ENDIF

C...  
C...  OUTPUT TO UNIT 16 INFORMATION ABOUT TIDAL POTENTIAL FORCING
C.... OUTPUT WILL VARY DEPENDING ON VALUES OF NTIP,NTIF AND NCOR
C...  
      WRITE(16,1112)
      WRITE(16,2102)
 2102 FORMAT(//,1X,'TIDAL POTENTIAL FORCING INFORMATION ',//)
      WRITE(16,22) NTIF
 22   FORMAT(/,1X,'TIDAL POTENTIAL IS FORCED FOR ',I5,
     &     ' CONSTITUENT(S) ')
      IF(NTIF.GT.0) WRITE(16,23)
 23   FORMAT(/,1X,'AMPLITUDE',4X,'FREQUENCY',5X,
     &     '    ETRF      ','NODAL FACTOR',2X,
     &     'EQU.ARG(DEG)',1X,'CONSTITUENT',/)
      DO I=1,NTIF
         WRITE(16,2107) TPK(I),AMIGT(I),ETRF(I),FFT(I),FACET(I),
     &        TIPOTAG(I)
 2107    FORMAT(1X,F9.7,1X,F15.12,2X,F10.7,5X,F10.7,1X,F10.3,7X,A5)
      END DO
C...  
C...  CONVERT FACET(I) VALUES FROM DEGREES TO RADIANS
C...  
      DO I=1,NTIF
         FACET(I)=FACET(I)*DEG2RAD
      END DO
C...  
C...  CHECK CONSISTENCY OF INPUT PARAMETERS NTIF AND NTIP
C...  
      IF(((NTIP.EQ.0).AND.(NTIF.NE.0)).OR.((NTIP.NE.0).AND.
     &     (NTIF.EQ.0))) THEN
         IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9961)
         WRITE(16,9961)
 9961    FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',
     &        'INPUT ERROR  !!!!!!!!!',
     &        //,1X,'YOUR SELECTION OF NTIF AND NTIP (UNIT 15 INPUT ',
     &        'PARAMETERS) IS INCONSISTENT',
     &        /,1X,'PLEASE CHECK THESE VALUES')
         IF(NFOVER.EQ.1) THEN
            IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9987)
            WRITE(16,9987)
 9987       FORMAT(/,1X,'PROGRAM WILL OVERRIDE THE SPECIFIED ',
     &           'INPUT AND NEGLECT TIDAL POTENTIAL TERMS',
     &           /,1X,' AND/OR RESET NTIP = 0',
     &           //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)
            NTIP=0
         ELSE
            IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9973)
            WRITE(16,9973)
            STOP
         ENDIF
         GOTO 1893
      ENDIF
C...  
C...  PRINT OUT LAT/LON VALUES TO BE USED IN COMPUTING TIDAL POTENTIAL
C.... IF NOT ALREADY DONE SO IN CORIOLIS SECTION AND TIDAL POTENTIAL
C.... IS ACTIVATED WITH NTIP=1
C...  

      IF(NTIP.GE.1) THEN
         IF(ICS.EQ.1) THEN
            WRITE(16,3605)
 3605       FORMAT(/,5X,'LONGITUDES AND LATITUDES ARE USED TO',
     &           ' COMPUTE THE TIDAL POTENTIAL FUNCTION',
     &           /,7X,'AND ARE BASED ON AN INVERSE CPP PROJECTION ',
     &           'OF THE INPUT COORDINATES',/)
         ELSE
            WRITE(16,2109)
 2109       FORMAT(/,5X,'LONGITUDES AND LATITUDES ARE USED TO',
     &           ' COMPUTE THE TIDAL POTENTIAL FUNCTION',
     &           /,7X,'AND ARE BASED ON INPUT NODAL COORDINATES ',/)
         ENDIF
      ENDIF

C...  
C...  INPUT FROM UNIT 15 THE TIDAL FORCING FREQUENCIES ON THE ELEVATION
C.... SPECIFIED BOUNDARIES: INCLUDING NBFR, FREQUENCIES, NODAL FACTORS,
C.... EQUILIBRIUM ARGUMENTS AND AN ELEVATION BOUNDARY CONDITION
C.... ALPHANUMERIC DESCRIPTOR
C...  
 1893 READ(15,*) NBFR
      MNBFR = NBFR
      IF (NBFR.EQ.0) MNBFR = 1

C     - Allocate arrays dimensioned by MNBFR

      call alloc_main5()

      WRITE(16,1112)
      WRITE(16,2106)
 2106 FORMAT(//,1X,'ELEVATION SPECIFIED BOUNDARY FORCING INFORMATION '
     &     ,//)
      WRITE(16,20) NBFR
 20   FORMAT(/,5X,'NUMBER OF PERIODIC, ELEVATION SPECIFIED ',
     &     'CONSTITUENTS =',I5)
      IF(NBFR.GE.1) WRITE(16,21)
 21   FORMAT(/,7X,'CONSTITUENT #',4X,'FREQUENCY',4X,'NODAL FACTOR',
     &     3X,'EQU.ARG (DEG)',2X,'CONSTITUENT',/)
      DO I=1,NBFR
         READ(15,'(A5)') BOUNTAG(I)
         READ(15,*) AMIG(I),FF(I),FACE(I)
         WRITE(16,1850) I,AMIG(I),FF(I),FACE(I),BOUNTAG(I)
 1850    FORMAT(12X,I2,6X,F16.12,2X,F10.7,2X,F10.3,10X,A5)
         FACE(I)=FACE(I)*DEG2RAD
         IF(AMIG(I).EQ.0.) THEN
            PER(I)=0.
         ELSE
            PER(I)=2.D0*PI/AMIG(I)
         ENDIF
      ENDDO
C...  
C...  INPUT ELEVATION BOUNDARY FORCING NODE NUMBER INFORMATION FROM UNIT 14 AND
C.... OUTPUT TO UNIT 16
C...  
C...  INPUT THE TOTAL NUMBER OF ELEVATION BOUNDARY SEGMENTS
C...  
      READ(14,*) NOPE

      WRITE(16,1852) NOPE
 1852 FORMAT(///,5X,'TOTAL NUMBER OF ELEVATION BOUNDARY FORCING',
     &     ' SEGMENTS ',' = ',I5)
C...  
C...  INPUT THE TOTAL NUMBER OF ELEVATION BOUNDARY NODES
C...  
      READ(14,*) NETA
      WRITE(16,1854) NETA
 1854 FORMAT(/,5X,'TOTAL NUMBER OF ELEVATION SPECIFIED BOUNDARY NODES ='
     &     ,I6)

C     allocate arrays dimensioned by NOPE and NETA
      MNOPE = NOPE
      IF (NOPE.EQ.0) MNOPE = 1
      MNETA = NETA
      IF (NETA.EQ.0) MNETA = 1

      call alloc_main2()     
C...  
C...  INPUT THE NODE NUMBERS ON EACH ELEVATION BOUNDARY FORCING SEGMENT
C...  
      MNEI=0
      JNMM=0
      DO K=1,NOPE
         READ(14,*) NVDLL(K)
         WRITE(16,281) K,NVDLL(K)
 281     FORMAT(//,5X,'TOTAL NUMBER OF NODES ON ELEVATION SPECIFIED ',
     &        'BOUNDARY SEGMENT ',2X,I2,2X,'=',1X,I5,/)
         DO I=1,NVDLL(K)
            READ(14,*) NBDV(K,I)
            WRITE(16,1855) NBDV(K,I)
 1855       FORMAT(7X,I7)
            IF (NNEIGH(NBDV(K,I)).NE.0) THEN
               NNEIGH(NBDV(K,I))=NNEIGH(NBDV(K,I))+1
               IF (NNEIGH(NBDV(K,I)).GT.MNEI) MNEI=NNEIGH(NBDV(K,I))
               NNEIGH(NBDV(K,I)) = 0
            ENDIF

            NBD(JNMM+I)=NBDV(K,I)
         ENDDO
         JNMM=JNMM+NVDLL(K)
      ENDDO
C...  
C...  CHECK TO MAKE SURE THAT JNMM EQUALS NETA
C...  
      IF(NETA.NE.JNMM) THEN
         IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9945)
         WRITE(16,9945)
 9945    FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL INPUT ERROR ',
     &        '!!!!!!!!!',
     &        //,1X,'THE INPUT PARAMETER NETA FROM UNIT 14 DOES NOT MATCH ',
     &        'THE TOTAL NUMBER OF BOUNDARY NODES',
     &        /,1X,' FROM ALL THE SPECIFIED SEGMENTS COMBINED')
         IF(NFOVER.EQ.1) THEN
            IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9989)
            WRITE(16,9989)
 9989       FORMAT(/,1X,'THE PROGRAM WILL NOW CORRECT THIS ERROR',
     &           /,1X,'PLEASE CHECK YOUR INPUT CAREFULLY !!!',
     &           //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)
            NETA=JNMM
         ELSE
            IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9973)
            WRITE(16,9973)
            STOP
         ENDIF
      ENDIF
C...  
C...  SET UP TO READ IN TIME SERIES ELEVATION SPECIFIED BOUNDARY CONDITIONS IF APPROPRIATE
C...  
      IF((NBFR.EQ.0).AND.(NOPE.GT.0)) THEN
         WRITE(16,1871)
 1871    FORMAT(/,5X,'TIME SERIES ELEVATION SPECIFIED VALUES WILL BE ',
     &        'READ FROM UNIT 19',
     &        /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE ',
     &        /,9X,'ELEVATION DATA WITH THE MODEL TIME STEP.')
      ENDIF

C...  
C...  INPUT FORCING CONDITIONS ON PERIODIC ELEVATION SPECIFIED BOUNDARIES FOR EACH
C...  OF THE ELEVATION FORCING FREQUENCIES FROM UNIT 15 AND OUTPUT TO UNIT 16
C...  

      DO I=1,NBFR
         WRITE(16,29) I,BOUNTAG(I)
 29      FORMAT(////,5X,'ELEVATION BOUNDARY TIDAL FORCING FOR',
     &        ' CONSTITUENT NUMBER',I4,1X,'DESIGNATED : ',A5)
         READ(15,'(A10)') ALPHA
         WRITE(16,31) ALPHA
 31      FORMAT(9X,'VERIFICATION OF CONSTITUENT : ',A10,/)
         WRITE(16,30)
 30      FORMAT(14X,'NODE',11X,'AMPL.',9X,'PHASE(DEG)',/)
         DO J=1,NETA
            READ(15,*) EMO(I,J),EFA(I,J)
            WRITE(16,1870) NBD(J),EMO(I,J),EFA(I,J)
 1870       FORMAT(10X,I8,4X,F14.5,4X,F12.3)
            EFA(I,J) = EFA(I,J)*DEG2RAD
         ENDDO
C     DO J=1,NETA
C     READ(15,*) UMO(I,J),UFA(I,J)
C     UFA(I,J) = UFA(I,J)*DEG2RAD
C     ENDDO
C     DO J=1,NETA
C     READ(15,*) VMO(I,J),VFA(I,J)
C     VFA(I,J) = VFA(I,J)*DEG2RAD
C     ENDDO
      ENDDO

C.....READ THE MINIMUM INNER ANGLE FOR WHICH VELOCITY AT FLOW BOUNDARY NODES
C.....WILL BE ZEROED IN THE TANGENTIAL DIRECTIONS WHEN NORMAL FLOW IS AN
C.....ESSENTIAL B.C.

      READ(15,*) ANGINN
      WRITE(16,1112)
      WRITE(16,7654) ANGINN
 7654 FORMAT(//,5X,'ANGINN = ',F8.2,' DEGREES',
     &     /,5X,'ALL FLOW BOUNDARY NODES WITH NORMAL FLOW AS AN ',
     &     'ESSENTIAL B.C. AND ',
     &     /,9X,'INNER ANGLES LESS THAN ANGINN WILL HAVE BOTH NORMAL ',
     &     /,9X,'AND TANGENTIAL VELOCITY COMPONENTS ZEROED',/)
      COSTSET=COS(ANGINN*DEG2RAD)

C...  
C...  INPUT FLOW BOUNDARY INFORMATION FROM UNIT 14 AND OUTPUT TO UNIT 16
C...  

C.....INTERIOR NODES, LBCODE=-1, COS=0, SIN=1
C.....BOUNDARY NODES, LBCODE=LBCODEI=IBTYPE,
C.....COS & SIN DETERMINED FROM NORMAL DIRECTION IN ALL CASES, ALTHOUGH THIS
C.......INFORMATION IS ONLY USED WHEN NORMAL FLOW IS AN ESSENTIAL B.C. AND
C.......FREE TANGENTIAL SLIP IS ALLOWED.

C.....INPUT THE TOTAL NUMBER OF FLOW BOUNDARY SEGMENTS

      WRITE(16,1112)
      WRITE(16,1878)
 1878 FORMAT(//,1X,'FLOW BOUNDARY INFORMATION ',/)
      READ(14,*) NBOU

      WRITE(16,1879) NBOU
 1879 FORMAT(//,5X,'THE TOTAL NUMBER OF FLOW BOUNDARY SEGMENTS = ',I5)

C.....INPUT THE TOTAL NUMBER OF FLOW BOUNDARY NODES

      READ(14,*) NVEL
      WRITE(16,1881) NVEL
 1881 FORMAT(/,5X,'THE TOTAL NUMBER OF FLOW BOUNDARY NODES = ',I5)

      MNBOU = NBOU
      IF (NBOU.EQ.0) MNBOU = 1
      MNVEL = NVEL*2            !Cvjp  -  11/28/99 -  upper bound guess for MNVEL

C.....Allocate space for nonperiodic zero and nonzero normal flow boundary arrays
C.....including barriers
      call alloc_main3()

C.....INPUT THE NUMBER OF NODES IN THE NEXT FLOW BOUNDARY SEGMENT
C.....AND THE BOUNDARY TYPE

      JGW=0
      JME=0
      NFLUXF=0
      NFLUXB=0
      NFLUXIB=0
      NFLUXIBP=0
      NVELEXT=0
      NIBSEG = 0
      NEBSEG = 0
      
      CALL ALLOC_EDGES0()

      DO K=1,NBOU
         READ(14,*) NVELL(K),IBTYPE
c     jcf dg - added variable SEGTYPE to record boundary segment IBTYPE
         SEGTYPE(K) = IBTYPE
C.......CHECK THAT IBTYPE PARAMETER HAS BEEN SET PROPERLY
         IF(    (IBTYPE.NE.0).AND.(IBTYPE.NE.10).AND.(IBTYPE.NE.20)
     &        .AND.(IBTYPE.NE.1).AND.(IBTYPE.NE.11).AND.(IBTYPE.NE.21)
     &        .AND.(IBTYPE.NE.2).AND.(IBTYPE.NE.12).AND.(IBTYPE.NE.22)
     &        .AND.(IBTYPE.NE.3).AND.(IBTYPE.NE.13).AND.(IBTYPE.NE.23)
cjj   wm001 - following 2 lines modified 
     &        .AND.(IBTYPE.NE.4).AND.(IBTYPE.NE.24)
     &        .AND.(IBTYPE.NE.5).AND.(IBTYPE.NE.25).AND.(IBTYPE.NE.30)) THEN
         IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9985) K
         WRITE(16,9985) K
 9985    FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &        //,1X,'THE FLOW BOUNDARY TYPE PARAMETER IBTYPE ',
     &        'HAS NOT BEEN CORRECTLY SET FOR ',
     &        /,1X,'FLOW BOUNDARY SEGMENT NO. ',I8,
     &        /,1X,'USER MUST CORRECT UNIT 14 INPUT FILE',
     &        //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
         STOP
      ENDIF
C.......WRITE OUT INFORMATION TO UNIT 16
      IF((IBTYPE.EQ.4).OR.(IBTYPE.EQ.24)) THEN
         WRITE(16,28) K,NVELL(K),K,2*NVELL(K)
 28      FORMAT(///,5X,'TOTAL NUMBER OF PAIRS FOR FLOW BOUNDARY',
     &        ' SEGMENT',2X,I2,2X,'=',2X,I5,/,
     &        5X,'TOTAL NUMBER OF NODES FOR FLOW BOUNDARY',
     &        ' SEGMENT',2X,I2,2X,'=',2X,I5)
      ELSE
         WRITE(16,128) K,NVELL(K)
 128     FORMAT(///,5X,'TOTAL NUMBER OF NODES FOR FLOW BOUNDARY',
     &        ' SEGMENT',2X,I2,2X,'=',2X,I5)
      ENDIF
C.......CONTINUE PROCESSING FLOW BOUNDARY INFORMATION
      IF(IBTYPE.EQ.0) THEN
         WRITE(16,2340)
 2340    FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,
     &        7X,'NO NORMAL FLOW AS AN ESSENTIAL B.C.',/,
     &        7X,'AND FREE TANGENTIAL SLIP',/)
      ENDIF
      IF(IBTYPE.EQ.1) THEN
         WRITE(16,2341)
 2341    FORMAT(5X,'THIS SEGMENT IS AN INTERNAL BOUNDARY WITH:',/,
     &        7X,'NO NORMAL FLOW AS AN ESSENTIAL B.C.',/,
     &        7X,'AND FREE TANGENTIAL SLIP',/)
      ENDIF
      IF(IBTYPE.EQ.2) THEN
         NFLUXF=1
         WRITE(16,2342)
 2342    FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,
     &        7X,'SPECIFIED NORMAL FLOW AS AN ESSENTIAL B.C.',/,
     &        7X,'AND FREE TANGENTIAL SLIP',/)
      ENDIF
      IF(IBTYPE.EQ.3) THEN
         NFLUXB=1
         WRITE(16,2344)
 2344    FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,
     &        7X,'A BARRIER WHICH ALLOWS FREE SURFACE',
     &        ' SUPERCRITICAL OUTFLOW',/,
     &        7X,'FROM THE DOMAIN ONCE THE BARRIER HAS BEEN',
     &        ' OVERTOPPED',/,
     &        7X,'AND FREE TANGENTIAL SLIP',/)
      ENDIF
      IF(IBTYPE.EQ.4) THEN
         NFLUXIB=1
         WRITE(16,2345)
 2345    FORMAT(5X,'THIS SEGMENT IS AN INTERNAL BARRIER BOUNDARY:',/,
     &        7X,'WITH CROSS BARRIER FLOW TREATED AS AN ESSENTIAL ',
     &        ' NORMAL FLOW BOUNDARY CONDITION',/,
     &        7X,'WHICH LEAVES/ENTERS THE DOMAIN ON ONE SIDE OF ',
     &        ' THE BARRIER AND ENTERS/LEAVES THE DOMAIN ON THE ',/,
     &        7X,'CORRESPONDING OPPOSITE SIDE OF THE BARRIER ',/,
     &        7X,'FLOW RATE AND DIRECTION ARE BASED ON BARRIER ',
     &        ' HEIGHT, SURFACE WATER ELEVATION',/,
     &        7X,'ON BOTH SIDES OF THE BARRIER, BARRIER COEFFICIENT',
     &        ' AND THE APPROPRIATE BARRIER FLOW FORMULA',/,
     &        7X,'FREE TANGENTIAL SLIP IS ALLOWED',/)
      ENDIF
      IF(IBTYPE.EQ.5) THEN
         NFLUXIB=1
         NFLUXIBP=1
         WRITE(16,2347)
 2347    FORMAT(5X,'THIS SEGMENT IS AN INTERNAL BARRIER BOUNDARY:',/,
     &        7X,'WITH ADDITIONAL CROSS BARRIER PIPES ',
     &        'LOCATED UNDER THE CROWN ',/,
     &        7X,'CROSS BARRIER FLOW IS TREATED AS AN ESSENTIAL',
     &        ' NORMAL FLOW BOUNDARY CONDITION',/,
     &        7X,'WHICH LEAVES/ENTERS THE DOMAIN ON ONE SIDE OF ',
     &        ' THE BARRIER AND ENTERS/LEAVES THE DOMAIN ON THE ',/,
     &        7X,'CORRESPONDING OPPOSITE SIDE OF THE BARRIER ',/,
     &        7X,'FLOW RATE AND DIRECTION ARE BASED ON BARRIER ',
     &        ' HEIGHT, SURFACE WATER ELEVATION',/,
     &        7X,'ON BOTH SIDES OF THE BARRIER, BARRIER COEFFICIENT',
     &        ' AND THE APPROPRIATE BARRIER FLOW FORMULA',/,
     &        7X,'IN ADDITION CROSS BARRIER PIPE FLOW RATE AND ',
     &        ' DIRECTION ARE BASED ON PIPE CROWN HEIGHT, ',/,
     &        7X,'SURFACE WATER ELEVATION ON BOTH SIDES OF THE ',
     &        'BARRIER, PIPE FRICTION COEFFICIENT, PIPE DIAMETER',/,
     &        7X,' AND THE APPROPRIATE PIPE FLOW FORMULA',/,
     &        7X,'FREE TANGENTIAL SLIP IS ALLOWED',/)
      ENDIF
      IF(IBTYPE.EQ.10) THEN
         WRITE(16,2350)
 2350    FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,
     &        7X,'NO NORMAL FLOW AS AN ESSENTIAL B.C.',/,
     &        7X,'AND NO TANGENTIAL SLIP',/)
      ENDIF
      IF(IBTYPE.EQ.11) THEN
         WRITE(16,2351)
 2351    FORMAT(5X,'THIS SEGMENT IS AN INTERNAL BOUNDARY WITH:',/,
     &        7X,'NO NORMAL FLOW AS AN ESSENTIAL B.C.',/,
     &        7X,'AND NO TANGENTIAL SLIP',/)
      ENDIF
      IF(IBTYPE.EQ.12) THEN
         NFLUXF=1
         WRITE(16,2352)
 2352    FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,
     &        7X,'SPECIFIED NORMAL FLOW AS AN ESSENTIAL B.C.',/,
     &        7X,'AND NO TANGENTIAL SLIP',/)
      ENDIF
      IF(IBTYPE.EQ.13) THEN
         NFLUXB=1
         WRITE(16,2354)
 2354    FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,
     &        7X,'A BARRIER WHICH ALLOWS FREE SURFACE',
     &        ' SUPERCRITICAL OUTFLOW',/,
     &        7X,'FROM THE DOMAIN ONCE THE BARRIER HAS BEEN',
     &        ' OVERTOPPED',/,
     &        7X,'AND NO TANGENTIAL SLIP',/)
      ENDIF
      IF(IBTYPE.EQ.20) THEN
         WRITE(16,2360)
 2360    FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,
     &        7X,'NO NORMAL FLOW AS A NATURAL B.C.',/,
     &        7X,'AND FREE TANGENTIAL SLIP',/)
      ENDIF
      IF(IBTYPE.EQ.21) THEN
         WRITE(16,2361)
 2361    FORMAT(5X,'THIS SEGMENT IS AN INTERNAL BOUNDARY WITH:',/,
     &        7X,'NO NORMAL FLOW AS A NATURAL B.C.',/,
     &        7X,'AND FREE TANGENTIAL SLIP',/)
      ENDIF
      IF(IBTYPE.EQ.22) THEN
         NFLUXF=1
         WRITE(16,2362)
 2362    FORMAT(5X,'THIS SEGMENT IS A EXTERNAL BOUNDARY WITH:',/,
     &        7X,'SPECIFIED NORMAL FLOW AS A NATURAL B.C.',/,
     &        7X,'AND FREE TANGENTIAL SLIP',/)
      ENDIF
      IF(IBTYPE.EQ.23) THEN
         NFLUXB=1
         WRITE(16,2356)
 2356    FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,
     &        7X,'A BARRIER WHICH ALLOWS FREE SURFACE',
     &        ' SUPERCRITICAL OUTFLOW',/,
     &        7X,'FROM THE DOMAIN ONCE THE BARRIER HAS BEEN',
     &        ' OVERTOPPED',/,
     &        7X,' IMPLEMENTED AS A NATURAL BOUNDARY CONDITION'
     &        ,7X,'FREE TANGENTIAL SLIP IS ALSO ALLOWED',/)
      ENDIF
      IF(IBTYPE.EQ.24) THEN
         NFLUXIB=1
         WRITE(16,2357)
 2357    FORMAT(5X,'THIS SEGMENT IS AN INTERNAL BARRIER BOUNDARY:',/,
     &        7X,'WITH CROSS BARRIER FLOW TREATED AS A NATURAL ',
     &        ' NORMAL FLOW BOUNDARY CONDITION',/,
     &        7X,'WHICH LEAVES/ENTERS THE DOMAIN ON ONE SIDE OF ',
     &        ' THE BARRIER AND ENTERS/LEAVES THE DOMAIN ON THE ',/,
     &        7X,'CORRESPONDING OPPOSITE SIDE OF THE BARRIER ',/,
     &        7X,'FLOW RATE AND DIRECTION ARE BASED ON BARRIER ',
     &        ' HEIGHT, SURFACE WATER ELEVATION',/,
     &        7X,'ON BOTH SIDES OF THE BARRIER, BARRIER COEFFICIENT',
     &        ' AND THE APPROPRIATE BARRIER FLOW FORMULA',/,
     &        7X,'FREE TANGENTIAL SLIP IS ALLOWED',/)
      ENDIF
      IF(IBTYPE.EQ.25) THEN
         NFLUXIB=1
         NFLUXIBP=1
         WRITE(16,2359)
 2359    FORMAT(5X,'THIS SEGMENT IS AN INTERNAL BARRIER BOUNDARY:',/,
     &        7X,'WITH ADDITIONAL CROSS BARRIER PIPES ',
     &        'LOCATED UNDER THE CROWN ',/,
     &        7X,'CROSS BARRIER FLOW IS TREATED AS A NATURAL',
     &        ' NORMAL FLOW BOUNDARY CONDITION',/,
     &        7X,'WHICH LEAVES/ENTERS THE DOMAIN ON ONE SIDE OF ',
     &        ' THE BARRIER AND ENTERS/LEAVES THE DOMAIN ON THE ',/,
     &        7X,'CORRESPONDING OPPOSITE SIDE OF THE BARRIER ',/,
     &        7X,'FLOW RATE AND DIRECTION ARE BASED ON BARRIER ',
     &        ' HEIGHT, SURFACE WATER ELEVATION',/,
     &        7X,'ON BOTH SIDES OF THE BARRIER, BARRIER COEFFICIENT',
     &        ' AND THE APPROPRIATE BARRIER FLOW FORMULA',/,
     &        7X,'IN ADDITION CROSS BARRIER PIPE FLOW RATE AND ',
     &        ' DIRECTION ARE BASED ON PIPE CROWN HEIGHT, ',/,
     &        7X,'SURFACE WATER ELEVATION ON BOTH SIDES OF THE ',
     &        'BARRIER, PIPE FRICTION COEFFICIENT, PIPE DIAMETER',/,
     &        7X,' AND THE APPROPRIATE PIPE FLOW FORMULA',/,
     &        7X,'FREE TANGENTIAL SLIP IS ALLOWED',/)
      ENDIF
      IF(IBTYPE.EQ.30) THEN
         NFLUXRBC=1
         WRITE(16,2355)
 2355    FORMAT(5X,'THIS SEGMENT IS AN OUTWARD RADIATING BOUNDARY:',/,
     &        7X,'NORMAL FLUX IS A NATURAL B.C. IN GWCE',/,
     &        7X,'NORMAL AND TANGENTIAL VELOCITY ARE COMPUTED FROM ',
     &        'THE MOMENTUM EQNS.',/)
      ENDIF


C.....INPUT INFORMATION FOR VARIOUS TYPES OF FLOW BOUNDARY SEGMENTS
C.......INPUT THE STANDARD NODE NUMBERS FOR THE Kth FLOW BOUNDARY SEGMENT
      IF((IBTYPE.NE.3).AND.(IBTYPE.NE.13).AND.(IBTYPE.NE.23).AND.
     &     (IBTYPE.NE.4).AND.(IBTYPE.NE.24).AND.
     &     (IBTYPE.NE.5).AND.(IBTYPE.NE.25)) THEN
         DO I=1,NVELL(K)
            READ(14,*) NBVV(K,I)
         END DO
         NPRBI=1
         NPIPE=0
      ENDIF
C.......INPUT THE NODE NUMBERS FOR THE Kth EXTERNAL BARRIER BOUNDARY SEGMENT
C........ALSO INPUT THE ELEVATION OF THE EXTERNAL BARRIER NODES ABOVE
C........THE GEOID AND THE COEFFICIENT OF FREE SURFACE SUPERCRITICAL
C........FLOW ALONG WITH EACH EXTERNAL BARRIER BOUNDARY NODE FROM UNIT 14
      IF((IBTYPE.EQ.3).OR.(IBTYPE.EQ.13).OR.(IBTYPE.EQ.23)) THEN
         DO I=1,NVELL(K)
            READ(14,*) NBVV(K,I),BARLANHTR(I),BARLANCFSPR(I)
         ENDDO
         NPRBI=1
         NPIPE=0
      ENDIF
C.......INPUT THE NODE NUMBERS FOR THE Kth INTERNAL BARRIER BOUNDARY SEGMENT
C........ALSO INPUT CONNECTION NODE NUMBER AND ELEVATION OF THE INTERNAL BARRIER
C........NODES ABOVE THE GEOID AND THE COEFFICIENTS OF FREE SURFACE SUPERCRITICAL
C........AND SUBCRITICAL FLOW ALONG WITH EACH INTERNAL BARRIER BOUNDARY NODE FROM
C........UNIT 14
      IF((IBTYPE.EQ.4).OR.(IBTYPE.EQ.24)) THEN
         DO I=1,NVELL(K)
            READ(14,*) NBVV(K,I), IBCONNR(I), BARINHTR(I), BARINCFSBR(I)
     &           , BARINCFSPR(I)
         ENDDO
         NPRBI=2
         NPIPE=0
      ENDIF
C.......INPUT THE NODE NUMBERS FOR THE Kth INTERNAL BARRIER BOUNDARY SEGMENT WITH
C........CROSS BARRIER PIPES; ALSO INPUT CONNECTION NODE NUMBER AND ELEVATION OF THE 
C.......INTERNAL BARRIER NODES ABOVE THE GEOID AND THE COEFFICIENTS OF FREE SURFACE 
C.......SUPERCRITICAL AND SUBCRITICAL FLOW ALONG WITH EACH INTERNAL BARRIER BOUNDARY
C.......NODE FROM UNIT 14; IN ADDITION INPUT THE CROSS BARRIER PIPE HEIGHT, CROSS
C.......BARRIER PIPE COEFFICIENT AND CROSS BARRIER PIPE DIAMETER
      IF((IBTYPE.EQ.5).OR.(IBTYPE.EQ.25)) THEN
         DO I=1,NVELL(K)
            READ(14,*) NBVV(K,I),IBCONNR(I),BARINHTR(I),BARINCFSBR(I),
     &           BARINCFSPR(I),PIPEHTR(I),PIPECOEFR(I),
     &           PIPEDIAMR(I)
         END DO
         NPRBI=2
         NPIPE=1
      ENDIF

C.....PROCESS INFORMATION FOR VARIOUS TYPES OF FLOW BOUNDARY SEGMENTS

      DO IPRBI=1,NPRBI

C.........LOAD PAIRED NODES INTO PRIMARY PROCESSING VECTORS AND RESET
C..........CONNECTING NODES FOR BACK FACE
C..........THUS BACK/CONNECTING NODES ARE BEING LOADED AS PRIMARY NODES
C..........AND FRONT NODES ARE RELOADED AS CONNECTING NODES
C..........NOTE THAT THE CLOCKWISE ORIENTATION OF ISLAND TYPE BOUNDARIES
C..........IS BEING MAINTAINED WHEN BACK NODES ARE RELOADED AS PRIMARY NODES
C..........ADDITIONAL INTERNAL BARRIER BOUNDARY INFORMATION IS ALSO RESET
         IF(IPRBI.EQ.2) THEN
            DO I=1,NVELL(K)
               NTRAN1(I)=NBVV(K,I)
               NTRAN2(I)=IBCONNR(I)
               BTRAN3(I)=BARINHTR(I)
               BTRAN4(I)=BARINCFSBR(I)
               BTRAN5(I)=BARINCFSPR(I)
               IF(NPIPE.EQ.1) THEN
                  BTRAN6(I)=PIPEHTR(I)
                  BTRAN7(I)=PIPECOEFR(I)
                  BTRAN8(I)=PIPEDIAMR(I)
               ENDIF
            END DO
            DO I=1,NVELL(K)
               NBVV(K,I)=NTRAN2(NVELL(K)+1-I)
               IBCONNR(I)=NTRAN1(NVELL(K)+1-I)
               BARINHTR(I)=BTRAN3(NVELL(K)+1-I)
               BARINCFSBR(I)=BTRAN4(NVELL(K)+1-I)
               BARINCFSPR(I)=BTRAN5(NVELL(K)+1-I)
               IF(NPIPE.EQ.1) THEN
                  PIPEHTR(I)=BTRAN6(NVELL(K)+1-I)
                  PIPECOEFR(I)=BTRAN7(NVELL(K)+1-I)
                  PIPEDIAMR(I)=BTRAN8(NVELL(K)+1-I)
               ENDIF
            ENDDO
         ENDIF

C.........WRITE OUT ADDITIONAL HEADER FOR INTERNAL BARRIER BOUNDARIES

         IF((IBTYPE.EQ.4).OR.(IBTYPE.EQ.24)) THEN
            IF(IPRBI.EQ.1) THEN
               WRITE(16,1842)
 1842          FORMAT(/,5X,'FRONT FACE OF INTERNAL BARRIER BOUNDARY',/)
            ELSE
               WRITE(16,1843)
 1843          FORMAT(/,5X,'BACK FACE OF INTERNAL BARRIER BOUNDARY',/)
            ENDIF
         ENDIF

C.........WRITE OUT ADDITIONAL HEADER FOR INTERNAL BARRIER BOUNDARIES 
C..........WITH CROSS BARRIER PIPES

         IF((IBTYPE.EQ.5).OR.(IBTYPE.EQ.25)) THEN
            IF(IPRBI.EQ.1) THEN
               WRITE(16,1844)
 1844          FORMAT(/,5X,'FRONT FACE OF INTERNAL BARRIER BOUNDARY',
     &              ' WITH CROSS BARRIER PIPES',/)
            ELSE
               WRITE(16,1845)
 1845          FORMAT(/,5X,'BACK FACE OF INTERNAL BARRIER BOUNDARY',
     &              ' WITH CROSS BARRIER PIPES',/)
            ENDIF
         ENDIF

C.........WRITE OUT GENERAL HEADER FOR BOUNDARY INFORMATION

         WRITE(16,1841)
 1841    FORMAT('    JGW    JME    ME2GW   NODE #  BNDRY CODE   INNER',
     &        ' ANGLE',7X,'COS',13X,'SIN',9X,'0.667*BNDRY LEN',/)

C.........COMPLETE THE BOUNDARY ARRAY FOR THE Kth FLOW BOUNDARY SEGMENT

         NBVV(K,0)=NBVV(K,1)    !UNCLOSED EXTERNAL
         IF((IBTYPE.EQ.1).OR.(IBTYPE.EQ.11).OR.(IBTYPE.EQ.21)) THEN
            IF(NBVV(K,NVELL(K)).NE.NBVV(K,1)) THEN !CLOSE AN UNCLOSED INTERNAL
               NVELL(K)=NVELL(K)+1
               NBVV(K,NVELL(K))=NBVV(K,1)
            ENDIF
         ENDIF
         IF(NBVV(K,NVELL(K)).EQ.NBVV(K,1)) THEN !CLOSED EXTERNAL OR INTERNAL
            NBVV(K,0)=NBVV(K,NVELL(K)-1)
         ENDIF
         NBVV(K,NVELL(K)+1)=NBVV(K,NVELL(K))

C.........PUT BOUNDARY INFORMATION INTO 2 TYPES OF ARRAYS, ONE FOR THE GWCE B.C.
C..........AND ONE FOR THE MOMENTUM EQUATION B.C.
C..........THE GWCE ARRAYS INCLUDE EVERY NODE IN THE UNIT 14 FILE, I.E., NODES
C..........ARE REPEATED WHERE SPECIFIED NORMAL FLOW AND NO NORMAL FLOW BOUNDARIES
C..........MEET AND AT THE BEGINNING AND END OF CLOSED EXTERNAL BOUNDARIES AND
C.........ISLANDS.
C..........THE MOMENTUM EQUATION ARRAYS ARE KEYED TO THE GWCE ARRAYS VIA THE
C..........ARRAY ME2GW WHICH INDICATES THE LOCATION IN THE GWCE ARRAYS THAT
C..........THE APPROPRIATE M.E. VALUE LIES.
C..........THE M.E. ARRAYS DO NOT REPEAT NODES THAT ARE DUPLICATED IN THE
C..........UNIT 14 FILE, I.E., WHEN SPECIFIED NORMAL FLOW AND NO NORMAL FLOW
C..........BOUNDARIES MEET, THE SPECIFIED NORMAL FLOW BOUNDARY CONDITION TAKES
C..........PRECEDENT.  ALSO THE BEGINNING AND ENDING NODES OF CLOSED EXTERNAL
C..........AND ISLAND BOUNDARIES ARE NOT REPEATED.

         DO I=1,NVELL(K)

C.........SET UP THE GWCE BOUNDARY ARRAYS WHICH CONSIST OF
C..........BOUNDARY NODE NUMBERS
C..........BOUNDARY CODES
C..........0.66667*LENGTH OF EACH BOUNDARY SEGMENT.  NOTE, THE LENGTH OF THE LAST
C..........BOUNDARY SEGMENT ON EACH BOUNDARY SHOULD BE ZERO

            JGW=JGW+1
            IF((IBTYPE.EQ.0).OR.(IBTYPE.EQ.10).OR.(IBTYPE.EQ.20).OR.
     &           (IBTYPE.EQ.2).OR.(IBTYPE.EQ.12).OR.(IBTYPE.EQ.22).OR.
     &           (IBTYPE.EQ.3).OR.(IBTYPE.EQ.13).OR.(IBTYPE.EQ.23).OR.
     &           (IBTYPE.EQ.30)) THEN
               NVELEXT=NVELEXT+1
            ENDIF
            NBV(JGW)=NBVV(K,I)
            NBVI=NBVV(K,I)
            NBVJ=NBVV(K,I+1)
            DELX=X(NBVJ)-X(NBVI)
            DELY=Y(NBVJ)-Y(NBVI)
            BNDLEN2O3(JGW)=2.D0*(SQRT(DELX*DELX+DELY*DELY))/3.D0

C...........COMPUTE THE INCLUDED ANGLE AND TEST TO DETERMINE WHETHER TO ZERO
C............TANGENTIAL VELOCITIES
C...........NOTE:.IMPLEMENTATION FOR ICS=2 REQUIRES COMPUTING ALL COORDINATES IN
C............A LOCALIZED SYSTEM (I.E. THE TRANSFORMATION IS CENTERED AT X0,Y0)

            IF(ICS.EQ.1) THEN
               XL0=X(NBVV(K,I))
               XL1=X(NBVV(K,I-1))
               XL2=X(NBVV(K,I+1))
               YL0=Y(NBVV(K,I))
               YL1=Y(NBVV(K,I-1))
               YL2=Y(NBVV(K,I+1))
            ELSE
               CALL CPP(XL0,YL0,SLAM(NBVV(K,I)),SFEA(NBVV(K,I)),
     &              SLAM(NBVV(K,I)),SFEA(NBVV(K,I)))
               CALL CPP(XL1,YL1,SLAM(NBVV(K,I-1)),SFEA(NBVV(K,I-1)),
     &              SLAM(NBVV(K,I)),SFEA(NBVV(K,I)))
               CALL CPP(XL2,YL2,SLAM(NBVV(K,I+1)),SFEA(NBVV(K,I+1)),
     &              SLAM(NBVV(K,I)),SFEA(NBVV(K,I)))
            ENDIF

C...........NOTE: INTERIOR ANGLE AT ENDS OF BOUNDARIES MUST BE EQUAL, EITHER:
C............A FICTICIOUSLY LARGE VALUE IF THE BOUNDARY IS NOT CLOSED OR
C............A TRUE VALUE IF THE BOUNDARY IS CLOSED

            THETA=0.
            IF((I.EQ.1).AND.(NBVV(K,I).EQ.NBVV(K,I-1))) THEN
               THETA1=-9999999.d0
               THETA=THETA1
               COSTHETA1=COSTSET
               COSTHETA=COSTHETA1
               CROSS1=0.d0
               CROSS=CROSS1
            ENDIF
            IF(I.EQ.NVELL(K)) THEN
               THETA=THETA1
               COSTHETA=COSTHETA1
               CROSS=CROSS1
            ENDIF
            IF(THETA.EQ.0.) THEN
               VL1X=XL1-XL0
               VL1Y=YL1-YL0
               VL2X=XL2-XL0
               VL2Y=YL2-YL0
               DOTVEC=VL1X*VL2X+VL1Y*VL2Y
               VECNORM=(SQRT(VL1X**2+VL1Y**2))*(SQRT(VL2X**2+VL2Y**2))
               COSTHETA=DOTVEC/VECNORM
               IF(COSTHETA.GT.1.0d0) COSTHETA=1.0d0
               IF(COSTHETA.LT.-1.0d0) COSTHETA=-1.0d0
               THETA=RAD2DEG*ACOS(COSTHETA)
               CROSS=-VL1X*VL2Y+VL2X*VL1Y
               IF(CROSS.LT.0) THETA=360.d0-THETA
               IF(I.EQ.1) THEN
                  THETA1=THETA
                  COSTHETA1=COSTHETA
                  CROSS1=CROSS
               ENDIF
            ENDIF

C...........CHECK WHETHER ANGLE IS LESS THAN MINIMUM ANGLE, IF SO CHANGE THE
C............BOUNDARY CODE TO ZERO TANGENTIAL VELOCITIES

            LBCODEI(JGW)=IBTYPE
            IF((COSTHETA.GT.COSTSET).AND.(CROSS.GT.0.0)) THEN
               IF(IBTYPE.EQ.0) LBCODEI(JGW)=10
               IF(IBTYPE.EQ.1) LBCODEI(JGW)=11
               IF(IBTYPE.EQ.2) LBCODEI(JGW)=12
               IF(IBTYPE.EQ.3) LBCODEI(JGW)=13
               IF((IBTYPE.GE.0).AND.(IBTYPE.LE.3)) THEN
                  WRITE(16,1856) NBVV(K,I),THETA
 1856             FORMAT(2X,I7,4X,'THE INNER ANGLE = ',F8.2,1X,
     &                 'TANGENTIAL SLIP WILL BE ZEROED')
               ENDIF
            ENDIF

C...........COMPUTE COS AND SIN OF OUTWARD NORMAL REGARDLESS OF BOUNDARY TYPE

            X1=X(NBVV(K,I-1))
            X2=X(NBVV(K,I+1))
            Y1=Y(NBVV(K,I-1))
            Y2=Y(NBVV(K,I+1))
            XL=SQRT((X1-X2)**2+(Y1-Y2)**2)
            CSII(JGW)=SFAC(NBVV(K,I))*(Y2-Y1)/XL
            SIII(JGW)=(X1-X2)/XL

C...........SET UP THE MOMENTUM EQUATION BOUNDARY ARRAY WHICH CONSISTS OF
C............A KEY TO THE GWCE BOUNDARY CONDITION ARRAY

            IF(I.EQ.1) THEN     !DEAL WITH FIRST NODE IN L.B. SEG
               IF(JGW.EQ.1) THEN !VERY FIRST L.B. SEG
                  JME=JME+1     !M.E. USES IT
                  ME2GW(JME)=JGW
               ENDIF
               IF(JGW.NE.1) THEN
                  IF(NBV(JGW).NE.NBV(JGW-1)) THEN !L.B. SEGS DON'T OVERLAP
                     JME=JME+1  !M.E. USES IT
                     ME2GW(JME)=JGW
                  ENDIF
                  IF(NBV(JGW).EQ.NBV(JGW-1)) THEN !L.B. SEGS OVERLAP
                     IF((LBCODEI(JGW).EQ.2) .OR. ! M.E. USES IT ONLY
     &                    (LBCODEI(JGW).EQ.12).OR. ! IF IT IS
     &                    (LBCODEI(JGW).EQ.22).OR. ! SPECIFIED FLOW,
     &                    (LBCODEI(JGW).EQ.3) .OR. ! AN OVERFLOW BARRIER
     &                    (LBCODEI(JGW).EQ.13).OR. ! OR A RADIATION
     &                    (LBCODEI(JGW).EQ.23).OR. ! BOUNDARY
     &                    (LBCODEI(JGW).EQ.30)) ME2GW(JME)=JGW
                  ENDIF                                                                
               ENDIF
            ENDIF
            IF((I.GT.1).AND.(I.LT.NVELL(K))) THEN !IF NOT FIRST OR
               JME=JME+1        !LAST NODE
               ME2GW(JME)=JGW   !M.E. USES IT
            ENDIF
            IF(I.EQ.NVELL(K)) THEN !DEAL WITH LAST NODE ON BOUNDARY
               IF((NBV(JGW).NE.NBVV(K,1)).AND. !IF UNCLOSED BOUNDARY
     &              (NBV(JGW).NE.NBV(1))) THEN !M.E. USES IT
                  JME=JME+1
                  ME2GW(JME)=JGW
               ENDIF
               IF(NBVV(K,I).EQ.NBV(1)) THEN !IF OVERLAPS WITH VERY FIRST
                  IF((LBCODEI(JGW).EQ.2) .OR. ! L.B. NODE
     &                 (LBCODEI(JGW).EQ.12).OR. ! M.E. USES IT ONLY IF IT IS
     &                 (LBCODEI(JGW).EQ.22).OR. ! SPECIFIED FLOW,
     &                 (LBCODEI(JGW).EQ.3) .OR. ! AN OVERFLOW BARRIER OR
     &                 (LBCODEI(JGW).EQ.13).OR. ! A RADIATION
     &                 (LBCODEI(JGW).EQ.23).OR. ! BOUNDARY
     &                 (LBCODEI(JGW).EQ.30)) ME2GW(1)=JGW
               ENDIF
            ENDIF

C...........LOAD EXTERNAL BARRIER BOUNDARY INFORMATION INTO THE CORRECT VECTORS
            IF((IBTYPE.EQ.3).OR.(IBTYPE.EQ.13).OR.(IBTYPE.EQ.23)) THEN
               BARLANHT(JGW)=BARLANHTR(I)
               BARLANCFSP(JGW)=BARLANCFSPR(I)
            ENDIF

C...........LOAD INTERNAL BARRIER BOUNDARY INFORMATION INTO THE CORRECT VECTORS
            IF((IBTYPE.EQ.4).OR.(IBTYPE.EQ.24)) THEN
               IBCONN(JGW)=IBCONNR(I)
               BARINHT(JGW)=BARINHTR(I)
               BARINCFSB(JGW)=BARINCFSBR(I)
               BARINCFSP(JGW)=BARINCFSPR(I)
            ENDIF

C...........LOAD INTERNAL BARRIER WITH PIPES BOUNDARY INFORMATION INTO 
C............THE CORRECT VECTORS
            IF((IBTYPE.EQ.5).OR.(IBTYPE.EQ.25)) THEN
               IBCONN(JGW)=IBCONNR(I)
               BARINHT(JGW)=BARINHTR(I)
               BARINCFSB(JGW)=BARINCFSBR(I)
               BARINCFSP(JGW)=BARINCFSPR(I)
               PIPEHT(JGW)=PIPEHTR(I)
               PIPECOEF(JGW)=PIPECOEFR(I)
               PIPEDIAM(JGW)=PIPEDIAMR(I)
            ENDIF
            
C...........WRITE OUT BOUNDARY CONDITION ARRAY INFORMATION

            WRITE(16,1857) JGW,JME,ME2GW(JME),NBV(JGW),LBCODEI(JGW),
     &           THETA,CSII(JGW),SIII(JGW),BNDLEN2O3(JGW)
 1857       FORMAT(1X,I6,1X,I6,1X,I6,3X,I6,3X,I4,9X,F8.2,2X,E16.8,1X,
     &           E16.8,2X,E16.8)

C...........CHECK EXTERNAL BARRIER HEIGHTS AGAINST DEPTHS
            IF((IBTYPE.EQ.3).OR.(IBTYPE.EQ.13).OR.(IBTYPE.EQ.23)) THEN
               IF(BARLANHT(JGW).LT.-DP(NBV(JGW))) THEN
                  IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,8367) 
     &                 JGW,NBV(JGW),BARLANHT(JGW),DP(NBV(JGW))
                  WRITE(16,8367) JGW,NBV(JGW),BARLANHT(JGW),DP(NBV(JGW))
 8367             FORMAT(////,1X,'!!!!!!!!!!  FATAL INPUT ERROR   !!!'
     &                 ,'!!!!!!!!!',//,
     &                 1X,'AT BOUNDARY NODE NO.',I6,' (GLOBAL NODE NO.',
     &                 I6, ' AND OF EXTERNAL BARRIER TYPE) ',/,
     &                 2X,'THE EXTERNAL BARRIER HEIGHT = ',E12.5,
     &                 2X,'IS EXCEEDED BY THE DEPTH SPECIFIED AT ',/,2X
     &                 ,'THE ASSOCIATED GLOBAL NODE = ',E12.5,/,2X,
     &                 'USER MUST SPECIFY CONSISTENT BARRIER HEIGHTS',
     &                 ' AND DEPTHS')
                  IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9973)
                  WRITE(16,9973)
                  STOP
               ENDIF
            ENDIF
C...........CHECK INTERNAL BARRIER HEIGHTS AGAINST DEPTHS
            IF((IBTYPE.EQ.4).OR.(IBTYPE.EQ.24)) THEN
               IF(BARINHT(JGW).LT.-DP(NBV(JGW))) THEN
                  IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,8368) 
     &                 JGW,NBV(JGW),BARINHT(JGW),DP(NBV(JGW))
                  WRITE(16,8368) JGW,NBV(JGW),BARINHT(JGW),DP(NBV(JGW))
 8368             FORMAT(////,1X,'!!!!!!!!!!  FATAL INPUT ERROR   !!!'
     &                 ,'!!!!!!!!!',//,
     &                 1X,'AT BOUNDARY NODE NO.',I6,' (GLOBAL NODE NO. ',
     &                 I6,' AND OF INTERNAL BARRIER TYPE) ',/,
     &                 2X,'THE INTERNAL BARRIER HEIGHT = ',E12.5,
     &                 2X,'IS EXCEEDED BY THE DEPTH SPECIFIED AT ',/,2X
     &                 ,'THE ASSOCIATED GLOBAL NODE = ',E12.5,/,2X,
     &                 'USER MUST SPECIFY CONSISTENT BARRIER HEIGHTS',
     &                 ' AND DEPTHS')
                  IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9973)
                  WRITE(16,9973)
                  STOP
               ENDIF
            ENDIF
            
C...........CHECK INTERNAL BARRIER WITH PIPES BARRIER HEIGHTS AGAINST DEPTHS
            IF((IBTYPE.EQ.5).OR.(IBTYPE.EQ.25)) THEN
               IF(BARINHT(JGW).LT.-DP(NBV(JGW))) THEN
                  IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,8370) 
     &                 JGW,NBV(JGW),BARINHT(JGW),DP(NBV(JGW))
                  WRITE(16,8370) JGW,NBV(JGW),BARINHT(JGW),DP(NBV(JGW))
 8370             FORMAT(////,1X,'!!!!!!!!!!  FATAL INPUT ERROR   !!!'
     &                 ,'!!!!!!!!!',//,
     &                 1X,'AT BOUNDARY NODE NO.',I6,' (GLOBAL NODE NO. ',
     &                 I6,' AND OF INTERNAL BARRIER TYPE) ',/,
     &                 2X,'THE INTERNAL BARRIER HEIGHT = ',E12.5,
     &                 2X,'IS EXCEEDED BY THE DEPTH SPECIFIED AT ',/,2X
     &                 ,'THE ASSOCIATED GLOBAL NODE = ',E12.5,/,2X,
     &                 'USER MUST SPECIFY CONSISTENT BARRIER HEIGHTS',
     &                 ' AND DEPTHS')
                  IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9973)
                  WRITE(16,9973)
                  STOP
               ENDIF
            ENDIF
C...........CHECK INTERNAL BARRIER WITH PIPES PIPE HEIGHTS AGAINST DEPTHS
            IF((IBTYPE.EQ.5).OR.(IBTYPE.EQ.25)) THEN
               IF(PIPEHT(JGW).LT.-DP(NBV(JGW))) THEN
                  IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,8372) 
     &                 JGW,NBV(JGW),BARINHT(JGW),DP(NBV(JGW))
                  WRITE(16,8372) JGW,NBV(JGW),BARINHT(JGW),DP(NBV(JGW))
 8372             FORMAT(////,1X,'!!!!!!!!!!  FATAL INPUT ERROR   !!!'
     &                 ,'!!!!!!!!!',//,
     &                 1X,'AT BOUNDARY NODE NO.',I6,' (GLOBAL NODE NO. ',
     &                 I6,' AND OF INTERNAL BARRIER TYPE) ',/,
     &                 2X,'THE BARRIER PIPE HEIGHT = ',E12.5,
     &                 2X,'IS EXCEEDED BY THE DEPTH SPECIFIED AT ',/,2X
     &                 ,'THE ASSOCIATED GLOBAL NODE = ',E12.5,/,2X,
     &                 'USER MUST SPECIFY CONSISTENT PIPE HEIGHTS',
     &                 ' AND DEPTHS')
                  IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9973)
                  WRITE(16,9973)
                  STOP
               ENDIF
            ENDIF

C...........CHECK FOR OVERLAPPING OF AN INTERNAL BARRIER BOUNDARY WITH
C............ANY EXTERNAL BARRIER BOUNDARY. IF THIS DOES OCCUR, TAKE
C............APPROPRIATE ACTION
            IF((IBTYPE.EQ.4).OR.(IBTYPE.EQ.24).OR.(IBTYPE.EQ.5)
     &           .OR.(IBTYPE.EQ.25)) THEN
               DO ICK=1,NVELEXT
C...............CHECK IF OVERLAP EXISTS
                  IF(NBV(ICK).EQ.NBV(JGW)) THEN
C.................CHECK FOR ILLEGAL OVERLAPS
                     IF((LBCODEI(ICK).EQ.2).OR.(LBCODEI(ICK).EQ.3).OR.
     &                    (LBCODEI(ICK).EQ.12).OR.(LBCODEI(ICK).EQ.13)) THEN 
                        IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,8567) 
     &                       JGW,NBV(JGW),ICK,NBV(ICK)
                        WRITE(16,8567) JGW,NBV(JGW),ICK,NBV(ICK)
 8567                   FORMAT(////,1X,'!!!!!!!!!!  FATAL INPUT ERROR   !!!'
     &                       ,'!!!!!!!!!',//,
     &                       1X,'BOUNDARY NODE NO. ',I6,' (GLOBAL NODE NO. ',
     &                       I9, 'AND OF INTERNAL BARRIER TYPE) ',/,
     &                       2X,'OVERLAPS BOUNDARY NODE NO.',I6,' (GLOBAL NODE'
     &                       ,' NO.',I6,' )',/,
     &                       2X,'THIS IS AN ILLEGAL TYPE OVERLAP !! - INTERNAL '
     &                       ,'BARRIER BOUNDARIES CAN ONLY OVERLAP WITH ',
     &                       'NO NORMAL FLOW EXTERNAL BOUNDARIES',/
     &                       2X,'(I.E. IBTYPE=0,10,20)')
                        IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9973)
                        WRITE(16,9973)
                        STOP
                     ENDIF
C.................CHECK FOR OVERLAPS WHICH REQUIRE ADJUSTMENTS OF BOUNDARY
C..................CODE ON THE EXTERNAL BOUNDARY
                     IF(((IBTYPE.EQ.4).AND.(LBCODEI(ICK).EQ.0))
     &                    .OR.((IBTYPE.EQ.5).AND.(LBCODEI(ICK).EQ.0))) THEN
                        WRITE(16,8568) JGW,ICK,ICK
 8568                   FORMAT(1X,'DUE TO LEGAL OVERLAPPING OF ',
     &                       'BOUNDARY NODE',I7,' (WHICH IS AN ESSENTIAL INTER'
     &                       ,'NAL BARRIER BOUNDARY NODE)', /,2X,
     &                       'AND BOUNDARY NODE',I7,' (WHICH IS AN ESSENTIAL ',
     &                       'EXTERNAL NO NORMAL FLOW WITH SLIP BOUNDARY',
     &                       ' NODE),',/,2X,
     &                       'THE BOUNDARY TYPE FOR BOUNDARY NODE ',I7,
     &                       ' IS BEING RESET TO IBTYPE=20',/,2X,
     &                       '(NATURAL NO NORMAL FLOW WITH SLIP BOUNDARY) ')
                        LBCODEI(ICK)=20
                     ENDIF
                     IF(((IBTYPE.EQ.4).AND.(LBCODEI(ICK).EQ.10)) 
     &                    .OR.((IBTYPE.EQ.5).AND.(LBCODEI(ICK).EQ.10))) THEN
                        WRITE(16,8569) JGW,ICK,ICK
 8569                   FORMAT(1X,'DUE TO LEGAL OVERLAPPING OF ',
     &                       'BOUNDARY NODE ',I7,' (WHICH IS AN ESSENTIAL INTER'
     &                       ,'NAL BARRIER BOUNDARY NODE)', /,2X,
     &                       'AND BOUNDARY NODE',I7,' (WHICH IS AN ESSENTIAL ',
     &                       'EXTERNAL NO NORMAL FLOW WITH NO SLIP BOUNDARY',
     &                       ' NODE),',/,2X,
     &                       'THE BOUNDARY TYPE FOR BOUNDARY NODE ',I7,
     &                       ' IS BEING RESET TO IBTYPE=20',/,2X,
     &                       '(NATURAL NO NORMAL FLOW WITH SLIP BOUNDARY) ')
                        LBCODEI(ICK)=20
                     ENDIF
                     IF(((IBTYPE.EQ.24).AND.(LBCODEI(ICK).EQ.10))  
     &                    .OR.((IBTYPE.EQ.25).AND.(LBCODEI(ICK).EQ.10))) THEN
                        WRITE(16,8570) JGW,ICK,ICK
 8570                   FORMAT(1X,'DUE TO LEGAL OVERLAPPING OF ',
     &                       'BOUNDARY NODE',I7,' (WHICH IS A NATURAL INTERNAL'
     &                       ,' BARRIER BOUNDARY NODE)', /,2X,
     &                       'AND BOUNDARY NODE',I7,' (WHICH IS AN ESSENTIAL ',
     &                       'EXTERNAL NO NORMAL FLOW WITH NO SLIP BOUNDARY',
     &                       ' NODE),',/,2X,
     &                       'THE BOUNDARY TYPE FOR BOUNDARY NODE',I7,
     &                       ' IS BEING RESET TO IBTYPE=0',/,2X,
     &                       '(ESSENTIAL NO NORMAL FLOW WITH SLIP BOUNDARY) ')
                        LBCODEI(ICK)=0
                     ENDIF
                  ENDIF
               END DO
            ENDIF

         ENDDO
      ENDDO
      
C.......Put barrier data into arrays more amiable to DG data structure

      IF ((IBTYPE.EQ.3).OR.(IBTYPE.EQ.13).OR.(IBTYPE.EQ.23)) THEN
         DO I = 1,NVELL(K)-1
            NEBSEG = NEBSEG + 1
            EBHT(NEBSEG)   = 0.5D0*( BARLANHTR(I)   + BARLANHTR(I+1)   )
            EBCFSP(NEBSEG) = 0.5D0*( BARLANCFSPR(I) + BARLANCFSPR(I+1) )
         ENDDO
      ENDIF
      
      IF ((IBTYPE.EQ.4).OR.(IBTYPE.EQ.24)) THEN
         DO I = 1,NVELL(K)-1
            NIBSEG = NIBSEG + 1
            IBHT(NIBSEG)   = 0.5D0*( BARINHTR(I)   + BARINHTR(I+1)   )
            IBCFSB(NIBSEG) = 0.5D0*( BARINCFSBR(I) + BARINCFSBR(I+1) )
            IBCFSP(NIBSEG) = 0.5D0*( BARINCFSPR(I) + BARINCFSPR(I+1) )
            BACKNODES(1,NIBSEG) = IBCONNR(I)
            BACKNODES(2,NIBSEG) = IBCONNR(I+1)
         ENDDO
      ENDIF
      ENDDO

C.....ONCE ALL FLOW BOUNDARY NODES HAVE BEEN PROCESSED, CHECK TO MAKE SURE
C.....THAT JGW LE MNVEL.  NOTE, JME MUST BE < JGW.

      IF(MNVEL.LT.JGW) THEN
         IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9947)
         WRITE(16,9947)
 9947    FORMAT(////,1X,'!!!!!!!!!!  FATAL INPUT ERROR   !!!!!!!!!!!!',
     &        //,1X,'THE DIMENSION PARAMETER MNVEL IS LESS THAN ',
     &        'THE TOTAL NUMBER OF FLOW BOUNDARY NODES',
     &        /,1X,'FROM ALL THE SPECIFIED FLOW SEGMENTS COMBINED',/)
         IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9973)
         WRITE(16,9973)
         STOP
      ENDIF

      NVEL=JGW
      NVELME=JME

C     
      DO IK=1,MNP               ! FINISH DET. MAX # NEIGHBORS
         IF(NNEIGH(IK).GT.MNEI) MNEI=NNEIGH(IK)
      ENDDO
      MNEI = MNEI+1

C     sb   Copied from read_input.F of v45.01
C     estimate the maximum array space needed for the neighbor
C     table by
C     increasing this number by 2, to provide array space for the
C     node itself
C     and in case the maximum number of nodes occurs at a boundary
C     node

      MNei = MNei+3


C.....TRANSFER FLOW BOUNDARY INFORMATION INTO NODAL ARRAYS

      DO I=1,NP
         LBCODE(I)=-1
         CSI(I)=0.
         SII(I)=1.
      END DO

      DO I=1,NVELME
         J=ME2GW(I)
         LBCODE(NBV(J))=LBCODEI(J)
         CSI(NBV(J))=CSII(J)
         SII(NBV(J))=SIII(J)
      END DO

C...  IF ANY NON ZERO NORMAL FLOW BOUNDARIES WERE SPECIFIED, (NFLUXF=1)
C.....READ FORCING INFORMATION FROM UNIT 15 FILE

      NFFR = 0
      IF(NFLUXF.EQ.1) THEN

C.....INPUT FROM THE NUMBER OF FREQUENCIES PRESENT IN NORMAL FLOW FORCING
C......DATA.  IF THIS = 0, NORMAL FLOW DATA IS READ IN FROM THE FORT.20 FILE.

         READ(15,*) NFFR
         MNFFR = NFFR
         IF (NFFR.EQ.0) MNFFR = 1

C.....Allocate space for periodic normal flow boundary conditions
         call alloc_main6()
C     
         DO I=1,NVELME
            J=ME2GW(I)
            QNAM(1,J)=0.
            QNPH(1,J)=0.
         ENDDO

C.....READ IN AND WRITE OUT INFO ON SPECIFIED NORMAL FLOW BOUNDARIES

         WRITE(16,1112)
         WRITE(16,2200)
 2200    FORMAT(//,1X,'NORMAL FLOW BOUNDARY FORCING INFORMATION ',//)
         IF(NFFR.EQ.0) THEN
            WRITE(16,2201)
 2201       FORMAT(/,5X,'NORMAL FLOW VALUES WILL BE READ FROM UNIT 20 ',
     &      /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE FLOW DATA ',
     &      /,9X,'WITH THE MODEL TIME STEP.')
         ENDIF
         IF(NFFR.NE.0) THEN
            WRITE(16,2202) NFFR
 2202       FORMAT(/,5X,'NUMBER OF PERIODIC NORMAL FLOW CONSTITUENTS =',
     &           I5)
            WRITE(16,2203)
 2203       FORMAT(/,7X,'CONSTITUENT #',4X,'FREQUENCY',4X,'NODAL FACTOR',
     &           3X,'EQU.ARG (DEG)',2X,'CONSTITUENT',/)
            DO I=1,NFFR
               READ(15,'(A5)') FBOUNTAG(I)
               READ(15,*) FAMIG(I),FFF(I),FFACE(I)
               WRITE(16,2204) I,FAMIG(I),FFF(I),FFACE(I),FBOUNTAG(I)
 2204          FORMAT(12X,I2,6X,F16.12,2X,F10.7,2X,F10.3,10X,A5)
               FFACE(I)=FFACE(I)*DEG2RAD
               IF(FAMIG(I).EQ.0.) THEN
                  FPER(I)=0.
               ELSE
                  FPER(I)=2.D0*PI/FAMIG(I)
               ENDIF
            END DO

C.......INPUT PERIODIC NORMAL FLOW FORCING CONDITIONS ON DESIGNATED FLOW BOUNDARIES
C........FOR EACH OF THE FORCING FREQUENCIES FROM UNIT 15 AND OUTPUT TO UNIT 16

            K = 1
            II = 1
            DO I=1,NFFR
               WRITE(16,2206) I,FBOUNTAG(I)
 2206          FORMAT(////,5X,'PERIODIC NORMAL FLOW CONSTITUENT ',
     &              'NUMBER',I4,1X,'DESIGNATED : ',A5)
               READ(15,'(A10)') ALPHA
               WRITE(16,31) ALPHA
               WRITE(16,30)
               DO J=1,NVEL
                  IF((LBCODEI(J).EQ.2).OR.(LBCODEI(J).EQ.12)
     &                 .OR.(LBCODEI(J).EQ.22)) THEN
                     
C.....Modified arrangement of QNAM and QNPH for DG

                     IF (DGSWE.EQ.1) THEN
                        READ(15,*) QNAM(I,II),QNPH(I,II)
                        WRITE(16,2205) NBV(J),QNAM(I,II),QNPH(I,II)
                        QNPH(I,II)=QNPH(I,II)*DEG2RAD
                        II = II + 1
                     ELSE
                        READ(15,*) QNAM(I,J),QNPH(I,J)
                        WRITE(16,2205) NBV(J),QNAM(I,J),QNPH(I,J)
 2205                   FORMAT(10X,I8,4X,F14.5,4X,F12.3)
                        QNPH(I,J)=QNPH(I,J)*DEG2RAD
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
            II = 1
         ENDIF
      ENDIF

C...  IF ANY EXTERNAL BARRIER BOUNDARIES WERE SPECIFIED, (NFLUXB=1)
C.....WRITE OUT EXTERNAL BARRIER BOUNDARY INFORMATION TO UNIT 16 FILE
C.....NOTE THAT THIS INFORMATION WAS READ IN FROM THE UNIT 14 FILE

      IF(NFLUXB.EQ.1) THEN

C.....WRITE OUT INFO ON SPECIFIED EXTERNAL BARRIER BOUNDARIES

         WRITE(16,1112)
         WRITE(16,2220)
 2220    FORMAT(//,1X,'EXTERNAL BARRIER BOUNDARY INFORMATION ',/)

C.......OUTPUT ELEVATION OF EXTERNAL BARRIER NODES ABOVE THE GEOID AND
C........THE COEFFICIENT OF FREE SURFACE SUPERCRITICAL FLOW AT
C........DESIGNATED EXTERNAL BARRIER BOUNDARY NODES TO UNIT 16
         WRITE(16,2224)
 2224    FORMAT(//,9X,'NODE',10X,'BARRIER HEIGHT',
     &        6X,'SUPER-CRIT. EXTERNAL BAR. COEF.',/)
         DO J=1,NVEL
            IF((LBCODEI(J).EQ.3).OR.(LBCODEI(J).EQ.13)
     &           .OR.(LBCODEI(J).EQ.23)) THEN
               WRITE(16,2225) NBV(J),BARLANHT(J),BARLANCFSP(J)
 2225          FORMAT(5X,I8,6X,F14.5,15X,F12.3)
            ENDIF
         END DO
      ENDIF

C...  IF ANY INTERNAL BARRIER BOUNDARIES WERE SPECIFIED, (NFLUXIB=1)
C.....WRITE INTERNAL BARRIER BOUNDARY INFORMATION TO UNIT 16 FILE

      IF(NFLUXIB.EQ.1) THEN

C.....WRITE OUT INFO ON SPECIFIED INTERNAL BARRIER BOUNDARIES

         WRITE(16,1112)
         WRITE(16,2320)
 2320    FORMAT(//,1X,'INTERNAL BARRIER BOUNDARY INFORMATION ',/)

C.......WRITE CONNECTION NODE NUMBER AND ELEVATION OF THE INTERNAL BARRIER
C........NODES ABOVE THE GEOID AND THE COEFFICIENTS OF FREE SURFACE SUPERCRITICAL
C........AND SUBCRITICAL FLOW AT DESIGNATED INTERNAL BARRIER BOUNDARY NODES
C........TO UNIT 16 (NOTE THAT THIS INFORMATION WAS INPUT FROM THE UNIT 14
C........FILE WITH BOUNDARY NODE INFORMATION)
         WRITE(16,2324)
 2324    FORMAT(//,9X,'NODE',6X,'CONNECTED NODE',6X,'BARRIER HEIGHT',
     &        4X,'SUB-CRIT. INT. BAR. COEF.',
     &        4X,'SUPER-CRIT. INT. BAR. COEF.',/)

         DO J=1,NVEL
            IF((LBCODEI(J).EQ.4).OR.(LBCODEI(J).EQ.24)) THEN

               WRITE(16,2325) NBV(J),IBCONN(J),BARINHT(J),
     &              BARINCFSB(J),BARINCFSP(J)
 2325          FORMAT(5X,I8,7X,I8,6X,F14.5,12X,F12.3,17X,F12.3)
            ENDIF
         END DO
      ENDIF

cjj   wm001 - begin add                        
C...  IF ANY INTERNAL BARRIER BOUNDARIES WITH CROSS BARRIER PIPES
C.....WERE SPECIFIED, (NFLUXIBP=1)
C.....WRITE INTERNAL BARRIER BOUNDARY INFORMATION WITH CROSS 
C.....BARRIER PIPE INFORMATION TO UNIT 16 FILE

      IF(NFLUXIBP.EQ.1) THEN

C.....WRITE OUT INFO ON SPECIFIED INTERNAL BARRIER BOUNDARIES

         WRITE(16,1112)
         WRITE(16,2326)
 2326    FORMAT(//,1X,'INTERNAL BARRIER BOUNDARY WITH CROSS BARRIER',
     &        ' PIPE INFORMATION ',/)

C.......WRITE CONNECTION NODE NUMBER AND ELEVATION OF THE INTERNAL BARRIER
C........NODES ABOVE THE GEOID AND THE COEFFICIENTS OF FREE SURFACE SUPERCRITICAL
C........AND SUBCRITICAL FLOW AT DESIGNATED INTERNAL BARRIER BOUNDARY NODES
C........IN ADDITION TO CROSS BARRIER PIPE CROWN HEIGHT, CROSS BARRIER PIPE
C........COEFFICIENT AND CROSS BARRIER PIPE DIAMETER TO UNIT 16
C........(NOTE THAT THIS INFORMATION WAS INPUT FROM THE UNIT 14 FILE WITH 
C........BOUNDARY NODE INFORMATION)
         WRITE(16,2327)
 2327    FORMAT(//,7X,'NODE',4X,'CONNECTED NODE',4X,'BARRIER HEIGHT',
     &        4X,'SUB-CRIT INT BAR COEF',
     &        4X,'SUPER-CRIT INT BAR COEF',
     &        4X,'PIPEHT  ',
     &        4X,'PIPECOEF',
     &        4X,'PIPEDIAM',/)
         DO J=1,NVEL
            IF((LBCODEI(J).EQ.5).OR.(LBCODEI(J).EQ.25)) THEN
               WRITE(16,2328) NBV(J),IBCONN(J),BARINHT(J),
     &              BARINCFSB(J),BARINCFSP(J),
     &              PIPEHT(J),PIPECOEF(J),PIPEDIAM(J)
 2328          FORMAT(3X,I8,5X,I8,4X,F14.5,8X,F12.3,12X,F12.3,
     &              2X,F10.5,2X,F10.5,2X,F10.5)
            ENDIF
         END DO
      ENDIF
cjj   wm001 - end add                        

C...  
C...  READ IN INFORMATION CONCERNING OUTPUT REQUIREMENTS FROM UNIT 15 AND
C...  OUTPUT THIS TO UNIT 16
C...  
      WRITE(16,1112)
      WRITE(16,3000)
 3000 FORMAT(//,1X,'OUTPUT INFORMATION WILL BE PROVIDED AS'
     &     ,' FOLLOWS :')

C...  
C...  INPUT INFORMATION FOR ELEVATION RECORDING STATIONS
C...  

C.... READ IN NOUTE,TOUTSE,TOUTFE,NSPOOLE : IF ABS(NOUTE)>0, INTERPOLATED
C.... ELEVATIONS AT ELEVATION STATIONS ARE SPOOLED TO UNIT 61 EVERY NSPOOLE
C.... TIME STEPS BETWEEN TIMES TOUTSE AND TOUTFE

      READ(15,*) NOUTE,TOUTSE,TOUTFE,NSPOOLE
      WRITE(16,3001) NOUTE
 3001 FORMAT(///,1X,'ELEVATION RECORDING STATION OUTPUT : ',
     &     //,5X,'NOUTE = ',I2)

C.... CHECK INPUT PARAMETER NOUTE

      IF(ABS(NOUTE).GT.2) THEN
         IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,3002)
         WRITE(16,3002)
 3002    FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &        //,1X,'YOUR SELECTION OF THE UNIT 15 INPUT PARAMETER',
     &        ' NOUTE',
     &        /,1X,'IS NOT AN ALLOWABLE VALUE.  CHECK YOUR INPUT!!')
         IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9973)
         WRITE(16,9973)
         STOP
      ENDIF

C.... IF STATION ELEVATION OUTPUT WILL NOT BE GENERATED

      IF(NOUTE.EQ.0) THEN
         WRITE(16,3003)
 3003    FORMAT(/,5X,'NO OUTPUT WILL BE SPOOLED AT ELEVATION ',
     &        'RECORDING STATIONS')
      ENDIF

C.... IF STATION ELEVATION OUTPUT WILL BE GENERATED

      IF(NOUTE.NE.0) THEN

C......COMPUTE NTCYSE, NTCYFE, WHICH = TOUTSE AND TOUTFE IN TIMESTEPS

         NTCYSE=INT((TOUTSE-STATIM)*(86400.D0/DTDP)+0.5d0)
         NTCYFE=INT((TOUTFE-STATIM)*(86400.D0/DTDP)+0.5d0)

         IF(NTCYFE.GT.NT) NTCYFE=NT

C......COMPUTE NTRSPE = THE NO. OF DATA SETS TO BE SPOOLED TO UNIT 61

         IF(NSPOOLE.EQ.0) NTRSPE=0
         IF(NSPOOLE.NE.0) NTRSPE=INT((NTCYFE-NTCYSE)/NSPOOLE)

C......WRITE TOUTSE,TOUTFE,NTCYSE,NTCYFE,NSPOOLE TO UNIT 16

         WRITE(16,3004) TOUTSE,NTCYSE,TOUTFE,NTCYFE,NSPOOLE
 3004    FORMAT(/,5X,'DATA RECORDS WILL START AFTER TOUTSE =',F8.3,
     &        ' DAY(S) RELATIVE',
     &        /,9X,'TO THE STARTING TIME OR',I9,
     &        ' TIME STEPS INTO THE SIMULATION',
     &        //,5X,'DATA RECORDS WILL STOP AFTER TOUTFE =',F8.3,
     &        ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &        I9,' TIME STEPS INTO THE SIMULATION',
     &        //,5X,'INFORMATION WILL BE SPOOLED TO UNIT 61 EVERY',
     &        ' NSPOOLE =',I8,' TIME STEPS')
         IF(ABS(NOUTE).EQ.1) WRITE(16,3005)
 3005    FORMAT(/,5X,'UNIT 61 FORMAT WILL BE ASCII')
         IF(ABS(NOUTE).EQ.2) WRITE(16,3006)
 3006    FORMAT(/,5X,'UNIT 61 FORMAT WILL BE BINARY')
      ENDIF

C.... REGARDLESS OF WHETHER NOUTE=0, READ IN THE NUMBER OF ELEVATION
C.... RECORDING STATIONS

      READ(15,*) NSTAE
      WRITE(16,3007) NSTAE
 3007 FORMAT(///,5X,'NUMBER OF INPUT ELEVATION RECORDING STATIONS = ',
     &     I5)

      IF(NSTAE.GT.0) THEN
         IF(ICS.EQ.1) WRITE(16,3008)
 3008    FORMAT(/,7X,'STATION #   ELEMENT',9X,'X',13X,'Y',/)
         IF(ICS.EQ.2) WRITE(16,3009)
 3009    FORMAT(/,5X,'STATION   ELEMENT',3X,'LAMBDA(DEG)',
     &        4X,'FEA(DEG)',10X,'XCP',12X,'YCP',/)
         MNSTAE = NSTAE
      ENDIF
      IF (NSTAE.EQ.0) MNSTAE = 1


C     Allocate arrays dimensioned by MNSTAE
      call alloc_main7()


C.... INPUT COORDINATES OF ELEVATION RECORDING STATIONS THEN COMPUTE
C.... THE ELEMENT NO. THE STATION LIES IN

      DO I=1,NSTAE
         NNE(I)=0
         IF(ICS.EQ.1) THEN
            READ(15,*) XEL(I),YEL(I)
         ELSE
            READ(15,*) SLEL(I),SFEL(I)
            SLEL(I)=SLEL(I)*DEG2RAD
            SFEL(I)=SFEL(I)*DEG2RAD
            CALL CPP(XEL(I),YEL(I),SLEL(I),SFEL(I),SLAM0,SFEA0)
         ENDIF
         AEMIN=1.0E+25
         KMIN=0
         DO K=1,NE
            N1=NM(K,1)
            N2=NM(K,2)
            N3=NM(K,3)
            X1=X(N1)
            X2=X(N2)
            X3=X(N3)
            X4=XEL(I)
            Y1=Y(N1)
            Y2=Y(N2)
            Y3=Y(N3)
            Y4=YEL(I)
            A1=(X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4)
            A2=(X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1)
            A3=(Y4-Y1)*(X2-X1)-(X4-X1)*(Y2-Y1)
            AA=ABS(A1)+ABS(A2)+ABS(A3)
            AE=ABS(AA-AREAS(K))/AREAS(K)
            IF(AE.LT.AEMIN) THEN
               AEMIN=AE
               KMIN=K
            ENDIF
            IF(AE.LT.1.0E-5) NNE(I)=K
         END DO

         IF(NNE(I).EQ.0) THEN
            IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9784) I
            WRITE(16,9784) I
 9784       FORMAT(///,1X,'!!!!!!!!!!  WARNING - NONFATAL ',
     &           'INPUT ERROR  !!!!!!!!!',//
     &           ,1X,'ELEVATION RECORDING STATION ',I6,' DOES NOT LIE',
     &           ' WITHIN ANY ELEMENT IN THE DEFINED',
     &           /,1X,'COMPUTATIONAL DOMAIN,   PLEASE CHECK THE INPUT',
     &           ' COORDINATES FOR THIS STATION')
            IF(NFOVER.EQ.1) THEN
               IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9790) AEMIN
               WRITE(16,9790) AEMIN
 9790          FORMAT(/,1X,'PROGRAM WILL ESTIMATE NEAREST ELEMENT',
     &              /,1X,'PROXIMITY INDEX FOR THIS STATION EQUALS ',E15.6,
     &              //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)
               NNE(I)=KMIN
            ELSE
               IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9791) AEMIN
               WRITE(16,9791) AEMIN
 9791          FORMAT(/,1X,'PROGRAM WILL NOT CORRECT ERROR ',
     &              'SINCE NON-FATAL ERROR OVERIDE OPTION, NFOVER,',
     &          /,1X,'HAS BEEN SELECTED EQUAL TO 0',
     &          /,1X,'PROXIMITY INDEX FOR THIS STATION EQUALS ',E15.6,
     &          //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',
     &              //)
               STOP
            ENDIF
         ENDIF

         IF(ICS.EQ.1) THEN
            WRITE(16,1880) I,NNE(I),XEL(I),YEL(I)
 1880       FORMAT(8X,I3,6X,I7,2(2X,F14.2))
         ELSE
            WRITE(16,1883) I,NNE(I),SLEL(I)*RAD2DEG,
     &           SFEL(I)*RAD2DEG,XEL(I),YEL(I)
 1883       FORMAT(6X,I3,4X,I7,2(2X,F13.8),2X,2(1X,F13.2))
         ENDIF

C.... PRE-COMPUTE INFORMATION REQUIRED TO INTERPOLATE AT ELEV. RECORDING STATIONS

         N1=NM(NNE(I),1)
         N2=NM(NNE(I),2)
         N3=NM(NNE(I),3)
         X1=X(N1)
         X2=X(N2)
         X3=X(N3)
         X4=XEL(I)
         Y1=Y(N1)
         Y2=Y(N2)
         Y3=Y(N3)
         Y4=YEL(I)
         STAIE1(I)=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))/AREAS(NNE(I))
         STAIE2(I)=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))/AREAS(NNE(I))
         STAIE3(I)=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))/AREAS(NNE(I))

      END DO

C...  
C...  INPUT INFORMATION FOR VELOCITY RECORDING STATIONS
C...  

C.... READ IN NOUTV,TOUTSV,TOUTFV,NSPOOLV : IF NOUTV<>0,INTERPOLATED VELOCITIES AT
C.... VELOCITY STATIONS ARE SPOOLED TO UNIT 62 EVERY NSPOOLV TIME STEPS BETWEEN
C.... TIMES TOUTSV AND TOUTFV; IF ABS(NOUTV)=2, OUTPUT WILL BE BINARY

      READ(15,*) NOUTV,TOUTSV,TOUTFV,NSPOOLV
      WRITE(16,3101) NOUTV
 3101 FORMAT(////,1X,'VELOCITY RECORDING STATION OUTPUT : ',
     &     //,5X,'NOUTV = ',I2)

C.... CHECK INPUT PARAMETER NOUTV

      IF(ABS(NOUTV).GT.2) THEN
         IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,3102)
         WRITE(16,3102)
 3102    FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &        //,1X,'YOUR SELECTION OF THE UNIT 15 INPUT PARAMETER',
     &        ' NOUTV',
     &        /,1X,'IS NOT AN ALLOWABLE VALUE.  CHECK YOUR INPUT!!')
         IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9973)
         WRITE(16,9973)
         STOP
      ENDIF

C.... IF STATION VELOCITY OUTPUT WILL NOT BE GENERATED

      IF(NOUTV.EQ.0) THEN
         WRITE(16,3103)
 3103    FORMAT(///,5X,'NO OUTPUT WILL BE SPOOLED AT VELOCITY',
     &        ' RECORDING STATIONS')
      ENDIF

C.... IF STATION VELOCITY OUTPUT WILL BE GENERATED

      IF(NOUTV.NE.0) THEN

C......COMPUTE NTCYSV, NTCYFV, WHICH = TOUTSV AND TOUTFV IN TIME STEPS

         NTCYSV=INT((TOUTSV-STATIM)*(86400.D0/DTDP) + 0.5d0)
         NTCYFV=INT((TOUTFV-STATIM)*(86400.D0/DTDP) + 0.5d0)
         IF(NTCYFV.GT.NT) NTCYFV=NT

C......CALCULATE NTRSPV = THE NO. OF DATA SETS TO BE SPOOLED TO UNIT 62

         IF(NSPOOLV.EQ.0) NTRSPV=0
         IF(NSPOOLV.NE.0) NTRSPV=INT((NTCYFV-NTCYSV)/NSPOOLV)

C......WRITE NOUTV,TOUTSV,TOUTFV,NTCYSV,NTCYFV,NSPOOLV TO UNIT 16

         WRITE(16,3104) TOUTSV,NTCYSV,TOUTFV,NTCYFV,NSPOOLV
 3104    FORMAT(/,5X,'DATA RECORDS WILL START AFTER TOUTSV =',F8.3,
     &        ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &        I9,' TIME STEPS INTO THE SIMULATION',
     &        //,5X,'DATA RECORDS WILL STOP AFTER TOUTFV =',F8.3,
     &        ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &        I9,' TIME STEPS INTO THE SIMULATION',
     &        //,5X,'INFORMATION WILL BE SPOOLED TO UNIT 62 EVERY ',
     &        ' NSPOOLV =',I8,' TIME STEPS')
         IF(ABS(NOUTV).EQ.1) WRITE(16,3105)
 3105    FORMAT(/,5X,'UNIT 62 FORMAT WILL BE ASCII')
         IF(ABS(NOUTV).EQ.2) WRITE(16,3106)
 3106    FORMAT(/,5X,'UNIT 62 FORMAT WILL BE BINARY')
      ENDIF

C.... REGARDLESS OF WHETHER NOUTV=0, READ IN THE NUMBER OF VELOCITY
C.... RECORDING STATIONS

      READ(15,*) NSTAV
      WRITE(16,3107) NSTAV
 3107 FORMAT(////,5X,'NUMBER OF INPUT VELOCITY RECORDING STATIONS = ',
     &     I5)

      IF(NSTAV.GT.0) THEN
         IF(ICS.EQ.1) WRITE(16,3108)
 3108    FORMAT(/,7X,'STATION #   ELEMENT',9X,'X',13X,'Y',/)
         IF(ICS.EQ.2) WRITE(16,3109)
 3109    FORMAT(/,5X,'STATION   ELEMENT',3X,'LAMBDA(DEG)',
     &        4X,'FEA(DEG)',10X,'XCP',12X,'YCP',/)
         MNSTAV = NSTAV
      ENDIF
      IF (NSTAV.EQ.0) MNSTAV = 1

C     Allocate arrays dimensioned by MNSTAV
      call alloc_main8()

C.... INPUT COORDINATES OF VELOCITY RECORDING STATIONS
C.... THEN COMPUTE ELEMENT NO. WITHIN WHICH STATION LIES

      DO I=1,NSTAV
         NNV(I)=0
         IF(ICS.EQ.1) THEN
            READ(15,*) XEV(I),YEV(I)
         ELSE
            READ(15,*) SLEV(I),SFEV(I)
            SLEV(I)=SLEV(I)*DEG2RAD
            SFEV(I)=SFEV(I)*DEG2RAD
            CALL CPP(XEV(I),YEV(I),SLEV(I),SFEV(I),SLAM0,SFEA0)
         ENDIF
         AEMIN=1.0E+25
         KMIN=0
         DO K=1,NE
            N1=NM(K,1)
            N2=NM(K,2)
            N3=NM(K,3)
            X1=X(N1)
            X2=X(N2)
            X3=X(N3)
            X4=XEV(I)
            Y1=Y(N1)
            Y2=Y(N2)
            Y3=Y(N3)
            Y4=YEV(I)
            A1=(X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4)
            A2=(X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1)
            A3=(Y4-Y1)*(X2-X1)-(X4-X1)*(Y2-Y1)
            AA=ABS(A1)+ABS(A2)+ABS(A3)
            AE=ABS(AA-AREAS(K))/AREAS(K)
            IF(AE.LT.AEMIN) THEN
               AEMIN=AE
               KMIN=K
            ENDIF
            IF(AE.LT.1.0E-5) NNV(I)=K
         END DO

         IF(NNV(I).EQ.0) THEN
            IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9786) I
            WRITE(16,9786) I
 9786       FORMAT(///,1X,'!!!!!!!!!!  WARNING - NONFATAL ',
     &           'INPUT ERROR  !!!!!!!!!',//
     &           ,1X,'VELOCITY RECORDING STATION ',I6,' DOES NOT LIE'
     &           ,' WITHIN ANY ELEMENT IN THE DEFINED',
     &           /,1X,'COMPUTATIONAL DOMAIN,   PLEASE CHECK THE INPUT'
     &           ,' COORDINATES FOR THIS STATION')
            IF(NFOVER.EQ.1) THEN
               IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9790) AEMIN
               WRITE(16,9790) AEMIN
               NNV(I)=KMIN
            ELSE
               IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9791) AEMIN
               WRITE(16,9791) AEMIN
               STOP
            ENDIF
         ENDIF

         IF(ICS.EQ.1) THEN
            WRITE(16,1880) I,NNV(I),XEV(I),YEV(I)
         ELSE
            WRITE(16,1883) I,NNV(I),SLEV(I)*RAD2DEG,SFEV(I)*RAD2DEG,
     &           XEV(I),YEV(I)
         ENDIF

C.... PRE-COMPUTE INFORMATION REQUIRED TO INTERPOLATE AT VEL. RECORDING STATIONS

         N1=NM(NNV(I),1)
         N2=NM(NNV(I),2)
         N3=NM(NNV(I),3)
         X1=X(N1)
         X2=X(N2)
         X3=X(N3)
         X4=XEV(I)
         Y1=Y(N1)
         Y2=Y(N2)
         Y3=Y(N3)
         Y4=YEV(I)
         STAIV1(I)=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))/AREAS(NNV(I))
         STAIV2(I)=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))/AREAS(NNV(I))
         STAIV3(I)=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))/AREAS(NNV(I))

      END DO

C...  
C...  IF TRANSPORT IS INCLUDED IN THE RUN, INPUT INFORMATION FOR CONCENTRATION
C...  RECORDING STATIONS
C...  
      NOUTC=0
      IF(IM.EQ.10) THEN

C.....READ IN NOUTC,TOUTSC,TOUTFC,NSPOOLC : IF NOUTC<>0,INTERPOLATED
C.....CONCENTRATIONS ARE SPOOLED TO UNIT 81 EVERY NSPOOLC TIME STEPS
C.....BETWEEN TIMES TOUTSC AND TOUTFC; IF ABS(NOUTC)=2, OUTPUT WILL BE BINARY

         READ(15,*) NOUTC,TOUTSC,TOUTFC,NSPOOLC
         WRITE(16,3201) NOUTC
 3201    FORMAT(///,1X,'CONCENTRATION RECORDING STATION OUTPUT : ',
     &        //,5X,'NOUTC = ',I2)

C.....CHECK INPUT PARAMETER NOUTC

         IF(ABS(NOUTC).GT.2) THEN
            IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,3202)
            WRITE(16,3202)
 3202       FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &           //,1X,'YOUR SELECTION OF THE UNIT 15 INPUT PARAMETER',
     &           ' NOUTC',
     &           /,1X,'IS NOT AN ALLOWABLE VALUE.  CHECK YOUR INPUT!!')
            IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9973)
            WRITE(16,9973)
            STOP
         ENDIF

C.....IF STATION CONCENTRATION OUTPUT WILL NOT BE GENERATED

         IF(NOUTC.EQ.0) THEN
            WRITE(16,3203)
 3203       FORMAT(/,5X,'NO OUTPUT WILL BE SPOOLED AT CONCENTRATION',
     &           ' RECORDING STATIONS')
         ENDIF

C.....IF STATION CONCENTRATION OUTPUT WILL BE GENERATED

         NSTAC = 0
         IF(NOUTC.NE.0) THEN

C.......COMPUTE NTCYSC, NTCYFC, WHICH = TOUTSC AND TOUTFC IN TIMESTEPS

            NTCYSC=INT((TOUTSC-STATIM)*(86400.D0/DTDP) + 0.5d0)
            NTCYFC=INT((TOUTFC-STATIM)*(86400.D0/DTDP) + 0.5d0)
            IF(NTCYFC.GT.NT) NTCYFC=NT

C.......COMPUTE NTRSPC = THE NO. OF DATA SETS TO BE SPOOLED TO UNIT 81

            IF(NSPOOLC.EQ.0) NTRSPC=0
            IF(NSPOOLC.NE.0) NTRSPC=INT((NTCYFC-NTCYSC)/NSPOOLC)

C.......WRITE TOUTSC,TOUTFC,NTCYSC,NTCYFC,NSPOOLC TO UNIT 16

            WRITE(16,3204) TOUTSC,NTCYSC,TOUTFC,NTCYFC,NSPOOLC
 3204       FORMAT(/,5X,'DATA RECORDS WILL START AFTER TOUTSC =',F8.3,
     &           ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &           I9,' TIME STEPS INTO THE SIMULATION',
     &           //,5X,'DATA RECORDS WILL STOP AFTER TOUTFC =',F8.3,
     &           ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &           I9,' TIME STEPS INTO THE SIMULATION',
     &           //,5X,'INFORMATION WILL BE SPOOLED TO UNIT 81 EVERY',
     &           ' NSPOOLC =',I8,' TIME STEPS')
            IF(ABS(NOUTC).EQ.1) WRITE(16,3205)
 3205       FORMAT(/,5X,'UNIT 81 FORMAT WILL BE ASCII')
            IF(ABS(NOUTC).EQ.2) WRITE(16,3206)
 3206       FORMAT(/,5X,'UNIT 81 FORMAT WILL BE BINARY')
         ENDIF

C.....REGARDLESS OF WHETHER NOUTC=0, READ IN THE NUMBER OF CONCENTRATION
C.....RECORDING STATIONS

         READ(15,*) NSTAC
         WRITE(16,3207) NSTAC
 3207    FORMAT(///,5X,'NUMBER OF INPUT CONCENTRATION RECORDING ',
     &        'STATIONS = ',I5)

         IF(NSTAC.GT.0) THEN
            IF(ICS.EQ.1) WRITE(16,3208)
 3208       FORMAT(/,7X,'STATION #   ELEMENT',9X,'X',13X,'Y',/)
            IF(ICS.EQ.2) WRITE(16,3209)
 3209       FORMAT(/,5X,'STATION   ELEMENT',3X,'LAMBDA(DEG)',
     &           4X,'FEA(DEG)',10X,'XCP',12X,'YCP',/)
            MNSTAC = NSTAC
         ENDIF

C     Allocate arrays dimensioned by MNSTAC
         call alloc_main9()

C.....INPUT COORDINATES OF CONCENTRATION RECORDING STATIONS
C.....THEN COMPUTE ELEMENT NO. WITHIN WHICH STATION LIES

         DO I=1,NSTAC
            NNC(I)=0
            IF(ICS.EQ.1) THEN
               READ(15,*) XEC(I),YEC(I)
            ELSE
               READ(15,*) SLEC(I),SFEC(I)
               SLEC(I)=SLEC(I)*DEG2RAD
               SFEC(I)=SFEC(I)*DEG2RAD
               CALL CPP(XEC(I),YEC(I),SLEC(I),SFEC(I),SLAM0,SFEA0)
            ENDIF
            AEMIN=1.0E+25
            KMIN=0
            DO K=1,NE
               N1=NM(K,1)
               N2=NM(K,2)
               N3=NM(K,3)
               X1=X(N1)
               X2=X(N2)
               X3=X(N3)
               X4=XEC(I)
               Y1=Y(N1)
               Y2=Y(N2)
               Y3=Y(N3)
               Y4=YEC(I)
               A1=(X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4)
               A2=(X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1)
               A3=(Y4-Y1)*(X2-X1)-(X4-X1)*(Y2-Y1)
               AA=ABS(A1)+ABS(A2)+ABS(A3)
               AE=ABS(AA-AREAS(K))/AREAS(K)
               IF(AE.LT.AEMIN) THEN
                  AEMIN=AE
                  KMIN=K
               ENDIF
               IF(AE.LT.1.0E-5) NNC(I)=K
            END DO

            IF(NNC(I).EQ.0) THEN
               IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9785) I
               WRITE(16,9785) I
 9785          FORMAT(///,1X,'!!!!!!!!!!  WARNING - NONFATAL INPUT ERROR ',
     &          '!!!!!!!!!',//,
     &          ' CONCENTRATION RECORDING STATION ',I6,' DOES NOT LIE'
     &          ,' WITHIN ANY ELEMENT IN THE DEFINED',/,
     &          ' COMPUTATIONAL DOMAIN,   PLEASE CHECK THE INPUT',
     &          ' COORDINATES FOR THIS STATION')
               IF(NFOVER.EQ.1) THEN
                  IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9790) AEMIN
                  WRITE(16,9790) AEMIN
                  NNC(I)=KMIN
               ELSE
                  IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9791) AEMIN
                  WRITE(16,9791) AEMIN
                  STOP
               ENDIF
            ENDIF

            IF(ICS.EQ.1) THEN
               WRITE(16,1880) I,NNC(I),XEC(I),YEC(I)
            ELSE
               WRITE(16,1883) I,NNC(I),SLEC(I)*RAD2DEG,
     &              SFEC(I)*RAD2DEG,XEC(I),YEC(I)
            ENDIF

C.....PRE-COMPUTE INFORMATION REQUIRED TO INTERPOLATE AT CONCENTRATION
C.....RECORDING STATIONS

            N1=NM(NNC(I),1)
            N2=NM(NNC(I),2)
            N3=NM(NNC(I),3)
            X1=X(N1)
            X2=X(N2)
            X3=X(N3)
            X4=XEL(I)
            Y1=Y(N1)
            Y2=Y(N2)
            Y3=Y(N3)
            Y4=YEL(I)
            STAIC1(I)=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))/AREAS(NNC(I))
            STAIC2(I)=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))/AREAS(NNC(I))
            STAIC3(I)=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))/AREAS(NNC(I))

         END DO
      ENDIF
      IF (NSTAC.EQ.0) MNSTAC = 1

C...  
C...  IF METEOROLOICAL FORCING IS INCLUDED IN THE RUN, INPUT INFORMATION FOR MET
C...  RECORDING STATIONS - OUTPUT
C...  
      NOUTM=0
      NSTAM = 0
C     
      IF(NWS.NE.0) THEN

C.....READ IN NOUTM,TOUTSM,TOUTFM,NSPOOLM : IF NOUTM<>0,INTERPOLATED
C.....MET DATA ARE SPOOLED TO UNITS 71&72 EVERY NSPOOLM TIME STEPS
C.....BETWEEN TIMES TOUTSM AND TOUTFM; IF ABS(NOUTM)=2, OUTPUT WILL BE BINARY

         READ(15,*) NOUTM,TOUTSM,TOUTFM,NSPOOLM
         WRITE(16,3211) NOUTM
 3211    FORMAT(///,1X,'METEOROLOGICAL RECORDING STATION OUTPUT : ',
     &        //,5X,'NOUTM = ',I2)

C.....CHECK INPUT PARAMETER NOUTM

         IF(ABS(NOUTM).GT.2) THEN
            IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,3212)
            WRITE(16,3202)
 3212       FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &           //,1X,'YOUR SELECTION OF THE UNIT 15 INPUT PARAMETER',
     &           ' NOUTC',
     &           /,1X,'IS NOT AN ALLOWABLE VALUE.  CHECK YOUR INPUT!!')
            IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9973)
            WRITE(16,9973)
            STOP
         ENDIF

C.....IF STATION METEOROLOGICAL OUTPUT WILL NOT BE GENERATED

         IF(NOUTM.EQ.0) THEN
            WRITE(16,3213)
 3213       FORMAT(/,5X,'NO OUTPUT WILL BE SPOOLED AT METEOROLOGICAL',
     &           ' RECORDING STATIONS')
         ENDIF

C.....IF STATION MET OUTPUT WILL BE GENERATED

         IF(NOUTM.NE.0) THEN

C.......COMPUTE NTCYSM, NTCYFM, WHICH = TOUTSM AND TOUTFM IN TIMESTEPS

            NTCYSM=INT((TOUTSM-STATIM)*(86400.D0/DTDP) + 0.5d0)
            NTCYFM=INT((TOUTFM-STATIM)*(86400.D0/DTDP) + 0.5d0)
            IF(NTCYFM.GT.NT) NTCYFM=NT

C.......COMPUTE NTRSPM = THE NO. OF DATA SETS TO BE SPOOLED TO UNITS 71&72

            IF(NSPOOLM.EQ.0) NTRSPM=0
            IF(NSPOOLM.NE.0) NTRSPM=INT((NTCYFM-NTCYSM)/NSPOOLM)

C.......WRITE TOUTSM,TOUTFM,NTCYSM,NTCYFM,NSPOOLM TO UNIT 16

            WRITE(16,3214) TOUTSM,NTCYSM,TOUTFM,NTCYFM,NSPOOLM
 3214       FORMAT(/,5X,'DATA RECORDS WILL START AFTER TOUTSM =',F8.3,
     &           ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &           I9,' TIME STEPS INTO THE SIMULATION',
     &           //,5X,'DATA RECORDS WILL STOP AFTER TOUTFM =',F8.3,
     &           ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &           I9,' TIME STEPS INTO THE SIMULATION',
     &           //,5X,'INFORMATION WILL BE SPOOLED TO UNITS 71&72',
     &           ' EVERY NSPOOLM =',I8,' TIME STEPS')
            IF(ABS(NOUTM).EQ.1) WRITE(16,3215)
 3215       FORMAT(/,5X,'UNITS 71&72 FORMAT WILL BE ASCII')
            IF(ABS(NOUTM).EQ.2) WRITE(16,3216)
 3216       FORMAT(/,5X,'UNITS 71&72 FORMAT WILL BE BINARY')
         ENDIF

C.....REGARDLESS OF WHETHER NOUTM=0, READ IN THE NUMBER OF METEOROLOGICAL
C.....RECORDING STATIONS

         READ(15,*) NSTAM
         WRITE(16,3217) NSTAM
 3217    FORMAT(///,5X,'NUMBER OF INPUT METEOROLOGICAL RECORDING ',
     &        'STATIONS = ',I5)

         IF(NSTAM.GT.0) THEN
            IF(ICS.EQ.1) WRITE(16,3218)
 3218       FORMAT(/,7X,'STATION #   ELEMENT',9X,'X',13X,'Y',/)
            IF(ICS.EQ.2) WRITE(16,3219)
 3219       FORMAT(/,5X,'STATION   ELEMENT',3X,'LAMBDA(DEG)',
     &           4X,'FEA(DEG)',10X,'XCP',12X,'YCP',/)
            MNSTAM = NSTAM
         ENDIF

C     Allocate arrays dimensioned by MNSTAM
         call alloc_main10()

C.....INPUT COORDINATES OF METEOROLOGICAL RECORDING STATIONS
C.....THEN COMPUTE ELEMENT NO. WITHIN WHICH STATION LIES

         DO I=1,NSTAM
            NNM(I)=0
            IF(ICS.EQ.1) THEN
               READ(15,*) XEM(I),YEM(I)
            ELSE
               READ(15,*) SLEM(I),SFEM(I)
               SLEM(I)=SLEM(I)*DEG2RAD
               SFEM(I)=SFEM(I)*DEG2RAD
               CALL CPP(XEM(I),YEM(I),SLEM(I),SFEM(I),SLAM0,SFEA0)
            ENDIF
            AEMIN=1.0E+25
            KMIN=0
            DO K=1,NE
               N1=NM(K,1)
               N2=NM(K,2)
               N3=NM(K,3)
               X1=X(N1)
               X2=X(N2)
               X3=X(N3)
               X4=XEM(I)
               Y1=Y(N1)
               Y2=Y(N2)
               Y3=Y(N3)
               Y4=YEM(I)
               A1=(X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4)
               A2=(X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1)
               A3=(Y4-Y1)*(X2-X1)-(X4-X1)*(Y2-Y1)
               AA=ABS(A1)+ABS(A2)+ABS(A3)
               AE=ABS(AA-AREAS(K))/AREAS(K)
               IF(AE.LT.AEMIN) THEN
                  AEMIN=AE
                  KMIN=K
               ENDIF
               IF(AE.LT.1.0E-5) NNM(I)=K
            END DO

            IF(NNM(I).EQ.0) THEN
               IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9942) I
               WRITE(16,9942) I
 9942          FORMAT(///,1X,'!!!!!!!!!!  WARNING - NONFATAL INPUT ERROR ',
     &          '!!!!!!!!!',//,
     &          ' METEOROLOGICAL RECORDING STATION ',I6,' DOES NOT LIE'
     &          ,' WITHIN ANY ELEMENT IN THE DEFINED',/,
     &          ' COMPUTATIONAL DOMAIN,   PLEASE CHECK THE INPUT',
     &          ' COORDINATES FOR THIS STATION')
               IF(NFOVER.EQ.1) THEN
                  IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9790) AEMIN
                  WRITE(16,9790) AEMIN
                  NNM(I)=KMIN
               ELSE
                  IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9791) AEMIN
                  WRITE(16,9791) AEMIN
                  STOP
               ENDIF
            ENDIF

            IF(ICS.EQ.1) THEN
               WRITE(16,1880) I,NNM(I),XEM(I),YEM(I)
            ELSE
               WRITE(16,1883) I,NNM(I),SLEM(I)*RAD2DEG,
     &              SFEM(I)*RAD2DEG,XEM(I),YEM(I)
            ENDIF

C.....PRE-COMPUTE INFORMATION REQUIRED TO INTERPOLATE AT METEOROLOGICAL
C.....RECORDING STATIONS

            N1=NM(NNM(I),1)
            N2=NM(NNM(I),2)
            N3=NM(NNM(I),3)
            X1=X(N1)
            X2=X(N2)
            X3=X(N3)
            X4=XEM(I)
            Y1=Y(N1)
            Y2=Y(N2)
            Y3=Y(N3)
            Y4=YEM(I)
            STAIM1(I)=((X4-X3)*(Y2-Y3)+(X2-X3)*(Y3-Y4))/AREAS(NNM(I))
            STAIM2(I)=((X4-X1)*(Y3-Y1)-(Y4-Y1)*(X3-X1))/AREAS(NNM(I))
            STAIM3(I)=(-(X4-X1)*(Y2-Y1)+(Y4-Y1)*(X2-X1))/AREAS(NNM(I))

         END DO
      ENDIF
      IF (NSTAM.EQ.0) MNSTAM = 1

C...  
C...  INPUT INFORMATION ABOUT GLOBAL ELEVATION DATA OUTPUT
C...  

C.... READ IN NOUTGE,TOUTSGE,TOUTFGE,NSPOOLGE : IF NOUTGE<>0, GLOBAL ELEV.
C.... OUTPUT IS SPOOLED TO UNIT 63 EVERY NSPOOLGE TIME STEPS BETWEEN
C.... TIMES TOUTSGE AND TOUTFGE; IF ABS(NOUTGE)=2, OUTPUT WILL BE BINARY

      READ(15,*) NOUTGE,TOUTSGE,TOUTFGE,NSPOOLGE
      WRITE(16,3301) NOUTGE
 3301 FORMAT(////,1X,'GLOBAL NODAL ELEVATION INFORMATION OUTPUT: ',
     &     //,5X,'NOUTGE = ',I2)

C.... CHECK INPUT PARAMETER NOUTGE

      IF(ABS(NOUTGE).GT.2) THEN
         IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,3302)
         WRITE(16,3302)
 3302    FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &        //,1X,'YOUR SELECTION OF THE UNIT 15 INPUT PARAMETER',
     &        ' NOUTGE',
     &        /,1X,'IS NOT AN ALLOWABLE VALUE.  CHECK YOUR INPUT!!')
         IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9973)
         WRITE(16,9973)
         STOP
      ENDIF

C.... IF GLOBAL ELEVATION OUTPUT WILL NOT BE GENERATED

      IF(NOUTGE.EQ.0) THEN
         WRITE(16,3303)
 3303    FORMAT(///,5X,'NO GLOBAL ELEVATION OUTPUT WILL BE SPOOLED')
      ENDIF

C.... IF GLOBAL ELEVATION OUTPUT WILL BE GENERATED

      IF(NOUTGE.NE.0) THEN

C......COMPUTE NTCYSGE, NTCYFGE, WHICH = TOUTSGE AND TOUTFGE IN TIMESTEPS

         NTCYSGE=INT((TOUTSGE-STATIM)*(86400.D0/DTDP) + 0.5d0)
         NTCYFGE=INT((TOUTFGE-STATIM)*(86400.D0/DTDP) + 0.5d0)
         IF(NTCYFGE.GT.NT) NTCYFGE=NT

C......CALCULATE NDSETSE = THE # OF DATA SETS TO BE SPOOLED TO UNIT 63

         IF(NSPOOLGE.EQ.0) NDSETSE=0
         IF(NSPOOLGE.NE.0) NDSETSE=INT((NTCYFGE-NTCYSGE)/NSPOOLGE)

C......WRITE NOUTGE,TOUTSGE,TOUTFGE,NTCYSGE,NTCYFGE,NSPOOLGE TO UNIT 16

         WRITE(16,3304) TOUTSGE,NTCYSGE,TOUTFGE,NTCYFGE,NSPOOLGE
 3304    FORMAT(/,5X,'DATA RECORDS WILL START AFTER TOUTSGE =',F8.3,
     &        ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &        I9,' TIME STEPS INTO THE SIMULATION',
     &        //,5X,'DATA RECORDS WILL STOP AFTER TOUTFGE =',F8.3,
     &        ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &        I9,' TIME STEPS INTO THE SIMULATION',
     &        //,5X,'INFORMATION WILL BE SPOOLED TO UNIT 63 EVERY ',
     &        'NSPOOLGE =',I8,' TIME STEPS')
         IF(ABS(NOUTGE).EQ.1) WRITE(16,3305)
 3305    FORMAT(/,5X,'UNIT 63 FORMAT WILL BE ASCII')
         IF(ABS(NOUTGE).EQ.2) WRITE(16,3306)
 3306    FORMAT(/,5X,'UNIT 63 FORMAT WILL BE BINARY')
      ENDIF

C...  
C...  INPUT INFORMATION ABOUT GLOBAL VELOCITY DATA OUTPUT
C...  

C.... READ IN NOUTGV,TOUTSGV,TOUTFGV,NSPOOLGV : IF NOUTGV<>0, GLOBAL VEL.
C.... OUTPUT IS SPOOLED TO UNIT 64 EVERY NSPOOLGV TIME STEPS BETWEEN
C.... TIMES TOUTSGV AND TOUTFGV; IF ABS(NOUTGV)=2, OUTPUT WILL BE BINARY

      READ(15,*) NOUTGV,TOUTSGV,TOUTFGV,NSPOOLGV
      WRITE(16,3351) NOUTGV
 3351 FORMAT(////,1X,'GLOBAL NODAL VELOCITY INFORMATION OUTPUT : ',
     &     //,5X,'NOUTGV = ',I2)

C.... CHECK INPUT PARAMETER NOUTGV

      IF(ABS(NOUTGV).GT.2) THEN
         IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,3352)
         WRITE(16,3352)
 3352    FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &        //,1X,'YOUR SELECTION OF THE UNIT 15 INPUT PARAMETER',
     &        ' NOUTGV',
     &        /,1X,'IS NOT AN ALLOWABLE VALUE.  CHECK YOUR INPUT!!')
         IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9973)
         WRITE(16,9973)
         STOP
      ENDIF

C.... IF GLOBAL VELOCITY OUTPUT WILL NOT BE GENERATED

      IF(NOUTGV.EQ.0) THEN
         WRITE(16,3353)
 3353    FORMAT(///,5X,'NO GLOBAL VELOCITY OUTPUT WILL BE SPOOLED')
      ENDIF

C.... IF GLOBAL VELOCITY OUTPUT WILL BE GENERATED

      IF(NOUTGV.NE.0) THEN

C......COMPUTE NTCYSGV, NTCYFGV, WHICH = TOUTSGV AND TOUTFGV IN TIMESTEPS

         NTCYSGV=INT((TOUTSGV-STATIM)*(86400.D0/DTDP) + 0.5d0)
         NTCYFGV=INT((TOUTFGV-STATIM)*(86400.D0/DTDP) + 0.5d0)
         IF(NTCYFGV.GT.NT) NTCYFGV=NT

C......CALCULATE NDSETSV = THE # OF DATA SETS TO BE SPOOLED TO UNIT 64

         IF(NSPOOLGV.EQ.0) NDSETSV=0
         IF(NSPOOLGV.NE.0) NDSETSV=INT((NTCYFGV-NTCYSGV)/NSPOOLGV)

C......WRITE NOUTGV,TOUTSGV,TOUTFGV,NTCYSGV,NTCYFGV,NSPOOLGV TO UNIT 16

         WRITE(16,3354) TOUTSGV,NTCYSGV,TOUTFGV,NTCYFGV,NSPOOLGV
 3354    FORMAT(/,5X,'DATA RECORDS WILL START AFTER TOUTSGV =',F8.3,
     &        ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &        I9,' TIME STEPS INTO THE SIMULATION',
     &        //,5X,'DATA RECORDS WILL STOP AFTER TOUTFGV =',F8.3,
     &        ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &        I9,' TIME STEPS INTO THE SIMULATION',
     &        //,5X,'INFORMATION WILL BE SPOOLED TO UNIT 64 EVERY ',
     &        'NSPOOLGV =',I8,' TIME STEPS')
         IF(ABS(NOUTGV).EQ.1) WRITE(16,3355)
 3355    FORMAT(/,5X,'UNIT 64 FORMAT WILL BE ASCII')
         IF(ABS(NOUTGV).EQ.2) WRITE(16,3356)
 3356    FORMAT(/,5X,'UNIT 64 FORMAT WILL BE BINARY')
      ENDIF

C...  
C...  IF TRANSPORT IS INCLUDED IN THE RUN, INPUT INFORMATION ABOUT GLOBAL
C...  CONCENTRATION DATA OUTPUT
C...  
      NOUTGC=0
      IF(IM.EQ.10) THEN

C.....READ IN NOUTGC,TOUTSGC,TOUTFGC,NSPOOLGC : IF NOUTGC<>0, GLOBAL
C.....CONCENTRATION OUTPUT IS SPOOLED TO UNIT 73 EVERY NSPOOLGC TIME STEPS
C.....BETWEEN TIMES TOUTSGC AND TOUTFGC; IF ABS(NOUTGC)=2, OUTPUT WILL BE BINARY

         READ(15,*) NOUTGC,TOUTSGC,TOUTFGC,NSPOOLGC
         WRITE(16,3401) NOUTGC
 3401    FORMAT(////,1X,'GLOBAL NODAL CONCENTRATION INFORMATION OUTPUT:',
     &        //,5X,'NOUTGC = ',I2)

C.....CHECK INPUT PARAMETER NOUTGC

         IF(ABS(NOUTGC).GT.2) THEN
            IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,3402)
            WRITE(16,3402)
 3402       FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &           //,1X,'YOUR SELECTION OF THE UNIT 15 INPUT PARAMETER',
     &           ' NOUTGC',
     &           /,1X,'IS NOT AN ALLOWABLE VALUE.  CHECK YOUR INPUT!!')
            IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9973)
            WRITE(16,9973)
            STOP
         ENDIF

C.....IF GLOBAL CONCENTRATION OUTPUT WILL NOT BE GENERATED

         IF(NOUTGC.EQ.0) THEN
            WRITE(16,3403)
 3403       FORMAT(///,5X,'NO GLOBAL CONCENTRATION OUTPUT WILL BE ',
     &           'SPOOLED')
         ENDIF

C.....IF GLOBAL CONCENTRATION OUTPUT WILL BE GENERATED

         IF(NOUTGC.NE.0) THEN

C.......COMPUTE NTCYSGC, NTCYFGC, WHICH = TOUTSGC AND TOUTFGC IN TIMESTEPS

            NTCYSGC=INT((TOUTSGC-STATIM)*(86400.D0/DTDP) + 0.5d0)
            NTCYFGC=INT((TOUTFGC-STATIM)*(86400.D0/DTDP) + 0.5d0)
            IF(NTCYFGC.GT.NT) NTCYFGC=NT

C.......CALCULATE NDSETSC = THE # OF DATA SETS TO BE SPOOLED TO UNIT 73

            IF(NSPOOLGC.EQ.0) NDSETSC=0
            IF(NSPOOLGC.NE.0) NDSETSC=INT((NTCYFGC-NTCYSGC)/NSPOOLGC)

C.......WRITE NOUTGC,TOUTSGC,TOUTFGC,NTCYSGC,NTCYFGC,NSPOOLGC TO UNIT 16

            WRITE(16,3404) TOUTSGC,NTCYSGC,TOUTFGC,NTCYFGC,NSPOOLGC
 3404       FORMAT(/,5X,'DATA RECORDS WILL START AFTER TOUTSGC =',F8.3,
     &           ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &           I9,' TIME STEPS INTO THE SIMULATION',
     &           //,5X,'DATA RECORDS WILL STOP AFTER TOUTFGC =',F8.3,
     &           ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &           I9,' TIME STEPS INTO THE SIMULATION',
     &           //,5X,'INFORMATION WILL BE SPOOLED TO UNIT 73 EVERY ',
     &           'NSPOOLGC =',I8,' TIME STEPS')
            IF(ABS(NOUTGC).EQ.1) WRITE(16,3405)
 3405       FORMAT(/,5X,'UNIT 73 FORMAT WILL BE ASCII')
            IF(ABS(NOUTGC).EQ.2) WRITE(16,3406)
 3406       FORMAT(/,5X,'UNIT 73 FORMAT WILL BE BINARY')
         ENDIF

      ENDIF

C...  
C...  IF NWS<>0   INPUT INFORMATION ABOUT GLOBAL WIND DATA OUTPUT
C...  
      IF(NWS.NE.0) THEN

C......READ IN NOUTGW,TOUTSGW,TOUTFGW,NSPOOLGW : IF NOUTGW<>0, GLOBAL WIND
C......OUTPUT IS SPOOLED TO UNIT 74 EVERY NSPOOLGW TIME STEPS BETWEEN
C......TIMES TOUTSGW AND TOUTFGW; IF ABS(NOUTGW)=2, OUTPUT WILL BE BINARY

         READ(15,*) NOUTGW,TOUTSGW,TOUTFGW,NSPOOLGW
         WRITE(16,3451) NOUTGW
 3451    FORMAT(////,1X,'GLOBAL WIND STRESS INFORMATION OUTPUT : ',
     &        //,5X,'NOUTGW = ',I2)

C......CHECK INPUT PARAMETER NOUTGW

         IF(ABS(NOUTGW).GT.2) THEN
            IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,3452)
            WRITE(16,3452)
 3452       FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &           //,1X,'YOUR SELECTION OF THE UNIT 15 INPUT PARAMETER',
     &           ' NOUTGW',
     &           /,1X,'IS NOT AN ALLOWABLE VALUE.  CHECK YOUR INPUT!!')
            IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,3453)
            WRITE(16,3453)
            STOP
         ENDIF

C......IF GLOBAL WIND STRESS OUTPUT WILL NOT BE GENERATED

         IF(NOUTGW.EQ.0) THEN
            WRITE(16,3453)
 3453       FORMAT(///,5X,'NO GLOBAL WIND STRESS OUTPUT WILL BE SPOOLED')
         ENDIF

C........IF GLOBAL WIND STRESS OUTPUT WILL BE GENERATED

         IF(NOUTGW.NE.0) THEN

C........COMPUTE NTCYSGW, NTCYFGW, WHICH = TOUTSGW AND TOUTFGW IN TIMESTEPS

            NTCYSGW=INT((TOUTSGW-STATIM)*(86400.D0/DTDP) + 0.5d0)
            NTCYFGW=INT((TOUTFGW-STATIM)*(86400.D0/DTDP) + 0.5d0)
            IF(NTCYFGW.GT.NT) NTCYFGW=NT

C........CALCULATE NDSETSW = THE # OF DATA SETS TO BE SPOOLED TO UNIT 74

            IF(NSPOOLGW.EQ.0) NDSETSW=0
            IF(NSPOOLGW.NE.0) NDSETSW=INT((NTCYFGW-NTCYSGW)/NSPOOLGW)

C........WRITE NOUTGW,TOUTSGW,TOUTFGW,NTCYSGW,NTCYFGW,NSPOOLGW TO UNIT 16

            WRITE(16,3454) TOUTSGW,NTCYSGW,TOUTFGW,NTCYFGW,NSPOOLGW
 3454       FORMAT(/,5X,'DATA RECORDS WILL START AFTER TOUTSGW =',F8.3,
     &           ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &           I9,' TIME STEPS INTO THE SIMULATION',
     &           //,5X,'DATA RECORDS WILL STOP AFTER TOUTFGW =',F8.3,
     &           ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',
     &           I9,' TIME STEPS INTO THE SIMULATION',
     &           //,5X,'INFORMATION WILL BE SPOOLED TO UNIT 74 EVERY ',
     &           'NSPOOLGW =',I8,' TIME STEPS')
            IF(ABS(NOUTGW).EQ.1) WRITE(16,3455)
 3455       FORMAT(/,5X,'UNIT 74 FORMAT WILL BE ASCII')
            IF(ABS(NOUTGW).EQ.2) WRITE(16,3456)
 3456       FORMAT(/,5X,'UNIT 74 FORMAT WILL BE BINARY')
         ENDIF

      ENDIF

C...  
C...  READ AND CHECK INFORMATION ABOUT HARMONIC ANALYSIS OF MODEL RESULTS
C...  
      READ(15,*) NFREQ 
      WRITE(16,99392) NFREQ  
99392 FORMAT(////,1X,'HARMONIC ANALYSIS INFORMATION OUTPUT : ',
     &     //,5X,'HARMONIC ANALYSIS PERFORMED FOR ',I4,' CONSTITUENTS',/)
      MNHARF = NFREQ

      IF (NFREQ.EQ.0) MNHARF = 1

C     Allocate harmonic analysis arrays

      IF (NFREQ.GT.0) THEN
         CALL ALLOC_HA()
         CALL ALLOC_MAIN14()
      ENDIF

      IF(NFREQ.LT.0) THEN
         IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,99391)
         WRITE(16,99391)
99391    FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
     &        //,1X,'YOUR SELECTION OF NHARFR (A UNIT 15 '
     &        ,'INPUT PARAMETER) IS NOT AN ALLOWABLE VALUE',/,1X,
     &        'PLEASE CHECK YOUR INPUT',
     &        //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
         STOP
      ENDIF
      IF(NFREQ.GT.0) WRITE(16,2330)
 2330 FORMAT(/,7X,'FREQUENCY',4X,'NODAL FACTOR',6X,'EQU.ARG(DEG)',
     +     1X,'CONSTITUENT',/)
      DO 1201 I=1,NFREQ  
         READ(15,'(A10)') NAMEFR(I)
         READ(15,*) HAFREQ(I),HAFF(I),HAFACE(I)
         WRITE(16,2331) HAFREQ(I),HAFF(I),HAFACE(I),NAMEFR(I)
 2331    FORMAT(4X,F15.12,2X,F10.7,5X,F10.3,7X,A10)
 1201 CONTINUE

C.... READ IN INTERVAL INFORMATION FOR HARMONIC ANALYSIS
C.... COMPUTE THAS AND THAF IN TERMS OF THE NUMBER OF TIME STEPS

      READ(15,*) THAS,THAF,NHAINC,FMV
      ITHAS=INT((THAS-STATIM)*(86400.D0/DTDP) + 0.5d0)
      THAS=ITHAS*DTDP/86400.D0 + STATIM
      ITHAF=INT((THAF-STATIM)*(86400.D0/DTDP) + 0.5d0)
      THAF=ITHAF*DTDP/86400.D0 + STATIM
      ITMV = ITHAF - (ITHAF-ITHAS)*FMV
      IF(NFREQ.GT.0) THEN
         WRITE(16,34634) THAS,ITHAS,THAF,ITHAF,NHAINC
34634    FORMAT(/,5X,'HARMONIC ANALYSIS WILL START AFTER THAS =',F8.3,
     &        ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',I9,
     &        ' TIME STEPS INTO THE SIMULATION',
     &        //,5X,'HARMONIC ANALYSIS WILL STOP AFTER THAF =',F8.3,
     &        ' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',I9,
     &        ' TIME STEPS INTO THE SIMULATION'
     &        ,//,5X,'INFORMATION WILL BE ANALYZED EVERY ',
     &        'NHAINC =',I8,' TIME STEPS.')
         WRITE(16,34639) FMV*100.,ITMV
34639    FORMAT(/,5X,'MEANS AND VARIANCES WILL BE COMPUTED FOR THE ',
     &        'FINAL ',F10.5,' %',/9X,'OF THE HARMONIC ANALYSIS ',
     &        'PERIOD OR AFTER ',I9,' TIME STEPS INTO THE ',
     &        'SIMULATION.',/9X,' RESULTS ARE WRITTEN TO UNIT 55.')

      ELSE
         WRITE(16,34645)
34645    FORMAT(///,5X,'NO HARMONIC ANALYSIS WILL BE DONE')
      ENDIF

      IF((FMV.GT.0.).AND.(NFREQ.GT.0).AND.(C2DDI)) CHARMV = .TRUE.

C.... READ IN AND WRITE OUT INFORMATION ON WHERE HARMONIC ANALYSIS WILL BE DONE

      READ(15,*) NHASE,NHASV,NHAGE,NHAGV
      IF((NHASE.LT.0).OR.(NHASE.GT.1)) THEN
         IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,99661)
         WRITE(16,99661)
99661    FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',
     &        'INPUT ERROR  !!!!!!!!!',//
     &        ,1X,'YOUR SELECTION OF NHASE (A UNIT 15 '
     &        ,'INPUT PARAMETER) IS NOT AN ALLOWABLE VALUE',/,1X,
     &        'PLEASE CHECK YOUR INPUT')
         IF(NFOVER.EQ.1) THEN
            IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,99671)
            WRITE(16,99671)
99671       FORMAT(/,1X,'PROGRAM WILL OVERRIDE SPECIFIED INPUT',
     &           ' AND SET NHASE EQUAL TO 0 ',
     &           //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)
            NHASE=0
         ELSE
            IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9973)
            WRITE(16,9973)
            STOP
         ENDIF
      ENDIF
      IF(NHASE.EQ.1) THEN
         WRITE(16,34641)
34641    FORMAT(///,5X,'STATION ELEVATION HARMONIC ANAL WILL BE ',
     &        'WRITTEN TO UNIT 51')
      ENDIF
      IF((NHASV.LT.0).OR.(NHASV.GT.1)) THEN
         IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,99662)
         WRITE(16,99662)
99662    FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',
     &        'INPUT ERROR  !!!!!!!!!',//
     &        ,1X,'YOUR SELECTION OF NHASV (A UNIT 15 '
     &        ,'INPUT PARAMETER) IS NOT AN ALLOWABLE VALUE',/,1X,
     &        'PLEASE CHECK YOUR INPUT')
         IF(NFOVER.EQ.1) THEN
            IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,99672)
            WRITE(16,99672)
99672       FORMAT(/,1X,'PROGRAM WILL OVERRIDE SPECIFIED INPUT',
     &           ' AND SET NHASV EQUAL TO 0 ',
     &           //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)
            NHASV=0
         ELSE
            IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9973)
            WRITE(16,9973)
            STOP
         ENDIF
      ENDIF
      IF(NHASV.EQ.1) THEN
         WRITE(16,34642)
34642    FORMAT(///,5X,'STATION VELOCITY HARMONIC ANAL WILL BE ',
     &        'WRITTEN TO UNIT 52')
      ENDIF
      IF((NHAGE.LT.0).OR.(NHAGE.GT.1)) THEN
         IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,99663)
         WRITE(16,99663)
99663    FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',
     &        'INPUT ERROR  !!!!!!!!!',//
     &        ,1X,'YOUR SELECTION OF NHAGE (A UNIT 15 '
     &        ,'INPUT PARAMETER) IS NOT AN ALLOWABLE VALUE',/,1X,
     &        'PLEASE CHECK YOUR INPUT')
         IF(NFOVER.EQ.1) THEN
            IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,99673)
            WRITE(16,99673)
99673       FORMAT(/,1X,'PROGRAM WILL OVERRIDE SPECIFIED INPUT',
     &           ' AND SET NHAGE EQUAL TO 0 ',
     &           //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)
            NHAGE=0
         ELSE
            IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9973)
            WRITE(16,9973)
            STOP
         ENDIF
      ENDIF
      IF(NHAGE.EQ.1) THEN
         WRITE(16,34643)
34643    FORMAT(///,5X,'GLOBAL ELEVATION HARMONIC ANAL WILL BE ',
     &        'WRITTEN TO UNIT 53')
      ENDIF
      IF((NHAGV.LT.0).OR.(NHAGV.GT.1)) THEN
         IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,99664)
         WRITE(16,99664)
99664    FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',
     &        'INPUT ERROR  !!!!!!!!!',//
     &        ,1X,'YOUR SELECTION OF NHAGV (A UNIT 15 '
     &        ,'INPUT PARAMETER) IS NOT AN ALLOWABLE VALUE',/,1X,
     &        'PLEASE CHECK YOUR INPUT')
         IF(NFOVER.EQ.1) THEN
            IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,99674)
            WRITE(16,99674)
99674       FORMAT(/,1X,'PROGRAM WILL OVERRIDE SPECIFIED INPUT',
     &           ' AND SET NHAGV EQUAL TO 0 ',
     &           //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)
            NHAGV=0
         ELSE
            IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9973)
            WRITE(16,9973)
            STOP
         ENDIF
      ENDIF
      IF(NHAGV.EQ.1) THEN
         WRITE(16,34644)
34644    FORMAT(///,5X,'GLOBAL VELOCITY HARMONIC ANAL WILL BE ',
     &        'WRITTEN TO UNIT 54')
      ENDIF

C.... ESTABLISH INDICATOR OF WHETHER ANY HARMONIC ANALYSIS WILL BE DONE

      IHARIND=NFREQ*(NHASE+NHASV+NHAGE+NHAGV)
      IF(IHARIND.GT.0) IHARIND=1

C...  
C...  INPUT INFORMATION ABOUT HOT START OUTPUT
C...  
      READ(15,*) NHSTAR,NHSINC
      WRITE(16,99655)
99655 FORMAT(////,1X,'HOT START OUTPUT INFORMATION OUTPUT : ')
      IF((NHSTAR.LT.0).OR.(NHSTAR.GT.1)) THEN
         IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,99665)
         WRITE(16,99665)
99665    FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',
     &        'INPUT ERROR  !!!!!!!!!',//
     &        ,1X,'YOUR SELECTION OF NHSTAR (A UNIT 15 '
     &        ,'INPUT PARAMETER) IS NOT AN ALLOWABLE VALUE',/,1X,
     &        'PLEASE CHECK YOUR INPUT')
         IF(NFOVER.EQ.1) THEN
            IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,99675)
            WRITE(16,99675)
99675       FORMAT(/,1X,'PROGRAM WILL OVERRIDE SPECIFIED INPUT',
     &           ' AND SET NHSTAR EQUAL TO 0 ',
     &           //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)
            NHSTAR=0
         ELSE
            IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9973)
            WRITE(16,9973)
            STOP
         ENDIF
      ENDIF
      IF(NHSTAR.EQ.1) THEN
         WRITE(16,34636) NHSINC
34636    FORMAT(/,5X,'HOT START OUTPUT WILL BE WRITTEN TO UNIT',
     &        ' 67 OR 68 EVERY ',I5,' TIME STEPS')
      ELSE
         WRITE(16,34646)
34646    FORMAT(///,5X,'NO HOT START OUTPUT WILL BE GENERATED')
      ENDIF
      IF((IHOT.EQ.0).OR.(IHOT.EQ.68)) IHSFIL=67
      IF(IHOT.EQ.67) IHSFIL=68
C...  
C...  INPUT INFORMATION ABOUT SOLVER
C...  

C...  THIS SECTION TO LUMP THE GWCE MATRIX
C     vjp 11/30/99 made lumping a compile time option

#ifdef LUMP
      CLUMP = .TRUE.
      ILUMP=1
#else
      CLUMP = .FALSE.
      ILUMP=0
#endif
      
      READ(15,*) ITITER,ISLDIA,CONVCR,ITMAX

      WRITE(16,99656)
99656 FORMAT(//,1X,'SOLVER INFORMATION OUTPUT : ')

C     - allocate arrays dimensioned by MNEI
      call alloc_main11()

C...  LINES TO USE THE ITERATIVE MATRIX SOLVER

      IF((ISLDIA.LT.0).OR.(ISLDIA.GT.5)) THEN
         IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9920)
         WRITE(16,9920)
 9920    FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL INPUT ERROR',
     &        ' !!!!!!!!!',//,1X,'ISLDIA (A UNIT 15 INPUT PARAMETER) ',
     &        'MUST BE 0-5',/,1X,'PLEASE CHECK YOUR INPUT')
         IF(NFOVER.EQ.1) THEN
            IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9921)
            WRITE(16,9921)
 9921       FORMAT(/,1X,'PROGRAM WILL OVERRIDE SPECIFIED INPUT',
     &           ' AND SET ISLDIA EQUAL TO 0 ',
     &           //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)
            ISLDIA=0
         ELSE
            IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,9973)
            WRITE(16,9973)
            STOP
         ENDIF
      ENDIF

#ifdef  CMPI
C     sb-PDG1 deleted
C     READ(15,*) MNPROC
C--   
#else
      MNPROC = 1
#endif

C.....INITIALIZE AVERAGING FOR INTERNAL BARRIER WATER LEVELS
C......BARAVGWT=0.000 -> NO AVERAGING PERFORMED
cjj   wm001 changed one line                
      BARAVGWT=0.000D0
      IBSTART=0
      DO I=1,NVEL
         RBARWL1AVG(I)=0.D0
         RBARWL2AVG(I)=0.D0
      END DO
C.....INITIALIZE NIBNODECODE(I)
      DO I=1,NP
         NIBNODECODE(I)=0
      END DO
C...  
C...  COMPUTE THE NEIGHBOR TABLE
C...  
      IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,1196)
      WRITE(16,1196)
 1196 FORMAT(/,1X,'THE NEIGHBOR TABLE IS BEING COMPUTED ',/)
C     
      CALL NEIGHB(NE,NP,NM,NNEIGH,NEITAB,NEIMIN,NEIMAX,X,Y,NSCREEN,NNEIGH_ELEM,NEIGH_ELEM)
C     
      IF(NSCREEN.EQ.1.AND.MYPROC.EQ.0) 
     &     WRITE(6,1195) NEIMIN,NEIMAX,NEIMAX
      WRITE(16,1195) NEIMIN,NEIMAX,NEIMAX
 1195 FORMAT(1X,'THE NEIGHBOR TABLE IS COMPLETED ',
     &     /,5X,'THE MINIMUM NUMBER OF NEIGHBORS FOR ANY NODE = ',I3,
     &     /,5X,'1+THE MAXIMUM NUMBER OF NEIGHBORS FOR ANY NODE = ',I3,
     &     /,5X,'THE PARAMETER MNEI CAN BE SET AS SMALL AS ',I3,/)


C.....Allocate arrays dealing with wind forcing

      CALL ALLOC_MAIN12()
      
C.....Allocate arrays for wave modified bottom friction (EJK)

      IF ((FRW.EQ.1).OR.(NRS.EQ.1)) THEN
         CALL ALLOC_MAIN15()
      ENDIF
      
C     sb-
C.....Count maximum of the number of the elements associated with a node

      DO I = 1,MNP
         NNDEL(I) = 0
      ENDDO
      DO IK=1,MNE
         NNDEL(NM(IK,1)) = NNDEL(NM(IK,1)) + 1
         NNDEL(NM(IK,2)) = NNDEL(NM(IK,2)) + 1
         NNDEL(NM(IK,3)) = NNDEL(NM(IK,3)) + 1
      ENDDO
      MNNDEL = 0
      DO IK=1,MNP
         IF(NNDEL(IK).GT.MNNDEL) MNNDEL = NNDEL(IK)
      ENDDO

C.....Allocate space for Arrays dimensioned by MNNDEL
      CALL ALLOC_MAIN16()
      
      DO I = 1,MNP
         NNDEL(I) = 0
      ENDDO

C.....Make node-to-elements table


      DO IK=1,MNE
         NM1 = NM(IK,1)
         NM2 = NM(IK,2)
         NM3 = NM(IK,3)
         
         NNDEL(NM1) = NNDEL(NM1) + 1
         NDEL(NM1,NNDEL(NM1)) = IK

         NNDEL(NM2) = NNDEL(NM2) + 1
         NDEL(NM2,NNDEL(NM2)) = IK

         NNDEL(NM3) = NNDEL(NM3) + 1
         NDEL(NM3,NNDEL(NM3)) = IK
      ENDDO
C--   

C.....Write table of parameter sizes (vjp 11/28/99)

      WRITE(16,4010) MNPROC,MNE,MNP,MNEI,MNOPE,MNETA,MNBOU,MNVEL,
     &     MNTIF,MNBFR,MNSTAE,MNSTAV,MNSTAC,MNSTAM,NWLAT,NWLON,MNHARF,MNFFR
      IF (NWS.EQ.0) WRITE(16,4011)
      IF (NWS.EQ.1) WRITE(16,4012)
      IF (ABS(NWS).EQ.2) WRITE(16,4013)
      IF (NWS.EQ.3) WRITE(16,4014)
      IF (ABS(NWS).EQ.4) WRITE(16,4015)
      IF (ABS(NWS).EQ.5) WRITE(16,4115)
      IF (NWS.EQ.10) WRITE(16,4016)
      IF (NWS.EQ.11) WRITE(16,4017)
      IF ((NFREQ.EQ.0).OR.(FMV.EQ.0.)) WRITE(16,4021)
      IF ((NFREQ.GE.1).AND.(FMV.NE.0.)) WRITE(16,4022)
      IF (ILUMP.EQ.0) WRITE(16,4031)
      IF (ILUMP.EQ.1) WRITE(16,4032)
      IF (IM.EQ.0) WRITE(16,4101)
      IF (IM.EQ.10) WRITE(16,4109)
      IF (IM.EQ.1) WRITE(16,4102)
      IF (IM.EQ.2) WRITE(16,4103)
      WRITE(16,4105)
      WRITE(16,4108)



 4010 FORMAT(' *****************************************************',/,
     &     ' *   Based on information extracted from the ADCIRC  *',/,
     &     ' *   UNIT 14 and 15 (grid and horiz run info) files  *',/,
     &     ' *   the following paramter values will be set:      *',/,
     &     ' *                                                   *',/,
     &     ' *       MNPROC = ',I5,'                              *',/,
     &     ' *       MNE = ',I8,1X,'     MNP  = ',I8,1X,'        *',/,
     &     ' *       MNEI = ',I7,2X,'                            *',/,
     &     ' *       MNOPE = ',I6,3X,'   MNETA = ',I6,3X,'       *',/,
     &     ' *       MNBOU = ',I6,3X,'   MNVEL = ',I6,3X,'       *',/,
     &     ' *       MNTIF = ',I6,3X,'   MNBFR = ',I6,3X,'       *',/,
     &     ' *       MNSTAE = ',I5,4X,'  MNSTAV = ',I5,4X,'      *',/,
     &     ' *       MNSTAC = ',I5,4X,'  MNSTAM = ',I5,4X,'      *',/,
     &     ' *       MNWLAT = ',I5,4X,'  MNWLON = ',I5,4X,'      *',/,
     &     ' *       MNHARF = ',I5,4X,'  MNFFR = ',I6,3X,'       *',/,
     &     ' *                                                   *')
 4011 FORMAT(' *   Also, NO wind forcing will be used,             *')
 4012 FORMAT(' *   Also, Standard wind stress and pres will be used,*')
 4013 FORMAT(' *   Also, Semi-standard wind forcing will be used,  *')
 4014 FORMAT(' *   Also, Fleet numeric wind forcing will be used,  *')
 4015 FORMAT(' *   Also, PBL/JAG wind forcing will be used,        *')
 4115 FORMAT(' *   Also, Standard wind vel and pres will be used,  *')
 4016 FORMAT(' *   Also, AVN wind & pressure forcing will be used, *')
 4017 FORMAT(' *   Also, ETA wind & pressure forcing will be used, *')
 4021 FORMAT(' *   means and variance calculation will NOT be made,*')
 4022 FORMAT(' *   means and variance calculation will be made,    *')
 4031 FORMAT(' *   the GWCE matrix will be left in consistent form *')
 4032 FORMAT(' *   the GWCE matrix will be LUMPED                  *')
 4101 FORMAT(' *   the model will be set up for a 2DDI run,        *')
 4109 FORMAT(' *   the model will be set up for a 2DDI run + transp*')
 4102 FORMAT(' *   the model will be set up for a 3D-VS run,       *')
 4103 FORMAT(' *   the model will be set up for a 3D-DSS run,      *')
 4105 FORMAT(' *   and an iterative solver will be used            *')
 4108 FORMAT(' *****************************************************',/)

C.....Close fort.15, fort.14, and fort.dg files

      CLOSE(14)
      CLOSE(15)

      RETURN 
      END SUBROUTINE
