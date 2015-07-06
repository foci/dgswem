!***********************************************************************
!     
!     SUBROUTINE READ_INPUT()
!     
!     Surprisingly enough this subroutine reads in the input
!     
!     Written by (too) many people
!     
!***********************************************************************

      SUBROUTINE READ_INPUT(s,dg_here,global_here,nodalattr_here)
      use sizes
!.....Use appropriate modules

      USE GLOBAL
#ifdef HARM
      USE HARM
#endif
      USE WIND
      USE DG
      USE NodalAttributes
      USE fort_dg, ONLY: read_keyword_fort_dg,read_fixed_fort_dg
#ifdef CMPI
      USE MESSENGER_ELEM
#endif
#ifdef SWAN
!asey 101118: Enable hot-start file generation by SWAN.
      USE Couple2Swan, ONLY: SwanHotStartUnit
#endif

      IMPLICIT NONE

      type (sizes_type) :: s
      type (dg_type) :: dg_here
      type (global_type) :: global_here
      type (nodalattr_type) :: nodalattr_here
      
!.....Declare local variables

      INTEGER NIBP, IBN1, IK, NDISC, NBBN, NVEL2, II, i,j,jj,k,IPRBI_here,ICK_here
      CHARACTER(256) LINE,LINE2
#ifndef HARM
      integer nfreq_dummy
#endif
      
!     sb-PDG1 added
#ifdef CMPI
      INTEGER IDUM80
      CHARACTER CDUM80
#endif

#ifdef HPX
      INTEGER IDUM80
      CHARACTER CDUM80
#endif
!--   
!     ek...Zero out all the variables in the Nodal Attributes Module
!     ek...Added from version 46

      CALL InitNAModule(nodalattr_here)

!     sb-PDG1 added
#ifdef CMPI
!     
!     When running in parallel check to make sure that the number of
!     processors (MNPROC - obtained from the job control script via a call
!     to the messenger module) is the same as that used in ADCPREP and written
!     in the fort.80 file
!     
      OPEN(80,FILE='fort.80')
      READ(80,'(A)') CDUM80     !Skip global_here%RUNDES
      READ(80,'(A)') CDUM80     !Skip global_here%RUNID
      READ(80,'(A)') CDUM80     !Skip global_here%AGRID
      READ(80,*) IDUM80         !Skip NELG & NNODG
      READ(80,*) IDUM80         !Read in NPROC
      CLOSE(80)
      IF(IDUM80.NE.s%MNPROC) THEN
         IF(s%MYPROC.EQ.0) THEN
            WRITE(*,'(A)') '*** ERROR IN PARALLEL SETUP!'
            WRITE(*,'(2A,I4,A)') '*** Number of CPUS for submitted job ',&
           '(NCPU = ',s%MNPROC,') is not equal to the'
            WRITE(*,'(2A,I4,A)') '*** number of CPUS specified during',&
           ' ADCPREP (see fort.80: NCPU = ',IDUM80,').'
            WRITE(*,'(A)') '*** dgswem will now quit!'
         ENDIF
         CALL MESSAGE_FINI(s)
         STOP
      ENDIF
#endif
!--   
#ifdef HPX
      OPEN(80,FILE='fort.80')
      READ(80,'(A)') CDUM80     !Skip global_here%RUNDES
      READ(80,'(A)') CDUM80     !Skip global_here%RUNID
      READ(80,'(A)') CDUM80     !Skip global_here%AGRID
      READ(80,*) IDUM80         !Skip NELG & NNODG
      READ(80,*) IDUM80         !Read in NPROC
      CLOSE(80)
      IF(IDUM80.NE.s%MNPROC) THEN
         IF(s%MYPROC.EQ.0) THEN
            WRITE(*,'(A)') '*** ERROR IN PARALLEL SETUP!'
            WRITE(*,'(2A,I4,A)') '*** Number of CPUS for submitted job ',&
           '(NCPU = ',s%MNPROC,') is not equal to the'
            WRITE(*,'(2A,I4,A)') '*** number of CPUS specified during',&
           ' ADCPREP (see fort.80: NCPU = ',IDUM80,').'
            WRITE(*,'(A)') '*** dgswem will now quit!'
         ENDIF
         STOP
      ENDIF
#endif
      global_here%ScreenUnit = 6

!.....Initialize all runtime option logicals to false

      s%C2DDI  = .FALSE.
      s%C3D    = .FALSE.
      global_here%vertexslope = .False.
      s%C3DDSS = .FALSE.
      s%C3DVS  = .FALSE.
      s%CLUMP  = .FALSE.
      s%CTIP   = .FALSE.
      s%CHARMV = .FALSE.

!.....Open statement for unit 14, 15, and 25 (fort.dg) input files

      OPEN(14,FILE=s%DIRNAME//'/'//'fort.14')
      OPEN(15,FILE=s%DIRNAME//'/'//'fort.15')
      OPEN(17,FILE=s%DIRNAME//'/'//'fort.17')
      OPEN(25,FILE=s%DIRNAME//'/'//'fort.dg')            

!.....Open statement for unit 16 output file
      
      OPEN(16,FILE=s%DIRNAME//'/'//'fort.16')

!.....General purpose format statements

 1112 FORMAT(/,1X,79('_'))
 9972 FORMAT(////,1X,'!!!!!!!!!! INPUT ERROR !!!!!!!!!',/)
 9973 FORMAT(/,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
 9974 FORMAT(/,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!!!',//)



!.....Print out header for output including version number and copyright

      WRITE(16,1112)
      WRITE(16,1112)
      WRITE(16,1114)
      WRITE(16,1112)
      IF (s%MYPROC.EQ.0) THEN
         WRITE(6,1112)
         WRITE(6,1114)
         WRITE(6,1112)
      ENDIF

 1114 FORMAT(//,19X,'dgswem.11.13 ',&
     //,5X,'A disconintuous Galerkin shallow water equation model' &
     /,10X,'for coastal seas and estuarine research',&
     ///,7X,'-  Based off of the ADCIRC model created by',&
     //,10X,'R.A. Luettich, Jr., University of North Carolina',&
     /,10X,'J.J. Westerink, University of Notre Dame',&
     //,7X,'-  The ADCIRC source code has been copyrighted by',&
     /,10X,'R.A. Luettich, Jr. & J.J. Westerink, 1994-2001',&
     /,10X,'No part of the adcirc base code may be reproduced',&
     /,10X,'or redistributed without the written permission of',&
     /,10X,'the above authors.',&
     ///,5X,'The DG version of the code was written largely ab initio,',   &
     /,5X,'though some data structures from the CG version were used.',&
     //,7X,'Developed by:',&
     /,10X,'Ethan J. Kubatko, The Ohio State University (EJK)',&
     /,10X,'Clint N. Dawson, UT ICES (cnd)',&
     /,10X,'Shintaro Bunya, The University of Tokyo (sb)',&
     /,10X,'Craig Michoski, UT ICES (cem)',&
     /,10X,'Christopher Mirabito, MIT',&
     /,10X,'Damrongsak Wirasaet, University of Notre Dame',&
     /,10X,'Casey Dietrich, North Carolina State University',//)


!.....Write out header information describing how code has been set up

      WRITE(16,1210)
 1210 FORMAT(//,1X,'THE SOURCE CODE HAS BEEN CONFIGURED ',&
     'BY THE PREPROCESSOR AS FOLLOWS:',/)

#ifdef C3DDSS
      WRITE(16,*) '      - 3D DSS MODEL OPTION'
#endif 

#ifdef C3DVS
      WRITE(16,*) '      - 3D VS MODEL OPTION'
#else
      WRITE(16,*) '      - 2D DEPTH INTEGRATED MODEL OPTION'
#endif 

#ifdef CMACHSUN
      WRITE(16,*) '      - CODE SETUP TO RUN ON SUN 4 OR SPARC ',&
     'COMPUTERS'
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

      WRITE(16,*) '      - NONVECTORIZABLE PARTS OF CODE OPTIMIZED FOR',&
     ' MEMORY'
      WRITE(16,*) '      - CODE WILL USE JCG ITERATIVE GWCE SOLVER'
      WRITE(16,1112)

      
      
!.....Read in the fort.dg file

      !srb - check fort.dg format (for backwards compatibility)
      READ(25,*) LINE
      CLOSE(25)
      LINE2 = ADJUSTL(LINE)
      
      IF (LINE2(1:1) .EQ. "1") THEN
        CALL read_fixed_fort_dg(s,dg_here,global_here)   ! first line of old fort.dg is a 1 for the global_here%dgswe option
      ELSE
        CALL read_keyword_fort_dg(s,dg_here,global_here) ! otherwise assume keyword format
      ENDIF     
      
      global_here%RHOWAT0 = 1000.D0
                  
      
!.....Input from unit 15 and output to unit 16 rundescription and run ID

      READ(15,'(A32)') global_here%RUNDES
      READ(15,'(A24)') global_here%RUNID
      WRITE(16,1) global_here%RUNDES
 1    FORMAT(//,1X,'RUN DESCRIPTION : ',A32)
      WRITE(16,209) global_here%RUNID
 209  FORMAT(/,1X,'RUN IDENTIFICATION : ',A24)

!.....Read and process global_here%NFOVER - nonfatal error override otion

      READ(15,*) global_here%NFOVER
      WRITE(16,1112)
      WRITE(16,1250)
 1250 FORMAT(//,1X,'GENERAL RUN INFORMATION',/)
      IF(global_here%NFOVER.EQ.1) THEN
         WRITE(16,1951) global_here%NFOVER
 1951    FORMAT(5X,'global_here%NFOVER = ',I2,&
        /,9X,'IF NON-FATAL ERRORS ARE DETECTED, THEY WILL BE ',&
        'CORRECTED AND EXECUTION CONTINUED')
      ELSE
         WRITE(16,1952) global_here%NFOVER
 1952    FORMAT(/,5X,'global_here%NFOVER = ',I3,&
        /,9X,'NON-FATAL ERRORS WILL STOP EXECUTION ',/)
      ENDIF

!.....Read and process global_here%NABOUT - abbreviated unit 16 output option

      READ(15,*) global_here%NABOUT
      IF (global_here%NABOUT.EQ.1) THEN
         WRITE(16,3501) global_here%NABOUT
 3501    FORMAT(5X,'global_here%NABOUT = ',I2,&
        /,9X,'ABREVIATED OUTPUT WILL BE PROVIDED TO UNIT 16',&
        /,9X,'UNIT 14, 21, 22 INPUT DATA WILL NOT BE ECHO PRINTED',/)
      ELSE
         WRITE(16,3502) global_here%NABOUT
 3502    FORMAT(/,5X,'global_here%NABOUT = ',I3,&
        /,9X,'DETAILED OUTPUT WILL BE PROVIDED TO UNIT 16',&
        /,9X,'UNIT 14, 15, 21, 22 INPUT DATA WILL BE ECHO PRINTED',/)
      ENDIF

!.....Read and process global_here%NSCREEN - screen ouput option

      READ(15,*) global_here%NSCREEN
      global_here%NSCREEN_INC = global_here%NSCREEN
      IF (global_here%NSCREEN.NE.0) global_here%NSCREEN = 1
      IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) THEN
         WRITE(16,3561) global_here%NSCREEN
 3561    FORMAT(5X,'global_here%NSCREEN = ',I2,&
        /,9X,'SCREEN OUTPUT WILL BE PROVIDED TO UNIT 6',/)
      ELSE
         WRITE(16,3562) global_here%NSCREEN
 3562    FORMAT(/,5X,'global_here%NSCREEN = ',I3,&
        /,9X,'SCREEN OUTPUT WILL NOT BE PROVIDED TO UNIT 6',/)
      ENDIF
      
!.....Read and process global_here%IHOT - hot start option

      READ(15,*) global_here%IHOT
      IF ((global_here%IHOT.NE.0).AND.(global_here%IHOT.NE.67).AND.(global_here%IHOT.NE.68)) THEN
         IF ((global_here%NSCREEN.EQ.1).AND.(s%MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'global_here%IHOT =',global_here%IHOT
            WRITE(6,9732)
            WRITE(6,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'global_here%IHOT =',global_here%IHOT
         WRITE(16,9732)
         WRITE(16,9973)       
 9732    FORMAT(/,1X,'Your selection of global_here%IHOT (a UNIT 15 input ',&
        'parameter) is not an allowable value')
         STOP
      ENDIF
      IF (global_here%IHOT.NE.0) THEN
         WRITE(16,9733) global_here%IHOT
 9733    FORMAT(/,5X,'dgswem will be hot started using information ',&
        'on UNIT ',I2)
      ELSE
         WRITE(16,9734)
 9734    FORMAT(/,5X,'dgswem will be cold started')
      ENDIF
#ifdef SWAN
!asey 100205: Enable hot-start file generation by SWAN.
      SwanHotStartUnit = global_here%IHOT
#endif     


!.....Read and process global_here%ICS - cartesian/spherical coordinate option

      READ(15,*) global_here%ICS
      IF ((global_here%ICS.NE.1).AND.(global_here%ICS.NE.2)) THEN
         IF ((global_here%NSCREEN.EQ.1).AND.(s%MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'global_here%ICS =',global_here%ICS
            WRITE(6,9735)
            WRITE(6,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'global_here%ICS =',global_here%ICS
         WRITE(16,9735)
         WRITE(16,9973)
 9735    FORMAT(/,1X,'Your selection of global_here%ICS (a UNIT 15 input ',&
        'parameter) is not an allowable value')
         STOP
      ENDIF
      IF (global_here%ICS.EQ.1) THEN
         WRITE(16,9736) global_here%ICS
 9736    FORMAT(/,5X,'global_here%ICS = ',I2,&
        /,9X,'Governing equations are in Cartesian coordinates')
      ELSE
         WRITE(16,9737) global_here%ICS
 9737    FORMAT(/,5X,'global_here%ICS = ',I2,&
        /,9X,'Governing equations are in Spherical coordinates',&
        /,9X,'mapped using a CPP projection')
      ENDIF

!.....Read and process global_here%IM - 2D/3D model option

      READ(15,*) global_here%IM
      IF (global_here%IM.EQ.0) THEN
         s%C2DDI = .TRUE.
      ELSEIF (global_here%IM.EQ.1) THEN
         s%C3D  = .TRUE.
         s%C3DVS  = .TRUE.
      ELSEIF (global_here%IM.EQ.2) THEN
         s%C3D  = .TRUE.
         s%C3DDSS = .TRUE.
      ELSEIF (global_here%IM.EQ.10) THEN
         s%C2DDI = .TRUE.
      ELSE
         IF ((global_here%NSCREEN.EQ.1).AND.(s%MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'global_here%IM =',global_here%IM
            WRITE(6,9721)
            WRITE(6,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'global_here%IM =',global_here%IM
         WRITE(16,9721)
         WRITE(16,9973)
 9721    FORMAT(/,1X,'Your selection of global_here%IM (a UNIT 15 input ',&
        'parameter) is not an allowable value')
         STOP
      ENDIF

!.....Read and process nodalattr_here%NOLIBF - nonlinear bottom friction option

      READ(15,*) nodalattr_here%NOLIBF
      IF ((nodalattr_here%NOLIBF.LT.0).OR.(nodalattr_here%NOLIBF.GT.2)) THEN
         IF ((global_here%NSCREEN.EQ.1).AND.(s%MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'nodalattr_here%NOLIBF =',nodalattr_here%NOLIBF
            WRITE(6,9722)
            WRITE(6,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'nodalattr_here%NOLIBF =',nodalattr_here%NOLIBF
         WRITE(16,9722)
         WRITE(16,9973)
 9722    FORMAT(/,1X,'Your selection of nodalattr_here%NOLIBF (a UNIT 15 input ',&
        'parameter) is not an allowable value')
         STOP
      ENDIF
      WRITE(16,9845) nodalattr_here%NOLIBF
 9845 FORMAT(/,5X,'nodalattr_here%NOLIBF = ',I3)
      IF (nodalattr_here%NOLIBF.EQ.0) WRITE(16,2050)
 2050 FORMAT(9X,'THE MODEL WILL USE LINEAR BOTTOM FRICTION')
      IF (nodalattr_here%NOLIBF.EQ.1) WRITE(16,2051)
 2051 FORMAT(9X,'THE MODEL WILL USE STANDARD QUADRATIC BOTTOM FRICTION')
      IF (nodalattr_here%NOLIBF.EQ.2) WRITE(16,2052)
 2052 FORMAT(9X,'THE MODEL WILL USE STANDARD QUADRATIC BOTTOM FRICTION',&
     'IN DEEP WATER ',&
     /,9X,'AND A FRICTION FACTOR THAT INCREASES AS THE DEPTH ',&
     'DECREASES IN SHALLOW WATER')

!.....Read and process global_here%NOLIFA - nonlinear finite amplitude option

      READ(15,*) global_here%NOLIFA
      IF ((global_here%NOLIFA.LT.0).OR.(global_here%NOLIFA.GT.3)) THEN
         IF ((global_here%NSCREEN.EQ.1).AND.(s%MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'global_here%NOLIFA =',global_here%NOLIFA
            WRITE(6,9723)
            WRITE(6,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'global_here%NOLIFA =',global_here%NOLIFA
         WRITE(16,9723)
         WRITE(16,9973)
 9723    FORMAT(/,1X,'Your selection of global_here%NOLIFA (a UNIT 15 input ',&
        'parameter) is not an allowable value')
         STOP
      ENDIF
      WRITE(16,9846) global_here%NOLIFA
 9846 FORMAT(/,5X,'global_here%NOLIFA = ',I3)
      IF (global_here%NOLIFA.EQ.0) WRITE(16,2053)
 2053 FORMAT(9X,'THE MODEL WILL NOT USE FINITE AMPLITUDE TERMS OR ',&
     'WETTING AND DRYING')
      IF (global_here%NOLIFA.EQ.1) WRITE(16,2054)
 2054 FORMAT(9X,'THE MODEL WILL USE FINITE AMPLITUDE TERMS BUT NO ',&
     'WETTING AND DRYING')
      IF (global_here%NOLIFA.EQ.2) WRITE(16,2049)
 2049 FORMAT(9X,'THE MODEL WILL USE FINITE AMPLITUDE TERMS AND ',&
     'WETTING AND DRYING')
      IF (global_here%NOLIFA.EQ.3) WRITE(16,2048)
 2048 FORMAT(9X,'THE MODEL WILL USE FINITE AMPLITUDE TERMS AND ',&
     'WETTING AND DRYING',/,10X,&
     'AND INCLUDES THE ABILITY TO INITIALIZE ',&
     'NODES WITH DEPTHS GREATER THAN global_here%H0 AS DRY')         
      global_here%NSTARTDRY = 0
      IF (global_here%NOLIFA.EQ.3) THEN
         global_here%NOLIFA = 2
         global_here%NSTARTDRY = 1
      ENDIF

!.....Read and process global_here%NOLICA - advective term spatial gradinet

      READ(15,*) global_here%NOLICA
      IF ((global_here%NOLICA.LT.0).OR.(global_here%NOLICA.GT.1)) THEN
         IF ((global_here%NSCREEN.EQ.1).AND.(s%MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'global_here%NOLICA =',global_here%NOLICA
            WRITE(6,9724)
            WRITE(6,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'global_here%NOLICA =',global_here%NOLICA
         WRITE(16,9724)
         WRITE(16,9973)
 9724    FORMAT(/,1X,'Your selection of global_here%NOLICA (a UNIT 15 input ',&
        'parameter) is not an allowable value')
         STOP
      ENDIF
      WRITE(16,9847) global_here%NOLICA
 9847 FORMAT(/,5X,'global_here%NOLICA = ',I3)
      IF (global_here%NOLICA.EQ.0) WRITE(16,2055)
 2055 FORMAT(9X,'THE MODEL WILL NOT USE SPATIAL DERIVATIVE ',&
     'COMPONENTS OF THE ADVECTIVE TERMS')
      IF (global_here%NOLICA.EQ.1) WRITE(16,2056)
 2056 FORMAT(9X,'THE MODEL WILL USE SPATIAL DERIVATIVE ',&
     'COMPONENTS OF THE ADVECTIVE TERMS')

!.....Read and process global_here%NOLICAT - GWCE advective term time derivative

      READ(15,*) global_here%NOLICAT
      IF ((global_here%NOLICAT.LT.0).OR.(global_here%NOLICAT.GT.1)) THEN
         IF ((global_here%NSCREEN.EQ.1).AND.(s%MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'global_here%NOLICAT =',global_here%NOLICAT
            WRITE(6,9725)
            WRITE(6,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'global_here%NOLICAT =',global_here%NOLICAT
         WRITE(16,9725)
         WRITE(16,9973)
 9725    FORMAT(/,1X,'Your selection of global_here%NOLICAT (a UNIT 15 input ',&
        'parameter) is not an allowable value')
         STOP
      ENDIF
      IF ((global_here%NOLIFA.GE.1).AND.(global_here%NOLICAT.EQ.0)) THEN
         IF ((global_here%NSCREEN.EQ.1).AND.(s%MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'global_here%NOLICAT =',global_here%NOLICAT
            WRITE(6,9726)
            IF (global_here%NFOVER.EQ.1) THEN
               WRITE(6,9974)
            ELSE
               WRITE(6,9973)
            ENDIF
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'global_here%NOLICAT =',global_here%NOLICAT
         WRITE(16,9726)
         WRITE(16,9974)
 9726    FORMAT(/,1X,'Your selection of global_here%NOLICAT (a UNIT 15 input ',&
        'parameter) is inconsistent with your ',&
        /,1X,'selection of global_here%NOLIFA and may lead to mass ',&
        'dg_here%balance problems')
         IF (global_here%NFOVER.EQ.1) THEN
            WRITE(6,9974)
         ELSE
            WRITE(6,9973)
            STOP
         ENDIF
      ENDIF
      IF ((global_here%NOLIFA.EQ.0).AND.(global_here%NOLICAT.EQ.1)) THEN
         IF ((global_here%NSCREEN.EQ.1).AND.(s%MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'global_here%NOLICAT =',global_here%NOLICAT
            WRITE(6,9726)
            IF (global_here%NFOVER.EQ.1) THEN
               WRITE(6,9974)
            ELSE
               WRITE(6,9973)
            ENDIF
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'global_here%NOLICAT =',global_here%NOLICAT
         WRITE(16,9726)
         WRITE(16,9974)
         IF (global_here%NFOVER.EQ.1) THEN
            WRITE(6,9974)
         ELSE
            WRITE(6,9973)
            STOP
         ENDIF
      ENDIF
      IF (global_here%NOLICA.NE.global_here%NOLICAT) THEN
         IF ((global_here%NSCREEN.EQ.1).AND.(s%MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'global_here%NOLICAT =',global_here%NOLICAT
            WRITE(6,9727)
            IF (global_here%NFOVER.EQ.1) THEN
               WRITE(6,9974)
            ELSE
               WRITE(6,9973)
            ENDIF
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'global_here%NOLICAT =',global_here%NOLICAT
         WRITE(16,9727)
         WRITE(16,9974)
 9727    FORMAT(/,1X,'Your selection of global_here%NOLICAT (a UNIT 15 input ',&
        'parameter) is inconsistent with your ',&
        /,1X,'selection of global_here%NOLICA and may lead to mass ',&
        'dg_here%balance problems')
         IF (global_here%NFOVER.EQ.1) THEN
            WRITE(6,9974)
         ELSE
            WRITE(6,9973)
            STOP
         ENDIF
      ENDIF
      WRITE(16,9848) global_here%NOLICAT
 9848 FORMAT(/,5X,'global_here%NOLICAT = ',I3)
      IF (global_here%NOLICAT.EQ.0) WRITE(16,2057)
 2057 FORMAT(9X,'THE MODEL WILL NOT USE TIME DERIVATIVE COMPONENTS ',&
     /,9X,'OF THE ADVECTIVE TERMS IN THE GWCE')
      IF (global_here%NOLICAT.EQ.1) WRITE(16,2058)
 2058 FORMAT(9X,'THE MODEL WILL USE TIME DERIVATIVE COMPONENTS ',&
     /,9X,'OF THE ADVECTIVE TERMS IN THE GWCE')

!.....Read and process nodalattr_here%NWP - spatially varying bottom friction

      READ(15,*) nodalattr_here%NWP
      CALL ReadNodalAttr(s, nodalattr_here, global_here%NSCREEN, global_here%ScreenUnit, s%MYPROC, global_here%NABOUT) ! Ek added call to nodalatt

!.....Read and process global_here%NCOR - spatially varying Coriolis parameter

      READ(15,*) global_here%NCOR
      IF ((global_here%NCOR.NE.0).AND.(global_here%NCOR.NE.1)) THEN
         IF ((global_here%NSCREEN.EQ.1).AND.(s%MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'global_here%NCOR =',global_here%NCOR
            WRITE(6,9729)
            WRITE(6,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'global_here%NCOR =',global_here%NCOR
         WRITE(16,9729)
         WRITE(16,9973)
 9729    FORMAT(/,1X,'Your selection of global_here%NCOR (a UNIT 15 input ',&
        'parameter) is not an allowable value')
         STOP
      ENDIF
      IF ((global_here%ICS.EQ.1).AND.(global_here%NCOR.EQ.1)) THEN
         IF ((global_here%NSCREEN.EQ.1).AND.(s%MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'global_here%NCOR =',global_here%NCOR
            WRITE(6,9730)
            WRITE(6,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'global_here%NCOR =',global_here%NCOR
         WRITE(16,9730)
         WRITE(16,9973)
 9730    FORMAT(/,1X,'Your selection of global_here%NCOR (a UNIT 15 input ',&
        'parameter) is inconsistent with your ',&
        /,1X,'selection of coordinate systems.  Spatially ',&
        'variable Coriolis should be used only with ',&
        /,1X,'Spherical coordinates')
         STOP
      ENDIF
      IF (global_here%NCOR.EQ.0) THEN
         WRITE(16,233) global_here%NCOR
 233     FORMAT(/,5X,'global_here%NCOR = ',I2,&
        /,9X,'A CONSTANT VALUE OF THE CORIOLIS PARAMETER WILL BE ',&
        /,9X,'USED THROUGHOUT THE DOMAIN')
      ELSE
         WRITE(16,234) global_here%NCOR
 234     FORMAT(/,5X,'global_here%NCOR = ',I2,&
        /,9X,'SPATIALLY VARYING CORIOLIS VALUES WILL BE COMPUTED ',&
        'FROM INPUT LATITUDES')
      ENDIF

!.....Read and process global_here%NTIP - tidal potential forcing

      READ(15,*) global_here%NTIP
      IF ((global_here%NTIP.LT.0).OR.(global_here%NTIP.GT.2)) THEN
         IF ((global_here%NSCREEN.EQ.1).AND.(s%MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'global_here%NTIP =',global_here%NTIP
            WRITE(6,9710)
            WRITE(6,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'global_here%NTIP =',global_here%NTIP
         WRITE(16,9710)
         WRITE(16,9973)
 9710    FORMAT(/,1X,'Your selection of global_here%NTIP (a UNIT 15 input ',&
        'parameter) is not an allowable value')
         STOP
      ENDIF

      IF ((global_here%ICS.EQ.1).AND.(global_here%NTIP.GE.1)) THEN
         IF ((global_here%NSCREEN.EQ.1).AND.(s%MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'global_here%NTIP =',global_here%NTIP
            WRITE(6,9711)
            WRITE(6,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'global_here%NTIP =',global_here%NTIP
         WRITE(16,9711)
         WRITE(16,9973)
 9711    FORMAT(/,1X,'Your selection of global_here%NTIP (a UNIT 15 input ',&
        'parameter) is inconsistent with your ',&
        /,1X,'selection of coordinate systems.  Tidal',&
        'potential forcing should be used only with ',&
        /,1X,'Spherical coordinates')
         STOP
      ENDIF
      IF (global_here%NTIP.NE.0) s%CTIP = .TRUE.
      IF (global_here%NTIP.EQ.0) THEN
         WRITE(16,235) global_here%NTIP
 235    FORMAT(/,5X,'global_here%NTIP = ',I2,&
    /,9X,'TIDAL POTENTIAL FORCING IS NOT USED IN THE COMPUTATION')
      ENDIF
      IF (global_here%NTIP.GE.1) THEN
        WRITE(16,236) global_here%NTIP
 236     FORMAT(/,5X,'global_here%NTIP = ',I2,&
        /,9X,'TIDAL POTENTIAL FORCING IS USED IN THE COMPUTATION ',&
        'BASED ON INPUT LONGITUDES/LATITUDES')
      ENDIF
      IF (global_here%NTIP.EQ.2) THEN
         WRITE(16,239)
 239     FORMAT(9X,'SELF ATTRACTION/LOAD TIDE FORCING IS ALSO USED ',&
        'IN THE COMPUTATION')
      ENDIF

!.....Read and process global_here%NWS - wind and pressure forcing & wave rad stress

      READ(15,*) global_here%NWS
      IF ( (global_here%NWS.NE.0)  .AND.(global_here%NWS.NE.1)       .AND.(ABS(global_here%NWS).NE.2)  .AND.&
     (global_here%NWS.NE.3)  .AND.(ABS(global_here%NWS).NE.4)  .AND.(ABS(global_here%NWS).NE.5)  .AND.&
     (global_here%NWS.NE.6)  .AND.(global_here%NWS.NE.10)      .AND.(global_here%NWS.NE.11)      .AND.&
     (ABS(global_here%NWS).NE.12 ).AND.&
     (global_here%NWS.NE.100).AND.(global_here%NWS.NE.101)     .AND.(ABS(global_here%NWS).NE.102).AND.&
     (global_here%NWS.NE.103).AND.(ABS(global_here%NWS).NE.104).AND.(ABS(global_here%NWS).NE.105).AND.&
     (global_here%NWS.NE.106).AND.(global_here%NWS.NE.110)     .AND.(global_here%NWS.NE.111)     .AND.&
     (global_here%NWS.NE.200).AND.(global_here%NWS.NE.201)     .AND.(ABS(global_here%NWS).NE.202).AND.&
     (global_here%NWS.NE.203).AND.(ABS(global_here%NWS).NE.204).AND.(ABS(global_here%NWS).NE.205).AND.&
     (global_here%NWS.NE.8).AND.&
     (global_here%NWS.NE.206).AND.(global_here%NWS.NE.210)     .AND.(global_here%NWS.NE.211) .AND.&
!asey 101118: Added the following cases for coupling to unstructured SWAN.
   (ABS(global_here%NWS).NE.300).AND.(ABS(global_here%NWS).NE.303).AND.(ABS(global_here%NWS).NE.304).AND.&
   (ABS(global_here%NWS).NE.305).AND.(ABS(global_here%NWS).NE.306).AND.(ABS(global_here%NWS).NE.308).AND.&
   (ABS(global_here%NWS).NE.309).AND.(ABS(global_here%NWS).NE.310).AND.(ABS(global_here%NWS).NE.311).AND.&
   (ABS(global_here%NWS).NE.312)) THEN
         IF ((global_here%NSCREEN.EQ.1).AND.(s%MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'global_here%NWS =',global_here%NWS
            WRITE(6,9712)
            WRITE(6,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'global_here%NWS =',global_here%NWS
         WRITE(16,9712)
         WRITE(16,9973)
 9712    FORMAT(/,1X,'Your selection of global_here%NWS (a UNIT 15 input ',&
        'parameter) is not an allowable value')
         STOP
      ENDIF

!.....Set wave radiation stress flag and adjust global_here%NWS accordingly
!asey 101118: Had to make some changes in this section.
      global_here%NRS = 0
      global_here%FRW = 0
      IF (ABS(global_here%NWS/100).EQ.1) THEN
        global_here%NRS=1
        global_here%NWS = (ABS(global_here%NWS) - 100)*(global_here%NWS/ABS(global_here%NWS))
      ENDIF
      IF (ABS(global_here%NWS/100).EQ.2) THEN
        global_here%NRS = 1
        global_here%FRW = 1
        global_here%NWS = (ABS(global_here%NWS) - 200)*(global_here%NWS/ABS(global_here%NWS))
      ENDIF
#ifdef SWAN
!asey 101118: Added the option for coupling directly to SWAN.
      IF(ABS(global_here%NWS/100).EQ.3) THEN
        global_here%NRS = 3
        global_here%NWS = (ABS(global_here%NWS) - 300)*(global_here%NWS/ABS(global_here%NWS))
      ENDIF
#endif

      IF (global_here%NWS.EQ.0) THEN
         WRITE(16,237) global_here%NWS
 237     FORMAT(/,5X,'global_here%NWS = ',I2,&
    /,9X,'WIND STRESS OR SURFACE PRESSURE ARE NOT USED TO FORCE',&
         'THE COMPUTATION')
      ENDIF
      IF (global_here%NWS.EQ.1) THEN
        WRITE(16,238) global_here%NWS
 238    FORMAT(/,5X,'global_here%NWS = ',I2,&
    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',&
    /,9X,' THE COMPUTATION',&
    /,9X,'VALUES ARE READ AT ADCIRC GRID NODES FROM UNIT 22',&
    /,9X,' EVERY MODEL TIME STEP')
      ENDIF
      IF (global_here%NWS.EQ.2) THEN
        WRITE(16,2381) global_here%NWS
 2381   FORMAT(/,5X,'global_here%NWS = ',I2,&
    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',&
    /,9X,' THE COMPUTATION',&
    /,9X,'VALUES ARE READ AT ADCIRC GRID NODES FROM UNIT 22',&
    /,9X,'THE UNIT 22 FILE BEGINS AT TIME=global_here%STATIM.',&
    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',&
    /,9X,'WITH THE MODEL TIME STEP.')
      ENDIF

      IF (global_here%NWS.EQ.-2) THEN
        WRITE(16,2380) global_here%NWS
 2380   FORMAT(/,5X,'global_here%NWS = ',I2,&
    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',&
    /,9X,' THE COMPUTATION',&
    /,9X,'VALUES ARE READ AT ADCIRC GRID NODES FROM UNIT 22',&
    /,9X,'THE UNIT 22 FILE BEGINS AT THE TIME OF THE HOT START.',&
    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',&
       /,9X,'WITH THE MODEL TIME STEP.')
      ENDIF

      IF (global_here%NWS.EQ.3) THEN
         WRITE(16,2382) global_here%NWS
 2382    FORMAT(/,5X,'global_here%NWS = ',I2,&
    /,9X,'WIND STRESS ONLY IS USED TO FORCE THE COMPUTATION.',&
    /,9X,'WIND SPEEDS AND DIRECTIONS ARE READ FROM A FLEET ',&
    /,9X,'NUMERIC FORMAT FILE AT UNIT 22 AND INTERPOLATED TO',&
    /,9X,'THE ADCIRC GRID. ',&
    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',&
    /,9X,'WITH THE MODEL TIME STEP.',&
    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',&
        'DRAG LAW.')
      ENDIF
      IF (global_here%NWS.EQ.4) THEN
         WRITE(16,2383) global_here%NWS
 2383    FORMAT(/,5X,'global_here%NWS = ',I2,&
    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',&
    /,9X,' THE COMPUTATION',&
    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ AT SELECTED',&
    /,9X,'ADCIRC GRID NODES FROM A PBL FILE AT UNIT 22.',&
    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',&
    /,9X,'WITH THE MODEL TIME STEP.',&
    /,9X,'THE UNIT 22 FILE BEGINS AT TIME=global_here%STATIM.',&
    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',&
        'DRAG LAW.')
      ENDIF
      IF (global_here%NWS.EQ.-4) THEN
         WRITE(16,2388) global_here%NWS
 2388    FORMAT(/,5X,'global_here%NWS = ',I2,&
    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',&
    /,9X,' THE COMPUTATION',&
    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ AT SELECTED',&
    /,9X,'ADCIRC GRID NODES FROM A PBL FILE AT UNIT 22.',&
    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',&
    /,9X,'WITH THE MODEL TIME STEP.',&
    /,9X,'THE UNIT 22 FILE BEGINS AT THE TIME OF THE HOT START.',&
    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',&
        'DRAG LAW.')
      ENDIF
      IF (global_here%NWS.EQ.5) THEN
         WRITE(16,2384) global_here%NWS
 2384    FORMAT(/,5X,'global_here%NWS = ',I2,&
    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',&
    /,9X,' THE COMPUTATION',&
    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ AT ADCIRC ',&
    /,9X,'GRID NODES FROM UNIT 22',&
    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',&
    /,9X,'WITH THE MODEL TIME STEP.',&
    /,9X,'THE UNIT 22 FILE BEGINS AT TIME=global_here%STATIM.',&
     /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',&
     'DRAG LAW.')
      ENDIF
      IF (global_here%NWS.EQ.-5) THEN
         WRITE(16,2389) global_here%NWS
 2389    FORMAT(/,5X,'global_here%NWS = ',I2,&
    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',&
    /,9X,' THE COMPUTATION',&
    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ AT ADCIRC ',&
    /,9X,'GRID NODES FROM UNIT 22',&
    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',&
    /,9X,'WITH THE MODEL TIME STEP.',&
    /,9X,'THE UNIT 22 FILE BEGINS AT THE TIME OF THE HOT START.',&
    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',&
        'DRAG LAW.')
      ENDIF
      IF (global_here%NWS.EQ.6) THEN
         WRITE(16,2385) global_here%NWS
 2385    FORMAT(/,5X,'global_here%NWS = ',I2,&
    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',&
    /,9X,' THE COMPUTATION',&
    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ FROM A ',&
    /,9X,'REGULARLY SPACED GRID FROM UNIT 22',&
    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',&
    /,9X,'WITH THE MODEL TIME STEP AND IN SPACE TO BRING THE ',&
    /,9X,'MET DATA FROM A REGULAR GRID TO THE ADCIRC GRID.'&
    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',&
        'DRAG LAW.')
      ENDIF
      IF (global_here%NWS.EQ.10) THEN
         WRITE(16,2386) global_here%NWS
 2386    FORMAT(/,5X,'global_here%NWS = ',I2,&
    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',&
    /,9X,' THE COMPUTATION',&
    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ EVERY N',&
    /,9X,' HOURS FROM A DIFFERENT FILE AT UNITS 200, 200+N,',&
    ' 200+2N, ETC.',&
    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',&
    /,9X,'WITH THE MODEL TIME STEP AND IN SPACE TO BRING THE ',&
    /,9X,'MET DATA FROM A GAUSSIAN GRID TO THE ADCRIC GRID.',&
    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',&
        'DRAG LAW.')
      ENDIF
      IF (global_here%NWS.EQ.11) THEN
         WRITE(16,2387) global_here%NWS
 2387    FORMAT(/,5X,'global_here%NWS = ',I2,&
    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',&
    /,9X,' THE COMPUTATION',&
    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ EVERY 3 ',&
    /,9X,'HOURS FROM ETA-29 FILES AT UNITS 200, 201, 202, ETC.',&
    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',&
    /,9X,'WITH THE MODEL TIME STEP AND IN SPACE TO BRING THE ',&
    /,9X,'WIND DATA FROM THE 29 KM E GRID TO THE ADCRIC GRID.',&
    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',&
        'DRAG LAW.')
      ENDIF
      IF (global_here%NRS.EQ.0) THEN
         WRITE(16,2390) global_here%NRS
 2390    FORMAT(/,5X,'global_here%NRS = ',I2,&
        /,9X,'WAVE RADIATION STRESS IS NOT USED TO FORCE THE ',&
        'COMPUTATION')
      ENDIF
      
!.....ek added for global_here%NWS=-12,12

      IF(global_here%NWS.EQ.12) THEN
         WRITE(16,12384) global_here%NWS
12384    FORMAT(/,5X,'global_here%NWS = ',I2,&
    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',&
    /,9X,' THE COMPUTATION',&
    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ FROM RAW ',&
    /,9X,'OWI DATA FILES (UNIT 221-224).',&
    /,9X,'META DATA IS READ FROM UNIT 220.',&
    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',&
    /,9X,'WITH THE MODEL TIME STEP.',&
    /,9X,'THE UNIT 22 FILE BEGINS AT TIME=global_here%STATIM.',&
    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',&
        'DRAG LAW.')
      ENDIF
      IF(global_here%NWS.EQ.-12) THEN
         WRITE(16,12389) global_here%NWS
12389    FORMAT(/,5X,'global_here%NWS = ',I3,&
    /,9X,'WIND STRESS AND SURFACE PRESSURE ARE USED TO FORCE',&
    /,9X,' THE COMPUTATION',&
    /,9X,'WIND VELOCITY AND PRESSURE VALUES ARE READ FROM RAW ',&
    /,9X,'OWI DATA FILES (UNIT 221-224).',&
    /,9X,'META DATA IS READ FROM UNIT 220.',&
    /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE WIND DATA ',&
    /,9X,'WITH THE MODEL TIME STEP.',&
    /,9X,'THE UNIT 22 FILE BEGINS AT THE TIME OF THE HOT START.',&
    /,9X,'WIND SPEEDS ARE CONVERTED TO STRESS USING THE GARRET ',&
        'DRAG LAW.')
      ENDIF
      IF (global_here%NRS.EQ.1) THEN
         WRITE(16,2391) global_here%NRS
 2391    FORMAT(/,5X,'global_here%NRS = ',I2,&
    /,9X,'WAVE RADIATION STRESS IS USED TO FORCE THE COMPUTATION',&
    /,9X,'STRESSES ARE READ AT SELECTED ADCIRC GRID NODES FROM A',&
    /,9X,'PBL TYPE FILE AT UNIT 23.  INTERPOLATION IN TIME IS ',&
    /,9X,'DONE TO SYNC THE STRESS DATA WITH THE MODEL TIME STEP.',&
    /,9X,'FOR A COLD START, THE UNIT 23 FILE BEGINS AT THE TIME ',&
    /,9X,'OF THE COLD START.  FOR A HOT START, THE UNIT 23 FILE ',&
    /,9X,'BEGINS AT THE TIME OF THE HOT START.')
      ENDIF
#ifdef SWAN
!asey 101118: Added the following lines.
      IF(global_here%NRS.EQ.3) THEN
         WRITE(16,2393) global_here%NRS
 2393    FORMAT(/,5X,'global_here%NRS = ',I2,&
    /,9X,'WAVES WILL BE COUPLED TO SWAN!')
      ENDIF
#endif

!.....Read and process global_here%NRAMP - whether a global_here%ramp function will be used

      READ(15,*) global_here%NRAMP
      IF ((global_here%NRAMP.NE.0).AND.(global_here%NRAMP.GT.7)) THEN
         IF ((global_here%NSCREEN.EQ.1).AND.(s%MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'global_here%NRAMP =',global_here%NRAMP
            WRITE(6,9713)
            WRITE(6,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'global_here%NRAMP =',global_here%NRAMP
         WRITE(16,9713)
         WRITE(16,9973)
 9713    FORMAT(/,1X,'Your selection of global_here%NRAMP (a UNIT 15 input ',&
        'parameter) is not an allowable value')
         STOP
      ENDIF
      IF (global_here%NRAMP.EQ.0) THEN
         WRITE(16,240) global_here%NRAMP
 240     FORMAT(/,5X,'global_here%NRAMP = ',I2,&
   /,9X,'NO global_here%RAMP FUNCTION IS USED IN THE COMPUTATION')
      ELSE
         WRITE(16,241) global_here%NRAMP
 241     FORMAT(/,5X,'global_here%NRAMP = ',I2,&
   /,9X,'A HYPERBOLIC TANGENT global_here%RAMP IS APPLIED TO THE FORCING ',&
        'FUNCTIONS')
      ENDIF

!.....Read and process global_here%G - gravity
!...  
      READ(15,*) global_here%G
      global_here%G2ROOT = SQRT(global_here%G/2.0d0)
      IF ((global_here%ICS.EQ.2).AND.(abs(global_here%G-9.81).gt.0.01)) THEN
         IF ((global_here%NSCREEN.EQ.1).AND.(s%MYPROC.EQ.0)) THEN
            WRITE(6,9972)
            WRITE(6,*) 'global_here%G =',global_here%G
            WRITE(6,9714)
            WRITE(6,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'global_here%G =',global_here%G
         WRITE(16,9714)
         WRITE(16,9973)
 9714    FORMAT(/,1X,'Your specification of the gravitational ',&
        'constant, global_here%G, (a UNIT 15 input) is not ',&
        /,1X,'consistant with the use of spherical coordinates.',&
        '  global_here%G must be in units of m/s^2')
         STOP
      ENDIF
      WRITE(16,5) global_here%G
    5 FORMAT(///,5X,'GRAVITATIONAL CONSTANT global_here%G =',F10.5,/)

!.....Read and process nodalattr_here%TAU0 - weighting coefficient in the GWCE

      READ(15,*) nodalattr_here%TAU0
      IF (nodalattr_here%TAU0.LT.0) THEN
         WRITE(16,6)
 6       FORMAT(/,5X,'WEIGHTING COEFFICIENT FOR THE GENERALIZED',&
        ' WAVE CONTINUITY EQUATION :',&
        /,5x,'THIS VALUE WILL BE  SELECTED BASED ON NODAL DEPTH',&
        ' ONCE DEPTHS HAVE BEEN PROCESSED',    &
        /,5X,' DEPTH > 10       -> nodalattr_here%TAU0 = 0.005  ',    &
        /,5X,' 10 >/ DEPTH        -> nodalattr_here%TAU0 = 0.020 ',  &
        /,5X,' nodalattr_here%STARTDRY VALUE = -77777  -> nodalattr_here%TAU0 = 0.020 ',  &
        /,5X,' nodalattr_here%STARTDRY VALUE = -88888  -> nodalattr_here%TAU0 = 0.020 ',/)  
      ELSE
         WRITE(16,7) nodalattr_here%TAU0
 7       FORMAT(/,5X,'WEIGHTING COEFFICIENT FOR THE GENERALIZED',&
        ' WAVE CONTINUITY EQUATION :',&
        /,5X, 'nodalattr_here%TAU0 = ',E15.8,2X,'1/sec',/)
      ENDIF

!.....Input from unit 15 and output to unit 16 time integration info

      WRITE(16,1112)
      WRITE(16,245)
 245  FORMAT(//,1X,'TIME INTEGRATION INFORMATION',//)

!.....Read and process global_here%DT - model time step

      READ(15,*) global_here%DTDP
      global_here%DT = global_here%DTDP
      WRITE(16,9) global_here%DTDP
    9 FORMAT(5X,'TIME STEP =',F12.6,5X,'SECONDS',/)

!.....Read and process global_here%STATIM - simulation starting time

      READ(15,*) global_here%STATIM
      WRITE(16,1113) global_here%STATIM
 1113 FORMAT(5X,'STARTING TIME FOR SIMULATION = ',F14.6,' DAYS',/)

!.....Read anf process global_here%REFTIM - harmonic reference time

      READ(15,*) global_here%REFTIM
      WRITE(16,1115) global_here%REFTIM
 1115 FORMAT(5X,'Harmonic REFERENCE TIME = ',F14.6,' DAYS',/)

!.....Read in and process additional timing information for wind.
!asey 101118: Changed global_here%NRS.EQ.1 to global_here%NRS.GE.1 throughout this section.
      IF ((global_here%NWS.EQ.0).AND.(global_here%NRS.GE.1)) READ(15,*) global_here%RSTIMINC
      IF ((global_here%NWS.EQ.1).AND.(global_here%NRS.GE.1)) READ(15,*) global_here%RSTIMINC
      IF (ABS(global_here%NWS).EQ.2) THEN
         IF (global_here%NRS.EQ.0) READ(15,*) global_here%WTIMINC
         IF (global_here%NRS.GE.1) READ(15,*) global_here%WTIMINC,global_here%RSTIMINC
      ENDIF
      IF (global_here%NWS.EQ.3) THEN
         READ(15,*) global_here%IREFYR,global_here%IREFMO,global_here%IREFDAY,global_here%IREFHR,global_here%IREFMIN,global_here%REFSEC
         WRITE(16,1116) global_here%IREFMO,global_here%IREFDAY,global_here%IREFYR,global_here%IREFHR,global_here%IREFMIN,global_here%REFSEC
 1116    FORMAT(5X,'WIND REFERENCE TIME FOR SIMULATION = ',&
        I2,'/',I2,'/',I2,'  ',I2,':',I2,':',f7.4,/)
         CALL TIMECONV(global_here%IREFYR,global_here%IREFMO,global_here%IREFDAY,global_here%IREFHR,global_here%IREFMIN,global_here%REFSEC,&
        global_here%WREFTIM, S%MYPROC, global_here%NScreen, global_here%ScreenUnit )
         IF (global_here%NRS.EQ.0) READ(15,*) global_here%NWLAT,global_here%NWLON,global_here%WLATMAX,global_here%WLONMIN,global_here%WLATINC,&
        global_here%WLONINC,global_here%WTIMINC
         IF (global_here%NRS.GE.1) READ(15,*) global_here%NWLAT,global_here%NWLON,global_here%WLATMAX,global_here%WLONMIN,global_here%WLATINC,&
        global_here%WLONINC,global_here%WTIMINC,global_here%RSTIMINC
      ENDIF
      IF (ABS(global_here%NWS).EQ.4) THEN
         IF (global_here%NRS.EQ.0) READ(15,*) global_here%WTIMINC
         IF (global_here%NRS.GE.1) READ(15,*) global_here%WTIMINC,global_here%RSTIMINC
      ENDIF

      IF (ABS(global_here%NWS).EQ.5) THEN
         IF (global_here%NRS.EQ.0) READ(15,*) global_here%WTIMINC
         IF (global_here%NRS.GE.1) READ(15,*) global_here%WTIMINC,global_here%RSTIMINC
      ENDIF
      IF (global_here%NWS.EQ.6) THEN
         IF (global_here%NRS.EQ.0) READ(15,*) global_here%NWLAT,global_here%NWLON,global_here%WLATMAX,global_here%WLONMIN,global_here%WLATINC,&
        global_here%WLONINC,global_here%WTIMINC
         IF (global_here%NRS.GE.1) READ(15,*) global_here%NWLAT,global_here%NWLON,global_here%WLATMAX,global_here%WLONMIN,global_here%WLATINC,&
        global_here%WLONINC,global_here%WTIMINC,global_here%RSTIMINC
      ENDIF
      IF(ABS(global_here%NWS).EQ.8) THEN
         IF(global_here%NRS.EQ.0) THEN
            READ(15,*) global_here%IREFYR,global_here%IREFMO,global_here%IREFDAY,global_here%IREFHR,StormNumber,BLAdj
         ELSEIF (global_here%NRS.GE.1) THEN
            READ(15,*) global_here%IREFYR,global_here%IREFMO,global_here%IREFDAY,global_here%IREFHR,StormNumber,BLAdj,&
           global_here%RSTIMINC
            WRITE(16,6111) global_here%IREFMO,global_here%IREFDAY,global_here%IREFYR,global_here%IREFHR
 6111       FORMAT(5X,'WIND REFERENCE TIME FOR SIMULATION = ',&
           I2,'/',I2,'/',I2,'  ',I2,'H',/)
         ENDIF
         CALL TIMECONV(global_here%IREFYR,global_here%IREFMO,global_here%IREFDAY,global_here%IREFHR,0,0.0d0,&
        WindRefTime, S%MYPROC, global_here%NScreen, global_here%ScreenUnit)
      ENDIF

      IF (global_here%NWS.EQ.10) THEN
         global_here%NWLAT=190
         global_here%NWLON=384
         IF (global_here%NRS.EQ.0) READ(15,*) global_here%WTIMINC
         IF (global_here%NRS.GE.1) READ(15,*) global_here%WTIMINC,global_here%RSTIMINC
      ENDIF
      IF (global_here%NWS.EQ.11) THEN
         global_here%NWLAT=271
         global_here%NWLON=181
         global_here%WTIMINC=10800.
         IF (global_here%NRS.GE.1) READ(15,*) global_here%RSTIMINC
      ENDIF
      IF (global_here%NWS.EQ.11) THEN
         global_here%NWLAT=271
         global_here%NWLON=181
         global_here%WTIMINC=10800.
         IF (global_here%NRS.GE.1) READ(15,*) global_here%RSTIMINC
      ENDIF
!.....ek added global_here%NWS=12 (OWI format) from version 46

      IF(ABS(global_here%NWS).EQ.12) THEN
         IF(global_here%NRS.EQ.0) READ(15,*) global_here%WTIMINC
         IF(global_here%NRS.GE.1) READ(15,*) global_here%WTIMINC,global_here%RSTIMINC ! sb46.28sb03
      ENDIF

      IF (global_here%NWS.NE.0) WRITE(16,1117) global_here%WTIMINC
 1117 FORMAT(5X,'WIND TIME INCREMENT (SEC) = ',F10.2,/)
      IF (global_here%NRS.NE.0) WRITE(16,1118) global_here%RSTIMINC
 1118 FORMAT(5X,'RADIATION STRESS TIME INCREMENT (SEC) = ',F10.2,/)

!.....Read and process global_here%RNDAY - simulation duration i days

      READ(15,*) global_here%RNDAY
      WRITE(16,10) global_here%RNDAY
 10   FORMAT(5X,'TOTAL LENGTH OF NUMERICAL SIMULATION =',F12.4,&
     5X,'DAYS',/)

!.....Compute total number of time steps global_here%NT

      global_here%NT = INT(global_here%RNDAY*(86400.D0/global_here%DTDP) + 0.5D0)
      WRITE(16,1920) global_here%NT
 1920 FORMAT(5X,'NUMBER OF TIME STEPS  =',I8,/)

!.....Read and process effective length of hyperbolic tangent global_here%ramp

!...  
!...  READ AND PROCESS EFFECTIVE LENGTH OF THE HYPERBOLIC TANGENT global_here%RAMP(S)
!...  IN DAYS
!...  
!     jgf46.08 Add fine-grained global_here%ramp functions.
!     jgf46.21 Add global_here%FluxSettlingTime for global_here%IBTYPE=52 to accomodate
!     MS river during Katrina, split ramps for flux b.c.s into internal
!     and external.
      global_here%FluxSettlingTime = 0.0d0
      global_here%DRamp = 1.0d0
      SELECT CASE(global_here%NRamp)
!     ---------
      CASE(0,1)                 ! Either no global_here%ramp, or same global_here%ramp for all forcings
!     ---------
         READ(15,*) global_here%DRamp
         global_here%DRampIntFlux = global_here%DRamp
         global_here%DRampExtFlux = global_here%DRamp
         global_here%DRampElev    = global_here%DRamp
         global_here%DRampTip     = global_here%DRamp
         global_here%DRampMete    = global_here%DRamp
         global_here%DRampWRad    = global_here%DRamp
!     -------
      CASE(2)                   ! global_here%Ramp for external flux boundary conditions.
!     -------
         READ(15,*) global_here%DRamp,global_here%DRampExtFlux,global_here%FluxSettlingTime
         global_here%DRampIntFlux = global_here%DRamp
         global_here%DRampElev    = global_here%DRamp
         global_here%DRampTip     = global_here%DRamp
         global_here%DRampMete    = global_here%DRamp
         global_here%DRampWRad    = global_here%DRamp
!     -------
      CASE(3)                   ! global_here%Ramp for internal flux boundary conditions.
!     -------
         READ(15,*) global_here%DRamp,global_here%DRampExtFlux,global_here%FluxSettlingTime,global_here%DRampIntFlux
         global_here%DRampElev = global_here%DRamp
         global_here%DRampTip  = global_here%DRamp
         global_here%DRampMete = global_here%DRamp
         global_here%DRampWRad = global_here%DRamp
!     -------
      CASE(4)                   ! global_here%Ramp for surface elevation specified boundary conditions.
!     -------
         READ(15,*) global_here%DRamp,global_here%DRampExtFlux,global_here%FluxSettlingTime,global_here%DRampIntFlux,&
        global_here%DRampElev
         global_here%DRampTip  = global_here%DRamp
         global_here%DRampMete = global_here%DRamp
         global_here%DRampWRad = global_here%DRamp
!     -------
      CASE(5)                   ! global_here%Ramp for tidal potential
!     -------
         READ(15,*) global_here%DRamp,global_here%DRampExtFlux,global_here%FluxSettlingTime,global_here%DRampIntFlux,&
        global_here%DRampElev,global_here%DRampTip
         global_here%DRampMete = global_here%DRamp
         global_here%DRampWRad = global_here%DRamp
!     -------
      CASE(6)                   ! global_here%Ramp for wind and atmospheric pressure
!     -------
         READ(15,*) global_here%DRamp,global_here%DRampExtFlux,global_here%FluxSettlingTime,global_here%DRampIntFlux,&
        global_here%DRampElev,global_here%DRampTip,global_here%DRampMete
         global_here%DRampWRad = global_here%DRamp
!     -------
      CASE(7)                   ! global_here%Ramp for wave radiation stress
!     -------
         READ(15,*) global_here%DRamp,global_here%DRampExtFlux,global_here%FluxSettlingTime,global_here%DRampIntFlux,&
        global_here%DRampElev,global_here%DRampTip,global_here%DRampMete,global_here%DRampWRad
!     ------------
      CASE DEFAULT              ! fall-through
!     ------------
         IF(global_here%NSCREEN.NE.0.AND.s%MYPROC.EQ.0) THEN
            WRITE(global_here%ScreenUnit,9972)
            WRITE(global_here%ScreenUnit,*) 'global_here%NRAMP =',global_here%NRAMP
            WRITE(global_here%ScreenUnit,9713)
            WRITE(global_here%ScreenUnit,9973)
         ENDIF
         WRITE(16,9972)
         WRITE(16,*) 'global_here%NRAMP =',global_here%NRAMP
         WRITE(16,9713)
         WRITE(16,9973)
#ifdef CMPI
         call MESSAGE_FINI(s)
#endif
         STOP
      END SELECT

!.....Compute the total number of timesteps the ramping function is used
 
      global_here%IRAMPING = INT(global_here%DRAMP*(86400.D0/global_here%DTDP) + 0.5D0)

!.....Read GWCE time weighting factors

      READ(15,*) global_here%A00,global_here%B00,global_here%C00
      WRITE(16,14)
 14   FORMAT(//,5X,'TIME WEIGHTING FACTORS IN THE WAVE EQUATION :'/)
      WRITE(16,15) global_here%A00,global_here%B00,global_here%C00
 15   FORMAT(9X,'AT TIME LEVEL K+1 : ',F8.5,&
     /,9X,'AT TIME LEVEL K   : ',F8.5,&
     /,9X,'AT TIME LEVEL K-1 : ',F8.5,/)

!.....Read minimum depth or wet/dry parameters from unit 15

      IF (global_here%NOLIFA.NE.2) THEN
         READ(15,*) global_here%H0
         WRITE(16,16) global_here%H0
 16      FORMAT(//,5X,'THE BATHYMETRIC DEPTH AT ALL NODES WILL BE ',&
        'INCREASED TO global_here%H0= ',F12.4,' IF NECESSARY'/)
      ENDIF
      IF (global_here%NOLIFA.EQ.2) THEN
         READ(15,*) global_here%H0,global_here%NODEDRYMIN,global_here%NODEWETMIN,global_here%VELMIN
         WRITE(16,17) global_here%H0,global_here%NODEWETMIN,global_here%VELMIN,global_here%NODEDRYMIN
 17      FORMAT(//,5X,'DRYING WILL OCCUR WHEN THE WATER DEPTH < global_here%H0',&
        /,5X,'global_here%H0 = ',E16.8,&
        /,5X,'AND global_here%NODEREP > global_here%NODEWETMIN = ',I6,' TIME STEPS',&
        /,5X,'global_here%NODEREP = NUMBER OF TIME STEPS SINCE A NODE ',&
        'CHANGED STATE (EITHER WETTED OR DRIED)',&
        //,5X,'WETTING WILL OCCUR WHEN THERE IS A FAVORABLE ',&
        'PRESSURE GRADIENT THAT',&
        /,5X,'WOULD DRIVE A STEADY VELOCITY TOWARDS A DRY NODE',&
        /,5X,'THAT IS GREATER THAN global_here%VELMIN = ',F10.5,&
        /,5X,'AND global_here%NODEREP > global_here%NODEDRYMIN = ',I6,' TIME STEPS',/)
      ENDIF

!.....Read grid information from units 14 and 15

      READ(14,'(A24)') global_here%AGRID
      READ(14,*) global_here%NE,global_here%NP
      s%MNP = global_here%NP
      s%MNE = global_here%NE
      READ(15,*) global_here%SLAM0,global_here%SFEA0
      global_here%SLAM0 = global_here%SLAM0*DEG2RAD
      global_here%SFEA0 = global_here%SFEA0*DEG2RAD
      WRITE(16,1112)
      WRITE(16,246)
 246  FORMAT(//,1X,'GRID INFORMATION',//)

!.....Allocate arrays dimensioned by MNP and MNE

      CALL ALLOC_MAIN1(s,global_here)

!.....If global_here%ICS = 1 input nodal coordinates and bathymetry from unit 14
!.....If either global_here%NTIP or global_here%NCOR = 1 compute the inverse CPP projection

      IF (global_here%ICS.EQ.1) THEN
         DO I = 1,global_here%NP
            READ(14,*) global_here%JKI,global_here%X(global_here%JKI),global_here%Y(global_here%JKI),global_here%DP(global_here%JKI)
            IF(global_here%JKI.NE.I) THEN
               IF ((global_here%NSCREEN.EQ.1).AND.(s%MYPROC.EQ.0)) WRITE(6,99801)
               WRITE(16,99801)
99801          FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',&
              'INPUT ERROR  !!!!!!!!!',&
              //,1X,'YOUR NODE NUMBERING IS NOT SEQUENTIAL ',&
              'CHECK YOUR UNIT 14 INPUT FILE CAREFULLY',//)
            ENDIF
!     IF ((global_here%NTIP.GE.1).OR.(global_here%NCOR.EQ.1))
            CALL INVCP(global_here%X(global_here%JKI),global_here%Y(global_here%JKI),global_here%SLAM(global_here%JKI),global_here%SFEA(global_here%JKI),global_here%SLAM0,global_here%SFEA0)
!     ENDIF
         ENDDO
      ENDIF
      
!.....If global_here%ICS = 2 input nodal coordinates and bathymetry from unit 14
!.....and compute CPP projection

      IF (global_here%ICS.EQ.2) THEN
         DO I = 1,global_here%NP
            READ(14,*) global_here%JKI,global_here%SLAM(global_here%JKI),global_here%SFEA(global_here%JKI),global_here%DP(global_here%JKI)
            IF (global_here%JKI.NE.I) THEN
               IF ((global_here%NSCREEN.EQ.1).AND.(s%MYPROC.EQ.0)) WRITE(6,99801)
               WRITE(16,99801)
            ENDIF
            global_here%SLAM(global_here%JKI) = DEG2RAD*global_here%SLAM(global_here%JKI)
            global_here%SFEA(global_here%JKI) = DEG2RAD*global_here%SFEA(global_here%JKI)
            CALL CPP(global_here%X(global_here%JKI),global_here%Y(global_here%JKI),global_here%SLAM(global_here%JKI),global_here%SFEA(global_here%JKI),global_here%SLAM0,global_here%SFEA0)
         ENDDO
      ENDIF

!.....If global_here%ICS = 1 set the global_here%SFAC equal to unity

      IF(global_here%ICS.EQ.1) THEN
         DO I=1,global_here%NP
            global_here%SFAC(I)=1.0d0
         ENDDO
      ENDIF
      
!.....If global_here%ICS = 2 compute global_here%SFAC to adjust equations to CPP coordinates

      IF (global_here%ICS.EQ.2) THEN
         DO I = 1,global_here%NP
            global_here%SFAC(I) = COS(global_here%SFEA0)/COS(global_here%SFEA(I))
         ENDDO
      ENDIF

!.....If wetting and drying will not be used make sure all bathymetric
!.....depths are > or = to global_here%H0

      IF ((global_here%NOLIFA.EQ.0).OR.(global_here%NOLIFA.EQ.1)) THEN
         DO I = 1,global_here%NP
            IF (global_here%DP(I).LT.global_here%H0) global_here%DP(I)=global_here%H0
         ENDDO
      ENDIF

!.....Read the global connectivity table from unit 14, compute element
!.....global_here%areas, check that sufficient accuracy is provided by the code to
!.....handle the input grid, check to make sure correct convention has
!.....been used for inputting the connectivity table

      DO I = 1, global_here%NP
         global_here%NNEIGH(I) = 0
      ENDDO

      DO I = 1,global_here%NE
         READ(14,*) global_here%JKI,global_here%NHY,global_here%NM(global_here%JKI,1),global_here%NM(global_here%JKI,2),global_here%NM(global_here%JKI,3)
         global_here%NNEIGH(global_here%NM(global_here%JKI,1)) = global_here%NNEIGH(global_here%NM(global_here%JKI,1)) + 1
         global_here%NNEIGH(global_here%NM(global_here%JKI,2)) = global_here%NNEIGH(global_here%NM(global_here%JKI,2)) + 1
         global_here%NNEIGH(global_here%NM(global_here%JKI,3)) = global_here%NNEIGH(global_here%NM(global_here%JKI,3)) + 1
         IF (global_here%JKI.NE.I) THEN
            IF ((global_here%NSCREEN.EQ.1).AND.(s%MYPROC.EQ.0)) WRITE(6,99802)
            WRITE(16,99802)
99802       FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',&
           'INPUT ERROR  !!!!!!!!!',&
           //,1X,'YOUR ELEMENT NUMBERING IS NOT SEQUENTIAL ',&
           /,1X,'CHECK YOUR UNIT 14 INPUT FILE CAREFULLY',//)
         ENDIF
         global_here%X1 = global_here%X(global_here%NM(global_here%JKI,1))
         global_here%X2 = global_here%X(global_here%NM(global_here%JKI,2))
         global_here%X3 = global_here%X(global_here%NM(global_here%JKI,3))
         global_here%Y1 = global_here%Y(global_here%NM(global_here%JKI,1))
         global_here%Y2 = global_here%Y(global_here%NM(global_here%JKI,2))
         global_here%Y3 = global_here%Y(global_here%NM(global_here%JKI,3))
         global_here%AVGXY = (ABS(global_here%X1)+ABS(global_here%X2)+ABS(global_here%X3)+ABS(global_here%Y1)+ABS(global_here%Y2)+ABS(global_here%Y3))/6.D0
         global_here%DIF1R = global_here%AVGXY/(((global_here%X2-global_here%X1)**2+(global_here%Y2-global_here%Y1)**2)**0.5d0)
         global_here%DIF2R = global_here%AVGXY/(((global_here%X3-global_here%X2)**2+(global_here%Y3-global_here%Y2)**2)**0.5d0)
         global_here%DIF3R = global_here%AVGXY/(((global_here%X3-global_here%X1)**2+(global_here%Y3-global_here%Y1)**2)**0.5d0)
         global_here%DIF1R = LOG10(global_here%DIF1R)
         global_here%DIF2R = LOG10(global_here%DIF2R)
         global_here%DIF3R = LOG10(global_here%DIF3R)
         IF((global_here%DIF1R.GT.NPREC).OR.(global_here%DIF2R.GT.NPREC).OR.(global_here%DIF3R.GT.NPREC))THEN
            IF ((global_here%NSCREEN.EQ.1).AND.(s%MYPROC.EQ.0)) WRITE(6,9898) global_here%JKI
            WRITE(16,9898) global_here%JKI
 9898       FORMAT(////,1X,'!!!!!!!!!!  WARNING  !!!!!!!!!',&
           //,1X,'IF THE GRID COORDINATES HAVE 32 BITS ',&
           '(APPROX 7 DIGITS) OF PRECISION',&
           /,1X,'A ROBUST MODEL SOLUTION CAN NOT BE GUARANTEED',&
           'AT ELEMENT NO. ',I10,&
           //,1X,'MORE PRECISION MUST BE USED IN THE GRID',//)
         ENDIF
         
!.....NOTE: This is 2 times the actual element area (why?)

         global_here%AREAS(global_here%JKI) = (global_here%X1 - global_here%X3)*(global_here%Y2 - global_here%Y3) + (global_here%X3 - global_here%X2)*(global_here%Y1 - global_here%Y3)
         IF (global_here%AREAS(global_here%JKI).LT.0.0) THEN
            IF((global_here%NSCREEN.EQ.1).AND.(s%MYPROC.EQ.0)) WRITE(6,9899) global_here%JKI
            WRITE(16,9899) global_here%JKI
 9899       FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',&
           //,1X,'THE CONNECTIVITY FOR ELEMENT ',I6,&
           '  HAS BEEN INCORRECTLY SPECIFIED ',&
           /,1X,'CHECK INPUT AND ENSURE THAT COUNTERCLOCKWISE',&
           ' CONVENTION HAS BEEN USED ',&
           //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
            STOP
         ENDIF
      ENDDO

!.....If baroclinic 2D run, read in initial density field

!     IF (global_here%IM.EQ.100) THEN
!     OPEN(11,FILE=DIRNAME//'/'//'fort.11')
!     READ(11,*)
!     READ(11,*)
!     READ(11,*) global_here%NP2
!     IF (global_here%NP2.global_here%NE.global_here%NP) THEN
!     IF ((global_here%NSCREEN.EQ.1).AND.(MYPROC.EQ.0)) WRITE(6,9943)
!     WRITE(16,9943)
!99   43     FORMAT(////,' !!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',
!     &              //,' THE NUMBER OF NODES (global_here%NP2) IN THE BAROCLINIC',
!     &                 ' INITIAL CONDITION FILE (UNIT 11) ',
!     &               /,' MUST EQUAL THE NUMBER OF NODES (global_here%NP) IN ',
!     &                 'THE ADCIRC GRID FILE (UNIT 14)'
!     &              //,' !!!!! EXECUTION WILL NOW BE TERMINATED !!!!!')
!     STOP
!     ENDIF
!     
!     DO I = 1,global_here%NP
!     READ(11,*) global_here%JKI,DASIGT(global_here%JKI),DATEMP(global_here%JKI),DASAL(global_here%JKI)
!     ENDDO
!     CLOSE(11)
!     ENDIF

!.....Process nodalattr_here%startdry info from unit 12 (if global_here%NOLIFA = 3 -> global_here%NSTARTDRY = 1
!.....nodalattr_here%STARTDRY now set in nodal attributes

      IF ((global_here%NSTARTDRY.EQ.1).AND.(nodalattr_here%NWP.EQ.0)) THEN
         
         ALLOCATE(nodalattr_here%STARTDRY(global_here%NP))

!.....Open unit 12 file

         OPEN(12,FILE=S%DIRNAME//'/'//'fort.12')

!.....Read nodalattr_here%startdry info from unit 12

         READ(12,'(A24)') global_here%AGRID2
         READ(12,*) global_here%NE2,global_here%NP2

!.....Check that global_here%NE2 and global_here%NP2 mathe with grid file

         IF ((global_here%NE2.NE.global_here%NE).OR.(global_here%NP2.NE.global_here%NP)) THEN
            IF ((global_here%NSCREEN.EQ.1).AND.(s%MYPROC.EQ.0)) WRITE(6,9900)
            WRITE(16,9900)
 9900       FORMAT(////,1X,'!!!!!!!!!!  FATAL ERROR  !!!!!!!!!',&
           //,1X,'THE PARAMETER global_here%NE2 AND global_here%NP2 MUST MATCH global_here%NE AND global_here%NP ',&
           /,1X,'USER MUST CHECK FORT.12 INPUT FILE ',&
           //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
            STOP
         ENDIF

!.....Read in nodalattr_here%startdry code values

         DO I = 1,global_here%NP
            READ(12,*) global_here%JKI,global_here%DUM1,global_here%DUM2,nodalattr_here%STARTDRY(global_here%JKI)
            IF (dg_here%MODAL_IC.NE.3) THEN
               IF (nodalattr_here%STARTDRY(global_here%JKI).EQ.-88888) THEN
                  nodalattr_here%STARTDRY(global_here%JKI) = 1
               ELSE
                  nodalattr_here%STARTDRY(global_here%JKI) = 0
               ENDIF
            ENDIF
            IF (global_here%JKI.NE.I) THEN
               IF ((global_here%NSCREEN.EQ.1).AND.(s%MYPROC.EQ.0)) WRITE(6,99805)
               WRITE(16,99805)
99805          FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',&
              'INPUT ERROR  !!!!!!!!!',&
              //,1X,'YOUR NODE NUMBERING IS NOT SEQUENTIAL ',&
              'CHECK YOUR UNIT 12 INPUT FILE CAREFULLY',//)
            ENDIF
         ENDDO

!.....Close unite 12 file

         CLOSE(12)

!.....If not using startup elevation file

      ENDIF


!.....Reset nodalattr_here%tau0var values based on input values of nodalattr_here%startdry and
!.....automatic selection of nodalattr_here%tau0 on a processor

!     DO I = 1,global_here%NP
!     IF (nodalattr_here%TAU0.LT.0.D0) THEN
!     IF (global_here%DP(I).LE.10.D0) nodalattr_here%TAU0VAR(I) = 0.02D0
!     IF (global_here%DP(I).GT.10.D0) nodalattr_here%TAU0VAR(I) = 0.005D0
!     IF (nodalattr_here%STARTDRY(I).EQ.-77777) nodalattr_here%TAU0VAR(I) = 0.02D0
!     IF (nodalattr_here%STARTDRY(I).EQ.-88888) nodalattr_here%TAU0VAR(I) = 0.02D0
!     WRITE(16,248) MYPROC,I,nodalattr_here%TAU0VAR(I)
!     248       FORMAT(/,' myproc = ',I6,' node = ',I8,
!     &             ' nodalattr_here%tau0 set to ',F12.6,/)
!     ELSE
!     C          nodalattr_here%TAU0VAR(I) = nodalattr_here%TAU0
!     ENDIF
!     ENDDO

!.....Output to unit 16 grid information including global_here%AGRID, NE, global_here%NP, global_here%H0,
!.....and nodal coordinates and bathymetry

      WRITE(16,2039) global_here%AGRID
 2039 FORMAT(/,5X,'GRID IDENTIFICATION : ',A24,/)
      IF(global_here%NSTARTDRY.EQ.1) WRITE(16,2038) global_here%AGRID2
 2038 FORMAT(5X,'nodalattr_here%STARTDRY FILE IDENTIFICATION : ',A24,/)
      WRITE(16,3) global_here%NP
    3 FORMAT(5X,'TOTAL NUMBER OF NODES =',I6,/)
      WRITE(16,4) global_here%NE
    4 FORMAT(5X,'TOTAL NUMBER OF ELEMENTS =',I6,/)
      IF(global_here%ICS.EQ.2) WRITE(16,13) global_here%SLAM0*RAD2DEG,global_here%SFEA0*RAD2DEG
 13   FORMAT(5X,'LONGITUDE ABOUT WHICH CPP PROJECTION IS CENTERED',&
     '  global_here%SLAM0 = ',F9.4,' DEGREES',&
     /,5X,'LATITUDE  ABOUT WHICH CPP PROJECTION IS CENTERED',&
     '  global_here%SFEA0 = ',F9.4,' DEGREES',/)
      IF (global_here%NSTARTDRY.EQ.0) THEN
         IF (global_here%NABOUT.NE.1) THEN
            WRITE(16,24)
 24         FORMAT(/,1X,'NODAL COORDINATES AND BATHYMETRY :')
            IF (global_here%ICS.EQ.1) THEN
               IF ((global_here%NTIP.EQ.0).AND.(global_here%NCOR.EQ.0)) THEN
                  WRITE(16,25)
 25               FORMAT(/,10X,'NODE NO.',10X,'global_here%X',20X,'global_here%Y',15X,'global_here%DP',/)
                  DO I = 1,global_here%NP
                     WRITE (16,2008) I,global_here%X(I),global_here%Y(I),global_here%DP(I)
 2008                FORMAT(5X,I6,2(2X,F20.2),2X,F12.2)
                  ENDDO
               ELSE
                  WRITE(16,9195)
 9195             FORMAT(/,1X,'   NODE ',7X,'global_here%X',14X,'global_here%Y',9X,&
                 'LAMBDA(DEG)',6X,'FEA(DEG)',9X,'global_here%DP',/)
                  DO I = 1,global_here%NP
                     WRITE (16,9197) I,global_here%X(I),global_here%Y(I),global_here%SLAM(I)*RAD2DEG,&
                    global_here%SFEA(I)*RAD2DEG,global_here%DP(I)
 9197                FORMAT(1X,I6,2(1X,F14.1),1X,2(1X,E15.7),1X,F8.2)
                  ENDDO
               ENDIF
            ELSE
               WRITE(16,9225)
 9225          FORMAT(/,1X,'   NODE ',2X,'LAMBDA(DEG)',5X,'FEA(DEG)',11X,&
              'XCP',14X,'YCP',11X,'global_here%DP',/)
               DO I = 1,global_here%NP
                  WRITE (16,9228) I,global_here%SLAM(I)*RAD2DEG,global_here%SFEA(I)*RAD2DEG,&
                 global_here%X(I),global_here%Y(I),global_here%DP(I)
 9228             FORMAT(1X,I6,2(1X,F14.8),2(1X,F15.1),1X,F10.2)
               ENDDO
            ENDIF
         ELSE
            WRITE(16,3511)
 3511       FORMAT(/,5X,'NODAL COORDINATES AND BATHYMETRY',&
           ' INFORMATION IS AVAILABLE IN THE',&
           /,6X,'UNIT 14 INPUT FILE')
         ENDIF
      ELSE
         IF (global_here%NABOUT.NE.1) THEN
            WRITE(16,24)
            IF (global_here%ICS.EQ.1) THEN
               IF ((global_here%NTIP.EQ.0).AND.(global_here%NCOR.EQ.0)) THEN
                  WRITE(16,3527)
 3527             FORMAT(/,10X,'NODE NO.',10X,'global_here%X',20X,'global_here%Y',15X,'global_here%DP',&
                 5X,'nodalattr_here%STARTDRY',/)
                  DO I = 1,global_here%NP
                     IF (nodalattr_here%STARTDRY(I).EQ.-88888.D0) THEN
                        WRITE (16,3529) I,global_here%X(I),global_here%Y(I),global_here%DP(I),nodalattr_here%STARTDRY(I)
 3529                   FORMAT(5X,I6,2(2X,F20.2),2X,F12.2,2X,F12.0)
                     ELSE
                        WRITE (16,2008) I,global_here%X(I),global_here%Y(I),global_here%DP(I)
                     ENDIF
                  ENDDO
               ELSE
                  WRITE(16,3530)
 3530             FORMAT(/,1X,'   NODE ',7X,'global_here%X',14X,'global_here%Y',9X,&
                 'LAMBDA(DEG)',6X,'FEA(DEG)',9X,'global_here%DP',&
                 5X,'nodalattr_here%STARTDRY',/)
                  DO I = 1,global_here%NP
                     IF (nodalattr_here%STARTDRY(I).EQ.-88888.D0) THEN
                        WRITE (16,3531) I,global_here%X(I),global_here%Y(I),global_here%SLAM(I)*RAD2DEG,&
                       global_here%SFEA(I)*RAD2DEG,global_here%DP(I),nodalattr_here%STARTDRY(I)
 3531                   FORMAT(1X,I6,2(1X,F14.1),1X,2(1X,E15.7),1X,F8.2,&
                       1X,F10.0)
                     ELSE
                        WRITE (16,9197) I,global_here%X(I),global_here%Y(I),global_here%SLAM(I)*RAD2DEG,&
                       global_here%SFEA(I)*RAD2DEG,global_here%DP(I)
                     ENDIF
                  ENDDO
               ENDIF
            ELSE
               WRITE(16,3535)
 3535          FORMAT(/,1X,'   NODE ',2X,'LAMBDA(DEG)',5X,'FEA(DEG)',11X,&
              'XCP',14X,'YCP',11X,'global_here%DP',&
              5X,'nodalattr_here%STARTDRY',/)
               DO I = 1,global_here%NP
                  IF (nodalattr_here%STARTDRY(I).EQ.-88888.D0) THEN
                     WRITE (16,3537) I,global_here%SLAM(I)*RAD2DEG,global_here%SFEA(I)*RAD2DEG,&
                    global_here%X(I),global_here%Y(I),global_here%DP(I),nodalattr_here%STARTDRY(I)
 3537                FORMAT(1X,I6,2(1X,F14.8),2(1X,F15.1),1X,F10.2,2X,F10.0)
                  ELSE
                     WRITE (16,9228) I,global_here%SLAM(I)*RAD2DEG,global_here%SFEA(I)*RAD2DEG,&
                    global_here%X(I),global_here%Y(I),global_here%DP(I)
                  ENDIF
               ENDDO
            ENDIF
         ELSE
            WRITE(16,3540)
 3540       FORMAT(/,5X,'NODAL COORDINATES AND BATHYMETRY',&
           ' INFORMATION IS AVAILABLE IN THE',&
           /,6X,'UNIT 14 AND 12 INPUT FILES')
         ENDIF
      ENDIF

!.....Output to unit 16 the global connectivity table (node numbers for
!.....elements

      IF (global_here%NABOUT.NE.1) THEN
         WRITE(16,26)
 26      FORMAT(//,5X,'GLOBAL NODE NUMBERS FOR EACH ELEMENT :')
         WRITE(16,27)
 27      FORMAT(/,9X,'ELEMENT',8X,'global_here%N1',9X,'global_here%N2',10X,'global_here%N3',/)
         DO I = 1,global_here%NE
            WRITE(16,2009) I,global_here%NM(I,1),global_here%NM(I,2),global_here%NM(I,3)
 2009       FORMAT(8X,4(I7,4X))
         ENDDO
      ELSE
         WRITE(16,3512)
 3512    FORMAT(/,5X,'THE GLOBAL CONNECTIVITY TABLE',&
        ' INFORMATION IS AVAILABLE IN THE',&
        /,6X,'UNIT 14 INPUT FILE')
      ENDIF
      
!.....Read information concerning bottom friction coefficient
!.....If nodalattr_here%NWP = 1, input nodal friction coefficients from unit 21
!.....If nodalattr_here%NWP = 2, set nodal friction coefficients equal to nodalattr_here%Cf
!.....If nodalattr_here%NWP = 3

!...  
!...  READ INFORMATION CONCERNING BOTTOM FRICTION COEFFICIENT
!...  IF nodalattr_here%NWP=1, INPUT NODAL FRICTION COEFFICIENTS FROM UNIT 21
!...  IF nodalattr_here%NWP=0, SET NODAL FRICTION COEFFICIENTS EQUAL TO nodalattr_here%CF
!...  IF nodalattr_here%NWP=2, READ ADDITIONAL FRICTIONAL PARAMETERS FOR BRIDGE PILINGS
!...  
      WRITE(16,1112)
      WRITE(16,2045)
 2045 FORMAT(//,' BOTTOM FRICTION INFORMATION',//)

      nodalattr_here%HBREAK=1.
      nodalattr_here%FTHETA=1.
      nodalattr_here%FGAMMA=1.
      IF(nodalattr_here%NOLIBF.EQ.0) READ(15,*) nodalattr_here%TAU
      nodalattr_here%CF=nodalattr_here%TAU
      IF(nodalattr_here%NOLIBF.EQ.1) READ(15,*) nodalattr_here%CF
      IF(nodalattr_here%NOLIBF.EQ.2) READ(15,*) nodalattr_here%CF,nodalattr_here%HBREAK,nodalattr_here%FTHETA,nodalattr_here%FGAMMA

      IF (nodalattr_here%NWP.EQ.0) THEN
         ALLOCATE(nodalattr_here%FRIC(global_here%NP))
         DO I=1,global_here%NP
            nodalattr_here%FRIC(I)=nodalattr_here%CF
         END DO
         IF(nodalattr_here%NOLIBF.EQ.2) THEN
            WRITE(16,101) nodalattr_here%CF,nodalattr_here%HBREAK,nodalattr_here%FTHETA,nodalattr_here%FGAMMA
 101        FORMAT(5X,'HYBRID FRICTION RELATIONSHIP PARAMTERS, CFMIN =',&
           F12.8,'  nodalattr_here%HBREAK = ',F8.2,&
           /,5X,'nodalattr_here%FTHETA = ',F8.2,'  nodalattr_here%FGAMMA = ',F10.4,//)
         ENDIF
         IF(nodalattr_here%NOLIBF.EQ.1) THEN
            WRITE(16,8) nodalattr_here%CF
 8          FORMAT(5X,'NONLINEAR FRICTION COEFFICIENT nodalattr_here%CF =',F12.8,/)
         ENDIF
         IF(nodalattr_here%NOLIBF.EQ.0) THEN
            WRITE(16,106) nodalattr_here%TAU
 106        FORMAT(5X,'LINEAR BOTTOM FRICTION nodalattr_here%TAU =',F12.8,5X,'1/sec'/)
            IF(nodalattr_here%TAU.NE.nodalattr_here%TAU0) THEN !CHECK nodalattr_here%TAU VALUE AGAINST nodalattr_here%TAU0
               IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9951)
               WRITE(16,9951)
 9951          FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',&
              'INPUT ERROR  !!!!!!!!!',&
              //,1X,'TYPICALLY YOUR INPUT VALUE FOR ',&
              'nodalattr_here%TAU0 SHOULD BE SET EQUAL TO nodalattr_here%TAU')
            ENDIF
         ENDIF
      ENDIF


!     IF(nodalattr_here%NWP.EQ.1) THEN
!     OPEN(21,FILE=DIRNAME//'/'//'fort.21')
!     READ(21,'(A20)') global_here%AFRIC
!     DO I=1,global_here%NP
!     READ(21,*) global_here%NHG,nodalattr_here%FRIC(global_here%NHG)
!     IF(global_here%NHG.global_here%NE.I) THEN
!     IF(global_here%NSCREEN.EQ.1.AND.MYPROC.EQ.0) WRITE(6,99803)
!     WRITE(16,99803)
!99   803       FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ',
!     &                     'INPUT ERROR  !!!!!!!!!',
!     &        //,1X,'YOUR NODAL FRICTION NUMBERING IS NOT SEQUENTIAL ',
!     &        /,1X,'CHECK YOUR UNIT 21 INPUT FILE CAREFULLY',
!     &        //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
!     STOP
!     ENDIF
!     END DO
!     WRITE(16,3601) global_here%AFRIC
!     3601   FORMAT(/,5X,'FRICTION FILE IDENTIFICATN : ',A20,/)
!     IF(global_here%NABOUT.global_here%NE.1) THEN
!     WRITE(16,2080)
!     2080     FORMAT(/,10X,'NODE',5X,'BOTTOM FRICTION nodalattr_here%FRIC',5X,/)
!     DO I=1,global_here%NP
!     WRITE(16,2087) I,nodalattr_here%FRIC(I)
!     2087       FORMAT(7X,I6,6X,E15.10)
!     END DO
!     ELSE
!     WRITE(16,3504)
!     3504     FORMAT(/,5X,'NODAL BOTTOM FRICTION VALUES ARE AVAILABLE',
!     &           /,6X,' IN UNIT 21 INPUT FILE')
!     ENDIF
!     ENDIF

!     IF(nodalattr_here%NWP.EQ.2) THEN
!     CALL ALLOC_MAIN13()   !allocate bridge piling arrays
!     DO I=1,global_here%NP
!     global_here%NBNNUM(I)=0
!     global_here%BK(I)=0.d0
!     global_here%BALPHA(I)=0.d0
!     global_here%BDELX(I)=1.d0
!     ENDDO
!     OPEN(21,FILE=DIRNAME//'/'//'fort.21')
!     READ(21,'(A20)') global_here%AFRIC
!     READ(21,*) global_here%NBPNODES
!     DO I=1,global_here%NBPNODES
!     READ(21,*) global_here%NBNNUM(I),global_here%BK(I),global_here%BALPHA(I),global_here%BDELX(I),global_here%POAN
!     global_here%BDELX(I)=4.d0*global_here%BDELX(I)/global_here%POAN
!     ENDDO
!     WRITE(16,3602) global_here%AFRIC
!     3602   FORMAT(/,5X,'BRIDGE PIER FRICTION FILE IDENTIFICATN : ',A20,/)
!     IF(global_here%NABOUT.global_here%NE.1) THEN
!     WRITE(16,2081)
!     2081     FORMAT(/,10X,'NODE',3X,'PIER SHAPE FACTOR',3X,
!     &                 'CONSTRICTION FRACTION',3X,'EFFECTIVE global_here%DELX'/)
!     DO I=1,global_here%NBPNODES
!     WRITE(16,2082) global_here%NBNNUM(I),global_here%BK(I),global_here%BALPHA(I),global_here%BDELX(I)
!     2082       FORMAT(5X,I8,10X,F7.3,12X,F7.3,13X,F9.3)
!     END DO
!     ELSE
!     WRITE(16,2083)
!     2083     FORMAT(/,5X,'BRIDGE PILING FRICTION VALUES ARE AVAILABLE',
!     &           /,6X,' IN UNIT 21 INPUT FILE')
!     ENDIF
!     ENDIF

!...  
!...  READ IN AND WRITE OUT EDDY VISCOSITY/DIFFUSIVITY COEFFICIENTS
!...  


      IF (global_here%IM.EQ.10) THEN
         READ(15,*) nodalattr_here%ESLM,nodalattr_here%ESLC
         DO I=1,global_here%NP
            nodalattr_here%EVM(I)=nodalattr_here%ESLM
            nodalattr_here%EVC(I)=nodalattr_here%ESLC
         END DO
         WRITE(16,111) nodalattr_here%ESLM,nodalattr_here%ESLC
 111     FORMAT(5X,'nodalattr_here%EVM, EDDY VISCOSITY COEFFICIENT =',E15.8,/,&
        5X,'nodalattr_here%EVC, EDDY DIFFUSIVITY COEFFICIENT =',E15.8,//)
      ELSE
         READ(15,*) nodalattr_here%ESLM
         IF (nodalattr_here%NWP.EQ.0) THEN
            ALLOCATE(nodalattr_here%EVM(global_here%NP))
            DO I=1,global_here%NP
               nodalattr_here%EVM(I)=nodalattr_here%ESLM
            END DO
            WRITE(16,11) nodalattr_here%ESLM
 11         FORMAT(5X,'nodalattr_here%EVM, EDDY VISCOSITY COEFFICIENT =',E15.8,//)
         ENDIF
      ENDIF


!...  02/19/2007 s.b.
      global_here%EVMSUM = 0.D0
      IF (nodalattr_here%NWP.EQ.0) THEN
         DO I=1,global_here%NP
            global_here%EVMSUM = global_here%EVMSUM + ABS(nodalattr_here%EVM(I))
         ENDDO
      ENDIF
      
!     
!     ek: Initialize nodal attributes, now that grid has been read
!     in from unit 14 file.

      IF (nodalattr_here%NWP.GT.0)&
     CALL InitNodalAttr(global_here,nodalattr_here,global_here%DP, global_here%NP, global_here%G, global_here%NScreen, global_here%ScreenUnit,S%MYPROC,global_here%NAbOut)
      


!...  
!...  READ CORIOLIS INFORMATION AND COMPUTE THE CORIOLIS VECTOR
!...  OUTPUT RESULTING CORIOLIS INFORMATION
!...  
      WRITE(16,1112)
      WRITE(16,2090)
 2090 FORMAT(//,1X,'CORIOLIS INFORMATION ',//)

      READ(15,*) global_here%CORI
      IF(global_here%NCOR.EQ.0) THEN
         DO I=1,global_here%NP
            global_here%CORIF(I)=global_here%CORI
         END DO
      ENDIF
      IF(global_here%NCOR.EQ.1) THEN
         DO I=1,global_here%NP
            global_here%CORIF(I)=2.0d0*7.29212d-5*SIN(global_here%SFEA(I))
         END DO
      ENDIF

      IF(global_here%NCOR.EQ.0) THEN
         WRITE(16,12) global_here%CORI
 12      FORMAT(5X,'CONSTANT CORIOLIS COEFFICIENT =',E15.8,5X,'1/SEC',/)
      ENDIF
      IF(global_here%NCOR.EQ.1) THEN
         WRITE(16,3604)
 3604    FORMAT(/,5X,'LATITUDES ARE USED TO COMPUTE VARIABLE CORIOLIS',&
        /,7X,'AND ARE BASED ON INPUT NODAL COORDINATES',/)
         IF(global_here%NABOUT.NE.1) THEN
            WRITE(16,2092)
 2092       FORMAT(/,10X,' NODE ',5X,'NODAL CORIOLIS global_here%CORIF',/)
            DO I=1,global_here%NP
               WRITE(16,2096) I,global_here%CORIF(I)
 2096          FORMAT(7X,I6,10X,E16.9)
            END DO
         ENDIF
      ENDIF

!...  
!...  READ AND PROCESS INFORMATION ABOUT THE TIDAL POTENTIAL CONSTITUENTS
!...  

      READ(15,*) global_here%NTIF
      s%MNTIF = global_here%ntif
      if (global_here%ntif .eq. 0) s%MNTIF = 1

!.... allocate tidal potential arrays

      call alloc_main4a(s,global_here)

!.... READ TIDAL POTENTIAL AMPLITUDE, FREQUENCIES, NODAL FACTORS,
!.... EQUILIBRIUM ARGUMENTS AND ALPHANUMERIC LABEL
!.... 
      DO I=1,global_here%NTIF
         READ(15,'(A5)')  global_here%TIPOTAG(I)
         READ(15,*)  global_here%TPK(I),global_here%AMIGT(I),global_here%ETRF(I),global_here%FFT(I),global_here%FACET(I)
         IF(global_here%AMIGT(I).EQ.0.) THEN
            global_here%PERT(I)=0.
         ELSE
            global_here%PERT(I)=2.D0*PI/global_here%AMIGT(I)
         ENDIF
      END DO

!...  LINES TO USE EARTH LOAD/SELF-ATTRACTION PART OF TIDAL POTENTIAL FORCING

      CALL ALLOC_MAIN4b(s,global_here)
      IF(global_here%NTIP.EQ.2) THEN
         OPEN(24,FILE='fort.24')
         DO I=1,global_here%NTIF
            READ(24,9930)
 9930       FORMAT(///)
            DO J=1,global_here%NP
               READ(24,*) JJ,global_here%SALTAMP(I,JJ),global_here%SALTPHA(I,JJ)
               global_here%SALTPHA(I,JJ)=global_here%SALTPHA(I,JJ)*DEG2RAD
            END DO
         END DO
      ELSE
         DO I=1,global_here%NTIF
            DO J=1,global_here%NP
               global_here%SALTAMP(I,J)=0.d0
               global_here%SALTPHA(I,J)=0.d0
            END DO
         END DO
         CLOSE(24)
      ENDIF

!...  
!...  OUTPUT TO UNIT 16 INFORMATION ABOUT TIDAL POTENTIAL FORCING
!.... OUTPUT WILL VARY DEPENDING ON VALUES OF global_here%NTIP,global_here%NTIF AND global_here%NCOR
!...  
      WRITE(16,1112)
      WRITE(16,2102)
 2102 FORMAT(//,1X,'TIDAL POTENTIAL FORCING INFORMATION ',//)
      WRITE(16,22) global_here%NTIF
 22   FORMAT(/,1X,'TIDAL POTENTIAL IS FORCED FOR ',I5,&
     ' CONSTITUENT(S) ')
      IF(global_here%NTIF.GT.0) WRITE(16,23)
 23   FORMAT(/,1X,'AMPLITUDE',4X,'FREQUENCY',5X,&
     '    global_here%ETRF      ','NODAL FACTOR',2X,&
     'EQU.global_here%ARG(DEG)',1X,'CONSTITUENT',/)
      DO I=1,global_here%NTIF
         WRITE(16,2107) global_here%TPK(I),global_here%AMIGT(I),global_here%ETRF(I),global_here%FFT(I),global_here%FACET(I),&
        global_here%TIPOTAG(I)
 2107    FORMAT(1X,F10.7,1X,F15.12,2X,F10.7,5X,F10.7,1X,F10.3,7X,A5)
      END DO
!...  
!...  CONVERT global_here%FACET(I) VALUES FROM DEGREES TO RADIANS
!...  
      DO I=1,global_here%NTIF
         global_here%FACET(I)=global_here%FACET(I)*DEG2RAD
      END DO
!...  
!...  CHECK CONSISTENCY OF INPUT PARAMETERS global_here%NTIF AND global_here%NTIP
!...  
      IF(((global_here%NTIP.EQ.0).AND.(global_here%NTIF.NE.0)).OR.((global_here%NTIP.NE.0).AND.&
     (global_here%NTIF.EQ.0))) THEN
         IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9961)
         WRITE(16,9961)
 9961    FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',&
        'INPUT ERROR  !!!!!!!!!',&
        //,1X,'YOUR SELECTION OF global_here%NTIF AND global_here%NTIP (UNIT 15 INPUT ',&
        'PARAMETERS) IS INCONSISTENT',&
        /,1X,'PLEASE CHECK THESE VALUES')
         IF(global_here%NFOVER.EQ.1) THEN
            IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9987)
            WRITE(16,9987)
 9987       FORMAT(/,1X,'PROGRAM WILL OVERRIDE THE SPECIFIED ',&
           'INPUT AND NEGLECT TIDAL POTENTIAL TERMS',&
           /,1X,' AND/OR RESET global_here%NTIP = 0',&
           //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)
            global_here%NTIP=0
         ELSE
            IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9973)
            WRITE(16,9973)
            STOP
         ENDIF
         GOTO 1893
      ENDIF
!...  
!...  PRINT OUT LAT/LON VALUES TO BE USED IN COMPUTING TIDAL POTENTIAL
!.... IF NOT ALREADY DONE SO IN CORIOLIS SECTION AND TIDAL POTENTIAL
!.... IS ACTIVATED WITH global_here%NTIP=1
!...  

      IF(global_here%NTIP.GE.1) THEN
         IF(global_here%ICS.EQ.1) THEN
            WRITE(16,3605)
 3605       FORMAT(/,5X,'LONGITUDES AND LATITUDES ARE USED TO',&
           ' COMPUTE THE TIDAL POTENTIAL FUNCTION',&
           /,7X,'AND ARE BASED ON AN INVERSE CPP PROJECTION ',&
           'OF THE INPUT COORDINATES',/)
         ELSE
            WRITE(16,2109)
 2109       FORMAT(/,5X,'LONGITUDES AND LATITUDES ARE USED TO',&
           ' COMPUTE THE TIDAL POTENTIAL FUNCTION',&
           /,7X,'AND ARE BASED ON INPUT NODAL COORDINATES ',/)
         ENDIF
      ENDIF

!...  
!...  INPUT FROM UNIT 15 THE TIDAL FORCING FREQUENCIES ON THE ELEVATION
!.... SPECIFIED BOUNDARIES: INCLUDING global_here%NBFR, FREQUENCIES, NODAL FACTORS,
!.... EQUILIBRIUM ARGUMENTS AND AN ELEVATION BOUNDARY CONDITION
!.... ALPHANUMERIC DESCRIPTOR
!...  
 1893 READ(15,*) global_here%NBFR
      s%MNBFR = global_here%NBFR
      IF (global_here%NBFR.EQ.0) s%MNBFR = 1

!     - Allocate arrays dimensioned by MNBFR

      call alloc_main5(s,global_here)

      WRITE(16,1112)
      WRITE(16,2106)
 2106 FORMAT(//,1X,'ELEVATION SPECIFIED BOUNDARY FORCING INFORMATION '&
     ,//)
      WRITE(16,20) global_here%NBFR
 20   FORMAT(/,5X,'NUMBER OF PERIODIC, ELEVATION SPECIFIED ',&
     'CONSTITUENTS =',I5)
      IF(global_here%NBFR.GE.1) WRITE(16,21)
 21   FORMAT(/,7X,'CONSTITUENT #',4X,'FREQUENCY',4X,'NODAL FACTOR',&
     3X,'EQU.global_here%ARG (DEG)',2X,'CONSTITUENT',/)
      DO I=1,global_here%NBFR
         READ(15,'(A5)') global_here%BOUNTAG(I)
         READ(15,*) global_here%AMIG(I),global_here%FF(I),global_here%FACE(I)
         WRITE(16,1850) I,global_here%AMIG(I),global_here%FF(I),global_here%FACE(I),global_here%BOUNTAG(I)
 1850    FORMAT(12X,I2,6X,F16.12,2X,F10.7,2X,F10.3,10X,A5)
         global_here%FACE(I)=global_here%FACE(I)*DEG2RAD
         IF(global_here%AMIG(I).EQ.0.) THEN
            global_here%PER(I)=0.
         ELSE
            global_here%PER(I)=2.D0*PI/global_here%AMIG(I)
         ENDIF
      ENDDO
!...  
!...  INPUT ELEVATION BOUNDARY FORCING NODE NUMBER INFORMATION FROM UNIT 14 AND
!.... OUTPUT TO UNIT 16
!...  
!...  INPUT THE TOTAL NUMBER OF ELEVATION BOUNDARY SEGMENTS
!...  
      READ(14,*) global_here%NOPE

      WRITE(16,1852) global_here%NOPE
 1852 FORMAT(///,5X,'TOTAL NUMBER OF ELEVATION BOUNDARY FORCING',&
     ' SEGMENTS ',' = ',I5)
!...  
!...  INPUT THE TOTAL NUMBER OF ELEVATION BOUNDARY NODES
!...  
      READ(14,*) global_here%NETA
      WRITE(16,1854) global_here%NETA
 1854 FORMAT(/,5X,'TOTAL NUMBER OF ELEVATION SPECIFIED BOUNDARY NODES ='&
     ,I6)

!     allocate arrays dimensioned by global_here%NOPE and global_here%NETA
      s%MNOPE = global_here%NOPE
      IF (global_here%NOPE.EQ.0) s%MNOPE = 1
      s%MNETA = global_here%NETA
      IF (global_here%NETA.EQ.0) s%MNETA = 1

      call alloc_main2(s,global_here)     
!...  
!...  INPUT THE NODE NUMBERS ON EACH ELEVATION BOUNDARY FORCING SEGMENT
!...  
      s%MNEI=0
      global_here%JNMM=0
      DO K=1,global_here%NOPE
         READ(14,*) global_here%NVDLL(K)
         WRITE(16,281) K,global_here%NVDLL(K)
 281     FORMAT(//,5X,'TOTAL NUMBER OF NODES ON ELEVATION SPECIFIED ',&
        'BOUNDARY SEGMENT ',2X,I2,2X,'=',1X,I5,/)
         DO I=1,global_here%NVDLL(K)
            READ(14,*) global_here%NBDV(K,I)
            WRITE(16,1855) global_here%NBDV(K,I)
 1855       FORMAT(7X,I7)
            IF (global_here%NNEIGH(global_here%NBDV(K,I)).NE.0) THEN
               global_here%NNEIGH(global_here%NBDV(K,I))=global_here%NNEIGH(global_here%NBDV(K,I))+1
               IF (global_here%NNEIGH(global_here%NBDV(K,I)).GT.s%MNEI) s%MNEI=global_here%NNEIGH(global_here%NBDV(K,I))
               global_here%NNEIGH(global_here%NBDV(K,I)) = 0
            ENDIF

            global_here%NBD(global_here%JNMM+I)=global_here%NBDV(K,I)
         ENDDO
         global_here%JNMM=global_here%JNMM+global_here%NVDLL(K)
      ENDDO
!...  
!...  CHECK TO MAKE SURE THAT global_here%JNMM EQUALS global_here%NETA
!...  
      IF(global_here%NETA.NE.global_here%JNMM) THEN
         IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9945)
         WRITE(16,9945)
 9945    FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL INPUT ERROR ',&
        '!!!!!!!!!',&
        //,1X,'THE INPUT PARAMETER global_here%NETA FROM UNIT 14 DOES NOT MATCH ',&
        'THE TOTAL NUMBER OF BOUNDARY NODES',&
        /,1X,' FROM ALL THE SPECIFIED SEGMENTS COMBINED')
         IF(global_here%NFOVER.EQ.1) THEN
            IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9989)
            WRITE(16,9989)
 9989       FORMAT(/,1X,'THE PROGRAM WILL NOW CORRECT THIS ERROR',&
           /,1X,'PLEASE CHECK YOUR INPUT CAREFULLY !!!',&
           //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)
            global_here%NETA=global_here%JNMM
         ELSE
            IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9973)
            WRITE(16,9973)
            STOP
         ENDIF
      ENDIF
!...  
!...  SET UP TO READ IN TIME SERIES ELEVATION SPECIFIED BOUNDARY CONDITIONS IF APPROPRIATE
!...  
      IF((global_here%NBFR.EQ.0).AND.(global_here%NOPE.GT.0)) THEN
         WRITE(16,1871)
 1871    FORMAT(/,5X,'TIME SERIES ELEVATION SPECIFIED VALUES WILL BE ',&
        'READ FROM UNIT 19',&
        /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE ',&
        /,9X,'ELEVATION DATA WITH THE MODEL TIME STEP.')
      ENDIF

!...  
!...  INPUT FORCING CONDITIONS ON PERIODIC ELEVATION SPECIFIED BOUNDARIES FOR EACH
!...  OF THE ELEVATION FORCING FREQUENCIES FROM UNIT 15 AND OUTPUT TO UNIT 16
!...  

      DO I=1,global_here%NBFR
         WRITE(16,29) I,global_here%BOUNTAG(I)
 29      FORMAT(////,5X,'ELEVATION BOUNDARY TIDAL FORCING FOR',&
        ' CONSTITUENT NUMBER',I4,1X,'DESIGNATED : ',A5)
         READ(15,'(A10)') global_here%ALPHA
         WRITE(16,31) global_here%ALPHA
 31      FORMAT(9X,'VERIFICATION OF CONSTITUENT : ',A10,/)
         WRITE(16,30)
 30      FORMAT(14X,'NODE',11X,'AMPL.',9X,'PHASE(DEG)',/)
         DO J=1,global_here%NETA
            READ(15,*) global_here%EMO(I,J),global_here%EFA(I,J)
            WRITE(16,1870) global_here%NBD(J),global_here%EMO(I,J),global_here%EFA(I,J)
 1870       FORMAT(10X,I8,4X,F14.5,4X,F12.3)
            global_here%EFA(I,J) = global_here%EFA(I,J)*DEG2RAD
         ENDDO
!     DO J=1,global_here%NETA
!     READ(15,*) global_here%UMO(I,J),global_here%UFA(I,J)
!     global_here%UFA(I,J) = global_here%UFA(I,J)*DEG2RAD
!     ENDDO
!     DO J=1,global_here%NETA
!     READ(15,*) global_here%VMO(I,J),global_here%VFA(I,J)
!     global_here%VFA(I,J) = global_here%VFA(I,J)*DEG2RAD
!     ENDDO
      ENDDO

!.....READ THE MINIMUM INNER ANGLE FOR WHICH VELOCITY AT FLOW BOUNDARY NODES
!.....WILL BE ZEROED IN THE TANGENTIAL DIRECTIONS WHEN NORMAL FLOW IS AN
!.....ESSENTIAL B.C.

      READ(15,*) global_here%ANGINN
      WRITE(16,1112)
      WRITE(16,7654) global_here%ANGINN
 7654 FORMAT(//,5X,'global_here%ANGINN = ',F8.2,' DEGREES',&
     /,5X,'ALL FLOW BOUNDARY NODES WITH NORMAL FLOW AS AN ',&
     'ESSENTIAL B.C. AND ',&
     /,9X,'INNER ANGLES LESS THAN global_here%ANGINN WILL HAVE BOTH NORMAL ',&
     /,9X,'AND TANGENTIAL VELOCITY COMPONENTS ZEROED',/)
      global_here%COSTSET=COS(global_here%ANGINN*DEG2RAD)

!...  
!...  INPUT FLOW BOUNDARY INFORMATION FROM UNIT 14 AND OUTPUT TO UNIT 16
!...  

!.....INTERIOR NODES, global_here%LBCODE=-1, COS=0, SIN=1
!.....BOUNDARY NODES, global_here%LBCODE=global_here%LBCODEI=global_here%IBTYPE,
!.....COS & SIN DETERMINED FROM NORMAL DIRECTION IN ALL CASES, ALTHOUGH THIS
!.......INFORMATION IS ONLY USED WHEN NORMAL FLOW IS AN ESSENTIAL B.C. AND
!.......FREE TANGENTIAL SLIP IS ALLOWED.

!.....INPUT THE TOTAL NUMBER OF FLOW BOUNDARY SEGMENTS

      WRITE(16,1112)
      WRITE(16,1878)
 1878 FORMAT(//,1X,'FLOW BOUNDARY INFORMATION ',/)
      READ(14,*) global_here%NBOU

      WRITE(16,1879) global_here%NBOU
 1879 FORMAT(//,5X,'THE TOTAL NUMBER OF FLOW BOUNDARY SEGMENTS = ',I5)

!.....INPUT THE TOTAL NUMBER OF FLOW BOUNDARY NODES

      READ(14,*) global_here%NVEL
      WRITE(16,1881) global_here%NVEL
 1881 FORMAT(/,5X,'THE TOTAL NUMBER OF FLOW BOUNDARY NODES = ',I5)

      s%MNBOU = global_here%NBOU
      IF (global_here%NBOU.EQ.0) s%MNBOU = 1
      s%MNVEL = global_here%NVEL*2            !Cvjp  -  11/28/99 -  upper bound guess for MNVEL

!.....Allocate space for nonperiodic zero and nonzero normal flow boundary arrays
!.....including barriers
      call alloc_main3(s,global_here)

!.....INPUT THE NUMBER OF NODES IN THE NEXT FLOW BOUNDARY SEGMENT
!.....AND THE BOUNDARY TYPE

      global_here%JGW=0
      global_here%JME=0
      global_here%NFLUXF=0
      global_here%NFLUXB=0
      global_here%NFLUXIB=0
      global_here%NFLUXIBP=0
      global_here%NVELEXT=0
      dg_here%NIBSEG = 0
      dg_here%NEBSEG = 0
      
      CALL ALLOC_EDGES0(s,dg_here)

      DO K=1,global_here%NBOU
         READ(14,*) global_here%NVELL(K),global_here%IBTYPE
!     jcf dg - added variable global_here%SEGTYPE to record boundary segment global_here%IBTYPE
         global_here%SEGTYPE(K) = global_here%IBTYPE
!.......CHECK THAT global_here%IBTYPE PARAMETER HAS BEEN SET PROPERLY
         IF(    (global_here%IBTYPE.NE.0).AND.(global_here%IBTYPE.NE.10).AND.(global_here%IBTYPE.NE.20)&
        .AND.(global_here%IBTYPE.NE.1).AND.(global_here%IBTYPE.NE.11).AND.(global_here%IBTYPE.NE.21)&
        .AND.(global_here%IBTYPE.NE.2).AND.(global_here%IBTYPE.NE.12).AND.(global_here%IBTYPE.NE.22)&
        .AND.(global_here%IBTYPE.NE.3).AND.(global_here%IBTYPE.NE.13).AND.(global_here%IBTYPE.NE.23)&
!jj   wm001 - following 2 lines modified 
        .AND.(global_here%IBTYPE.NE.4).AND.(global_here%IBTYPE.NE.24)&
        .AND.(global_here%IBTYPE.NE.5).AND.(global_here%IBTYPE.NE.25).AND.(global_here%IBTYPE.NE.30)) THEN
         IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9985) K
         WRITE(16,9985) K
 9985    FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',&
        //,1X,'THE FLOW BOUNDARY TYPE PARAMETER global_here%IBTYPE ',&
        'HAS NOT BEEN CORRECTLY SET FOR ',&
        /,1X,'FLOW BOUNDARY SEGMENT NO. ',I8,&
        /,1X,'USER MUST CORRECT UNIT 14 INPUT FILE',&
        //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
         STOP
      ENDIF
!.......WRITE OUT INFORMATION TO UNIT 16
      IF((global_here%IBTYPE.EQ.4).OR.(global_here%IBTYPE.EQ.24)) THEN
         WRITE(16,28) K,global_here%NVELL(K),K,2*global_here%NVELL(K)
 28      FORMAT(///,5X,'TOTAL NUMBER OF PAIRS FOR FLOW BOUNDARY',&
        ' SEGMENT',2X,I2,2X,'=',2X,I5,/,&
        5X,'TOTAL NUMBER OF NODES FOR FLOW BOUNDARY',&
        ' SEGMENT',2X,I2,2X,'=',2X,I5)
      ELSE
         WRITE(16,128) K,global_here%NVELL(K)
 128     FORMAT(///,5X,'TOTAL NUMBER OF NODES FOR FLOW BOUNDARY',&
        ' SEGMENT',2X,I2,2X,'=',2X,I5)
      ENDIF
!.......CONTINUE PROCESSING FLOW BOUNDARY INFORMATION
      IF(global_here%IBTYPE.EQ.0) THEN
         WRITE(16,2340)
 2340    FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,&
        7X,'NO NORMAL FLOW AS AN ESSENTIAL B.C.',/,&
        7X,'AND FREE TANGENTIAL SLIP',/)
      ENDIF
      IF(global_here%IBTYPE.EQ.1) THEN
         WRITE(16,2341)
 2341    FORMAT(5X,'THIS SEGMENT IS AN INTERNAL BOUNDARY WITH:',/,&
        7X,'NO NORMAL FLOW AS AN ESSENTIAL B.C.',/,&
        7X,'AND FREE TANGENTIAL SLIP',/)
      ENDIF
      IF(global_here%IBTYPE.EQ.2) THEN
         global_here%NFLUXF=1
         WRITE(16,2342)
 2342    FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,&
        7X,'SPECIFIED NORMAL FLOW AS AN ESSENTIAL B.C.',/,&
        7X,'AND FREE TANGENTIAL SLIP',/)
      ENDIF
      IF(global_here%IBTYPE.EQ.3) THEN
         global_here%NFLUXB=1
         WRITE(16,2344)
 2344    FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,&
        7X,'A BARRIER WHICH ALLOWS FREE SURFACE',&
        ' SUPERCRITICAL OUTFLOW',/,&
        7X,'FROM THE DOMAIN ONCE THE BARRIER HAS BEEN',&
        ' OVERTOPPED',/,&
        7X,'AND FREE TANGENTIAL SLIP',/)
      ENDIF
      IF(global_here%IBTYPE.EQ.4) THEN
         global_here%NFLUXIB=1
         WRITE(16,2345)
 2345    FORMAT(5X,'THIS SEGMENT IS AN INTERNAL BARRIER BOUNDARY:',/,&
        7X,'WITH global_here%CROSS BARRIER FLOW TREATED AS AN ESSENTIAL ',&
        ' NORMAL FLOW BOUNDARY CONDITION',/,&
        7X,'WHICH LEAVES/ENTERS THE DOMAIN ON ONE global_here%SIDE OF ',&
        ' THE BARRIER AND ENTERS/LEAVES THE DOMAIN ON THE ',/,&
        7X,'CORRESPONDING OPPOSITE global_here%SIDE OF THE BARRIER ',/,&
        7X,'FLOW RATE AND DIRECTION ARE BASED ON BARRIER ',&
        ' HEIGHT, SURFACE WATER ELEVATION',/,&
        7X,'ON BOTH SIDES OF THE BARRIER, BARRIER COEFFICIENT',&
        ' AND THE APPROPRIATE BARRIER FLOW FORMULA',/,&
        7X,'FREE TANGENTIAL SLIP IS ALLOWED',/)
      ENDIF
      IF(global_here%IBTYPE.EQ.5) THEN
         global_here%NFLUXIB=1
         global_here%NFLUXIBP=1
         WRITE(16,2347)
 2347    FORMAT(5X,'THIS SEGMENT IS AN INTERNAL BARRIER BOUNDARY:',/,&
        7X,'WITH ADDITIONAL global_here%CROSS BARRIER PIPES ',&
        'LOCATED UNDER THE CROWN ',/,&
        7X,'global_here%CROSS BARRIER FLOW IS TREATED AS AN ESSENTIAL',&
        ' NORMAL FLOW BOUNDARY CONDITION',/,&
        7X,'WHICH LEAVES/ENTERS THE DOMAIN ON ONE global_here%SIDE OF ',&
        ' THE BARRIER AND ENTERS/LEAVES THE DOMAIN ON THE ',/,&
        7X,'CORRESPONDING OPPOSITE global_here%SIDE OF THE BARRIER ',/,&
        7X,'FLOW RATE AND DIRECTION ARE BASED ON BARRIER ',&
        ' HEIGHT, SURFACE WATER ELEVATION',/,&
        7X,'ON BOTH SIDES OF THE BARRIER, BARRIER COEFFICIENT',&
        ' AND THE APPROPRIATE BARRIER FLOW FORMULA',/,&
        7X,'IN ADDITION global_here%CROSS BARRIER PIPE FLOW RATE AND ',&
        ' DIRECTION ARE BASED ON PIPE CROWN HEIGHT, ',/,&
        7X,'SURFACE WATER ELEVATION ON BOTH SIDES OF THE ',&
        'BARRIER, PIPE FRICTION COEFFICIENT, PIPE DIAMETER',/,&
        7X,' AND THE APPROPRIATE PIPE FLOW FORMULA',/,&
        7X,'FREE TANGENTIAL SLIP IS ALLOWED',/)
      ENDIF
      IF(global_here%IBTYPE.EQ.10) THEN
         WRITE(16,2350)
 2350    FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,&
        7X,'NO NORMAL FLOW AS AN ESSENTIAL B.C.',/,&
        7X,'AND NO TANGENTIAL SLIP',/)
      ENDIF
      IF(global_here%IBTYPE.EQ.11) THEN
         WRITE(16,2351)
 2351    FORMAT(5X,'THIS SEGMENT IS AN INTERNAL BOUNDARY WITH:',/,&
        7X,'NO NORMAL FLOW AS AN ESSENTIAL B.C.',/,&
        7X,'AND NO TANGENTIAL SLIP',/)
      ENDIF
      IF(global_here%IBTYPE.EQ.12) THEN
         global_here%NFLUXF=1
         WRITE(16,2352)
 2352    FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,&
        7X,'SPECIFIED NORMAL FLOW AS AN ESSENTIAL B.C.',/,&
        7X,'AND NO TANGENTIAL SLIP',/)
      ENDIF
      IF(global_here%IBTYPE.EQ.13) THEN
         global_here%NFLUXB=1
         WRITE(16,2354)
 2354    FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,&
        7X,'A BARRIER WHICH ALLOWS FREE SURFACE',&
        ' SUPERCRITICAL OUTFLOW',/,&
        7X,'FROM THE DOMAIN ONCE THE BARRIER HAS BEEN',&
        ' OVERTOPPED',/,&
        7X,'AND NO TANGENTIAL SLIP',/)
      ENDIF
      IF(global_here%IBTYPE.EQ.20) THEN
         WRITE(16,2360)
 2360    FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,&
        7X,'NO NORMAL FLOW AS A NATURAL B.C.',/,&
        7X,'AND FREE TANGENTIAL SLIP',/)
      ENDIF
      IF(global_here%IBTYPE.EQ.21) THEN
         WRITE(16,2361)
 2361    FORMAT(5X,'THIS SEGMENT IS AN INTERNAL BOUNDARY WITH:',/,&
        7X,'NO NORMAL FLOW AS A NATURAL B.C.',/,&
        7X,'AND FREE TANGENTIAL SLIP',/)
      ENDIF
      IF(global_here%IBTYPE.EQ.22) THEN
         global_here%NFLUXF=1
         WRITE(16,2362)
 2362    FORMAT(5X,'THIS SEGMENT IS A EXTERNAL BOUNDARY WITH:',/,&
        7X,'SPECIFIED NORMAL FLOW AS A NATURAL B.C.',/,&
        7X,'AND FREE TANGENTIAL SLIP',/)
      ENDIF
      IF(global_here%IBTYPE.EQ.23) THEN
         global_here%NFLUXB=1
         WRITE(16,2356)
 2356    FORMAT(5X,'THIS SEGMENT IS AN EXTERNAL BOUNDARY WITH:',/,&
        7X,'A BARRIER WHICH ALLOWS FREE SURFACE',&
        ' SUPERCRITICAL OUTFLOW',/,&
        7X,'FROM THE DOMAIN ONCE THE BARRIER HAS BEEN',&
        ' OVERTOPPED',/,&
        7X,' IMPLEMENTED AS A NATURAL BOUNDARY CONDITION'&
        ,7X,'FREE TANGENTIAL SLIP IS ALSO ALLOWED',/)
      ENDIF
      IF(global_here%IBTYPE.EQ.24) THEN
         global_here%NFLUXIB=1
         WRITE(16,2357)
 2357    FORMAT(5X,'THIS SEGMENT IS AN INTERNAL BARRIER BOUNDARY:',/,&
        7X,'WITH global_here%CROSS BARRIER FLOW TREATED AS A NATURAL ',&
        ' NORMAL FLOW BOUNDARY CONDITION',/,&
        7X,'WHICH LEAVES/ENTERS THE DOMAIN ON ONE global_here%SIDE OF ',&
        ' THE BARRIER AND ENTERS/LEAVES THE DOMAIN ON THE ',/,&
        7X,'CORRESPONDING OPPOSITE global_here%SIDE OF THE BARRIER ',/,&
        7X,'FLOW RATE AND DIRECTION ARE BASED ON BARRIER ',&
        ' HEIGHT, SURFACE WATER ELEVATION',/,&
        7X,'ON BOTH SIDES OF THE BARRIER, BARRIER COEFFICIENT',&
        ' AND THE APPROPRIATE BARRIER FLOW FORMULA',/,&
        7X,'FREE TANGENTIAL SLIP IS ALLOWED',/)
      ENDIF
      IF(global_here%IBTYPE.EQ.25) THEN
         global_here%NFLUXIB=1
         global_here%NFLUXIBP=1
         WRITE(16,2359)
 2359    FORMAT(5X,'THIS SEGMENT IS AN INTERNAL BARRIER BOUNDARY:',/,&
        7X,'WITH ADDITIONAL global_here%CROSS BARRIER PIPES ',&
        'LOCATED UNDER THE CROWN ',/,&
        7X,'global_here%CROSS BARRIER FLOW IS TREATED AS A NATURAL',&
        ' NORMAL FLOW BOUNDARY CONDITION',/,&
        7X,'WHICH LEAVES/ENTERS THE DOMAIN ON ONE global_here%SIDE OF ',&
        ' THE BARRIER AND ENTERS/LEAVES THE DOMAIN ON THE ',/,&
        7X,'CORRESPONDING OPPOSITE global_here%SIDE OF THE BARRIER ',/,&
        7X,'FLOW RATE AND DIRECTION ARE BASED ON BARRIER ',&
        ' HEIGHT, SURFACE WATER ELEVATION',/,&
        7X,'ON BOTH SIDES OF THE BARRIER, BARRIER COEFFICIENT',&
        ' AND THE APPROPRIATE BARRIER FLOW FORMULA',/,&
        7X,'IN ADDITION global_here%CROSS BARRIER PIPE FLOW RATE AND ',&
        ' DIRECTION ARE BASED ON PIPE CROWN HEIGHT, ',/,&
        7X,'SURFACE WATER ELEVATION ON BOTH SIDES OF THE ',&
        'BARRIER, PIPE FRICTION COEFFICIENT, PIPE DIAMETER',/,&
        7X,' AND THE APPROPRIATE PIPE FLOW FORMULA',/,&
        7X,'FREE TANGENTIAL SLIP IS ALLOWED',/)
      ENDIF
      IF(global_here%IBTYPE.EQ.30) THEN
         global_here%NFLUXRBC=1
         WRITE(16,2355)
 2355    FORMAT(5X,'THIS SEGMENT IS AN OUTWARD RADIATING BOUNDARY:',/,&
        7X,'NORMAL FLUX IS A NATURAL B.C. IN GWCE',/,&
        7X,'NORMAL AND TANGENTIAL VELOCITY ARE COMPUTED FROM ',&
        'THE MOMENTUM EQNS.',/)
      ENDIF


!.....INPUT INFORMATION FOR VARIOUS TYPES OF FLOW BOUNDARY SEGMENTS
!.......INPUT THE STANDARD NODE NUMBERS FOR THE Kth FLOW BOUNDARY SEGMENT
      IF((global_here%IBTYPE.NE.3).AND.(global_here%IBTYPE.NE.13).AND.(global_here%IBTYPE.NE.23).AND.&
     (global_here%IBTYPE.NE.4).AND.(global_here%IBTYPE.NE.24).AND.&
     (global_here%IBTYPE.NE.5).AND.(global_here%IBTYPE.NE.25)) THEN
         DO I=1,global_here%NVELL(K)
            READ(14,*) global_here%NBVV(K,I)
         END DO
         global_here%NPRBI=1
         global_here%NPIPE=0
      ENDIF
!.......INPUT THE NODE NUMBERS FOR THE Kth EXTERNAL BARRIER BOUNDARY SEGMENT
!........ALSO INPUT THE ELEVATION OF THE EXTERNAL BARRIER NODES ABOVE
!........THE GEOID AND THE COEFFICIENT OF FREE SURFACE SUPERCRITICAL
!........FLOW ALONG WITH EACH EXTERNAL BARRIER BOUNDARY NODE FROM UNIT 14
      IF((global_here%IBTYPE.EQ.3).OR.(global_here%IBTYPE.EQ.13).OR.(global_here%IBTYPE.EQ.23)) THEN
         DO I=1,global_here%NVELL(K)
            READ(14,*) global_here%NBVV(K,I),global_here%BARLANHTR(I),global_here%BARLANCFSPR(I)
         ENDDO
         global_here%NPRBI=1
         global_here%NPIPE=0
      ENDIF
!.......INPUT THE NODE NUMBERS FOR THE Kth INTERNAL BARRIER BOUNDARY SEGMENT
!........ALSO INPUT CONNECTION NODE NUMBER AND ELEVATION OF THE INTERNAL BARRIER
!........NODES ABOVE THE GEOID AND THE COEFFICIENTS OF FREE SURFACE SUPERCRITICAL
!........AND SUBCRITICAL FLOW ALONG WITH EACH INTERNAL BARRIER BOUNDARY NODE FROM
!........UNIT 14
      IF((global_here%IBTYPE.EQ.4).OR.(global_here%IBTYPE.EQ.24)) THEN
         DO I=1,global_here%NVELL(K)
            READ(14,*) global_here%NBVV(K,I), global_here%IBCONNR(I), global_here%BARINHTR(I), global_here%BARINCFSBR(I)&
           , global_here%BARINCFSPR(I)
         ENDDO
         global_here%NPRBI=2
         global_here%NPIPE=0
      ENDIF
!.......INPUT THE NODE NUMBERS FOR THE Kth INTERNAL BARRIER BOUNDARY SEGMENT WITH
!........global_here%CROSS BARRIER PIPES; ALSO INPUT CONNECTION NODE NUMBER AND ELEVATION OF THE 
!.......INTERNAL BARRIER NODES ABOVE THE GEOID AND THE COEFFICIENTS OF FREE SURFACE 
!.......SUPERCRITICAL AND SUBCRITICAL FLOW ALONG WITH EACH INTERNAL BARRIER BOUNDARY
!.......NODE FROM UNIT 14; IN ADDITION INPUT THE global_here%CROSS BARRIER PIPE HEIGHT, global_here%CROSS
!.......BARRIER PIPE COEFFICIENT AND global_here%CROSS BARRIER PIPE DIAMETER
      IF((global_here%IBTYPE.EQ.5).OR.(global_here%IBTYPE.EQ.25)) THEN
         DO I=1,global_here%NVELL(K)
            READ(14,*) global_here%NBVV(K,I),global_here%IBCONNR(I),global_here%BARINHTR(I),global_here%BARINCFSBR(I),&
           global_here%BARINCFSPR(I),global_here%PIPEHTR(I),global_here%PIPECOEFR(I),&
           global_here%PIPEDIAMR(I)
         END DO
         global_here%NPRBI=2
         global_here%NPIPE=1
      ENDIF

!.....PROCESS INFORMATION FOR VARIOUS TYPES OF FLOW BOUNDARY SEGMENTS

      DO IPRBI_here=1,global_here%NPRBI

         global_here%IPRBI = IPRBI_here
!.........LOAD PAIRED NODES INTO PRIMARY PROCESSING VECTORS AND RESET
!..........CONNECTING NODES FOR BACK global_here%FACE
!..........THUS BACK/CONNECTING NODES ARE BEING LOADED AS PRIMARY NODES
!..........AND FRONT NODES ARE RELOADED AS CONNECTING NODES
!..........NOTE THAT THE CLOCKWISE ORIENTATION OF ISLAND TYPE BOUNDARIES
!..........IS BEING MAINTAINED WHEN BACK NODES ARE RELOADED AS PRIMARY NODES
!..........ADDITIONAL INTERNAL BARRIER BOUNDARY INFORMATION IS ALSO RESET
         IF(global_here%IPRBI.EQ.2) THEN
            DO I=1,global_here%NVELL(K)
               global_here%NTRAN1(I)=global_here%NBVV(K,I)
               global_here%NTRAN2(I)=global_here%IBCONNR(I)
               global_here%BTRAN3(I)=global_here%BARINHTR(I)
               global_here%BTRAN4(I)=global_here%BARINCFSBR(I)
               global_here%BTRAN5(I)=global_here%BARINCFSPR(I)
               IF(global_here%NPIPE.EQ.1) THEN
                  global_here%BTRAN6(I)=global_here%PIPEHTR(I)
                  global_here%BTRAN7(I)=global_here%PIPECOEFR(I)
                  global_here%BTRAN8(I)=global_here%PIPEDIAMR(I)
               ENDIF
            END DO
            DO I=1,global_here%NVELL(K)
               global_here%NBVV(K,I)=global_here%NTRAN2(global_here%NVELL(K)+1-I)
               global_here%IBCONNR(I)=global_here%NTRAN1(global_here%NVELL(K)+1-I)
               global_here%BARINHTR(I)=global_here%BTRAN3(global_here%NVELL(K)+1-I)
               global_here%BARINCFSBR(I)=global_here%BTRAN4(global_here%NVELL(K)+1-I)
               global_here%BARINCFSPR(I)=global_here%BTRAN5(global_here%NVELL(K)+1-I)
               IF(global_here%NPIPE.EQ.1) THEN
                  global_here%PIPEHTR(I)=global_here%BTRAN6(global_here%NVELL(K)+1-I)
                  global_here%PIPECOEFR(I)=global_here%BTRAN7(global_here%NVELL(K)+1-I)
                  global_here%PIPEDIAMR(I)=global_here%BTRAN8(global_here%NVELL(K)+1-I)
               ENDIF
            ENDDO
         ENDIF

!.........WRITE OUT ADDITIONAL HEADER FOR INTERNAL BARRIER BOUNDARIES

         IF((global_here%IBTYPE.EQ.4).OR.(global_here%IBTYPE.EQ.24)) THEN
            IF(global_here%IPRBI.EQ.1) THEN
               WRITE(16,1842)
 1842          FORMAT(/,5X,'FRONT global_here%FACE OF INTERNAL BARRIER BOUNDARY',/)
            ELSE
               WRITE(16,1843)
 1843          FORMAT(/,5X,'BACK global_here%FACE OF INTERNAL BARRIER BOUNDARY',/)
            ENDIF
         ENDIF

!.........WRITE OUT ADDITIONAL HEADER FOR INTERNAL BARRIER BOUNDARIES 
!..........WITH global_here%CROSS BARRIER PIPES

         IF((global_here%IBTYPE.EQ.5).OR.(global_here%IBTYPE.EQ.25)) THEN
            IF(global_here%IPRBI.EQ.1) THEN
               WRITE(16,1844)
 1844          FORMAT(/,5X,'FRONT global_here%FACE OF INTERNAL BARRIER BOUNDARY',&
              ' WITH global_here%CROSS BARRIER PIPES',/)
            ELSE
               WRITE(16,1845)
 1845          FORMAT(/,5X,'BACK global_here%FACE OF INTERNAL BARRIER BOUNDARY',&
              ' WITH global_here%CROSS BARRIER PIPES',/)
            ENDIF
         ENDIF

!.........WRITE OUT GENERAL HEADER FOR BOUNDARY INFORMATION

         WRITE(16,1841)
 1841    FORMAT('    global_here%JGW    global_here%JME    global_here%ME2GW   NODE #  BNDRY CODE   INNER',&
        ' ANGLE',7X,'COS',13X,'SIN',9X,'0.667*BNDRY LEN',/)

!.........COMPLETE THE BOUNDARY ARRAY FOR THE Kth FLOW BOUNDARY SEGMENT

         global_here%NBVV(K,0)=global_here%NBVV(K,1)    !UNCLOSED EXTERNAL
         IF((global_here%IBTYPE.EQ.1).OR.(global_here%IBTYPE.EQ.11).OR.(global_here%IBTYPE.EQ.21)) THEN
            IF(global_here%NBVV(K,global_here%NVELL(K)).NE.global_here%NBVV(K,1)) THEN !CLOSE AN UNCLOSED INTERNAL
               global_here%NVELL(K)=global_here%NVELL(K)+1
               global_here%NBVV(K,global_here%NVELL(K))=global_here%NBVV(K,1)
            ENDIF
         ENDIF
         IF(global_here%NBVV(K,global_here%NVELL(K)).EQ.global_here%NBVV(K,1)) THEN !CLOSED EXTERNAL OR INTERNAL
            global_here%NBVV(K,0)=global_here%NBVV(K,global_here%NVELL(K)-1)
         ENDIF
         global_here%NBVV(K,global_here%NVELL(K)+1)=global_here%NBVV(K,global_here%NVELL(K))

!.........PUT BOUNDARY INFORMATION INTO 2 TYPES OF ARRAYS, ONE FOR THE GWCE B.C.
!..........AND ONE FOR THE MOMENTUM EQUATION B.C.
!..........THE GWCE ARRAYS INCLUDE EVERY NODE IN THE UNIT 14 FILE, I.E., NODES
!..........ARE REPEATED WHERE SPECIFIED NORMAL FLOW AND NO NORMAL FLOW BOUNDARIES
!..........MEET AND AT THE BEGINNING AND END OF CLOSED EXTERNAL BOUNDARIES AND
!.........ISLANDS.
!..........THE MOMENTUM EQUATION ARRAYS ARE KEYED TO THE GWCE ARRAYS VIA THE
!..........ARRAY global_here%ME2GW WHICH INDICATES THE LOCATION IN THE GWCE ARRAYS THAT
!..........THE APPROPRIATE M.E. VALUE LIES.
!..........THE M.E. ARRAYS DO NOT REPEAT NODES THAT ARE DUPLICATED IN THE
!..........UNIT 14 FILE, I.E., WHEN SPECIFIED NORMAL FLOW AND NO NORMAL FLOW
!..........BOUNDARIES MEET, THE SPECIFIED NORMAL FLOW BOUNDARY CONDITION TAKES
!..........PRECEDENT.  ALSO THE BEGINNING AND ENDING NODES OF CLOSED EXTERNAL
!..........AND ISLAND BOUNDARIES ARE NOT REPEATED.

         DO I=1,global_here%NVELL(K)

!.........SET UP THE GWCE BOUNDARY ARRAYS WHICH CONSIST OF
!..........BOUNDARY NODE NUMBERS
!..........BOUNDARY CODES
!..........0.66667*LENGTH OF EACH BOUNDARY SEGMENT.  NOTE, THE LENGTH OF THE LAST
!..........BOUNDARY SEGMENT ON EACH BOUNDARY SHOULD BE ZERO

            global_here%JGW=global_here%JGW+1
            IF((global_here%IBTYPE.EQ.0).OR.(global_here%IBTYPE.EQ.10).OR.(global_here%IBTYPE.EQ.20).OR.&
           (global_here%IBTYPE.EQ.2).OR.(global_here%IBTYPE.EQ.12).OR.(global_here%IBTYPE.EQ.22).OR.&
           (global_here%IBTYPE.EQ.3).OR.(global_here%IBTYPE.EQ.13).OR.(global_here%IBTYPE.EQ.23).OR.&
           (global_here%IBTYPE.EQ.30)) THEN
               global_here%NVELEXT=global_here%NVELEXT+1
            ENDIF
            global_here%NBV(global_here%JGW)=global_here%NBVV(K,I)
            global_here%NBVI=global_here%NBVV(K,I)
            global_here%NBVJ=global_here%NBVV(K,I+1)
            global_here%DELX=global_here%X(global_here%NBVJ)-global_here%X(global_here%NBVI)
            global_here%DELY=global_here%Y(global_here%NBVJ)-global_here%Y(global_here%NBVI)
            global_here%BNDLEN2O3(global_here%JGW)=2.D0*(SQRT(global_here%DELX*global_here%DELX+global_here%DELY*global_here%DELY))/3.D0

!...........COMPUTE THE INCLUDED ANGLE AND TEST TO DETERMINE WHETHER TO ZERO
!............TANGENTIAL VELOCITIES
!...........NOTE:.IMPLEMENTATION FOR global_here%ICS=2 REQUIRES COMPUTING ALL COORDINATES IN
!............A LOCALIZED SYSTEM (I.E. THE TRANSFORMATION IS CENTERED AT X0,Y0)

            IF(global_here%ICS.EQ.1) THEN
               global_here%XL0=global_here%X(global_here%NBVV(K,I))
               global_here%XL1=global_here%X(global_here%NBVV(K,I-1))
               global_here%XL2=global_here%X(global_here%NBVV(K,I+1))
               global_here%YL0=global_here%Y(global_here%NBVV(K,I))
               global_here%YL1=global_here%Y(global_here%NBVV(K,I-1))
               global_here%YL2=global_here%Y(global_here%NBVV(K,I+1))
            ELSE
               CALL CPP(global_here%XL0,global_here%YL0,global_here%SLAM(global_here%NBVV(K,I)),global_here%SFEA(global_here%NBVV(K,I)),&
              global_here%SLAM(global_here%NBVV(K,I)),global_here%SFEA(global_here%NBVV(K,I)))
               CALL CPP(global_here%XL1,global_here%YL1,global_here%SLAM(global_here%NBVV(K,I-1)),global_here%SFEA(global_here%NBVV(K,I-1)),&
              global_here%SLAM(global_here%NBVV(K,I)),global_here%SFEA(global_here%NBVV(K,I)))
               CALL CPP(global_here%XL2,global_here%YL2,global_here%SLAM(global_here%NBVV(K,I+1)),global_here%SFEA(global_here%NBVV(K,I+1)),&
              global_here%SLAM(global_here%NBVV(K,I)),global_here%SFEA(global_here%NBVV(K,I)))
            ENDIF

!...........NOTE: INTERIOR ANGLE AT ENDS OF BOUNDARIES MUST BE EQUAL, EITHER:
!............A FICTICIOUSLY LARGE VALUE IF THE BOUNDARY IS NOT CLOSED OR
!............A TRUE VALUE IF THE BOUNDARY IS CLOSED

            global_here%THETA=0.
            IF((I.EQ.1).AND.(global_here%NBVV(K,I).EQ.global_here%NBVV(K,I-1))) THEN
               global_here%THETA1=-9999999.d0
               global_here%THETA=global_here%THETA1
               global_here%COSTHETA1=global_here%COSTSET
               global_here%COSTHETA=global_here%COSTHETA1
               global_here%CROSS1=0.d0
               global_here%CROSS=global_here%CROSS1
            ENDIF
            IF(I.EQ.global_here%NVELL(K)) THEN
               global_here%THETA=global_here%THETA1
               global_here%COSTHETA=global_here%COSTHETA1
               global_here%CROSS=global_here%CROSS1
            ENDIF
            IF(global_here%THETA.EQ.0.) THEN
               global_here%VL1X=global_here%XL1-global_here%XL0
               global_here%VL1Y=global_here%YL1-global_here%YL0
               global_here%VL2X=global_here%XL2-global_here%XL0
               global_here%VL2Y=global_here%YL2-global_here%YL0
               global_here%DOTVEC=global_here%VL1X*global_here%VL2X+global_here%VL1Y*global_here%VL2Y
               global_here%VECNORM=(SQRT(global_here%VL1X**2+global_here%VL1Y**2))*(SQRT(global_here%VL2X**2+global_here%VL2Y**2))
               global_here%COSTHETA=global_here%DOTVEC/global_here%VECNORM
               IF(global_here%COSTHETA.GT.1.0d0) global_here%COSTHETA=1.0d0
               IF(global_here%COSTHETA.LT.-1.0d0) global_here%COSTHETA=-1.0d0
               global_here%THETA=RAD2DEG*ACOS(global_here%COSTHETA)
               global_here%CROSS=-global_here%VL1X*global_here%VL2Y+global_here%VL2X*global_here%VL1Y
               IF(global_here%CROSS.LT.0) global_here%THETA=360.d0-global_here%THETA
               IF(I.EQ.1) THEN
                  global_here%THETA1=global_here%THETA
                  global_here%COSTHETA1=global_here%COSTHETA
                  global_here%CROSS1=global_here%CROSS
               ENDIF
            ENDIF

!...........CHECK WHETHER ANGLE IS LESS THAN MINIMUM ANGLE, IF SO CHANGE THE
!............BOUNDARY CODE TO ZERO TANGENTIAL VELOCITIES

            global_here%LBCODEI(global_here%JGW)=global_here%IBTYPE
            IF((global_here%COSTHETA.GT.global_here%COSTSET).AND.(global_here%CROSS.GT.0.0)) THEN
               IF(global_here%IBTYPE.EQ.0) global_here%LBCODEI(global_here%JGW)=10
               IF(global_here%IBTYPE.EQ.1) global_here%LBCODEI(global_here%JGW)=11
               IF(global_here%IBTYPE.EQ.2) global_here%LBCODEI(global_here%JGW)=12
               IF(global_here%IBTYPE.EQ.3) global_here%LBCODEI(global_here%JGW)=13
               IF((global_here%IBTYPE.GE.0).AND.(global_here%IBTYPE.LE.3)) THEN
                  WRITE(16,1856) global_here%NBVV(K,I),global_here%THETA
 1856             FORMAT(2X,I7,4X,'THE INNER ANGLE = ',F8.2,1X,&
                 'TANGENTIAL SLIP WILL BE ZEROED')
               ENDIF
            ENDIF

!...........COMPUTE COS AND SIN OF OUTWARD NORMAL REGARDLESS OF BOUNDARY TYPE

            global_here%X1=global_here%X(global_here%NBVV(K,I-1))
            global_here%X2=global_here%X(global_here%NBVV(K,I+1))
            global_here%Y1=global_here%Y(global_here%NBVV(K,I-1))
            global_here%Y2=global_here%Y(global_here%NBVV(K,I+1))
            global_here%XL=SQRT((global_here%X1-global_here%X2)**2+(global_here%Y1-global_here%Y2)**2)
            global_here%CSII(global_here%JGW)=global_here%SFAC(global_here%NBVV(K,I))*(global_here%Y2-global_here%Y1)/global_here%XL
            global_here%SIII(global_here%JGW)=(global_here%X1-global_here%X2)/global_here%XL

!...........SET UP THE MOMENTUM EQUATION BOUNDARY ARRAY WHICH CONSISTS OF
!............A KEY TO THE GWCE BOUNDARY CONDITION ARRAY

            IF(I.EQ.1) THEN     !DEAL WITH FIRST NODE IN L.B. SEG
               IF(global_here%JGW.EQ.1) THEN !VERY FIRST L.B. SEG
                  global_here%JME=global_here%JME+1     !M.E. USES IT
                  global_here%ME2GW(global_here%JME)=global_here%JGW
               ENDIF
               IF(global_here%JGW.NE.1) THEN
                  IF(global_here%NBV(global_here%JGW).NE.global_here%NBV(global_here%JGW-1)) THEN !L.B. SEGS DON'T OVERLAP
                     global_here%JME=global_here%JME+1  !M.E. USES IT
                     global_here%ME2GW(global_here%JME)=global_here%JGW
                  ENDIF
                  IF(global_here%NBV(global_here%JGW).EQ.global_here%NBV(global_here%JGW-1)) THEN !L.B. SEGS OVERLAP
                     IF((global_here%LBCODEI(global_here%JGW).EQ.2) .OR. &! M.E. USES IT ONLY&
                    (global_here%LBCODEI(global_here%JGW).EQ.12).OR. &! IF IT IS&
                    (global_here%LBCODEI(global_here%JGW).EQ.22).OR. &! SPECIFIED FLOW,&
                    (global_here%LBCODEI(global_here%JGW).EQ.3) .OR. &! AN OVERFLOW BARRIER&
                    (global_here%LBCODEI(global_here%JGW).EQ.13).OR. &! OR A RADIATION&
                    (global_here%LBCODEI(global_here%JGW).EQ.23).OR. &! BOUNDARY&
                    (global_here%LBCODEI(global_here%JGW).EQ.30)) global_here%ME2GW(global_here%JME)=global_here%JGW
                  ENDIF                                                                
               ENDIF
            ENDIF
            IF((I.GT.1).AND.(I.LT.global_here%NVELL(K))) THEN !IF NOT FIRST OR
               global_here%JME=global_here%JME+1        !LAST NODE
               global_here%ME2GW(global_here%JME)=global_here%JGW   !M.E. USES IT
            ENDIF
            IF(I.EQ.global_here%NVELL(K)) THEN !DEAL WITH LAST NODE ON BOUNDARY
               IF((global_here%NBV(global_here%JGW).NE.global_here%NBVV(K,1)).AND.& !IF UNCLOSED BOUNDARY&
              (global_here%NBV(global_here%JGW).NE.global_here%NBV(1))) THEN !M.E. USES IT
                  global_here%JME=global_here%JME+1
                  global_here%ME2GW(global_here%JME)=global_here%JGW
               ENDIF
               IF(global_here%NBVV(K,I).EQ.global_here%NBV(1)) THEN !IF OVERLAPS WITH VERY FIRST
                  IF((global_here%LBCODEI(global_here%JGW).EQ.2) .OR.& ! L.B. NODE&
                 (global_here%LBCODEI(global_here%JGW).EQ.12).OR.& ! M.E. USES IT ONLY IF IT IS&
                 (global_here%LBCODEI(global_here%JGW).EQ.22).OR.& ! SPECIFIED FLOW,&
                 (global_here%LBCODEI(global_here%JGW).EQ.3) .OR.& ! AN OVERFLOW BARRIER OR&
                 (global_here%LBCODEI(global_here%JGW).EQ.13).OR.& ! A RADIATION&
                 (global_here%LBCODEI(global_here%JGW).EQ.23).OR.& ! BOUNDARY&
                 (global_here%LBCODEI(global_here%JGW).EQ.30)) global_here%ME2GW(1)=global_here%JGW
               ENDIF
            ENDIF

!...........LOAD EXTERNAL BARRIER BOUNDARY INFORMATION INTO THE CORRECT VECTORS
            IF((global_here%IBTYPE.EQ.3).OR.(global_here%IBTYPE.EQ.13).OR.(global_here%IBTYPE.EQ.23)) THEN
               global_here%BARLANHT(global_here%JGW)=global_here%BARLANHTR(I)
               global_here%BARLANCFSP(global_here%JGW)=global_here%BARLANCFSPR(I)
            ENDIF

!...........LOAD INTERNAL BARRIER BOUNDARY INFORMATION INTO THE CORRECT VECTORS
            IF((global_here%IBTYPE.EQ.4).OR.(global_here%IBTYPE.EQ.24)) THEN
               global_here%IBCONN(global_here%JGW)=global_here%IBCONNR(I)
               global_here%BARINHT(global_here%JGW)=global_here%BARINHTR(I)
               global_here%BARINCFSB(global_here%JGW)=global_here%BARINCFSBR(I)
               global_here%BARINCFSP(global_here%JGW)=global_here%BARINCFSPR(I)
            ENDIF

!...........LOAD INTERNAL BARRIER WITH PIPES BOUNDARY INFORMATION INTO 
!............THE CORRECT VECTORS
            IF((global_here%IBTYPE.EQ.5).OR.(global_here%IBTYPE.EQ.25)) THEN
               global_here%IBCONN(global_here%JGW)=global_here%IBCONNR(I)
               global_here%BARINHT(global_here%JGW)=global_here%BARINHTR(I)
               global_here%BARINCFSB(global_here%JGW)=global_here%BARINCFSBR(I)
               global_here%BARINCFSP(global_here%JGW)=global_here%BARINCFSPR(I)
               global_here%PIPEHT(global_here%JGW)=global_here%PIPEHTR(I)
               global_here%PIPECOEF(global_here%JGW)=global_here%PIPECOEFR(I)
               global_here%PIPEDIAM(global_here%JGW)=global_here%PIPEDIAMR(I)
            ENDIF
            
!...........WRITE OUT BOUNDARY CONDITION ARRAY INFORMATION

            WRITE(16,1857) global_here%JGW,global_here%JME,global_here%ME2GW(global_here%JME),global_here%NBV(global_here%JGW),global_here%LBCODEI(global_here%JGW),&
           global_here%THETA,global_here%CSII(global_here%JGW),global_here%SIII(global_here%JGW),global_here%BNDLEN2O3(global_here%JGW)
 1857       FORMAT(1X,I6,1X,I6,1X,I6,3X,I6,3X,I4,9X,F8.2,2X,E16.8,1X,&
           E16.8,2X,E16.8)

!...........CHECK EXTERNAL BARRIER HEIGHTS AGAINST DEPTHS
            IF((global_here%IBTYPE.EQ.3).OR.(global_here%IBTYPE.EQ.13).OR.(global_here%IBTYPE.EQ.23)) THEN
               IF(global_here%BARLANHT(global_here%JGW).LT.-global_here%DP(global_here%NBV(global_here%JGW))) THEN
                  IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,8367) &
                 global_here%JGW,global_here%NBV(global_here%JGW),global_here%BARLANHT(global_here%JGW),global_here%DP(global_here%NBV(global_here%JGW))
                  WRITE(16,8367) global_here%JGW,global_here%NBV(global_here%JGW),global_here%BARLANHT(global_here%JGW),global_here%DP(global_here%NBV(global_here%JGW))
 8367             FORMAT(////,1X,'!!!!!!!!!!  FATAL INPUT ERROR   !!!'&
                 ,'!!!!!!!!!',//,&
                 1X,'AT BOUNDARY NODE NO.',I6,' (GLOBAL NODE NO.',&
                 I6, ' AND OF EXTERNAL BARRIER TYPE) ',/,&
                 2X,'THE EXTERNAL BARRIER HEIGHT = ',E12.5,&
                 2X,'IS EXCEEDED BY THE DEPTH SPECIFIED AT ',/,2X&
                 ,'THE ASSOCIATED GLOBAL NODE = ',E12.5,/,2X,&
                 'USER MUST SPECIFY CONSISTENT BARRIER HEIGHTS',&
                 ' AND DEPTHS')
                  IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9973)
                  WRITE(16,9973)
                  STOP
               ENDIF
            ENDIF
!...........CHECK INTERNAL BARRIER HEIGHTS AGAINST DEPTHS
            IF((global_here%IBTYPE.EQ.4).OR.(global_here%IBTYPE.EQ.24)) THEN
               IF(global_here%BARINHT(global_here%JGW).LT.-global_here%DP(global_here%NBV(global_here%JGW))) THEN
                  IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,8368) &
                 global_here%JGW,global_here%NBV(global_here%JGW),global_here%BARINHT(global_here%JGW),global_here%DP(global_here%NBV(global_here%JGW))
                  WRITE(16,8368) global_here%JGW,global_here%NBV(global_here%JGW),global_here%BARINHT(global_here%JGW),global_here%DP(global_here%NBV(global_here%JGW))
 8368             FORMAT(////,1X,'!!!!!!!!!!  FATAL INPUT ERROR   !!!'&
                 ,'!!!!!!!!!',//,&
                 1X,'AT BOUNDARY NODE NO.',I6,' (GLOBAL NODE NO. ',&
                 I6,' AND OF INTERNAL BARRIER TYPE) ',/,&
                 2X,'THE INTERNAL BARRIER HEIGHT = ',E12.5,&
                 2X,'IS EXCEEDED BY THE DEPTH SPECIFIED AT ',/,2X&
                 ,'THE ASSOCIATED GLOBAL NODE = ',E12.5,/,2X,&
                 'USER MUST SPECIFY CONSISTENT BARRIER HEIGHTS',&
                 ' AND DEPTHS')
                  IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9973)
                  WRITE(16,9973)
                  STOP
               ENDIF
            ENDIF
            
!...........CHECK INTERNAL BARRIER WITH PIPES BARRIER HEIGHTS AGAINST DEPTHS
            IF((global_here%IBTYPE.EQ.5).OR.(global_here%IBTYPE.EQ.25)) THEN
               IF(global_here%BARINHT(global_here%JGW).LT.-global_here%DP(global_here%NBV(global_here%JGW))) THEN
                  IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,8370) &
                 global_here%JGW,global_here%NBV(global_here%JGW),global_here%BARINHT(global_here%JGW),global_here%DP(global_here%NBV(global_here%JGW))
                  WRITE(16,8370) global_here%JGW,global_here%NBV(global_here%JGW),global_here%BARINHT(global_here%JGW),global_here%DP(global_here%NBV(global_here%JGW))
 8370             FORMAT(////,1X,'!!!!!!!!!!  FATAL INPUT ERROR   !!!'&
                 ,'!!!!!!!!!',//,&
                 1X,'AT BOUNDARY NODE NO.',I6,' (GLOBAL NODE NO. ',&
                 I6,' AND OF INTERNAL BARRIER TYPE) ',/,&
                 2X,'THE INTERNAL BARRIER HEIGHT = ',E12.5,&
                 2X,'IS EXCEEDED BY THE DEPTH SPECIFIED AT ',/,2X&
                 ,'THE ASSOCIATED GLOBAL NODE = ',E12.5,/,2X,&
                 'USER MUST SPECIFY CONSISTENT BARRIER HEIGHTS',&
                 ' AND DEPTHS')
                  IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9973)
                  WRITE(16,9973)
                  STOP
               ENDIF
            ENDIF
!...........CHECK INTERNAL BARRIER WITH PIPES PIPE HEIGHTS AGAINST DEPTHS
            IF((global_here%IBTYPE.EQ.5).OR.(global_here%IBTYPE.EQ.25)) THEN
               IF(global_here%PIPEHT(global_here%JGW).LT.-global_here%DP(global_here%NBV(global_here%JGW))) THEN
                  IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,8372) &
                 global_here%JGW,global_here%NBV(global_here%JGW),global_here%BARINHT(global_here%JGW),global_here%DP(global_here%NBV(global_here%JGW))
                  WRITE(16,8372) global_here%JGW,global_here%NBV(global_here%JGW),global_here%BARINHT(global_here%JGW),global_here%DP(global_here%NBV(global_here%JGW))
 8372             FORMAT(////,1X,'!!!!!!!!!!  FATAL INPUT ERROR   !!!'&
                 ,'!!!!!!!!!',//,&
                 1X,'AT BOUNDARY NODE NO.',I6,' (GLOBAL NODE NO. ',&
                 I6,' AND OF INTERNAL BARRIER TYPE) ',/,&
                 2X,'THE BARRIER PIPE HEIGHT = ',E12.5,&
                 2X,'IS EXCEEDED BY THE DEPTH SPECIFIED AT ',/,2X&
                 ,'THE ASSOCIATED GLOBAL NODE = ',E12.5,/,2X,&
                 'USER MUST SPECIFY CONSISTENT PIPE HEIGHTS',&
                 ' AND DEPTHS')
                  IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9973)
                  WRITE(16,9973)
                  STOP
               ENDIF
            ENDIF

!...........CHECK FOR OVERLAPPING OF AN INTERNAL BARRIER BOUNDARY WITH
!............ANY EXTERNAL BARRIER BOUNDARY. IF THIS DOES OCCUR, TAKE
!............APPROPRIATE ACTION
            IF((global_here%IBTYPE.EQ.4).OR.(global_here%IBTYPE.EQ.24).OR.(global_here%IBTYPE.EQ.5)&
           .OR.(global_here%IBTYPE.EQ.25)) THEN
               DO ICK_here=1,global_here%NVELEXT
                  global_here%ICK=ICK_here
!...............CHECK IF OVERLAP EXISTS
                  IF(global_here%NBV(global_here%ICK).EQ.global_here%NBV(global_here%JGW)) THEN
!.................CHECK FOR ILLEGAL OVERLAPS
                     IF((global_here%LBCODEI(global_here%ICK).EQ.2).OR.(global_here%LBCODEI(global_here%ICK).EQ.3).OR.&
                    (global_here%LBCODEI(global_here%ICK).EQ.12).OR.(global_here%LBCODEI(global_here%ICK).EQ.13)) THEN 
                        IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,8567) &
                       global_here%JGW,global_here%NBV(global_here%JGW),global_here%ICK,global_here%NBV(global_here%ICK)
                        WRITE(16,8567) global_here%JGW,global_here%NBV(global_here%JGW),global_here%ICK,global_here%NBV(global_here%ICK)
 8567                   FORMAT(////,1X,'!!!!!!!!!!  FATAL INPUT ERROR   !!!'&
                       ,'!!!!!!!!!',//,&
                       1X,'BOUNDARY NODE NO. ',I6,' (GLOBAL NODE NO. ',&
                       I9, 'AND OF INTERNAL BARRIER TYPE) ',/,&
                       2X,'OVERLAPS BOUNDARY NODE NO.',I6,' (GLOBAL NODE'&
                       ,' NO.',I6,' )',/,&
                       2X,'THIS IS AN ILLEGAL TYPE OVERLAP !! - INTERNAL '&
                       ,'BARRIER BOUNDARIES CAN ONLY OVERLAP WITH ',&
                       'NO NORMAL FLOW EXTERNAL BOUNDARIES',/&
                       2X,'(I.E. global_here%IBTYPE=0,10,20)')
                        IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9973)
                        WRITE(16,9973)
                        STOP
                     ENDIF
!.................CHECK FOR OVERLAPS WHICH REQUIRE ADJUSTMENTS OF BOUNDARY
!..................CODE ON THE EXTERNAL BOUNDARY
                     IF(((global_here%IBTYPE.EQ.4).AND.(global_here%LBCODEI(global_here%ICK).EQ.0))&
                    .OR.((global_here%IBTYPE.EQ.5).AND.(global_here%LBCODEI(global_here%ICK).EQ.0))) THEN
                        WRITE(16,8568) global_here%JGW,global_here%ICK,global_here%ICK
 8568                   FORMAT(1X,'DUE TO LEGAL OVERLAPPING OF ',&
                       'BOUNDARY NODE',I7,' (WHICH IS AN ESSENTIAL INTER'&
                       ,'NAL BARRIER BOUNDARY NODE)', /,2X,&
                       'AND BOUNDARY NODE',I7,' (WHICH IS AN ESSENTIAL ',&
                       'EXTERNAL NO NORMAL FLOW WITH SLIP BOUNDARY',&
                       ' NODE),',/,2X,&
                       'THE BOUNDARY TYPE FOR BOUNDARY NODE ',I7,&
                       ' IS BEING RESET TO global_here%IBTYPE=20',/,2X,&
                       '(NATURAL NO NORMAL FLOW WITH SLIP BOUNDARY) ')
                        global_here%LBCODEI(global_here%ICK)=20
                     ENDIF
                     IF(((global_here%IBTYPE.EQ.4).AND.(global_here%LBCODEI(global_here%ICK).EQ.10)) &
                    .OR.((global_here%IBTYPE.EQ.5).AND.(global_here%LBCODEI(global_here%ICK).EQ.10))) THEN
                        WRITE(16,8569) global_here%JGW,global_here%ICK,global_here%ICK
 8569                   FORMAT(1X,'DUE TO LEGAL OVERLAPPING OF ',&
                       'BOUNDARY NODE ',I7,' (WHICH IS AN ESSENTIAL INTER'&
                       ,'NAL BARRIER BOUNDARY NODE)', /,2X,&
                       'AND BOUNDARY NODE',I7,' (WHICH IS AN ESSENTIAL ',&
                       'EXTERNAL NO NORMAL FLOW WITH NO SLIP BOUNDARY',&
                       ' NODE),',/,2X,&
                       'THE BOUNDARY TYPE FOR BOUNDARY NODE ',I7,&
                       ' IS BEING RESET TO global_here%IBTYPE=20',/,2X,&
                       '(NATURAL NO NORMAL FLOW WITH SLIP BOUNDARY) ')
                        global_here%LBCODEI(global_here%ICK)=20
                     ENDIF
                     IF(((global_here%IBTYPE.EQ.24).AND.(global_here%LBCODEI(global_here%ICK).EQ.10))  &
                    .OR.((global_here%IBTYPE.EQ.25).AND.(global_here%LBCODEI(global_here%ICK).EQ.10))) THEN
                        WRITE(16,8570) global_here%JGW,global_here%ICK,global_here%ICK
 8570                   FORMAT(1X,'DUE TO LEGAL OVERLAPPING OF ',&
                       'BOUNDARY NODE',I7,' (WHICH IS A NATURAL INTERNAL'&
                       ,' BARRIER BOUNDARY NODE)', /,2X,&
                       'AND BOUNDARY NODE',I7,' (WHICH IS AN ESSENTIAL ',&
                       'EXTERNAL NO NORMAL FLOW WITH NO SLIP BOUNDARY',&
                       ' NODE),',/,2X,&
                       'THE BOUNDARY TYPE FOR BOUNDARY NODE',I7,&
                       ' IS BEING RESET TO global_here%IBTYPE=0',/,2X,&
                       '(ESSENTIAL NO NORMAL FLOW WITH SLIP BOUNDARY) ')
                        global_here%LBCODEI(global_here%ICK)=0
                     ENDIF
                  ENDIF
               END DO
            ENDIF

         ENDDO
      ENDDO
      
!.......Put barrier data into arrays more amiable to DG data structure

      IF ((global_here%IBTYPE.EQ.3).OR.(global_here%IBTYPE.EQ.13).OR.(global_here%IBTYPE.EQ.23)) THEN
         DO I = 1,global_here%NVELL(K)-1
            dg_here%NEBSEG = dg_here%NEBSEG + 1
            dg_here%EBHT(dg_here%NEBSEG)   = 0.5D0*( global_here%BARLANHTR(I)   + global_here%BARLANHTR(I+1)   )
            dg_here%EBCFSP(dg_here%NEBSEG) = 0.5D0*( global_here%BARLANCFSPR(I) + global_here%BARLANCFSPR(I+1) )
         ENDDO
      ENDIF
      
      IF ((global_here%IBTYPE.EQ.4).OR.(global_here%IBTYPE.EQ.24)) THEN
         DO I = 1,global_here%NVELL(K)-1
            dg_here%NIBSEG = dg_here%NIBSEG + 1
            dg_here%IBHT(dg_here%NIBSEG)   = 0.5D0*( global_here%BARINHTR(I)   + global_here%BARINHTR(I+1)   )
            dg_here%IBCFSB(dg_here%NIBSEG) = 0.5D0*( global_here%BARINCFSBR(I) + global_here%BARINCFSBR(I+1) )
            dg_here%IBCFSP(dg_here%NIBSEG) = 0.5D0*( global_here%BARINCFSPR(I) + global_here%BARINCFSPR(I+1) )
            dg_here%BACKNODES(1,dg_here%NIBSEG) = global_here%IBCONNR(I)
            dg_here%BACKNODES(2,dg_here%NIBSEG) = global_here%IBCONNR(I+1)
         ENDDO
      ENDIF
      ENDDO

!.....ONCE ALL FLOW BOUNDARY NODES HAVE BEEN PROCESSED, CHECK TO MAKE SURE
!.....THAT global_here%JGW LE MNVEL.  NOTE, global_here%JME MUST BE < global_here%JGW.

      IF(s%MNVEL.LT.global_here%JGW) THEN
         IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9947)
         WRITE(16,9947)
 9947    FORMAT(////,1X,'!!!!!!!!!!  FATAL INPUT ERROR   !!!!!!!!!!!!',&
        //,1X,'THE DIMENSION PARAMETER MNVEL IS LESS THAN ',&
        'THE TOTAL NUMBER OF FLOW BOUNDARY NODES',&
        /,1X,'FROM ALL THE SPECIFIED FLOW SEGMENTS COMBINED',/)
         IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9973)
         WRITE(16,9973)
         STOP
      ENDIF

      global_here%NVEL=global_here%JGW
      global_here%NVELME=global_here%JME

!     
      DO IK=1,s%MNP               ! FINISH DET. MAX # NEIGHBORS
         IF(global_here%NNEIGH(IK).GT.s%MNEI) s%MNEI=global_here%NNEIGH(IK)
      ENDDO
      s%MNEI = s%MNEI+1

!     sb   Copied from read_input.dg%F of v45.01
!     estimate the maximum array space needed for the neighbor
!     table by
!     increasing this number by 2, to provide array space for the
!     node itself
!     and in case the maximum number of nodes occurs at a boundary
!     node

      S%MNEI = S%MNEI+3


!.....TRANSFER FLOW BOUNDARY INFORMATION INTO NODAL ARRAYS

      DO I=1,global_here%NP
         global_here%LBCODE(I)=-1
         global_here%CSI(I)=0.
         global_here%SII(I)=1.
      END DO

      DO I=1,global_here%NVELME
         J=global_here%ME2GW(I)
         global_here%LBCODE(global_here%NBV(J))=global_here%LBCODEI(J)
         global_here%CSI(global_here%NBV(J))=global_here%CSII(J)
         global_here%SII(global_here%NBV(J))=global_here%SIII(J)
      END DO

!...  IF ANY NON ZERO NORMAL FLOW BOUNDARIES WERE SPECIFIED, (global_here%NFLUXF=1)
!.....READ FORCING INFORMATION FROM UNIT 15 FILE

      global_here%NFFR = 0
      IF(global_here%NFLUXF.EQ.1) THEN

!.....INPUT FROM THE NUMBER OF FREQUENCIES PRESENT IN NORMAL FLOW FORCING
!......DATA.  IF THIS = 0, NORMAL FLOW DATA IS READ IN FROM THE FORT.20 FILE.

         READ(15,*) global_here%NFFR
         s%MNFFR = global_here%NFFR
         IF (global_here%NFFR.EQ.0) s%MNFFR = 1

!.....Allocate space for periodic normal flow boundary conditions
         call alloc_main6(s,global_here)
!     
         DO I=1,global_here%NVELME
            J=global_here%ME2GW(I)
            global_here%QNAM(1,J)=0.
            global_here%QNPH(1,J)=0.
         ENDDO

!.....READ IN AND WRITE OUT INFO ON SPECIFIED NORMAL FLOW BOUNDARIES

         WRITE(16,1112)
         WRITE(16,2200)
 2200    FORMAT(//,1X,'NORMAL FLOW BOUNDARY FORCING INFORMATION ',//)
         IF(global_here%NFFR.EQ.0) THEN
            WRITE(16,2201)
 2201       FORMAT(/,5X,'NORMAL FLOW VALUES WILL BE READ FROM UNIT 20 ',&
      /,9X,'INTERPOLATION IN TIME IS DONE TO SYNC THE FLOW DATA ',&
      /,9X,'WITH THE MODEL TIME STEP.')
         ENDIF
         IF(global_here%NFFR.NE.0) THEN
            WRITE(16,2202) global_here%NFFR
 2202       FORMAT(/,5X,'NUMBER OF PERIODIC NORMAL FLOW CONSTITUENTS =',&
           I5)
            WRITE(16,2203)
 2203       FORMAT(/,7X,'CONSTITUENT #',4X,'FREQUENCY',4X,'NODAL FACTOR',&
           3X,'EQU.global_here%ARG (DEG)',2X,'CONSTITUENT',/)
            DO I=1,global_here%NFFR
               READ(15,'(A5)') global_here%FBOUNTAG(I)
               READ(15,*) global_here%FAMIG(I),global_here%FFF(I),global_here%FFACE(I)
               WRITE(16,2204) I,global_here%FAMIG(I),global_here%FFF(I),global_here%FFACE(I),global_here%FBOUNTAG(I)
 2204          FORMAT(12X,I2,6X,F16.12,2X,F10.7,2X,F10.3,10X,A5)
               global_here%FFACE(I)=global_here%FFACE(I)*DEG2RAD
               IF(global_here%FAMIG(I).EQ.0.) THEN
                  global_here%FPER(I)=0.
               ELSE
                  global_here%FPER(I)=2.D0*PI/global_here%FAMIG(I)
               ENDIF
            END DO

!.......INPUT PERIODIC NORMAL FLOW FORCING CONDITIONS ON DESIGNATED FLOW BOUNDARIES
!........FOR EACH OF THE FORCING FREQUENCIES FROM UNIT 15 AND OUTPUT TO UNIT 16

            K = 1
            II = 1
            DO I=1,global_here%NFFR
               WRITE(16,2206) I,global_here%FBOUNTAG(I)
 2206          FORMAT(////,5X,'PERIODIC NORMAL FLOW CONSTITUENT ',&
              'NUMBER',I4,1X,'DESIGNATED : ',A5)
               READ(15,'(A10)') global_here%ALPHA
               WRITE(16,31) global_here%ALPHA
               WRITE(16,30)
               DO J=1,global_here%NVEL
                  IF((global_here%LBCODEI(J).EQ.2).OR.(global_here%LBCODEI(J).EQ.12)&
                 .OR.(global_here%LBCODEI(J).EQ.22)) THEN
                     
!.....Modified arrangement of global_here%QNAM and global_here%QNPH for DG

                     IF (global_here%DGSWE.EQ.1) THEN
                        READ(15,*) global_here%QNAM(I,II),global_here%QNPH(I,II)
                        WRITE(16,2205) global_here%NBV(J),global_here%QNAM(I,II),global_here%QNPH(I,II)
                        global_here%QNPH(I,II)=global_here%QNPH(I,II)*DEG2RAD
                        II = II + 1
                     ELSE
                        READ(15,*) global_here%QNAM(I,J),global_here%QNPH(I,J)
                        WRITE(16,2205) global_here%NBV(J),global_here%QNAM(I,J),global_here%QNPH(I,J)
 2205                   FORMAT(10X,I8,4X,F14.5,4X,F12.3)
                        global_here%QNPH(I,J)=global_here%QNPH(I,J)*DEG2RAD
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
            II = 1
         ENDIF
      ENDIF

!...  IF ANY EXTERNAL BARRIER BOUNDARIES WERE SPECIFIED, (global_here%NFLUXB=1)
!.....WRITE OUT EXTERNAL BARRIER BOUNDARY INFORMATION TO UNIT 16 FILE
!.....NOTE THAT THIS INFORMATION WAS READ IN FROM THE UNIT 14 FILE

      IF(global_here%NFLUXB.EQ.1) THEN

!.....WRITE OUT INFO ON SPECIFIED EXTERNAL BARRIER BOUNDARIES

         WRITE(16,1112)
         WRITE(16,2220)
 2220    FORMAT(//,1X,'EXTERNAL BARRIER BOUNDARY INFORMATION ',/)

!.......OUTPUT ELEVATION OF EXTERNAL BARRIER NODES ABOVE THE GEOID AND
!........THE COEFFICIENT OF FREE SURFACE SUPERCRITICAL FLOW AT
!........DESIGNATED EXTERNAL BARRIER BOUNDARY NODES TO UNIT 16
         WRITE(16,2224)
 2224    FORMAT(//,9X,'NODE',10X,'BARRIER HEIGHT',&
        6X,'SUPER-CRIT. EXTERNAL BAR. global_here%COEF.',/)
         DO J=1,global_here%NVEL
            IF((global_here%LBCODEI(J).EQ.3).OR.(global_here%LBCODEI(J).EQ.13)&
           .OR.(global_here%LBCODEI(J).EQ.23)) THEN
               WRITE(16,2225) global_here%NBV(J),global_here%BARLANHT(J),global_here%BARLANCFSP(J)
 2225          FORMAT(5X,I8,6X,F14.5,15X,F12.3)
            ENDIF
         END DO
      ENDIF

!...  IF ANY INTERNAL BARRIER BOUNDARIES WERE SPECIFIED, (global_here%NFLUXIB=1)
!.....WRITE INTERNAL BARRIER BOUNDARY INFORMATION TO UNIT 16 FILE

      IF(global_here%NFLUXIB.EQ.1) THEN

!.....WRITE OUT INFO ON SPECIFIED INTERNAL BARRIER BOUNDARIES

         WRITE(16,1112)
         WRITE(16,2320)
 2320    FORMAT(//,1X,'INTERNAL BARRIER BOUNDARY INFORMATION ',/)

!.......WRITE CONNECTION NODE NUMBER AND ELEVATION OF THE INTERNAL BARRIER
!........NODES ABOVE THE GEOID AND THE COEFFICIENTS OF FREE SURFACE SUPERCRITICAL
!........AND SUBCRITICAL FLOW AT DESIGNATED INTERNAL BARRIER BOUNDARY NODES
!........TO UNIT 16 (NOTE THAT THIS INFORMATION WAS INPUT FROM THE UNIT 14
!........FILE WITH BOUNDARY NODE INFORMATION)
         WRITE(16,2324)
 2324    FORMAT(//,9X,'NODE',6X,'CONNECTED NODE',6X,'BARRIER HEIGHT',&
        4X,'SUB-CRIT. INT. BAR. global_here%COEF.',&
        4X,'SUPER-CRIT. INT. BAR. global_here%COEF.',/)

         DO J=1,global_here%NVEL
            IF((global_here%LBCODEI(J).EQ.4).OR.(global_here%LBCODEI(J).EQ.24)) THEN

               WRITE(16,2325) global_here%NBV(J),global_here%IBCONN(J),global_here%BARINHT(J),&
              global_here%BARINCFSB(J),global_here%BARINCFSP(J)
 2325          FORMAT(5X,I8,7X,I8,6X,F14.5,12X,F12.3,17X,F12.3)
            ENDIF
         END DO
      ENDIF

!jj   wm001 - begin add                        
!...  IF ANY INTERNAL BARRIER BOUNDARIES WITH global_here%CROSS BARRIER PIPES
!.....WERE SPECIFIED, (global_here%NFLUXIBP=1)
!.....WRITE INTERNAL BARRIER BOUNDARY INFORMATION WITH global_here%CROSS 
!.....BARRIER PIPE INFORMATION TO UNIT 16 FILE

      IF(global_here%NFLUXIBP.EQ.1) THEN

!.....WRITE OUT INFO ON SPECIFIED INTERNAL BARRIER BOUNDARIES

         WRITE(16,1112)
         WRITE(16,2326)
 2326    FORMAT(//,1X,'INTERNAL BARRIER BOUNDARY WITH global_here%CROSS BARRIER',&
        ' PIPE INFORMATION ',/)

!.......WRITE CONNECTION NODE NUMBER AND ELEVATION OF THE INTERNAL BARRIER
!........NODES ABOVE THE GEOID AND THE COEFFICIENTS OF FREE SURFACE SUPERCRITICAL
!........AND SUBCRITICAL FLOW AT DESIGNATED INTERNAL BARRIER BOUNDARY NODES
!........IN ADDITION TO global_here%CROSS BARRIER PIPE CROWN HEIGHT, global_here%CROSS BARRIER PIPE
!........COEFFICIENT AND global_here%CROSS BARRIER PIPE DIAMETER TO UNIT 16
!........(NOTE THAT THIS INFORMATION WAS INPUT FROM THE UNIT 14 FILE WITH 
!........BOUNDARY NODE INFORMATION)
         WRITE(16,2327)
 2327    FORMAT(//,7X,'NODE',4X,'CONNECTED NODE',4X,'BARRIER HEIGHT',&
        4X,'SUB-CRIT INT BAR global_here%COEF',&
        4X,'SUPER-CRIT INT BAR global_here%COEF',&
        4X,'global_here%PIPEHT  ',&
        4X,'global_here%PIPECOEF',&
        4X,'global_here%PIPEDIAM',/)
         DO J=1,global_here%NVEL
            IF((global_here%LBCODEI(J).EQ.5).OR.(global_here%LBCODEI(J).EQ.25)) THEN
               WRITE(16,2328) global_here%NBV(J),global_here%IBCONN(J),global_here%BARINHT(J),&
              global_here%BARINCFSB(J),global_here%BARINCFSP(J),&
              global_here%PIPEHT(J),global_here%PIPECOEF(J),global_here%PIPEDIAM(J)
 2328          FORMAT(3X,I8,5X,I8,4X,F14.5,8X,F12.3,12X,F12.3,&
              2X,F10.5,2X,F10.5,2X,F10.5)
            ENDIF
         END DO
      ENDIF
!jj   wm001 - end add                        

!...  
!...  READ IN INFORMATION CONCERNING OUTPUT REQUIREMENTS FROM UNIT 15 AND
!...  OUTPUT THIS TO UNIT 16
!...  
      WRITE(16,1112)
      WRITE(16,3000)
 3000 FORMAT(//,1X,'OUTPUT INFORMATION WILL BE PROVIDED AS'&
     ,' FOLLOWS :')

!...  
!...  INPUT INFORMATION FOR ELEVATION RECORDING STATIONS
!...  

!.... READ IN global_here%NOUTE,global_here%TOUTSE,global_here%TOUTFE,global_here%NSPOOLE : IF ABS(global_here%NOUTE)>0, INTERPOLATED
!.... ELEVATIONS AT ELEVATION STATIONS ARE SPOOLED TO UNIT 61 EVERY global_here%NSPOOLE
!.... TIME STEPS BETWEEN TIMES global_here%TOUTSE AND global_here%TOUTFE

      READ(15,*) global_here%NOUTE,global_here%TOUTSE,global_here%TOUTFE,global_here%NSPOOLE
      WRITE(16,3001) global_here%NOUTE
 3001 FORMAT(///,1X,'ELEVATION RECORDING STATION OUTPUT : ',&
     //,5X,'global_here%NOUTE = ',I2)

!.... CHECK INPUT PARAMETER global_here%NOUTE

      IF(ABS(global_here%NOUTE).GT.2) THEN
         IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,3002)
         WRITE(16,3002)
 3002    FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',&
        //,1X,'YOUR SELECTION OF THE UNIT 15 INPUT PARAMETER',&
        ' global_here%NOUTE',&
        /,1X,'IS NOT AN ALLOWABLE VALUE.  CHECK YOUR INPUT!!')
         IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9973)
         WRITE(16,9973)
         STOP
      ENDIF

!.... IF STATION ELEVATION OUTPUT WILL NOT BE GENERATED

      IF(global_here%NOUTE.EQ.0) THEN
         WRITE(16,3003)
 3003    FORMAT(/,5X,'NO OUTPUT WILL BE SPOOLED AT ELEVATION ',&
        'RECORDING STATIONS')
      ENDIF

!.... IF STATION ELEVATION OUTPUT WILL BE GENERATED

      IF(global_here%NOUTE.NE.0) THEN

!......COMPUTE global_here%NTCYSE, global_here%NTCYFE, WHICH = global_here%TOUTSE AND global_here%TOUTFE IN TIMESTEPS

         global_here%NTCYSE=INT((global_here%TOUTSE-global_here%STATIM)*(86400.D0/global_here%DTDP)+0.5d0)
         global_here%NTCYFE=INT((global_here%TOUTFE-global_here%STATIM)*(86400.D0/global_here%DTDP)+0.5d0)

         IF(global_here%NTCYFE.GT.global_here%NT) global_here%NTCYFE=global_here%NT

!......COMPUTE global_here%NTRSPE = THE NO. OF DATA SETS TO BE SPOOLED TO UNIT 61

         IF(global_here%NSPOOLE.EQ.0) global_here%NTRSPE=0
         IF(global_here%NSPOOLE.NE.0) global_here%NTRSPE=INT((global_here%NTCYFE-global_here%NTCYSE)/global_here%NSPOOLE)

!......WRITE global_here%TOUTSE,global_here%TOUTFE,global_here%NTCYSE,global_here%NTCYFE,global_here%NSPOOLE TO UNIT 16

         WRITE(16,3004) global_here%TOUTSE,global_here%NTCYSE,global_here%TOUTFE,global_here%NTCYFE,global_here%NSPOOLE
 3004    FORMAT(/,5X,'DATA RECORDS WILL START AFTER global_here%TOUTSE =',F8.3,&
        ' global_here%DAY(S) RELATIVE',&
        /,9X,'TO THE STARTING TIME OR',I9,&
        ' TIME STEPS INTO THE SIMULATION',&
        //,5X,'DATA RECORDS WILL STOP AFTER global_here%TOUTFE =',F8.3,&
        ' global_here%DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',&
        I9,' TIME STEPS INTO THE SIMULATION',&
        //,5X,'INFORMATION WILL BE SPOOLED TO UNIT 61 EVERY',&
        ' global_here%NSPOOLE =',I8,' TIME STEPS')
         IF(ABS(global_here%NOUTE).EQ.1) WRITE(16,3005)
 3005    FORMAT(/,5X,'UNIT 61 FORMAT WILL BE ASCII')
         IF(ABS(global_here%NOUTE).EQ.2) WRITE(16,3006)
 3006    FORMAT(/,5X,'UNIT 61 FORMAT WILL BE BINARY')
      ENDIF

!.... REGARDLESS OF WHETHER global_here%NOUTE=0, READ IN THE NUMBER OF ELEVATION
!.... RECORDING STATIONS

      READ(15,*) global_here%NSTAE
      WRITE(16,3007) global_here%NSTAE
 3007 FORMAT(///,5X,'NUMBER OF INPUT ELEVATION RECORDING STATIONS = ',&
     I5)

      IF(global_here%NSTAE.GT.0) THEN
         IF(global_here%ICS.EQ.1) WRITE(16,3008)
 3008    FORMAT(/,7X,'STATION #   ELEMENT',9X,'global_here%X',13X,'global_here%Y',/)
         IF(global_here%ICS.EQ.2) WRITE(16,3009)
 3009    FORMAT(/,5X,'STATION   ELEMENT',3X,'LAMBDA(DEG)',&
        4X,'FEA(DEG)',10X,'XCP',12X,'YCP',/)
         s%MNSTAE = global_here%NSTAE
      ENDIF
      IF (global_here%NSTAE.EQ.0) s%MNSTAE = 1


!     Allocate arrays dimensioned by MNSTAE
      call alloc_main7(s,global_here)


!.... INPUT COORDINATES OF ELEVATION RECORDING STATIONS THEN COMPUTE
!.... THE ELEMENT NO. THE STATION LIES IN

      DO I=1,global_here%NSTAE
         global_here%NNE(I)=0
         IF(global_here%ICS.EQ.1) THEN
            READ(15,*) global_here%XEL(I),global_here%YEL(I)
         ELSE
            READ(15,*) global_here%SLEL(I),global_here%SFEL(I)
            global_here%SLEL(I)=global_here%SLEL(I)*DEG2RAD
            global_here%SFEL(I)=global_here%SFEL(I)*DEG2RAD
            CALL CPP(global_here%XEL(I),global_here%YEL(I),global_here%SLEL(I),global_here%SFEL(I),global_here%SLAM0,global_here%SFEA0)
         ENDIF
         global_here%AEMIN=1.0E+25
         global_here%KMIN=0
         DO K=1,global_here%NE
            global_here%N1=global_here%NM(K,1)
            global_here%N2=global_here%NM(K,2)
            global_here%N3=global_here%NM(K,3)
            global_here%X1=global_here%X(global_here%N1)
            global_here%X2=global_here%X(global_here%N2)
            global_here%X3=global_here%X(global_here%N3)
            global_here%X4=global_here%XEL(I)
            global_here%Y1=global_here%Y(global_here%N1)
            global_here%Y2=global_here%Y(global_here%N2)
            global_here%Y3=global_here%Y(global_here%N3)
            global_here%Y4=global_here%YEL(I)
            global_here%A1=(global_here%X4-global_here%X3)*(global_here%Y2-global_here%Y3)+(global_here%X2-global_here%X3)*(global_here%Y3-global_here%Y4)
            global_here%A2=(global_here%X4-global_here%X1)*(global_here%Y3-global_here%Y1)-(global_here%Y4-global_here%Y1)*(global_here%X3-global_here%X1)
            global_here%A3=(global_here%Y4-global_here%Y1)*(global_here%X2-global_here%X1)-(global_here%X4-global_here%X1)*(global_here%Y2-global_here%Y1)
            global_here%AA=ABS(global_here%A1)+ABS(global_here%A2)+ABS(global_here%A3)
            global_here%AE=ABS(global_here%AA-global_here%AREAS(K))/global_here%AREAS(K)
            IF(global_here%AE.LT.global_here%AEMIN) THEN
               global_here%AEMIN=global_here%AE
               global_here%KMIN=K
            ENDIF
            IF(global_here%AE.LT.1.0E-5) global_here%NNE(I)=K
         END DO

         IF(global_here%NNE(I).EQ.0) THEN
            IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9784) I
            WRITE(16,9784) I
 9784       FORMAT(///,1X,'!!!!!!!!!!  WARNING - NONFATAL ',&
           'INPUT ERROR  !!!!!!!!!',//&
           ,1X,'ELEVATION RECORDING STATION ',I6,' DOES NOT LIE',&
           ' WITHIN ANY ELEMENT IN THE DEFINED',&
           /,1X,'COMPUTATIONAL DOMAIN,   PLEASE CHECK THE INPUT',&
           ' COORDINATES FOR THIS STATION')
            IF(global_here%NFOVER.EQ.1) THEN
               IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9790) global_here%AEMIN
               WRITE(16,9790) global_here%AEMIN
 9790          FORMAT(/,1X,'PROGRAM WILL ESTIMATE NEAREST ELEMENT',&
              /,1X,'PROXIMITY INDEX FOR THIS STATION EQUALS ',E15.6,&
              //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)
               global_here%NNE(I)=global_here%KMIN
            ELSE
               IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9791) global_here%AEMIN
               WRITE(16,9791) global_here%AEMIN
 9791          FORMAT(/,1X,'PROGRAM WILL NOT CORRECT ERROR ',&
              'SINCE NON-FATAL ERROR OVERIDE OPTION, global_here%NFOVER,',&
          /,1X,'HAS BEEN SELECTED EQUAL TO 0',&
          /,1X,'PROXIMITY INDEX FOR THIS STATION EQUALS ',E15.6,&
          //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',&
              //)
               STOP
            ENDIF
         ENDIF

         IF(global_here%ICS.EQ.1) THEN
            WRITE(16,1880) I,global_here%NNE(I),global_here%XEL(I),global_here%YEL(I)
 1880       FORMAT(8X,I3,6X,I7,2(2X,F14.2))
         ELSE
            WRITE(16,1883) I,global_here%NNE(I),global_here%SLEL(I)*RAD2DEG,&
           global_here%SFEL(I)*RAD2DEG,global_here%XEL(I),global_here%YEL(I)
 1883       FORMAT(6X,I3,4X,I7,2(2X,F13.8),2X,2(1X,F13.2))
         ENDIF

!.... PRE-COMPUTE INFORMATION REQUIRED TO INTERPOLATE AT ELEV. RECORDING STATIONS

         global_here%N1=global_here%NM(global_here%NNE(I),1)
         global_here%N2=global_here%NM(global_here%NNE(I),2)
         global_here%N3=global_here%NM(global_here%NNE(I),3)
         global_here%X1=global_here%X(global_here%N1)
         global_here%X2=global_here%X(global_here%N2)
         global_here%X3=global_here%X(global_here%N3)
         global_here%X4=global_here%XEL(I)
         global_here%Y1=global_here%Y(global_here%N1)
         global_here%Y2=global_here%Y(global_here%N2)
         global_here%Y3=global_here%Y(global_here%N3)
         global_here%Y4=global_here%YEL(I)
         global_here%STAIE1(I)=((global_here%X4-global_here%X3)*(global_here%Y2-global_here%Y3)+(global_here%X2-global_here%X3)*(global_here%Y3-global_here%Y4))/global_here%AREAS(global_here%NNE(I))
         global_here%STAIE2(I)=((global_here%X4-global_here%X1)*(global_here%Y3-global_here%Y1)-(global_here%Y4-global_here%Y1)*(global_here%X3-global_here%X1))/global_here%AREAS(global_here%NNE(I))
         global_here%STAIE3(I)=(-(global_here%X4-global_here%X1)*(global_here%Y2-global_here%Y1)+(global_here%Y4-global_here%Y1)*(global_here%X2-global_here%X1))/global_here%AREAS(global_here%NNE(I))

      END DO

!...  
!...  INPUT INFORMATION FOR VELOCITY RECORDING STATIONS
!...  

!.... READ IN global_here%NOUTV,global_here%TOUTSV,global_here%TOUTFV,global_here%NSPOOLV : IF global_here%NOUTV<>0,INTERPOLATED VELOCITIES AT
!.... VELOCITY STATIONS ARE SPOOLED TO UNIT 62 EVERY global_here%NSPOOLV TIME STEPS BETWEEN
!.... TIMES global_here%TOUTSV AND global_here%TOUTFV; IF ABS(global_here%NOUTV)=2, OUTPUT WILL BE BINARY

      READ(15,*) global_here%NOUTV,global_here%TOUTSV,global_here%TOUTFV,global_here%NSPOOLV
      WRITE(16,3101) global_here%NOUTV
 3101 FORMAT(////,1X,'VELOCITY RECORDING STATION OUTPUT : ',&
     //,5X,'global_here%NOUTV = ',I2)

!.... CHECK INPUT PARAMETER global_here%NOUTV

      IF(ABS(global_here%NOUTV).GT.2) THEN
         IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,3102)
         WRITE(16,3102)
 3102    FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',&
        //,1X,'YOUR SELECTION OF THE UNIT 15 INPUT PARAMETER',&
        ' global_here%NOUTV',&
        /,1X,'IS NOT AN ALLOWABLE VALUE.  CHECK YOUR INPUT!!')
         IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9973)
         WRITE(16,9973)
         STOP
      ENDIF

!.... IF STATION VELOCITY OUTPUT WILL NOT BE GENERATED

      IF(global_here%NOUTV.EQ.0) THEN
         WRITE(16,3103)
 3103    FORMAT(///,5X,'NO OUTPUT WILL BE SPOOLED AT VELOCITY',&
        ' RECORDING STATIONS')
      ENDIF

!.... IF STATION VELOCITY OUTPUT WILL BE GENERATED

      IF(global_here%NOUTV.NE.0) THEN

!......COMPUTE global_here%NTCYSV, global_here%NTCYFV, WHICH = global_here%TOUTSV AND global_here%TOUTFV IN TIME STEPS

         global_here%NTCYSV=INT((global_here%TOUTSV-global_here%STATIM)*(86400.D0/global_here%DTDP) + 0.5d0)
         global_here%NTCYFV=INT((global_here%TOUTFV-global_here%STATIM)*(86400.D0/global_here%DTDP) + 0.5d0)
         IF(global_here%NTCYFV.GT.global_here%NT) global_here%NTCYFV=global_here%NT

!......CALCULATE global_here%NTRSPV = THE NO. OF DATA SETS TO BE SPOOLED TO UNIT 62

         IF(global_here%NSPOOLV.EQ.0) global_here%NTRSPV=0
         IF(global_here%NSPOOLV.NE.0) global_here%NTRSPV=INT((global_here%NTCYFV-global_here%NTCYSV)/global_here%NSPOOLV)

!......WRITE global_here%NOUTV,global_here%TOUTSV,global_here%TOUTFV,global_here%NTCYSV,global_here%NTCYFV,global_here%NSPOOLV TO UNIT 16

         WRITE(16,3104) global_here%TOUTSV,global_here%NTCYSV,global_here%TOUTFV,global_here%NTCYFV,global_here%NSPOOLV
 3104    FORMAT(/,5X,'DATA RECORDS WILL START AFTER global_here%TOUTSV =',F8.3,&
        ' global_here%DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',&
        I9,' TIME STEPS INTO THE SIMULATION',&
        //,5X,'DATA RECORDS WILL STOP AFTER global_here%TOUTFV =',F8.3,&
        ' global_here%DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',&
        I9,' TIME STEPS INTO THE SIMULATION',&
        //,5X,'INFORMATION WILL BE SPOOLED TO UNIT 62 EVERY ',&
        ' global_here%NSPOOLV =',I8,' TIME STEPS')
         IF(ABS(global_here%NOUTV).EQ.1) WRITE(16,3105)
 3105    FORMAT(/,5X,'UNIT 62 FORMAT WILL BE ASCII')
         IF(ABS(global_here%NOUTV).EQ.2) WRITE(16,3106)
 3106    FORMAT(/,5X,'UNIT 62 FORMAT WILL BE BINARY')
      ENDIF

!.... REGARDLESS OF WHETHER global_here%NOUTV=0, READ IN THE NUMBER OF VELOCITY
!.... RECORDING STATIONS

      READ(15,*) global_here%NSTAV
      WRITE(16,3107) global_here%NSTAV
 3107 FORMAT(////,5X,'NUMBER OF INPUT VELOCITY RECORDING STATIONS = ',&
     I5)

      IF(global_here%NSTAV.GT.0) THEN
         IF(global_here%ICS.EQ.1) WRITE(16,3108)
 3108    FORMAT(/,7X,'STATION #   ELEMENT',9X,'global_here%X',13X,'global_here%Y',/)
         IF(global_here%ICS.EQ.2) WRITE(16,3109)
 3109    FORMAT(/,5X,'STATION   ELEMENT',3X,'LAMBDA(DEG)',&
        4X,'FEA(DEG)',10X,'XCP',12X,'YCP',/)
         s%MNSTAV = global_here%NSTAV
      ENDIF
      IF (global_here%NSTAV.EQ.0) s%MNSTAV = 1

!     Allocate arrays dimensioned by MNSTAV
      call alloc_main8(s,global_here)

!.... INPUT COORDINATES OF VELOCITY RECORDING STATIONS
!.... THEN COMPUTE ELEMENT NO. WITHIN WHICH STATION LIES

      DO I=1,global_here%NSTAV
         global_here%NNV(I)=0
         IF(global_here%ICS.EQ.1) THEN
            READ(15,*) global_here%XEV(I),global_here%YEV(I)
         ELSE
            READ(15,*) global_here%SLEV(I),global_here%SFEV(I)
            global_here%SLEV(I)=global_here%SLEV(I)*DEG2RAD
            global_here%SFEV(I)=global_here%SFEV(I)*DEG2RAD
            CALL CPP(global_here%XEV(I),global_here%YEV(I),global_here%SLEV(I),global_here%SFEV(I),global_here%SLAM0,global_here%SFEA0)
         ENDIF
         global_here%AEMIN=1.0E+25
         global_here%KMIN=0
         DO K=1,global_here%NE
            global_here%N1=global_here%NM(K,1)
            global_here%N2=global_here%NM(K,2)
            global_here%N3=global_here%NM(K,3)
            global_here%X1=global_here%X(global_here%N1)
            global_here%X2=global_here%X(global_here%N2)
            global_here%X3=global_here%X(global_here%N3)
            global_here%X4=global_here%XEV(I)
            global_here%Y1=global_here%Y(global_here%N1)
            global_here%Y2=global_here%Y(global_here%N2)
            global_here%Y3=global_here%Y(global_here%N3)
            global_here%Y4=global_here%YEV(I)
            global_here%A1=(global_here%X4-global_here%X3)*(global_here%Y2-global_here%Y3)+(global_here%X2-global_here%X3)*(global_here%Y3-global_here%Y4)
            global_here%A2=(global_here%X4-global_here%X1)*(global_here%Y3-global_here%Y1)-(global_here%Y4-global_here%Y1)*(global_here%X3-global_here%X1)
            global_here%A3=(global_here%Y4-global_here%Y1)*(global_here%X2-global_here%X1)-(global_here%X4-global_here%X1)*(global_here%Y2-global_here%Y1)
            global_here%AA=ABS(global_here%A1)+ABS(global_here%A2)+ABS(global_here%A3)
            global_here%AE=ABS(global_here%AA-global_here%AREAS(K))/global_here%AREAS(K)
            IF(global_here%AE.LT.global_here%AEMIN) THEN
               global_here%AEMIN=global_here%AE
               global_here%KMIN=K
            ENDIF
            IF(global_here%AE.LT.1.0E-5) global_here%NNV(I)=K
         END DO

         IF(global_here%NNV(I).EQ.0) THEN
            IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9786) I
            WRITE(16,9786) I
 9786       FORMAT(///,1X,'!!!!!!!!!!  WARNING - NONFATAL ',&
           'INPUT ERROR  !!!!!!!!!',//&
           ,1X,'VELOCITY RECORDING STATION ',I6,' DOES NOT LIE'&
           ,' WITHIN ANY ELEMENT IN THE DEFINED',&
           /,1X,'COMPUTATIONAL DOMAIN,   PLEASE CHECK THE INPUT'&
           ,' COORDINATES FOR THIS STATION')
            IF(global_here%NFOVER.EQ.1) THEN
               IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9790) global_here%AEMIN
               WRITE(16,9790) global_here%AEMIN
               global_here%NNV(I)=global_here%KMIN
            ELSE
               IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9791) global_here%AEMIN
               WRITE(16,9791) global_here%AEMIN
               STOP
            ENDIF
         ENDIF

         IF(global_here%ICS.EQ.1) THEN
            WRITE(16,1880) I,global_here%NNV(I),global_here%XEV(I),global_here%YEV(I)
         ELSE
            WRITE(16,1883) I,global_here%NNV(I),global_here%SLEV(I)*RAD2DEG,global_here%SFEV(I)*RAD2DEG,&
           global_here%XEV(I),global_here%YEV(I)
         ENDIF

!.... PRE-COMPUTE INFORMATION REQUIRED TO INTERPOLATE AT global_here%VEL. RECORDING STATIONS

         global_here%N1=global_here%NM(global_here%NNV(I),1)
         global_here%N2=global_here%NM(global_here%NNV(I),2)
         global_here%N3=global_here%NM(global_here%NNV(I),3)
         global_here%X1=global_here%X(global_here%N1)
         global_here%X2=global_here%X(global_here%N2)
         global_here%X3=global_here%X(global_here%N3)
         global_here%X4=global_here%XEV(I)
         global_here%Y1=global_here%Y(global_here%N1)
         global_here%Y2=global_here%Y(global_here%N2)
         global_here%Y3=global_here%Y(global_here%N3)
         global_here%Y4=global_here%YEV(I)
         global_here%STAIV1(I)=((global_here%X4-global_here%X3)*(global_here%Y2-global_here%Y3)+(global_here%X2-global_here%X3)*(global_here%Y3-global_here%Y4))/global_here%AREAS(global_here%NNV(I))
         global_here%STAIV2(I)=((global_here%X4-global_here%X1)*(global_here%Y3-global_here%Y1)-(global_here%Y4-global_here%Y1)*(global_here%X3-global_here%X1))/global_here%AREAS(global_here%NNV(I))
         global_here%STAIV3(I)=(-(global_here%X4-global_here%X1)*(global_here%Y2-global_here%Y1)+(global_here%Y4-global_here%Y1)*(global_here%X2-global_here%X1))/global_here%AREAS(global_here%NNV(I))

      END DO

!...  
!...  IF TRANSPORT IS INCLUDED IN THE RUN, INPUT INFORMATION FOR CONCENTRATION
!...  RECORDING STATIONS
!...  
      global_here%NOUTC=0
      IF(global_here%IM.EQ.10) THEN

!.....READ IN global_here%NOUTC,global_here%TOUTSC,global_here%TOUTFC,global_here%NSPOOLC : IF global_here%NOUTC<>0,INTERPOLATED
!.....CONCENTRATIONS ARE SPOOLED TO UNIT 81 EVERY global_here%NSPOOLC TIME STEPS
!.....BETWEEN TIMES global_here%TOUTSC AND global_here%TOUTFC; IF ABS(global_here%NOUTC)=2, OUTPUT WILL BE BINARY

         READ(15,*) global_here%NOUTC,global_here%TOUTSC,global_here%TOUTFC,global_here%NSPOOLC
         WRITE(16,3201) global_here%NOUTC
 3201    FORMAT(///,1X,'CONCENTRATION RECORDING STATION OUTPUT : ',&
        //,5X,'global_here%NOUTC = ',I2)

!.....CHECK INPUT PARAMETER global_here%NOUTC

         IF(ABS(global_here%NOUTC).GT.2) THEN
            IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,3202)
            WRITE(16,3202)
 3202       FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',&
           //,1X,'YOUR SELECTION OF THE UNIT 15 INPUT PARAMETER',&
           ' global_here%NOUTC',&
           /,1X,'IS NOT AN ALLOWABLE VALUE.  CHECK YOUR INPUT!!')
            IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9973)
            WRITE(16,9973)
            STOP
         ENDIF

!.....IF STATION CONCENTRATION OUTPUT WILL NOT BE GENERATED

         IF(global_here%NOUTC.EQ.0) THEN
            WRITE(16,3203)
 3203       FORMAT(/,5X,'NO OUTPUT WILL BE SPOOLED AT CONCENTRATION',&
           ' RECORDING STATIONS')
         ENDIF

!.....IF STATION CONCENTRATION OUTPUT WILL BE GENERATED

         global_here%NSTAC = 0
         IF(global_here%NOUTC.NE.0) THEN

!.......COMPUTE global_here%NTCYSC, global_here%NTCYFC, WHICH = global_here%TOUTSC AND global_here%TOUTFC IN TIMESTEPS

            global_here%NTCYSC=INT((global_here%TOUTSC-global_here%STATIM)*(86400.D0/global_here%DTDP) + 0.5d0)
            global_here%NTCYFC=INT((global_here%TOUTFC-global_here%STATIM)*(86400.D0/global_here%DTDP) + 0.5d0)
            IF(global_here%NTCYFC.GT.global_here%NT) global_here%NTCYFC=global_here%NT

!.......COMPUTE global_here%NTRSPC = THE NO. OF DATA SETS TO BE SPOOLED TO UNIT 81

            IF(global_here%NSPOOLC.EQ.0) global_here%NTRSPC=0
            IF(global_here%NSPOOLC.NE.0) global_here%NTRSPC=INT((global_here%NTCYFC-global_here%NTCYSC)/global_here%NSPOOLC)

!.......WRITE global_here%TOUTSC,global_here%TOUTFC,global_here%NTCYSC,global_here%NTCYFC,global_here%NSPOOLC TO UNIT 16

            WRITE(16,3204) global_here%TOUTSC,global_here%NTCYSC,global_here%TOUTFC,global_here%NTCYFC,global_here%NSPOOLC
 3204       FORMAT(/,5X,'DATA RECORDS WILL START AFTER global_here%TOUTSC =',F8.3,&
           ' global_here%DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',&
           I9,' TIME STEPS INTO THE SIMULATION',&
           //,5X,'DATA RECORDS WILL STOP AFTER global_here%TOUTFC =',F8.3,&
           ' global_here%DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',&
           I9,' TIME STEPS INTO THE SIMULATION',&
           //,5X,'INFORMATION WILL BE SPOOLED TO UNIT 81 EVERY',&
           ' global_here%NSPOOLC =',I8,' TIME STEPS')
            IF(ABS(global_here%NOUTC).EQ.1) WRITE(16,3205)
 3205       FORMAT(/,5X,'UNIT 81 FORMAT WILL BE ASCII')
            IF(ABS(global_here%NOUTC).EQ.2) WRITE(16,3206)
 3206       FORMAT(/,5X,'UNIT 81 FORMAT WILL BE BINARY')
         ENDIF

!.....REGARDLESS OF WHETHER global_here%NOUTC=0, READ IN THE NUMBER OF CONCENTRATION
!.....RECORDING STATIONS

         READ(15,*) global_here%NSTAC
         WRITE(16,3207) global_here%NSTAC
 3207    FORMAT(///,5X,'NUMBER OF INPUT CONCENTRATION RECORDING ',&
        'STATIONS = ',I5)

         IF(global_here%NSTAC.GT.0) THEN
            IF(global_here%ICS.EQ.1) WRITE(16,3208)
 3208       FORMAT(/,7X,'STATION #   ELEMENT',9X,'global_here%X',13X,'global_here%Y',/)
            IF(global_here%ICS.EQ.2) WRITE(16,3209)
 3209       FORMAT(/,5X,'STATION   ELEMENT',3X,'LAMBDA(DEG)',&
           4X,'FEA(DEG)',10X,'XCP',12X,'YCP',/)
            s%MNSTAC = global_here%NSTAC
         ENDIF

!     Allocate arrays dimensioned by MNSTAC
         call alloc_main9(s,global_here)

!.....INPUT COORDINATES OF CONCENTRATION RECORDING STATIONS
!.....THEN COMPUTE ELEMENT NO. WITHIN WHICH STATION LIES

         DO I=1,global_here%NSTAC
            global_here%NNC(I)=0
            IF(global_here%ICS.EQ.1) THEN
               READ(15,*) global_here%XEC(I),global_here%YEC(I)
            ELSE
               READ(15,*) global_here%SLEC(I),global_here%SFEC(I)
               global_here%SLEC(I)=global_here%SLEC(I)*DEG2RAD
               global_here%SFEC(I)=global_here%SFEC(I)*DEG2RAD
               CALL CPP(global_here%XEC(I),global_here%YEC(I),global_here%SLEC(I),global_here%SFEC(I),global_here%SLAM0,global_here%SFEA0)
            ENDIF
            global_here%AEMIN=1.0E+25
            global_here%KMIN=0
            DO K=1,global_here%NE
               global_here%N1=global_here%NM(K,1)
               global_here%N2=global_here%NM(K,2)
               global_here%N3=global_here%NM(K,3)
               global_here%X1=global_here%X(global_here%N1)
               global_here%X2=global_here%X(global_here%N2)
               global_here%X3=global_here%X(global_here%N3)
               global_here%X4=global_here%XEC(I)
               global_here%Y1=global_here%Y(global_here%N1)
               global_here%Y2=global_here%Y(global_here%N2)
               global_here%Y3=global_here%Y(global_here%N3)
               global_here%Y4=global_here%YEC(I)
               global_here%A1=(global_here%X4-global_here%X3)*(global_here%Y2-global_here%Y3)+(global_here%X2-global_here%X3)*(global_here%Y3-global_here%Y4)
               global_here%A2=(global_here%X4-global_here%X1)*(global_here%Y3-global_here%Y1)-(global_here%Y4-global_here%Y1)*(global_here%X3-global_here%X1)
               global_here%A3=(global_here%Y4-global_here%Y1)*(global_here%X2-global_here%X1)-(global_here%X4-global_here%X1)*(global_here%Y2-global_here%Y1)
               global_here%AA=ABS(global_here%A1)+ABS(global_here%A2)+ABS(global_here%A3)
               global_here%AE=ABS(global_here%AA-global_here%AREAS(K))/global_here%AREAS(K)
               IF(global_here%AE.LT.global_here%AEMIN) THEN
                  global_here%AEMIN=global_here%AE
                  global_here%KMIN=K
               ENDIF
               IF(global_here%AE.LT.1.0E-5) global_here%NNC(I)=K
            END DO

            IF(global_here%NNC(I).EQ.0) THEN
               IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9785) I
               WRITE(16,9785) I
 9785          FORMAT(///,1X,'!!!!!!!!!!  WARNING - NONFATAL INPUT ERROR ',&
          '!!!!!!!!!',//,&
          ' CONCENTRATION RECORDING STATION ',I6,' DOES NOT LIE'&
          ,' WITHIN ANY ELEMENT IN THE DEFINED',/,&
          ' COMPUTATIONAL DOMAIN,   PLEASE CHECK THE INPUT',&
          ' COORDINATES FOR THIS STATION')
               IF(global_here%NFOVER.EQ.1) THEN
                  IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9790) global_here%AEMIN
                  WRITE(16,9790) global_here%AEMIN
                  global_here%NNC(I)=global_here%KMIN
               ELSE
                  IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9791) global_here%AEMIN
                  WRITE(16,9791) global_here%AEMIN
                  STOP
               ENDIF
            ENDIF

            IF(global_here%ICS.EQ.1) THEN
               WRITE(16,1880) I,global_here%NNC(I),global_here%XEC(I),global_here%YEC(I)
            ELSE
               WRITE(16,1883) I,global_here%NNC(I),global_here%SLEC(I)*RAD2DEG,&
              global_here%SFEC(I)*RAD2DEG,global_here%XEC(I),global_here%YEC(I)
            ENDIF

!.....PRE-COMPUTE INFORMATION REQUIRED TO INTERPOLATE AT CONCENTRATION
!.....RECORDING STATIONS

            global_here%N1=global_here%NM(global_here%NNC(I),1)
            global_here%N2=global_here%NM(global_here%NNC(I),2)
            global_here%N3=global_here%NM(global_here%NNC(I),3)
            global_here%X1=global_here%X(global_here%N1)
            global_here%X2=global_here%X(global_here%N2)
            global_here%X3=global_here%X(global_here%N3)
            global_here%X4=global_here%XEL(I)
            global_here%Y1=global_here%Y(global_here%N1)
            global_here%Y2=global_here%Y(global_here%N2)
            global_here%Y3=global_here%Y(global_here%N3)
            global_here%Y4=global_here%YEL(I)
            global_here%STAIC1(I)=((global_here%X4-global_here%X3)*(global_here%Y2-global_here%Y3)+(global_here%X2-global_here%X3)*(global_here%Y3-global_here%Y4))/global_here%AREAS(global_here%NNC(I))
            global_here%STAIC2(I)=((global_here%X4-global_here%X1)*(global_here%Y3-global_here%Y1)-(global_here%Y4-global_here%Y1)*(global_here%X3-global_here%X1))/global_here%AREAS(global_here%NNC(I))
            global_here%STAIC3(I)=(-(global_here%X4-global_here%X1)*(global_here%Y2-global_here%Y1)+(global_here%Y4-global_here%Y1)*(global_here%X2-global_here%X1))/global_here%AREAS(global_here%NNC(I))

         END DO
      ENDIF
      IF (global_here%NSTAC.EQ.0) s%MNSTAC = 1

!...  
!...  IF METEOROLOICAL FORCING IS INCLUDED IN THE RUN, INPUT INFORMATION FOR MET
!...  RECORDING STATIONS - OUTPUT
!...  
      global_here%NOUTM=0
      global_here%NSTAM = 0
!     
      IF(global_here%NWS.NE.0) THEN

!.....READ IN global_here%NOUTM,global_here%TOUTSM,global_here%TOUTFM,global_here%NSPOOLM : IF global_here%NOUTM<>0,INTERPOLATED
!.....MET DATA ARE SPOOLED TO UNITS 71&72 EVERY global_here%NSPOOLM TIME STEPS
!.....BETWEEN TIMES global_here%TOUTSM AND global_here%TOUTFM; IF ABS(global_here%NOUTM)=2, OUTPUT WILL BE BINARY

         READ(15,*) global_here%NOUTM,global_here%TOUTSM,global_here%TOUTFM,global_here%NSPOOLM
         WRITE(16,3211) global_here%NOUTM
 3211    FORMAT(///,1X,'METEOROLOGICAL RECORDING STATION OUTPUT : ',&
        //,5X,'global_here%NOUTM = ',I2)

!.....CHECK INPUT PARAMETER global_here%NOUTM

         IF(ABS(global_here%NOUTM).GT.2) THEN
            IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,3212)
            WRITE(16,3202)
 3212       FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',&
           //,1X,'YOUR SELECTION OF THE UNIT 15 INPUT PARAMETER',&
           ' global_here%NOUTC',&
           /,1X,'IS NOT AN ALLOWABLE VALUE.  CHECK YOUR INPUT!!')
            IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9973)
            WRITE(16,9973)
            STOP
         ENDIF

!.....IF STATION METEOROLOGICAL OUTPUT WILL NOT BE GENERATED

         IF(global_here%NOUTM.EQ.0) THEN
            WRITE(16,3213)
 3213       FORMAT(/,5X,'NO OUTPUT WILL BE SPOOLED AT METEOROLOGICAL',&
           ' RECORDING STATIONS')
         ENDIF

!.....IF STATION MET OUTPUT WILL BE GENERATED

         IF(global_here%NOUTM.NE.0) THEN

!.......COMPUTE global_here%NTCYSM, global_here%NTCYFM, WHICH = global_here%TOUTSM AND global_here%TOUTFM IN TIMESTEPS

            global_here%NTCYSM=INT((global_here%TOUTSM-global_here%STATIM)*(86400.D0/global_here%DTDP) + 0.5d0)
            global_here%NTCYFM=INT((global_here%TOUTFM-global_here%STATIM)*(86400.D0/global_here%DTDP) + 0.5d0)
            IF(global_here%NTCYFM.GT.global_here%NT) global_here%NTCYFM=global_here%NT

!.......COMPUTE global_here%NTRSPM = THE NO. OF DATA SETS TO BE SPOOLED TO UNITS 71&72

            IF(global_here%NSPOOLM.EQ.0) global_here%NTRSPM=0
            IF(global_here%NSPOOLM.NE.0) global_here%NTRSPM=INT((global_here%NTCYFM-global_here%NTCYSM)/global_here%NSPOOLM)

!.......WRITE global_here%TOUTSM,global_here%TOUTFM,global_here%NTCYSM,global_here%NTCYFM,global_here%NSPOOLM TO UNIT 16

            WRITE(16,3214) global_here%TOUTSM,global_here%NTCYSM,global_here%TOUTFM,global_here%NTCYFM,global_here%NSPOOLM
 3214       FORMAT(/,5X,'DATA RECORDS WILL START AFTER global_here%TOUTSM =',F8.3,&
           ' global_here%DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',&
           I9,' TIME STEPS INTO THE SIMULATION',&
           //,5X,'DATA RECORDS WILL STOP AFTER global_here%TOUTFM =',F8.3,&
           ' global_here%DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',&
           I9,' TIME STEPS INTO THE SIMULATION',&
           //,5X,'INFORMATION WILL BE SPOOLED TO UNITS 71&72',&
           ' EVERY global_here%NSPOOLM =',I8,' TIME STEPS')
            IF(ABS(global_here%NOUTM).EQ.1) WRITE(16,3215)
 3215       FORMAT(/,5X,'UNITS 71&72 FORMAT WILL BE ASCII')
            IF(ABS(global_here%NOUTM).EQ.2) WRITE(16,3216)
 3216       FORMAT(/,5X,'UNITS 71&72 FORMAT WILL BE BINARY')
         ENDIF

!.....REGARDLESS OF WHETHER global_here%NOUTM=0, READ IN THE NUMBER OF METEOROLOGICAL
!.....RECORDING STATIONS

         READ(15,*) global_here%NSTAM
         WRITE(16,3217) global_here%NSTAM
 3217    FORMAT(///,5X,'NUMBER OF INPUT METEOROLOGICAL RECORDING ',&
        'STATIONS = ',I5)

         IF(global_here%NSTAM.GT.0) THEN
            IF(global_here%ICS.EQ.1) WRITE(16,3218)
 3218       FORMAT(/,7X,'STATION #   ELEMENT',9X,'global_here%X',13X,'global_here%Y',/)
            IF(global_here%ICS.EQ.2) WRITE(16,3219)
 3219       FORMAT(/,5X,'STATION   ELEMENT',3X,'LAMBDA(DEG)',&
           4X,'FEA(DEG)',10X,'XCP',12X,'YCP',/)
            s%MNSTAM = global_here%NSTAM
         ENDIF

!     Allocate arrays dimensioned by MNSTAM
         call alloc_main10(s,global_here)

!.....INPUT COORDINATES OF METEOROLOGICAL RECORDING STATIONS
!.....THEN COMPUTE ELEMENT NO. WITHIN WHICH STATION LIES

         DO I=1,global_here%NSTAM
            global_here%NNM(I)=0
            IF(global_here%ICS.EQ.1) THEN
               READ(15,*) global_here%XEM(I),global_here%YEM(I)
            ELSE
               READ(15,*) global_here%SLEM(I),global_here%SFEM(I)
               global_here%SLEM(I)=global_here%SLEM(I)*DEG2RAD
               global_here%SFEM(I)=global_here%SFEM(I)*DEG2RAD
               CALL CPP(global_here%XEM(I),global_here%YEM(I),global_here%SLEM(I),global_here%SFEM(I),global_here%SLAM0,global_here%SFEA0)
            ENDIF
            global_here%AEMIN=1.0E+25
            global_here%KMIN=0
            DO K=1,global_here%NE
               global_here%N1=global_here%NM(K,1)
               global_here%N2=global_here%NM(K,2)
               global_here%N3=global_here%NM(K,3)
               global_here%X1=global_here%X(global_here%N1)
               global_here%X2=global_here%X(global_here%N2)
               global_here%X3=global_here%X(global_here%N3)
               global_here%X4=global_here%XEM(I)
               global_here%Y1=global_here%Y(global_here%N1)
               global_here%Y2=global_here%Y(global_here%N2)
               global_here%Y3=global_here%Y(global_here%N3)
               global_here%Y4=global_here%YEM(I)
               global_here%A1=(global_here%X4-global_here%X3)*(global_here%Y2-global_here%Y3)+(global_here%X2-global_here%X3)*(global_here%Y3-global_here%Y4)
               global_here%A2=(global_here%X4-global_here%X1)*(global_here%Y3-global_here%Y1)-(global_here%Y4-global_here%Y1)*(global_here%X3-global_here%X1)
               global_here%A3=(global_here%Y4-global_here%Y1)*(global_here%X2-global_here%X1)-(global_here%X4-global_here%X1)*(global_here%Y2-global_here%Y1)
               global_here%AA=ABS(global_here%A1)+ABS(global_here%A2)+ABS(global_here%A3)
               global_here%AE=ABS(global_here%AA-global_here%AREAS(K))/global_here%AREAS(K)
               IF(global_here%AE.LT.global_here%AEMIN) THEN
                  global_here%AEMIN=global_here%AE
                  global_here%KMIN=K
               ENDIF
               IF(global_here%AE.LT.1.0E-5) global_here%NNM(I)=K
            END DO

            IF(global_here%NNM(I).EQ.0) THEN
               IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9942) I
               WRITE(16,9942) I
 9942          FORMAT(///,1X,'!!!!!!!!!!  WARNING - NONFATAL INPUT ERROR ',&
          '!!!!!!!!!',//,&
          ' METEOROLOGICAL RECORDING STATION ',I6,' DOES NOT LIE'&
          ,' WITHIN ANY ELEMENT IN THE DEFINED',/,&
          ' COMPUTATIONAL DOMAIN,   PLEASE CHECK THE INPUT',&
          ' COORDINATES FOR THIS STATION')
               IF(global_here%NFOVER.EQ.1) THEN
                  IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9790) global_here%AEMIN
                  WRITE(16,9790) global_here%AEMIN
                  global_here%NNM(I)=global_here%KMIN
               ELSE
                  IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9791) global_here%AEMIN
                  WRITE(16,9791) global_here%AEMIN
                  STOP
               ENDIF
            ENDIF

            IF(global_here%ICS.EQ.1) THEN
               WRITE(16,1880) I,global_here%NNM(I),global_here%XEM(I),global_here%YEM(I)
            ELSE
               WRITE(16,1883) I,global_here%NNM(I),global_here%SLEM(I)*RAD2DEG,&
              global_here%SFEM(I)*RAD2DEG,global_here%XEM(I),global_here%YEM(I)
            ENDIF

!.....PRE-COMPUTE INFORMATION REQUIRED TO INTERPOLATE AT METEOROLOGICAL
!.....RECORDING STATIONS

            global_here%N1=global_here%NM(global_here%NNM(I),1)
            global_here%N2=global_here%NM(global_here%NNM(I),2)
            global_here%N3=global_here%NM(global_here%NNM(I),3)
            global_here%X1=global_here%X(global_here%N1)
            global_here%X2=global_here%X(global_here%N2)
            global_here%X3=global_here%X(global_here%N3)
            global_here%X4=global_here%XEM(I)
            global_here%Y1=global_here%Y(global_here%N1)
            global_here%Y2=global_here%Y(global_here%N2)
            global_here%Y3=global_here%Y(global_here%N3)
            global_here%Y4=global_here%YEM(I)
            global_here%STAIM1(I)=((global_here%X4-global_here%X3)*(global_here%Y2-global_here%Y3)+(global_here%X2-global_here%X3)*(global_here%Y3-global_here%Y4))/global_here%AREAS(global_here%NNM(I))
            global_here%STAIM2(I)=((global_here%X4-global_here%X1)*(global_here%Y3-global_here%Y1)-(global_here%Y4-global_here%Y1)*(global_here%X3-global_here%X1))/global_here%AREAS(global_here%NNM(I))
            global_here%STAIM3(I)=(-(global_here%X4-global_here%X1)*(global_here%Y2-global_here%Y1)+(global_here%Y4-global_here%Y1)*(global_here%X2-global_here%X1))/global_here%AREAS(global_here%NNM(I))

         END DO
      ENDIF
      IF (global_here%NSTAM.EQ.0) s%MNSTAM = 1

!...  
!...  INPUT INFORMATION ABOUT GLOBAL ELEVATION DATA OUTPUT
!...  

!.... READ IN global_here%NOUTGE,global_here%TOUTSGE,global_here%TOUTFGE,global_here%NSPOOLGE : IF global_here%NOUTGE<>0, GLOBAL ELEV.
!.... OUTPUT IS SPOOLED TO UNIT 63 EVERY global_here%NSPOOLGE TIME STEPS BETWEEN
!.... TIMES global_here%TOUTSGE AND global_here%TOUTFGE; IF ABS(global_here%NOUTGE)=2, OUTPUT WILL BE BINARY

      READ(15,*) global_here%NOUTGE,global_here%TOUTSGE,global_here%TOUTFGE,global_here%NSPOOLGE
      WRITE(16,3301) global_here%NOUTGE
 3301 FORMAT(////,1X,'GLOBAL NODAL ELEVATION INFORMATION OUTPUT: ',&
     //,5X,'global_here%NOUTGE = ',I2)

!.... CHECK INPUT PARAMETER global_here%NOUTGE

      IF(ABS(global_here%NOUTGE).GT.2) THEN
         IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,3302)
         WRITE(16,3302)
 3302    FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',&
        //,1X,'YOUR SELECTION OF THE UNIT 15 INPUT PARAMETER',&
        ' global_here%NOUTGE',&
        /,1X,'IS NOT AN ALLOWABLE VALUE.  CHECK YOUR INPUT!!')
         IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9973)
         WRITE(16,9973)
         STOP
      ENDIF

!.... IF GLOBAL ELEVATION OUTPUT WILL NOT BE GENERATED

      IF(global_here%NOUTGE.EQ.0) THEN
         WRITE(16,3303)
 3303    FORMAT(///,5X,'NO GLOBAL ELEVATION OUTPUT WILL BE SPOOLED')
      ENDIF

!.... IF GLOBAL ELEVATION OUTPUT WILL BE GENERATED

      IF(global_here%NOUTGE.NE.0) THEN

!......COMPUTE global_here%NTCYSGE, global_here%NTCYFGE, WHICH = global_here%TOUTSGE AND global_here%TOUTFGE IN TIMESTEPS

         global_here%NTCYSGE=INT((global_here%TOUTSGE-global_here%STATIM)*(86400.D0/global_here%DTDP) + 0.5d0)
         global_here%NTCYFGE=INT((global_here%TOUTFGE-global_here%STATIM)*(86400.D0/global_here%DTDP) + 0.5d0)
         IF(global_here%NTCYFGE.GT.global_here%NT) global_here%NTCYFGE=global_here%NT

!......CALCULATE global_here%NDSETSE = THE # OF DATA SETS TO BE SPOOLED TO UNIT 63

         IF(global_here%NSPOOLGE.EQ.0) global_here%NDSETSE=0
         IF(global_here%NSPOOLGE.NE.0) global_here%NDSETSE=INT((global_here%NTCYFGE-global_here%NTCYSGE)/global_here%NSPOOLGE)

!......WRITE global_here%NOUTGE,global_here%TOUTSGE,global_here%TOUTFGE,global_here%NTCYSGE,global_here%NTCYFGE,global_here%NSPOOLGE TO UNIT 16

         WRITE(16,3304) global_here%TOUTSGE,global_here%NTCYSGE,global_here%TOUTFGE,global_here%NTCYFGE,global_here%NSPOOLGE
 3304    FORMAT(/,5X,'DATA RECORDS WILL START AFTER global_here%TOUTSGE =',F8.3,&
        ' global_here%DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',&
        I9,' TIME STEPS INTO THE SIMULATION',&
        //,5X,'DATA RECORDS WILL STOP AFTER global_here%TOUTFGE =',F8.3,&
        ' global_here%DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',&
        I9,' TIME STEPS INTO THE SIMULATION',&
        //,5X,'INFORMATION WILL BE SPOOLED TO UNIT 63 EVERY ',&
        'global_here%NSPOOLGE =',I8,' TIME STEPS')
         IF(ABS(global_here%NOUTGE).EQ.1) WRITE(16,3305)
 3305    FORMAT(/,5X,'UNIT 63 FORMAT WILL BE ASCII')
         IF(ABS(global_here%NOUTGE).EQ.2) WRITE(16,3306)
 3306    FORMAT(/,5X,'UNIT 63 FORMAT WILL BE BINARY')
      ENDIF

!...  
!...  INPUT INFORMATION ABOUT GLOBAL VELOCITY DATA OUTPUT
!...  

!.... READ IN global_here%NOUTGV,global_here%TOUTSGV,global_here%TOUTFGV,global_here%NSPOOLGV : IF global_here%NOUTGV<>0, GLOBAL global_here%VEL.
!.... OUTPUT IS SPOOLED TO UNIT 64 EVERY global_here%NSPOOLGV TIME STEPS BETWEEN
!.... TIMES global_here%TOUTSGV AND global_here%TOUTFGV; IF ABS(global_here%NOUTGV)=2, OUTPUT WILL BE BINARY

      READ(15,*) global_here%NOUTGV,global_here%TOUTSGV,global_here%TOUTFGV,global_here%NSPOOLGV
      WRITE(16,3351) global_here%NOUTGV
 3351 FORMAT(////,1X,'GLOBAL NODAL VELOCITY INFORMATION OUTPUT : ',&
     //,5X,'global_here%NOUTGV = ',I2)

!.... CHECK INPUT PARAMETER global_here%NOUTGV

      IF(ABS(global_here%NOUTGV).GT.2) THEN
         IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,3352)
         WRITE(16,3352)
 3352    FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',&
        //,1X,'YOUR SELECTION OF THE UNIT 15 INPUT PARAMETER',&
        ' global_here%NOUTGV',&
        /,1X,'IS NOT AN ALLOWABLE VALUE.  CHECK YOUR INPUT!!')
         IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9973)
         WRITE(16,9973)
         STOP
      ENDIF

!.... IF GLOBAL VELOCITY OUTPUT WILL NOT BE GENERATED

      IF(global_here%NOUTGV.EQ.0) THEN
         WRITE(16,3353)
 3353    FORMAT(///,5X,'NO GLOBAL VELOCITY OUTPUT WILL BE SPOOLED')
      ENDIF

!.... IF GLOBAL VELOCITY OUTPUT WILL BE GENERATED

      IF(global_here%NOUTGV.NE.0) THEN

!......COMPUTE global_here%NTCYSGV, global_here%NTCYFGV, WHICH = global_here%TOUTSGV AND global_here%TOUTFGV IN TIMESTEPS

         global_here%NTCYSGV=INT((global_here%TOUTSGV-global_here%STATIM)*(86400.D0/global_here%DTDP) + 0.5d0)
         global_here%NTCYFGV=INT((global_here%TOUTFGV-global_here%STATIM)*(86400.D0/global_here%DTDP) + 0.5d0)
         IF(global_here%NTCYFGV.GT.global_here%NT) global_here%NTCYFGV=global_here%NT

!......CALCULATE global_here%NDSETSV = THE # OF DATA SETS TO BE SPOOLED TO UNIT 64

         IF(global_here%NSPOOLGV.EQ.0) global_here%NDSETSV=0
         IF(global_here%NSPOOLGV.NE.0) global_here%NDSETSV=INT((global_here%NTCYFGV-global_here%NTCYSGV)/global_here%NSPOOLGV)

!......WRITE global_here%NOUTGV,global_here%TOUTSGV,global_here%TOUTFGV,global_here%NTCYSGV,global_here%NTCYFGV,global_here%NSPOOLGV TO UNIT 16

         WRITE(16,3354) global_here%TOUTSGV,global_here%NTCYSGV,global_here%TOUTFGV,global_here%NTCYFGV,global_here%NSPOOLGV
 3354    FORMAT(/,5X,'DATA RECORDS WILL START AFTER global_here%TOUTSGV =',F8.3,&
        ' global_here%DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',&
        I9,' TIME STEPS INTO THE SIMULATION',&
        //,5X,'DATA RECORDS WILL STOP AFTER global_here%TOUTFGV =',F8.3,&
        ' global_here%DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',&
        I9,' TIME STEPS INTO THE SIMULATION',&
        //,5X,'INFORMATION WILL BE SPOOLED TO UNIT 64 EVERY ',&
        'global_here%NSPOOLGV =',I8,' TIME STEPS')
         IF(ABS(global_here%NOUTGV).EQ.1) WRITE(16,3355)
 3355    FORMAT(/,5X,'UNIT 64 FORMAT WILL BE ASCII')
         IF(ABS(global_here%NOUTGV).EQ.2) WRITE(16,3356)
 3356    FORMAT(/,5X,'UNIT 64 FORMAT WILL BE BINARY')
      ENDIF

!...  
!...  IF TRANSPORT IS INCLUDED IN THE RUN, INPUT INFORMATION ABOUT GLOBAL
!...  CONCENTRATION DATA OUTPUT
!...  
      global_here%NOUTGC=0
      IF(global_here%IM.EQ.10) THEN

!.....READ IN global_here%NOUTGC,global_here%TOUTSGC,global_here%TOUTFGC,global_here%NSPOOLGC : IF global_here%NOUTGC<>0, GLOBAL
!.....CONCENTRATION OUTPUT IS SPOOLED TO UNIT 73 EVERY global_here%NSPOOLGC TIME STEPS
!.....BETWEEN TIMES global_here%TOUTSGC AND global_here%TOUTFGC; IF ABS(global_here%NOUTGC)=2, OUTPUT WILL BE BINARY

         READ(15,*) global_here%NOUTGC,global_here%TOUTSGC,global_here%TOUTFGC,global_here%NSPOOLGC
         WRITE(16,3401) global_here%NOUTGC
 3401    FORMAT(////,1X,'GLOBAL NODAL CONCENTRATION INFORMATION OUTPUT:',&
        //,5X,'global_here%NOUTGC = ',I2)

!.....CHECK INPUT PARAMETER global_here%NOUTGC

         IF(ABS(global_here%NOUTGC).GT.2) THEN
            IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,3402)
            WRITE(16,3402)
 3402       FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',&
           //,1X,'YOUR SELECTION OF THE UNIT 15 INPUT PARAMETER',&
           ' global_here%NOUTGC',&
           /,1X,'IS NOT AN ALLOWABLE VALUE.  CHECK YOUR INPUT!!')
            IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9973)
            WRITE(16,9973)
            STOP
         ENDIF

!.....IF GLOBAL CONCENTRATION OUTPUT WILL NOT BE GENERATED

         IF(global_here%NOUTGC.EQ.0) THEN
            WRITE(16,3403)
 3403       FORMAT(///,5X,'NO GLOBAL CONCENTRATION OUTPUT WILL BE ',&
           'SPOOLED')
         ENDIF

!.....IF GLOBAL CONCENTRATION OUTPUT WILL BE GENERATED

         IF(global_here%NOUTGC.NE.0) THEN

!.......COMPUTE global_here%NTCYSGC, global_here%NTCYFGC, WHICH = global_here%TOUTSGC AND global_here%TOUTFGC IN TIMESTEPS

            global_here%NTCYSGC=INT((global_here%TOUTSGC-global_here%STATIM)*(86400.D0/global_here%DTDP) + 0.5d0)
            global_here%NTCYFGC=INT((global_here%TOUTFGC-global_here%STATIM)*(86400.D0/global_here%DTDP) + 0.5d0)
            IF(global_here%NTCYFGC.GT.global_here%NT) global_here%NTCYFGC=global_here%NT

!.......CALCULATE global_here%NDSETSC = THE # OF DATA SETS TO BE SPOOLED TO UNIT 73

            IF(global_here%NSPOOLGC.EQ.0) global_here%NDSETSC=0
            IF(global_here%NSPOOLGC.NE.0) global_here%NDSETSC=INT((global_here%NTCYFGC-global_here%NTCYSGC)/global_here%NSPOOLGC)

!.......WRITE global_here%NOUTGC,global_here%TOUTSGC,global_here%TOUTFGC,global_here%NTCYSGC,global_here%NTCYFGC,global_here%NSPOOLGC TO UNIT 16

            WRITE(16,3404) global_here%TOUTSGC,global_here%NTCYSGC,global_here%TOUTFGC,global_here%NTCYFGC,global_here%NSPOOLGC
 3404       FORMAT(/,5X,'DATA RECORDS WILL START AFTER global_here%TOUTSGC =',F8.3,&
           ' global_here%DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',&
           I9,' TIME STEPS INTO THE SIMULATION',&
           //,5X,'DATA RECORDS WILL STOP AFTER global_here%TOUTFGC =',F8.3,&
           ' global_here%DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',&
           I9,' TIME STEPS INTO THE SIMULATION',&
           //,5X,'INFORMATION WILL BE SPOOLED TO UNIT 73 EVERY ',&
           'global_here%NSPOOLGC =',I8,' TIME STEPS')
            IF(ABS(global_here%NOUTGC).EQ.1) WRITE(16,3405)
 3405       FORMAT(/,5X,'UNIT 73 FORMAT WILL BE ASCII')
            IF(ABS(global_here%NOUTGC).EQ.2) WRITE(16,3406)
 3406       FORMAT(/,5X,'UNIT 73 FORMAT WILL BE BINARY')
         ENDIF

      ENDIF

!...  
!...  IF global_here%NWS<>0   INPUT INFORMATION ABOUT GLOBAL WIND DATA OUTPUT
!...  
      IF(global_here%NWS.NE.0) THEN

!......READ IN global_here%NOUTGW,global_here%TOUTSGW,global_here%TOUTFGW,global_here%NSPOOLGW : IF global_here%NOUTGW<>0, GLOBAL WIND
!......OUTPUT IS SPOOLED TO UNIT 74 EVERY global_here%NSPOOLGW TIME STEPS BETWEEN
!......TIMES global_here%TOUTSGW AND global_here%TOUTFGW; IF ABS(global_here%NOUTGW)=2, OUTPUT WILL BE BINARY

         READ(15,*) global_here%NOUTGW,global_here%TOUTSGW,global_here%TOUTFGW,global_here%NSPOOLGW
         WRITE(16,3451) global_here%NOUTGW
 3451    FORMAT(////,1X,'GLOBAL WIND STRESS INFORMATION OUTPUT : ',&
        //,5X,'global_here%NOUTGW = ',I2)

!......CHECK INPUT PARAMETER global_here%NOUTGW

         IF(ABS(global_here%NOUTGW).GT.2) THEN
            IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,3452)
            WRITE(16,3452)
 3452       FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',&
           //,1X,'YOUR SELECTION OF THE UNIT 15 INPUT PARAMETER',&
           ' global_here%NOUTGW',&
           /,1X,'IS NOT AN ALLOWABLE VALUE.  CHECK YOUR INPUT!!')
            IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,3453)
            WRITE(16,3453)
            STOP
         ENDIF

!......IF GLOBAL WIND STRESS OUTPUT WILL NOT BE GENERATED

         IF(global_here%NOUTGW.EQ.0) THEN
            WRITE(16,3453)
 3453       FORMAT(///,5X,'NO GLOBAL WIND STRESS OUTPUT WILL BE SPOOLED')
         ENDIF

!........IF GLOBAL WIND STRESS OUTPUT WILL BE GENERATED

         IF(global_here%NOUTGW.NE.0) THEN

!........COMPUTE global_here%NTCYSGW, global_here%NTCYFGW, WHICH = global_here%TOUTSGW AND global_here%TOUTFGW IN TIMESTEPS

            global_here%NTCYSGW=INT((global_here%TOUTSGW-global_here%STATIM)*(86400.D0/global_here%DTDP) + 0.5d0)
            global_here%NTCYFGW=INT((global_here%TOUTFGW-global_here%STATIM)*(86400.D0/global_here%DTDP) + 0.5d0)
            IF(global_here%NTCYFGW.GT.global_here%NT) global_here%NTCYFGW=global_here%NT

!........CALCULATE global_here%NDSETSW = THE # OF DATA SETS TO BE SPOOLED TO UNIT 74

            IF(global_here%NSPOOLGW.EQ.0) global_here%NDSETSW=0
            IF(global_here%NSPOOLGW.NE.0) global_here%NDSETSW=INT((global_here%NTCYFGW-global_here%NTCYSGW)/global_here%NSPOOLGW)

!........WRITE global_here%NOUTGW,global_here%TOUTSGW,global_here%TOUTFGW,global_here%NTCYSGW,global_here%NTCYFGW,global_here%NSPOOLGW TO UNIT 16

            WRITE(16,3454) global_here%TOUTSGW,global_here%NTCYSGW,global_here%TOUTFGW,global_here%NTCYFGW,global_here%NSPOOLGW
 3454       FORMAT(/,5X,'DATA RECORDS WILL START AFTER global_here%TOUTSGW =',F8.3,&
           ' global_here%DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',&
           I9,' TIME STEPS INTO THE SIMULATION',&
           //,5X,'DATA RECORDS WILL STOP AFTER global_here%TOUTFGW =',F8.3,&
           ' global_here%DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',&
           I9,' TIME STEPS INTO THE SIMULATION',&
           //,5X,'INFORMATION WILL BE SPOOLED TO UNIT 74 EVERY ',&
           'global_here%NSPOOLGW =',I8,' TIME STEPS')
            IF(ABS(global_here%NOUTGW).EQ.1) WRITE(16,3455)
 3455       FORMAT(/,5X,'UNIT 74 FORMAT WILL BE ASCII')
            IF(ABS(global_here%NOUTGW).EQ.2) WRITE(16,3456)
 3456       FORMAT(/,5X,'UNIT 74 FORMAT WILL BE BINARY')
         ENDIF

      ENDIF

!...  
!...  READ AND CHECK INFORMATION ABOUT HARMONIC ANALYSIS OF MODEL RESULTS
!...  
#ifdef harm
      READ(15,*) NFREQ 
      WRITE(16,99392) NFREQ  
99392 FORMAT(////,1X,'HARMONIC ANALYSIS INFORMATION OUTPUT : ',&
     //,5X,'HARMONIC ANALYSIS PERFORMED FOR ',I4,' CONSTITUENTS',/)
      s%MNHARF = NFREQ

      IF (NFREQ.EQ.0) s%MNHARF = 1

!     Allocate harmonic analysis arrays

      IF (NFREQ.GT.0) THEN
         CALL ALLOC_HA()
         CALL ALLOC_MAIN14(s,global_here)
      ENDIF

      IF(NFREQ.LT.0) THEN
         IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,99391)
         WRITE(16,99391)
99391    FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',&
        //,1X,'YOUR SELECTION OF global_here%NHARFR (A UNIT 15 '&
        ,'INPUT PARAMETER) IS NOT AN ALLOWABLE VALUE',/,1X,&
        'PLEASE CHECK YOUR INPUT',&
        //,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
         STOP
      ENDIF      
      IF(NFREQ.GT.0) WRITE(16,2330)
 2330 FORMAT(/,7X,'FREQUENCY',4X,'NODAL FACTOR',6X,'EQU.global_here%ARG(DEG)',&
     1X,'CONSTITUENT',/)
      DO 1201 I=1,NFREQ  
         READ(15,'(A10)') NAMEFR(I)
         READ(15,*) HAFREQ(I),HAFF(I),HAFACE(I)
         WRITE(16,2331) HAFREQ(I),HAFF(I),HAFACE(I),NAMEFR(I)
#else
      READ(15,*) nfreq_dummy
      IF (nfreq_dummy.EQ.0) s%MNHARF = 1
      if (nfreq_dummy.ne.0) then
         write(16,*) "HARM not supported. Stopping execution."
         stop
      end if
#endif 
 2331    FORMAT(4X,F15.12,2X,F10.7,5X,F10.3,7X,A10)
 1201 CONTINUE

!.... READ IN INTERVAL INFORMATION FOR HARMONIC ANALYSIS
!.... COMPUTE global_here%THAS AND global_here%THAF IN TERMS OF THE NUMBER OF TIME STEPS
      READ(15,*) global_here%THAS,global_here%THAF,global_here%NHAINC,global_here%FMV
      global_here%ITHAS=INT((global_here%THAS-global_here%STATIM)*(86400.D0/global_here%DTDP) + 0.5d0)
      global_here%THAS=global_here%ITHAS*global_here%DTDP/86400.D0 + global_here%STATIM
      global_here%ITHAF=INT((global_here%THAF-global_here%STATIM)*(86400.D0/global_here%DTDP) + 0.5d0)
      global_here%THAF=global_here%ITHAF*global_here%DTDP/86400.D0 + global_here%STATIM
      global_here%ITMV = global_here%ITHAF - (global_here%ITHAF-global_here%ITHAS)*global_here%FMV
#ifdef HARM
      IF(NFREQ.GT.0) THEN
         WRITE(16,34634) global_here%THAS,global_here%ITHAS,global_here%THAF,global_here%ITHAF,global_here%NHAINC
34634    FORMAT(/,5X,'HARMONIC ANALYSIS WILL START AFTER global_here%THAS =',F8.3,&
        ' global_here%DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',I9,&
        ' TIME STEPS INTO THE SIMULATION',&
        //,5X,'HARMONIC ANALYSIS WILL STOP AFTER global_here%THAF =',F8.3,&
        ' global_here%DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',I9,&
        ' TIME STEPS INTO THE SIMULATION'&
        ,//,5X,'INFORMATION WILL BE ANALYZED EVERY ',&
        'global_here%NHAINC =',I8,' TIME STEPS.')
         WRITE(16,34639) global_here%FMV*100.,global_here%ITMV
34639    FORMAT(/,5X,'MEANS AND VARIANCES WILL BE COMPUTED FOR THE ',&
        'FINAL ',F10.5,' %',/9X,'OF THE HARMONIC ANALYSIS ',&
        'PERIOD OR AFTER ',I9,' TIME STEPS INTO THE ',&
        'SIMULATION.',/9X,' RESULTS ARE WRITTEN TO UNIT 55.')

      ELSE
#endif
         WRITE(16,34645)
34645    FORMAT(///,5X,'NO HARMONIC ANALYSIS WILL BE DONE')
#ifdef HARM
      ENDIF
      IF((global_here%FMV.GT.0.).AND.(NFREQ.GT.0).AND.(s%C2DDI)) s%CHARMV = .TRUE.
#endif

!.... READ IN AND WRITE OUT INFORMATION ON WHERE HARMONIC ANALYSIS WILL BE DONE

      READ(15,*) global_here%NHASE,global_here%NHASV,global_here%NHAGE,global_here%NHAGV
      IF((global_here%NHASE.LT.0).OR.(global_here%NHASE.GT.1)) THEN
         IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,99661)
         WRITE(16,99661)
99661    FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',&
        'INPUT ERROR  !!!!!!!!!',//&
        ,1X,'YOUR SELECTION OF global_here%NHASE (A UNIT 15 '&
        ,'INPUT PARAMETER) IS NOT AN ALLOWABLE VALUE',/,1X,&
        'PLEASE CHECK YOUR INPUT')
         IF(global_here%NFOVER.EQ.1) THEN
            IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,99671)
            WRITE(16,99671)
99671       FORMAT(/,1X,'PROGRAM WILL OVERRIDE SPECIFIED INPUT',&
           ' AND SET global_here%NHASE EQUAL TO 0 ',&
           //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)
            global_here%NHASE=0
         ELSE
            IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9973)
            WRITE(16,9973)
            STOP
         ENDIF
      ENDIF
      IF(global_here%NHASE.EQ.1) THEN
         WRITE(16,34641)
34641    FORMAT(///,5X,'STATION ELEVATION HARMONIC ANAL WILL BE ',&
        'WRITTEN TO UNIT 51')
      ENDIF
      IF((global_here%NHASV.LT.0).OR.(global_here%NHASV.GT.1)) THEN
         IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,99662)
         WRITE(16,99662)
99662    FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',&
        'INPUT ERROR  !!!!!!!!!',//&
        ,1X,'YOUR SELECTION OF global_here%NHASV (A UNIT 15 '&
        ,'INPUT PARAMETER) IS NOT AN ALLOWABLE VALUE',/,1X,&
        'PLEASE CHECK YOUR INPUT')
         IF(global_here%NFOVER.EQ.1) THEN
            IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,99672)
            WRITE(16,99672)
99672       FORMAT(/,1X,'PROGRAM WILL OVERRIDE SPECIFIED INPUT',&
           ' AND SET global_here%NHASV EQUAL TO 0 ',&
           //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)
            global_here%NHASV=0
         ELSE
            IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9973)
            WRITE(16,9973)
            STOP
         ENDIF
      ENDIF
      IF(global_here%NHASV.EQ.1) THEN
         WRITE(16,34642)
34642    FORMAT(///,5X,'STATION VELOCITY HARMONIC ANAL WILL BE ',&
        'WRITTEN TO UNIT 52')
      ENDIF
      IF((global_here%NHAGE.LT.0).OR.(global_here%NHAGE.GT.1)) THEN
         IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,99663)
         WRITE(16,99663)
99663    FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',&
        'INPUT ERROR  !!!!!!!!!',//&
        ,1X,'YOUR SELECTION OF global_here%NHAGE (A UNIT 15 '&
        ,'INPUT PARAMETER) IS NOT AN ALLOWABLE VALUE',/,1X,&
        'PLEASE CHECK YOUR INPUT')
         IF(global_here%NFOVER.EQ.1) THEN
            IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,99673)
            WRITE(16,99673)
99673       FORMAT(/,1X,'PROGRAM WILL OVERRIDE SPECIFIED INPUT',&
           ' AND SET global_here%NHAGE EQUAL TO 0 ',&
           //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)
            global_here%NHAGE=0
         ELSE
            IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9973)
            WRITE(16,9973)
            STOP
         ENDIF
      ENDIF
      IF(global_here%NHAGE.EQ.1) THEN
         WRITE(16,34643)
34643    FORMAT(///,5X,'GLOBAL ELEVATION HARMONIC ANAL WILL BE ',&
              'WRITTEN TO UNIT 53')
      ENDIF
      IF((global_here%NHAGV.LT.0).OR.(global_here%NHAGV.GT.1)) THEN
         IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,99664)
         WRITE(16,99664)
99664    FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',&
        'INPUT ERROR  !!!!!!!!!',//&
        ,1X,'YOUR SELECTION OF global_here%NHAGV (A UNIT 15 '&
        ,'INPUT PARAMETER) IS NOT AN ALLOWABLE VALUE',/,1X,&
        'PLEASE CHECK YOUR INPUT')
         IF(global_here%NFOVER.EQ.1) THEN
            IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,99674)
            WRITE(16,99674)
99674       FORMAT(/,1X,'PROGRAM WILL OVERRIDE SPECIFIED INPUT',&
           ' AND SET global_here%NHAGV EQUAL TO 0 ',&
           //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)
            global_here%NHAGV=0
         ELSE
            IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9973)
            WRITE(16,9973)
            STOP
         ENDIF
      ENDIF
      IF(global_here%NHAGV.EQ.1) THEN
         WRITE(16,34644)
34644    FORMAT(///,5X,'GLOBAL VELOCITY HARMONIC ANAL WILL BE ',&
        'WRITTEN TO UNIT 54')
      ENDIF
      
!.... ESTABLISH INDICATOR OF WHETHER ANY HARMONIC ANALYSIS WILL BE DONE

#ifdef HARM
      global_here%IHARIND=NFREQ*(global_here%NHASE+global_here%NHASV+global_here%NHAGE+global_here%NHAGV)
      IF(global_here%IHARIND.GT.0) global_here%IHARIND=1
#endif

!...  
!...  INPUT INFORMATION ABOUT HOT START OUTPUT
!...  
      READ(15,*) global_here%NHSTAR,global_here%NHSINC
      WRITE(16,99655)
99655 FORMAT(////,1X,'HOT START OUTPUT INFORMATION OUTPUT : ')
      IF((global_here%NHSTAR.LT.0).OR.(global_here%NHSTAR.GT.1)) THEN
         IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,99665)
         WRITE(16,99665)
99665    FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL ',&
        'INPUT ERROR  !!!!!!!!!',//&
        ,1X,'YOUR SELECTION OF global_here%NHSTAR (A UNIT 15 '&
        ,'INPUT PARAMETER) IS NOT AN ALLOWABLE VALUE',/,1X,&
        'PLEASE CHECK YOUR INPUT')
         IF(global_here%NFOVER.EQ.1) THEN
            IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,99675)
            WRITE(16,99675)
99675       FORMAT(/,1X,'PROGRAM WILL OVERRIDE SPECIFIED INPUT',&
           ' AND SET global_here%NHSTAR EQUAL TO 0 ',&
           //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)
            global_here%NHSTAR=0
         ELSE
            IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9973)
            WRITE(16,9973)
            STOP
         ENDIF
      ENDIF
      IF(global_here%NHSTAR.EQ.1) THEN
         WRITE(16,34636) global_here%NHSINC
34636    FORMAT(/,5X,'HOT START OUTPUT WILL BE WRITTEN TO UNIT',&
        ' 67 OR 68 EVERY ',I5,' TIME STEPS')
      ELSE
         WRITE(16,34646)
34646    FORMAT(///,5X,'NO HOT START OUTPUT WILL BE GENERATED')
      ENDIF
      IF((global_here%IHOT.EQ.0).OR.(global_here%IHOT.EQ.68)) global_here%IHSFIL=67
      IF(global_here%IHOT.EQ.67) global_here%IHSFIL=68
!...  
!...  INPUT INFORMATION ABOUT SOLVER
!...  

!...  THIS SECTION TO LUMP THE GWCE MATRIX
!     vjp 11/30/99 made lumping a compile time option

#ifdef LUMP
      s%CLUMP = .TRUE.
      global_here%ILUMP=1
#else
      s%CLUMP = .FALSE.
      global_here%ILUMP=0
#endif
      
      READ(15,*) global_here%ITITER,global_here%ISLDIA,global_here%CONVCR,global_here%ITMAX

      WRITE(16,99656)
99656 FORMAT(//,1X,'SOLVER INFORMATION OUTPUT : ')

!     - allocate arrays dimensioned by MNEI
      call alloc_main11(s,global_here)

!...  LINES TO USE THE ITERATIVE MATRIX SOLVER

      IF((global_here%ISLDIA.LT.0).OR.(global_here%ISLDIA.GT.5)) THEN
         IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9920)
         WRITE(16,9920)
 9920    FORMAT(////,1X,'!!!!!!!!!!  WARNING - NONFATAL INPUT ERROR',&
        ' !!!!!!!!!',//,1X,'global_here%ISLDIA (A UNIT 15 INPUT PARAMETER) ',&
        'MUST BE 0-5',/,1X,'PLEASE CHECK YOUR INPUT')
         IF(global_here%NFOVER.EQ.1) THEN
            IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9921)
            WRITE(16,9921)
 9921       FORMAT(/,1X,'PROGRAM WILL OVERRIDE SPECIFIED INPUT',&
           ' AND SET global_here%ISLDIA EQUAL TO 0 ',&
           //,1X,'!!!!!! EXECUTION WILL CONTINUE !!!!!!',//)
            global_here%ISLDIA=0
         ELSE
            IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,9973)
            WRITE(16,9973)
            STOP
         ENDIF
      ENDIF

#ifdef  CMPI
!     sb-PDG1 deleted
!     READ(15,*) MNPROC
!--   
#else
      s%MNPROC = 1
#endif

!.....INITIALIZE AVERAGING FOR INTERNAL BARRIER WATER LEVELS
!......global_here%BARAVGWT=0.000 -> NO AVERAGING PERFORMED
!jj   wm001 changed one line                
      global_here%BARAVGWT=0.000D0
      global_here%IBSTART=0
      DO I=1,global_here%NVEL
         global_here%RBARWL1AVG(I)=0.D0
         global_here%RBARWL2AVG(I)=0.D0
      END DO
!.....INITIALIZE global_here%NIBNODECODE(I)
      DO I=1,global_here%NP
         global_here%NIBNODECODE(I)=0
      END DO
!...  
!...  COMPUTE THE NEIGHBOR TABLE
!...  
      IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) WRITE(6,1196)
      WRITE(16,1196)
 1196 FORMAT(/,1X,'THE NEIGHBOR TABLE IS BEING COMPUTED ',/)
!     
      CALL NEIGHB(s,global_here%NE,global_here%NP,global_here%NM,global_here%NNEIGH,global_here%NEITAB,global_here%NEIMIN,global_here%NEIMAX,global_here%X,global_here%Y,global_here%NSCREEN,global_here%NNEIGH_ELEM,global_here%NEIGH_ELEM)
!     
      IF(global_here%NSCREEN.EQ.1.AND.s%MYPROC.EQ.0) &
     WRITE(6,1195) global_here%NEIMIN,global_here%NEIMAX,global_here%NEIMAX
      WRITE(16,1195) global_here%NEIMIN,global_here%NEIMAX,global_here%NEIMAX
 1195 FORMAT(1X,'THE NEIGHBOR TABLE IS COMPLETED ',&
     /,5X,'THE MINIMUM NUMBER OF NEIGHBORS FOR ANY NODE = ',I3,&
     /,5X,'1+THE MAXIMUM NUMBER OF NEIGHBORS FOR ANY NODE = ',I3,&
     /,5X,'THE PARAMETER MNEI CAN BE SET AS SMALL AS ',I3,/)


!.....Allocate arrays dealing with wind forcing

      CALL ALLOC_MAIN12(s,global_here)
      
!.....Allocate arrays for wave modified bottom friction (EJK)

      IF ((global_here%FRW.EQ.1).OR.(global_here%NRS.EQ.1)) THEN
         CALL ALLOC_MAIN15(s,global_here)
      ENDIF
      
!     sb-
!.....Count maximum of the number of the elements associated with a node

      DO I = 1,s%MNP
         global_here%NNDEL(I) = 0
      ENDDO
      DO IK=1,s%MNE
         global_here%NNDEL(global_here%NM(IK,1)) = global_here%NNDEL(global_here%NM(IK,1)) + 1
         global_here%NNDEL(global_here%NM(IK,2)) = global_here%NNDEL(global_here%NM(IK,2)) + 1
         global_here%NNDEL(global_here%NM(IK,3)) = global_here%NNDEL(global_here%NM(IK,3)) + 1
      ENDDO
      s%MNNDEL = 0
      DO IK=1,s%MNP
         IF(global_here%NNDEL(IK).GT.s%MNNDEL) s%MNNDEL = global_here%NNDEL(IK)
      ENDDO

!.....Allocate space for Arrays dimensioned by MNNDEL
      CALL ALLOC_MAIN16(s,global_here)
      
      DO I = 1,s%MNP
         global_here%NNDEL(I) = 0
      ENDDO

!.....Make node-to-elements table


      DO IK=1,s%MNE
         global_here%NM1 = global_here%NM(IK,1)
         global_here%NM2 = global_here%NM(IK,2)
         global_here%NM3 = global_here%NM(IK,3)
         
         global_here%NNDEL(global_here%NM1) = global_here%NNDEL(global_here%NM1) + 1
         global_here%NDEL(global_here%NM1,global_here%NNDEL(global_here%NM1)) = IK

         global_here%NNDEL(global_here%NM2) = global_here%NNDEL(global_here%NM2) + 1
         global_here%NDEL(global_here%NM2,global_here%NNDEL(global_here%NM2)) = IK

         global_here%NNDEL(global_here%NM3) = global_here%NNDEL(global_here%NM3) + 1
         global_here%NDEL(global_here%NM3,global_here%NNDEL(global_here%NM3)) = IK
      ENDDO
!--   

!.....Write table of parameter sizes (vjp 11/28/99)

      WRITE(16,4010) s%MNPROC,s%MNE,s%MNP,s%MNEI,s%MNOPE,s%MNETA,s%MNBOU,s%MNVEL,&
     s%MNTIF,s%MNBFR,s%MNSTAE,s%MNSTAV,s%MNSTAC,s%MNSTAM,global_here%NWLAT,global_here%NWLON,s%MNHARF,s%MNFFR
      IF (global_here%NWS.EQ.0) WRITE(16,4011)
      IF (global_here%NWS.EQ.1) WRITE(16,4012)
      IF (ABS(global_here%NWS).EQ.2) WRITE(16,4013)
      IF (global_here%NWS.EQ.3) WRITE(16,4014)
      IF (ABS(global_here%NWS).EQ.4) WRITE(16,4015)
      IF (ABS(global_here%NWS).EQ.5) WRITE(16,4115)
      IF (global_here%NWS.EQ.10) WRITE(16,4016)
      IF (global_here%NWS.EQ.11) WRITE(16,4017)
#ifdef HARM
      IF ((NFREQ.EQ.0).OR.(global_here%FMV.EQ.0.)) WRITE(16,4021)
      IF ((NFREQ.GE.1).AND.(global_here%FMV.NE.0.)) WRITE(16,4022)
#else
      WRITE(16,4021)
#endif      
      IF (global_here%ILUMP.EQ.0) WRITE(16,4031)
      IF (global_here%ILUMP.EQ.1) WRITE(16,4032)
      IF (global_here%IM.EQ.0) WRITE(16,4101)
      IF (global_here%IM.EQ.10) WRITE(16,4109)
      IF (global_here%IM.EQ.1) WRITE(16,4102)
      IF (global_here%IM.EQ.2) WRITE(16,4103)
      WRITE(16,4105)
      WRITE(16,4108)



 4010 FORMAT(' *****************************************************',/,&
     ' *   Based on information extracted from the ADCIRC  *',/,&
     ' *   UNIT 14 and 15 (grid and horiz run info) files  *',/,&
     ' *   the following paramter values will be set:      *',/,&
     ' *                                                   *',/,&
     ' *       MNPROC = ',I5,'                              *',/,&
     ' *       MNE = ',I8,1X,'     MNP  = ',I8,1X,'        *',/,&
     ' *       MNEI = ',I7,2X,'                            *',/,&
     ' *       MNOPE = ',I6,3X,'   MNETA = ',I6,3X,'       *',/,&
     ' *       MNBOU = ',I6,3X,'   MNVEL = ',I6,3X,'       *',/,&
     ' *       MNTIF = ',I6,3X,'   MNBFR = ',I6,3X,'       *',/,&
     ' *       MNSTAE = ',I5,4X,'  MNSTAV = ',I5,4X,'      *',/,&
     ' *       MNSTAC = ',I5,4X,'  MNSTAM = ',I5,4X,'      *',/,&
     ' *       MNWLAT = ',I5,4X,'  MNWLON = ',I5,4X,'      *',/,&
     ' *       MNHARF = ',I5,4X,'  MNFFR = ',I6,3X,'       *',/,&
     ' *                                                   *')
 4011 FORMAT(' *   Also, NO wind forcing will be used,             *')
 4012 FORMAT(' *   Also, Standard wind stress and pres will be used,*')
 4013 FORMAT(' *   Also, Semi-standard wind forcing will be used,  *')
 4014 FORMAT(' *   Also, Fleet numeric wind forcing will be used,  *')
 4015 FORMAT(' *   Also, PBL/JAG wind forcing will be used,        *')
 4115 FORMAT(' *   Also, Standard wind global_here%vel and pres will be used,  *')
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

!.....Close fort.15, fort.14, and fort.dg files

      CLOSE(14)
      CLOSE(15)

      RETURN 
      END SUBROUTINE
