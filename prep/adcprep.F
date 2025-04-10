C----------------------------------------------------------------------------
C
C                        ADCPREP Version 10_sb4_1 (  9/15/2005 )
C
C                     current for hp_DG_ADCIRC v10_sb4   9/15/2005
C----------------------------------------------------------------------------
C
C  Program Development History
C  ---------------------------
C   Written for ADCIRC_v24.05    ( S. Chippada 1996 )
C   Updated for ADCIRC_v24.05    ( M. Martinez 1997 )
C   Hilbert Space Filling Curve  ( C. Edwards 1997 )
C   Updated for ADCIRC_v34.04    ( V. Parr 1998 )
C   METIS 4.0 Interface          ( V. Parr 1999 )
C   Added UNIT 12 for polders    ( V. Parr 1999 )
C   hotstarts for Windstress     ( V. Parr 1999 )
C   Updated for ADCIRC_v34.08    ( V. Parr 1999 )
C   Modified for PADCIRC90       ( V. Parr & J. Romo 1999 )
C   Parallel version for 2 PEs   ( V. Parr & J. Romo 2000 )
C   Add hotprep1 command         ( V. Parr & J. Romo 2001 )
C   Add hotprep2 command         ( R. Luettich, 10/2001 )
C   Expanded to include 3D       ( T. Campbell 2002)
C   Bug fix                      ( J. Fleming/R. Luettich 2/2002)
C   Added wave stress capability ( M. Brown 3/2003)
C   NWS=3,6 fixes                ( M. Brown 3/2003)
C   NWS=+-4 revised to not convert to NWS = 5  ( R. Luettich 3/2003)
C   Wave stress routines revise  ( R. Luettich 3/2003)
C   Islands kept as islands      ( R. Luettich 5/2004)
C   Vic Parr's bug fix for NWS=+-4, Wave stress  ( R. Luettich 5/2004)
C   Updated for hp_DG_ADCIRC     ( S. Bunya 7/2005)
C   Updated for hp_DG_ADCIRC     ( S. Bunya 9/2005)
C
C----------------------------------------------------------------------------
C
C  ADCPREP performs 3 operations:
C
C    prep     =  localize grid and input files for a cold start
C    hotprep1 =  assumes all directories and files already in place from previous run segment
C                uses fort.15 file in master directory and creates new localized fort.15 files
C                for a hot start in the local processor directories.  In doing so it 
C                converts NWS = 2,4,5 to NWS = -2,-4,-5.  This forces 1st entry in met file to
C                correspond to the hot start time.
C    hotprep2 =  localize grid and input files (including hotstart file) for a hot start
C                NOTE: this skips the harmonic analysis/means & variance and 3D part of the
C                      hotstart file if these are present.
C
C  Pre-processor
C  -------------
C  The Pre-processor "prep" performs data decomposition of the ADCIRC grid 
C  and localizes global input files to subdomains of the decomposition for 
C  parallel runs of ADCIRC.
C   
C  Subdirectories "PExxxx" of the working directory ( where "xxxx" is a 
C  4-digit integer representing the MPI_rank of a Processing Element ) 
C  are created and "localized" input files are placed in the appropriate 
C  subdirectory with the same filename as the global file in the working 
C  directory.  A unit 14 and 15 file, called respectively, setup.14 and
C  setup.15  are written to the working directory to allow setup to
C  prepare the source code for a parallel run.
C
C  The global input files currently localized by ADCPREP are:
C
C   fort.14    Global Mesh File
C   fort.15    Global Input File
C   fort.12    Global StartDry Input File
C   fort.19    Global Time-Variant Elevation Boundary Conditions
C   fort.20    Global Time-Variant Normal Flow Boundary Conditions
C   fort.21    Global Nodal Friction and Eddy Viscosity
C   fort.22    Global Wind Stress and Atmospheric Pressure
C   fort.23    Global Wave Stress Forcing
C   fort.67    Global Hot Start file
C   fort.68    Global Hot Start file
C   fort.150   METIS input information
C   fort.151   METIS output information
C 
C  When a parallel ADCIRC job is run, each PExxxx takes its input from its 
C  subdirectory PExxxx and writes local output files to the same directory.
C
C----------------------------------------------------------------------------
C
      PROGRAM ADCPREP
      USE PRE_GLOBAL 
C
      INTEGER I,PE,NWTYPE
      REAL TOL
      LOGICAL FOUND
      CHARACTER DIRCMD*72,CMD*6,PENUM*6,PHASE*7,ISET*6,DIR1*80,DIR2*80,
     &    CONVDIR*6, JOB*80
C
#ifdef CMPI
#ifdef CRAY
      INTEGER*8 NPROC_PRE,MYRANK_PRE
#else
      INTEGER NPROC_PRE,MYRANK_PRE
#endif
      DATA NPROC_PRE/1/,MYRANK_PRE/0/
      INTEGER MNPROC2,MNP2,MNE2,MNEI2,NWEIR2
      INTEGER,ALLOCATABLE :: NNEG2(:,:), WEIR2(:), WEIRD2(:)
#endif

C
C-- Initialize constants for case ICS = 2 inputs
C
      R  =  6378206.4D0
      DEG2RAD = 3.141592653589793D0/180.0D0
      RAD2DEG = 180.0D0/3.141592653589793D0
      
C
C     Initialize flags.
C
      USE_DEFAULT = .FALSE.
      APERIODIC_FLOW_BC = .FALSE.
      PARTITION  = .FALSE.
      PREP_ALL   = .FALSE.
      PREP_13    = .FALSE.
      PREP_15    = .FALSE.
      HOT_LOCAL   = .FALSE.
      HOT_GLOBAL  = .FALSE.

#ifdef CMPI
C
C-- If running parallel start up message-passing interface
C   Note: parallel execution is restricted to 2 processors
C
      CALL MESSAGE_INIT(NPROC_PRE,MYRANK_PRE)
      IF (NPROC_PRE.NE.2) THEN
        PRINT*,"ERROR: Parallel execution of ADCIRC Preprocessor"
        PRINT*,"ERROR: must be on 2 processors.  Now quiting..."
        CALL SYNCHRONIZE()
        CALL MESSAGE_FINI()
        STOP
      ENDIF
      IF (MYRANK_PRE.EQ.1) GO TO  99
#endif

C
C-- Say Hello Gracie
C
      print *," *****************************************"
      print *," ADCPREP Fortran90  Version 1.8  5/18/2004"
      print *," Parallel version of ADCIRC Pre-processor "
      print *," *****************************************"
      print *, " "
C
C-- Prompt for user input
C
      print *, "Input number of processors for parallel ADCIRC run:"
      READ(*,*) MNPROC
C
C-- Copy MNPROC to NPROC since they are to be the same
C
      NPROC=MNPROC

      print *, "prepare for coldstart/hotstart type 1/hotstart type 2 ",
     &         "( prep / hotprep1 / hotprep2 )"
      READ(*,*)  JOB
C
C-- Short circuit rest of code to simply relocalize a fort.15 file for
C     a hotstart assuming all other files are already in place.
C
      IF(TRIM(JOB).EQ."hotprep1") THEN
        CALL HOTPREP()
        GO TO 999 
      ENDIF

C
C---------------------------------------------------------------------------
C  Complete Parallel ADCIRC Pre-Processor Starts here
C---------------------------------------------------------------------------
C
      CALL MAKE_DIRNAME()
      CALL SIZEUP()
      CALL ALLOC_MAIN1()
C
C-- Read the Global Grid File ( Unit 14 )
C
        CALL READ14() 
        print *, 'Global Grid file read successfully'


#ifdef  CMPI
        WRITE(150,*) MNPROC
        WRITE(150,*) MNP
        WRITE(150,*) MNE
        WRITE(150,*) MNEI
        WRITE(150,*) NWEIR
        DO I=1, MNE
           WRITE(150,'(3I8)') NNEG(1,I),NNEG(2,I),NNEG(3,I)
        ENDDO
        DO I=1, NWEIR
           WRITE(150,'(2I8)') WEIR(I), WEIRD(I)
        ENDDO
        CLOSE(150) 
#endif
C
C-- Read the Global Input File ( Unit 15 )
C
        CALL READ15()
        print *, 'Global Input file read successfully'

C
C-- Read Hard Bottom File ( fort.hb )
C
        IF(HBFLAG.EQ.1) THEN
          CALL READHB() 
          print *, 'Hard Bottom file read successfully'
        ENDIF

C
C-- Partiton Nodes with METIS 4.0 Graph Partition Library
C

#ifdef  CMPI
 99     CONTINUE 
        CALL SYNCHRONIZE()
        IF (MYRANK_PRE.EQ.1) THEN
          READ(150,*) MNPROC2
          READ(150,*) MNP2
          READ(150,*) MNE2
          READ(150,*) MNEI2
          READ(150,*) NWEIR2
          ALLOCATE(NNEG2(3,MNE2)) 
          DO I=1, MNE2
             READ(150,'(3I8)') NNEG2(1,I),NNEG2(2,I),NNEG2(3,I)
          ENDDO
          ALLOCATE(WEIR2(NWEIR2))
          ALLOCATE(WEIRD2(NWEIR2))
          DO I=1, NWEIR2
             READ(150,'(2I8)') WEIR2(I), WEIRD2(I)
          ENDDO
          CLOSE(150) 
          CALL METIS2(MNPROC2,MNP2,MNE2,MNEI2,NWEIR2,
     &                NNEG2,WEIR2,WEIRD2)
        ENDIF
        CALL SYNCHRONIZE()
        IF (MYRANK_PRE.EQ.0) THEN
          DO I=1, MNP
            READ(151,*) PROC(I)
          ENDDO
        ELSE
          GO TO 999
        ENDIF
#else
        CALL METIS()
#endif

        print *, " "
        print *, 'METIS has partitioned nodes'
C
C-- Decompose the ADCIRC grid into MNPROC subdomains
C
        print *, " "
        print *, "Determine the parameters MNPP and MNEP"
        CALL DOMSIZE()
C
        print *, "Allocate arrays dimensioned by MNPP and MNEP"
        CALL ALLOC_MAIN2()

        print *, " "
        print *, "Decomposition of grid begins"
        CALL DECOMP()
        print *, "Decomposition successful"
C
C-- Create MNPROC sub-directories of the working directory
C
        DO PE=0, MNPROC-1
           CMD = 'mkdir '
           PENUM  = 'PE0000'
           CALL IWRITE(PENUM,3,6,PE)
           DIRCMD= CMD//PENUM  
           CALL ISHELL(DIRCMD)
        ENDDO
C
C-- Write Message-Passing File for each PE
C
        print *, "Writing Message-Passing Info Files for each PE"
        CALL PREP18()
Csb-DG1
C
C-- Write Message-Passing File for each PE
C
        print *, "Writing DG Message-Passing Info Files for each PE"
        CALL PREPDG18()
C--

C
C-- Write a Local Grid File ( fort.14 ) for each PE
C
        print *, "Writing Local UNIT 14 (Grid) File and"
        print *, "fort.hb (Hard Bottom) File for each PE"
        CALL PREP14_AND_HB()
C
C-- Write a Local Input file ( fort.15 ) for each PE
C
        print *, "Writing Local UNIT 15 (Run Info) File for each PE"
        CALL PREP15()
C
Casey 101118: Need to copy SWAN control file into local sub-directories.
#ifdef SWAN
        CALL PREPUNSWAN()
#endif
C
C-- Write a Local StartDry Input file ( fort.12 ) for each PE

        if (NOLIFA.EQ.3) then
          print *, "Writing Local UNIT 12 (StartDry) File for each PE"
          CALL PREP12()
        endif

        CALL PREP13()
C
C-- If required write a Local fort.19 file for each PE
C
        IF ((NOPE.GT.0).AND.(NBFR.EQ.0)) THEN
          print *, "Ready to Write Local UNIT 19 File for each PE"
          CALL PREP19()
        ENDIF
C

C-- Namo 7/12/21: Write a local fort.20 file for each PE
        print *, "Number of flow nodes is: ", NFLBN
        if ((NFLBN .ne. 0) .and. (APERIODIC_FLOW_BC)) then
          print *, "Ready to write Local UNIT 20 File for each PE"
          call PREP20()
        endif

C-- If required write a Local fort.21 file for each PE
C
        IF (NWP.EQ.1)  THEN
          print *, "Ready to Write Local UNIT 21 File for each PE"
          CALL PREP21()
        ENDIF
C
C-- If required write a Local Wind Stress file for each PE
C
        IF (NWS.NE.0)  THEN
          print *, "Ready to Write Local UNIT 22 File for each PE"
          IF((NWS.EQ.1).OR.(ABS(NWS).EQ.2).OR.(ABS(NWS).EQ.4).OR.    
     &       (NWS.EQ.3).OR.(ABS(NWS).EQ.5).OR.(    NWS .EQ.6) ) THEN
            CALL PREP22()
c         ELSEIF (NWS.EQ.10.OR.NWS.EQ.11) THEN
c           CALL PREP200()
          ENDIF
        ENDIF
C
C-- If required write a Local Wave Stress file for each PE
C
        IF(NRS.EQ.1) CALL PREP23()
C
C-- If required write a Local Hot Start file for each PE
C
        IF((IHOT.EQ.67).OR.(IHOT.EQ.68)) THEN
          print *, "Ready to write Local Hot Start Files for each PE"
          CALL PREP67_68()
          ENDIF
C
C-- Save domain-decomposition information for post-processor
C
        print *, ""
        print *, "Writing domain-decomposition file for post-processor"
        CALL PREP80()
Csb--
C-- Save domain-decomposition information for post-processor
C
        print *, ""
        print *, "Writing fort.dg file for each PE"
        CALL PREP_DG()
C
        if (MIC.EQ.1) then
          CALL READIC()
          print *, 'Initial conditions files read successfully'
          print *, ""
          print *, "Writing initial conditions files for each PE"
          CALL PREP_IC()
        endif
C--
C
        print *, ""
        print *, "Finished pre-processing input files"
C
 999    CONTINUE
#ifdef CMPI
        CALL SYNCHRONIZE()
        CALL MESSAGE_FINI()    
#endif
C
        STOP
      END





      SUBROUTINE GETMSG( STRING, MSG )
      INTEGER I, I1
      CHARACTER  STRING*(*),MSG*(*), TARGET
c
      I1 = 0
      TARGET = "!"
c
c-- Find beginning of message
c
      DO I=1, 80
         IF (STRING(I:I).EQ.TARGET) THEN
           I1 = I
           GOTO 100
         ENDIF
      ENDDO
c
 100  CONTINUE
c
c--Copy message to ouput string
c
      DO I=1, I1-1
         MSG(I:I) = " "
      ENDDO
      MSG(I1:80)  = STRING(I1:80)
c
      RETURN
      END



      SUBROUTINE NEWINDEX ( ISTRING, OSTRING, INDEX )
      INTEGER I,I1,I2,I3,I4,INDEX
      CHARACTER  ISTRING*(*),OSTRING*(*),TARGET
      CHARACTER TEMP1*80, TEMP2*100
c
      I1 = 0
      I2 = 0
      I3 = 0
      I4 = 0
      TARGET = " "
c
c-- Find first non-blank character of String
c
      DO I=1, 80
         IF (ISTRING(I:I).NE.TARGET) THEN
           I1 = I
           GOTO 100
         ENDIF
      ENDDO
c
c-- Find next blank character of String
c
 100  CONTINUE
      DO I=I1+1,80
         IF (ISTRING(I:I).EQ.TARGET) THEN
           I2 = I
           GOTO 200
         ENDIF
      ENDDO
c
 200  CONTINUE
c
c-- Create a temporary string containing new index
c
      WRITE(TEMP1(1:80),'(I8)') INDEX
c
c-- Find first non-blank character of String
c
      DO I=1, 80
         IF (TEMP1(I:I).NE.TARGET) THEN
           I3 = I
           GOTO 300
         ENDIF
      ENDDO
c
c-- Find next blank character of String
c
 300  CONTINUE
      DO I=I3+1,80
         IF (TEMP1(I:I).EQ.TARGET) THEN
           I4 = I
           GOTO 400
         ENDIF
      ENDDO
c
 400  CONTINUE
c
Cdebug print *, "i1 i2 i3 i4 ",I1, I2, I3 , I4
      TEMP2(1:100) = TEMP1(I3-1:I4-1)//ISTRING(I2:80)
c
c-- Write out first 80 characters of concatenated strings
c
      OSTRING(1:80) = TEMP2(1:80)
c
      RETURN
      END


      SUBROUTINE INSERT( ISTRING, OSTRING, NUMS, N )
      INTEGER I,J,I1,N,NUMS(N)
      CHARACTER  ISTRING*80,OSTRING*80,BLANK
      CHARACTER  TEMP1*80, TEMP2*160
c
      I1 = 0
      BLANK = " "
c
c-- Create Tempoarary String TEMP1 containing NUMS
c
      IF (N.EQ.1) THEN
        WRITE(TEMP1(1:80),'(I8)') NUMS(1)
Casey 101118: Changed this section to allow N = 3.
      ELSEIF (N.EQ.2) THEN
        WRITE(TEMP1(1:80),'(2I8)') NUMS(1),NUMS(2)
      ELSE
        WRITE(TEMP1(1:80),'(3I8)') NUMS(1),NUMS(2),NUMS(3)
      ENDIF
c
c-- Find length of TEMP1 string
c
      DO I=80,1,-1
         IF (TEMP1(I:I).NE.BLANK) THEN
           LEN1 = I
           GOTO 10
         ENDIF
      ENDDO
 10   CONTINUE  
c
c-- Scan input string for character after old number list
c
      I = 1
      DO NUM=1, N
         DO J=I,80
            IF (ISTRING(J:J).NE.BLANK) THEN
              I = J
              GOTO 100
            ENDIF
         ENDDO
 100     CONTINUE
         DO J=I,80
           IF (ISTRING(J:J).EQ.BLANK) THEN
             I = J
             GOTO 200
           ENDIF
         ENDDO
 200     CONTINUE
      ENDDO
      I1 = MAX(0,I)
c
c-- Insert Integer List into Message
c
      IF (I1.NE.0) THEN
C-- if there is a message
        TEMP2(1:160) = TEMP1(1:LEN1+1)//ISTRING(I1:80)
      ELSE
        TEMP2(1:160) = TEMP1(1:LEN1+1)
      ENDIF
c
c-- Write out first 80 characters of concatenated string
c
      OSTRING(1:80) = TEMP2(1:80)
c
      RETURN
      END



      SUBROUTINE HOTPREP()
      USE PRE_GLOBAL
C
C---------------------------------------------------------------------------C
C                                                                           C
C  This routine modifies the local fort.15 files for running a hotstart job.C
C                                                                           C
C  The Decomposition Variables are defined in the include file pre_global.F C
C                                                                           C
C---------------------------------------------------------------------------C
C
        IMPLICIT NONE
        INTEGER I
        INTEGER IHOTS,NWSS,IPROC
        CHARACTER LOCFN*14,CMD*21,PE*6
        CHARACTER*80 BUFIN
C
        print *, "Enter restart file you wish to use: ( 67 / 68 )"
        READ(*,*) IHOTS
        OPEN(15,FILE="fort.15")
        DO I=1, 16
           READ(15,*) BUFIN
        ENDDO
        READ(BUFIN,*)  NWSS
        IF (NWSS.EQ.2.OR.NWSS.EQ.4.OR.NWSS.EQ.5) NWSS = -NWSS
        CLOSE(15)
C
C-- Modify each local  ( fort.15 ) file
C
        DO IPROC = 1,MNPROC
C
           LOCFN(1:14) = 'PE0000/fort.15'
           CALL IWRITE(LOCFN,3,6,IPROC-1)
           OPEN (15+IPROC,FILE=LOCFN)
           OPEN (10,FILE='tmp')
C
           DO I=1, 5
              READ(15+IPROC,'(A80)') BUFIN
              WRITE(10,'(A80)') BUFIN
           ENDDO
           READ(15+IPROC,'(A80)') BUFIN
           WRITE(10,*) IHOTS
C
           DO I=1, 9
              READ(15+IPROC,'(A80)') BUFIN
              WRITE(10,'(A80)') BUFIN
           ENDDO
           READ(15+IPROC,'(A80)') BUFIN
           WRITE(10,*) NWSS
C
199        READ(15+IPROC,'(A80)',END=200) BUFIN
           WRITE(10,'(A80)') BUFIN
           GO TO 199
200        CONTINUE
           CLOSE(15+IPROC)
           CLOSE(10)
C
           CMD(1:21) = "mv tmp "//LOCFN(1:14)
           CALL ISHELL(CMD)    
        ENDDO
C
        RETURN
        END
