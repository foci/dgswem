C**************************************************************************
C    last changes in this file VERSION 41.11
C
C  mod history
C  v41.11 - 09/14/01 - RL - eliminated MNWLAT, MNWLON
C  v41.06 - 04/02/01 - RL - eliminated MNWP
C  v41.06mxxx - date - programmer - describe change 
C                    - mark change in code with  cinitials-mxxx
C-----------------------------------------------------------------------
C
C     
C     v10_sb5    - 10/11/05 - sb - Table that associates nodes to elements (NDEL)
C
C**************************************************************************
C
      MODULE SIZES
      IMPLICIT NONE
C
C...SET NUMBER OF BYTES "SZ" IN REAL(SZ) DECLARATIONS       
C...SET "NBYTE" FOR PROCESSING INPUT DATA RECORD LENGTH

#ifdef  REAL4
      INTEGER, PARAMETER :: SZ = 4
      INTEGER, PARAMETER :: NBYTE=4
#endif
#ifdef REAL16
      INTEGER, PARAMETER :: SZ = 16
      INTEGER, PARAMETER :: NBYTE=16
#endif
#ifndef REAL4
#ifndef REAL16
      INTEGER, PARAMETER :: SZ = 8
      INTEGER, PARAMETER :: NBYTE=8
#endif
#endif

C...SET MAX OF DIGITS OF PRECISION "NPREC" THE GRID CAN BE EXPECTED TO HAVE
C...NOTE: IF THE GRID WAS BUILT ON A 32 BIT COMPUTER, IT SHOULD BE
C   ACCURATE TO ABOUT 7 DIGITS.  THUS, IF THE NODAL SPACING REQUIRES MORE
C   THAN 5 DIGITS OF PRECISION, THE MODEL RESULTS MAY NOT BE TRUSTWORTHY.
 
      INTEGER, PARAMETER ::  NPREC=7
C
      INTEGER ::  MNPROC,MNE,MNP,MNEI,MNOPE,MNETA,MNBOU,MNVEL,
     &  MNTIF,MNBFR,MNFFR,MNSTAE,MNSTAV,MNSTAC,MNSTAM,MNHARF,layers
Csb-
      INTEGER MNNDEL
C--

      LOGICAL C2DDI,C3D,C3DDSS,C3DVS,CLUMP,CTIP,CHARMV
C
C For Definition of Working Directory
C
      INTEGER,SAVE :: MYPROC

#ifdef CMPI
      INTEGER,SAVE :: LNAME = 6
      CHARACTER*6,SAVE :: DIRNAME = 'PE0000'
#else
      INTEGER,SAVE :: LNAME = 1
      CHARACTER*1,SAVE :: DIRNAME = '.'
#endif

      LOGICAL      :: WRITE_LOCAL_FILES
      LOGICAL      :: WRITE_LOCAL_HOT_START_FILES
      CHARACTER(256),  SAVE :: ROOTDIR
      CHARACTER(2048), SAVE :: INPUTDIR, GLOBALDIR, LOCALDIR
      CHARACTER(2048), SAVE :: GBLINPUTDIR, HOTSTARTDIR
      
C---------------------end of data declarations--------------------------------C


      CONTAINS

      SUBROUTINE MAKE_DIRNAME()
      INTEGER :: LNAME, IARGC, ARGCOUNT, I, iprefix
      CHARACTER(2048) :: CMDLINEARG
      CHARACTER(8)    :: PREFIX(2) = (/ '/PE0000 ', '/DOM0000' /)
      logical         :: fileFound

      INPUTDIR  = ""
      GLOBALDIR = ""
      LOCALDIR  = ""
      ARGCOUNT  = IARGC()
      
!      WRITE_LOCAL_FILES = MNPROC == 1
Casey 121128: The DG output will be written to local files anyway.
C             This setting only affects the SWAN files.
C     WRITE_LOCAL_FILES = 1
      WRITE_LOCAL_FILES = .FALSE.
      WRITE_LOCAL_HOT_START_FILES = .TRUE.
      
#ifdef CMPI
      WRITE(DIRNAME(3:6),'(I4.4)') MYPROC
#else
      MYPROC=0
#endif

      IF (ARGCOUNT > 0) THEN
        I = 0
        DO WHILE (I < ARGCOUNT)
          I = I + 1
          CALL GETARG( I, CMDLINEARG )
          IF (CMDLINEARG(1:2) == "-I") THEN
            I = I + 1
            CALL GETARG( I, INPUTDIR )
          ELSEIF (CMDLINEARG(1:2) == "-O") THEN
            I = I + 1
            CALL GETARG(I,GLOBALDIR)
          ELSEIF (CMDLINEARG(1:2) == "-L") THEN
            WRITE_LOCAL_FILES = .TRUE.
          ENDIF
        ENDDO
      ENDIF

C.....Default root working directory

      IF (INPUTDIR == "") THEN
        ROOTDIR = '.'
        INPUTDIR = '.'
      ELSE
        ROOTDIR = INPUTDIR
      ENDIF
      
C.....Set the global input directory

      GBLINPUTDIR = ROOTDIR
      
#ifdef CMPI
C      iprefix = 0
C      do i = 1, 2
C        INQUIRE(file=TRIM(ROOTDIR)//TRIM(PREFIX(i))//'/'//'fort.14',
C     &    exist=fileFound)
C       if (fileFound) then
C          iprefix = i
C          exit
C        end if
C      end do
C      if (.not. fileFound) then
C        print *, "Failed to find prefix directory"
C        call msg_abort()
C      end if

C      WRITE(INPUTDIR,'(2A)') TRIM(ROOTDIR),PREFIX(iprefix)
C      LNAME = LEN_TRIM(INPUTDIR)
C      WRITE(INPUTDIR(LNAME-3:LNAME),'(I4.4)') MYPROC
#else
C      WRITE(INPUTDIR,'(A)') TRIM(ROOTDIR)
#endif

Casey 121128: Uncommented these lines.
       if (GLOBALDIR == "") then
         ROOTDIR = '.'
       else
         ROOTDIR = GLOBALDIR
       endif

       WRITE(GLOBALDIR,'(A)') TRIM(ROOTDIR)


#ifdef CMPI
C      WRITE(LOCALDIR,'(2A)') TRIM(ROOTDIR),TRIM(PREFIX(iprefix))
C      LNAME = LEN_TRIM(LOCALDIR)
C      WRITE(LOCALDIR(LNAME-3:LNAME),'(I4.4)') MYPROC
C      call MAKEDIR(trim(LOCALDIR))
#else
C      WRITE(LOCALDIR,'(A)') TRIM(ROOTDIR)
#endif
C      if (WRITE_LOCAL_FILES) GLOBALDIR = LOCALDIR
C      HOTSTARTDIR = LOCALDIR

      RETURN
      END SUBROUTINE

      END MODULE SIZES
