!**************************************************************************
!    last changes in this file VERSION 41.11
!
!  mod history
!  v41.11 - 09/14/01 - RL - eliminated MNWLAT, MNWLON
!  v41.06 - 04/02/01 - RL - eliminated MNWP
!  v41.06mxxx - date - programmer - describe change 
!                    - mark change in code with  cinitials-mxxx
!-----------------------------------------------------------------------
!
!     
!     v10_sb5    - 10/11/05 - sb - Table that associates nodes to elements (NDEL)
!
!**************************************************************************
!
      MODULE SIZES
      IMPLICIT NONE
!
!...SET NUMBER OF BYTES "SZ" IN REAL(SZ) DECLARATIONS       
!...SET "NBYTE" FOR PROCESSING INPUT DATA RECORD LENGTH

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

!...SET MAX OF DIGITS OF PRECISION "NPREC" THE GRID CAN BE EXPECTED TO HAVE
!...NOTE: IF THE GRID WAS BUILT ON A 32 BIT COMPUTER, IT SHOULD BE
!   ACCURATE TO ABOUT 7 DIGITS.  THUS, IF THE NODAL SPACING REQUIRES MORE
!   THAN 5 DIGITS OF PRECISION, THE MODEL RESULTS MAY NOT BE TRUSTWORTHY.
 
      INTEGER, PARAMETER ::  NPREC=7
!
      INTEGER ::  MNPROC,MNE,MNP,MNEI,MNOPE,MNETA,MNBOU,MNVEL,&
       MNTIF,MNBFR,MNFFR,MNSTAE,MNSTAV,MNSTAC,MNSTAM,MNHARF
      INTEGER, TARGET :: layers
!sb-
      INTEGER MNNDEL
!--

      LOGICAL C2DDI,C3D,C3DDSS,C3DVS,CLUMP,CTIP,CHARMV
!
! For Definition of Working Directory
!
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
      
!---------------------end of data declarations--------------------------------C


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
!asey 121128: The DG output will be written to local files anyway.
!             This setting only affects the SWAN files.
!     WRITE_LOCAL_FILES = 1
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

!.....Default root working directory

      IF (INPUTDIR == "") THEN
        ROOTDIR = '.'
        INPUTDIR = '.'
      ELSE
        ROOTDIR = INPUTDIR
      ENDIF
      
!.....Set the global input directory

      GBLINPUTDIR = ROOTDIR
      
#ifdef CMPI
!      iprefix = 0
!      do i = 1, 2
!        INQUIRE(file=TRIM(ROOTDIR)//TRIM(PREFIX(i))//'/'//'fort.14',
!     &    exist=fileFound)
!       if (fileFound) then
!          iprefix = i
!          exit
!        end if
!      end do
!      if (.not. fileFound) then
!        print *, "Failed to find prefix directory"
!        call msg_abort()
!      end if

!      WRITE(INPUTDIR,'(2A)') TRIM(ROOTDIR),PREFIX(iprefix)
!      LNAME = LEN_TRIM(INPUTDIR)
!      WRITE(INPUTDIR(LNAME-3:LNAME),'(I4.4)') MYPROC
#else
!      WRITE(INPUTDIR,'(A)') TRIM(ROOTDIR)
#endif

!asey 121128: Uncommented these lines.
       if (GLOBALDIR == "") then
         ROOTDIR = '.'
       else
         ROOTDIR = GLOBALDIR
       endif

       WRITE(GLOBALDIR,'(A)') TRIM(ROOTDIR)


#ifdef CMPI
!      WRITE(LOCALDIR,'(2A)') TRIM(ROOTDIR),TRIM(PREFIX(iprefix))
!      LNAME = LEN_TRIM(LOCALDIR)
!      WRITE(LOCALDIR(LNAME-3:LNAME),'(I4.4)') MYPROC
!      call MAKEDIR(trim(LOCALDIR))
#else
!      WRITE(LOCALDIR,'(A)') TRIM(ROOTDIR)
#endif
!      if (WRITE_LOCAL_FILES) GLOBALDIR = LOCALDIR
!      HOTSTARTDIR = LOCALDIR

      RETURN
      END SUBROUTINE

      END MODULE SIZES
