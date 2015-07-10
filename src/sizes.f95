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

  type sizes_type
     
     INTEGER ::  MNPROC,MNE,MNP,MNEI,MNOPE,MNETA,MNBOU,MNVEL,&
          MNTIF,MNBFR,MNFFR,MNSTAE,MNSTAV,MNSTAC,MNSTAM,MNHARF
     INTEGER :: layers
     !sb-
     INTEGER MNNDEL
     !--
     
     LOGICAL C2DDI,C3D,C3DDSS,C3DVS,CLUMP,CTIP,CHARMV
     !
     ! For Definition of Working Directory
     !
     INTEGER :: MYPROC
     
#ifdef CMPI
     INTEGER :: LNAME = 6
     CHARACTER*6 :: DIRNAME = 'PE0000'
#else
#ifdef HPX
     INTEGER :: LNAME = 6
     CHARACTER*6 :: DIRNAME = 'PE0000'
#else
     INTEGER :: LNAME = 1
     CHARACTER*1 :: DIRNAME = '.'
#endif
#endif
     
     LOGICAL      :: WRITE_LOCAL_FILES
     LOGICAL      :: WRITE_LOCAL_HOT_START_FILES
     CHARACTER(256) :: ROOTDIR
     CHARACTER(2048) :: INPUTDIR, GLOBALDIR, LOCALDIR
     CHARACTER(2048) :: GBLINPUTDIR, HOTSTARTDIR



     ! file unit numbers
     ! 
     integer :: fort12unit 
     integer :: fort13unit
     integer :: fort14unit
     integer :: fort15unit
     integer :: fort16unit
     integer :: fort17unit
     integer :: fort19unit

     integer :: fort20unit
     integer :: fort22unit
     integer :: fort23unit
     integer :: fort24unit

     integer :: fort41unit !4441

     integer :: fort55unit

     integer :: fort61unit
     integer :: fort62unit
     integer :: fort63unit
     integer :: fort64unit

     integer :: fort67unit
     integer :: fort68unit


     integer :: fort71unit
     integer :: fort72unit
     integer :: fort73unit
     integer :: fort74unit

     
     integer :: fort81unit
     integer :: fort82unit
     integer :: fort83unit
     integer :: fort84unit
     integer :: fort85unit
     integer :: fort88unit
     integer :: fort89unit

     integer :: fort94unit
     integer :: fort96unit

     integer :: rads64unit ! 164 or 716?
     
     integer :: fort199unit
     
     integer :: fortdgunit 

     integer :: dg18unit ! 18??

     integer :: tecfileunit ! 777
     integer :: tecfile_maxunit ! 778


  end type sizes_type
  
END MODULE SIZES
