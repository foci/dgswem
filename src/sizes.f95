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
     
  end type sizes_type
  
END MODULE SIZES
