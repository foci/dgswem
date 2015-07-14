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
     ! 6 <- screen unit?

     integer :: fort12unit
     integer :: fort13unit
     integer :: fort14unit
     integer :: fort15unit
     integer :: fort16unit ! output
     integer :: fort17unit ! output
     integer :: fort19unit ! input cstart

     integer :: fort20unit ! input cstart
     integer :: fort22unit ! input cstart
     integer :: fort23unit
     integer :: fort24unit

     integer :: fort4lunit
!     integer :: fort41unit !4441
     ! 42 slopelimiter?

     integer :: fort55unit ! output

     integer :: fort61unit ! output
     integer :: fort62unit ! output
     integer :: fort63unit ! output
     integer :: fort64unit ! output

     integer :: fort67unit
     integer :: fort68unit


     integer :: fort71unit ! output
     integer :: fort72unit ! output
     integer :: fort73unit ! output
     integer :: fort74unit ! output

     integer :: fort80unit ! input

     integer :: fort81unit ! output
     integer :: fort82unit ! output
     integer :: fort83unit ! output
     integer :: fort84unit ! output
     integer :: fort85unit ! output
     integer :: fort88unit ! output
     integer :: fort89unit ! output

     integer :: fort94unit ! output
     integer :: fort96unit ! output

     integer :: rads64unit ! 164 or 716?  ! output
     
     integer :: fort199unit
     
     integer :: fortdgunit 

     integer :: dg18unit ! 18??

     ! 263, 264, 214, 288, 289, 290, 291 ! Hot starts, opened in write_results

     ! 333 wind.f95

     ! 440 rhs_dg_hydro
     ! 444
     ! 445

     ! 632, 642, 652 opened in write_results, DG.**.IC

     ! 667 opened in wind.f95, fort.22.adc

     integer :: tecfileunit ! 777
     integer :: tecfile_maxunit ! 778

     ! 963 ! maxele.63 opened in write_results
     
     ! 4441 ?

  end type sizes_type
  
END MODULE SIZES
