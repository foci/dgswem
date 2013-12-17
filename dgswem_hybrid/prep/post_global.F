C----------------------------------------------------------------------
C
C  updated by S.Bunya Aug 2005 to read elemental decomposition data         
C
C----------------------------------------------------------------------
      MODULE POST_GLOBAL
      IMPLICIT NONE
      INTEGER MNPROC,MNP,MNE,MNPP,MNSTAE,MNSTAV,MNHARF,MNWLAT,MNWLON
Csb
      INTEGER MNEP
C--      
C
C   PARM14   section
      INTEGER NELG,NNODG
C
C   STRING14 section
      CHARACTER*80 AGRID
C
C   PARM15  section
      INTEGER NOUTE,NSPOOLE,NSTAE
      INTEGER NOUTV,NSPOOLV,NSTAV
      INTEGER NOUTC,NSPOOLC,NSTAC
      INTEGER NOUTGE,NSPOOLGE
      INTEGER NOUTGV,NSPOOLGV
      INTEGER NOUTGC,NSPOOLGC
      INTEGER NOUTGW,NSPOOLGW
      INTEGER NHASE,NHASV,NHAGE,NHAGV
C
C
C  STRING15 section
      CHARACTER*80 RUNDES,RUNID
C
C  PARM15-3DVS section
      INTEGER :: I3DSD,NSPO3DSD,NHN3DSD
      INTEGER :: I3DSV,NSPO3DSV,NHN3DSV
      INTEGER :: I3DST,NSPO3DST,NHN3DST
      REAL(4) :: TO3DSSD,TO3DFSD
      REAL(4) :: TO3DSSV,TO3DFSV
      REAL(4) :: TO3DSST,TO3DFST
      INTEGER :: I3DGD,NSPO3DGD
      INTEGER :: I3DGV,NSPO3DGV
      INTEGER :: I3DGT,NSPO3DGT
      REAL(4) :: TO3DSGD,TO3DFGD
      REAL(4) :: TO3DSGV,TO3DFGV
      REAL(4) :: TO3DSGT,TO3DFGT
C
C  ELESTAT section
      REAL TOUTSE,TOUTFE,TOUTSGE,TOUTFGE
C
C  VELSTAT section
      REAL TOUTSGV,TOUTFGV,TOUTSV,TOUTFV
C
C  CONSTAT  section
      REAL TOUTSC,TOUTFC,TOUTSGC,TOUTFGC
      REAL TOUTSGW,TOUTFGW
C
C--Degress-to-Radians and Radians-to-Degrees
C
C  CONVERT section
      REAL*8 DEG2RAD,RAD2DEG,R
C
C
C--Used to DASD file operations
C
C DASD
      INTEGER NBYTE
C
C--------------------------------------------------------------------------C
C                                                                          C
C              DATA DECOMPOSITION DECLARATIONS BEGIN HERE                  C
C                                                                          C
C--------------------------------------------------------------------------C
C
C--Local Map Variable Declarations
C
Cvjp Note:  NSTACP  and IMAP_STAC_LG are not used yet
C
C LOCALI section
      INTEGER NPROC
      INTEGER, ALLOCATABLE ::  NNODP(:)
      INTEGER, ALLOCATABLE ::  NOD_RES_TOT(:)
      INTEGER, ALLOCATABLE ::  NSTAEP(:), NSTAVP(:), NSTACP(:)
Csb
      INTEGER, ALLOCATABLE ::  NELP(:)
C--
C
C
C--Local-to-Global Mapping Variables
C
C  LOC2G section
      INTEGER, ALLOCATABLE :: IMAP_NOD_LG(:,:)
      INTEGER, ALLOCATABLE :: IMAP_STAE_LG(:,:)
      INTEGER, ALLOCATABLE :: IMAP_STAV_LG(:,:)
      INTEGER, ALLOCATABLE :: IMAP_STAC_LG(:,:)
Csb
      INTEGER, ALLOCATABLE :: IMAP_EL_LG(:,:)
C--
C
C--Global-to-Local Mapping Variables
C
C  GLOB2L section
      INTEGER, ALLOCATABLE :: IMAP_NOD_GL(:,:)
C
C--3DVS section
C
      INTEGER,ALLOCATABLE  ::  N3DSDP(:),N3DSVP(:),N3DSTP(:)
      INTEGER,ALLOCATABLE  ::  IMAP_3DSD_LG(:,:)
      INTEGER,ALLOCATABLE  ::  IMAP_3DSV_LG(:,:)
      INTEGER,ALLOCATABLE  ::  IMAP_3DST_LG(:,:)
C
C--------------------------------------------------------------------------C
c   END OF DECLARATIONS
C--------------------------------------------------------------------------C

       CONTAINS


       SUBROUTINE ALLOC_MAIN1()
C                     
      ALLOCATE ( NNODP(MNPROC))
      ALLOCATE ( NOD_RES_TOT(MNPROC))
      ALLOCATE ( NSTAEP(MNPROC), NSTAVP(MNPROC), NSTACP(MNPROC))
      ALLOCATE ( IMAP_NOD_LG(MNPP,MNPROC))
      ALLOCATE ( IMAP_NOD_GL(2,MNP))

      ALLOCATE ( IMAP_STAE_LG(MNSTAE,MNPROC))
      ALLOCATE ( IMAP_STAV_LG(MNSTAV,MNPROC))
      ALLOCATE ( IMAP_STAC_LG(MNSTAV,MNPROC))
Csb
      ALLOCATE ( NELP(MNPROC))
      ALLOCATE ( IMAP_EL_LG(MNEP,MNPROC))
C--
C
      END SUBROUTINE ALLOC_MAIN1
      END MODULE POST_GLOBAL

