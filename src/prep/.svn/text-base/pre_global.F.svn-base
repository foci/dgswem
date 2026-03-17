C----------------------------------------------------------------------------
C
C                           MODULE PRE_GLOBAL
C
C----------------------------------------------------------------------------
C
C                  For use with ADCPREP Version 1.7 (  5/24/04 )
C
C                     current for hp_DG_ADCIRC v10_sb4   9/15/2005
C----------------------------------------------------------------------------

      MODULE PRE_GLOBAL
      USE PRESIZES
      IMPLICIT NONE
C
C  PARM14 section
      INTEGER  NELG,NNODG,NOPE,NETA,NBOU,NVEL
      INTEGER,ALLOCATABLE ::  NNEG(:,:)
      INTEGER,ALLOCATABLE ::  NVDLL(:),NBDV(:,:)
      INTEGER,ALLOCATABLE ::  NVELL(:),NBVV(:,:)
      INTEGER,ALLOCATABLE ::  IBTYPE(:)
      INTEGER,ALLOCATABLE ::  IBCONNR(:,:),LBCODE(:)
      INTEGER,ALLOCATABLE ::  WEIR(:),WEIRD(:)

C   METIS interface section
      INTEGER,ALLOCATABLE ::  PROC(:)
C
C  GRID14 section
      REAL(8)  ,ALLOCATABLE ::    X(:),Y(:),SLAM(:),SFEA(:)
      REAL(SZ) ,ALLOCATABLE ::   DP(:)
C
C   BARRIER14 section
      REAL(SZ) ,ALLOCATABLE ::  BAR1(:,:),BAR2(:,:),BAR3(:,:)
 
C   FLOWBC section
      INTEGER NFLBN,NFLBNP
      INTEGER,ALLOCATABLE   ::  FLBN(:),FLBNX(:),FLBNXP(:)
C
C   STRING14 section
      CHARACTER(80) AGRID,NSIZES,SIZEMSG
      CHARACTER(80) NOPEMSG,NETAMSG
      CHARACTER(80) NBOUMSG,NVELMSG
      CHARACTER(80),ALLOCATABLE   ::  NVELLMSG(:),NVDLLMSG(:)
C
C  PARM15  section
      INTEGER  NOUTE,NSPOOLE,NSTAE,NSTAE_MAX
      INTEGER  NOUTV,NSPOOLV,NSTAV,NSTAV_MAX
      INTEGER  NOUTM,NSPOOLM,NSTAM,NSTAM_MAX
      INTEGER  NOUTC,NSPOOLC,NSTAC,NSTAC_MAX
      INTEGER  NFOVER,NABOUT,NSCREEN
      INTEGER  IHOT,ICS,IM
      INTEGER  NOLIBF,NOLIFA,NOLICA,NOLICAT
      INTEGER  NWP,NCOR,NTIP,NWS,NRAMP
      INTEGER  NRS
      INTEGER  NTIF,NBFR,NFFR
      INTEGER  NOUTGE,NSPOOLGE
      INTEGER  NOUTGV,NSPOOLGV
      INTEGER  NOUTGC,NSPOOLGC
      INTEGER  NOUTGW,NSPOOLGW
      INTEGER  NHARFR,NHAINC,NHSINC
      INTEGER  NHASE,NHASV,NHAGE,NHAGV,NHSTAR
      INTEGER  IREFYR,IREFMO,IREFDAY,IREFHR,IREFMIN
      INTEGER  NWLAT,NWLON,NFLUXF
      INTEGER IGC                      ! vertical grid code in fort.15
cvjp 3/23/06 added to limit iteration info in log file
      INTEGER STDOUT
      INTEGER,ALLOCATABLE ::   NNSEG(:),NNSVG(:),NNSCG(:),NNSMG(:)
      
      LOGICAL APERIODIC_FLOW_BC     ! .true. if EXIST_FLUX is .true. .AND.
                                    ! NFFR=0 in fort.15
      
C     jgf46.00 Nodal attributes section.
      CHARACTER(len=80), ALLOCATABLE :: NodalAttributes(:) ! labels
C
C  STRING15 section
C
      CHARACTER(80) RUNDES,RUNID,OVERMSG,ABOUTMSG,SCREENMSG,HOTMSG
      CHARACTER(80) ICSMSG,IMMSG,IBFMSG,IFAMSG,ICAMSG,ICATMSG,NWPMSG,
     &  NCORMSG
      CHARACTER(80) NTIPMSG,NWSMSG,RAMPMSG,GMSG,TAU0MSG,DTMSG,
     &  STATMSG,REFTMSG
      CHARACTER(80) RNDAYMSG,DRAMPMSG,COEFMSG,H0MSG,SLMSG,TAUMSG,
     &  ESLMSG,CORIMSG
      CHARACTER(80)  NTIFMSG,NBFRMSG,ANGMSG,NFFRMSG,STAEMSG,
     &  NSTAEMSG,STAVMSG,RSTIMMSG
      CHARACTER(80)  NSTAVMSG,STACMSG,NSTACMSG,STAMMSG,NSTAMMSG
      CHARACTER(80)  OUTGEMSG,OUTGVMSG,OUTGCMSG,OUTGWMSG
      CHARACTER(80)  HARFRMSG,HARPARMSG,OUTHARMSG,HSTARMSG
      CHARACTER(80)  SOLVMSG,WSMSG1,WSMSG2
      CHARACTER*80,ALLOCATABLE :: STAELOC(:),STAVLOC(:),STACLOC(:)
      CHARACTER*80,ALLOCATABLE :: STAMLOC(:),HAFNAM(:),HAFREMSG(:)
      CHARACTER*80,ALLOCATABLE :: TIPOTAG(:),BOUNTAG(:),FBOUNTAG(:)
      CHARACTER*80,ALLOCATABLE :: ALPHA1(:),ALPHA2(:),FREQMSG(:),
     &  QNMSG(:,:)
      CHARACTER*80,ALLOCATABLE :: AMIGMSG(:),EMOMSG(:,:),TPKMSG(:)
C
C   PARM15-3DVS section
      INTEGER :: IDIAG,ISLIP,ICG,NFEN,IEVC
      INTEGER :: I3DSD,NSPO3DSD       ! NHN3DSD is in presizes
      INTEGER :: I3DSV,NSPO3DSV       ! NHN3DSV is in presizes
      INTEGER :: I3DST,NSPO3DST       ! NHN3DST is in presizes
      REAL(SZ) :: TO3DSSD,TO3DFSD
      REAL(SZ) :: TO3DSSV,TO3DFSV
      REAL(SZ) :: TO3DSST,TO3DFST
      INTEGER,ALLOCATABLE :: ISDHOUT(:),ISVHOUT(:),ISTHOUT(:)
      INTEGER :: I3DGD,NSPO3DGD
      INTEGER :: I3DGV,NSPO3DGV
      INTEGER :: I3DGT,NSPO3DGT
      REAL(SZ) :: TO3DSGD,TO3DFGD
      REAL(SZ) :: TO3DSGV,TO3DFGV
      REAL(SZ) :: TO3DSGT,TO3DFGT
      REAL(SZ) :: KP,Z0S,Z0B,ALP1,ALP2,ALP3,EVMIN,EVCON,THETA1,THETA2
C
C   STRING15-3DVS section
      CHARACTER(80) :: IDIAGMSG,IDENMSG,SLIPMSG,Z0MSG,ALPMSG,FEMSG
      CHARACTER(80) :: EVCMSG,THETAMSG
      CHARACTER(80) :: DSDMSG,DSVMSG,DSTMSG,DGDMSG,DGVMSG,DGTMSG
Csb
C  FORT.HB section
      REAL(SZ) ,ALLOCATABLE ::   DP_HB(:)
C
C   ELESTAT section
      REAL(SZ) TOUTSE,TOUTFE,TOUTSGE,TOUTFGE
      REAL(8),ALLOCATABLE :: XEL(:),YEL(:),SLEL(:),SFEL(:)
C
C   VELSTAT section
      REAL(SZ) TOUTSGV,TOUTFGV,TOUTSV,TOUTFV
      REAL(8),ALLOCATABLE :: XEV(:),YEV(:),SLEV(:),SFEV(:)
C
C   CONSTAT section
      REAL(SZ) TOUTSC,TOUTFC,TOUTSGC,TOUTFGC
      REAL(8),ALLOCATABLE :: XEC(:),YEC(:),SLEC(:),SFEC(:)
C
C   METSTAT section
      REAL(SZ) TOUTSM,TOUTFM,TOUTSGM,TOUTFGM
      REAL(8),ALLOCATABLE ::  XEM(:),YEM(:),SLEM(:),SFEM(:)
      REAL(SZ) TOUTSGW,TOUTFGW
C
C   INPUT15 section
      INTEGER    NODEDRYMIN,NODEWETRMP
      REAL(SZ)   G
      REAL(SZ)   TAU0,RNDAY,DRAMP
      REAL(SZ)   A00,B00,C00,H0,VELMIN
      REAL(8)    SLAM0,SFEA0
      REAL(SZ)   TAU,CF,ESL,CORI,ANGINN
      REAL(SZ)   ESLM,ESLC,HBREAK,FTHETA,FGAMMA
      REAL(8)    REFSEC
      REAL(SZ)   WLATMAX,WLONMIN,WLATINC,WLONINC
      REAL(SZ),ALLOCATABLE ::  FF(:),FACE(:)
      REAL(8) ,ALLOCATABLE ::  EMO(:,:),EFA(:,:)
      REAL(SZ),ALLOCATABLE ::  TPK(:),ETRF(:),FFT(:),FACET(:)
      REAL(SZ),ALLOCATABLE ::  QNAM(:,:),QNPH(:,:)
      REAL(SZ),ALLOCATABLE ::  FFF(:),FFACE(:)
      REAL(SZ),ALLOCATABLE :: Sigma(:) !(NFEN) User specified vertical grid
      REAL(SZ),ALLOCATABLE :: EVTot(:) !(NFEN) Vertical eddy viscosity
C
C
C  INPUT15D section
      REAL(SZ)   DT
      REAL(8)    STATIM,REFTIM,WTIMINC
      REAL(8),ALLOCATABLE ::   AMIG(:),AMIGT(:),FAMIG(:)
C
C
C   HARMON section
      REAL(SZ) THAS,THAF,FMV
      REAL(SZ),ALLOCATABLE ::   HAFREQ(:),HAFF(:),HAFACE(:)
C
C
C   SOLVER  section
      INTEGER ITITER,ISLDIA,ITMAX
      REAL(SZ) CONVCR
C
C--Degress-to-Radians and Radians-to-Degrees
C
C  CONVERT section
      REAL(8) DEG2RAD,RAD2DEG,R
C
C
C--------------------------------------------------------------------------C
C                                                                          C
C              DATA DECOMPOSITION DECLARATIONS BEGIN HERE                  C
C                                                                          C
C--------------------------------------------------------------------------C
C
C--Local Map Variable Declarations
C
C   LOCALI section
      INTEGER NPROC,MSHAR,NWEIR 
      INTEGER,ALLOCATABLE  ::  NELP(:), NNODP(:), NNEP(:,:,:)
Csb-DG1
      INTEGER,ALLOCATABLE  ::  NIEL_SEND(:),NIEL_RECV(:)
      INTEGER,ALLOCATABLE  ::  IEL_SEND(:,:,:),IEL_RECV(:,:,:)
      INTEGER,ALLOCATABLE  ::  NUM_COMM_PE_SEND(:),NUM_COMM_PE_RECV(:)
      INTEGER,ALLOCATABLE  ::  IMAP_NEIGHEDG(:,:,:) ! Maps (IEDG,IEL) --> (IEDG,IEL) of the neighbor edge
C--
Csb-02/26/06
      INTEGEr              ::  MNNDEL
      INTEGER,ALLOCATABLE  ::  NNDEL(:),NDEL(:,:)
C--
      INTEGER,ALLOCATABLE  ::  NOD_RES_TOT(:),NSTACP(:), NSTAMP(:)
      INTEGER,ALLOCATABLE  ::  NWEIRP(:),NSTAEP(:),NSTAVP(:)
      INTEGER,ALLOCATABLE  ::  NOPEP(:),NETAP(:),NVDLLP(:)
      INTEGER,ALLOCATABLE  ::  NBDVP(:,:)
      INTEGER,ALLOCATABLE  ::  NBOUP(:),NVELP(:),NVELLP(:)
      INTEGER,ALLOCATABLE  ::  NBVVP(:,:), IBTYPEP(:,:)
      INTEGER,ALLOCATABLE  ::  LBCODEP(:,:)
      INTEGER,ALLOCATABLE  ::  ITOTPROC(:),NFLUXFP(:)
      INTEGER,ALLOCATABLE  ::  IBCONNRP(:,:)
C
C
C DIAGS section
      REAL(SZ),ALLOCATABLE ::   PROC_SV(:)
C
C--Local-to-Global Mapping Variables
C
C   LOC2G section
      INTEGER,ALLOCATABLE  ::   IMAP_NOD_LG(:,:),IMAP_EL_LG(:,:)
      INTEGER,ALLOCATABLE  ::   EL_SHARE(:)
      INTEGER,ALLOCATABLE  ::   IMAP_STAE_LG(:,:), STAE_SHARE(:)
      INTEGER,ALLOCATABLE  ::   IMAP_STAV_LG(:,:), STAV_SHARE(:)
      INTEGER,ALLOCATABLE  ::   IMAP_STAC_LG(:,:), STAC_SHARE(:)
      INTEGER,ALLOCATABLE  ::   IMAP_STAM_LG(:,:), STAM_SHARE(:)
      INTEGER,ALLOCATABLE  ::   OBNODE_LG(:,:)
      INTEGER,ALLOCATABLE  ::   WEIRP_LG(:,:)
      INTEGER,ALLOCATABLE  ::   WEIRDP_LG(:,:)
      INTEGER,ALLOCATABLE  ::   LBINDEX_LG(:,:)
C
C--Global-to-Local Mapping Variables
C
C   GLOB2L section
      INTEGER,ALLOCATABLE  ::   IMAP_NOD_GL(:,:),IMAP_NOD_GL2(:,:)

C
C--Message-Passing Variables
C
C   MSGTAB section
      INTEGER,ALLOCATABLE  ::   NUM_COMM_PE(:), COMM_PE_NUM(:,:)
      INTEGER,ALLOCATABLE  ::   IRECV_TOT(:,:),IRECV(:)
      INTEGER,ALLOCATABLE  ::   ISEND_TOT(:,:),ISEND(:)
C
C--3DVS section
C
      INTEGER,ALLOCATABLE  ::  N3DSDP(:),N3DSVP(:),N3DSTP(:)
      INTEGER,ALLOCATABLE  ::  IMAP_3DSD_LG(:,:)
      INTEGER,ALLOCATABLE  ::  IMAP_3DSV_LG(:,:)
      INTEGER,ALLOCATABLE  ::  IMAP_3DST_LG(:,:)

C-------------------end of data declarations------------------------C


      CONTAINS

      SUBROUTINE ALLOC_MAIN1()
C
C     Allocate space for Arrays except those dimensioned by MNPP and MNEP
C
      ALLOCATE ( PROC(MNP) )
      ALLOCATE ( X(MNP),Y(MNP),DP(MNP),SLAM(MNP),SFEA(MNP) )
Csb-For hard bottom
      IF(HBFLAG.EQ.1) ALLOCATE ( DP_HB(MNP) )
C--
      ALLOCATE ( NNEG(3,MNE) )
      ALLOCATE ( NVDLL(MNOPE),NBDV(MNOPE,MNETA) )
      ALLOCATE ( NVELL(MNBOU),NBVV(MNBOU,0:MNVEL),IBTYPE(MNBOU) )
      ALLOCATE ( IBCONNR(MNBOU,MNVEL),LBCODE(MNVEL) )
      ALLOCATE ( WEIR(MNVEL),WEIRD(MNVEL) )
      ALLOCATE ( BAR1(MNBOU,MNVEL),BAR2(MNBOU,MNVEL),BAR3(MNBOU,MNVEL) )
      ALLOCATE ( FLBN(MNVEL),FLBNX(MNVEL),FLBNXP(MNVEL) )
      ALLOCATE ( NVDLLMSG(MNOPE),NVELLMSG(MNBOU+1) )
      ALLOCATE ( NNSEG(MNSTAE),NNSVG(MNSTAV),NNSCG(MNSTAC),
     &           NNSMG(MNSTAM) )
      ALLOCATE ( STAELOC(MNSTAE),STAVLOC(MNSTAV),STACLOC(MNSTAC) )
      ALLOCATE ( STAMLOC(MNSTAM) )
      ALLOCATE ( HAFNAM(MNHARF),HAFREMSG(MNHARF) )
      ALLOCATE ( TIPOTAG(MNTIF),BOUNTAG(MNBFR),FBOUNTAG(MNFFR) )
      ALLOCATE ( ALPHA1(MNBFR),ALPHA2(MNFFR),FREQMSG(MNFFR),
     &           QNMSG(MNFFR,MNVEL) )
      ALLOCATE ( AMIGMSG(MNBFR),EMOMSG(MNBFR,MNETA),TPKMSG(MNTIF) )
      ALLOCATE ( XEL(MNSTAE),YEL(MNSTAE),SLEL(MNSTAE),SFEL(MNSTAE) )
      ALLOCATE ( XEV(MNSTAV),YEV(MNSTAV),SLEV(MNSTAV),SFEV(MNSTAV) )
      ALLOCATE ( XEC(MNSTAC),YEC(MNSTAC),SLEC(MNSTAC),SFEC(MNSTAC) )
      ALLOCATE ( XEM(MNSTAM),YEM(MNSTAM),SLEM(MNSTAM),SFEM(MNSTAM) )
      ALLOCATE ( FF(MNBFR),FACE(MNBFR) )
      ALLOCATE ( EMO(MNBFR,MNETA),EFA(MNBFR,MNETA) )
      ALLOCATE ( TPK(MNTIF),ETRF(MNTIF),FFT(MNTIF),FACET(MNTIF) )
      ALLOCATE ( QNAM(MNFFR,MNVEL),QNPH(MNFFR,MNVEL) )
      ALLOCATE ( FFF(MNFFR),FFACE(MNFFR) )
      ALLOCATE ( AMIG(MNBFR),AMIGT(MNTIF), FAMIG(MNFFR) )
      ALLOCATE ( HAFREQ(MNHARF),HAFF(MNHARF),HAFACE(MNHARF) )
      ALLOCATE ( NELP(MNPROC), NNODP(MNPROC) )
      ALLOCATE ( NOD_RES_TOT(MNPROC) )
      ALLOCATE ( NWEIRP(MNPROC) )
      ALLOCATE ( NSTAEP(MNPROC),NSTAVP(MNPROC),NSTACP(MNPROC) )
      ALLOCATE ( NSTAMP(MNPROC) )
      ALLOCATE ( NOPEP(MNPROC),NETAP(MNPROC),NVDLLP(MNOPE) )
      ALLOCATE ( NBDVP(MNOPE,MNETA) )
      ALLOCATE ( NBOUP(MNPROC),NVELP(MNPROC),NVELLP(MNBOU) )
      ALLOCATE ( NBVVP(MNBOU,0:MNVEL), IBTYPEP(MNBOU,MNPROC) )
      ALLOCATE ( LBCODEP(MNVEL,MNPROC) )
      ALLOCATE ( ITOTPROC(MNP),NFLUXFP(MNPROC) )
      ALLOCATE ( IBCONNRP(MNBOU,MNVEL) )
      ALLOCATE ( PROC_SV(MNPROC) )
      ALLOCATE ( IMAP_STAE_LG(MNSTAE,MNPROC) )
      ALLOCATE ( IMAP_STAV_LG(MNSTAV,MNPROC) )
      ALLOCATE ( IMAP_STAC_LG(MNSTAC,MNPROC) )
      ALLOCATE ( IMAP_STAM_LG(MNSTAM,MNPROC) )
      ALLOCATE ( OBNODE_LG(MNETA,MNPROC) )
      ALLOCATE ( WEIRP_LG(MNVEL,MNPROC) )
      ALLOCATE ( WEIRDP_LG(MNVEL,MNPROC) )
      ALLOCATE ( LBINDEX_LG(MNBOU,MNVEL) )
      ALLOCATE ( IMAP_NOD_GL(2,MNP),IMAP_NOD_GL2(2*MNEI,MNP) )
      ALLOCATE ( NUM_COMM_PE(MNPROC), COMM_PE_NUM(MNPROC,MNPROC) )
      ALLOCATE ( IRECV_TOT(MNPROC,MNPROC),IRECV(MNP) )
      ALLOCATE ( ISEND_TOT(MNPROC,MNPROC),ISEND(MNP) )
      ALLOCATE ( STAE_SHARE(MNSTAE), STAV_SHARE(MNSTAV) )
      ALLOCATE ( STAM_SHARE(MNSTAM) )
      ALLOCATE ( STAC_SHARE(MNSTAC))
      ALLOCATE(EL_SHARE(MNEP))
C
C--Allocate 3DVS stuff
      IF(NHN3DSD.GT.0)
     &  ALLOCATE ( N3DSDP(MNPROC), IMAP_3DSD_LG(NHN3DSD,MNPROC) )
      IF(NHN3DSV.GT.0)
     &  ALLOCATE ( N3DSVP(MNPROC), IMAP_3DSV_LG(NHN3DSV,MNPROC) )
      IF(NHN3DST.GT.0)
     &  ALLOCATE ( N3DSTP(MNPROC), IMAP_3DST_LG(NHN3DST,MNPROC) )
C
C
      RETURN
      END SUBROUTINE


      SUBROUTINE ALLOC_MAIN2()
C
C     Allocate space for Arrays dimensioned by MNPP and MNEP
C
      ALLOCATE ( IMAP_NOD_LG(MNPP,MNPROC),IMAP_EL_LG(MNEP,MNPROC) )
      ALLOCATE ( NNEP(3,MNEP,MNPROC) )
C
      RETURN
      END SUBROUTINE


Csb-DG1
      SUBROUTINE ALLOC_DG1()

      ALLOCATE ( NIEL_SEND(MNPROC) )
      ALLOCATE ( NIEL_RECV(MNPROC) )
      ALLOCATE ( NUM_COMM_PE_SEND(MNPROC) )
      ALLOCATE ( NUM_COMM_PE_RECV(MNPROC) )
      ALLOCATE ( IMAP_NEIGHEDG(2,3,NELG) )
Csb-02/26/06
      ALLOCATE ( NNDEL(NELG) )
C--
      RETURN
      END SUBROUTINE
C--


Csb-02/26/06
      SUBROUTINE ALLOC_DG1B()

      ALLOCATE ( NDEL(MNP,MNNDEL) )

      RETURN
      END SUBROUTINE
C--


Csb-DG1
      SUBROUTINE ALLOC_DG2()

      INTEGER IPROC, MNIEL_SEND, MNIEL_RECV
      
      MNIEL_SEND = 0
      MNIEL_RECV = 0

      DO IPROC = 1,MNPROC
        IF(NIEL_SEND(IPROC).GT.MNIEL_SEND) THEN
          MNIEL_SEND = NIEL_SEND(IPROC)
        ENDIF
        IF(NIEL_RECV(IPROC).GT.MNIEL_RECV) THEN
          MNIEL_RECV = NIEL_RECV(IPROC)
        ENDIF
      ENDDO

      ALLOCATE ( IEL_SEND(2,MNIEL_SEND,MNPROC) )
      ALLOCATE ( IEL_RECV(2,MNIEL_RECV,MNPROC) )

      RETURN
      END SUBROUTINE
C--

C
C
C--------------------------------------------------------------------------C
C                                                                          C
C                DEFINITIONS OF DOMAIN DECOMPOSITION VARIABLES             C
C                                                                          C
C--------------------------------------------------------------------------C
C
C Processing Element Definitions:
C
C   MNPROC               = Maximum Number of PEs
C   NPROC                = Actual Number of PEs - this should eventually be
C                          dropped and MNPROC used throughout preprocessor
C                          routines - RL
C   MSHAR                = Max. Number PEs assigned to any Global Node
C
C Nodal and Element Definitions:
C
C   X(I)                 = X-coordinate of Global Node I
C   Y(I)                 = Y-coordinate of Global Node I
C   DP(I)                = Bathymetry of Global Node I
C
C   NELG                 = Number of Global Elements
C   NELP(PE)             = Number of Elements Assigned to PE
C
C   NNODG                = Number of Global Nodes
C   NNODP(PE)            = Number of Nodes Assigned to PE
C
C   NOD_RES_TOT(PE)      = Number of Resident Nodes on PE
C   ITOTPROC(I)          = Number of PEs assigned to Global Node I
C
C   NNEG(3,I)            = Three Nodes of Global Element I
C   NNEP(3,I,PE)         = Three Nodes of Element I on PE
C 
C   IMAP_NOD_GL(1,I)     = PE assigned to Global Node I
C   IMAP_NOD_GL(2,I)     = Local Node Number of Global Node I
C
C   IMAP_NOD_LG(I,PE)    = Global Node Number of Local Node I on PE
C   IMAP_EL_LG(I,PE)     = Global Element Number of Local Element I on PE
C
C   IMAP_NOD_GL2(2(PE-1)+1,I)  = PE assigned to Global Node I
C   IMAP_NOD_GL2(2(PE-1)+2,I)  = Local Node Number of Global Node I on PE
C
C Open Boundary Nodes and Segment Definitions:
C
C   NETA                 = Number of Global Open Boundary Nodes 
C   NETAP(PE)            = Number of Open Boundary Nodes on PE
C
C   NOPE                 = Number of Global Open Boundary Segments 
C   NOPEP(PE)            = Number of Open Boundary Segments on PE
C
C   NVDLL(K)             = Number of Nodes on Global Open Boundary Segment K 
C   NVDLLP(K)            = Number of Nodes on Open Boundary Segment K on PE
C
C   NBDV(K,I)            = Global Node Number of I-th Node on Open Boundary
C                          Segment K
C   NBDVP(K,I)           = Local Node Number of I-th Node on Open Boundary 
C                          Segment K on PE
C
C   OBNODE_LG(I,PE)      = Global Open Boundary Node Number of Local
C                          Open Boundary Node I on PE
C 
C Land Boundary Nodes and Segment Definitions:
C
C   NVEL                 = Total Number of Global Land Boundary Nodes 
C   NVELP(PE)            = Total Number of Land Boundary Nodes on PE
C
C   NBOU                 = Number of Global Land Boundary Segments
C   NBOUP(PE)            = Number of Land Boundary Segments on PE
C    
C   NVELL(K)             = Number of Global Land Boundary Nodes of Segment K
C   NVELLP(K)            = Number of Land Boundary Nodes of Segment K on PE
C
C   NBVV(K,I)            = Global Node Number of I-th Node on Land Boundary
C                          Segment K
C   NBVVP(K,I)           = Local Node Number of I-th Node on Land Boundary
C                          Segment K on PE
C
C   IBTYPE(K)            = Type Global Land Boundary Segment K 
C   IBTYPEP(K,PE)        = Type Land Boundary Segment K on PE     
C                          ( 0=mainland, 1=island )
C
C   LBCODE(I)            = Boundary Type of Land Boundary Node I
C   LBCODEP(I,PE)        = Boundary Type of Land Boundary Node I on PE
C
C   NFLBN                = Number of Global Flow Boundary Nodes
C   FLBN(I)              = Global Node number of Ith Flow Boundary Node
C   FLBNX(I)             = Index of Ith Flow Boundary Node
C
C   NFLBNP(PE)           = Number of Flow boundary Nodes on PE
C   FLBNXP(I)            = Index of Ith Flow Boundary Node
C
Cvjp modified array to drop last dimension to save memory space
C   LBINDEX_LG(K,I,PE)   = Global Index of I-th Node on Land Boundary Segment
C                          K on PE 
C
C   NWEIR                = Total Number of Global Weir Land Boundary Pairs 
C   NWEIRP(PE)           = Total Number of Land Boundary Nodes on PE
C
C   WEIRP_LG(I,PE)       = Global Node Number of I-th Weir Node on PE       
C   WEIRDP_LG(I,PE)      = Global Node Number of I-th Dual Weir Node on PE    
C
C Elevation Station Definitions:
C
C   NSTAEP(PE)           = Number of Elevation Stations on PE
C   IMAP_STAE_LG(I,PE)   = Global Number of Local Elevation Station I on PE
C
C Velocity Station Definitions:
C
C   NSTAVP(PE)           = Number of Velocity Stations on PE
C   IMAP_STAV_LG(I,PE)   = Global Number of Local Velocity  Station I on PE
C
C Concentration Station Definitions:
C
C   NSTACP(PE)           = Number of Concentration Stations on PE
C   IMAP_STAC_LG(I,PE)   = Global Number of Local Concentration Station I on PE
C
C Meterological Station Definitions:
C
C   NSTAMP(PE)           = Number of Meterological Stations on PE
C   IMAP_STAM_LG(I,PE)   = Global Number of Local Meterological Station I on PE
C
C
C Message-Passing Definitions:
C
C   NUM_COMM_PE(PE)      = Number of PEs communicating with PE
C   COMM_PE_NUM(IPE,PE)  = IPE-th PE communicating with PE
C
C   IRECV_TOT(IPE,PE)    = Number of Nodes Received by PE from IPE
C   ISEND_TOT(IPE,PE)    = Number of Nodes Sent by PE to IPE
C
C   PROC_SV(PE)          = Surface-to-Volume Ratio on PE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      END MODULE PRE_GLOBAL 
