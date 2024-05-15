C**************************************************************************
C  mod history
C  v41.06mxxx  - date - programmer - describe change 
C                    - mark change in code with  cinitials-mxxx
C  v41.10      - 07/25/01 - rl - from 41.09 - bug fix in GWCE lateral viscosity term
C  v41.09      - 06/30/01 - jw - from 41.08 - made minor mods as per vp version 41.05
C  v41.06      - 04/02/01 - rl - changed MNWP to MNP in wind forcing 
C                                ALLOCATION statements
C  v41.02      - 09/04 - rl
C  v40.02m004b - 05/17 - vjp - corrected dimensioning problem cvjpm004b
C  v40.02m001  - 12/21 - jjw - add cross barrier pipes cjjwm001
C-----------------------------------------------------------------------
C
C     v9_sb2.2.1 - 08/16/05 - sb - Hard bottom implementation
C     v10_sb5    - 10/11/05 - sb - Table that associates nodes to elements (NODEELEM)
C                - 10/27/05 - sb - H0L and H0H are added
C
C**************************************************************************
C 
      MODULE GLOBAL
      USE SIZES
C...
C...SET GLOBAL PARAMETER CONSTANTS
C...

#ifdef SWAN
Casey 121128: Added variables for the output of global files from SWAN.
      INTEGER,ALLOCATABLE,TARGET :: NODES_LG(:)
      INTEGER,PARAMETER :: BUFSIZE_MAX = 131072
      INTEGER :: float_type
      INTEGER :: integerBuffer(BUFSIZE_MAX)
      INTEGER :: integerResultBuffer(BUFSIZE_MAX)
      INTEGER :: NP_G
      REAL(SZ) :: buf(BUFSIZE_MAX)
      REAL(SZ) :: resultBuf(BUFSIZE_MAX)
      TYPE OutputDataDescript_t
        CHARACTER(12) :: file_name
        CHARACTER(20) :: field_name
        INTEGER,POINTER :: iarray(:)
        INTEGER,POINTER :: iarray_g(:)
        INTEGER,POINTER :: imap(:)
        INTEGER :: int_initial_value
        INTEGER :: num_fd_records
        INTEGER :: num_items_per_record
        INTEGER :: num_records_this
        INTEGER :: specifier
        LOGICAL :: ConsiderWetDry
        REAL(SZ),POINTER :: array(:)
        REAL(SZ),POINTER :: array2(:)
        REAL(SZ),POINTER :: array3(:)
        REAL(SZ),POINTER :: array_g(:)
        REAL(SZ),POINTER :: array2_g(:)
        REAL(SZ),POINTER :: array3_g(:)
        REAL(SZ),POINTER :: hotstart(:)
        REAL(SZ),POINTER :: hotstart_g(:)
        REAL(SZ) :: initial_value
        real(SZ) :: alternate_value
      ENDTYPE OutputDataDescript_t
#endif

      Logical vertexslope

C.....Background Atmospheric Pressure in mb
      REAL(SZ), PARAMETER :: PRBCKGRND = 1013.0d0   !tcm 20100617 v49.16 Added
      REAL(8) :: rhoAir = 1.15d0         ! kg/m^3
      REAL(8) :: windReduction = 0.9d0
      REAL(8) :: one2ten = 0.8928d0
      REAL(8) :: ten2one = 1.d0/0.8928d0
Casey 121019: Added multiplication factor to be used before sending winds to coupled wave models.
      REAL(SZ) :: WaveWindMultiplier = 1.D0
C.....PI etc
      REAL(8), PARAMETER ::  PI=3.141592653589793D0
      REAL(8), PARAMETER ::  DEG2RAD = PI/180.D0  ! degrees to radians
      REAL(8), PARAMETER ::  RAD2DEG = 180.D0/PI  ! radians to degrees
      REAL(8), PARAMETER ::  e=2.718281828459045d0
C.....Days to seconds conversion
      REAL(8), PARAMETER :: day2sec=24.d0*3600.d0
      REAL(8), PARAMETER :: sec2day=1.d0/day2sec
      ! jgf50.38.03: Additional useful constants
      REAL(8), PARAMETER :: Rearth = 6378206.4d0 ! radius of earth (m)
      REAL(8), PARAMETER :: nm2m = 1852.d0 ! nautical miles to meters
      REAL(8), PARAMETER :: m2nm = 1.d0/nm2m ! meters to nautical miles
      REAL(8), PARAMETER :: kt2ms = nm2m / 3600.0 ! knots to m/s
      REAL(8), PARAMETER :: ms2kt = 1.d0 / kt2ms ! m/s to knots
      REAL(8), PARAMETER :: omega = 2.0d0*pi / 86164.2d0
      REAL(8), PARAMETER :: mb2pa = 100.d0

C.....parameters used in barrier overflow
      REAL(SZ), PARAMETER ::  BARMIN=0.01D0
      REAL(SZ) DEPAVG,DEPMAX,DEPMIN
C
C.....Sediment Transport stuff added (Ethan Kubatko 8-1-2003)
C
C.....Declare variables used in sediment transport section
C
      INTEGER, TARGET :: SEDFLAG
      INTEGER MAXEL, ELEM_ED, NBOR_ED, NBOR_EL,ITDG,ModetoNode
      Integer tracer_flag, chem_flag
      INTEGER N1,N2,NO_NBORS,NBOR,SEDFLAG_W, OPEN_INDEX
      INTEGER, TARGET :: DG_TO_CG
      INTEGER NSCREEN_INC
      
      INTEGER ScreenUnit
      
      REAL(SZ) AREA_SUM, CEN_SUM, NLEQ, LEQ, NLEQG
      REAL(SZ), TARGET :: reaction_rate
      Character*100, TARGET :: sed_equationX,sed_equationY
      
      REAL(SZ) FluxSettlingTime
      INTEGER  FluxSettlingIT
C     jgf46.08 Fine grained ramp functions (jgf46.21 split flux b.c.s)
      REAL(SZ) RampExtFlux,DRampExtFlux ! Ramp for external flux b.c.s
      REAL(SZ) RampIntFlux,DRampIntFlux ! Ramp for internal flux b.c.s
      REAL(SZ) RampElev,DRampElev    ! Ramp for elevation boundary conditions.
      REAL(SZ) RampTip,DRampTip      ! Ramp for tidal potential
      REAL(SZ) RampMete,DRampMete    ! Ramp for wind and atmospheric pressure
      REAL(SZ) RampWRad,DRampWRad    ! Ramp for wave radiation stress

C.....Create allocatable arrays

      REAL(SZ),ALLOCATABLE :: ANGTAB(:,:),CENTAB(:,:),ELETAB(:,:)
      REAL(SZ),ALLOCATABLE :: DG_ANG(:),DP_DG(:)
      INTEGER, ALLOCATABLE :: EL_COUNT(:)
      INTEGER, ALLOCATABLE :: NNOEL(:,:)

Csb...Table that maps nodes to elements
      INTEGER, ALLOCATABLE :: NNDEL(:)
      INTEGER, ALLOCATABLE :: NDEL(:,:)
C--
      
C.....Declare variables for DG SW

      INTEGER, TARGET :: DGSWE
      INTEGER EL_IN, EL_EX, SD_IN, SD_EX, EDGE(3)
      INTEGER SIDE(2),TESTPROBLEM
      REAL(SZ) FX_IN,FY_IN,GX_IN,GY_IN,HX_IN,HY_IN     
      REAL(SZ) FX_EX,FY_EX,GX_EX,GY_EX,HX_EX,HY_EX     
      REAL(SZ) F_AVG,G_AVG,H_AVG,JUMP(4),HT_IN,HT_EX
      REAL(SZ) C_ROE,U_ROE,V_ROE,EIGVAL(4),RI(4,4),LE(4,4),A_ROE(4,4)
      Real(SZ) UMag_IN, UMag_EX, ZE_ROE,QX_ROE,QY_ROE,bed_ROE
      REAL(SZ) Q_N,Q_T,U_N,U_T,U_IN,U_EX,V_IN,V_EX
      REAL(SZ) ZE_SUM,QX_SUM,QY_SUM,DG_MAX,DG_MIN,U_N_EXT,U_T_EXT
      REAL(SZ) Q_N_INT,Q_T_INT,U_N_INT,U_T_INT,Q_N_EXT,Q_T_EXT

      REAL(SZ) BX_INT,BY_INT,SOURCE_1,SOURCE_2,SOURCE_SUM,k_hat
      REAL(SZ) FRIC_AVG, DP_MID, F_HAT, G_HAT, H_HAT,i_hat,j_hat
      REAL(SZ) INFLOW_ZE,INFLOW_QX,INFLOW_QY,H_LEN,INFLOW_LEN
      REAL(SZ) ZE_NORM,QX_NORM,QY_NORM,ZE_DECT,QX_DECT,QY_DECT
      REAL(SZ),ALLOCATABLE :: FX_MID(:,:),GX_MID(:,:),HX_MID(:,:)
      REAL(SZ),ALLOCATABLE :: FY_MID(:,:),GY_MID(:,:),HY_MID(:,:)
      REAL(SZ),ALLOCATABLE :: ZE_C(:),QX_C(:),QY_C(:),dynP_DG(:)
      REAL(SZ),ALLOCATABLE :: ZE_DG(:),QX_DG(:),QY_DG(:),HB_DG(:)
      REAL(SZ),Allocatable :: iota2_DG(:),iota_DG(:),iotaa_DG(:)
      REAL(SZ),Allocatable :: bed_DG(:,:), bed_N_int(:), bed_N_ext(:)
C...
C...DECLARE ALL ARRAYS
C...
#ifdef CMPI
      INTEGER, ALLOCATABLE ::   IDUMY(:)
      REAL(SZ),ALLOCATABLE ::   DUMY1(:),DUMY2(:)
      REAL(SZ),ALLOCATABLE ::   DGDUMY1(:,:,:),DGDUMY2(:,:,:)
#endif

      INTEGER, ALLOCATABLE ::   pdg_el(:)
      REAL(SZ),ALLOCATABLE ::   ETAS(:),ETA1(:),ETA2(:),ETAMAX(:)
      REAL(SZ),ALLOCATABLE ::   entrop(:,:)
      Real(SZ),Allocatable::tracer(:),tracer2(:),MassMax(:),bed_int(:,:)
      REAL(SZ),ALLOCATABLE ::   UU1(:),UU2(:),VV1(:),VV2(:), dyn_P(:)
      REAL(SZ),ALLOCATABLE ::   DP(:),DP0(:),DPe(:),SFAC(:) !,STARTDRY(:)
      REAL(SZ),ALLOCATABLE ::   QU(:),QV(:),QW(:)
      REAL(SZ),ALLOCATABLE ::   CORIF(:) !,EVM(:)
      REAL(SZ),ALLOCATABLE ::   TPK(:),FFT(:)
      REAL(SZ),ALLOCATABLE ::   FACET(:),ETRF(:)
      REAL(SZ),ALLOCATABLE ::   ESBIN1(:),ESBIN2(:)
      REAL(SZ),ALLOCATABLE ::   QTEMA(:,:),QTEMB(:,:),QN2(:)
      REAL(SZ),ALLOCATABLE ::   BNDLEN2O3(:)
      REAL(SZ),ALLOCATABLE ::   CSII(:),SIII(:)
      REAL(SZ),ALLOCATABLE ::   QNAM(:,:),QNPH(:,:)
      REAL(SZ),ALLOCATABLE ::   QNIN1(:),QNIN2(:)
      REAL(SZ),ALLOCATABLE ::   CSI(:),SII(:)
C      REAL(SZ),ALLOCATABLE ::   TAU0VAR(:)
      REAL(SZ),ALLOCATABLE ::   ET00(:),BT00(:)
      REAL(SZ),ALLOCATABLE ::   STAIE1(:),STAIE2(:),STAIE3(:)
      REAL(8),ALLOCATABLE ::    XEV(:),YEV(:),SLEV(:),SFEV(:)
      REAL(SZ),ALLOCATABLE ::   UU00(:),VV00(:)
      REAL(SZ),ALLOCATABLE ::   STAIV1(:),STAIV2(:),STAIV3(:)
      REAL(8),ALLOCATABLE ::    XEC(:),YEC(:),SLEC(:),SFEC(:)
      REAL(SZ),ALLOCATABLE ::   CC00(:)
      REAL(SZ),ALLOCATABLE ::   STAIC1(:),STAIC2(:),STAIC3(:)
      REAL(8),ALLOCATABLE ::    XEM(:),YEM(:),SLEM(:),SFEM(:)
      REAL(SZ),ALLOCATABLE ::   RMU00(:),RMV00(:),RMP00(:),RMR00(:)
      REAL(SZ),ALLOCATABLE ::   STAIM1(:),STAIM2(:),STAIM3(:)
      REAL(SZ),ALLOCATABLE ::   CH1(:),QB(:),QA(:),SOURSIN(:)!,EVC(:)
      REAL(SZ),ALLOCATABLE ::   WSX1(:),WSY1(:),PR1(:)
      REAL(SZ),ALLOCATABLE ::   WSX2(:),WSY2(:),PR2(:)
      REAL(SZ),ALLOCATABLE ::   WVNX1(:),WVNY1(:),PRN1(:)
      REAL(SZ),ALLOCATABLE ::   WVNX2(:),WVNY2(:),PRN2(:)
      REAL(SZ),ALLOCATABLE ::   RSNX1(:),RSNY1(:),RSNX2(:),RSNY2(:)
#ifdef SWAN
Casey 101118: Added the following arrays for output of radiation stress gradients.
      REAL(SZ),ALLOCATABLE,TARGET :: RSNXOUT(:), RSNYOUT(:)
#endif
C.....Allocate additional variables for wave friction
      REAL(SZ),ALLOCATABLE ::WAVE_T1(:),WAVE_H1(:),WAVE_A1(:),WAVE_D1(:)
      REAL(SZ),ALLOCATABLE ::WAVE_T2(:),WAVE_H2(:),WAVE_A2(:),WAVE_D2(:)
      REAL(SZ),ALLOCATABLE ::WAVE_T(:), WAVE_H(:), WAVE_A(:), WAVE_D(:)
      REAL(SZ),ALLOCATABLE ::   WB(:)
      REAL(SZ),ALLOCATABLE ::   WVNXOUT(:),WVNYOUT(:)
      REAL(SZ),ALLOCATABLE ::   TKXX(:),TKYY(:),TKXY(:)
      REAL(SZ),ALLOCATABLE ::    EMO(:,:),EFA(:,:)
      REAL(SZ),ALLOCATABLE ::    UMO(:,:),UFA(:,:),VMO(:,:),VFA(:,:)
      REAL(SZ),ALLOCATABLE ::    XEL(:),YEL(:),SLEL(:),SFEL(:)
      REAL(SZ),ALLOCATABLE ::    AREAS(:),SFACDUB(:,:),YDUB(:,:,:)

      !temp
      REAL(8) ,ALLOCATABLE ::   RTEMP2(:)


      REAL(SZ),ALLOCATABLE ::   XVELAV(:),YVELAV(:),XVELVA(:),YVELVA(:)
      REAL(SZ),ALLOCATABLE ::   ELAV(:),ELVA(:)
      REAL(SZ),ALLOCATABLE ::   AUV11(:),AUV12(:),AUV13(:),AUV14(:)
      REAL(SZ),ALLOCATABLE ::   AUVXX(:), AUVYY(:), AUVXY(:), AUVYX(:)
      REAL(SZ),ALLOCATABLE ::   DUU1(:),DUV1(:),DVV1(:),BSX1(:),BSY1(:)
      REAL(SZ),ALLOCATABLE ::   TIP1(:),TIP2(:)
      REAL(SZ),ALLOCATABLE ::   SALTAMP(:,:),SALTPHA(:,:)
      REAL(SZ),ALLOCATABLE ::   OBCCOEF(:,:),COEF(:,:)
      REAL(SZ),ALLOCATABLE ::   WKSP(:),RPARM(:)
      REAL(SZ),ALLOCATABLE ::   ABD(:,:),ZX(:)
      REAL(SZ),ALLOCATABLE ::   GRAVX(:),GRAVY(:)

      INTEGER,ALLOCATABLE ::    ME2GW(:)
      INTEGER,ALLOCATABLE ::    NBV(:),LBCODEI(:)
      INTEGER,ALLOCATABLE ::    NNODECODE(:),NODECODE(:),NODEREP(:)
      INTEGER,ALLOCATABLE ::    NIBCNT(:)
      INTEGER,ALLOCATABLE ::    NM(:,:)
      INTEGER,ALLOCATABLE ::    NNEIGH(:),MJU(:),NODELE(:)
      INTEGER,ALLOCATABLE ::    NEITAB(:,:),NNEIGH_ELEM(:)
      INTEGER,ALLOCATABLE ::    NIBNODECODE(:),NEIGH_ELEM(:,:)
      INTEGER,ALLOCATABLE ::    LBCODE(:)
      INTEGER,ALLOCATABLE ::    NNC(:)
      INTEGER,ALLOCATABLE ::    NNE(:)
      INTEGER,ALLOCATABLE ::    NNV(:)
      INTEGER,ALLOCATABLE ::    NNM(:)
      INTEGER,ALLOCATABLE ::    IWKSP(:),IPARM(:),IPV(:)
      INTEGER,ALLOCATABLE ::    NVDLL(:),NBD(:)
      INTEGER,ALLOCATABLE ::    NBDV(:,:)
      INTEGER,ALLOCATABLE ::    NVELL(:)
      INTEGER,ALLOCATABLE ::    NBVV(:,:)
      INTEGER,ALLOCATABLE ::    NELED(:,:)
      INTEGER,ALLOCATABLE ::    SEGTYPE(:)
      INTEGER,ALLOCATABLE ::    NOT_AN_EDGE(:)
      INTEGER,ALLOCATABLE ::    WEIR_BUDDY_NODE(:,:)
      INTEGER,ALLOCATABLE ::    ONE_OR_TWO(:)
      INTEGER,ALLOCATABLE ::    EDFLG(:,:)

C.....for internal barrier boundaries with flowthrough pipes
      REAL(SZ),ALLOCATABLE ::   BARLANHTR(:),BARLANCFSPR(:)
      REAL(SZ),ALLOCATABLE ::   BARINHTR(:),BARINCFSBR(:),BARINCFSPR(:)
      REAL(SZ),ALLOCATABLE ::   PIPEHTR(:),PIPECOEFR(:),PIPEDIAMR(:)
      REAL(SZ),ALLOCATABLE ::   BARLANHT(:),BARLANCFSP(:)
      REAL(SZ),ALLOCATABLE ::   FFF(:),FFACE(:)
      REAL(SZ),ALLOCATABLE ::   BTRAN3(:),BTRAN4(:),BTRAN5(:)
      REAL(SZ),ALLOCATABLE ::   BTRAN6(:),BTRAN7(:),BTRAN8(:)
      REAL(SZ),ALLOCATABLE ::   BARINHT(:),BARINCFSB(:),BARINCFSP(:)
      REAL(SZ),ALLOCATABLE ::   PIPEHT(:),PIPECOEF(:),PIPEDIAM(:)
      REAL(SZ),ALLOCATABLE ::   RBARWL1AVG(:),RBARWL2AVG(:)
      REAL(SZ),ALLOCATABLE ::   ELEXLEN(:,:)
      INTEGER, ALLOCATABLE ::   IBCONN(:),IBCONNR(:),NTRAN1(:),NTRAN2(:)

C.....for bridge pilings
      REAL(SZ)                  POAN,Fr,FRICBP
      INTEGER                   NBPNODES
      REAL(SZ),ALLOCATABLE ::   BK(:),BALPHA(:),BDELX(:)
      INTEGER, ALLOCATABLE ::   NBNNUM(:)

C.....Younghun added:
      INTEGER, TARGET :: DGSAKE
      INTEGER, TARGET :: RAINFALL
      REAL(SZ),ALLOCATABLE :: H00(:), HQX00(:), HQY00(:)
      REAL(SZ),ALLOCATABLE :: AA00(:), AQ00(:) ! AA00 is for station records for 1D KWE
      REAL(SZ),ALLOCATABLE :: H_DG(:), A_DG(:)
      REAL(SZ),ALLOCATABLE :: H2(:)
C...
C...DECLARE COMMON BLOCKS
C...
      INTEGER NTSTEPS,ITMV
      REAL(SZ) DT,FMV
      REAL(8) TIMEBEG
      COMMON /MEANSQ/ TIMEBEG,DT,FMV,NTSTEPS,ITMV
C
      INTEGER NHARFR
      COMMON /LSQFREQS/ NHARFR
C
      INTEGER NP,NOLICA,NOLIFA,NSCREEN,IHOT,ICS
      COMMON /EXTMODE5/ NP,NOLICA,NOLIFA,NSCREEN,IHOT,ICS

C...
C...DECLARE REAL(8) AND CHAR VARIABLES, EQUIVALENCES
C...
      REAL(8) STATIM,REFTIM,TIME_A,DTDP,TIMEH,vdtdp,cfl_max
      REAL(8) AVGXY,DIF1R,DIF2R,DIF3R
      REAL(8) AEMIN,AE,AA,A1,A2,A3,X1,X2,X3,X4,Y1,Y2,Y3,Y4
      REAL(8) FDX1,FDX2,FDX3,FDY1,FDY2,FDY3
      REAL(8) FDX1OA,FDX2OA,FDX3OA,FDY1OA,FDY2OA,FDY3OA,AREAIE
      REAL(8) DDX1,DDX2,DDX3,DDY1,DDY2,DDY3
      REAL(8) DXX11,DXX12,DXX13,DXX21,DXX22,DXX23,DXX31,DXX32,DXX33
      REAL(8) DYY11,DYY12,DYY13,DYY21,DYY22,DYY23,DYY31,DYY32,DYY33
      REAL(8) DXY11,DXY12,DXY13,DXY21,DXY22,DXY23,DXY31,DXY32,DXY33
      REAL(8) XL0,XL1,XL2,YL0,YL1,YL2,SLAM0,SFEA0
      REAL(8) WREFTIM,WTIMED,WTIME2,WTIME1,WTIMINC,QTIME1,QTIME2
      REAL(8) FTIMINC,ETIMINC,RSTIME1,RSTIME2,RSTIMINC
      REAL(8) DELX,DELY,DIST,DELDIST,DELETA
      REAL(8),ALLOCATABLE :: AMIG(:),AMIGT(:),FAMIG(:)
      REAL(8),ALLOCATABLE :: PER(:),PERT(:),FPER(:)
      REAL(8),ALLOCATABLE :: FREQ(:),FF(:),FACE(:)
      REAL(8),ALLOCATABLE :: SLAM(:),SFEA(:),X(:),Y(:)


      CHARACTER*32 RUNDES
      CHARACTER*24 RUNID,AGRID,AGRID2,AFRIC
      CHARACTER*4  RDES4(8),RID4(6),AID4(6)
      CHARACTER*8  RDES8(4),RID8(3),AID8(3)
      CHARACTER*10 ALPHA
      CHARACTER*5,ALLOCATABLE :: TIPOTAG(:),BOUNTAG(:),FBOUNTAG(:)
      EQUIVALENCE (RDES4(1),RDES8(1),RUNDES), (RID4(1),RID8(1),RUNID),
     &            (AID4(1),AID8(1),AGRID)

C...
C...EXPLICITLY DECLARE ADDITIONAL VARIABLES
C...
      INTEGER  FRW,NODEDRYMIN,NODEWETMIN,IBTYPE,ICK
      INTEGER  IDR,IM,IPRBI,JGW,JKI,JME
      INTEGER  JNMM,KMIN,N3,NABOUT
      INTEGER  NBFR,NBOU,NBVI,NBVJ,NCOR,NE,NE2,NP2
      INTEGER  NEIMIN,NEIMAX,NETA,NFFR,NFLUXB,NFLUXF,NFLUXIB,NFLUXRBC
      INTEGER  NFLUXIBP,NPIPE
      INTEGER  NFOVER,NHG,NHY,NOLICAT,NOPE,NOUTC
      INTEGER  NOUTE,NSPOOLE,NOUTV,NSPOOLV,NPRBI
      INTEGER  NRAMP,NRS,NSTAE,NSTARTDRY,NSTAV,NT,NTCYFE
      INTEGER  NTCYFV,NTCYSE,NTCYSV,NTIF,NTIP,NTRSPE,NTRSPV
      INTEGER  NVEL,NVELEXT,NVELME,NWLAT,NWLON,NWS
      INTEGER  IBSTART, ICHA, ICSTP, IDSETFLG, IE, IER
      INTEGER  IESTP
      INTEGER  IFNLCAT, IFNLCT, IFNLFA
      INTEGER  IFWIND, IGCP, IGEP, IGPP, IGVP, IGWP, IHABEG
#ifdef SWAN
Casey 101118: Added a variable for output of radiation stress gradients.
      INTEGER  IGRadS
#endif
      INTEGER  IHARIND, IHOTSTP, IHSFIL, IJ, ILUMP, IMHS, IPSTP
      INTEGER  IREFYR,IREFMO,IREFDAY,IREFHR,IREFMIN, ISLDIA
      INTEGER  ITIME_A, ITEMPSTP, ITEST, ITHAF, ITHAS
      INTEGER  ITHS, ITITER, ITMAX, IVSTP, IWSTP, IWTIME, IWTIMEP
      INTEGER  IWYR, J12, J13, J21, J23, J31, J32
      INTEGER  JN, KEMAX, KVMAX, LRC, LUMPT, MMAX
      INTEGER  MBW, MDF, MMIN, NA, NBDI, NBDJ, NBNCTOT
      INTEGER  NBW, NC1, NC2, NC3, NCBND
      INTEGER  NCELE, NCI, NCJ, NCTOT, NCYC, NDRY, NDSETSC
      INTEGER  NDSETSE, NDSETSV, NDSETSW, NHAGE, NHAGV
      INTEGER  NHAINC, NHASE
      INTEGER  NHASV, NHSINC, NHSTAR, NM1, NM123, NM2, NM3
      INTEGER  NMI1, NMI2, NMI3, NMJ1, NMJ2, NMJ3, NNBB
      INTEGER  NNBB1, NNBB2, NOUTGC, NOUTGE, NOUTGV, NOUTGW, NOUTM
      INTEGER  NSCOUC, NSCOUE, NSCOUGC, NSCOUGE, NSCOUGV
      INTEGER  NSCOUGW, NSCOUM
      INTEGER  NSCOUV, NSPOOLC, NSPOOLGC, NSPOOLGE
      INTEGER  NSPOOLGV, NSPOOLGW, NSPOOLM
      INTEGER  NSTAC, NSTAM, NTCYFC, NTCYFGC, NTCYFGE
      INTEGER  NTCYFGV, NTCYFGW
      INTEGER  NTCYFM, NTCYSC, NTCYSGC, NTCYSGE
      INTEGER  NTCYSGV, NTCYSGW, NTCYSM
      INTEGER  NTRSPC, NTRSPM, NUMITR, NW, NWET, NWSEGWI, NWSGGWI
      INTEGER  NCCHANGE
      INTEGER  IRAMPING
C
      REAL(SZ) ADVECX, ADVECY, AGIRD, AH, AO12, AO6, ARG
      REAL(SZ) ARG1, ARG2, ARGJ, ARGJ1, ARGJ2, ARGSALT, ARGT
      REAL(SZ) ARGTP, AUV21, AUV22, BARAVGWT
      REAL(SZ) BEDSTR, BNDLEN2O3NC, BSXN1
      REAL(SZ) BSXN2, BSXN3, BSXPP3, BSYN1, BSYN2, BSYN3, BSYPP3
      REAL(SZ) C1, C2, C3, CBEDSTRD, CBEDSTRE, CCRITD, CCSFEA
      REAL(SZ) CELERITY, CH1N1, CH1N2, CH1N3, CHSUM, COND, CONVCR
      REAL(SZ) CORIFPP, DDU, DHDX, DHDY
      REAL(SZ) DISPERX, DISPERY, DT2, DTO2, DTOHPP, DUU1N1, DUU1N2
      REAL(SZ) DUU1N3, DUV1N1, DUV1N2, DUV1N3, DVV1N1, DVV1N2, DVV1N3
      REAL(SZ) DXXYY11, DXXYY12, DXXYY13, DXXYY21
      REAL(SZ) DXXYY22, DXXYY23, DXXYY31
      REAL(SZ) DXXYY32, DXXYY33, DXYH11, DXYH12
      REAL(SZ) DXYH13, DXYH21, DXYH22
      REAL(SZ) DXYH23, DXYH31, DXYH32, DXYH33, E0N1, E0N2, E0N3
      REAL(SZ) E1N1, E1N1SQ, E1N2, E1N2SQ, E1N3, E1N3SQ, ECONST
      REAL(SZ) EE1, EE2, EE3, ELMAX, EP, ESN1, ESN2, BE1, BE2, BE3
      REAL(SZ) ESN3, ETIME1, ETIME2, ETRATIO, EVC1, EVC2, EVC3
      REAL(SZ) EVCEA, EVMPPODT, EVMPPDT, FDDD, FDDDODT, FDDOD, FDDODODT
      REAL(SZ) FIIN, G, GA00, GB00A00, GC00, GDTO2, GFAO2
      REAL(SZ) GHPP, GO3, HABSMIN, HEA, HH1, HH1N1, HH1N2
      REAL(SZ) HH1N3, HH2, HH2N1, HH2N2, HH2N3, HHU1N1, HHU1N2
      REAL(SZ) HHU1N3, HHV1N1, HHV1N2, HHV1N3, HPP, HSD, HSE
      REAL(SZ) HTOT
      REAL(SZ) G2ROOT    !  added for optimization of roe_flux

      REAL(SZ) P11, P22, P33, PR1N1, PR1N2, PR1N3, QFORCEI
      REAL(SZ) QFORCEJ, QTEMA1, QTEMA2, QTEMA3, QTEMB1, QTEMB2, QTEMB3
      REAL(SZ) QTRATIO, QUNORM, QVNORM, RAMP1, RAMP2, RBARWL, RBARWL1
      REAL(SZ) RBARWL1F, RBARWL2, RBARWL2F, RFF, RFF1, RFF2, RHO0
      REAL(SZ) RSTRATIO, RSX, RSY, S2SFEA, SADVDTO3, SALTMUL, SFACPP
      REAL(SZ) SS1N1, SS1N2, SS1N3, T0N1, T0N2, T0N3, T0XN1
      REAL(SZ) T0XN2, T0XN3, T0XPP3, T0YN1, T0YN2, T0YN3, T0YPP3
      REAL(SZ) TADVODT, TAU0AVG, THAF, THAS, THENALLDSSSTUP
      REAL(SZ) TIMEIT, TIPN1,TOUTFC
      REAL(SZ) TIPN2, TIPN3, TKWET, TOUTFGC, TOUTFGE, TOUTFGV, TOUTFGW
      REAL(SZ) TOUTFM, TOUTSGC, TOUTSGE, TOUTSGV, TOUTSGW, TOUTSM, TPMUL
      REAL(SZ) TT0L, TT0R, U11, U1N1, U1N2, U1N3, U22
      REAL(SZ) U33, UEA, UHPP, UHPP3, UN1, UPEA, UPP
      REAL(SZ) UPPDT, UPPDTDDX1, UPPDTDDX2, UPPDTDDX3, UV1, V11, V1N1
      REAL(SZ) R11,R22,R33
C.....Add additional variable declarations for wave friction (EK)
      REAL(SZ) UVW1, UVW2, COSA1, SINA1, WD1, WD, WDXX, WDYY,WDXY,CHYBR
      REAL(SZ) VCOEFXX, VCOEFYY, VCOEFXY, VCOEFYX, VCOEF2
C.....End (EK)
      REAL(SZ) V1N2, V1N3, V22, V33, VCOEF3N1
      REAL(SZ) VCOEF3N2, VCOEF3N3, VCOEF3X, VCOEF3Y, VEA, VEL, VELABS
      REAL(SZ) VELMAX, VELNORM, VELTAN, VHPP, VHPP3, VPEA, VPP
      REAL(SZ) VPPDT, VPPDTDDY1, VPPDTDDY2, VPPDTDDY3
      REAL(SZ) WDRAGCO, WINDMAG, WINDX
      REAL(SZ) WINDY, WS, WSMOD, WSX, WSXN1, WSXN2, WSXN3
      REAL(SZ) WSY, WSYN1, WSYN2, WSYN3, WTRATIO,A00,B00,C00,ANGINN
      REAL(SZ) CORI,COSTHETA,COSTHETA1,COSTSET,CROSS,CROSS1
      REAL(SZ) DAY,DOTVEC,DRAMP,DUM1,DUM2,EVMSUM
      REAL(SZ) H0,H0L,H0H ! H0L and H0H are added by sb 10/27/05
      REAL(SZ) RNDAY,THETA,THETA1,TOUTSC,RAMP,RHOWAT0
      REAL(SZ) TOUTSE,TOUTFE
      REAL(SZ) TOUTSV,TOUTFV,XL
      REAL(SZ) VECNORM,VL1X,VL1Y,VL2X,VL2Y,WLATMAX
      REAL(SZ) WLONMIN,WLATINC,WLONINC
      REAL(SZ) VELMIN


      REAL(8)  RNP_GLOBAL
      REAL(8)  REFSEC   ! required to run in either 32-bit or 64-bit

C...Younghun added follow arrays for open-channels

      REAL(SZ), ALLOCATABLE :: CDEPTHR(:,:),CTWIDTHR(:,:),CBWIDTHR(:,:)
      REAL(SZ), ALLOCATABLE :: CDEPTH(:,:),CTWIDTH(:,:),CBWIDTH(:,:),
     & CSSLOPE(:,:)
      LOGICAL :: METONLY
      

C-------------------end of data declarations----------------------------------C


      CONTAINS

C     ----------------------------------------------------------------
C      F U N C T I O N   S P H E R I C A L   D I S T A N C E
C     ----------------------------------------------------------------
C     jgf49.1001 Function to get the distance along the surface of
C     a sphere (the earth's surface in this case).
C     ----------------------------------------------------------------
      REAL(SZ) FUNCTION sphericalDistance(dx, dy, y1, y2)
      IMPLICIT NONE
      REAL(SZ), intent(in) :: DX    ! longitude distance in radians
      REAL(SZ), intent(in) :: DY    ! latitude distance in radians
      REAL(SZ), intent(in) :: y1    ! degrees latitude of starting point
      REAL(SZ), intent(in) :: y2    ! degrees latitude of ending point
C
C     compute the distances based on haversine formula for
C     distance along a sphere
      sphericalDistance = Rearth *
     &  ( 2.0d0*ASIN(
     &     sqrt(sin(DY/2.0d0)**2.0d0
     &      + cos(y1*DEG2RAD)*cos(y2*DEG2RAD)*sin(DX/2.0d0)**2.0d0
     &   )))

      RETURN
C     ----------------------------------------------------------------
      END FUNCTION sphericalDistance
C     ----------------------------------------------------------------

      SUBROUTINE ALLOC_MAIN1()
C
C     Allocate space for Arrays dimensioned by MNE and MNP
C
      ALLOCATE ( SLAM(MNP),SFEA(MNP),X(MNP),Y(MNP))
      ALLOCATE ( DP(MNP),DP0(MNP),DPe(MNE),SFAC(MNP))!,STARTDRY(MNP))
      ALLOCATE ( NM(MNE,3))
      ALLOCATE ( NOT_AN_EDGE(MNP) )
      ALLOCATE ( WEIR_BUDDY_NODE(MNP,2) )
      ALLOCATE ( ONE_OR_TWO(MNP) )
      ALLOCATE ( NELED(3,MNE))
      ALLOCATE ( ETAS(MNP))
      ALLOCATE ( QW(MNP))
C.....ek now allocated in nodal attributes
C      ALLOCATE ( FRIC(MNP),EVM(MNP))
      ALLOCATE ( UU1(MNP),VV1(MNP))
      ALLOCATE ( NNODECODE(MNP),NODEREP(MNP))
      ALLOCATE ( NNEIGH(MNP),MJU(MNP),NODELE(MNP))
      ALLOCATE ( NIBNODECODE(MNP),NNEIGH_ELEM(MNP))
C      ALLOCATE ( TAU0VAR(MNP))
      ALLOCATE ( CH1(MNP),QB(MNP),QA(MNP),SOURSIN(MNP))!,EVC(MNP))
      ALLOCATE ( TKXX(MNP),TKYY(MNP),TKXY(MNP))
      ALLOCATE ( AREAS(MNE),SFACDUB(3,MNE),RTEMP2(MNE))
      ALLOCATE ( YDUB(36,MNE,8))
      ALLOCATE ( UU2(max(MNP,MNE)))
      ALLOCATE ( VV2(max(MNP,MNE)))
      ALLOCATE ( ETA1(max(MNP,MNE)),ETA2(max(MNP,MNE)), ETAMAX(MNP),
     &  MassMax(MNE) )
      Allocate ( tracer(max(MNP,MNE)), tracer2(max(MNP,MNE)),
     &  bed_int(max(mnp,mne),layers), entrop(5,MNE) )
      ALLOCATE ( CORIF(MNP), dyn_P(max(MNP,MNE)), pdg_el(MNE) )
      ALLOCATE ( QU(MNP),QV(MNP))
      ALLOCATE ( LBCODE(MNP))
      ALLOCATE ( CSI(MNP),SII(MNP))
      ALLOCATE ( NODECODE(MNP))
      ALLOCATE (GRAVX(MNP),GRAVY(MNP))
      ALLOCATE ( NIBCNT(MNP) )    !  added 7/31/2000 to fix wetdry bug
Csb-added by sb 10/11/2005
      ALLOCATE( NNDEL(MNP) )
      ALLOCATE( EDFLG(3,MNE) )
C.....Younghun added for DG-SAKE
      ALLOCATE ( H2(max(MNP,MNE)) )

C--

#ifdef CMPI
      ALLOCATE ( IDUMY(1), DUMY1(1), DUMY2(1) )
#endif

#ifdef CVEC
      ALLOCATE ( QTEMA(MNE,3),QTEMB(MNE,3))
#endif

      IF ( C3DDSS) THEN
        ALLOCATE( AUV11(MNP),AUV12(MNP),AUV13(MNP),AUV14(MNP))
      ELSEIF (C2DDI) THEN
        ALLOCATE( AUV11(MNP),AUV12(MNP))
        ALLOCATE( AUVXX(MNP),AUVYY(MNP),AUVXY(MNP),AUVYX(MNP))
      ENDIF

      IF ( C3D) THEN
       ALLOCATE( DUU1(MNP),DUV1(MNP),DVV1(MNP),BSX1(MNP),BSY1(MNP))
       endif
C
      RETURN
      END SUBROUTINE


      SUBROUTINE ALLOC_MAIN2()
C
C     Allocate space for Arrays dimensioned by MNOPE and MNETA
C
      ALLOCATE ( ESBIN1(MNETA),ESBIN2(MNETA))
      ALLOCATE ( NBDV(MNOPE,MNETA))
      ALLOCATE ( NVDLL(MNOPE),NBD(MNETA))
      ALLOCATE ( EMO(MNBFR,MNETA),EFA(MNBFR,MNETA))
      ALLOCATE ( UMO(MNBFR,MNETA),UFA(MNBFR,MNETA))
      ALLOCATE ( VMO(MNBFR,MNETA),VFA(MNBFR,MNETA))

      RETURN
      END SUBROUTINE


C
C     Allocate space for nonperiodic zero and nonzero normal flow boundary arrays
C     including barriers
C
      SUBROUTINE ALLOC_MAIN3()
      ALLOCATE ( NBV(MNVEL),LBCODEI(MNVEL))
      ALLOCATE ( BNDLEN2O3(MNVEL))
      ALLOCATE ( ME2GW(MNVEL))
      ALLOCATE ( CSII(MNVEL),SIII(MNVEL))
      ALLOCATE ( BARLANHT(MNVEL),BARLANCFSP(MNVEL))
      ALLOCATE ( BARLANHTR(MNVEL),BARLANCFSPR(MNVEL))
      ALLOCATE ( BARINHT(MNVEL),BARINCFSB(MNVEL),BARINCFSP(MNVEL))
      ALLOCATE ( PIPEHT(MNVEL),PIPECOEF(MNVEL),PIPEDIAM(MNVEL))
      ALLOCATE ( IBCONN(MNVEL))
      ALLOCATE ( BARINHTR(MNVEL),BARINCFSBR(MNVEL),BARINCFSPR(MNVEL))
      ALLOCATE ( PIPEHTR(MNVEL),PIPECOEFR(MNVEL),PIPEDIAMR(MNVEL))
      ALLOCATE ( IBCONNR(MNVEL),NTRAN1(MNVEL),NTRAN2(MNVEL))
      ALLOCATE ( BTRAN3(MNVEL),BTRAN4(MNVEL),BTRAN5(MNVEL))
      ALLOCATE ( BTRAN6(MNVEL),BTRAN7(MNVEL),BTRAN8(MNVEL))
      ALLOCATE ( RBARWL1AVG(MNVEL),RBARWL2AVG(MNVEL))
      ALLOCATE ( QNIN1(MNVEL),QNIN2(MNVEL))
      ALLOCATE ( NBVV(MNBOU,0:MNVEL))
      ALLOCATE ( NVELL(MNBOU))
      ALLOCATE ( SEGTYPE(MNBOU))
      
C...Younghun added
      ALLOCATE (CDEPTH(2,MNVEL),CTWIDTH(2,MNVEL),CBWIDTH(2,MNVEL),
     & CSSLOPE(2,MNVEL))
      ALLOCATE (CDEPTHR(MNBOU,MNVEL),CTWIDTHR(MNBOU,MNVEL),CBWIDTHR(MNBOU,MNVEL))
      RETURN
      END SUBROUTINE


C
C     Allocate space for tidal potential terms 
C
      SUBROUTINE ALLOC_MAIN4a()
      ALLOCATE ( TPK(MNTIF),AMIGT(MNTIF),FFT(MNTIF) )
      ALLOCATE ( FACET(MNTIF),PERT(MNTIF),ETRF(MNTIF) )
      ALLOCATE ( TIPOTAG(MNTIF) )

      IF ( CTIP ) THEN
        ALLOCATE( TIP1(MNP),TIP2(MNP))
        ENDIF  

      RETURN
      END SUBROUTINE


C
C     Allocate space for Earth load/self-attraction tide 
C
      SUBROUTINE ALLOC_MAIN4b()
      ALLOCATE ( SALTAMP(MNTIF,MNP),SALTPHA(MNTIF,MNP) )
      RETURN
      END SUBROUTINE


C
C     Allocate space for Arrays dimensioned by MNBFR      
C
      SUBROUTINE ALLOC_MAIN5()
      ALLOCATE ( AMIG(MNBFR),PER(MNBFR))
      ALLOCATE ( FF(MNBFR),FACE(MNBFR))
      ALLOCATE ( BOUNTAG(MNBFR) )

      RETURN
      END SUBROUTINE



C
C     Allocate space for periodic normal flow boundary conditions
C
      SUBROUTINE ALLOC_MAIN6()
      ALLOCATE ( QNAM(MNFFR,MNVEL),QNPH(MNFFR,MNVEL))
      ALLOCATE ( FBOUNTAG(MNFFR))
      ALLOCATE ( FAMIG(MNFFR), FFF(MNFFR), FFACE(MNFFR), FPER(MNFFR) )

      RETURN
      END SUBROUTINE


C
C     Allocate space for arrays used for station elevation output
C
      SUBROUTINE ALLOC_MAIN7()
      ALLOCATE ( NNE(MNSTAE),ET00(MNSTAE),BT00(MNSTAE))
      ALLOCATE ( STAIE1(MNSTAE),STAIE2(MNSTAE),STAIE3(MNSTAE))
      ALLOCATE ( XEL(MNSTAE),YEL(MNSTAE),SLEL(MNSTAE),SFEL(MNSTAE))
      ALLOCATE ( H00(MNSTAE))  ! Younghun added
      ALLOCATE ( AA00(MNSTAE))
      RETURN
      END SUBROUTINE


C
C     Allocate space for arrays used for station velocity output
C
      SUBROUTINE ALLOC_MAIN8()
      ALLOCATE ( XEV(MNSTAV),YEV(MNSTAV),SLEV(MNSTAV),SFEV(MNSTAV))
      ALLOCATE ( NNV(MNSTAV))
      ALLOCATE ( UU00(MNSTAV),VV00(MNSTAV))
      ALLOCATE ( STAIV1(MNSTAV),STAIV2(MNSTAV),STAIV3(MNSTAV))
      ALLOCATE ( HQX00(MNSTAV), HQY00(MNSTAV))
      ALLOCATE ( AQ00(MNSTAV) )
      RETURN
      END SUBROUTINE


C
C     Allocate space for arrays used for station concentration output
C
      SUBROUTINE ALLOC_MAIN9()
      ALLOCATE ( XEC(MNSTAC),YEC(MNSTAC),SLEC(MNSTAC),SFEC(MNSTAC))
      ALLOCATE ( NNC(MNSTAC))
      ALLOCATE ( CC00(MNSTAC))
      ALLOCATE ( STAIC1(MNSTAC),STAIC2(MNSTAC),STAIC3(MNSTAC))

      RETURN
      END SUBROUTINE


C
C     Allocate space for arrays used for station meteorological output
C
      SUBROUTINE ALLOC_MAIN10()
      ALLOCATE ( XEM(MNSTAM),YEM(MNSTAM),SLEM(MNSTAM),SFEM(MNSTAM))
      ALLOCATE ( NNM(MNSTAM))
      ALLOCATE ( RMU00(MNSTAM),RMV00(MNSTAM),RMP00(MNSTAM),RMR00(MNSTAM))
      ALLOCATE ( STAIM1(MNSTAM),STAIM2(MNSTAM),STAIM3(MNSTAM))

      RETURN
      END SUBROUTINE


C
C     Allocate space for Arrays dimensioned by MNEI   
C
      SUBROUTINE ALLOC_MAIN11()

C  Arrays used by JCG iterative solver

      ALLOCATE( OBCCOEF(MNETA,MNEI-1),COEF(MNP,MNEI))
      ALLOCATE( IWKSP(3*MNP),WKSP(4*MNP+400) )
      ALLOCATE( IPARM(12),RPARM(12) )

C  Neighbor Table

      ALLOCATE ( NEITAB(MNP,MNEI),NEIGH_ELEM(MNP,MNEI))

      RETURN
      END SUBROUTINE


C
C     Allocate space for wind forcing   
C
      SUBROUTINE ALLOC_MAIN12()
      ALLOCATE ( WSX1(MNP),WSY1(MNP),PR1(MNP) )
      ALLOCATE ( WSX2(MNP),WSY2(MNP),PR2(MNP) )
      ALLOCATE ( WVNX1(MNP),WVNY1(MNP),PRN1(MNP) )
      ALLOCATE ( WVNX2(MNP),WVNY2(MNP),PRN2(MNP) )
      ALLOCATE ( RSNX1(MNP),RSNY1(MNP),RSNX2(MNP),RSNY2(MNP) )
      ALLOCATE ( WVNXOUT(MNP),WVNYOUT(MNP) )
#ifdef SWAN
Casey 101118: Added the next line for output of radiation stress gradients.
      ALLOCATE ( RSNXOUT(MNP),RSNYOUT(MNP) )
#endif
      RETURN
      END SUBROUTINE


C
C     Allocate space for bridge piling friction arrays   
C
      SUBROUTINE ALLOC_MAIN13()
      ALLOCATE ( NBNNUM(MNP),BK(MNP),BALPHA(MNP),BDELX(MNP) )
      RETURN
      END SUBROUTINE
C
C     Allocate space for harmonic analysis means and variance calculations, this is in 
C     global data because the variables are used in main source, outside of HA analysis
C     subroutines.  This should probably be changed.
C
      SUBROUTINE ALLOC_MAIN14()
      ALLOCATE ( XVELAV(MNP),YVELAV(MNP),XVELVA(MNP),YVELVA(MNP) )
      ALLOCATE ( ELAV(MNP),ELVA(MNP) )
      RETURN
      END SUBROUTINE
      
C.....Allocate for wave friction

      SUBROUTINE ALLOC_MAIN15()
      ALLOCATE ( WAVE_T1(MNP),WAVE_H1(MNP),WAVE_A1(MNP),WAVE_D1(MNP) )
      ALLOCATE ( WAVE_T2(MNP),WAVE_H2(MNP),WAVE_A2(MNP),WAVE_D2(MNP) )
      ALLOCATE ( WAVE_T(MNP),WAVE_H(MNP),WAVE_A(MNP),WAVE_D(MNP) )
      ALLOCATE ( WB(MNP) )
      END SUBROUTINE

C.....Allocate arrays used for node to element table
C
      SUBROUTINE ALLOC_NNOEL1(MNP)
      ALLOCATE ( EL_COUNT(MNP) )
      RETURN
      END SUBROUTINE
C
      SUBROUTINE ALLOC_NNOEL2(MNP,MAXEL)
      ALLOCATE ( DP_DG(MAXEL),DG_ANG(MAXEL) )
      ALLOCATE ( NNOEL(MNP,MAXEL),CENTAB(MNP,MAXEL+1) )
      ALLOCATE (ELETAB(MNP,MAXEL+1),ANGTAB(MNP,MAXEL+1),dynP_DG(MAXEL) )
      ALLOCATE ( ZE_DG(MAXEL), QX_DG(MAXEL), QY_DG(MAXEL), HB_DG(MAXEL))
      Allocate ( iota2_DG(MAXEL), iota_DG(MAXEL), iotaa_DG(MAXEL) )
      Allocate(bed_DG(MAXEL,layers),bed_N_int(layers),bed_N_ext(layers))
      ALLOCATE ( H_DG(MAXEL), A_DG(MAXEL) )
      RETURN
      END SUBROUTINE

Csb
C     Allocate space for Arrays dimensioned by MNNDEL
C
      SUBROUTINE ALLOC_MAIN16()

C.....Node-to-elements table      
      ALLOCATE( NDEL(MNP,MNNDEL) )

      RETURN
      END SUBROUTINE


      END MODULE GLOBAL
