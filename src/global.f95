!**************************************************************************
!  mod history
!  v41.06mxxx  - date - programmer - describe change 
!                    - mark change in code with  cinitials-mxxx
!  v41.10      - 07/25/01 - rl - from 41.09 - bug fix in GWCE lateral viscosity term
!  v41.09      - 06/30/01 - jw - from 41.08 - made minor mods as per vp version 41.05
!  v41.06      - 04/02/01 - rl - changed MNWP to MNP in wind forcing 
!                                ALLOCATION statements
!  v41.02      - 09/04 - rl
!  v40.02m004b - 05/17 - vjp - corrected dimensioning problem cvjpm004b
!  v40.02m001  - 12/21 - jjw - add cross barrier pipes cjjwm001
!-----------------------------------------------------------------------
!
!     v9_sb2.2.1 - 08/16/05 - sb - Hard bottom implementation
!     v10_sb5    - 10/11/05 - sb - Table that associates nodes to elements (NODEELEM)
!                - 10/27/05 - sb - H0L and H0H are added
!
!**************************************************************************
! 
      MODULE GLOBAL
      USE SIZES
!...
!...SET GLOBAL PARAMETER CONSTANTS
!...

#ifdef SWAN
!asey 121128: Added variables for the output of global files from SWAN.
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

!.....PI and degrees to radians conversions
      REAL(8), PARAMETER  ::  PI=3.141592653589793D0
      REAL(8), PARAMETER  ::  DEG2RAD = PI/180.D0
      REAL(8), PARAMETER  ::  RAD2DEG = 180.D0/PI

!.....parameters used in barrier overflow 
      REAL(SZ), PARAMETER ::  BARMIN=0.01D0
      REAL(SZ) DEPAVG,DEPMAX,DEPMIN
!
!.....Sediment Transport stuff added (Ethan Kubatko 8-1-2003)
!
!.....Declare variables used in sediment transport section
!
      INTEGER, TARGET :: SEDFLAG
      INTEGER MAXEL, ELEM_ED, NBOR_ED, NBOR_EL,ITDG
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
!     jgf46.08 Fine grained ramp functions (jgf46.21 split flux b.c.s)
      REAL(SZ) RampExtFlux,DRampExtFlux ! Ramp for external flux b.c.s
      REAL(SZ) RampIntFlux,DRampIntFlux ! Ramp for internal flux b.c.s
      REAL(SZ) RampElev,DRampElev    ! Ramp for elevation boundary conditions.
      REAL(SZ) RampTip,DRampTip      ! Ramp for tidal potential
      REAL(SZ) RampMete,DRampMete    ! Ramp for wind and atmospheric pressure
      REAL(SZ) RampWRad,DRampWRad    ! Ramp for wave radiation stress

!.....Create allocatable arrays

      REAL(SZ),ALLOCATABLE :: ANGTAB(:,:),CENTAB(:,:),ELETAB(:,:)
      REAL(SZ),ALLOCATABLE :: DG_ANG(:),DP_DG(:)
      INTEGER, ALLOCATABLE :: EL_COUNT(:)
      INTEGER, ALLOCATABLE :: NNOEL(:,:)

!sb...Table that maps nodes to elements
      INTEGER, ALLOCATABLE :: NNDEL(:)
      INTEGER, ALLOCATABLE :: NDEL(:,:)
!--
      
!.....Declare variables for DG SW

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
!...
!...DECLARE ALL ARRAYS
!...
#ifdef CMPI
      INTEGER, ALLOCATABLE ::   IDUMY(:)
      REAL(SZ),ALLOCATABLE ::   DUMY1(:),DUMY2(:)
      REAL(SZ),ALLOCATABLE ::   DGDUMY1(:,:,:),DGDUMY2(:,:,:)
#endif

      INTEGER, ALLOCATABLE ::   pdg_el(:)
      REAL(SZ),ALLOCATABLE ::   ETAS(:),ETA1(:),ETA2(:),ETAMAX(:),entrop(:,:)
      Real(SZ),Allocatable ::   tracer(:),tracer2(:),MassMax(:),bed_int(:,:)
      REAL(SZ),ALLOCATABLE ::   UU1(:),UU2(:),VV1(:),VV2(:), dyn_P(:)
      REAL(SZ),ALLOCATABLE ::   DP(:),DP0(:),DPe(:),SFAC(:) !,STARTDRY(:)
      REAL(SZ),ALLOCATABLE ::   QU(:),QV(:),QW(:)
      REAL(SZ),ALLOCATABLE ::   CORIF(:) !,EVM(:)
      REAL(SZ),ALLOCATABLE ::   TPK(:),FFT(:)
      REAL(SZ),ALLOCATABLE ::   FACET(:),ETRF(:)
      REAL(SZ),ALLOCATABLE ::   ESBIN1(:),ESBIN2(:)
      REAL(SZ),ALLOCATABLE ::   QTEMA(:,:),QTEMB(:,:)
      REAL(SZ),ALLOCATABLE ::   QN0(:),QN1(:),QN2(:)
      REAL(SZ),ALLOCATABLE ::   BNDLEN2O3(:)
      REAL(SZ),ALLOCATABLE ::   CSII(:),SIII(:)
      REAL(SZ),ALLOCATABLE ::   QNAM(:,:),QNPH(:,:)
      REAL(SZ),ALLOCATABLE ::   QNIN1(:),QNIN2(:)
      REAL(SZ),ALLOCATABLE ::   CSI(:),SII(:)
!      REAL(SZ),ALLOCATABLE ::   TAU0VAR(:)
      REAL(SZ),ALLOCATABLE ::   ET00(:),BT00(:)
      REAL(SZ),ALLOCATABLE ::   STAIE1(:),STAIE2(:),STAIE3(:)
      REAL(8),ALLOCATABLE ::    XEV(:),YEV(:),SLEV(:),SFEV(:)
      REAL(SZ),ALLOCATABLE ::   UU00(:),VV00(:)
      REAL(SZ),ALLOCATABLE ::   STAIV1(:),STAIV2(:),STAIV3(:)
      REAL(8),ALLOCATABLE ::    XEC(:),YEC(:),SLEC(:),SFEC(:)
      REAL(SZ),ALLOCATABLE ::   CC00(:)
      REAL(SZ),ALLOCATABLE ::   STAIC1(:),STAIC2(:),STAIC3(:)
      REAL(8),ALLOCATABLE ::    XEM(:),YEM(:),SLEM(:),SFEM(:)
      REAL(SZ),ALLOCATABLE ::   RMU00(:),RMV00(:),RMP00(:)
      REAL(SZ),ALLOCATABLE ::   STAIM1(:),STAIM2(:),STAIM3(:)
      REAL(SZ),ALLOCATABLE ::   CH1(:),QB(:),QA(:),SOURSIN(:)!,EVC(:)
      REAL(SZ),ALLOCATABLE ::   WSX1(:),WSY1(:),PR1(:)
      REAL(SZ),ALLOCATABLE ::   WSX2(:),WSY2(:),PR2(:)
      REAL(SZ),ALLOCATABLE ::   WVNX1(:),WVNY1(:),PRN1(:)
      REAL(SZ),ALLOCATABLE ::   WVNX2(:),WVNY2(:),PRN2(:)
      REAL(SZ),ALLOCATABLE ::   RSNX1(:),RSNY1(:),RSNX2(:),RSNY2(:)
#ifdef SWAN
!asey 101118: Added the following arrays for output of radiation stress gradients.
      REAL(SZ),ALLOCATABLE,TARGET :: RSNXOUT(:), RSNYOUT(:)
#endif
!.....Allocate additional variables for wave friction
      REAL(SZ),ALLOCATABLE ::   WAVE_T1(:),WAVE_H1(:),WAVE_A1(:),WAVE_D1(:)
      REAL(SZ),ALLOCATABLE ::   WAVE_T2(:),WAVE_H2(:),WAVE_A2(:),WAVE_D2(:)
      REAL(SZ),ALLOCATABLE ::   WAVE_T(:), WAVE_H(:), WAVE_A(:), WAVE_D(:)
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

!.....for internal barrier boundaries with flowthrough pipes
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
      REAL(SZ),ALLOCATABLE ::   RPIPEWL1AVG(:),RPIPEWL2AVG(:)
      REAL(SZ),ALLOCATABLE ::   ELEXLEN(:,:)
      INTEGER, ALLOCATABLE ::   IBCONN(:),IBCONNR(:),NTRAN1(:),NTRAN2(:)

!.....for bridge pilings
      REAL(SZ)                  POAN,Fr,FRICBP
      INTEGER                   NBPNODES
      REAL(SZ),ALLOCATABLE ::   BK(:),BALPHA(:),BDELX(:)
      INTEGER, ALLOCATABLE ::   NBNNUM(:)

!...
!...DECLARE COMMON BLOCKS
!...
      INTEGER NTSTEPS,ITMV
      REAL(SZ) DT,FMV
      REAL(8) TIMEBEG
      COMMON /MEANSQ/ TIMEBEG,DT,FMV,NTSTEPS,ITMV
!
      INTEGER NHARFR
      COMMON /LSQFREQS/ NHARFR
!
      INTEGER NP,NOLICA,NOLIFA,NSCREEN,IHOT,ICS
      COMMON /EXTMODE5/ NP,NOLICA,NOLIFA,NSCREEN,IHOT,ICS

!...
!...DECLARE REAL(8) AND CHAR VARIABLES, EQUIVALENCES
!...
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
      EQUIVALENCE (RDES4(1),RDES8(1),RUNDES), (RID4(1),RID8(1),RUNID),&
                 (AID4(1),AID8(1),AGRID)

!...
!...EXPLICITLY DECLARE ADDITIONAL VARIABLES
!...
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
!asey 101118: Added a variable for output of radiation stress gradients.
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
!
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
!.....Add additional variable declarations for wave friction (EK)
      REAL(SZ) UVW1, UVW2, COSA1, SINA1, WD1, WD, WDXX, WDYY,WDXY,CHYBR
      REAL(SZ) VCOEFXX, VCOEFYY, VCOEFXY, VCOEFYX, VCOEF2
!.....End (EK)
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
      


!-------------------end of data declarations----------------------------------C


      CONTAINS

      SUBROUTINE ALLOC_MAIN1(s)
      use sizes
      implicit none
      type (sizes_type) :: s
!
!     Allocate space for Arrays dimensioned by MNE and MNP
!
      ALLOCATE ( SLAM(s%MNP),SFEA(s%MNP),X(s%MNP),Y(s%MNP))
      ALLOCATE ( DP(s%MNP),DP0(s%MNP),DPe(s%MNE),SFAC(s%MNP))!,STARTDRY(s%MNP))
      ALLOCATE ( NM(s%MNE,3))
      ALLOCATE ( NOT_AN_EDGE(s%MNP) )
      ALLOCATE ( WEIR_BUDDY_NODE(s%MNP,2) )
      ALLOCATE ( ONE_OR_TWO(s%MNP) )
      ALLOCATE ( NELED(3,s%MNE))
      ALLOCATE ( ETAS(s%MNP))
      ALLOCATE ( QW(s%MNP))
!.....ek now allocated in nodal attributes
!      ALLOCATE ( FRIC(s%MNP),EVM(s%MNP))
      ALLOCATE ( UU1(s%MNP),VV1(s%MNP))
      ALLOCATE ( NNODECODE(s%MNP),NODEREP(s%MNP))
      ALLOCATE ( NNEIGH(s%MNP),MJU(s%MNP),NODELE(s%MNP))
      ALLOCATE ( NIBNODECODE(s%MNP),NNEIGH_ELEM(s%MNP))
!      ALLOCATE ( TAU0VAR(s%MNP))
      ALLOCATE ( CH1(s%MNP),QB(s%MNP),QA(s%MNP),SOURSIN(s%MNP))!,EVC(s%MNP))
      ALLOCATE ( TKXX(s%MNP),TKYY(s%MNP),TKXY(s%MNP))
      ALLOCATE ( AREAS(s%MNE),SFACDUB(3,s%MNE),RTEMP2(s%MNE))
      ALLOCATE ( YDUB(36,s%MNE,8))
      ALLOCATE ( UU2(max(s%MNP,s%MNE)))
      ALLOCATE ( VV2(max(s%MNP,s%MNE)))
      ALLOCATE ( ETA1(max(s%MNP,s%MNE)),ETA2(max(s%MNP,s%MNE)), ETAMAX(s%MNP), MassMax(s%MNE) )
      Allocate ( tracer(max(s%MNP,s%MNE)), tracer2(max(s%MNP,s%MNE)), bed_int(max(s%mnp,s%mne),s%layers), entrop(5,s%MNE) )
      ALLOCATE ( CORIF(s%MNP), dyn_P(max(s%MNP,s%MNE)), pdg_el(s%MNE) )
      ALLOCATE ( QU(s%MNP),QV(s%MNP))
      ALLOCATE ( LBCODE(s%MNP))
      ALLOCATE ( CSI(s%MNP),SII(s%MNP))
      ALLOCATE ( NODECODE(s%MNP))
      ALLOCATE (GRAVX(s%MNP),GRAVY(s%MNP))
      ALLOCATE ( NIBCNT(s%MNP) )    !  added 7/31/2000 to fix wetdry bug
!sb-added by sb 10/11/2005
      ALLOCATE( NNDEL(s%MNP) )
      ALLOCATE( EDFLG(3,s%MNE) )

!--

#ifdef CMPI
      ALLOCATE ( IDUMY(1), DUMY1(1), DUMY2(1) )
#endif

#ifdef CVEC
      ALLOCATE ( QTEMA(s%MNE,3),QTEMB(s%MNE,3))
#endif

      IF ( s%C3DDSS) THEN
        ALLOCATE( AUV11(s%MNP),AUV12(s%MNP),AUV13(s%MNP),AUV14(s%MNP))
      ELSEIF (s%C2DDI) THEN
        ALLOCATE( AUV11(s%MNP),AUV12(s%MNP))
        ALLOCATE( AUVXX(s%MNP),AUVYY(s%MNP),AUVXY(s%MNP),AUVYX(s%MNP))
      ENDIF

      IF ( s%C3D) THEN
       ALLOCATE( DUU1(s%MNP),DUV1(s%MNP),DVV1(s%MNP),BSX1(s%MNP),BSY1(s%MNP))
       endif
!
      RETURN
      END SUBROUTINE


      SUBROUTINE ALLOC_MAIN2(s)
      use sizes
      implicit none
      type (sizes_type) :: s
!
!     Allocate space for Arrays dimensioned by MNOPE and MNETA
!
      ALLOCATE ( ESBIN1(s%MNETA),ESBIN2(s%MNETA))
      ALLOCATE ( NBDV(s%MNOPE,s%MNETA))
      ALLOCATE ( NVDLL(s%MNOPE),NBD(s%MNETA))
      ALLOCATE ( EMO(s%MNBFR,s%MNETA),EFA(s%MNBFR,s%MNETA))
      ALLOCATE ( UMO(s%MNBFR,s%MNETA),UFA(s%MNBFR,s%MNETA))
      ALLOCATE ( VMO(s%MNBFR,s%MNETA),VFA(s%MNBFR,s%MNETA))

      RETURN
      END SUBROUTINE


!
!     Allocate space for nonperiodic zero and nonzero normal flow boundary arrays
!     including barriers
!
      SUBROUTINE ALLOC_MAIN3(s)
      use sizes
      implicit none
      type (sizes_type) :: s
      ALLOCATE ( QN0(s%MNVEL),QN1(s%MNVEL),QN2(s%MNVEL))
      ALLOCATE ( NBV(s%MNVEL),LBCODEI(s%MNVEL))
      ALLOCATE ( BNDLEN2O3(s%MNVEL))
      ALLOCATE ( ME2GW(s%MNVEL))
      ALLOCATE ( CSII(s%MNVEL),SIII(s%MNVEL))
      ALLOCATE ( BARLANHT(s%MNVEL),BARLANCFSP(s%MNVEL))
      ALLOCATE ( BARLANHTR(s%MNVEL),BARLANCFSPR(s%MNVEL))
      ALLOCATE ( BARINHT(s%MNVEL),BARINCFSB(s%MNVEL),BARINCFSP(s%MNVEL))
      ALLOCATE ( PIPEHT(s%MNVEL),PIPECOEF(s%MNVEL),PIPEDIAM(s%MNVEL))
      ALLOCATE ( IBCONN(s%MNVEL))
      ALLOCATE ( BARINHTR(s%MNVEL),BARINCFSBR(s%MNVEL),BARINCFSPR(s%MNVEL))
      ALLOCATE ( PIPEHTR(s%MNVEL),PIPECOEFR(s%MNVEL),PIPEDIAMR(s%MNVEL))
      ALLOCATE ( IBCONNR(s%MNVEL),NTRAN1(s%MNVEL),NTRAN2(s%MNVEL))
      ALLOCATE ( BTRAN3(s%MNVEL),BTRAN4(s%MNVEL),BTRAN5(s%MNVEL))
      ALLOCATE ( BTRAN6(s%MNVEL),BTRAN7(s%MNVEL),BTRAN8(s%MNVEL))
      ALLOCATE ( RBARWL1AVG(s%MNVEL),RBARWL2AVG(s%MNVEL))
      ALLOCATE ( RPIPEWL1AVG(s%MNVEL),RPIPEWL2AVG(s%MNVEL))
      ALLOCATE ( QNIN1(s%MNVEL),QNIN2(s%MNVEL))
      ALLOCATE ( NBVV(s%MNBOU,0:s%MNVEL))
      ALLOCATE ( NVELL(s%MNBOU))
      ALLOCATE ( SEGTYPE(s%MNBOU))

      RETURN
      END SUBROUTINE


!
!     Allocate space for tidal potential terms 
!
      SUBROUTINE ALLOC_MAIN4a(s)
      use sizes
      implicit none
      type (sizes_type) :: s
      ALLOCATE ( TPK(s%MNTIF),AMIGT(s%MNTIF),FFT(s%MNTIF) )
      ALLOCATE ( FACET(s%MNTIF),PERT(s%MNTIF),ETRF(s%MNTIF) )
      ALLOCATE ( TIPOTAG(s%MNTIF) )

      IF ( s%CTIP ) THEN
        ALLOCATE( TIP1(s%MNP),TIP2(s%MNP))
        ENDIF  

      RETURN
      END SUBROUTINE


!
!     Allocate space for Earth load/self-attraction tide 
!
      SUBROUTINE ALLOC_MAIN4b(s)
      use sizes
      implicit none
      type (sizes_type) :: s
      ALLOCATE ( SALTAMP(s%MNTIF,s%MNP),SALTPHA(s%MNTIF,s%MNP) )
      RETURN
      END SUBROUTINE


!
!     Allocate space for Arrays dimensioned by s%MNBFR      
!
      SUBROUTINE ALLOC_MAIN5(s)
      use sizes
      implicit none
      type (sizes_type) :: s
      ALLOCATE ( AMIG(s%MNBFR),PER(s%MNBFR))
      ALLOCATE ( FF(s%MNBFR),FACE(s%MNBFR))
      ALLOCATE ( BOUNTAG(s%MNBFR) )

      RETURN
      END SUBROUTINE



!
!     Allocate space for periodic normal flow boundary conditions
!
      SUBROUTINE ALLOC_MAIN6(s)
      use sizes
      implicit none
      type (sizes_type) :: s
      ALLOCATE ( QNAM(s%MNFFR,s%MNVEL),QNPH(s%MNFFR,s%MNVEL))
      ALLOCATE ( FBOUNTAG(s%MNFFR))
      ALLOCATE ( FAMIG(s%MNFFR), FFF(s%MNFFR), FFACE(s%MNFFR), FPER(s%MNFFR) )

      RETURN
      END SUBROUTINE


!
!     Allocate space for arrays used for station elevation output
!
      SUBROUTINE ALLOC_MAIN7(s)
      use sizes
      implicit none
      type (sizes_type) :: s
      ALLOCATE ( NNE(s%MNSTAE),ET00(s%MNSTAE),BT00(s%MNSTAE))
      ALLOCATE ( STAIE1(s%MNSTAE),STAIE2(s%MNSTAE),STAIE3(s%MNSTAE))
      ALLOCATE ( XEL(s%MNSTAE),YEL(s%MNSTAE),SLEL(s%MNSTAE),SFEL(s%MNSTAE))
      RETURN
      END SUBROUTINE


!
!     Allocate space for arrays used for station velocity output
!
      SUBROUTINE ALLOC_MAIN8(s)
      use sizes
      implicit none
      type (sizes_type) :: s
      ALLOCATE ( XEV(s%MNSTAV),YEV(s%MNSTAV),SLEV(s%MNSTAV),SFEV(s%MNSTAV))
      ALLOCATE ( NNV(s%MNSTAV))
      ALLOCATE ( UU00(s%MNSTAV),VV00(s%MNSTAV))
      ALLOCATE ( STAIV1(s%MNSTAV),STAIV2(s%MNSTAV),STAIV3(s%MNSTAV))

      RETURN
      END SUBROUTINE


!
!     Allocate space for arrays used for station concentration output
!
      SUBROUTINE ALLOC_MAIN9(s)
      use sizes
      implicit none
      type (sizes_type) :: s
      ALLOCATE ( XEC(s%MNSTAC),YEC(s%MNSTAC),SLEC(s%MNSTAC),SFEC(s%MNSTAC))
      ALLOCATE ( NNC(s%MNSTAC))
      ALLOCATE ( CC00(s%MNSTAC))
      ALLOCATE ( STAIC1(s%MNSTAC),STAIC2(s%MNSTAC),STAIC3(s%MNSTAC))

      RETURN
      END SUBROUTINE


!
!     Allocate space for arrays used for station meteorological output
!
      SUBROUTINE ALLOC_MAIN10(s)
      use sizes
      implicit none
      type (sizes_type) :: s
      ALLOCATE ( XEM(s%MNSTAM),YEM(s%MNSTAM),SLEM(s%MNSTAM),SFEM(s%MNSTAM))
      ALLOCATE ( NNM(s%MNSTAM))
      ALLOCATE ( RMU00(s%MNSTAM),RMV00(s%MNSTAM),RMP00(s%MNSTAM))
      ALLOCATE ( STAIM1(s%MNSTAM),STAIM2(s%MNSTAM),STAIM3(s%MNSTAM))

      RETURN
      END SUBROUTINE


!
!     Allocate space for Arrays dimensioned by s%s%MNEI   
!
      SUBROUTINE ALLOC_MAIN11(s)
      use sizes
      implicit none
      type (sizes_type) :: s

!  Arrays used by JCG iterative solver

      ALLOCATE( OBCCOEF(s%MNETA,s%MNEI-1),COEF(s%MNP,s%MNEI))
      ALLOCATE( IWKSP(3*s%MNP),WKSP(4*s%MNP+400) )
      ALLOCATE( IPARM(12),RPARM(12) )

!  Neighbor Table

      ALLOCATE ( NEITAB(s%MNP,s%MNEI),NEIGH_ELEM(s%MNP,s%MNEI))

      RETURN
      END SUBROUTINE


!
!     Allocate space for wind forcing   
!
      SUBROUTINE ALLOC_MAIN12(s)
      use sizes
      implicit none
      type (sizes_type) :: s
      ALLOCATE ( WSX1(s%MNP),WSY1(s%MNP),PR1(s%MNP) )
      ALLOCATE ( WSX2(s%MNP),WSY2(s%MNP),PR2(s%MNP) )
      ALLOCATE ( WVNX1(s%MNP),WVNY1(s%MNP),PRN1(s%MNP) )
      ALLOCATE ( WVNX2(s%MNP),WVNY2(s%MNP),PRN2(s%MNP) )
      ALLOCATE ( RSNX1(s%MNP),RSNY1(s%MNP),RSNX2(s%MNP),RSNY2(s%MNP) )
      ALLOCATE ( WVNXOUT(s%MNP),WVNYOUT(s%MNP) )
#ifdef SWAN
!asey 101118: Added the next line for output of radiation stress gradients.
      ALLOCATE ( RSNXOUT(s%MNP),RSNYOUT(s%MNP) )
#endif
      RETURN
      END SUBROUTINE


!
!     Allocate space for bridge piling friction arrays   
!
      SUBROUTINE ALLOC_MAIN13(s)
      use sizes
      implicit none
      type (sizes_type) :: s
      ALLOCATE ( NBNNUM(s%MNP),BK(s%MNP),BALPHA(s%MNP),BDELX(s%MNP) )
      RETURN
      END SUBROUTINE
!
!     Allocate space for harmonic analysis means and variance calculations, this is in 
!     global data because the variables are used in main source, outside of HA analysis
!     subroutines.  This should probably be changed.
!
      SUBROUTINE ALLOC_MAIN14(s)
      use sizes
      implicit none
      type (sizes_type) :: s
      ALLOCATE ( XVELAV(s%MNP),YVELAV(s%MNP),XVELVA(s%MNP),YVELVA(s%MNP) )
      ALLOCATE ( ELAV(s%MNP),ELVA(s%MNP) )
      RETURN
      END SUBROUTINE
      
!.....Allocate for wave friction

      SUBROUTINE ALLOC_MAIN15(s)
      use sizes
      implicit none
      type (sizes_type) :: s
      ALLOCATE ( WAVE_T1(s%MNP),WAVE_H1(s%MNP),WAVE_A1(s%MNP),WAVE_D1(s%MNP) )
      ALLOCATE ( WAVE_T2(s%MNP),WAVE_H2(s%MNP),WAVE_A2(s%MNP),WAVE_D2(s%MNP) )
      ALLOCATE ( WAVE_T(s%MNP),WAVE_H(s%MNP),WAVE_A(s%MNP),WAVE_D(s%MNP) )
      ALLOCATE ( WB(s%MNP) )
      END SUBROUTINE

!.....Allocate arrays used for node to element table
!
      SUBROUTINE ALLOC_NNOEL1(s)
      use sizes
      implicit none
      type (sizes_type) :: s
      ALLOCATE ( EL_COUNT(s%MNP) )
      RETURN
      END SUBROUTINE
!
      SUBROUTINE ALLOC_NNOEL2(s,MAXEL)
      use sizes
      implicit none
      type (sizes_type) :: s
      integer :: MAXEL
      ALLOCATE ( DP_DG(MAXEL),DG_ANG(MAXEL) )
      ALLOCATE ( NNOEL(s%MNP,MAXEL),CENTAB(s%MNP,MAXEL+1) )
      ALLOCATE ( ELETAB(s%MNP,MAXEL+1),ANGTAB(s%MNP,MAXEL+1),dynP_DG(MAXEL) )
      ALLOCATE ( ZE_DG(MAXEL), QX_DG(MAXEL), QY_DG(MAXEL), HB_DG(MAXEL))
      Allocate ( iota2_DG(MAXEL), iota_DG(MAXEL), iotaa_DG(MAXEL) )
      Allocate ( bed_DG(MAXEL,s%layers),bed_N_int(s%layers), bed_N_ext(s%layers) )
      RETURN
      END SUBROUTINE

!sb
!     Allocate space for Arrays dimensioned by   MNNDEL
!
      SUBROUTINE ALLOC_MAIN16(s)
      use sizes
      implicit none
      type (sizes_type) :: s

!.....Node-to-elements table      
      ALLOCATE( NDEL(s%MNP,s%MNNDEL) )

      RETURN
      END SUBROUTINE


      END MODULE GLOBAL
