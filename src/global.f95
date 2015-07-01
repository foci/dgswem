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

!.....PI and degrees to radians conversions
      REAL(8), PARAMETER  ::  PI=3.141592653589793D0
      REAL(8), PARAMETER  ::  DEG2RAD = PI/180.D0
      REAL(8), PARAMETER  ::  RAD2DEG = 180.D0/PI

!.....parameters used in barrier overflow 
      REAL(SZ), PARAMETER ::  BARMIN=0.01D0
      REAL(SZ) DEPAVG,DEPMAX,DEPMIN

      type global_type

      Logical vertexslope
!
!.....Sediment Transport stuff added (Ethan Kubatko 8-1-2003)
!
!.....Declare variables used in sediment transport section
!
      INTEGER :: SEDFLAG
      INTEGER MAXEL, ELEM_ED, NBOR_ED, NBOR_EL,ITDG
      Integer tracer_flag, chem_flag
      INTEGER N1,N2,NO_NBORS,NBOR,SEDFLAG_W, OPEN_INDEX
      INTEGER :: DG_TO_CG
      INTEGER NSCREEN_INC
      
      INTEGER ScreenUnit
      
      REAL(SZ) AREA_SUM, CEN_SUM, NLEQ, LEQ, NLEQG

      REAL(SZ) :: reaction_rate
      Character*100 :: sed_equationX,sed_equationY
      
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

      INTEGER :: DGSWE
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
!      COMMON /MEANSQ/ TIMEBEG,DT,FMV,NTSTEPS,ITMV
!
      INTEGER NHARFR
!      COMMON /LSQFREQS/ NHARFR
!
      INTEGER NP,NOLICA,NOLIFA,NSCREEN,IHOT,ICS
!      COMMON /EXTMODE5/ NP,NOLICA,NOLIFA,NSCREEN,IHOT,ICS

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
!      EQUIVALENCE (RDES4(1),RDES8(1),RUNDES), (RID4(1),RID8(1),RUNID),&
!                 (AID4(1),AID8(1),AGRID)

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
      
   end type global_type

!-------------------end of data declarations----------------------------------C


      CONTAINS

      SUBROUTINE ALLOC_MAIN1(s,global_here)
      use sizes
      implicit none
      type (sizes_type) :: s
      type (global_type) :: global_here
      
!
!     Allocate space for Arrays dimensioned by MNE and MNP
!
      ALLOCATE ( global_here%SLAM(s%MNP),global_here%SFEA(s%MNP),global_here%X(s%MNP),global_here%Y(s%MNP))
      ALLOCATE ( global_here%DP(s%MNP),global_here%DP0(s%MNP),global_here%DPe(s%MNE),global_here%SFAC(s%MNP))!,STARTDRY(s%MNP))
      ALLOCATE ( global_here%NM(s%MNE,3))
      ALLOCATE ( global_here%NOT_AN_EDGE(s%MNP) )
      ALLOCATE ( global_here%WEIR_BUDDY_NODE(s%MNP,2) )
      ALLOCATE ( global_here%ONE_OR_TWO(s%MNP) )
      ALLOCATE ( global_here%NELED(3,s%MNE))
      ALLOCATE ( global_here%ETAS(s%MNP))
      ALLOCATE ( global_here%QW(s%MNP))
!.....ek now allocated in nodal attributes
!      ALLOCATE ( FRIC(s%MNP),EVM(s%MNP))
      ALLOCATE ( global_here%UU1(s%MNP),global_here%VV1(s%MNP))
      ALLOCATE ( global_here%NNODECODE(s%MNP),global_here%NODEREP(s%MNP))
      ALLOCATE ( global_here%NNEIGH(s%MNP),global_here%MJU(s%MNP),global_here%NODELE(s%MNP))
      ALLOCATE ( global_here%NIBNODECODE(s%MNP),global_here%NNEIGH_ELEM(s%MNP))
!      ALLOCATE ( TAU0VAR(s%MNP))
      ALLOCATE ( global_here%CH1(s%MNP),global_here%QB(s%MNP),global_here%QA(s%MNP),global_here%SOURSIN(s%MNP))!,EVC(s%MNP))
      ALLOCATE ( global_here%TKXX(s%MNP),global_here%TKYY(s%MNP),global_here%TKXY(s%MNP))
      ALLOCATE ( global_here%AREAS(s%MNE),global_here%SFACDUB(3,s%MNE),global_here%RTEMP2(s%MNE))
      ALLOCATE ( global_here%YDUB(36,s%MNE,8))
      ALLOCATE ( global_here%UU2(max(s%MNP,s%MNE)))
      ALLOCATE ( global_here%VV2(max(s%MNP,s%MNE)))
      ALLOCATE ( global_here%ETA1(max(s%MNP,s%MNE)),global_here%ETA2(max(s%MNP,s%MNE)), global_here%ETAMAX(s%MNP), global_here%MassMax(s%MNE) )
      Allocate ( global_here%tracer(max(s%MNP,s%MNE)), global_here%tracer2(max(s%MNP,s%MNE)), global_here%bed_int(max(s%mnp,s%mne),s%layers), global_here%entrop(5,s%MNE) )
      ALLOCATE ( global_here%CORIF(s%MNP), global_here%dyn_P(max(s%MNP,s%MNE)), global_here%pdg_el(s%MNE) )
      ALLOCATE ( global_here%QU(s%MNP),global_here%QV(s%MNP))
      ALLOCATE ( global_here%LBCODE(s%MNP))
      ALLOCATE ( global_here%CSI(s%MNP),global_here%SII(s%MNP))
      ALLOCATE ( global_here%NODECODE(s%MNP))
      ALLOCATE (global_here%GRAVX(s%MNP),global_here%GRAVY(s%MNP))
      ALLOCATE ( global_here%NIBCNT(s%MNP) )    !  added 7/31/2000 to fix wetdry bug
!sb-added by sb 10/11/2005
      ALLOCATE( global_here%NNDEL(s%MNP) )
      ALLOCATE( global_here%EDFLG(3,s%MNE) )

!--

#ifdef CMPI
      ALLOCATE ( IDUMY(1), DUMY1(1), DUMY2(1) )
#endif

#ifdef CVEC
      ALLOCATE ( global_here%QTEMA(s%MNE,3),global_here%QTEMB(s%MNE,3))
#endif

      IF ( s%C3DDSS) THEN
        ALLOCATE( global_here%AUV11(s%MNP),global_here%AUV12(s%MNP),global_here%AUV13(s%MNP),global_here%AUV14(s%MNP))
      ELSEIF (s%C2DDI) THEN
        ALLOCATE( global_here%AUV11(s%MNP),global_here%AUV12(s%MNP))
        ALLOCATE( global_here%AUVXX(s%MNP),global_here%AUVYY(s%MNP),global_here%AUVXY(s%MNP),global_here%AUVYX(s%MNP))
      ENDIF

      IF ( s%C3D) THEN
       ALLOCATE( global_here%DUU1(s%MNP),global_here%DUV1(s%MNP),global_here%DVV1(s%MNP),global_here%BSX1(s%MNP),global_here%BSY1(s%MNP))
       endif
!
      RETURN
      END SUBROUTINE


      SUBROUTINE ALLOC_MAIN2(s,global_here)
      use sizes
      implicit none
      type (sizes_type) :: s
      type (global_type) :: global_here

!
!     Allocate space for Arrays dimensioned by MNOPE and MNETA
!
      ALLOCATE ( global_here%ESBIN1(s%MNETA),global_here%ESBIN2(s%MNETA))
      ALLOCATE ( global_here%NBDV(s%MNOPE,s%MNETA))
      ALLOCATE ( global_here%NVDLL(s%MNOPE),global_here%NBD(s%MNETA))
      ALLOCATE ( global_here%EMO(s%MNBFR,s%MNETA),global_here%EFA(s%MNBFR,s%MNETA))
      ALLOCATE ( global_here%UMO(s%MNBFR,s%MNETA),global_here%UFA(s%MNBFR,s%MNETA))
      ALLOCATE ( global_here%VMO(s%MNBFR,s%MNETA),global_here%VFA(s%MNBFR,s%MNETA))

      RETURN
      END SUBROUTINE


!
!     Allocate space for nonperiodic zero and nonzero normal flow boundary arrays
!     including barriers
!
      SUBROUTINE ALLOC_MAIN3(s,global_here)
      use sizes
      implicit none
      type (sizes_type) :: s
      type (global_type) :: global_here

      ALLOCATE ( global_here%QN0(s%MNVEL),global_here%QN1(s%MNVEL),global_here%QN2(s%MNVEL))
      ALLOCATE ( global_here%NBV(s%MNVEL),global_here%LBCODEI(s%MNVEL))
      ALLOCATE ( global_here%BNDLEN2O3(s%MNVEL))
      ALLOCATE ( global_here%ME2GW(s%MNVEL))
      ALLOCATE ( global_here%CSII(s%MNVEL),global_here%SIII(s%MNVEL))
      ALLOCATE ( global_here%BARLANHT(s%MNVEL),global_here%BARLANCFSP(s%MNVEL))
      ALLOCATE ( global_here%BARLANHTR(s%MNVEL),global_here%BARLANCFSPR(s%MNVEL))
      ALLOCATE ( global_here%BARINHT(s%MNVEL),global_here%BARINCFSB(s%MNVEL),global_here%BARINCFSP(s%MNVEL))
      ALLOCATE ( global_here%PIPEHT(s%MNVEL),global_here%PIPECOEF(s%MNVEL),global_here%PIPEDIAM(s%MNVEL))
      ALLOCATE ( global_here%IBCONN(s%MNVEL))
      ALLOCATE ( global_here%BARINHTR(s%MNVEL),global_here%BARINCFSBR(s%MNVEL),global_here%BARINCFSPR(s%MNVEL))
      ALLOCATE ( global_here%PIPEHTR(s%MNVEL),global_here%PIPECOEFR(s%MNVEL),global_here%PIPEDIAMR(s%MNVEL))
      ALLOCATE ( global_here%IBCONNR(s%MNVEL),global_here%NTRAN1(s%MNVEL),global_here%NTRAN2(s%MNVEL))
      ALLOCATE ( global_here%BTRAN3(s%MNVEL),global_here%BTRAN4(s%MNVEL),global_here%BTRAN5(s%MNVEL))
      ALLOCATE ( global_here%BTRAN6(s%MNVEL),global_here%BTRAN7(s%MNVEL),global_here%BTRAN8(s%MNVEL))
      ALLOCATE ( global_here%RBARWL1AVG(s%MNVEL),global_here%RBARWL2AVG(s%MNVEL))
      ALLOCATE ( global_here%RPIPEWL1AVG(s%MNVEL),global_here%RPIPEWL2AVG(s%MNVEL))
      ALLOCATE ( global_here%QNIN1(s%MNVEL),global_here%QNIN2(s%MNVEL))
      ALLOCATE ( global_here%NBVV(s%MNBOU,0:s%MNVEL))
      ALLOCATE ( global_here%NVELL(s%MNBOU))
      ALLOCATE ( global_here%SEGTYPE(s%MNBOU))

      RETURN
      END SUBROUTINE


!
!     Allocate space for tidal potential terms 
!
      SUBROUTINE ALLOC_MAIN4a(s,global_here)
      use sizes
      implicit none
      type (sizes_type) :: s
      type (global_type) :: global_here

      ALLOCATE ( global_here%TPK(s%MNTIF),global_here%AMIGT(s%MNTIF),global_here%FFT(s%MNTIF) )
      ALLOCATE ( global_here%FACET(s%MNTIF),global_here%PERT(s%MNTIF),global_here%ETRF(s%MNTIF) )
      ALLOCATE ( global_here%TIPOTAG(s%MNTIF) )

      IF ( s%CTIP ) THEN
        ALLOCATE( global_here%TIP1(s%MNP),global_here%TIP2(s%MNP))
        ENDIF  

      RETURN
      END SUBROUTINE


!
!     Allocate space for Earth load/self-attraction tide 
!
      SUBROUTINE ALLOC_MAIN4b(s,global_here)
      use sizes
      implicit none
      type (sizes_type) :: s
      type (global_type) :: global_here

      ALLOCATE ( global_here%SALTAMP(s%MNTIF,s%MNP),global_here%SALTPHA(s%MNTIF,s%MNP) )
      RETURN
      END SUBROUTINE


!
!     Allocate space for Arrays dimensioned by s%MNBFR      
!
      SUBROUTINE ALLOC_MAIN5(s,global_here)
      use sizes
      implicit none
      type (sizes_type) :: s
       type (global_type) :: global_here

     ALLOCATE ( global_here%AMIG(s%MNBFR),global_here%PER(s%MNBFR))
      ALLOCATE ( global_here%FF(s%MNBFR),global_here%FACE(s%MNBFR))
      ALLOCATE ( global_here%BOUNTAG(s%MNBFR) )

      RETURN
      END SUBROUTINE



!
!     Allocate space for periodic normal flow boundary conditions
!
      SUBROUTINE ALLOC_MAIN6(s,global_here)
      use sizes
      implicit none
      type (sizes_type) :: s
      type (global_type) :: global_here

      ALLOCATE ( global_here%QNAM(s%MNFFR,s%MNVEL),global_here%QNPH(s%MNFFR,s%MNVEL))
      ALLOCATE ( global_here%FBOUNTAG(s%MNFFR))
      ALLOCATE ( global_here%FAMIG(s%MNFFR), global_here%FFF(s%MNFFR), global_here%FFACE(s%MNFFR), global_here%FPER(s%MNFFR) )

      RETURN
      END SUBROUTINE


!
!     Allocate space for arrays used for station elevation output
!
      SUBROUTINE ALLOC_MAIN7(s,global_here)
      use sizes
      implicit none
      type (sizes_type) :: s
      type (global_type) :: global_here

      ALLOCATE ( global_here%NNE(s%MNSTAE),global_here%ET00(s%MNSTAE),global_here%BT00(s%MNSTAE))
      ALLOCATE ( global_here%STAIE1(s%MNSTAE),global_here%STAIE2(s%MNSTAE),global_here%STAIE3(s%MNSTAE))
      ALLOCATE ( global_here%XEL(s%MNSTAE),global_here%YEL(s%MNSTAE),global_here%SLEL(s%MNSTAE),global_here%SFEL(s%MNSTAE))
      RETURN
      END SUBROUTINE


!
!     Allocate space for arrays used for station velocity output
!
      SUBROUTINE ALLOC_MAIN8(s,global_here)
      use sizes
      implicit none
      type (sizes_type) :: s
      type (global_type) :: global_here

      ALLOCATE ( global_here%XEV(s%MNSTAV),global_here%YEV(s%MNSTAV),global_here%SLEV(s%MNSTAV),global_here%SFEV(s%MNSTAV))
      ALLOCATE ( global_here%NNV(s%MNSTAV))
      ALLOCATE ( global_here%UU00(s%MNSTAV),global_here%VV00(s%MNSTAV))
      ALLOCATE ( global_here%STAIV1(s%MNSTAV),global_here%STAIV2(s%MNSTAV),global_here%STAIV3(s%MNSTAV))

      RETURN
      END SUBROUTINE


!
!     Allocate space for arrays used for station concentration output
!
      SUBROUTINE ALLOC_MAIN9(s,global_here)
      use sizes
      implicit none
      type (sizes_type) :: s
      type (global_type) :: global_here

      ALLOCATE ( global_here%XEC(s%MNSTAC),global_here%YEC(s%MNSTAC),global_here%SLEC(s%MNSTAC),global_here%SFEC(s%MNSTAC))
      ALLOCATE ( global_here%NNC(s%MNSTAC))
      ALLOCATE ( global_here%CC00(s%MNSTAC))
      ALLOCATE ( global_here%STAIC1(s%MNSTAC),global_here%STAIC2(s%MNSTAC),global_here%STAIC3(s%MNSTAC))

      RETURN
      END SUBROUTINE


!
!     Allocate space for arrays used for station meteorological output
!
      SUBROUTINE ALLOC_MAIN10(s,global_here)
      use sizes
      implicit none
      type (sizes_type) :: s
      type (global_type) :: global_here

      ALLOCATE ( global_here%XEM(s%MNSTAM),global_here%YEM(s%MNSTAM),global_here%SLEM(s%MNSTAM),global_here%SFEM(s%MNSTAM))
      ALLOCATE ( global_here%NNM(s%MNSTAM))
      ALLOCATE ( global_here%RMU00(s%MNSTAM),global_here%RMV00(s%MNSTAM),global_here%RMP00(s%MNSTAM))
      ALLOCATE ( global_here%STAIM1(s%MNSTAM),global_here%STAIM2(s%MNSTAM),global_here%STAIM3(s%MNSTAM))

      RETURN
      END SUBROUTINE


!
!     Allocate space for Arrays dimensioned by s%s%MNEI   
!
      SUBROUTINE ALLOC_MAIN11(s,global_here)
      use sizes
      implicit none
      type (sizes_type) :: s
      type (global_type) :: global_here

!  Arrays used by JCG iterative solver

      ALLOCATE( global_here%OBCCOEF(s%MNETA,s%MNEI-1),global_here%COEF(s%MNP,s%MNEI))
      ALLOCATE( global_here%IWKSP(3*s%MNP),global_here%WKSP(4*s%MNP+400) )
      ALLOCATE( global_here%IPARM(12),global_here%RPARM(12) )

!  Neighbor Table

      ALLOCATE ( global_here%NEITAB(s%MNP,s%MNEI),global_here%NEIGH_ELEM(s%MNP,s%MNEI))

      RETURN
      END SUBROUTINE


!
!     Allocate space for wind forcing   
!
      SUBROUTINE ALLOC_MAIN12(s,global_here)
      use sizes
      implicit none
      type (sizes_type) :: s
      type (global_type) :: global_here

      ALLOCATE ( global_here%WSX1(s%MNP),global_here%WSY1(s%MNP),global_here%PR1(s%MNP) )
      ALLOCATE ( global_here%WSX2(s%MNP),global_here%WSY2(s%MNP),global_here%PR2(s%MNP) )
      ALLOCATE ( global_here%WVNX1(s%MNP),global_here%WVNY1(s%MNP),global_here%PRN1(s%MNP) )
      ALLOCATE ( global_here%WVNX2(s%MNP),global_here%WVNY2(s%MNP),global_here%PRN2(s%MNP) )
      ALLOCATE ( global_here%RSNX1(s%MNP),global_here%RSNY1(s%MNP),global_here%RSNX2(s%MNP),global_here%RSNY2(s%MNP) )
      ALLOCATE ( global_here%WVNXOUT(s%MNP),global_here%WVNYOUT(s%MNP) )
#ifdef SWAN
!asey 101118: Added the next line for output of radiation stress gradients.
      ALLOCATE ( RSNXOUT(s%MNP),RSNYOUT(s%MNP) )
#endif
      RETURN
      END SUBROUTINE


!
!     Allocate space for bridge piling friction arrays   
!
      SUBROUTINE ALLOC_MAIN13(s,global_here)
      use sizes
      implicit none
      type (sizes_type) :: s
      type (global_type) :: global_here

      ALLOCATE ( global_here%NBNNUM(s%MNP),global_here%BK(s%MNP),global_here%BALPHA(s%MNP),global_here%BDELX(s%MNP) )
      RETURN
      END SUBROUTINE
!
!     Allocate space for harmonic analysis means and variance calculations, this is in 
!     global data because the variables are used in main source, outside of HA analysis
!     subroutines.  This should probably be changed.
!
      SUBROUTINE ALLOC_MAIN14(s,global_here)
      use sizes
      implicit none
      type (sizes_type) :: s
      type (global_type) :: global_here

      ALLOCATE ( global_here%XVELAV(s%MNP),global_here%YVELAV(s%MNP),global_here%XVELVA(s%MNP),global_here%YVELVA(s%MNP) )
      ALLOCATE ( global_here%ELAV(s%MNP),global_here%ELVA(s%MNP) )
      RETURN
      END SUBROUTINE
      
!.....Allocate for wave friction

      SUBROUTINE ALLOC_MAIN15(s,global_here)
      use sizes
      implicit none
      type (sizes_type) :: s
      type (global_type) :: global_here

      ALLOCATE ( global_here%WAVE_T1(s%MNP),global_here%WAVE_H1(s%MNP),global_here%WAVE_A1(s%MNP),global_here%WAVE_D1(s%MNP) )
      ALLOCATE ( global_here%WAVE_T2(s%MNP),global_here%WAVE_H2(s%MNP),global_here%WAVE_A2(s%MNP),global_here%WAVE_D2(s%MNP) )
      ALLOCATE ( global_here%WAVE_T(s%MNP),global_here%WAVE_H(s%MNP),global_here%WAVE_A(s%MNP),global_here%WAVE_D(s%MNP) )
      ALLOCATE ( global_here%WB(s%MNP) )
      END SUBROUTINE

!.....Allocate arrays used for node to element table
!
      SUBROUTINE ALLOC_NNOEL1(s,global_here)
      use sizes
      implicit none
      type (sizes_type) :: s
      type (global_type) :: global_here

      ALLOCATE ( global_here%EL_COUNT(s%MNP) )
      RETURN
      END SUBROUTINE
!
      SUBROUTINE ALLOC_NNOEL2(s,global_here,MAXEL)
      use sizes
      implicit none
      type (sizes_type) :: s
      type (global_type) :: global_here

      integer :: maxel
      ALLOCATE ( global_here%DP_DG(maxel),global_here%DG_ANG(maxel) )
      ALLOCATE ( global_here%NNOEL(s%MNP,maxel),global_here%CENTAB(s%MNP,maxel+1) )
      ALLOCATE ( global_here%ELETAB(s%MNP,maxel+1),global_here%ANGTAB(s%MNP,maxel+1),global_here%dynP_DG(maxel) )
      ALLOCATE ( global_here%ZE_DG(maxel), global_here%QX_DG(maxel), global_here%QY_DG(maxel), global_here%HB_DG(maxel))
      Allocate ( global_here%iota2_DG(maxel), global_here%iota_DG(maxel), global_here%iotaa_DG(maxel) )
      Allocate ( global_here%bed_DG(maxel,s%layers),global_here%bed_N_int(s%layers), global_here%bed_N_ext(s%layers) )
      RETURN
      END SUBROUTINE

!sb
!     Allocate space for Arrays dimensioned by   MNNDEL
!
      SUBROUTINE ALLOC_MAIN16(s,global_here)
      use sizes
      implicit none
      type (sizes_type) :: s
      type (global_type) :: global_here

!.....Node-to-elements table      
      ALLOCATE( global_here%NDEL(s%MNP,s%MNNDEL) )

      RETURN
      END SUBROUTINE


      END MODULE GLOBAL
