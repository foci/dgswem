!***********************************************************************
!     
!     MODULE DG
!     
!     This module declares all the necessary variables for DG
!     
!     Written by Ethan Kubatko (3-10-05)
!     
!-----------------------------------------------------------------------
!     
!     Modification history after v9_sb1
!     
!     v9_sb1     - 08/05    - sb - Parallel runs
!     v9_sb2     - 08/05    - sb - Wetting and drying
!     v9_sb2.2.1 - 08/16/05 - sb - Hard bottom implementation
!     v10_sb5      - Oct 03 - sb - Consolidate Ethan's slope limiter
!     Moving to versioning repo now
!     
!***********************************************************************

      MODULE DG
      
!.....Use appropriate modules

      USE SIZES

!.....Declare integer variables

      INTEGER, TARGET :: DGFLAG,DGHOT,DGHOTSPOOL
      INTEGER DOF,dofh,dofl,dofx
      INTEGER EL
      INTEGER, TARGET :: MNES,artdif,tune_by_hand
      INTEGER IRK
      INTEGER J1, J2, J3,negp_fixed,nagp_fixed
      INTEGER NAGP(8),NCHECK(8), NEGP(8), NEDGES, NRK !42 hardwires for ph=7
      INTEGER NIEDS, NLEDS, NEEDS, NFEDS, NREDS, NEBEDS, NIBEDS
      INTEGER NIBSEG, NEBSEG
      INTEGER MNED, MNLED, MNSED, MNRAED, MNRIED
      INTEGER, TARGET :: MODAL_IC
      INTEGER P_READ
      INTEGER, TARGET :: SLOPEFLAG
      INTEGER test_el
      INTEGER, TARGET :: FLUXTYPE
      INTEGER, TARGET :: RK_STAGE, RK_ORDER
      Integer, TARGET :: padapt,pflag,pl,ph,px,lebesgueP, gflag
      INTEGER pa
      logical init_parser,stblzr
!     
      integer iwrite

!.....Declare real variables

      REAL(SZ) C13, C16
      REAL(SZ), TARGET :: diorism, porosity, SEVDM
      REAL(SZ) DOT, DHB_X, DHB_Y, DPHIDX, DPHIDY
      Real(SZ), TARGET :: slimit,plimit,pflag2con1,pflag2con2
      REAL(SZ) EFA_GP, EMO_GP,slimit1,slimit2,slimit3
      REAL(SZ) EL_ANG,slimit4, bg_dif,trc_dif,slimit5
      REAL(SZ) FG_L,l2er_global,temperg
      REAL(SZ), TARGET :: slope_weight
      REAL(SZ) HB_IN, HB_EX, H_TRI
      REAL(SZ) MAG1, MAG2
      REAL(SZ) NX, NY 
      REAL(SZ), TARGET :: kappa,s0,uniform_dif
      REAL(SZ) SFAC_IN,SFAC_EX
      REAL(SZ) RAMPDG
      REAL(SZ) S1, S2, SAV, SOURCE_X, SOURCE_Y
      REAL(SZ) TIMEDG, TIMEH_DG, TK
      REAL(SZ) QX_EX, QX_IN, QY_EX, QY_IN
      REAL(SZ) QNAM_GP, QNPH_GP
      REAL(SZ) SL2_M, SL2_NYU
      REAL(SZ) SL3_MD, EVMAvg, SEVDMAvg
      REAL(SZ) UMAG
      REAL(SZ) WSX_GP, WSY_GP
      REAL(SZ) ZE_EX, ZE_IN, QMag_IN, QMag_EX
      Real(SZ) subphi_IN,subphi_EX
      Real(SZ) iota_EX, iota_IN,iota2_EX, iota2_IN

!     sb...Wetting and drying
      INTEGER, ALLOCATABLE :: WDFLG(:) ! = 1 if wet, =0 if dry
      INTEGER, ALLOCATABLE :: WDFLG_TMP(:)
      INTEGER, ALLOCATABLE :: DOFW(:)
      INTEGER, ALLOCATABLE :: EL_UPDATED(:)
      REAL(SZ), ALLOCATABLE :: LEDGE_NVEC(:,:,:)

                                ! LEDGE_NVEC(1,1:3,1:NELEM) = whether a node is on the land boundary
                                ! LEDGE_NVEC(2,1:3,1:NELEM) = x component of the normal vector
                                ! LEDGE_NVEC(3,1:3,1:NELEM) = y component of the normal vector

!Declare some stuff for function parsing for bed load

      CHARACTER (LEN=7), DIMENSION(4), PARAMETER :: varx = [ CHARACTER(7) :: &
           'ZE_ROE', 'QX_ROE','QY_ROE','bed_ROE']
      CHARACTER (LEN=7), DIMENSION(4), PARAMETER :: vary = [ CHARACTER(7) :: &
           'ZE_ROE', 'QX_ROE','QY_ROE','bed_ROE']

!      CHARACTER (LEN=*), DIMENSION(4),  PARAMETER :: varx  = (/ 'ZE_ROE', 'QX_ROE','QY_ROE','bed_ROE' /)
!      CHARACTER (LEN=*), DIMENSION(4),  PARAMETER :: vary  = (/ 'ZE_ROE', 'QX_ROE','QY_ROE','bed_ROE' /)

      CHARACTER (LEN=200) funcx(4), funcy(4)
      Real(sz)  valx(4), valy(4)

!.....Declare real variable arrays

      REAL(SZ) DRPSI(3), DSPSI(3)
      REAL(SZ) VEC1(2), VEC2(2)
      
!.....Declare allocatable integer arrays

      INTEGER, ALLOCATABLE :: DOFS(:), PCOUNT(:), PDG(:)
      INTEGER, ALLOCATABLE :: NCOUNT(:)
      INTEGER, ALLOCATABLE :: NEDEL(:,:), NEDNO(:,:), NEDSD(:,:)
      INTEGER, ALLOCATABLE :: NIEDN(:), NLEDN(:), NEEDN(:)
      INTEGER, ALLOCATABLE :: NFEDN(:), NREDN(:), NEBEDN(:), NIBEDN(:)
      INTEGER, ALLOCATABLE :: NIBSEGN(:,:)
      INTEGER, ALLOCATABLE :: NEBSEGN(:)
      INTEGER, ALLOCATABLE :: EL_NBORS(:,:)
      INTEGER, ALLOCATABLE :: BACKNODES(:,:)
      INTEGER, ALLOCATABLE :: MARK(:)


!.....Declare allocatable real arrays
      
      REAL(SZ), ALLOCATABLE :: ATVD(:,:), BTVD(:,:), CTVD(:,:)
      REAL(SZ), ALLOCATABLE :: DTVD(:), MAX_BOA_DT(:),e1(:),balance(:)
      Real(SZ), Allocatable :: RKC_T(:),RKC_U(:),RKC_Tprime(:)
      Real(SZ), Allocatable :: RKC_Tdprime(:),RKC_a(:),RKC_b(:),RKC_c(:)
      Real(SZ), Allocatable :: RKC_mu(:),RKC_tildemu(:),RKC_nu(:),RKC_gamma(:)
      REAL(SZ), ALLOCATABLE :: BATH(:,:,:), DBATHDX(:,:,:), DBATHDY(:,:,:)
      REAL(SZ), ALLOCATABLE :: SFAC_ELEM(:,:,:)
      REAL(SZ), ALLOCATABLE :: BATHED(:,:,:,:),SFACED(:,:,:,:)
      REAL(SZ), ALLOCATABLE :: COSNX(:), SINNX(:)
      REAL(SZ), ALLOCATABLE :: DP_NODE(:,:,:)
      REAL(SZ), ALLOCATABLE :: DP_VOL(:,:)
      REAL(SZ), ALLOCATABLE :: DRPHI(:,:,:), DSPHI(:,:,:)
      REAL(SZ), ALLOCATABLE :: DRDX(:), DSDX(:), DRDY(:), DSDY(:)
      REAL(SZ), ALLOCATABLE :: DXPHI2(:,:,:), DYPHI2(:,:,:), PHI2(:,:,:)
      REAL(SZ), ALLOCATABLE :: EFA_DG(:,:,:), EMO_DG(:,:,:)
      REAL(SZ), ALLOCATABLE :: UFA_DG(:,:,:), UMO_DG(:,:,:)
      REAL(SZ), ALLOCATABLE :: VFA_DG(:,:,:), VMO_DG(:,:,:)
      REAL(SZ), ALLOCATABLE :: XLEN(:)
      REAL(SZ), ALLOCATABLE :: HB(:,:,:)
      REAL(SZ), ALLOCATABLE :: MANN(:,:)
      REAL(SZ), ALLOCATABLE :: IBHT(:), EBHT(:)
      REAL(SZ), ALLOCATABLE :: EBCFSP(:), IBCFSP(:), IBCFSB(:)
      REAL(SZ), ALLOCATABLE :: JACOBI(:,:,:,:)
      REAL(SZ), ALLOCATABLE :: M_INV(:,:),phi_edge_fixed(:,:,:)
      REAL(SZ), ALLOCATABLE :: PHI_AREA(:,:,:), PHI_EDGE(:,:,:,:)
      REAL(SZ), ALLOCATABLE :: PHI_CENTER(:,:), PHI_CORNER(:,:,:)
      REAL(SZ), ALLOCATABLE :: PHI_CHECK(:,:,:)
      REAL(SZ), ALLOCATABLE :: PHI_CORNER1(:,:,:,:)
      REAL(SZ), ALLOCATABLE :: PHI_MID(:,:,:)
      REAL(SZ), ALLOCATABLE :: PHI_INTEGRATED(:,:)
      REAL(SZ), ALLOCATABLE :: PSI_CHECK(:,:)
      REAL(SZ), ALLOCATABLE :: PSI1(:,:), PSI2(:,:), PSI3(:,:)
      REAL(SZ), ALLOCATABLE :: Q_HAT(:)
      REAL(SZ), ALLOCATABLE :: QIB(:)
      REAL(SZ), ALLOCATABLE :: QX(:,:,:), QY(:,:,:), ZE(:,:,:)
      Real(SZ), ALLOCATABLE, target:: bed(:,:,:,:)
      Real(SZ), Allocatable :: dynP(:,:,:),dynP_MAX(:), dynP_MIN(:)
      Real(SZ), Allocatable :: iota(:,:,:),iotaa(:,:,:),iota2(:,:,:)
      Real(SZ), Allocatable :: iota_MAX(:),iota_MIN(:),iotaa2(:,:,:),iotaa3(:,:,:)
      Real(SZ), Allocatable :: iota2_MAX(:),iota2_MIN(:)
      Real(SZ), pointer :: arrayfix(:,:,:)
      REAL(SZ), ALLOCATABLE :: CORI_EL(:), FRIC_EL(:)
      REAL(SZ), ALLOCATABLE :: ZE_MAX(:),ZE_MIN(:),DPE_MIN(:)
      REAL(SZ), ALLOCATABLE :: WATER_DEPTH_OLD(:,:),WATER_DEPTH(:,:)
      REAL(SZ), ALLOCATABLE :: ADVECTQX(:),ADVECTQY(:)
      REAL(SZ), ALLOCATABLE :: SOURCEQX(:),SOURCEQY(:)
      REAL(SZ), ALLOCATABLE :: LZ(:,:,:,:),MZ(:,:,:,:)
      Real(SZ), Allocatable :: HZ(:,:,:,:),TZ(:,:,:,:)
      REAL(SZ), ALLOCATABLE :: QNAM_DG(:,:,:), QNPH_DG(:,:,:)
      REAL(SZ), ALLOCATABLE :: RHS_ZE(:,:,:), RHS_bed(:,:,:,:)
      REAL(SZ), ALLOCATABLE :: RHS_QX(:,:,:),RHS_QY(:,:,:)
      Real(SZ), Allocatable :: RHS_iota(:,:,:),RHS_iota2(:,:,:), RHS_dynP(:,:,:)
      Real(SZ), Allocatable :: RHS_bed_IN(:,:),RHS_bed_EX(:,:), bed_HAT_O(:)
      REAL(SZ), ALLOCATABLE :: XAGP(:,:), YAGP(:,:), WAGP(:,:)
      REAL(SZ), ALLOCATABLE :: XEGP(:,:), YEGP(:,:), WEGP(:,:)
      REAL(SZ), ALLOCATABLE :: SL3(:,:)
      REAL(SZ), ALLOCATABLE :: XBC(:), YBC(:)
      REAL(SZ), ALLOCATABLE :: XFAC(:,:,:,:), YFAC(:,:,:,:), SRFAC(:,:,:,:)
      REAL(SZ), ALLOCATABLE :: EDGEQ(:,:,:,:)
      REAL(SZ), ALLOCATABLE :: PHI(:), DPHIDZ1(:), DPHIDZ2(:)
      REAL(SZ), ALLOCATABLE :: PHI_STAE(:,:), PHI_STAV(:,:)
      Real(SZ), Allocatable :: bed_IN(:),bed_EX(:),bed_HAT(:)

!.....These (below) are defined in prep_slopelim.F

      Integer :: lim_count,lim_count_roll
      Integer,Allocatable :: fact(:),focal_neigh(:,:),focal_up(:),bi(:),bj(:)

      Real(SZ),Allocatable :: XBCb(:),YBCb(:),xi1(:,:),xi2(:,:)
      Real(SZ),Allocatable :: xtransform(:,:),ytransform(:,:)
      Real(SZ),Allocatable :: xi1BCb(:),xi2BCb(:),xi1vert(:,:)
      Real(SZ),Allocatable :: xi2vert(:,:),xtransformv(:,:),ytransformv(:,:)
      Real(SZ),Allocatable :: XBCv(:,:),YBCv(:,:)
      Real(SZ),Allocatable :: xi1BCv(:,:),xi2BCv(:,:)
      Real(SZ),Allocatable :: Area_integral(:,:,:)
      Real(SZ),Allocatable :: f(:,:,:,:),g0(:,:,:,:),varsigma0(:,:,:,:)
      Real(SZ),Allocatable :: fv(:,:,:,:),g0v(:,:,:,:),varsigma0v(:,:,:,:)
      Real(SZ),Allocatable :: var2sigmag(:,:,:),var2sigmav(:,:,:)
      Real(SZ),Allocatable :: Nmatrix(:,:,:,:),NmatrixInv(:,:,:,:)
      Real(SZ),Allocatable :: deltx(:),delty(:),pmatrix(:,:,:)

!.....These (below) are defined in slopelimiter.F 


      Real(SZ),Allocatable :: ZEmin(:,:),ZEmax(:,:),QXmin(:,:)
      Real(SZ),Allocatable :: QXmax(:,:),QYmin(:,:),QYmax(:,:)
      Real(SZ),Allocatable :: iotamin(:,:),iotamax(:,:)
      Real(SZ),Allocatable :: iota2min(:,:),iota2max(:,:)

#ifdef SLOPEALL
      Real(SZ),Allocatable :: ZEtaylor(:,:,:),QXtaylor(:,:,:),QYtaylor(:,:,:)
      Real(SZ),Allocatable :: iotataylor(:,:,:),iota2taylor(:,:,:)
      Real(SZ),Allocatable :: ZEtaylorvert(:,:,:),QXtaylorvert(:,:,:)
      Real(SZ),Allocatable :: QYtaylorvert(:,:,:),iotataylorvert(:,:,:)
      Real(SZ),Allocatable :: iota2taylorvert(:,:,:)
      Real(SZ),Allocatable :: alphaZE0(:,:,:),alphaQX0(:,:,:)
      Real(SZ),Allocatable :: alphaQY0(:,:,:),alphaiota0(:,:,:)
      Real(SZ),Allocatable :: alphaiota20(:,:,:)
      Real(SZ),Allocatable :: alphaZE(:,:),alphaQX(:,:),alphaQY(:,:)
      Real(SZ),Allocatable :: alphaiota(:,:),alphaiota2(:,:)
      Real(SZ),Allocatable :: alphaZEm(:,:),alphaQXm(:,:),alphaQYm(:,:)
      Real(SZ),Allocatable :: alphaiotam(:,:),alphaiota2m(:,:)
      Real(SZ),Allocatable :: alphaZE_max(:,:),alphaQX_max(:,:)
      Real(SZ),Allocatable :: alphaQY_max(:,:)
      Real(SZ),Allocatable :: alphaiota_max(:,:),alphaiota2_max(:,:)
      Real(SZ),Allocatable :: limitZE(:,:),limitQX(:,:),limitQY(:,:)
      Real(SZ),Allocatable :: limitiota(:,:),limitiota2(:,:)
      Real(SZ),Allocatable :: ZEconst(:,:),QXconst(:,:),QYconst(:,:)
      Real(SZ),Allocatable :: iotaconst(:,:),iota2const(:,:)
#endif

    
!**********************END OF DATA DECLARATIONS ***********************

      CONTAINS
      
!.....Set edge array sizes

      SUBROUTINE ALLOC_EDGES0(s)
        type (sizes_type) :: s
      ALLOCATE ( IBHT(3*s%MNE), EBHT(3*s%MNE) )
      ALLOCATE ( EBCFSP(3*s%MNE), IBCFSP(3*s%MNE), IBCFSB(3*s%MNE) )
      ALLOCATE ( BACKNODES(2,3*s%MNE) )
      RETURN
      END SUBROUTINE

      SUBROUTINE ALLOC_EDGES1()
      ALLOCATE ( NEDNO(2,MNED), NEDEL(2,MNED), NEDSD(2,MNED) )
      ALLOCATE ( NIBSEGN(2,MNED) )
      ALLOCATE ( NEBSEGN(MNED) )
      ALLOCATE ( NIEDN(MNED), NLEDN(MNED), NEEDN(MNED) )
      ALLOCATE ( NFEDN(MNED),  NREDN(MNED), NIBEDN(MNED), NEBEDN(MNED) )
      ALLOCATE ( NCOUNT(MNED) )
      ALLOCATE ( COSNX(MNED), SINNX(MNED), XLEN(MNED) )
      ALLOCATE ( Q_HAT(MNED) )
      RETURN
      END SUBROUTINE
      
!.....Set DG SWE array sizes

      SUBROUTINE ALLOC_DG1(MNBFR)
        integer :: MNBFR
      ALLOCATE ( EFA_DG(MNBFR,NEEDS+2,2), EMO_DG(MNBFR,NEEDS+2,2) )
      ALLOCATE ( UFA_DG(MNBFR,NEEDS+2,2), UMO_DG(MNBFR,NEEDS+2,2) )
      ALLOCATE ( VFA_DG(MNBFR,NEEDS+2,2), VMO_DG(MNBFR,NEEDS+2,2) )
      RETURN
      END SUBROUTINE
      
      SUBROUTINE ALLOC_DG2(MNFFR)
        integer :: MNFFR
      ALLOCATE ( QNAM_DG(MNFFR,NFEDS,2), QNPH_DG(MNFFR,NFEDS,2) )
      RETURN
      END SUBROUTINE
      
      SUBROUTINE ALLOC_DG3(MNP)
        integer :: MNP
      ALLOCATE ( QIB(MNP) )
      RETURN
      END SUBROUTINE

      SUBROUTINE ALLOC_DG4(s)
        type (sizes_type) :: s
!sb-20070228 NRK+1-->NRK+2 --- XX(:,:,NRK+2) will be used by slope limiter
      ALLOCATE ( HB(DOFH,S%MNE,NRK+2),e1(s%layers+5),balance(s%layers+5) )
      ALLOCATE ( MANN(DOFH,S%MNE) ) !,arrayfix(DOFH,S%MNE,NRK+2) )
      ALLOCATE ( QY(DOFH,S%MNE,NRK+2), QX(DOFH,S%MNE,NRK+2) )
      ALLOCATE ( ZE(DOFH,S%MNE,NRK+2), bed(DOFH,S%MNE,NRK+2,s%layers) )
      Allocate ( iota(DOFH,S%MNE,NRK+2),iota2(DOFH,S%MNE,NRK+2) ) 
      Allocate ( dynP(DOFH,S%MNE,NRK+2) ) 
      Allocate ( iotaa(DOFH,S%MNE,NRK+2),iotaa2(DOFH,S%MNE,NRK+2) )
      Allocate ( iotaa3(DOFH,S%MNE,NRK+2) )
!sb-20060711 For wet/dry
      ALLOCATE ( ZE_MAX(S%MNE),ZE_MIN(S%MNE),DPE_MIN(S%MNE) )
      Allocate ( iota_MAX(S%MNE),iota_MIN(S%MNE),iota2_MAX(S%MNE),iota2_MIN(S%MNE) )
      Allocate ( dynP_MAX(S%MNE),dynP_MIN(S%MNE) )
      ALLOCATE ( WATER_DEPTH(S%MNE,3), WATER_DEPTH_OLD(S%MNE,3))
      ALLOCATE ( ADVECTQX(S%MNE), ADVECTQY(S%MNE))
      ALLOCATE ( SOURCEQX(S%MNE),SOURCEQY(S%MNE))
      ALLOCATE ( MARK(S%MNE))
!em-2012 for sediment
      Allocate ( bed_IN(s%layers),bed_EX(s%layers),bed_HAT(s%layers) )

!--
!sb-20070101
      ALLOCATE ( LZ(DOFH,2,2,S%MNE),MZ(DOFH,2,s%layers,S%MNE) )
      Allocate ( HZ(DOFH,2,2,S%MNE),TZ(DOFH,2,2,S%MNE) )
!--
      ALLOCATE ( RHS_QX(DOFH,S%MNE,NRK), RHS_QY(DOFH,S%MNE,NRK) )
      ALLOCATE ( RHS_ZE(DOFH,S%MNE,NRK), RHS_bed(DOFH,S%MNE,NRK,s%layers) )
      Allocate ( RHS_iota(DOFH,S%MNE,NRK),RHS_iota2(DOFH,S%MNE,NRK) )
      Allocate ( RHS_dynP(DOFH,S%MNE,NRK) )
      Allocate ( RHS_bed_IN(dofh,s%layers), RHS_bed_EX(dofh,s%layers) )
      Allocate ( bed_HAT_O(s%layers) )
      ALLOCATE ( DRDX(S%MNE), DSDX(S%MNE), DRDY(S%MNE), DSDY(S%MNE) )
      ALLOCATE ( CORI_EL(S%MNE), FRIC_EL(S%MNE) )
      ALLOCATE ( PHI(DOFH), DPHIDZ1(DOFH), DPHIDZ2(DOFH) )
      ALLOCATE ( DOFS(S%MNE), PCOUNT(S%MNE) )
      ALLOCATE ( PDG(S%MNP) )
      RETURN
      END SUBROUTINE

!.....Set RK time scheme parameters array sizes

      SUBROUTINE ALLOC_RK()
      ALLOCATE( ATVD(NRK,NRK), BTVD(NRK,NRK), CTVD(NRK,NRK)  )
      ALLOCATE( DTVD(NRK), MAX_BOA_DT(NRK) )
      Allocate( RKC_T(0:nrk),RKC_U(0:nrk),RKC_Tprime(0:nrk) )
      Allocate( RKC_Tdprime(0:nrk),RKC_a(0:nrk),RKC_b(0:nrk),RKC_c(0:nrk) ) 
      Allocate( RKC_mu(0:nrk),RKC_tildemu(0:nrk),RKC_nu(0:nrk),RKC_gamma(0:nrk) )
      RETURN
      END SUBROUTINE ALLOC_RK

!.....Set sizes for arrays used in orthobasis

      SUBROUTINE ALLOC_JACOBI()
      ALLOCATE ( JACOBI(ph+1,2*ph+3,2,NAGP(ph)+1) )
      ALLOCATE ( DXPHI2(ph+1,ph+1,NAGP(ph)+1),DYPHI2(ph+1,ph+1,NAGP(ph)+1) )
      ALLOCATE ( PHI2(ph+1,ph+1,NAGP(ph)+1) )
      ALLOCATE ( PHI_CORNER1(ph+1,ph+1,3,ph) )
      RETURN
      END SUBROUTINE
      
!.....Set sizes for arrays for area integrals
      
      SUBROUTINE ALLOC_AREA_GAUSS(s)
        type (sizes_type) :: s
      ALLOCATE ( XAGP(NAGP(ph),ph),YAGP(NAGP(ph),ph),WAGP(NAGP(ph),ph) ) 
      ALLOCATE ( PHI_AREA(DOFH,NAGP(ph)+1,ph ) )
      ALLOCATE ( DSPHI(DOFH,NAGP(ph)+1,ph),DRPHI(DOFH,NAGP(ph)+1,ph) )
      ALLOCATE ( PHI_CORNER(DOFH,3,ph),PHI_MID(DOFH,3,ph) )
      ALLOCATE ( PHI_CENTER(DOFH,DOFH) )
      ALLOCATE ( PSI1(NAGP(ph),ph),PSI2(NAGP(ph),ph),PSI3(NAGP(ph),ph) )
      ALLOCATE ( BATH(NAGP(ph),S%MNE,ph),DBATHDX(NAGP(ph),S%MNE,ph) )
      Allocate ( DBATHDY(NAGP(ph),S%MNE,ph) )
      ALLOCATE ( SFAC_ELEM(NAGP(ph),S%MNE,ph) )
      ALLOCATE ( XFAC(DOFH,NAGP(ph),S%MNE,ph), YFAC(DOFH,NAGP(ph),S%MNE,ph) )
      ALLOCATE ( SRFAC(DOFH,NAGP(ph),S%MNE,ph) )
      RETURN
      END SUBROUTINE
      
!.....Set sizes for arrays for edge integrals

      SUBROUTINE ALLOC_EDGE_GAUSS(s)
        type (sizes_type) :: s
      ALLOCATE ( XEGP(NEGP(ph),ph), WEGP(NEGP(ph),ph) )
      ALLOCATE ( PHI_EDGE(DOFH,NEGP(ph)+1,3,ph) )
      ALLOCATE ( M_INV(DOFH,ph) )
      ALLOCATE ( BATHED(NEGP(ph),3,S%MNE,ph),SFACED(NEGP(ph),3,S%MNE,ph) )
      ALLOCATE ( EDGEQ(DOFH,NEGP(ph),3,ph) )
      RETURN
      END SUBROUTINE

!.....Set sizes for the arrays for the slope limiter
!.....slopelim arrays

      SUBROUTINE ALLOC_SLOPELIM(s)
        type (sizes_type) :: s
      ALLOCATE ( XBC(S%MNE), YBC(S%MNE) )
      ALLOCATE ( EL_NBORS(4,S%MNE) )
      ALLOCATE ( SL3(3,S%MNE) )

!.....These are defined in prep_slopelim.F

      Allocate ( fact(0:ph) ,focal_neigh(S%MNE,3*S%MNEI),focal_up(S%MNE),bi(dofh),bj(dofh) )

      Allocate ( XBCb(S%MNE),YBCb(S%MNE),xi1(S%MNE,NAGP(ph)),xi2(S%MNE,NAGP(ph)) )
      Allocate ( xtransform(S%MNE,NAGP(ph)),ytransform(S%MNE,NAGP(ph)) )
      Allocate ( xi1BCb(S%MNE),xi2BCb(S%MNE),xi1vert(S%MNE,3) )
      Allocate ( xi2vert(S%MNE,3),xtransformv(S%MNE,3),ytransformv(S%MNE,3) )
      Allocate ( XBCv(S%MNE,S%MNE),YBCv(S%MNE,S%MNE) )
      Allocate ( xi1BCv(S%MNE,S%MNE),xi2BCv(S%MNE,S%MNE),Area_integral(S%MNE,0:ph,0:ph) )
      Allocate ( f(S%MNE,NAGP(ph),0:ph,0:ph),g0(S%MNE,NAGP(ph),0:ph,0:ph) )
      Allocate ( fv(S%MNE,3,0:ph,0:ph),g0v(S%MNE,3,0:ph,0:ph) )
      Allocate ( varsigma0(S%MNE,NAGP(ph),0:ph,0:ph) ) 
      Allocate ( varsigma0v(S%MNE,3,0:ph,0:ph) )
      Allocate ( pmatrix(S%MNE,dofh,dofh), var2sigmag(S%MNE,NAGP(ph),dofh) )
      Allocate ( Nmatrix(S%MNE,dofh,dofh,dofh),NmatrixInv(S%MNE,dofh,dofh,dofh) )
      Allocate ( deltx(S%MNE),delty(S%MNE),var2sigmav(S%MNE,3,dofh))

!.....These (below) are defined in slopelimiter.F (slopelimiter4)
      
      Allocate ( ZEmin(S%MNP,dofh),ZEmax(S%MNP,dofh),QXmin(S%MNP,dofh) )
      Allocate ( QXmax(S%MNP,dofh),QYmin(S%MNP,dofh),QYmax(S%MNP,dofh) )
      Allocate ( iotamin(S%MNP,dofh),iotamax(S%MNP,dofh) )
      Allocate ( iota2min(S%MNP,dofh),iota2max(S%MNP,dofh) )

#ifdef SLOPEALL
      Allocate ( ZEtaylor(S%MNE,dofh,1),QXtaylor(S%MNE,dofh,1) )
      Allocate ( iotataylor(S%MNE,dofh,1),iota2taylor(S%MNE,dofh,1) 
      Allocate ( ZEtaylorvert(S%MNE,dofh,3),QXtaylorvert(S%MNE,dofh,3) )
      Allocate ( QYtaylorvert(S%MNE,dofh,3),iotataylorvert(S%MNE,dofh,3) )
      Allocate ( iota2taylorvert(S%MNE,dofh,3), QYtaylor(S%MNE,dofh,1) )
      Allocate ( alphaZE0(S%MNE,dofh,3),alphaQX0(S%MNE,dofh,3) )
      Allocate ( alphaQY0(S%MNE,dofh,3),alphaiota0(S%MNE,dofh,3) )
      Allocate ( alphaiota20(S%MNE,dofh,3) )
      Allocate ( alphaZE(S%MNE,dofh),alphaQX(S%MNE,dofh),alphaQY(S%MNE,dofh) )
      Allocate ( alphaiota(S%MNE,dofh),alphaiota2(S%MNE,dofh) )
      Allocate ( alphaZEm(S%MNE,dofh),alphaQXm(S%MNE,dofh),alphaQYm(S%MNE,dofh) )
      Allocate ( alphaiotam(S%MNE,dofh),alphaiota2m(S%MNE,dofh) )
      Allocate ( alphaZE_max(S%MNE,dofh),alphaQX_max(S%MNE,dofh) )
      Allocate ( alphaQY_max(S%MNE,dofh) )
      Allocate ( alphaiota_max(S%MNE,dofh),alphaiota2_max(S%MNE,dofh) )
      Allocate ( limitZE(S%MNE,dofh),limitQX(S%MNE,dofh),limitQY(S%MNE,dofh) )
      Allocate ( limitiota(S%MNE,dofh),limitiota2(S%MNE,dofh) )
      Allocate ( ZEconst(S%MNE,dofh),QXconst(S%MNE,dofh),QYconst(S%MNE,dofh) )
      Allocate ( iotaconst(S%MNE,dofh),iota2const(S%MNE,dofh) )
#endif

      RETURN
      END SUBROUTINE

!sb...Set sizes for arrays for wetting and drying
      SUBROUTINE ALLOC_DG_WETDRY(s)
        type (sizes_type) :: s
      ALLOCATE ( WDFLG(S%MNE), DOFW(S%MNE) )
      ALLOCATE ( EL_UPDATED(S%MNE) )
      ALLOCATE ( WDFLG_TMP(S%MNE) )
      ALLOCATE ( LEDGE_NVEC(3,3,S%MNE) )
      ALLOCATE ( DP_VOL(S%MNE,ph) )
      ALLOCATE ( PHI_INTEGRATED(DOFH,ph) )
      ALLOCATE ( PHI_CHECK(DOFH,NCHECK(ph),ph) )
      ALLOCATE ( DP_NODE(NCHECK(ph),S%MNE,ph) )
      ALLOCATE ( PSI_CHECK(3,12*3) )
      RETURN
      END SUBROUTINE
      

      SUBROUTINE ALLOC_STAE(L)
        integer :: L
      ALLOCATE ( PHI_STAE(DOFH,L) )
      RETURN
      END SUBROUTINE
      
      SUBROUTINE ALLOC_STAV(L)
        integer :: L
      ALLOCATE ( PHI_STAV(DOFH,L) )
      RETURN
      END SUBROUTINE
      

      END MODULE DG


