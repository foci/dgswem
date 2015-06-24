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

      CHARACTER (LEN=7), DIMENSION(4), PARAMETER :: varx = [ CHARACTER(7) :: &
           'ZE_ROE', 'QX_ROE','QY_ROE','bed_ROE']
      CHARACTER (LEN=7), DIMENSION(4), PARAMETER :: vary = [ CHARACTER(7) :: &
           'ZE_ROE', 'QX_ROE','QY_ROE','bed_ROE']

      type dg_type

!.....Declare integer variables
         
      INTEGER :: DGFLAG,DGHOT,DGHOTSPOOL
      INTEGER :: MNES,artdif,tune_by_hand
      INTEGER :: MODAL_IC
      INTEGER :: SLOPEFLAG
      INTEGER :: FLUXTYPE
      INTEGER :: RK_STAGE, RK_ORDER
      Integer :: padapt,pflag,pl,ph,px,lebesgueP, gflag
      INTEGER DOF,dofh,dofl,dofx
      INTEGER EL
      INTEGER IRK
      INTEGER J1, J2, J3,negp_fixed,nagp_fixed
      INTEGER NAGP(8),NCHECK(8), NEGP(8), NEDGES, NRK !42 hardwires for ph=7
      INTEGER NIEDS, NLEDS, NEEDS, NFEDS, NREDS, NEBEDS, NIBEDS
      INTEGER NIBSEG, NEBSEG
      INTEGER MNED, MNLED, MNSED, MNRAED, MNRIED
      INTEGER P_READ
      INTEGER test_el
      INTEGER pa
      logical init_parser,stblzr
!     
      integer iwrite

!.....Declare real variables

      REAL(SZ) :: diorism, porosity, SEVDM
      Real(SZ) :: slimit,plimit,pflag2con1,pflag2con2
      REAL(SZ) :: slope_weight
      REAL(SZ) :: kappa,s0,uniform_dif

      REAL(SZ) C13, C16
      REAL(SZ) DOT, DHB_X, DHB_Y, DPHIDX, DPHIDY
      REAL(SZ) EFA_GP, EMO_GP,slimit1,slimit2,slimit3
      REAL(SZ) EL_ANG,slimit4, bg_dif,trc_dif,slimit5
      REAL(SZ) FG_L,l2er_global,temperg
      REAL(SZ) HB_IN, HB_EX, H_TRI
      REAL(SZ) MAG1, MAG2
      REAL(SZ) NX, NY 
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
#ifdef SEDLAY
      Real(SZ), ALLOCATABLE, target:: bed(:,:,:,:) !TODO: change this
#endif
      Real(SZ), Allocatable :: dynP(:,:,:),dynP_MAX(:), dynP_MIN(:)
      Real(SZ), Allocatable :: iota(:,:,:),iotaa(:,:,:),iota2(:,:,:)
      Real(SZ), Allocatable :: iota_MAX(:),iota_MIN(:),iotaa2(:,:,:),iotaa3(:,:,:)
      Real(SZ), Allocatable :: iota2_MAX(:),iota2_MIN(:)
      Real(SZ), pointer :: arrayfix(:,:,:) !TODO change this
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

   end type dg_type
    
!**********************END OF DATA DECLARATIONS ***********************

      CONTAINS
      
!.....Set edge array sizes

      SUBROUTINE ALLOC_EDGES0(s,dg)
        type (sizes_type) :: s
        type (dg_type) :: dg
      ALLOCATE ( dg%IBHT(3*s%MNE), dg%EBHT(3*s%MNE) )
      ALLOCATE ( dg%EBCFSP(3*s%MNE), dg%IBCFSP(3*s%MNE), dg%IBCFSB(3*s%MNE) )
      ALLOCATE ( dg%BACKNODES(2,3*s%MNE) )
      RETURN
      END SUBROUTINE

      SUBROUTINE ALLOC_EDGES1(dg)
        type (dg_type) :: dg
      ALLOCATE ( dg%NEDNO(2,dg%MNED), dg%NEDEL(2,dg%MNED), dg%NEDSD(2,dg%MNED) )
      ALLOCATE ( dg%NIBSEGN(2,dg%MNED) )
      ALLOCATE ( dg%NEBSEGN(dg%MNED) )
      ALLOCATE ( dg%NIEDN(dg%MNED), dg%NLEDN(dg%MNED), dg%NEEDN(dg%MNED) )
      ALLOCATE ( dg%NFEDN(dg%MNED),  dg%NREDN(dg%MNED), dg%NIBEDN(dg%MNED), dg%NEBEDN(dg%MNED) )
      ALLOCATE ( dg%NCOUNT(dg%MNED) )
      ALLOCATE ( dg%COSNX(dg%MNED), dg%SINNX(dg%MNED), dg%XLEN(dg%MNED) )
      ALLOCATE ( dg%Q_HAT(dg%MNED) )
      RETURN
      END SUBROUTINE
      
!.....Set DG SWE array sizes

      SUBROUTINE ALLOC_DG1(dg,MNBFR)
        type (dg_type) :: dg
        integer :: MNBFR
      ALLOCATE ( dg%EFA_DG(MNBFR,dg%NEEDS+2,2), dg%EMO_DG(MNBFR,dg%NEEDS+2,2) )
      ALLOCATE ( dg%UFA_DG(MNBFR,dg%NEEDS+2,2), dg%UMO_DG(MNBFR,dg%NEEDS+2,2) )
      ALLOCATE ( dg%VFA_DG(MNBFR,dg%NEEDS+2,2), dg%VMO_DG(MNBFR,dg%NEEDS+2,2) )
      RETURN
      END SUBROUTINE
      
      SUBROUTINE ALLOC_DG2(dg,MNFFR)
        type (dg_type) :: dg
        integer :: MNFFR
      ALLOCATE ( dg%QNAM_DG(MNFFR,dg%NFEDS,2), dg%QNPH_DG(MNFFR,dg%NFEDS,2) )
      RETURN
      END SUBROUTINE
      
      SUBROUTINE ALLOC_DG3(dg,MNP)
        type (dg_type) :: dg
        integer :: MNP
      ALLOCATE ( dg%QIB(MNP) )
      RETURN
      END SUBROUTINE

      SUBROUTINE ALLOC_DG4(s,dg)
        type (sizes_type) :: s
        type (dg_type) :: dg
!sb-20070228 NRK+1-->NRK+2 --- XX(:,:,NRK+2) will be used by slope limiter
      ALLOCATE ( dg%HB(dg%DOFH,S%MNE,dg%NRK+2),dg%e1(s%layers+5),dg%balance(s%layers+5) )
      ALLOCATE ( dg%MANN(dg%DOFH,S%MNE) ) !,dg%arrayfix(dg%DOFH,S%MNE,dg%NRK+2) )
      ALLOCATE ( dg%QY(dg%DOFH,S%MNE,dg%NRK+2), dg%QX(dg%DOFH,S%MNE,dg%NRK+2) )
      ALLOCATE ( dg%ZE(dg%DOFH,S%MNE,dg%NRK+2) )
#ifdef SEDLAY
      ALLOCATE ( dg%bed(dg%DOFH,S%MNE,dg%NRK+2,s%layers) )
#endif
      Allocate ( dg%iota(dg%DOFH,S%MNE,dg%NRK+2),dg%iota2(dg%DOFH,S%MNE,dg%NRK+2) ) 
      Allocate ( dg%dynP(dg%DOFH,S%MNE,dg%NRK+2) ) 
      Allocate ( dg%iotaa(dg%DOFH,S%MNE,dg%NRK+2),dg%iotaa2(dg%DOFH,S%MNE,dg%NRK+2) )
      Allocate ( dg%iotaa3(dg%DOFH,S%MNE,dg%NRK+2) )
!sb-20060711 For wet/dry
      ALLOCATE ( dg%ZE_MAX(S%MNE),dg%ZE_MIN(S%MNE),dg%DPE_MIN(S%MNE) )
      Allocate ( dg%iota_MAX(S%MNE),dg%iota_MIN(S%MNE),dg%iota2_MAX(S%MNE),dg%iota2_MIN(S%MNE) )
      Allocate ( dg%dynP_MAX(S%MNE),dg%dynP_MIN(S%MNE) )
      ALLOCATE ( dg%WATER_DEPTH(S%MNE,3), dg%WATER_DEPTH_OLD(S%MNE,3))
      ALLOCATE ( dg%ADVECTQX(S%MNE), dg%ADVECTQY(S%MNE))
      ALLOCATE ( dg%SOURCEQX(S%MNE),dg%SOURCEQY(S%MNE))
      ALLOCATE ( dg%MARK(S%MNE))
!em-2012 for sediment
      Allocate ( dg%bed_IN(s%layers),dg%bed_EX(s%layers),dg%bed_HAT(s%layers) )

!--
!sb-20070101
      ALLOCATE ( dg%LZ(dg%DOFH,2,2,S%MNE),dg%MZ(dg%DOFH,2,s%layers,S%MNE) )
      Allocate ( dg%HZ(dg%DOFH,2,2,S%MNE),dg%TZ(dg%DOFH,2,2,S%MNE) )
!--
      ALLOCATE ( dg%RHS_QX(dg%DOFH,S%MNE,dg%NRK), dg%RHS_QY(dg%DOFH,S%MNE,dg%NRK) )
      ALLOCATE ( dg%RHS_ZE(dg%DOFH,S%MNE,dg%NRK), dg%RHS_bed(dg%DOFH,S%MNE,dg%NRK,s%layers) )
      Allocate ( dg%RHS_iota(dg%DOFH,S%MNE,dg%NRK),dg%RHS_iota2(dg%DOFH,S%MNE,dg%NRK) )
      Allocate ( dg%RHS_dynP(dg%DOFH,S%MNE,dg%NRK) )
      Allocate ( dg%RHS_bed_IN(dg%dofh,s%layers), dg%RHS_bed_EX(dg%dofh,s%layers) )
      Allocate ( dg%bed_HAT_O(s%layers) )
      ALLOCATE ( dg%DRDX(S%MNE), dg%DSDX(S%MNE), dg%DRDY(S%MNE), dg%DSDY(S%MNE) )
      ALLOCATE ( dg%CORI_EL(S%MNE), dg%FRIC_EL(S%MNE) )
      ALLOCATE ( dg%PHI(dg%DOFH), dg%DPHIDZ1(dg%DOFH), dg%DPHIDZ2(dg%DOFH) )
      ALLOCATE ( dg%DOFS(S%MNE), dg%PCOUNT(S%MNE) )
      ALLOCATE ( dg%PDG(S%MNP) )
      RETURN
      END SUBROUTINE

!.....Set RK time scheme parameters array sizes

      SUBROUTINE ALLOC_RK(dg)
        type (dg_type) :: dg
      ALLOCATE( dg%ATVD(dg%NRK,dg%NRK), dg%BTVD(dg%NRK,dg%NRK), dg%CTVD(dg%NRK,dg%NRK)  )
      ALLOCATE( dg%DTVD(dg%NRK), dg%MAX_BOA_DT(dg%NRK) )
      Allocate( dg%RKC_T(0:dg%nrk),dg%RKC_U(0:dg%nrk),dg%RKC_Tprime(0:dg%nrk) )
      Allocate( dg%RKC_Tdprime(0:dg%nrk),dg%RKC_a(0:dg%nrk),dg%RKC_b(0:dg%nrk),dg%RKC_c(0:dg%nrk) ) 
      Allocate( dg%RKC_mu(0:dg%nrk),dg%RKC_tildemu(0:dg%nrk),dg%RKC_nu(0:dg%nrk),dg%RKC_gamma(0:dg%nrk) )
      RETURN
      END SUBROUTINE ALLOC_RK

!.....Set sizes for arrays used in orthobasis

      SUBROUTINE ALLOC_JACOBI(dg)
        type (dg_type) :: dg
      ALLOCATE ( dg%JACOBI(dg%ph+1,2*dg%ph+3,2,dg%NAGP(dg%ph)+1) )
      ALLOCATE ( dg%DXPHI2(dg%ph+1,dg%ph+1,dg%NAGP(dg%ph)+1),dg%DYPHI2(dg%ph+1,dg%ph+1,dg%NAGP(dg%ph)+1) )
      ALLOCATE ( dg%PHI2(dg%ph+1,dg%ph+1,dg%NAGP(dg%ph)+1) )
      ALLOCATE ( dg%PHI_CORNER1(dg%ph+1,dg%ph+1,3,dg%ph) )
      RETURN
      END SUBROUTINE
      
!.....Set sizes for arrays for area integrals
      
      SUBROUTINE ALLOC_AREA_GAUSS(s,dg)
        type (sizes_type) :: s
        type (dg_type) :: dg
      ALLOCATE ( dg%XAGP(dg%NAGP(dg%ph),dg%ph),dg%YAGP(dg%NAGP(dg%ph),dg%ph),dg%WAGP(dg%NAGP(dg%ph),dg%ph) ) 
      ALLOCATE ( dg%PHI_AREA(dg%DOFH,dg%NAGP(dg%ph)+1,dg%ph ) )
      ALLOCATE ( dg%DSPHI(dg%DOFH,dg%NAGP(dg%ph)+1,dg%ph),dg%DRPHI(dg%DOFH,dg%NAGP(dg%ph)+1,dg%ph) )
      ALLOCATE ( dg%PHI_CORNER(dg%DOFH,3,dg%ph),dg%PHI_MID(dg%DOFH,3,dg%ph) )
      ALLOCATE ( dg%PHI_CENTER(dg%DOFH,dg%DOFH) )
      ALLOCATE ( dg%PSI1(dg%NAGP(dg%ph),dg%ph),dg%PSI2(dg%NAGP(dg%ph),dg%ph),dg%PSI3(dg%NAGP(dg%ph),dg%ph) )
      ALLOCATE ( dg%BATH(dg%NAGP(dg%ph),S%MNE,dg%ph),dg%DBATHDX(dg%NAGP(dg%ph),S%MNE,dg%ph) )
      Allocate ( dg%DBATHDY(dg%NAGP(dg%ph),S%MNE,dg%ph) )
      ALLOCATE ( dg%SFAC_ELEM(dg%NAGP(dg%ph),S%MNE,dg%ph) )
      ALLOCATE ( dg%XFAC(dg%DOFH,dg%NAGP(dg%ph),S%MNE,dg%ph), dg%YFAC(dg%DOFH,dg%NAGP(dg%ph),S%MNE,dg%ph) )
      ALLOCATE ( dg%SRFAC(dg%DOFH,dg%NAGP(dg%ph),S%MNE,dg%ph) )
      RETURN
      END SUBROUTINE
      
!.....Set sizes for arrays for edge integrals

      SUBROUTINE ALLOC_EDGE_GAUSS(s,dg)
        type (sizes_type) :: s
        type (dg_type) :: dg
      ALLOCATE ( dg%XEGP(dg%NEGP(dg%ph),dg%ph), dg%WEGP(dg%NEGP(dg%ph),dg%ph) )
      ALLOCATE ( dg%PHI_EDGE(dg%DOFH,dg%NEGP(dg%ph)+1,3,dg%ph) )
      ALLOCATE ( dg%M_INV(dg%DOFH,dg%ph) )
      ALLOCATE ( dg%BATHED(dg%NEGP(dg%ph),3,S%MNE,dg%ph),dg%SFACED(dg%NEGP(dg%ph),3,S%MNE,dg%ph) )
      ALLOCATE ( dg%EDGEQ(dg%DOFH,dg%NEGP(dg%ph),3,dg%ph) )
      RETURN
      END SUBROUTINE

!.....Set sizes for the arrays for the slope limiter
!.....slopelim arrays

      SUBROUTINE ALLOC_SLOPELIM(s,dg)
        type (sizes_type) :: s
        type (dg_type) :: dg
      ALLOCATE ( dg%XBC(S%MNE), dg%YBC(S%MNE) )
      ALLOCATE ( dg%EL_NBORS(4,S%MNE) )
      ALLOCATE ( dg%SL3(3,S%MNE) )

!.....These are defined in prep_slopelim.dg%F

      Allocate ( dg%fact(0:dg%ph) ,dg%focal_neigh(S%MNE,3*S%MNEI),dg%focal_up(S%MNE),dg%bi(dg%dofh),dg%bj(dg%dofh) )

      Allocate ( dg%XBCb(S%MNE),dg%YBCb(S%MNE),dg%xi1(S%MNE,dg%NAGP(dg%ph)),dg%xi2(S%MNE,dg%NAGP(dg%ph)) )
      Allocate ( dg%xtransform(S%MNE,dg%NAGP(dg%ph)),dg%ytransform(S%MNE,dg%NAGP(dg%ph)) )
      Allocate ( dg%xi1BCb(S%MNE),dg%xi2BCb(S%MNE),dg%xi1vert(S%MNE,3) )
      Allocate ( dg%xi2vert(S%MNE,3),dg%xtransformv(S%MNE,3),dg%ytransformv(S%MNE,3) )
      Allocate ( dg%XBCv(S%MNE,S%MNE),dg%YBCv(S%MNE,S%MNE) )
      Allocate ( dg%xi1BCv(S%MNE,S%MNE),dg%xi2BCv(S%MNE,S%MNE),dg%Area_integral(S%MNE,0:dg%ph,0:dg%ph) )
      Allocate ( dg%f(S%MNE,dg%NAGP(dg%ph),0:dg%ph,0:dg%ph),dg%g0(S%MNE,dg%NAGP(dg%ph),0:dg%ph,0:dg%ph) )
      Allocate ( dg%fv(S%MNE,3,0:dg%ph,0:dg%ph),dg%g0v(S%MNE,3,0:dg%ph,0:dg%ph) )
      Allocate ( dg%varsigma0(S%MNE,dg%NAGP(dg%ph),0:dg%ph,0:dg%ph) ) 
      Allocate ( dg%varsigma0v(S%MNE,3,0:dg%ph,0:dg%ph) )
      Allocate ( dg%pmatrix(S%MNE,dg%dofh,dg%dofh), dg%var2sigmag(S%MNE,dg%NAGP(dg%ph),dg%dofh) )
      Allocate ( dg%Nmatrix(S%MNE,dg%dofh,dg%dofh,dg%dofh),dg%NmatrixInv(S%MNE,dg%dofh,dg%dofh,dg%dofh) )
      Allocate ( dg%deltx(S%MNE),dg%delty(S%MNE),dg%var2sigmav(S%MNE,3,dg%dofh))

!.....These (below) are defined in slopelimiter.dg%F (slopelimiter4)
      
      Allocate ( dg%ZEmin(S%MNP,dg%dofh),dg%ZEmax(S%MNP,dg%dofh),dg%QXmin(S%MNP,dg%dofh) )
      Allocate ( dg%QXmax(S%MNP,dg%dofh),dg%QYmin(S%MNP,dg%dofh),dg%QYmax(S%MNP,dg%dofh) )
      Allocate ( dg%iotamin(S%MNP,dg%dofh),dg%iotamax(S%MNP,dg%dofh) )
      Allocate ( dg%iota2min(S%MNP,dg%dofh),dg%iota2max(S%MNP,dg%dofh) )

#ifdef SLOPEALL
      Allocate ( dg%ZEtaylor(S%MNE,dg%dofh,1),dg%QXtaylor(S%MNE,dg%dofh,1) )
      Allocate ( dg%iotataylor(S%MNE,dg%dofh,1),dg%iota2taylor(S%MNE,dg%dofh,1) 
      Allocate ( dg%ZEtaylorvert(S%MNE,dg%dofh,3),dg%QXtaylorvert(S%MNE,dg%dofh,3) )
      Allocate ( dg%QYtaylorvert(S%MNE,dg%dofh,3),dg%iotataylorvert(S%MNE,dg%dofh,3) )
      Allocate ( dg%iota2taylorvert(S%MNE,dg%dofh,3), dg%QYtaylor(S%MNE,dg%dofh,1) )
      Allocate ( dg%alphaZE0(S%MNE,dg%dofh,3),dg%alphaQX0(S%MNE,dg%dofh,3) )
      Allocate ( dg%alphaQY0(S%MNE,dg%dofh,3),dg%alphaiota0(S%MNE,dg%dofh,3) )
      Allocate ( dg%alphaiota20(S%MNE,dg%dofh,3) )
      Allocate ( dg%alphaZE(S%MNE,dg%dofh),dg%alphaQX(S%MNE,dg%dofh),dg%alphaQY(S%MNE,dg%dofh) )
      Allocate ( dg%alphaiota(S%MNE,dg%dofh),dg%alphaiota2(S%MNE,dg%dofh) )
      Allocate ( dg%alphaZEm(S%MNE,dg%dofh),dg%alphaQXm(S%MNE,dg%dofh),dg%alphaQYm(S%MNE,dg%dofh) )
      Allocate ( dg%alphaiotam(S%MNE,dg%dofh),dg%alphaiota2m(S%MNE,dg%dofh) )
      Allocate ( dg%alphaZE_max(S%MNE,dg%dofh),dg%alphaQX_max(S%MNE,dg%dofh) )
      Allocate ( dg%alphaQY_max(S%MNE,dg%dofh) )
      Allocate ( dg%alphaiota_max(S%MNE,dg%dofh),dg%alphaiota2_max(S%MNE,dg%dofh) )
      Allocate ( dg%limitZE(S%MNE,dg%dofh),dg%limitQX(S%MNE,dg%dofh),dg%limitQY(S%MNE,dg%dofh) )
      Allocate ( dg%limitiota(S%MNE,dg%dofh),dg%limitiota2(S%MNE,dg%dofh) )
      Allocate ( dg%ZEconst(S%MNE,dg%dofh),dg%QXconst(S%MNE,dg%dofh),dg%QYconst(S%MNE,dg%dofh) )
      Allocate ( dg%iotaconst(S%MNE,dg%dofh),dg%iota2const(S%MNE,dg%dofh) )
#endif

      RETURN
      END SUBROUTINE

!sb...Set sizes for arrays for wetting and drying
      SUBROUTINE ALLOC_DG_WETDRY(s,dg)
        type (sizes_type) :: s
        type (dg_type) :: dg
      ALLOCATE ( dg%WDFLG(S%MNE), dg%DOFW(S%MNE) )
      ALLOCATE ( dg%EL_UPDATED(S%MNE) )
      ALLOCATE ( dg%WDFLG_TMP(S%MNE) )
      ALLOCATE ( dg%LEDGE_NVEC(3,3,S%MNE) )
      ALLOCATE ( dg%DP_VOL(S%MNE,dg%ph) )
      ALLOCATE ( dg%PHI_INTEGRATED(dg%DOFH,dg%ph) )
      ALLOCATE ( dg%PHI_CHECK(dg%DOFH,dg%NCHECK(dg%ph),dg%ph) )
      ALLOCATE ( dg%DP_NODE(dg%NCHECK(dg%ph),S%MNE,dg%ph) )
      ALLOCATE ( dg%PSI_CHECK(3,12*3) )
      RETURN
      END SUBROUTINE
      

      SUBROUTINE ALLOC_STAE(dg,L)
        type (dg_type) :: dg
        integer :: L
      ALLOCATE ( dg%PHI_STAE(dg%DOFH,L) )
      RETURN
      END SUBROUTINE
      
      SUBROUTINE ALLOC_STAV(dg,L)
        type (dg_type) :: dg
        integer :: L
      ALLOCATE ( dg%PHI_STAV(dg%DOFH,L) )
      RETURN
      END SUBROUTINE
      

      END MODULE DG


