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

#ifdef HPX
      INTEGER,ALLOCATABLE :: NELEMLOC(:)
      INTEGER,ALLOCATABLE :: IBELONGTO(:)
      LOGICAL,ALLOCATABLE :: RESELEM(:)
      INTEGER :: NEIGHPROC_R, NEIGHPROC_S
      INTEGER,ALLOCATABLE :: IPROC_R(:), NELEMRECV(:), IRECVLOC(:,:)
      INTEGER,ALLOCATABLE :: IPROC_S(:), NELEMSEND(:), ISENDLOC(:,:)
      INTEGER,ALLOCATABLE :: INDEX(:)
      INTEGER :: RDIM
#endif

   end type dg_type
    
!**********************END OF DATA DECLARATIONS ***********************

      CONTAINS
      
!.....Set edge array sizes

      SUBROUTINE ALLOC_EDGES0(s,dg_here)
        type (sizes_type) :: s
        type (dg_type) :: dg_here
      ALLOCATE ( dg_here%IBHT(3*s%MNE), dg_here%EBHT(3*s%MNE) )
      ALLOCATE ( dg_here%EBCFSP(3*s%MNE), dg_here%IBCFSP(3*s%MNE), dg_here%IBCFSB(3*s%MNE) )
      ALLOCATE ( dg_here%BACKNODES(2,3*s%MNE) )
      RETURN
      END SUBROUTINE

      SUBROUTINE ALLOC_EDGES1(dg_here)
        type (dg_type) :: dg_here
      ALLOCATE ( dg_here%NEDNO(2,dg_here%MNED), dg_here%NEDEL(2,dg_here%MNED), dg_here%NEDSD(2,dg_here%MNED) )
      ALLOCATE ( dg_here%NIBSEGN(2,dg_here%MNED) )
      ALLOCATE ( dg_here%NEBSEGN(dg_here%MNED) )
      ALLOCATE ( dg_here%NIEDN(dg_here%MNED), dg_here%NLEDN(dg_here%MNED), dg_here%NEEDN(dg_here%MNED) )
      ALLOCATE ( dg_here%NFEDN(dg_here%MNED),  dg_here%NREDN(dg_here%MNED), dg_here%NIBEDN(dg_here%MNED), dg_here%NEBEDN(dg_here%MNED) )
      ALLOCATE ( dg_here%NCOUNT(dg_here%MNED) )
      ALLOCATE ( dg_here%COSNX(dg_here%MNED), dg_here%SINNX(dg_here%MNED), dg_here%XLEN(dg_here%MNED) )
      ALLOCATE ( dg_here%Q_HAT(dg_here%MNED) )
      RETURN
      END SUBROUTINE
      
!.....Set DG SWE array sizes

      SUBROUTINE ALLOC_DG1(dg_here,MNBFR)
        type (dg_type) :: dg_here
        integer :: MNBFR
      ALLOCATE ( dg_here%EFA_DG(MNBFR,dg_here%NEEDS+2,2), dg_here%EMO_DG(MNBFR,dg_here%NEEDS+2,2) )
      ALLOCATE ( dg_here%UFA_DG(MNBFR,dg_here%NEEDS+2,2), dg_here%UMO_DG(MNBFR,dg_here%NEEDS+2,2) )
      ALLOCATE ( dg_here%VFA_DG(MNBFR,dg_here%NEEDS+2,2), dg_here%VMO_DG(MNBFR,dg_here%NEEDS+2,2) )
      RETURN
      END SUBROUTINE
      
      SUBROUTINE ALLOC_DG2(dg_here,MNFFR)
        type (dg_type) :: dg_here
        integer :: MNFFR
      ALLOCATE ( dg_here%QNAM_DG(MNFFR,dg_here%NFEDS,2), dg_here%QNPH_DG(MNFFR,dg_here%NFEDS,2) )
      RETURN
      END SUBROUTINE
      
      SUBROUTINE ALLOC_DG3(dg_here,MNP)
        type (dg_type) :: dg_here
        integer :: MNP
      ALLOCATE ( dg_here%QIB(MNP) )
      RETURN
      END SUBROUTINE

      SUBROUTINE ALLOC_DG4(s,dg_here)
        type (sizes_type) :: s
        type (dg_type) :: dg_here
!sb-20070228 NRK+1-->NRK+2 --- XX(:,:,NRK+2) will be used by slope limiter
      ALLOCATE ( dg_here%HB(dg_here%DOFH,S%MNE,dg_here%NRK+2),dg_here%e1(s%layers+5),dg_here%balance(s%layers+5) )
      ALLOCATE ( dg_here%MANN(dg_here%DOFH,S%MNE) ) !,dg_here%arrayfix(dg_here%DOFH,S%MNE,dg_here%NRK+2) )
      ALLOCATE ( dg_here%QY(dg_here%DOFH,S%MNE,dg_here%NRK+2), dg_here%QX(dg_here%DOFH,S%MNE,dg_here%NRK+2) )
      ALLOCATE ( dg_here%ZE(dg_here%DOFH,S%MNE,dg_here%NRK+2) )
#ifdef SEDLAY
      ALLOCATE ( dg_here%bed(dg_here%DOFH,S%MNE,dg_here%NRK+2,s%layers) )
#endif
      Allocate ( dg_here%iota(dg_here%DOFH,S%MNE,dg_here%NRK+2),dg_here%iota2(dg_here%DOFH,S%MNE,dg_here%NRK+2) ) 
      Allocate ( dg_here%dynP(dg_here%DOFH,S%MNE,dg_here%NRK+2) ) 
      Allocate ( dg_here%iotaa(dg_here%DOFH,S%MNE,dg_here%NRK+2),dg_here%iotaa2(dg_here%DOFH,S%MNE,dg_here%NRK+2) )
      Allocate ( dg_here%iotaa3(dg_here%DOFH,S%MNE,dg_here%NRK+2) )
!sb-20060711 For wet/dry
      ALLOCATE ( dg_here%ZE_MAX(S%MNE),dg_here%ZE_MIN(S%MNE),dg_here%DPE_MIN(S%MNE) )
      Allocate ( dg_here%iota_MAX(S%MNE),dg_here%iota_MIN(S%MNE),dg_here%iota2_MAX(S%MNE),dg_here%iota2_MIN(S%MNE) )
      Allocate ( dg_here%dynP_MAX(S%MNE),dg_here%dynP_MIN(S%MNE) )
      ALLOCATE ( dg_here%WATER_DEPTH(S%MNE,3), dg_here%WATER_DEPTH_OLD(S%MNE,3))
      ALLOCATE ( dg_here%ADVECTQX(S%MNE), dg_here%ADVECTQY(S%MNE))
      ALLOCATE ( dg_here%SOURCEQX(S%MNE),dg_here%SOURCEQY(S%MNE))
      ALLOCATE ( dg_here%MARK(S%MNE))
!em-2012 for sediment
      Allocate ( dg_here%bed_IN(s%layers),dg_here%bed_EX(s%layers),dg_here%bed_HAT(s%layers) )

!--
!sb-20070101
      ALLOCATE ( dg_here%LZ(dg_here%DOFH,2,2,S%MNE),dg_here%MZ(dg_here%DOFH,2,s%layers,S%MNE) )
      Allocate ( dg_here%HZ(dg_here%DOFH,2,2,S%MNE),dg_here%TZ(dg_here%DOFH,2,2,S%MNE) )
!--
      ALLOCATE ( dg_here%RHS_QX(dg_here%DOFH,S%MNE,dg_here%NRK), dg_here%RHS_QY(dg_here%DOFH,S%MNE,dg_here%NRK) )
      ALLOCATE ( dg_here%RHS_ZE(dg_here%DOFH,S%MNE,dg_here%NRK), dg_here%RHS_bed(dg_here%DOFH,S%MNE,dg_here%NRK,s%layers) )
      Allocate ( dg_here%RHS_iota(dg_here%DOFH,S%MNE,dg_here%NRK),dg_here%RHS_iota2(dg_here%DOFH,S%MNE,dg_here%NRK) )
      Allocate ( dg_here%RHS_dynP(dg_here%DOFH,S%MNE,dg_here%NRK) )
      Allocate ( dg_here%RHS_bed_IN(dg_here%dofh,s%layers), dg_here%RHS_bed_EX(dg_here%dofh,s%layers) )
      Allocate ( dg_here%bed_HAT_O(s%layers) )
      ALLOCATE ( dg_here%DRDX(S%MNE), dg_here%DSDX(S%MNE), dg_here%DRDY(S%MNE), dg_here%DSDY(S%MNE) )
      ALLOCATE ( dg_here%CORI_EL(S%MNE), dg_here%FRIC_EL(S%MNE) )
      ALLOCATE ( dg_here%PHI(dg_here%DOFH), dg_here%DPHIDZ1(dg_here%DOFH), dg_here%DPHIDZ2(dg_here%DOFH) )
      ALLOCATE ( dg_here%DOFS(S%MNE), dg_here%PCOUNT(S%MNE) )
      ALLOCATE ( dg_here%PDG(S%MNP) )
      RETURN
      END SUBROUTINE

!.....Set RK time scheme parameters array sizes

      SUBROUTINE ALLOC_RK(dg_here)
        type (dg_type) :: dg_here
      ALLOCATE( dg_here%ATVD(dg_here%NRK,dg_here%NRK), dg_here%BTVD(dg_here%NRK,dg_here%NRK), dg_here%CTVD(dg_here%NRK,dg_here%NRK)  )
      ALLOCATE( dg_here%DTVD(dg_here%NRK), dg_here%MAX_BOA_DT(dg_here%NRK) )
      Allocate( dg_here%RKC_T(0:dg_here%nrk),dg_here%RKC_U(0:dg_here%nrk),dg_here%RKC_Tprime(0:dg_here%nrk) )
      Allocate( dg_here%RKC_Tdprime(0:dg_here%nrk),dg_here%RKC_a(0:dg_here%nrk),dg_here%RKC_b(0:dg_here%nrk),dg_here%RKC_c(0:dg_here%nrk) ) 
      Allocate( dg_here%RKC_mu(0:dg_here%nrk),dg_here%RKC_tildemu(0:dg_here%nrk),dg_here%RKC_nu(0:dg_here%nrk),dg_here%RKC_gamma(0:dg_here%nrk) )
      RETURN
      END SUBROUTINE ALLOC_RK

!.....Set sizes for arrays used in orthobasis

      SUBROUTINE ALLOC_JACOBI(dg_here)
        type (dg_type) :: dg_here
      ALLOCATE ( dg_here%JACOBI(dg_here%ph+1,2*dg_here%ph+3,2,dg_here%NAGP(dg_here%ph)+1) )
      ALLOCATE ( dg_here%DXPHI2(dg_here%ph+1,dg_here%ph+1,dg_here%NAGP(dg_here%ph)+1),dg_here%DYPHI2(dg_here%ph+1,dg_here%ph+1,dg_here%NAGP(dg_here%ph)+1) )
      ALLOCATE ( dg_here%PHI2(dg_here%ph+1,dg_here%ph+1,dg_here%NAGP(dg_here%ph)+1) )
      ALLOCATE ( dg_here%PHI_CORNER1(dg_here%ph+1,dg_here%ph+1,3,dg_here%ph) )
      RETURN
      END SUBROUTINE
      
!.....Set sizes for arrays for area integrals
      
      SUBROUTINE ALLOC_AREA_GAUSS(s,dg_here)
        type (sizes_type) :: s
        type (dg_type) :: dg_here
      ALLOCATE ( dg_here%XAGP(dg_here%NAGP(dg_here%ph),dg_here%ph),dg_here%YAGP(dg_here%NAGP(dg_here%ph),dg_here%ph),dg_here%WAGP(dg_here%NAGP(dg_here%ph),dg_here%ph) ) 
      ALLOCATE ( dg_here%PHI_AREA(dg_here%DOFH,dg_here%NAGP(dg_here%ph)+1,dg_here%ph ) )
      ALLOCATE ( dg_here%DSPHI(dg_here%DOFH,dg_here%NAGP(dg_here%ph)+1,dg_here%ph),dg_here%DRPHI(dg_here%DOFH,dg_here%NAGP(dg_here%ph)+1,dg_here%ph) )
      ALLOCATE ( dg_here%PHI_CORNER(dg_here%DOFH,3,dg_here%ph),dg_here%PHI_MID(dg_here%DOFH,3,dg_here%ph) )
      ALLOCATE ( dg_here%PHI_CENTER(dg_here%DOFH,dg_here%DOFH) )
      ALLOCATE ( dg_here%PSI1(dg_here%NAGP(dg_here%ph),dg_here%ph),dg_here%PSI2(dg_here%NAGP(dg_here%ph),dg_here%ph),dg_here%PSI3(dg_here%NAGP(dg_here%ph),dg_here%ph) )
      ALLOCATE ( dg_here%BATH(dg_here%NAGP(dg_here%ph),S%MNE,dg_here%ph),dg_here%DBATHDX(dg_here%NAGP(dg_here%ph),S%MNE,dg_here%ph) )
      Allocate ( dg_here%DBATHDY(dg_here%NAGP(dg_here%ph),S%MNE,dg_here%ph) )
      ALLOCATE ( dg_here%SFAC_ELEM(dg_here%NAGP(dg_here%ph),S%MNE,dg_here%ph) )
      ALLOCATE ( dg_here%XFAC(dg_here%DOFH,dg_here%NAGP(dg_here%ph),S%MNE,dg_here%ph), dg_here%YFAC(dg_here%DOFH,dg_here%NAGP(dg_here%ph),S%MNE,dg_here%ph) )
      ALLOCATE ( dg_here%SRFAC(dg_here%DOFH,dg_here%NAGP(dg_here%ph),S%MNE,dg_here%ph) )
      RETURN
      END SUBROUTINE
      
!.....Set sizes for arrays for edge integrals

      SUBROUTINE ALLOC_EDGE_GAUSS(s,dg_here)
        type (sizes_type) :: s
        type (dg_type) :: dg_here
      ALLOCATE ( dg_here%XEGP(dg_here%NEGP(dg_here%ph),dg_here%ph), dg_here%WEGP(dg_here%NEGP(dg_here%ph),dg_here%ph) )
      ALLOCATE ( dg_here%PHI_EDGE(dg_here%DOFH,dg_here%NEGP(dg_here%ph)+1,3,dg_here%ph) )
      ALLOCATE ( dg_here%M_INV(dg_here%DOFH,dg_here%ph) )
      ALLOCATE ( dg_here%BATHED(dg_here%NEGP(dg_here%ph),3,S%MNE,dg_here%ph),dg_here%SFACED(dg_here%NEGP(dg_here%ph),3,S%MNE,dg_here%ph) )
      ALLOCATE ( dg_here%EDGEQ(dg_here%DOFH,dg_here%NEGP(dg_here%ph),3,dg_here%ph) )
      RETURN
      END SUBROUTINE

!.....Set sizes for the arrays for the slope limiter
!.....slopelim arrays

      SUBROUTINE ALLOC_SLOPELIM(s,dg_here)
        type (sizes_type) :: s
        type (dg_type) :: dg_here
      ALLOCATE ( dg_here%XBC(S%MNE), dg_here%YBC(S%MNE) )
      ALLOCATE ( dg_here%EL_NBORS(4,S%MNE) )
      ALLOCATE ( dg_here%SL3(3,S%MNE) )

!.....These are defined in prep_slopelim.dg_here%F

      Allocate ( dg_here%fact(0:dg_here%ph) ,dg_here%focal_neigh(S%MNE,3*S%MNEI),dg_here%focal_up(S%MNE),dg_here%bi(dg_here%dofh),dg_here%bj(dg_here%dofh) )

      Allocate ( dg_here%XBCb(S%MNE),dg_here%YBCb(S%MNE),dg_here%xi1(S%MNE,dg_here%NAGP(dg_here%ph)),dg_here%xi2(S%MNE,dg_here%NAGP(dg_here%ph)) )
      Allocate ( dg_here%xtransform(S%MNE,dg_here%NAGP(dg_here%ph)),dg_here%ytransform(S%MNE,dg_here%NAGP(dg_here%ph)) )
      Allocate ( dg_here%xi1BCb(S%MNE),dg_here%xi2BCb(S%MNE),dg_here%xi1vert(S%MNE,3) )
      Allocate ( dg_here%xi2vert(S%MNE,3),dg_here%xtransformv(S%MNE,3),dg_here%ytransformv(S%MNE,3) )
      Allocate ( dg_here%XBCv(S%MNE,S%MNE),dg_here%YBCv(S%MNE,S%MNE) )
      Allocate ( dg_here%xi1BCv(S%MNE,S%MNE),dg_here%xi2BCv(S%MNE,S%MNE),dg_here%Area_integral(S%MNE,0:dg_here%ph,0:dg_here%ph) )
      Allocate ( dg_here%f(S%MNE,dg_here%NAGP(dg_here%ph),0:dg_here%ph,0:dg_here%ph),dg_here%g0(S%MNE,dg_here%NAGP(dg_here%ph),0:dg_here%ph,0:dg_here%ph) )
      Allocate ( dg_here%fv(S%MNE,3,0:dg_here%ph,0:dg_here%ph),dg_here%g0v(S%MNE,3,0:dg_here%ph,0:dg_here%ph) )
      Allocate ( dg_here%varsigma0(S%MNE,dg_here%NAGP(dg_here%ph),0:dg_here%ph,0:dg_here%ph) ) 
      Allocate ( dg_here%varsigma0v(S%MNE,3,0:dg_here%ph,0:dg_here%ph) )
      Allocate ( dg_here%pmatrix(S%MNE,dg_here%dofh,dg_here%dofh), dg_here%var2sigmag(S%MNE,dg_here%NAGP(dg_here%ph),dg_here%dofh) )
      Allocate ( dg_here%Nmatrix(S%MNE,dg_here%dofh,dg_here%dofh,dg_here%dofh),dg_here%NmatrixInv(S%MNE,dg_here%dofh,dg_here%dofh,dg_here%dofh) )
      Allocate ( dg_here%deltx(S%MNE),dg_here%delty(S%MNE),dg_here%var2sigmav(S%MNE,3,dg_here%dofh))

!.....These (below) are defined in slopelimiter.dg_here%F (slopelimiter4)
      
      Allocate ( dg_here%ZEmin(S%MNP,dg_here%dofh),dg_here%ZEmax(S%MNP,dg_here%dofh),dg_here%QXmin(S%MNP,dg_here%dofh) )
      Allocate ( dg_here%QXmax(S%MNP,dg_here%dofh),dg_here%QYmin(S%MNP,dg_here%dofh),dg_here%QYmax(S%MNP,dg_here%dofh) )
      Allocate ( dg_here%iotamin(S%MNP,dg_here%dofh),dg_here%iotamax(S%MNP,dg_here%dofh) )
      Allocate ( dg_here%iota2min(S%MNP,dg_here%dofh),dg_here%iota2max(S%MNP,dg_here%dofh) )

#ifdef SLOPEALL
      Allocate ( dg_here%ZEtaylor(S%MNE,dg_here%dofh,1),dg_here%QXtaylor(S%MNE,dg_here%dofh,1) )
      Allocate ( dg_here%iotataylor(S%MNE,dg_here%dofh,1),dg_here%iota2taylor(S%MNE,dg_here%dofh,1) 
      Allocate ( dg_here%ZEtaylorvert(S%MNE,dg_here%dofh,3),dg_here%QXtaylorvert(S%MNE,dg_here%dofh,3) )
      Allocate ( dg_here%QYtaylorvert(S%MNE,dg_here%dofh,3),dg_here%iotataylorvert(S%MNE,dg_here%dofh,3) )
      Allocate ( dg_here%iota2taylorvert(S%MNE,dg_here%dofh,3), dg_here%QYtaylor(S%MNE,dg_here%dofh,1) )
      Allocate ( dg_here%alphaZE0(S%MNE,dg_here%dofh,3),dg_here%alphaQX0(S%MNE,dg_here%dofh,3) )
      Allocate ( dg_here%alphaQY0(S%MNE,dg_here%dofh,3),dg_here%alphaiota0(S%MNE,dg_here%dofh,3) )
      Allocate ( dg_here%alphaiota20(S%MNE,dg_here%dofh,3) )
      Allocate ( dg_here%alphaZE(S%MNE,dg_here%dofh),dg_here%alphaQX(S%MNE,dg_here%dofh),dg_here%alphaQY(S%MNE,dg_here%dofh) )
      Allocate ( dg_here%alphaiota(S%MNE,dg_here%dofh),dg_here%alphaiota2(S%MNE,dg_here%dofh) )
      Allocate ( dg_here%alphaZEm(S%MNE,dg_here%dofh),dg_here%alphaQXm(S%MNE,dg_here%dofh),dg_here%alphaQYm(S%MNE,dg_here%dofh) )
      Allocate ( dg_here%alphaiotam(S%MNE,dg_here%dofh),dg_here%alphaiota2m(S%MNE,dg_here%dofh) )
      Allocate ( dg_here%alphaZE_max(S%MNE,dg_here%dofh),dg_here%alphaQX_max(S%MNE,dg_here%dofh) )
      Allocate ( dg_here%alphaQY_max(S%MNE,dg_here%dofh) )
      Allocate ( dg_here%alphaiota_max(S%MNE,dg_here%dofh),dg_here%alphaiota2_max(S%MNE,dg_here%dofh) )
      Allocate ( dg_here%limitZE(S%MNE,dg_here%dofh),dg_here%limitQX(S%MNE,dg_here%dofh),dg_here%limitQY(S%MNE,dg_here%dofh) )
      Allocate ( dg_here%limitiota(S%MNE,dg_here%dofh),dg_here%limitiota2(S%MNE,dg_here%dofh) )
      Allocate ( dg_here%ZEconst(S%MNE,dg_here%dofh),dg_here%QXconst(S%MNE,dg_here%dofh),dg_here%QYconst(S%MNE,dg_here%dofh) )
      Allocate ( dg_here%iotaconst(S%MNE,dg_here%dofh),dg_here%iota2const(S%MNE,dg_here%dofh) )
#endif

      RETURN
      END SUBROUTINE

!sb...Set sizes for arrays for wetting and drying
      SUBROUTINE ALLOC_DG_WETDRY(s,dg_here)
        type (sizes_type) :: s
        type (dg_type) :: dg_here
      ALLOCATE ( dg_here%WDFLG(S%MNE), dg_here%DOFW(S%MNE) )
      ALLOCATE ( dg_here%EL_UPDATED(S%MNE) )
      ALLOCATE ( dg_here%WDFLG_TMP(S%MNE) )
      ALLOCATE ( dg_here%LEDGE_NVEC(3,3,S%MNE) )
      ALLOCATE ( dg_here%DP_VOL(S%MNE,dg_here%ph) )
      ALLOCATE ( dg_here%PHI_INTEGRATED(dg_here%DOFH,dg_here%ph) )
      ALLOCATE ( dg_here%PHI_CHECK(dg_here%DOFH,dg_here%NCHECK(dg_here%ph),dg_here%ph) )
      ALLOCATE ( dg_here%DP_NODE(dg_here%NCHECK(dg_here%ph),S%MNE,dg_here%ph) )
      ALLOCATE ( dg_here%PSI_CHECK(3,12*3) )
      RETURN
      END SUBROUTINE
      

      SUBROUTINE ALLOC_STAE(dg_here,L)
        type (dg_type) :: dg_here
        integer :: L
      ALLOCATE ( dg_here%PHI_STAE(dg_here%DOFH,L) )
      RETURN
      END SUBROUTINE
      
      SUBROUTINE ALLOC_STAV(dg_here,L)
        type (dg_type) :: dg_here
        integer :: L
      ALLOCATE ( dg_here%PHI_STAV(dg_here%DOFH,L) )
      RETURN
      END SUBROUTINE
      

      END MODULE DG


