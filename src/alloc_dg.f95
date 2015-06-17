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
!.....Set edge array sizes

      SUBROUTINE ALLOC_EDGES0()
      ALLOCATE ( IBHT(3*MNE), EBHT(3*MNE) )
      ALLOCATE ( EBCFSP(3*MNE), IBCFSP(3*MNE), IBCFSB(3*MNE) )
      ALLOCATE ( BACKNODES(2,3*MNE) )
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
      ALLOCATE ( EFA_DG(MNBFR,NEEDS+2,2), EMO_DG(MNBFR,NEEDS+2,2) )
      ALLOCATE ( UFA_DG(MNBFR,NEEDS+2,2), UMO_DG(MNBFR,NEEDS+2,2) )
      ALLOCATE ( VFA_DG(MNBFR,NEEDS+2,2), VMO_DG(MNBFR,NEEDS+2,2) )
      RETURN
      END SUBROUTINE
      
      SUBROUTINE ALLOC_DG2(MNFFR)
      ALLOCATE ( QNAM_DG(MNFFR,NFEDS,2), QNPH_DG(MNFFR,NFEDS,2) )
      RETURN
      END SUBROUTINE
      
      SUBROUTINE ALLOC_DG3(MNP)
      ALLOCATE ( QIB(MNP) )
      RETURN
      END SUBROUTINE

      SUBROUTINE ALLOC_DG4()
!sb-20070228 NRK+1-->NRK+2 --- XX(:,:,NRK+2) will be used by slope limiter
      ALLOCATE ( HB(DOFH,MNE,NRK+2),e1(layers+5),balance(layers+5) )
      ALLOCATE ( MANN(DOFH,MNE) ) !,arrayfix(DOFH,MNE,NRK+2) )
      ALLOCATE ( QY(DOFH,MNE,NRK+2), QX(DOFH,MNE,NRK+2) )
      ALLOCATE ( ZE(DOFH,MNE,NRK+2), bed(DOFH,MNE,NRK+2,layers) )
      Allocate ( iota(DOFH,MNE,NRK+2),iota2(DOFH,MNE,NRK+2) ) 
      Allocate ( dynP(DOFH,MNE,NRK+2) ) 
      Allocate ( iotaa(DOFH,MNE,NRK+2),iotaa2(DOFH,MNE,NRK+2) )
      Allocate ( iotaa3(DOFH,MNE,NRK+2) )
!sb-20060711 For wet/dry
      ALLOCATE ( ZE_MAX(MNE),ZE_MIN(MNE),DPE_MIN(MNE) )
      Allocate ( iota_MAX(MNE),iota_MIN(MNE),iota2_MAX(MNE),iota2_MIN(MNE) )
      Allocate ( dynP_MAX(MNE),dynP_MIN(MNE) )
      ALLOCATE ( WATER_DEPTH(MNE,3), WATER_DEPTH_OLD(MNE,3))
      ALLOCATE ( ADVECTQX(MNE), ADVECTQY(MNE))
      ALLOCATE ( SOURCEQX(MNE),SOURCEQY(MNE))
      ALLOCATE ( MARK(MNE))
!em-2012 for sediment
      Allocate ( bed_IN(layers),bed_EX(layers),bed_HAT(layers) )

!--
!sb-20070101
      ALLOCATE ( LZ(DOFH,2,2,MNE),MZ(DOFH,2,layers,MNE) )
      Allocate ( HZ(DOFH,2,2,MNE),TZ(DOFH,2,2,MNE) )
!--
      ALLOCATE ( RHS_QX(DOFH,MNE,NRK), RHS_QY(DOFH,MNE,NRK) )
      ALLOCATE ( RHS_ZE(DOFH,MNE,NRK), RHS_bed(DOFH,MNE,NRK,layers) )
      Allocate ( RHS_iota(DOFH,MNE,NRK),RHS_iota2(DOFH,MNE,NRK) )
      Allocate ( RHS_dynP(DOFH,MNE,NRK) )
      Allocate ( RHS_bed_IN(dofh,layers), RHS_bed_EX(dofh,layers) )
      Allocate ( bed_HAT_O(layers) )
      ALLOCATE ( DRDX(MNE), DSDX(MNE), DRDY(MNE), DSDY(MNE) )
      ALLOCATE ( CORI_EL(MNE), FRIC_EL(MNE) )
      ALLOCATE ( PHI(DOFH), DPHIDZ1(DOFH), DPHIDZ2(DOFH) )
      ALLOCATE ( DOFS(MNE), PCOUNT(MNE) )
      ALLOCATE ( PDG(MNP) )
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
      ALLOCATE ( JACOBI(dg%ph+1,2*dg%ph+3,2,NAGP(dg%ph)+1) )
      ALLOCATE ( DXPHI2(dg%ph+1,dg%ph+1,NAGP(dg%ph)+1),DYPHI2(dg%ph+1,dg%ph+1,NAGP(dg%ph)+1) )
      ALLOCATE ( PHI2(dg%ph+1,dg%ph+1,NAGP(dg%ph)+1) )
      ALLOCATE ( PHI_CORNER1(dg%ph+1,dg%ph+1,3,dg%ph) )
      RETURN
      END SUBROUTINE
      
!.....Set sizes for arrays for area integrals
      
      SUBROUTINE ALLOC_AREA_GAUSS()
      ALLOCATE ( XAGP(NAGP(dg%ph),dg%ph),YAGP(NAGP(dg%ph),dg%ph),WAGP(NAGP(dg%ph),dg%ph) ) 
      ALLOCATE ( PHI_AREA(DOFH,NAGP(dg%ph)+1,dg%ph ) )
      ALLOCATE ( DSPHI(DOFH,NAGP(dg%ph)+1,dg%ph),DRPHI(DOFH,NAGP(dg%ph)+1,dg%ph) )
      ALLOCATE ( PHI_CORNER(DOFH,3,dg%ph),PHI_MID(DOFH,3,dg%ph) )
      ALLOCATE ( PHI_CENTER(DOFH,DOFH) )
      ALLOCATE ( PSI1(NAGP(dg%ph),dg%ph),PSI2(NAGP(dg%ph),dg%ph),PSI3(NAGP(dg%ph),dg%ph) )
      ALLOCATE ( BATH(NAGP(dg%ph),MNE,dg%ph),DBATHDX(NAGP(dg%ph),MNE,dg%ph) )
      Allocate ( DBATHDY(NAGP(dg%ph),MNE,dg%ph) )
      ALLOCATE ( SFAC_ELEM(NAGP(dg%ph),MNE,dg%ph) )
      ALLOCATE ( XFAC(DOFH,NAGP(dg%ph),MNE,dg%ph), YFAC(DOFH,NAGP(dg%ph),MNE,dg%ph) )
      ALLOCATE ( SRFAC(DOFH,NAGP(dg%ph),MNE,dg%ph) )
      RETURN
      END SUBROUTINE
      
!.....Set sizes for arrays for edge integrals

      SUBROUTINE ALLOC_EDGE_GAUSS()
      ALLOCATE ( XEGP(NEGP(dg%ph),dg%ph), WEGP(NEGP(dg%ph),dg%ph) )
      ALLOCATE ( PHI_EDGE(DOFH,NEGP(dg%ph)+1,3,dg%ph) )
      ALLOCATE ( M_INV(DOFH,dg%ph) )
      ALLOCATE ( BATHED(NEGP(dg%ph),3,MNE,dg%ph),SFACED(NEGP(dg%ph),3,MNE,dg%ph) )
      ALLOCATE ( EDGEQ(DOFH,NEGP(dg%ph),3,dg%ph) )
      RETURN
      END SUBROUTINE

!.....Set sizes for the arrays for the slope limiter
!.....slopelim arrays

      SUBROUTINE ALLOC_SLOPELIM()
      ALLOCATE ( XBC(MNE), YBC(MNE) )
      ALLOCATE ( EL_NBORS(4,MNE) )
      ALLOCATE ( SL3(3,MNE) )

!.....These are defined in prep_slopelim.F

      Allocate ( fact(0:dg%ph) ,focal_neigh(MNE,3*MNEI),focal_up(MNE),bi(dofh),bj(dofh) )

      Allocate ( XBCb(MNE),YBCb(MNE),xi1(MNE,NAGP(dg%ph)),xi2(MNE,NAGP(dg%ph)) )
      Allocate ( xtransform(MNE,NAGP(dg%ph)),ytransform(MNE,NAGP(dg%ph)) )
      Allocate ( xi1BCb(MNE),xi2BCb(MNE),xi1vert(MNE,3) )
      Allocate ( xi2vert(MNE,3),xtransformv(MNE,3),ytransformv(MNE,3) )
      Allocate ( XBCv(MNE,MNE),YBCv(MNE,MNE) )
      Allocate ( xi1BCv(MNE,MNE),xi2BCv(MNE,MNE),Area_integral(MNE,0:dg%ph,0:dg%ph) )
      Allocate ( f(MNE,NAGP(dg%ph),0:dg%ph,0:dg%ph),g0(MNE,NAGP(dg%ph),0:dg%ph,0:dg%ph) )
      Allocate ( fv(MNE,3,0:dg%ph,0:dg%ph),g0v(MNE,3,0:dg%ph,0:dg%ph) )
      Allocate ( varsigma0(MNE,NAGP(dg%ph),0:dg%ph,0:dg%ph) ) 
      Allocate ( varsigma0v(MNE,3,0:dg%ph,0:dg%ph) )
      Allocate ( pmatrix(MNE,dofh,dofh), var2sigmag(MNE,NAGP(dg%ph),dofh) )
      Allocate ( Nmatrix(MNE,dofh,dofh,dofh),NmatrixInv(MNE,dofh,dofh,dofh) )
      Allocate ( deltx(MNE),delty(MNE),var2sigmav(MNE,3,dofh))

!.....These (below) are defined in slopelimiter.F (slopelimiter4)
      
      Allocate ( ZEmin(MNP,dofh),ZEmax(MNP,dofh),QXmin(MNP,dofh) )
      Allocate ( QXmax(MNP,dofh),QYmin(MNP,dofh),QYmax(MNP,dofh) )
      Allocate ( iotamin(MNP,dofh),iotamax(MNP,dofh) )
      Allocate ( iota2min(MNP,dofh),iota2max(MNP,dofh) )

#ifdef SLOPEALL
      Allocate ( ZEtaylor(MNE,dofh,1),QXtaylor(MNE,dofh,1) )
      Allocate ( iotataylor(MNE,dofh,1),iota2taylor(MNE,dofh,1) 
      Allocate ( ZEtaylorvert(MNE,dofh,3),QXtaylorvert(MNE,dofh,3) )
      Allocate ( QYtaylorvert(MNE,dofh,3),iotataylorvert(MNE,dofh,3) )
      Allocate ( iota2taylorvert(MNE,dofh,3), QYtaylor(MNE,dofh,1) )
      Allocate ( alphaZE0(MNE,dofh,3),alphaQX0(MNE,dofh,3) )
      Allocate ( alphaQY0(MNE,dofh,3),alphaiota0(MNE,dofh,3) )
      Allocate ( alphaiota20(MNE,dofh,3) )
      Allocate ( alphaZE(MNE,dofh),alphaQX(MNE,dofh),alphaQY(MNE,dofh) )
      Allocate ( alphaiota(MNE,dofh),alphaiota2(MNE,dofh) )
      Allocate ( alphaZEm(MNE,dofh),alphaQXm(MNE,dofh),alphaQYm(MNE,dofh) )
      Allocate ( alphaiotam(MNE,dofh),alphaiota2m(MNE,dofh) )
      Allocate ( alphaZE_max(MNE,dofh),alphaQX_max(MNE,dofh) )
      Allocate ( alphaQY_max(MNE,dofh) )
      Allocate ( alphaiota_max(MNE,dofh),alphaiota2_max(MNE,dofh) )
      Allocate ( limitZE(MNE,dofh),limitQX(MNE,dofh),limitQY(MNE,dofh) )
      Allocate ( limitiota(MNE,dofh),limitiota2(MNE,dofh) )
      Allocate ( ZEconst(MNE,dofh),QXconst(MNE,dofh),QYconst(MNE,dofh) )
      Allocate ( iotaconst(MNE,dofh),iota2const(MNE,dofh) )
#endif

      RETURN
      END SUBROUTINE

!sb...Set sizes for arrays for wetting and drying
      SUBROUTINE ALLOC_DG_WETDRY()
      ALLOCATE ( WDFLG(MNE), DOFW(MNE) )
      ALLOCATE ( EL_UPDATED(MNE) )
      ALLOCATE ( WDFLG_TMP(MNE) )
      ALLOCATE ( LEDGE_NVEC(3,3,MNE) )
      ALLOCATE ( DP_VOL(MNE,dg%ph) )
      ALLOCATE ( PHI_INTEGRATED(DOFH,dg%ph) )
      ALLOCATE ( PHI_CHECK(DOFH,NCHECK(dg%ph),dg%ph) )
      ALLOCATE ( DP_NODE(NCHECK(dg%ph),MNE,dg%ph) )
      ALLOCATE ( PSI_CHECK(3,12*3) )
      RETURN
      END SUBROUTINE
      

      SUBROUTINE ALLOC_STAE(L)
      ALLOCATE ( PHI_STAE(DOFH,L) )
      RETURN
      END SUBROUTINE
      
      SUBROUTINE ALLOC_STAV(L)
      ALLOCATE ( PHI_STAV(DOFH,L) )
      RETURN
      END SUBROUTINE
