!***********************************************************************
!     
!     SUBROUTINE INTERNAL_EDGE_HYDRO( )
!     
!     This subroutine does the following:
!     
!     1.  Calculates the values of the necessary variables at the edge
!     gauss points for INTERNAL edges
!     2.  Calls the appropriate subroutine to compute the flux at
!     these points.
!     3.  Calls the appropriate subroutine to compute the boundary
!     integrals.
!     
!     Written by Ethan Kubatko (06-11-2004)
!     
!     01-10-2011 - cem - adapted for p_enrichment and multicomponent
!     11-11-2011 - cem - adapted for layered sediment
!     
!-----------------------------------------------------------------------
!     
!     01-02-2007, sb, Modified for LDG
!     C***********************************************************************

      SUBROUTINE INTERNAL_EDGE_HYDRO(s,dg_here,global_here,IT)

!.....Use appropriate modules

      USE GLOBAL
      USE DG
!      USE NodalAttributes, ONLY :  ESLM
      use sizes

      IMPLICIT NONE

      type (sizes_type) :: s
      type (dg_type) :: dg_here
      type (global_type) :: global_here

!.....Declare local variables

      INTEGER L, LED_IN, LED_EX, GED, GP_IN, GP_EX,k,i,ll,IT
      REAL(SZ), PARAMETER :: ZERO = 1.D-12
      REAL(SZ) TX, TY, W_IN, W_EX
      REAL(SZ) LZ_XX_IN, LZ_XY_IN, LZ_YX_IN, LZ_YY_IN
      REAL(SZ) LZ_XX_EX, LZ_XY_EX, LZ_YX_EX, LZ_YY_EX
      REAL(SZ) HZ_X_EX, HZ_Y_EX, HZ_X_IN, HZ_Y_IN
      REAL(SZ) TZ_X_EX, TZ_Y_EX, TZ_X_IN, TZ_Y_IN
      REAL(SZ) EDFAC_IN, EDFAC_EX, DEN
      REAL(SZ) XLEN_EL_IN, XLEN_EL_EX,MZ_X_EX(s%layers),MZ_Y_EX(s%layers)
      REAL(SZ) MASS_EL_IN, MASS_EL_EX,MZ_X_IN(s%layers),MZ_Y_IN(s%layers)

      REAL(SZ) :: RHS_ZE_IN(MAX_DOFH), RHS_QX_IN(MAX_DOFH), RHS_QY_IN(MAX_DOFH)
      REAL(SZ) :: RHS_ZE_EX(MAX_DOFH), RHS_QX_EX(MAX_DOFH), RHS_QY_EX(MAX_DOFH)

#ifdef TRACE     
      REAL(SZ) :: RHS_iota_IN(MAX_DOFH), RHS_iota_EX(MAX_DOFH)
#endif
#ifdef CHEM
      REAL(SZ) :: RHS_iota_IN(MAX_DOFH), RHS_iota_EX(MAX_DOFH)
      REAL(SZ) :: RHS_iota2_IN(MAX_DOFH), RHS_iota2_EX(MAX_DOFH)
#endif
#ifdef dynp     
      REAL(SZ) :: RHS_dynP_IN(MAX_DOFH), RHS_dynP_EX(MAX_DOFH)
#endif

      REAL(SZ) ARK, BRK
      REAL(SZ) MAX_BOA          ! Maximum of beta_il/alpha_il for all l
      REAL(SZ) NLEQG_TMP, G_TMP
      REAL(SZ) F_HAT_O, G_HAT_O, H_HAT_O,i_hat_o,j_hat_o
      REAL(SZ) G_HAT_IN, H_HAT_IN
      REAL(SZ) G_HAT_EX, H_HAT_EX
      REAL(SZ) K_HAT_O
      Real(SZ) bed_HAT_IN(s%layers),bed_HAT_EX(s%layers)
!     

      if (dg_here%DOFH.gt.MAX_DOFH) then
         !FIXME: add exception 
      endif

      dg_here%test_el = 0
      DO 1000 L = 1,dg_here%NIEDS

!.......Retrieve the global and local global_here%edge number

         GED = dg_here%NIEDN(L)
         LED_IN = dg_here%NEDSD(1,GED)
         LED_EX = dg_here%NEDSD(2,GED)

!.......Retrieve the elements which share the global_here%edge

         global_here%EL_IN = dg_here%NEDEL(1,GED)
         global_here%EL_EX = dg_here%NEDEL(2,GED)

         dg_here%EL = global_here%EL_EX

         if (dg_here%DOFS(global_here%EL_EX).LT.dg_here%DOFS(global_here%EL_IN)) then
            dg_here%EL = global_here%EL_IN
         endif
         
         dg_here%pa = global_here%PDG_EL(dg_here%EL)

#ifdef P0
         if (dg_here%pa.eq.0) then
            dg_here%pa = 1
         endif
#endif

!.......If both elements on either global_here%side of global_here%edge are dry then skip

         IF((dg_here%WDFLG(global_here%EL_IN).EQ.0).AND.(dg_here%WDFLG(global_here%EL_EX).EQ.0)) GOTO 1000

         !dg_here%test_el = dg_here%test_el+1

!.....Save the current RHS

         DO K = 1,dg_here%DOFS(dg_here%EL)

            RHS_ZE_IN(K) = dg_here%RHS_ZE(K,global_here%EL_IN,dg_here%IRK)
            RHS_QX_IN(K) = dg_here%RHS_QX(K,global_here%EL_IN,dg_here%IRK)
            RHS_QY_IN(K) = dg_here%RHS_QY(K,global_here%EL_IN,dg_here%IRK)

            RHS_ZE_EX(K) = dg_here%RHS_ZE(K,global_here%EL_EX,dg_here%IRK)
            RHS_QX_EX(K) = dg_here%RHS_QX(K,global_here%EL_EX,dg_here%IRK)
            RHS_QY_EX(K) = dg_here%RHS_QY(K,global_here%EL_EX,dg_here%IRK)

#ifdef TRACE
            RHS_iota_IN(K) = dg_here%RHS_iota(K,global_here%EL_IN,dg_here%IRK)
            RHS_iota_EX(K) = dg_here%RHS_iota(K,global_here%EL_EX,dg_here%IRK)
#endif

#ifdef CHEM
            RHS_iota_IN(K) = dg_here%RHS_iota(K,global_here%EL_IN,dg_here%IRK)
            RHS_iota2_IN(K) = dg_here%RHS_iota2(K,global_here%EL_IN,dg_here%IRK)
            RHS_iota_EX(K) = dg_here%RHS_iota(K,global_here%EL_EX,dg_here%IRK)
            RHS_iota2_EX(K) = dg_here%RHS_iota2(K,global_here%EL_EX,dg_here%IRK)
#endif

#ifdef dynp
            RHS_dynP_IN(K) = dg_here%RHS_dynP(K,global_here%EL_IN,dg_here%IRK)
            RHS_dynP_EX(K) = dg_here%RHS_dynP(K,global_here%EL_EX,dg_here%IRK)
#endif

#ifdef SED_LAY
            do ll=1,s%layers
               dg_here%RHS_bed_IN(K,ll) = dg_here%RHS_bed(k,global_here%EL_IN,dg_here%IRK,ll)
               dg_here%RHS_bed_EX(K,ll) = dg_here%RHS_bed(k,global_here%EL_EX,dg_here%IRK,ll)
            enddo
#endif

         ENDDO

!.....Compute the sum of the lengths of three edges

         XLEN_EL_IN = dg_here%XLEN(global_here%NELED(1,global_here%EL_IN))
         XLEN_EL_IN = XLEN_EL_IN + dg_here%XLEN(global_here%NELED(2,global_here%EL_IN))
         XLEN_EL_IN = XLEN_EL_IN + dg_here%XLEN(global_here%NELED(3,global_here%EL_IN))

         XLEN_EL_EX = dg_here%XLEN(global_here%NELED(1,global_here%EL_EX))
         XLEN_EL_EX = XLEN_EL_EX + dg_here%XLEN(global_here%NELED(2,global_here%EL_EX))
         XLEN_EL_EX = XLEN_EL_EX + dg_here%XLEN(global_here%NELED(3,global_here%EL_EX))

!.....Compute the total mass in the elements

         MASS_EL_IN = (dg_here%ZE(1,global_here%EL_IN,dg_here%IRK)+dg_here%HB(1,global_here%EL_IN,1))*global_here%AREAS(global_here%EL_IN)*0.5D0
         MASS_EL_EX = (dg_here%ZE(1,global_here%EL_EX,dg_here%IRK)+dg_here%HB(1,global_here%EL_EX,1))*global_here%AREAS(global_here%EL_EX)*0.5D0

                  !Must add up the layers as they shift
#ifdef SED_LAY
         dg_here%HB(1,global_here%EL_IN,dg_here%irk) = 0.D0
         dg_here%HB(1,global_here%EL_EX,dg_here%irk) = 0.D0
         do ll=1,s%layers
            dg_here%HB(1,global_here%EL_IN,dg_here%irk) = dg_here%HB(1,global_here%EL_IN,dg_here%irk) + dg_here%bed(1,global_here%EL_IN,dg_here%irk,ll)
            dg_here%HB(1,global_here%EL_EX,dg_here%irk) = dg_here%HB(1,global_here%EL_EX,dg_here%irk) + dg_here%bed(1,global_here%EL_EX,dg_here%irk,ll)
         enddo
         MASS_EL_IN = (dg_here%ZE(1,global_here%EL_IN,dg_here%IRK)+dg_here%HB(1,global_here%EL_IN,1))*global_here%AREAS(global_here%EL_IN)*0.5D0
         MASS_EL_EX = (dg_here%ZE(1,global_here%EL_EX,dg_here%IRK)+dg_here%HB(1,global_here%EL_EX,1))*global_here%AREAS(global_here%EL_EX)*0.5D0
#endif
         
!.....Retrieve the components of the normal vector to the global_here%edge

         dg_here%NX = dg_here%COSNX(GED)
         dg_here%NY = dg_here%SINNX(GED)
         
!.....Set the components for the tangential vector to the global_here%edge

         TX = -dg_here%NY
         TY =  dg_here%NX
         
         EDFAC_IN = dg_here%XLEN(GED)/global_here%AREAS(global_here%EL_IN)
         EDFAC_EX = dg_here%XLEN(GED)/global_here%AREAS(global_here%EL_EX)

!.....Compute dg_here%ZE, dg_here%QX, dg_here%QY, and dg_here%HB at each global_here%edge Gauss quadrature point

         DO I = 1,dg_here%NEGP(dg_here%pa)

            GP_IN = I
            GP_EX = dg_here%NEGP(dg_here%pa) - I + 1

            dg_here%ZE_IN = dg_here%ZE(1,global_here%EL_IN,dg_here%IRK)
            dg_here%QX_IN = dg_here%QX(1,global_here%EL_IN,dg_here%IRK)
            dg_here%QY_IN = dg_here%QY(1,global_here%EL_IN,dg_here%IRK)
            dg_here%HB_IN = dg_here%BATHED(GP_IN,LED_IN,global_here%EL_IN,dg_here%pa)

            !When layered, these change
#ifdef SED_LAY
            do ll = 1,s%layers
               dg_here%bed_IN(ll) = dg_here%bed(1,global_here%EL_IN,dg_here%irk,ll)
            enddo
            dg_here%HB_IN = dg_here%HB(1,global_here%EL_IN,dg_here%irk)
#endif

            dg_here%SFAC_IN = dg_here%SFACED(GP_IN,LED_IN,global_here%EL_IN,dg_here%pa)
            
            dg_here%ZE_EX = dg_here%ZE(1,global_here%EL_EX,dg_here%IRK)
            dg_here%QX_EX = dg_here%QX(1,global_here%EL_EX,dg_here%IRK)
            dg_here%QY_EX = dg_here%QY(1,global_here%EL_EX,dg_here%IRK)
            dg_here%HB_EX = dg_here%BATHED(GP_EX,LED_EX,global_here%EL_EX,dg_here%pa)

#ifdef SED_LAY 
            do ll = 1,s%layers
               dg_here%bed_EX(ll) = dg_here%bed(1,global_here%EL_EX,dg_here%irk,ll)
            enddo
            dg_here%HB_EX = dg_here%HB(1,global_here%EL_EX,dg_here%irk)
#endif

            dg_here%SFAC_EX = dg_here%SFACED(GP_EX,LED_EX,global_here%EL_EX,dg_here%pa)


#ifdef TRACE
            dg_here%iota_IN = dg_here%iota(1,global_here%EL_IN,dg_here%IRK)
            dg_here%iota_EX = dg_here%iota(1,global_here%EL_EX,dg_here%IRK)

            TZ_X_IN = dg_here%TZ(1,1,1,global_here%EL_IN)
            TZ_Y_IN = dg_here%TZ(1,2,2,global_here%EL_IN)
            TZ_X_EX = dg_here%TZ(1,1,1,global_here%EL_EX)
            TZ_Y_EX = dg_here%TZ(1,2,2,global_here%EL_EX)
#endif

#ifdef CHEM
            dg_here%iota_IN = dg_here%iota(1,global_here%EL_IN,dg_here%IRK)
            dg_here%iota2_IN = dg_here%iota2(1,global_here%EL_IN,dg_here%IRK)
            dg_here%iota_EX = dg_here%iota(1,global_here%EL_EX,dg_here%IRK)
            dg_here%iota2_EX = dg_here%iota2(1,global_here%EL_EX,dg_here%IRK)
#endif

#ifdef dynp
            dynP_IN = dg_here%dynP(1,global_here%EL_IN,dg_here%IRK)
            dynP_EX = dg_here%dynP(1,global_here%EL_EX,dg_here%IRK)
#endif


#ifdef SED_LAY
            do ll = 1,s%layers
               dg_here%bed_IN(ll) = dg_here%bed(1,global_here%EL_IN,dg_here%IRK,ll)
               dg_here%bed_EX(ll) = dg_here%bed(1,global_here%EL_EX,dg_here%IRK,ll)
               
               MZ_X_IN(ll) = dg_here%MZ(1,1,ll,global_here%EL_IN)
               MZ_Y_IN(ll) = dg_here%MZ(1,2,ll,global_here%EL_IN)  

               MZ_X_EX(ll) = dg_here%MZ(1,1,ll,global_here%EL_EX)
               MZ_Y_EX(ll) = dg_here%MZ(1,2,ll,global_here%EL_EX)  
            enddo
#endif

                                !LDG terms
#ifdef WAVE_DIF
            HZ_X_IN = dg_here%HZ(1,1,1,global_here%EL_IN)
            HZ_Y_IN = dg_here%HZ(1,2,2,global_here%EL_IN)
            HZ_X_EX = dg_here%HZ(1,1,1,global_here%EL_EX)
            HZ_Y_EX = dg_here%HZ(1,2,2,global_here%EL_EX)
#endif
            LZ_XX_IN = dg_here%LZ(1,1,1,global_here%EL_IN)
            LZ_XY_IN = dg_here%LZ(1,1,2,global_here%EL_IN)
            LZ_YX_IN = dg_here%LZ(1,2,1,global_here%EL_IN)
            LZ_YY_IN = dg_here%LZ(1,2,2,global_here%EL_IN)

            LZ_XX_EX = dg_here%LZ(1,1,1,global_here%EL_EX)
            LZ_XY_EX = dg_here%LZ(1,1,2,global_here%EL_EX)
            LZ_YX_EX = dg_here%LZ(1,2,1,global_here%EL_EX)
            LZ_YY_EX = dg_here%LZ(1,2,2,global_here%EL_EX)

            DO K = 2,dg_here%DOFS(dg_here%EL)

               dg_here%ZE_IN = dg_here%ZE_IN + dg_here%ZE(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
               dg_here%QX_IN = dg_here%QX_IN + dg_here%QX(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
               dg_here%QY_IN = dg_here%QY_IN + dg_here%QY(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)

               dg_here%ZE_EX = dg_here%ZE_EX + dg_here%ZE(K,global_here%EL_EX,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
               dg_here%QX_EX = dg_here%QX_EX + dg_here%QX(K,global_here%EL_EX,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
               dg_here%QY_EX = dg_here%QY_EX + dg_here%QY(K,global_here%EL_EX,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)

#ifdef TRACE
               dg_here%iota_IN = dg_here%iota_IN + dg_here%iota(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
               dg_here%iota_EX = dg_here%iota_EX + dg_here%iota(K,global_here%EL_EX,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
#endif

#ifdef CHEM
               dg_here%iota_IN  = dg_here%iota_IN  + dg_here%iota(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
               dg_here%iota_EX  = dg_here%iota_EX  + dg_here%iota(K,global_here%EL_EX,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)

               dg_here%iota2_IN = dg_here%iota2_IN + dg_here%iota2(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
               dg_here%iota2_EX = dg_here%iota2_EX + dg_here%iota2(K,global_here%EL_EX,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
#endif

#ifdef dynp
               dynP_IN = dynP_IN + dg_here%dynP(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
               dynP_EX = dynP_EX + dg_here%dynP(K,global_here%EL_EX,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
#endif

#ifdef SED_LAY
               do ll=1,s%layers
                  dg_here%bed_IN(ll) = dg_here%bed_IN(ll) + dg_here%bed(K,global_here%EL_IN,dg_here%IRK,ll)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
                  dg_here%bed_EX(ll) = dg_here%bed_EX(ll) + dg_here%bed(K,global_here%EL_EX,dg_here%IRK,ll)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
                  dg_here%HB_IN = dg_here%HB_IN + dg_here%bed(k,global_here%EL_IN,dg_here%irk,ll)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
                  dg_here%HB_EX = dg_here%HB_EX + dg_here%bed(k,global_here%EL_EX,dg_here%irk,ll)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)

                  MZ_X_IN(ll) = MZ_X_IN(ll) + dg_here%MZ(K,1,ll,global_here%EL_IN)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
                  MZ_Y_IN(ll) = MZ_Y_IN(ll) + dg_here%MZ(K,2,ll,global_here%EL_IN)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)  
                  
                  MZ_X_EX(ll) = MZ_X_EX(ll) + dg_here%MZ(K,1,ll,global_here%EL_EX)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
                  MZ_Y_EX(ll) = MZ_Y_EX(ll) + dg_here%MZ(K,2,ll,global_here%EL_EX)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
               enddo
#endif

                                ! LDG terms
#ifdef WAVE_DIF
               HZ_X_IN = HZ_X_IN &
              + dg_here%HZ(K,1,1,global_here%EL_IN)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
               HZ_Y_IN = HZ_Y_IN &
              + dg_here%HZ(K,2,2,global_here%EL_IN)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
               HZ_X_EX = HZ_X_EX &
              + dg_here%HZ(K,1,1,global_here%EL_EX)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
               HZ_Y_EX = HZ_Y_EX &
              + dg_here%HZ(K,2,2,global_here%EL_EX)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
#endif
               LZ_XX_IN = LZ_XX_IN &
              + dg_here%LZ(K,1,1,global_here%EL_IN)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
               LZ_XY_IN = LZ_XY_IN &
              + dg_here%LZ(K,1,2,global_here%EL_IN)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
               LZ_YX_IN = LZ_YX_IN &
              + dg_here%LZ(K,2,1,global_here%EL_IN)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
               LZ_YY_IN = LZ_YY_IN &
              + dg_here%LZ(K,2,2,global_here%EL_IN)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)

               LZ_XX_EX = LZ_XX_EX &
              + dg_here%LZ(K,1,1,global_here%EL_EX)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
               LZ_XY_EX = LZ_XY_EX &
              + dg_here%LZ(K,1,2,global_here%EL_EX)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
               LZ_YX_EX = LZ_YX_EX &
              + dg_here%LZ(K,2,1,global_here%EL_EX)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
               LZ_YY_EX = LZ_YY_EX &
              + dg_here%LZ(K,2,2,global_here%EL_EX)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)

#ifdef TRACE
               TZ_X_IN = TZ_X_IN &
              + dg_here%TZ(K,1,1,global_here%EL_IN)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
               TZ_Y_IN = TZ_Y_IN &
              + dg_here%TZ(K,2,2,global_here%EL_IN)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
               TZ_X_EX = TZ_X_EX &
              + dg_here%TZ(K,1,1,global_here%EL_EX)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
               TZ_Y_EX = TZ_Y_EX &
              + dg_here%TZ(K,2,2,global_here%EL_EX)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)

#endif

            ENDDO

!.....Compute the numerical flux
            
            CALL NUMERICAL_FLUX(s,dg_here,global_here,IT,dg_here%test_el)

!......... dummy
            F_HAT_O  = global_here%F_HAT

#ifdef TRACE
            I_HAT_O  = global_here%I_HAT
#endif

#ifdef CHEM
            I_HAT_O  = global_here%I_HAT
            J_HAT_O  = global_here%J_HAT
#endif

#ifdef dynp
            K_HAT_O  = global_here%K_HAT
#endif

#ifdef SED_LAY
            bed_HAT_IN(:) = dg_here%bed_HAT(:)
            bed_HAT_EX(:) = dg_here%bed_HAT(:)
#endif

            G_HAT_IN = global_here%G_HAT
            H_HAT_IN = global_here%H_HAT
            G_HAT_EX = global_here%G_HAT
            H_HAT_EX = global_here%H_HAT

!.....Check if the flux is large enough to dry up the elements
!.....1.01D0 is a safty factor.

            IF ( (1.01*global_here%F_HAT*XLEN_EL_IN*dg_here%MAX_BOA_DT(dg_here%IRK).GE.MASS_EL_IN).OR.&
           (1.01*global_here%F_HAT*XLEN_EL_EX*dg_here%MAX_BOA_DT(dg_here%IRK)*(-1.D0).GE.&
           MASS_EL_EX)) THEN

!...........Put back the saved RHS

               DO K = 1,dg_here%DOFS(dg_here%EL)
                  dg_here%RHS_ZE(K,global_here%EL_IN,dg_here%IRK) = RHS_ZE_IN(K)
                  dg_here%RHS_QX(K,global_here%EL_IN,dg_here%IRK) = RHS_QX_IN(K)
                  dg_here%RHS_QY(K,global_here%EL_IN,dg_here%IRK) = RHS_QY_IN(K)

                  dg_here%RHS_ZE(K,global_here%EL_EX,dg_here%IRK) = RHS_ZE_EX(K)
                  dg_here%RHS_QX(K,global_here%EL_EX,dg_here%IRK) = RHS_QX_EX(K)
                  dg_here%RHS_QY(K,global_here%EL_EX,dg_here%IRK) = RHS_QY_EX(K)

#ifdef TRACE
                  dg_here%RHS_iota(K,global_here%EL_IN,dg_here%IRK) = RHS_iota_IN(K)
                  dg_here%RHS_iota(K,global_here%EL_EX,dg_here%IRK) = RHS_iota_EX(K)
#endif

#ifdef CHEM
                  dg_here%RHS_iota(K,global_here%EL_IN,dg_here%IRK) = RHS_iota_IN(K)
                  dg_here%RHS_iota2(K,global_here%EL_IN,dg_here%IRK) = RHS_iota2_IN(K)
                  dg_here%RHS_iota(K,global_here%EL_EX,dg_here%IRK) = RHS_iota_EX(K)
                  dg_here%RHS_iota2(K,global_here%EL_EX,dg_here%IRK) = RHS_iota2_EX(K)
#endif

#ifdef dynp
                  dg_here%RHS_dynP(K,global_here%EL_IN,dg_here%IRK) = RHS_dynP_IN(K)
                  dg_here%RHS_dynP(K,global_here%EL_EX,dg_here%IRK) = RHS_dynP_EX(K)
#endif

#ifdef SED_LAY
                  do ll = 1,s%layers
                     dg_here%RHS_bed(K,global_here%EL_IN,dg_here%IRK,ll) = dg_here%RHS_bed_IN(K,ll)
                     dg_here%RHS_bed(K,global_here%EL_EX,dg_here%IRK,ll) = dg_here%RHS_bed_EX(K,ll)
                  enddo
#endif

               ENDDO

               GOTO 100
            ENDIF
            
!........Check to make sure mass flux is not coming from a dry element

            IF (abs(global_here%F_HAT).gt.1.d-12) THEN
               !
               IF (dg_here%WDFLG(global_here%EL_IN).EQ.0) THEN
                  ! global_here%EL_IN is dry !
                  IF (global_here%F_HAT.GT.0) THEN
                     ! Flux going from the dry element (in)
                     ! On the wet global_here%side (ex): reflect boundary
                     global_here%Q_N_EXT = dg_here%QX_EX*dg_here%NX + dg_here%QY_EX*dg_here%NY
                     global_here%Q_T_EXT = dg_here%QX_EX*TX + dg_here%QY_EX*TY
                     global_here%Q_N_INT = -global_here%Q_N_EXT
                     global_here%Q_T_INT =  global_here%Q_T_EXT
                     DEN = 1.D0/(dg_here%NX*TY - dg_here%NY*TX)
                     dg_here%QX_IN = ( TY*global_here%Q_N_INT - dg_here%NY*global_here%Q_T_INT)*DEN
                     dg_here%QY_IN = (-TX*global_here%Q_N_INT + dg_here%NX*global_here%Q_T_INT)*DEN
                     dg_here%ZE_IN = dg_here%ZE_EX

#ifdef TRACE
                     dg_here%iota_IN = dg_here%iota_EX
#endif

#ifdef CHEM
                     dg_here%iota_IN = dg_here%iota_EX
                     dg_here%iota2_IN = dg_here%iota2_EX
#endif

#ifdef dynp
                     dynP_IN = dynP_EX
#endif

#ifdef SED_LAY
                     dg_here%bed_IN(:) = dg_here%bed_EX(:)
#endif

                     dg_here%HB_IN = dg_here%HB_EX
                     dg_here%SFAC_IN = dg_here%SFAC_EX
                     CALL NUMERICAL_FLUX(s,dg_here,global_here,IT,dg_here%test_el)

                     F_HAT_O  = global_here%F_HAT
                     G_HAT_EX = global_here%G_HAT
                     H_HAT_EX = global_here%H_HAT

                     ! on the dry global_here%side (in): do nothing
                     G_HAT_IN = global_here%G_HAT
                     H_HAT_IN = global_here%H_HAT
                     G_HAT_IN = 0.D0 
                     H_HAT_IN = 0.D0 

#ifdef TRACE
                     I_HAT_O  = global_here%I_HAT
#endif

#ifdef CHEM
                     I_HAT_O  = global_here%I_HAT
                     J_HAT_O  = global_here%J_HAT
#endif

#ifdef dynp
                     K_HAT_O  = global_here%K_HAT
#endif

#ifdef SED_LAY
                     bed_HAT_IN(:) = dg_here%bed_HAT(:)
                     bed_HAT_EX(:) = dg_here%bed_HAT(:)
#endif

                  ELSE
                     ! Flux coming from wet global_here%side (ex)
                     ! Leave fluxes on ex global_here%side intact and use zero gravtity flux for dry global_here%side      
                     NLEQG_TMP = global_here%NLEQG
                     global_here%NLEQG = 0.D0
                     G_TMP = global_here%G
                     global_here%G = 0.D0
                     CALL NUMERICAL_FLUX(s,dg_here,global_here,IT,dg_here%test_el)
                     global_here%NLEQG = NLEQG_TMP
                     global_here%G = G_TMP
                     G_HAT_IN = global_here%G_HAT
                     H_HAT_IN = global_here%H_HAT
#ifdef SED_LAY
                     bed_HAT_IN(:) = dg_here%bed_HAT(:)
#endif
                  ENDIF
                  
               ELSEIF (dg_here%WDFLG(global_here%EL_EX).EQ.0) THEN
                  
                  ! global_here%EL_EX is dry
                  IF (global_here%F_HAT.LT.0) THEN
                     ! Flux comming from dry size (ex)
                     ! On the wet global_here%side (in): reflect boundary
                     global_here%Q_N_INT = dg_here%QX_IN*dg_here%NX + dg_here%QY_IN*dg_here%NY
                     global_here%Q_T_INT = dg_here%QX_IN*TX + dg_here%QY_IN*TY
                     global_here%Q_N_EXT = -global_here%Q_N_INT
                     global_here%Q_T_EXT =  global_here%Q_T_INT
                     DEN = 1.D0/(dg_here%NX*TY - dg_here%NY*TX)
                     dg_here%QX_EX = ( TY*global_here%Q_N_EXT - dg_here%NY*global_here%Q_T_EXT)*DEN
                     dg_here%QY_EX = (-TX*global_here%Q_N_EXT + dg_here%NX*global_here%Q_T_EXT)*DEN
                     dg_here%ZE_EX = dg_here%ZE_IN
                     dg_here%HB_EX = dg_here%HB_IN
                     dg_here%SFAC_EX = dg_here%SFAC_IN

#ifdef TRACE
                     dg_here%iota_EX = dg_here%iota_IN
#endif

#ifdef CHEM
                     dg_here%iota_EX = dg_here%iota_IN
                     dg_here%iota2_EX = dg_here%iota2_IN
#endif

#ifdef dynp
                     dynP_EX = dynP_IN
#endif

#ifdef SED_LAY
                     dg_here%bed_EX(:) = dg_here%bed_IN(:)
#endif

                     CALL NUMERICAL_FLUX(s,dg_here,global_here,IT,dg_here%test_el)
                     F_HAT_O  = global_here%F_HAT

#ifdef TRACE
                     I_HAT_O  = global_here%I_HAT
#endif

#ifdef CHEM
                     I_HAT_O  = global_here%I_HAT
                     J_HAT_O  = global_here%J_HAT
#endif

#ifdef dynp
                     K_HAT_O  = global_here%K_HAT
#endif

#ifdef SED_LAY
                     bed_HAT_IN(:) = dg_here%bed_HAT(:)
                     bed_HAT_EX(:) = dg_here%bed_HAT(:)
#endif
                     
                     G_HAT_IN = global_here%G_HAT
                     H_HAT_IN = global_here%H_HAT

                     ! zero out mentum flux on the dry global_here%side (ex)
                     G_HAT_EX = global_here%G_HAT
                     H_HAT_EX = global_here%H_HAT
                     G_HAT_EX = 0.D0 
                     H_HAT_EX = 0.D0 
                  ELSE
                     ! Flux comming from the wet global_here%side (in)
                     ! Leave fluxes on in global_here%side intact and use zero gravtity for flux for dry (ex) global_here%side     
                     NLEQG_TMP = global_here%NLEQG
                     global_here%NLEQG = 0.D0
                     G_TMP = global_here%G
                     global_here%G = 0.D0
                     CALL NUMERICAL_FLUX(s,dg_here,global_here,IT,dg_here%test_el)
                     global_here%NLEQG = NLEQG_TMP
                     global_here%G = G_TMP
                     G_HAT_EX = global_here%G_HAT
                     H_HAT_EX = global_here%H_HAT
#ifdef TRACE
                     I_HAT_O  = global_here%I_HAT
#endif
#ifdef CHEM
                     I_HAT_O  = global_here%I_HAT
                     J_HAT_O  = global_here%J_HAT
#endif

#ifdef SED_LAY
                     bed_HAT_EX(:) =  dg_here%bed_HAT(:)
#endif
                  ENDIF
               ENDIF
            ENDIF

                                !Let us add some diffusive penalties
#ifdef WAVE_DIF
            F_HAT_O = F_HAT_O + 0.5D0*(HZ_X_IN*dg_here%SFAC_IN+HZ_X_EX*dg_here%SFAC_EX)*dg_here%NX + &
           0.5D0*(HZ_Y_IN+HZ_Y_EX)*dg_here%NY
            !F_HAT_EX = F_HAT_EX + 0.5D0*(HZ_X_IN+HZ_X_EX)*dg_here%SFAC_IN*dg_here%NX + 
            !&           0.5D0*(HZ_Y_EX+HZ_Y_EX)*dg_here%NY
#endif
            G_HAT_IN = G_HAT_IN + 0.5D0*(LZ_XX_IN*dg_here%SFAC_IN + &
           LZ_XX_EX*dg_here%SFAC_EX)*dg_here%NX&
           + 0.5D0*(LZ_XY_IN + LZ_XY_EX)*dg_here%NY
!nd try adding the penalty term in the LDG method
      !$           +ESLM*.01D0*(dg_here%QX_IN-dg_here%QX_EX)
            H_HAT_IN = H_HAT_IN + 0.5D0*(LZ_YX_IN*dg_here%SFAC_IN + &
           LZ_YX_EX*dg_here%SFAC_EX)*dg_here%NX&
           + 0.5D0*(LZ_YY_IN + LZ_YY_EX)*dg_here%NY
!nd try adding the penalty term in the LDG method
      !$           +ESLM*.01D0*(dg_here%QY_IN-dg_here%QY_EX)
            G_HAT_EX = G_HAT_EX + 0.5D0*(LZ_XX_IN*dg_here%SFAC_IN&
           + LZ_XX_EX*dg_here%SFAC_EX)*dg_here%NX&
           + 0.5D0*(LZ_XY_IN + LZ_XY_EX)*dg_here%NY
!nd try adding the penalty term in the LDG method
      !$           +ESLM*.01D0*(dg_here%QX_IN-dg_here%QX_EX)
            H_HAT_EX = H_HAT_EX + 0.5D0*(LZ_YX_IN*dg_here%SFAC_IN&
           + LZ_YX_EX*dg_here%SFAC_EX)*dg_here%NX&
           + 0.5D0*(LZ_YY_IN + LZ_YY_EX)*dg_here%NY
!nd try adding the penalty term in the LDG method
      !$           +ESLM*.01D0*(dg_here%QY_IN-dg_here%QY_EX) 

#ifdef TRACE

            I_HAT_O = I_HAT_O + 0.5D0*(TZ_X_IN*dg_here%SFAC_IN+TZ_X_EX*dg_here%SFAC_EX)*dg_here%NX + &
           0.5D0*(TZ_Y_IN+TZ_Y_EX)*dg_here%NY
      !&           +ESLM*.01D0*(dg_here%iota_IN-dg_here%iota_EX)

#endif
     

                                !Let us add some sediment diffusion
#ifdef SED_LAY
            do ll=1,s%layers
               bed_HAT_IN(ll) = bed_HAT_IN(ll) + 0.5D0 * ( (MZ_X_IN(ll)+&
              MZ_X_EX(ll))*dg_here%SFAC_IN*dg_here%NX+(MZ_Y_IN(ll)+MZ_Y_EX(ll))*dg_here%NY)

               bed_HAT_EX(ll) = bed_HAT_EX(ll) + 0.5D0 * ( (MZ_X_IN(ll)+&
              MZ_X_EX(ll))*dg_here%SFAC_IN*dg_here%NX+(MZ_Y_IN(ll)+MZ_Y_EX(ll))*dg_here%NY)
            enddo
#endif

!.....Compute the global_here%edge integral

            DO K = 1,dg_here%DOFS(dg_here%EL)
               
               W_IN = EDFAC_IN*dg_here%EDGEQ(K,GP_IN,LED_IN,dg_here%pa)
               W_EX = EDFAC_EX*dg_here%EDGEQ(K,GP_EX,LED_EX,dg_here%pa)

               dg_here%RHS_ZE(K,global_here%EL_IN,dg_here%IRK) = dg_here%RHS_ZE(K,global_here%EL_IN,dg_here%IRK) - W_IN*F_HAT_O
               dg_here%RHS_QX(K,global_here%EL_IN,dg_here%IRK) = dg_here%RHS_QX(K,global_here%EL_IN,dg_here%IRK) - W_IN*G_HAT_IN
               dg_here%RHS_QY(K,global_here%EL_IN,dg_here%IRK) = dg_here%RHS_QY(K,global_here%EL_IN,dg_here%IRK) - W_IN*H_HAT_IN

               dg_here%RHS_ZE(K,global_here%EL_EX,dg_here%IRK) = dg_here%RHS_ZE(K,global_here%EL_EX,dg_here%IRK) + W_EX*F_HAT_O
               dg_here%RHS_QX(K,global_here%EL_EX,dg_here%IRK) = dg_here%RHS_QX(K,global_here%EL_EX,dg_here%IRK) + W_EX*G_HAT_EX
               dg_here%RHS_QY(K,global_here%EL_EX,dg_here%IRK) = dg_here%RHS_QY(K,global_here%EL_EX,dg_here%IRK) + W_EX*H_HAT_EX

#ifdef TRACE

               dg_here%RHS_iota(K,global_here%EL_IN,dg_here%IRK) = dg_here%RHS_iota(K,global_here%EL_IN,dg_here%IRK) - W_IN*I_HAT_O 
               dg_here%RHS_iota(K,global_here%EL_EX,dg_here%IRK) = dg_here%RHS_iota(K,global_here%EL_EX,dg_here%IRK) + W_EX*I_HAT_O
#endif

#ifdef CHEM
               dg_here%RHS_iota(K,global_here%EL_IN,dg_here%IRK)  = dg_here%RHS_iota(K,global_here%EL_IN,dg_here%IRK)  - W_IN*I_HAT_O 
               dg_here%RHS_iota2(K,global_here%EL_IN,dg_here%IRK) = dg_here%RHS_iota2(K,global_here%EL_IN,dg_here%IRK) - W_IN*J_HAT_O 
               dg_here%RHS_iota(K,global_here%EL_EX,dg_here%IRK)  = dg_here%RHS_iota(K,global_here%EL_EX,dg_here%IRK)  + W_EX*I_HAT_O
               dg_here%RHS_iota2(K,global_here%EL_EX,dg_here%IRK) = dg_here%RHS_iota2(K,global_here%EL_EX,dg_here%IRK) + W_EX*J_HAT_O
#endif

#ifdef dynp
               dg_here%RHS_dynP(K,global_here%EL_IN,dg_here%IRK) = dg_here%RHS_dynP(K,global_here%EL_IN,dg_here%IRK) - W_IN*K_HAT_O 
               dg_here%RHS_dynP(K,global_here%EL_EX,dg_here%IRK) = dg_here%RHS_dynP(K,global_here%EL_EX,dg_here%IRK) + W_EX*K_HAT_O
#endif


#ifdef SED_LAY
               do ll = 1,s%layers
                  dg_here%RHS_bed(K,global_here%EL_IN,dg_here%IRK,ll) = dg_here%RHS_bed(K,global_here%EL_IN,dg_here%IRK,ll) - W_IN*bed_HAT_IN(ll) 
                  dg_here%RHS_bed(K,global_here%EL_EX,dg_here%IRK,ll) = dg_here%RHS_bed(K,global_here%EL_EX,dg_here%IRK,ll) + W_EX*bed_HAT_EX(ll)
               enddo
#endif

            ENDDO

         ENDDO

         GOTO 1000

 100     CONTINUE

!--------------------------------------------------------------------
!     Compute reflection flux
!--------------------------------------------------------------------

!.....Set the components for the tangential vector to the global_here%edge

         TX = -dg_here%NY
         TY =  dg_here%NX

!.....Compute dg_here%ZE, dg_here%QX, dg_here%QY, dg_here%bed, and dg_here%HB at each global_here%edge Gauss quadrature point

         DO I = 1,dg_here%NEGP(dg_here%pa)

            GP_IN = I
            GP_EX = dg_here%NEGP(dg_here%pa) - I + 1

!     ----------------------------- global_here%EL_IN ------------------------------

            dg_here%ZE_IN = dg_here%ZE(1,global_here%EL_IN,dg_here%IRK)
            dg_here%QX_IN = dg_here%QX(1,global_here%EL_IN,dg_here%IRK)
            dg_here%QY_IN = dg_here%QY(1,global_here%EL_IN,dg_here%IRK)
            dg_here%HB_IN = dg_here%BATHED(GP_IN,LED_IN,global_here%EL_IN,dg_here%pa)

            dg_here%SFAC_IN = dg_here%SFACED(GP_IN,LED_IN,global_here%EL_IN,dg_here%pa)

#ifdef TRACE
            dg_here%iota_IN = dg_here%iota(1,global_here%EL_IN,dg_here%IRK)
#endif

#ifdef CHEM
            dg_here%iota_IN = dg_here%iota(1,global_here%EL_IN,dg_here%IRK)
            dg_here%iota2_IN = dg_here%iota2(1,global_here%EL_IN,dg_here%IRK)
#endif

#ifdef dynp
            dynP_IN = dg_here%dynP(1,global_here%EL_IN,dg_here%IRK)
#endif
            
            !Must add up the layers as they shift
#ifdef SED_LAY
            dg_here%HB(1,global_here%EL_IN,dg_here%irk) = 0.D0
            do ll=1,s%layers
               dg_here%bed_IN(ll) = dg_here%bed(1,global_here%EL_IN,dg_here%irk,ll)
               dg_here%HB(1,global_here%EL_IN,dg_here%irk) = dg_here%HB(1,global_here%EL_IN,dg_here%irk) + dg_here%bed(1,global_here%EL_IN,dg_here%irk,ll)

               !sediment diffusion
               MZ_X_IN(ll) = dg_here%MZ(1,1,ll,global_here%EL_IN)
               MZ_Y_IN(ll) = dg_here%MZ(1,2,ll,global_here%EL_IN)  
            enddo
            dg_here%HB_IN = dg_here%HB(1,global_here%EL_IN,dg_here%irk)
#endif

            !diffusion
#ifdef WAVE_DIF
            HZ_X_IN = dg_here%HZ(1,1,1,global_here%EL_IN)
            HZ_Y_IN = dg_here%HZ(1,2,2,global_here%EL_IN)
#endif
            LZ_XX_IN = dg_here%LZ(1,1,1,global_here%EL_IN)
            LZ_XY_IN = dg_here%LZ(1,1,2,global_here%EL_IN)
            LZ_YX_IN = dg_here%LZ(1,2,1,global_here%EL_IN)
            LZ_YY_IN = dg_here%LZ(1,2,2,global_here%EL_IN)

#ifdef TRACE

            TZ_X_IN = dg_here%TZ(1,1,1,global_here%EL_IN)
            TZ_Y_IN = dg_here%TZ(1,2,2,global_here%EL_IN)

#endif

            DO K = 2,dg_here%DOFS(dg_here%EL)

               dg_here%ZE_IN = dg_here%ZE_IN + dg_here%ZE(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
               dg_here%QX_IN = dg_here%QX_IN + dg_here%QX(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
               dg_here%QY_IN = dg_here%QY_IN + dg_here%QY(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)

#ifdef TRACE
               dg_here%iota_IN = dg_here%iota_IN + dg_here%iota(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%PA)
#endif

#ifdef CHEM
               dg_here%iota_IN  = dg_here%iota_IN  + dg_here%iota(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%PA)
               dg_here%iota2_IN = dg_here%iota2_IN + dg_here%iota2(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%PA)
#endif

#ifdef dynp
               dynP_IN = dynP_IN + dg_here%dynP(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%PA)
#endif

#ifdef SED_LAY
               do ll=1,s%layers
                  dg_here%bed_IN(ll) = dg_here%bed_IN(ll) + dg_here%bed(k,global_here%EL_IN,dg_here%irk,ll)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
                  dg_here%HB_IN = dg_here%HB_IN + dg_here%bed(k,global_here%EL_IN,dg_here%irk,ll)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
                  
                  MZ_X_IN(ll) = MZ_X_IN(ll) + dg_here%MZ(K,1,ll,global_here%EL_IN)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
                  MZ_Y_IN(ll) = MZ_Y_IN(ll) + dg_here%MZ(K,2,ll,global_here%EL_IN)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)  
               enddo
#endif
                                ! LDG terms
#ifdef WAVE_DIF
               HZ_X_IN = HZ_X_IN&
              + dg_here%HZ(K,1,1,global_here%EL_IN)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
               HZ_Y_IN = HZ_Y_IN&
              + dg_here%HZ(K,2,2,global_here%EL_IN)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
#endif
               LZ_XX_IN = LZ_XX_IN&
              + dg_here%LZ(K,1,1,global_here%EL_IN)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
               LZ_XY_IN = LZ_XY_IN&
              + dg_here%LZ(K,1,2,global_here%EL_IN)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
               LZ_YX_IN = LZ_YX_IN&
              + dg_here%LZ(K,2,1,global_here%EL_IN)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
               LZ_YY_IN = LZ_YY_IN&
              + dg_here%LZ(K,2,2,global_here%EL_IN)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)

#ifdef TRACE
               TZ_X_IN = TZ_X_IN&
              + dg_here%TZ(K,1,1,global_here%EL_IN)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
               TZ_Y_IN = TZ_Y_IN&
              + dg_here%TZ(K,2,2,global_here%EL_IN)*dg_here%PHI_EDGE(K,GP_IN,LED_IN,dg_here%pa)
#endif

            ENDDO

!.....Compute the velocity in the normal and tangental direction

            global_here%Q_N_INT = dg_here%QX_IN*dg_here%NX + dg_here%QY_IN*dg_here%NY
            global_here%Q_T_INT = dg_here%QX_IN*TX + dg_here%QY_IN*TY

!.....Reflect the velocity in the normal direction

            global_here%Q_N_EXT = -global_here%Q_N_INT
            global_here%Q_T_EXT =  global_here%Q_T_INT

!.....Compute the global_here%x and global_here%y components of the external state flow

            DEN = 1.D0/(dg_here%NX*TY - dg_here%NY*TX)
            dg_here%QX_EX = ( TY*global_here%Q_N_EXT - dg_here%NY*global_here%Q_T_EXT)*DEN
            dg_here%QY_EX = (-TX*global_here%Q_N_EXT + dg_here%NX*global_here%Q_T_EXT)*DEN

            dg_here%ZE_EX = dg_here%ZE_IN
            dg_here%HB_EX = dg_here%HB_IN
            dg_here%SFAC_EX = dg_here%SFAC_IN

#ifdef TRACE
            dg_here%iota_EX = dg_here%iota_IN
#endif

#ifdef CHEM
            dg_here%iota_EX = dg_here%iota_IN
            dg_here%iota2_EX = dg_here%iota2_IN
#endif

#ifdef dynp
            dynP_EX = dynP_IN
#endif

#ifdef SED_LAY
            dg_here%bed_EX(:) = dg_here%bed_IN(:)
#endif

!.....Compute the numerical flux

!     CALL NUMERICAL_FLUX(IT)
            IF(dg_here%WDFLG(global_here%EL_IN).EQ.0) THEN
               NLEQG_TMP = global_here%NLEQG
               global_here%NLEQG = 0.D0
               G_TMP = global_here%G
               global_here%G = 0.D0
               CALL NUMERICAL_FLUX(s,dg_here,global_here,IT,dg_here%test_el)
               global_here%NLEQG = NLEQG_TMP
               global_here%G = G_TMP
            ELSE
               CALL NUMERICAL_FLUX(s,dg_here,global_here,IT,dg_here%test_el)
            ENDIF
            global_here%F_HAT = 0.D0        ! Ensure zero mass flux

!.....Add LDG terms for eddy viscosity

#ifdef WAVE_DIF
            global_here%F_HAT = global_here%F_HAT + HZ_X_IN*dg_here%NX*dg_here%SFAC_IN&
           + HZ_Y_IN*dg_here%NY
#endif
            global_here%G_HAT = global_here%G_HAT + LZ_XX_IN*dg_here%NX*dg_here%SFAC_IN&
           + LZ_XY_IN*dg_here%NY
            global_here%H_HAT = global_here%H_HAT + LZ_YX_IN*dg_here%NX*dg_here%SFAC_IN&
           + LZ_YY_IN*dg_here%NY

#ifdef TRACE

            global_here%I_HAT = global_here%I_HAT + TZ_X_IN*dg_here%NX*dg_here%SFAC_IN&
           + TZ_Y_IN*dg_here%NY

#endif

!.....Add LDG terms for sediment diffusion

#ifdef SED_LAY
            do ll=1,s%layers
               dg_here%bed_HAT(ll) = dg_here%bed_HAT(ll) + MZ_X_IN(ll)*dg_here%SFAC_IN*dg_here%NX &
              + MZ_Y_IN(ll)*dg_here%NY
            enddo
#endif

!.....Compute the global_here%edge integral
            DO K = 1,dg_here%DOFS(dg_here%EL)

               W_IN = EDFAC_IN*dg_here%EDGEQ(K,GP_IN,LED_IN,dg_here%pa)

               dg_here%RHS_ZE(K,global_here%EL_IN,dg_here%IRK) = dg_here%RHS_ZE(K,global_here%EL_IN,dg_here%IRK) - W_IN*global_here%F_HAT
               dg_here%RHS_QX(K,global_here%EL_IN,dg_here%IRK) = dg_here%RHS_QX(K,global_here%EL_IN,dg_here%IRK) - W_IN*global_here%G_HAT
               dg_here%RHS_QY(K,global_here%EL_IN,dg_here%IRK) = dg_here%RHS_QY(K,global_here%EL_IN,dg_here%IRK) - W_IN*global_here%H_HAT

#ifdef TRACE
               dg_here%RHS_iota(K,global_here%EL_IN,dg_here%IRK) = dg_here%RHS_iota(K,global_here%EL_IN,dg_here%IRK) - W_IN*global_here%I_HAT
#endif

#ifdef CHEM
               dg_here%RHS_iota(K,global_here%EL_IN,dg_here%IRK)  = dg_here%RHS_iota(K,global_here%EL_IN,dg_here%IRK)  - W_IN*global_here%I_HAT
               dg_here%RHS_iota2(K,global_here%EL_IN,dg_here%IRK) = dg_here%RHS_iota2(K,global_here%EL_IN,dg_here%IRK) - W_IN*global_here%J_HAT
#endif

#ifdef dynp
               dg_here%RHS_dynP(K,global_here%EL_IN,dg_here%IRK) = dg_here%RHS_dynP(K,global_here%EL_IN,dg_here%IRK) - W_IN*global_here%K_HAT
#endif

#ifdef SED_LAY
               do ll = 1,s%layers
                  dg_here%RHS_bed(K,global_here%EL_IN,dg_here%IRK,ll) = dg_here%RHS_bed(K,global_here%EL_IN,dg_here%IRK,ll) - W_IN*dg_here%bed_HAT(ll) 
               enddo
#endif

            ENDDO

!     ----------------------------- global_here%EL_EX ------------------------------

            dg_here%ZE_EX = dg_here%ZE(1,global_here%EL_EX,dg_here%IRK)
            dg_here%QX_EX = dg_here%QX(1,global_here%EL_EX,dg_here%IRK)
            dg_here%QY_EX = dg_here%QY(1,global_here%EL_EX,dg_here%IRK)
            dg_here%HB_EX = dg_here%BATHED(GP_EX,LED_EX,global_here%EL_EX,dg_here%pa)

            dg_here%SFAC_EX = dg_here%SFACED(GP_EX,LED_EX,global_here%EL_EX,dg_here%pa)

#ifdef TRACE
            dg_here%iota_EX = dg_here%iota(1,global_here%EL_EX,dg_here%IRK)
#endif

#ifdef CHEM
            dg_here%iota_EX = dg_here%iota(1,global_here%EL_EX,dg_here%IRK)
            dg_here%iota2_EX = dg_here%iota2(1,global_here%EL_EX,dg_here%IRK)
#endif

#ifdef dynp
            dynP_EX = dg_here%dynP(1,global_here%EL_EX,dg_here%IRK)
#endif

            !Must add up the layers as they shift
#ifdef SED_LAY
            dg_here%HB(1,global_here%EL_EX,dg_here%irk) = 0.D0
            do ll=1,s%layers
               dg_here%bed_EX(ll) = dg_here%bed(1,global_here%EL_EX,dg_here%irk,ll)
               dg_here%HB(1,global_here%EL_EX,dg_here%irk) = dg_here%HB(1,global_here%EL_EX,dg_here%irk) + dg_here%bed(1,global_here%EL_EX,dg_here%irk,ll)
               !sediment diffusion
               MZ_X_EX(ll) = dg_here%MZ(1,1,ll,global_here%EL_EX)
               MZ_Y_EX(ll) = dg_here%MZ(1,2,ll,global_here%EL_EX)  
            enddo
            dg_here%HB_EX = dg_here%HB(1,global_here%EL_EX,dg_here%irk)
#endif

            !LDG
#ifdef WAVE_DIF
            HZ_X_EX = dg_here%HZ(1,1,1,global_here%EL_EX)
            HZ_Y_EX = dg_here%HZ(1,2,2,global_here%EL_EX)
#endif
            LZ_XX_EX = dg_here%LZ(1,1,1,global_here%EL_EX)
            LZ_XY_EX = dg_here%LZ(1,1,2,global_here%EL_EX)
            LZ_YX_EX = dg_here%LZ(1,2,1,global_here%EL_EX)
            LZ_YY_EX = dg_here%LZ(1,2,2,global_here%EL_EX)

#ifdef TRACE

            TZ_X_EX = dg_here%TZ(1,1,1,global_here%EL_EX)
            TZ_Y_EX = dg_here%TZ(1,2,2,global_here%EL_EX)

#endif

            DO K = 2,dg_here%DOFS(dg_here%EL)

               dg_here%ZE_EX = dg_here%ZE_EX + dg_here%ZE(K,global_here%EL_EX,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
               dg_here%QX_EX = dg_here%QX_EX + dg_here%QX(K,global_here%EL_EX,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
               dg_here%QY_EX = dg_here%QY_EX + dg_here%QY(K,global_here%EL_EX,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)

#ifdef TRACE
               dg_here%iota_EX = dg_here%iota_EX + dg_here%iota(K,global_here%EL_EX,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
#endif

#ifdef CHEM
               dg_here%iota_EX  = dg_here%iota_EX  + dg_here%iota(K,global_here%EL_EX,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
               dg_here%iota2_EX = dg_here%iota2_EX + dg_here%iota2(K,global_here%EL_EX,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
#endif

#ifdef dynp
               dynP_EX = dynP_EX + dg_here%dynP(K,global_here%EL_EX,dg_here%IRK)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
#endif

#ifdef SED_LAY
               do ll=1,s%layers
                  dg_here%bed_EX(ll) = dg_here%bed_EX(ll) + dg_here%bed(k,global_here%EL_EX,dg_here%irk,ll)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
                  dg_here%HB_EX = dg_here%HB_EX + dg_here%bed(k,global_here%EL_EX,dg_here%irk,ll)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
                  
                  MZ_X_EX(ll) = MZ_X_EX(ll) + dg_here%MZ(K,1,ll,global_here%EL_EX)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
                  MZ_Y_EX(ll) = MZ_Y_EX(ll) + dg_here%MZ(K,2,ll,global_here%EL_EX)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa) 
               enddo
#endif

                                ! LDG terms
#ifdef WAVE_DIF
               HZ_X_EX = HZ_X_EX&
              + dg_here%HZ(K,1,1,global_here%EL_EX)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
               HZ_Y_EX = HZ_Y_EX&
              + dg_here%HZ(K,2,2,global_here%EL_EX)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
#endif
               LZ_XX_EX = LZ_XX_EX&
              + dg_here%LZ(K,1,1,global_here%EL_EX)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
               LZ_XY_EX = LZ_XY_EX&
              + dg_here%LZ(K,1,2,global_here%EL_EX)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
               LZ_YX_EX = LZ_YX_EX&
              + dg_here%LZ(K,2,1,global_here%EL_EX)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
               LZ_YY_EX = LZ_YY_EX&
              + dg_here%LZ(K,2,2,global_here%EL_EX)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)

#ifdef TRACE

               TZ_X_EX = TZ_X_EX&
              + dg_here%TZ(K,1,1,global_here%EL_EX)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)
               TZ_Y_EX = TZ_Y_EX&
              + dg_here%TZ(K,2,2,global_here%EL_EX)*dg_here%PHI_EDGE(K,GP_EX,LED_EX,dg_here%pa)

#endif

                                ! LDG terms for sediment diffusion               

            ENDDO

!.....Compute the velocity in the normal and tangental direction

            global_here%Q_N_EXT = dg_here%QX_EX*dg_here%NX + dg_here%QY_EX*dg_here%NY
            global_here%Q_T_EXT = dg_here%QX_EX*TX + dg_here%QY_EX*TY

!.....Reflect the velocity in the normal direction

            global_here%Q_N_INT = -global_here%Q_N_EXT
            global_here%Q_T_INT =  global_here%Q_T_EXT

!.....Compute the global_here%x and global_here%y components of the external state flow

            DEN = 1.D0/(dg_here%NX*TY - dg_here%NY*TX)
            dg_here%QX_IN = ( TY*global_here%Q_N_INT - dg_here%NY*global_here%Q_T_INT)*DEN
            dg_here%QY_IN = (-TX*global_here%Q_N_INT + dg_here%NX*global_here%Q_T_INT)*DEN

            dg_here%ZE_IN = dg_here%ZE_EX
            dg_here%HB_IN = dg_here%HB_EX

            dg_here%SFAC_IN = dg_here%SFAC_EX

#ifdef TRACE
            dg_here%iota_IN = dg_here%iota_EX
#endif

#ifdef CHEM
            dg_here%iota_IN = dg_here%iota_EX
            dg_here%iota2_IN = dg_here%iota2_EX
#endif

#ifdef dynp
            dynP_IN = dynP_EX
#endif

#ifdef SED_LAY
            dg_here%bed_IN(:) = dg_here%bed_EX(:)
#endif

!.....Compute the numerical flux

!     CALL NUMERICAL_FLUX(IT,L)
            IF(dg_here%WDFLG(global_here%EL_EX).EQ.0) THEN
               NLEQG_TMP = global_here%NLEQG
               global_here%NLEQG = 0.D0
               G_TMP = global_here%G
               global_here%G = 0.D0
               CALL NUMERICAL_FLUX(s,dg_here,global_here,IT,dg_here%test_el)
               global_here%NLEQG = NLEQG_TMP
               global_here%G = G_TMP
            ELSE
               CALL NUMERICAL_FLUX(s,dg_here,global_here,IT,dg_here%test_el)
            ENDIF
            global_here%F_HAT = 0.D0        ! Ensure zero mass flux

!.....Add LDG terms for eddy viscosity
#ifdef WAVE_DIF
            global_here%F_HAT = global_here%F_HAT + HZ_X_EX*dg_here%NX*dg_here%SFAC_EX&
           + HZ_Y_EX*dg_here%NY
#endif
            global_here%G_HAT = global_here%G_HAT + LZ_XX_EX*dg_here%NX*dg_here%SFAC_EX&
           + LZ_XY_EX*dg_here%NY
            global_here%H_HAT = global_here%H_HAT + LZ_YX_EX*dg_here%NX*dg_here%SFAC_EX&
           + LZ_YY_EX*dg_here%NY

#ifdef TRACE

            global_here%I_HAT = global_here%I_HAT + TZ_X_EX*dg_here%NX*dg_here%SFAC_EX&
           + TZ_Y_EX*dg_here%NY

#endif

!.....Add LDG terms for sediment diffusion

#ifdef SED_LAY
            do ll=1,s%layers
               dg_here%bed_HAT(ll) = dg_here%bed_HAT(ll) + MZ_X_EX(ll)*dg_here%SFAC_EX*dg_here%NX &
              + MZ_Y_EX(ll)*dg_here%NY
            enddo
#endif

!.....Compute the global_here%edge integral
            DO K = 1,dg_here%DOFS(dg_here%EL)

               W_EX = EDFAC_EX*dg_here%EDGEQ(K,GP_EX,LED_EX,dg_here%pa)

               dg_here%RHS_ZE(K,global_here%EL_EX,dg_here%IRK) = dg_here%RHS_ZE(K,global_here%EL_EX,dg_here%IRK) + W_EX*global_here%F_HAT
               dg_here%RHS_QX(K,global_here%EL_EX,dg_here%IRK) = dg_here%RHS_QX(K,global_here%EL_EX,dg_here%IRK) + W_EX*global_here%G_HAT
               dg_here%RHS_QY(K,global_here%EL_EX,dg_here%IRK) = dg_here%RHS_QY(K,global_here%EL_EX,dg_here%IRK) + W_EX*global_here%H_HAT

#ifdef TRACE
               dg_here%RHS_iota(K,global_here%EL_EX,dg_here%IRK) = dg_here%RHS_iota(K,global_here%EL_EX,dg_here%IRK) + W_EX*global_here%I_HAT
#endif

#ifdef CHEM
               dg_here%RHS_iota(K,global_here%EL_EX,dg_here%IRK)  = dg_here%RHS_iota(K,global_here%EL_EX,dg_here%IRK)  + W_EX*global_here%I_HAT
               dg_here%RHS_iota2(K,global_here%EL_EX,dg_here%IRK) = dg_here%RHS_iota2(K,global_here%EL_EX,dg_here%IRK) + W_EX*global_here%J_HAT
#endif

#ifdef dynp
               dg_here%RHS_dynP(K,global_here%EL_EX,dg_here%IRK) = dg_here%RHS_dynP(K,global_here%EL_EX,dg_here%IRK) + W_EX*global_here%K_HAT
#endif

#ifdef SED_LAY
               do ll = 1,s%layers
                  dg_here%RHS_bed(K,global_here%EL_EX,dg_here%IRK,ll) = dg_here%RHS_bed(K,global_here%EL_EX,dg_here%IRK,ll) + W_EX*dg_here%bed_HAT(ll) 
               enddo
#endif

            ENDDO

         ENDDO
         
 1000 CONTINUE

      RETURN
      END SUBROUTINE
