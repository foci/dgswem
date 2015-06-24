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

      SUBROUTINE INTERNAL_EDGE_HYDRO(s,IT)

!.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE NodalAttributes, ONLY :  ESLM
      use sizes

      IMPLICIT NONE

      type (sizes_type) :: s

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
      REAL(SZ), SAVE, ALLOCATABLE :: &
     RHS_ZE_IN(:), RHS_QX_IN(:), RHS_QY_IN(:),&
     RHS_ZE_EX(:), RHS_QX_EX(:), RHS_QY_EX(:)
#ifdef TRACE     
      REAL(SZ), SAVE, ALLOCATABLE :: &
     ,RHS_iota_IN(:), RHS_iota_EX(:)
#endif
#ifdef CHEM
      REAL(SZ), SAVE, ALLOCATABLE :: &
     ,RHS_iota_IN(:), RHS_iota_EX(:)&
     ,RHS_iota2_IN(:), RHS_iota2_EX(:)
#endif
#ifdef dynp     
      REAL(SZ), SAVE, ALLOCATABLE :: &
     ,RHS_dynP_IN(:), RHS_dynP_EX(:)
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

      IF(.NOT.ALLOCATED(RHS_ZE_IN)) THEN
         ALLOCATE ( RHS_ZE_IN(dg%DOFH),RHS_QX_IN(dg%DOFH),RHS_QY_IN(dg%DOFH) )
         ALLOCATE ( RHS_ZE_EX(dg%DOFH),RHS_QX_EX(dg%DOFH),RHS_QY_EX(dg%DOFH) )
#ifdef TRACE
         Allocate ( rhs_iota_in(dg%dofh),rhs_iota_ex(dg%dofh) )
#endif
#ifdef CHEM
         Allocate ( rhs_iota_in(dg%dofh),rhs_iota2_in(dg%dofh) )
         Allocate ( rhs_iota_ex(dg%dofh),rhs_iota2_ex(dg%dofh) )
#endif
#ifdef dynp
         Allocate ( rhs_dynP_IN(dg%dofh),rhs_dynP_EX(dg%dofh) )
#endif
      ENDIF

      dg%test_el = 0
      DO 1000 L = 1,dg%NIEDS

!.......Retrieve the global and local edge number

         GED = dg%NIEDN(L)
         LED_IN = dg%NEDSD(1,GED)
         LED_EX = dg%NEDSD(2,GED)

!.......Retrieve the elements which share the edge

         EL_IN = dg%NEDEL(1,GED)
         EL_EX = dg%NEDEL(2,GED)

         dg%EL = EL_EX

         if (dg%DOFS(EL_EX).LT.dg%DOFS(EL_IN)) then
            dg%EL = EL_IN
         endif
         
         dg%pa = PDG_EL(dg%EL)

#ifdef P0
         if (dg%pa.eq.0) then
            dg%pa = 1
         endif
#endif

!.......If both elements on either side of edge are dry then skip

         IF((dg%WDFLG(EL_IN).EQ.0).AND.(dg%WDFLG(EL_EX).EQ.0)) GOTO 1000

         !dg%test_el = dg%test_el+1

!.....Save the current RHS

         DO K = 1,dg%DOFS(dg%EL)

            RHS_ZE_IN(K) = dg%RHS_ZE(K,EL_IN,dg%IRK)
            RHS_QX_IN(K) = dg%RHS_QX(K,EL_IN,dg%IRK)
            RHS_QY_IN(K) = dg%RHS_QY(K,EL_IN,dg%IRK)

            RHS_ZE_EX(K) = dg%RHS_ZE(K,EL_EX,dg%IRK)
            RHS_QX_EX(K) = dg%RHS_QX(K,EL_EX,dg%IRK)
            RHS_QY_EX(K) = dg%RHS_QY(K,EL_EX,dg%IRK)

#ifdef TRACE
            RHS_iota_IN(K) = dg%RHS_iota(K,EL_IN,dg%IRK)
            RHS_iota_EX(K) = dg%RHS_iota(K,EL_EX,dg%IRK)
#endif

#ifdef CHEM
            RHS_iota_IN(K) = dg%RHS_iota(K,EL_IN,dg%IRK)
            RHS_iota2_IN(K) = dg%RHS_iota2(K,EL_IN,dg%IRK)
            RHS_iota_EX(K) = dg%RHS_iota(K,EL_EX,dg%IRK)
            RHS_iota2_EX(K) = dg%RHS_iota2(K,EL_EX,dg%IRK)
#endif

#ifdef dynp
            RHS_dynP_IN(K) = dg%RHS_dynP(K,EL_IN,dg%IRK)
            RHS_dynP_EX(K) = dg%RHS_dynP(K,EL_EX,dg%IRK)
#endif

#ifdef SED_LAY
            do ll=1,s%layers
               dg%RHS_bed_IN(K,ll) = dg%RHS_bed(k,EL_IN,dg%IRK,ll)
               dg%RHS_bed_EX(K,ll) = dg%RHS_bed(k,EL_EX,dg%IRK,ll)
            enddo
#endif

         ENDDO

!.....Compute the sum of the lengths of three edges

         XLEN_EL_IN = dg%XLEN(NELED(1,EL_IN))
         XLEN_EL_IN = XLEN_EL_IN + dg%XLEN(NELED(2,EL_IN))
         XLEN_EL_IN = XLEN_EL_IN + dg%XLEN(NELED(3,EL_IN))

         XLEN_EL_EX = dg%XLEN(NELED(1,EL_EX))
         XLEN_EL_EX = XLEN_EL_EX + dg%XLEN(NELED(2,EL_EX))
         XLEN_EL_EX = XLEN_EL_EX + dg%XLEN(NELED(3,EL_EX))

!.....Compute the total mass in the elements

         MASS_EL_IN = (dg%ZE(1,EL_IN,dg%IRK)+dg%HB(1,EL_IN,1))*AREAS(EL_IN)*0.5D0
         MASS_EL_EX = (dg%ZE(1,EL_EX,dg%IRK)+dg%HB(1,EL_EX,1))*AREAS(EL_EX)*0.5D0

                  !Must add up the layers as they shift
#ifdef SED_LAY
         dg%HB(1,EL_IN,dg%irk) = 0.D0
         dg%HB(1,EL_EX,dg%irk) = 0.D0
         do ll=1,s%layers
            dg%HB(1,EL_IN,dg%irk) = dg%HB(1,EL_IN,dg%irk) + dg%bed(1,EL_IN,dg%irk,ll)
            dg%HB(1,EL_EX,dg%irk) = dg%HB(1,EL_EX,dg%irk) + dg%bed(1,EL_EX,dg%irk,ll)
         enddo
         MASS_EL_IN = (dg%ZE(1,EL_IN,dg%IRK)+dg%HB(1,EL_IN,1))*AREAS(EL_IN)*0.5D0
         MASS_EL_EX = (dg%ZE(1,EL_EX,dg%IRK)+dg%HB(1,EL_EX,1))*AREAS(EL_EX)*0.5D0
#endif
         
!.....Retrieve the components of the normal vector to the edge

         dg%NX = dg%COSNX(GED)
         dg%NY = dg%SINNX(GED)
         
!.....Set the components for the tangential vector to the edge

         TX = -dg%NY
         TY =  dg%NX
         
         EDFAC_IN = dg%XLEN(GED)/AREAS(EL_IN)
         EDFAC_EX = dg%XLEN(GED)/AREAS(EL_EX)

!.....Compute dg%ZE, dg%QX, dg%QY, and dg%HB at each edge Gauss quadrature point

         DO I = 1,dg%NEGP(dg%pa)

            GP_IN = I
            GP_EX = dg%NEGP(dg%pa) - I + 1

            dg%ZE_IN = dg%ZE(1,EL_IN,dg%IRK)
            dg%QX_IN = dg%QX(1,EL_IN,dg%IRK)
            dg%QY_IN = dg%QY(1,EL_IN,dg%IRK)
            dg%HB_IN = dg%BATHED(GP_IN,LED_IN,EL_IN,dg%pa)

            !When layered, these change
#ifdef SED_LAY
            do ll = 1,s%layers
               dg%bed_IN(ll) = dg%bed(1,EL_IN,dg%irk,ll)
            enddo
            dg%HB_IN = dg%HB(1,EL_IN,dg%irk)
#endif

            dg%SFAC_IN = dg%SFACED(GP_IN,LED_IN,EL_IN,dg%pa)
            
            dg%ZE_EX = dg%ZE(1,EL_EX,dg%IRK)
            dg%QX_EX = dg%QX(1,EL_EX,dg%IRK)
            dg%QY_EX = dg%QY(1,EL_EX,dg%IRK)
            dg%HB_EX = dg%BATHED(GP_EX,LED_EX,EL_EX,dg%pa)

#ifdef SED_LAY 
            do ll = 1,s%layers
               dg%bed_EX(ll) = dg%bed(1,EL_EX,dg%irk,ll)
            enddo
            dg%HB_EX = dg%HB(1,EL_EX,dg%irk)
#endif

            dg%SFAC_EX = dg%SFACED(GP_EX,LED_EX,EL_EX,dg%pa)


#ifdef TRACE
            dg%iota_IN = dg%iota(1,EL_IN,dg%IRK)
            dg%iota_EX = dg%iota(1,EL_EX,dg%IRK)

            TZ_X_IN = dg%TZ(1,1,1,EL_IN)
            TZ_Y_IN = dg%TZ(1,2,2,EL_IN)
            TZ_X_EX = dg%TZ(1,1,1,EL_EX)
            TZ_Y_EX = dg%TZ(1,2,2,EL_EX)
#endif

#ifdef CHEM
            dg%iota_IN = dg%iota(1,EL_IN,dg%IRK)
            dg%iota2_IN = dg%iota2(1,EL_IN,dg%IRK)
            dg%iota_EX = dg%iota(1,EL_EX,dg%IRK)
            dg%iota2_EX = dg%iota2(1,EL_EX,dg%IRK)
#endif

#ifdef dynp
            dynP_IN = dg%dynP(1,EL_IN,dg%IRK)
            dynP_EX = dg%dynP(1,EL_EX,dg%IRK)
#endif


#ifdef SED_LAY
            do ll = 1,s%layers
               dg%bed_IN(ll) = dg%bed(1,EL_IN,dg%IRK,ll)
               dg%bed_EX(ll) = dg%bed(1,EL_EX,dg%IRK,ll)
               
               MZ_X_IN(ll) = dg%MZ(1,1,ll,EL_IN)
               MZ_Y_IN(ll) = dg%MZ(1,2,ll,EL_IN)  

               MZ_X_EX(ll) = dg%MZ(1,1,ll,EL_EX)
               MZ_Y_EX(ll) = dg%MZ(1,2,ll,EL_EX)  
            enddo
#endif

                                !LDG terms
#ifdef WAVE_DIF
            HZ_X_IN = dg%HZ(1,1,1,EL_IN)
            HZ_Y_IN = dg%HZ(1,2,2,EL_IN)
            HZ_X_EX = dg%HZ(1,1,1,EL_EX)
            HZ_Y_EX = dg%HZ(1,2,2,EL_EX)
#endif
            LZ_XX_IN = dg%LZ(1,1,1,EL_IN)
            LZ_XY_IN = dg%LZ(1,1,2,EL_IN)
            LZ_YX_IN = dg%LZ(1,2,1,EL_IN)
            LZ_YY_IN = dg%LZ(1,2,2,EL_IN)

            LZ_XX_EX = dg%LZ(1,1,1,EL_EX)
            LZ_XY_EX = dg%LZ(1,1,2,EL_EX)
            LZ_YX_EX = dg%LZ(1,2,1,EL_EX)
            LZ_YY_EX = dg%LZ(1,2,2,EL_EX)

            DO K = 2,dg%DOFS(dg%EL)

               dg%ZE_IN = dg%ZE_IN + dg%ZE(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)
               dg%QX_IN = dg%QX_IN + dg%QX(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)
               dg%QY_IN = dg%QY_IN + dg%QY(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)

               dg%ZE_EX = dg%ZE_EX + dg%ZE(K,EL_EX,dg%IRK)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)
               dg%QX_EX = dg%QX_EX + dg%QX(K,EL_EX,dg%IRK)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)
               dg%QY_EX = dg%QY_EX + dg%QY(K,EL_EX,dg%IRK)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)

#ifdef TRACE
               dg%iota_IN = dg%iota_IN + dg%iota(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)
               dg%iota_EX = dg%iota_EX + dg%iota(K,EL_EX,dg%IRK)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)
#endif

#ifdef CHEM
               dg%iota_IN  = dg%iota_IN  + dg%iota(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)
               dg%iota_EX  = dg%iota_EX  + dg%iota(K,EL_EX,dg%IRK)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)

               dg%iota2_IN = dg%iota2_IN + dg%iota2(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)
               dg%iota2_EX = dg%iota2_EX + dg%iota2(K,EL_EX,dg%IRK)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)
#endif

#ifdef dynp
               dynP_IN = dynP_IN + dg%dynP(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)
               dynP_EX = dynP_EX + dg%dynP(K,EL_EX,dg%IRK)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)
#endif

#ifdef SED_LAY
               do ll=1,s%layers
                  dg%bed_IN(ll) = dg%bed_IN(ll) + dg%bed(K,EL_IN,dg%IRK,ll)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)
                  dg%bed_EX(ll) = dg%bed_EX(ll) + dg%bed(K,EL_EX,dg%IRK,ll)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)
                  dg%HB_IN = dg%HB_IN + dg%bed(k,EL_IN,dg%irk,ll)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)
                  dg%HB_EX = dg%HB_EX + dg%bed(k,EL_EX,dg%irk,ll)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)

                  MZ_X_IN(ll) = MZ_X_IN(ll) + dg%MZ(K,1,ll,EL_IN)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)
                  MZ_Y_IN(ll) = MZ_Y_IN(ll) + dg%MZ(K,2,ll,EL_IN)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)  
                  
                  MZ_X_EX(ll) = MZ_X_EX(ll) + dg%MZ(K,1,ll,EL_EX)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)
                  MZ_Y_EX(ll) = MZ_Y_EX(ll) + dg%MZ(K,2,ll,EL_EX)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)
               enddo
#endif

                                ! LDG terms
#ifdef WAVE_DIF
               HZ_X_IN = HZ_X_IN &
              + dg%HZ(K,1,1,EL_IN)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)
               HZ_Y_IN = HZ_Y_IN &
              + dg%HZ(K,2,2,EL_IN)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)
               HZ_X_EX = HZ_X_EX &
              + dg%HZ(K,1,1,EL_EX)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)
               HZ_Y_EX = HZ_Y_EX &
              + dg%HZ(K,2,2,EL_EX)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)
#endif
               LZ_XX_IN = LZ_XX_IN &
              + dg%LZ(K,1,1,EL_IN)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)
               LZ_XY_IN = LZ_XY_IN &
              + dg%LZ(K,1,2,EL_IN)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)
               LZ_YX_IN = LZ_YX_IN &
              + dg%LZ(K,2,1,EL_IN)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)
               LZ_YY_IN = LZ_YY_IN &
              + dg%LZ(K,2,2,EL_IN)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)

               LZ_XX_EX = LZ_XX_EX &
              + dg%LZ(K,1,1,EL_EX)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)
               LZ_XY_EX = LZ_XY_EX &
              + dg%LZ(K,1,2,EL_EX)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)
               LZ_YX_EX = LZ_YX_EX &
              + dg%LZ(K,2,1,EL_EX)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)
               LZ_YY_EX = LZ_YY_EX &
              + dg%LZ(K,2,2,EL_EX)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)

#ifdef TRACE
               TZ_X_IN = TZ_X_IN &
              + dg%TZ(K,1,1,EL_IN)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)
               TZ_Y_IN = TZ_Y_IN &
              + dg%TZ(K,2,2,EL_IN)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)
               TZ_X_EX = TZ_X_EX &
              + dg%TZ(K,1,1,EL_EX)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)
               TZ_Y_EX = TZ_Y_EX &
              + dg%TZ(K,2,2,EL_EX)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)

#endif

            ENDDO

!.....Compute the numerical flux
            
            CALL NUMERICAL_FLUX(s,IT,dg%test_el)

!......... dummy
            F_HAT_O  = F_HAT

#ifdef TRACE
            I_HAT_O  = I_HAT
#endif

#ifdef CHEM
            I_HAT_O  = I_HAT
            J_HAT_O  = J_HAT
#endif

#ifdef dynp
            K_HAT_O  = K_HAT
#endif

#ifdef SED_LAY
            bed_HAT_IN(:) = dg%bed_HAT(:)
            bed_HAT_EX(:) = dg%bed_HAT(:)
#endif

            G_HAT_IN = G_HAT
            H_HAT_IN = H_HAT
            G_HAT_EX = G_HAT
            H_HAT_EX = H_HAT

!.....Check if the flux is large enough to dry up the elements
!.....1.01D0 is a safty factor.

            IF ( (1.01*F_HAT*XLEN_EL_IN*dg%MAX_BOA_DT(dg%IRK).GE.MASS_EL_IN).OR.&
           (1.01*F_HAT*XLEN_EL_EX*dg%MAX_BOA_DT(dg%IRK)*(-1.D0).GE.&
           MASS_EL_EX)) THEN

!...........Put back the saved RHS

               DO K = 1,dg%DOFS(dg%EL)
                  dg%RHS_ZE(K,EL_IN,dg%IRK) = RHS_ZE_IN(K)
                  dg%RHS_QX(K,EL_IN,dg%IRK) = RHS_QX_IN(K)
                  dg%RHS_QY(K,EL_IN,dg%IRK) = RHS_QY_IN(K)

                  dg%RHS_ZE(K,EL_EX,dg%IRK) = RHS_ZE_EX(K)
                  dg%RHS_QX(K,EL_EX,dg%IRK) = RHS_QX_EX(K)
                  dg%RHS_QY(K,EL_EX,dg%IRK) = RHS_QY_EX(K)

#ifdef TRACE
                  dg%RHS_iota(K,EL_IN,dg%IRK) = RHS_iota_IN(K)
                  dg%RHS_iota(K,EL_EX,dg%IRK) = RHS_iota_EX(K)
#endif

#ifdef CHEM
                  dg%RHS_iota(K,EL_IN,dg%IRK) = RHS_iota_IN(K)
                  dg%RHS_iota2(K,EL_IN,dg%IRK) = RHS_iota2_IN(K)
                  dg%RHS_iota(K,EL_EX,dg%IRK) = RHS_iota_EX(K)
                  dg%RHS_iota2(K,EL_EX,dg%IRK) = RHS_iota2_EX(K)
#endif

#ifdef dynp
                  dg%RHS_dynP(K,EL_IN,dg%IRK) = RHS_dynP_IN(K)
                  dg%RHS_dynP(K,EL_EX,dg%IRK) = RHS_dynP_EX(K)
#endif

#ifdef SED_LAY
                  do ll = 1,s%layers
                     dg%RHS_bed(K,EL_IN,dg%IRK,ll) = dg%RHS_bed_IN(K,ll)
                     dg%RHS_bed(K,EL_EX,dg%IRK,ll) = dg%RHS_bed_EX(K,ll)
                  enddo
#endif

               ENDDO

               GOTO 100
            ENDIF
            
!........Check to make sure mass flux is not coming from a dry element

            IF (abs(F_HAT).gt.1.d-12) THEN
               !
               IF (dg%WDFLG(EL_IN).EQ.0) THEN
                  ! EL_IN is dry !
                  IF (F_HAT.GT.0) THEN
                     ! Flux going from the dry element (in)
                     ! On the wet side (ex): reflect boundary
                     Q_N_EXT = dg%QX_EX*dg%NX + dg%QY_EX*dg%NY
                     Q_T_EXT = dg%QX_EX*TX + dg%QY_EX*TY
                     Q_N_INT = -Q_N_EXT
                     Q_T_INT =  Q_T_EXT
                     DEN = 1.D0/(dg%NX*TY - dg%NY*TX)
                     dg%QX_IN = ( TY*Q_N_INT - dg%NY*Q_T_INT)*DEN
                     dg%QY_IN = (-TX*Q_N_INT + dg%NX*Q_T_INT)*DEN
                     dg%ZE_IN = dg%ZE_EX

#ifdef TRACE
                     dg%iota_IN = dg%iota_EX
#endif

#ifdef CHEM
                     dg%iota_IN = dg%iota_EX
                     dg%iota2_IN = dg%iota2_EX
#endif

#ifdef dynp
                     dynP_IN = dynP_EX
#endif

#ifdef SED_LAY
                     dg%bed_IN(:) = dg%bed_EX(:)
#endif

                     dg%HB_IN = dg%HB_EX
                     dg%SFAC_IN = dg%SFAC_EX
                     CALL NUMERICAL_FLUX(s,IT,dg%test_el)

                     F_HAT_O  = F_HAT
                     G_HAT_EX = G_HAT
                     H_HAT_EX = H_HAT

                     ! on the dry side (in): do nothing
                     G_HAT_IN = G_HAT
                     H_HAT_IN = H_HAT
                     G_HAT_IN = 0.D0 
                     H_HAT_IN = 0.D0 

#ifdef TRACE
                     I_HAT_O  = I_HAT
#endif

#ifdef CHEM
                     I_HAT_O  = I_HAT
                     J_HAT_O  = J_HAT
#endif

#ifdef dynp
                     K_HAT_O  = K_HAT
#endif

#ifdef SED_LAY
                     bed_HAT_IN(:) = dg%bed_HAT(:)
                     bed_HAT_EX(:) = dg%bed_HAT(:)
#endif

                  ELSE
                     ! Flux coming from wet side (ex)
                     ! Leave fluxes on ex side intact and use zero gravtity flux for dry side      
                     NLEQG_TMP = NLEQG
                     NLEQG = 0.D0
                     G_TMP = G
                     G = 0.D0
                     CALL NUMERICAL_FLUX(s,IT,dg%test_el)
                     NLEQG = NLEQG_TMP
                     G = G_TMP
                     G_HAT_IN = G_HAT
                     H_HAT_IN = H_HAT
#ifdef SED_LAY
                     bed_HAT_IN(:) = dg%bed_HAT(:)
#endif
                  ENDIF
                  
               ELSEIF (dg%WDFLG(EL_EX).EQ.0) THEN
                  
                  ! EL_EX is dry
                  IF (F_HAT.LT.0) THEN
                     ! Flux comming from dry size (ex)
                     ! On the wet side (in): reflect boundary
                     Q_N_INT = dg%QX_IN*dg%NX + dg%QY_IN*dg%NY
                     Q_T_INT = dg%QX_IN*TX + dg%QY_IN*TY
                     Q_N_EXT = -Q_N_INT
                     Q_T_EXT =  Q_T_INT
                     DEN = 1.D0/(dg%NX*TY - dg%NY*TX)
                     dg%QX_EX = ( TY*Q_N_EXT - dg%NY*Q_T_EXT)*DEN
                     dg%QY_EX = (-TX*Q_N_EXT + dg%NX*Q_T_EXT)*DEN
                     dg%ZE_EX = dg%ZE_IN
                     dg%HB_EX = dg%HB_IN
                     dg%SFAC_EX = dg%SFAC_IN

#ifdef TRACE
                     dg%iota_EX = dg%iota_IN
#endif

#ifdef CHEM
                     dg%iota_EX = dg%iota_IN
                     dg%iota2_EX = dg%iota2_IN
#endif

#ifdef dynp
                     dynP_EX = dynP_IN
#endif

#ifdef SED_LAY
                     dg%bed_EX(:) = dg%bed_IN(:)
#endif

                     CALL NUMERICAL_FLUX(s,IT,dg%test_el)
                     F_HAT_O  = F_HAT

#ifdef TRACE
                     I_HAT_O  = I_HAT
#endif

#ifdef CHEM
                     I_HAT_O  = I_HAT
                     J_HAT_O  = J_HAT
#endif

#ifdef dynp
                     K_HAT_O  = K_HAT
#endif

#ifdef SED_LAY
                     bed_HAT_IN(:) = dg%bed_HAT(:)
                     bed_HAT_EX(:) = dg%bed_HAT(:)
#endif
                     
                     G_HAT_IN = G_HAT
                     H_HAT_IN = H_HAT

                     ! zero out mentum flux on the dry side (ex)
                     G_HAT_EX = G_HAT
                     H_HAT_EX = H_HAT
                     G_HAT_EX = 0.D0 
                     H_HAT_EX = 0.D0 
                  ELSE
                     ! Flux comming from the wet side (in)
                     ! Leave fluxes on in side intact and use zero gravtity for flux for dry (ex) side     
                     NLEQG_TMP = NLEQG
                     NLEQG = 0.D0
                     G_TMP = G
                     G = 0.D0
                     CALL NUMERICAL_FLUX(s,IT,dg%test_el)
                     NLEQG = NLEQG_TMP
                     G = G_TMP
                     G_HAT_EX = G_HAT
                     H_HAT_EX = H_HAT
#ifdef TRACE
                     I_HAT_O  = I_HAT
#endif
#ifdef CHEM
                     I_HAT_O  = I_HAT
                     J_HAT_O  = J_HAT
#endif

#ifdef SED_LAY
                     bed_HAT_EX(:) =  dg%bed_HAT(:)
#endif
                  ENDIF
               ENDIF
            ENDIF

                                !Let us add some diffusive penalties
#ifdef WAVE_DIF
            F_HAT_O = F_HAT_O + 0.5D0*(HZ_X_IN*dg%SFAC_IN+HZ_X_EX*dg%SFAC_EX)*dg%NX + &
           0.5D0*(HZ_Y_IN+HZ_Y_EX)*dg%NY
            !F_HAT_EX = F_HAT_EX + 0.5D0*(HZ_X_IN+HZ_X_EX)*dg%SFAC_IN*dg%NX + 
            !&           0.5D0*(HZ_Y_EX+HZ_Y_EX)*dg%NY
#endif
            G_HAT_IN = G_HAT_IN + 0.5D0*(LZ_XX_IN*dg%SFAC_IN + &
           LZ_XX_EX*dg%SFAC_EX)*dg%NX&
           + 0.5D0*(LZ_XY_IN + LZ_XY_EX)*dg%NY
!nd try adding the penalty term in the LDG method
      !$           +ESLM*.01D0*(dg%QX_IN-dg%QX_EX)
            H_HAT_IN = H_HAT_IN + 0.5D0*(LZ_YX_IN*dg%SFAC_IN + &
           LZ_YX_EX*dg%SFAC_EX)*dg%NX&
           + 0.5D0*(LZ_YY_IN + LZ_YY_EX)*dg%NY
!nd try adding the penalty term in the LDG method
      !$           +ESLM*.01D0*(dg%QY_IN-dg%QY_EX)
            G_HAT_EX = G_HAT_EX + 0.5D0*(LZ_XX_IN*dg%SFAC_IN&
           + LZ_XX_EX*dg%SFAC_EX)*dg%NX&
           + 0.5D0*(LZ_XY_IN + LZ_XY_EX)*dg%NY
!nd try adding the penalty term in the LDG method
      !$           +ESLM*.01D0*(dg%QX_IN-dg%QX_EX)
            H_HAT_EX = H_HAT_EX + 0.5D0*(LZ_YX_IN*dg%SFAC_IN&
           + LZ_YX_EX*dg%SFAC_EX)*dg%NX&
           + 0.5D0*(LZ_YY_IN + LZ_YY_EX)*dg%NY
!nd try adding the penalty term in the LDG method
      !$           +ESLM*.01D0*(dg%QY_IN-dg%QY_EX) 

#ifdef TRACE

            I_HAT_O = I_HAT_O + 0.5D0*(TZ_X_IN*dg%SFAC_IN+TZ_X_EX*dg%SFAC_EX)*dg%NX + &
           0.5D0*(TZ_Y_IN+TZ_Y_EX)*dg%NY
      !&           +ESLM*.01D0*(dg%iota_IN-dg%iota_EX)

#endif
     

                                !Let us add some sediment diffusion
#ifdef SED_LAY
            do ll=1,s%layers
               bed_HAT_IN(ll) = bed_HAT_IN(ll) + 0.5D0 * ( (MZ_X_IN(ll)+&
              MZ_X_EX(ll))*dg%SFAC_IN*dg%NX+(MZ_Y_IN(ll)+MZ_Y_EX(ll))*dg%NY)

               bed_HAT_EX(ll) = bed_HAT_EX(ll) + 0.5D0 * ( (MZ_X_IN(ll)+&
              MZ_X_EX(ll))*dg%SFAC_IN*dg%NX+(MZ_Y_IN(ll)+MZ_Y_EX(ll))*dg%NY)
            enddo
#endif

!.....Compute the edge integral

            DO K = 1,dg%DOFS(dg%EL)
               
               W_IN = EDFAC_IN*dg%EDGEQ(K,GP_IN,LED_IN,dg%pa)
               W_EX = EDFAC_EX*dg%EDGEQ(K,GP_EX,LED_EX,dg%pa)

               dg%RHS_ZE(K,EL_IN,dg%IRK) = dg%RHS_ZE(K,EL_IN,dg%IRK) - W_IN*F_HAT_O
               dg%RHS_QX(K,EL_IN,dg%IRK) = dg%RHS_QX(K,EL_IN,dg%IRK) - W_IN*G_HAT_IN
               dg%RHS_QY(K,EL_IN,dg%IRK) = dg%RHS_QY(K,EL_IN,dg%IRK) - W_IN*H_HAT_IN

               dg%RHS_ZE(K,EL_EX,dg%IRK) = dg%RHS_ZE(K,EL_EX,dg%IRK) + W_EX*F_HAT_O
               dg%RHS_QX(K,EL_EX,dg%IRK) = dg%RHS_QX(K,EL_EX,dg%IRK) + W_EX*G_HAT_EX
               dg%RHS_QY(K,EL_EX,dg%IRK) = dg%RHS_QY(K,EL_EX,dg%IRK) + W_EX*H_HAT_EX

#ifdef TRACE

               dg%RHS_iota(K,EL_IN,dg%IRK) = dg%RHS_iota(K,EL_IN,dg%IRK) - W_IN*I_HAT_O 
               dg%RHS_iota(K,EL_EX,dg%IRK) = dg%RHS_iota(K,EL_EX,dg%IRK) + W_EX*I_HAT_O
#endif

#ifdef CHEM
               dg%RHS_iota(K,EL_IN,dg%IRK)  = dg%RHS_iota(K,EL_IN,dg%IRK)  - W_IN*I_HAT_O 
               dg%RHS_iota2(K,EL_IN,dg%IRK) = dg%RHS_iota2(K,EL_IN,dg%IRK) - W_IN*J_HAT_O 
               dg%RHS_iota(K,EL_EX,dg%IRK)  = dg%RHS_iota(K,EL_EX,dg%IRK)  + W_EX*I_HAT_O
               dg%RHS_iota2(K,EL_EX,dg%IRK) = dg%RHS_iota2(K,EL_EX,dg%IRK) + W_EX*J_HAT_O
#endif

#ifdef dynp
               dg%RHS_dynP(K,EL_IN,dg%IRK) = dg%RHS_dynP(K,EL_IN,dg%IRK) - W_IN*K_HAT_O 
               dg%RHS_dynP(K,EL_EX,dg%IRK) = dg%RHS_dynP(K,EL_EX,dg%IRK) + W_EX*K_HAT_O
#endif


#ifdef SED_LAY
               do ll = 1,s%layers
                  dg%RHS_bed(K,EL_IN,dg%IRK,ll) = dg%RHS_bed(K,EL_IN,dg%IRK,ll) - W_IN*bed_HAT_IN(ll) 
                  dg%RHS_bed(K,EL_EX,dg%IRK,ll) = dg%RHS_bed(K,EL_EX,dg%IRK,ll) + W_EX*bed_HAT_EX(ll)
               enddo
#endif

            ENDDO

         ENDDO

         GOTO 1000

 100     CONTINUE

!--------------------------------------------------------------------
!     Compute reflection flux
!--------------------------------------------------------------------

!.....Set the components for the tangential vector to the edge

         TX = -dg%NY
         TY =  dg%NX

!.....Compute dg%ZE, dg%QX, dg%QY, dg%bed, and dg%HB at each edge Gauss quadrature point

         DO I = 1,dg%NEGP(dg%pa)

            GP_IN = I
            GP_EX = dg%NEGP(dg%pa) - I + 1

!     ----------------------------- EL_IN ------------------------------

            dg%ZE_IN = dg%ZE(1,EL_IN,dg%IRK)
            dg%QX_IN = dg%QX(1,EL_IN,dg%IRK)
            dg%QY_IN = dg%QY(1,EL_IN,dg%IRK)
            dg%HB_IN = dg%BATHED(GP_IN,LED_IN,EL_IN,dg%pa)

            dg%SFAC_IN = dg%SFACED(GP_IN,LED_IN,EL_IN,dg%pa)

#ifdef TRACE
            dg%iota_IN = dg%iota(1,EL_IN,dg%IRK)
#endif

#ifdef CHEM
            dg%iota_IN = dg%iota(1,EL_IN,dg%IRK)
            dg%iota2_IN = dg%iota2(1,EL_IN,dg%IRK)
#endif

#ifdef dynp
            dynP_IN = dg%dynP(1,EL_IN,dg%IRK)
#endif
            
            !Must add up the layers as they shift
#ifdef SED_LAY
            dg%HB(1,EL_IN,dg%irk) = 0.D0
            do ll=1,s%layers
               dg%bed_IN(ll) = dg%bed(1,EL_IN,dg%irk,ll)
               dg%HB(1,EL_IN,dg%irk) = dg%HB(1,EL_IN,dg%irk) + dg%bed(1,EL_IN,dg%irk,ll)

               !sediment diffusion
               MZ_X_IN(ll) = dg%MZ(1,1,ll,EL_IN)
               MZ_Y_IN(ll) = dg%MZ(1,2,ll,EL_IN)  
            enddo
            dg%HB_IN = dg%HB(1,EL_IN,dg%irk)
#endif

            !diffusion
#ifdef WAVE_DIF
            HZ_X_IN = dg%HZ(1,1,1,EL_IN)
            HZ_Y_IN = dg%HZ(1,2,2,EL_IN)
#endif
            LZ_XX_IN = dg%LZ(1,1,1,EL_IN)
            LZ_XY_IN = dg%LZ(1,1,2,EL_IN)
            LZ_YX_IN = dg%LZ(1,2,1,EL_IN)
            LZ_YY_IN = dg%LZ(1,2,2,EL_IN)

#ifdef TRACE

            TZ_X_IN = dg%TZ(1,1,1,EL_IN)
            TZ_Y_IN = dg%TZ(1,2,2,EL_IN)

#endif

            DO K = 2,dg%DOFS(dg%EL)

               dg%ZE_IN = dg%ZE_IN + dg%ZE(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)
               dg%QX_IN = dg%QX_IN + dg%QX(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)
               dg%QY_IN = dg%QY_IN + dg%QY(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)

#ifdef TRACE
               dg%iota_IN = dg%iota_IN + dg%iota(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%PA)
#endif

#ifdef CHEM
               dg%iota_IN  = dg%iota_IN  + dg%iota(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%PA)
               dg%iota2_IN = dg%iota2_IN + dg%iota2(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%PA)
#endif

#ifdef dynp
               dynP_IN = dynP_IN + dg%dynP(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%PA)
#endif

#ifdef SED_LAY
               do ll=1,s%layers
                  dg%bed_IN(ll) = dg%bed_IN(ll) + dg%bed(k,EL_IN,dg%irk,ll)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)
                  dg%HB_IN = dg%HB_IN + dg%bed(k,EL_IN,dg%irk,ll)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)
                  
                  MZ_X_IN(ll) = MZ_X_IN(ll) + dg%MZ(K,1,ll,EL_IN)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)
                  MZ_Y_IN(ll) = MZ_Y_IN(ll) + dg%MZ(K,2,ll,EL_IN)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)  
               enddo
#endif
                                ! LDG terms
#ifdef WAVE_DIF
               HZ_X_IN = HZ_X_IN&
              + dg%HZ(K,1,1,EL_IN)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)
               HZ_Y_IN = HZ_Y_IN&
              + dg%HZ(K,2,2,EL_IN)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)
#endif
               LZ_XX_IN = LZ_XX_IN&
              + dg%LZ(K,1,1,EL_IN)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)
               LZ_XY_IN = LZ_XY_IN&
              + dg%LZ(K,1,2,EL_IN)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)
               LZ_YX_IN = LZ_YX_IN&
              + dg%LZ(K,2,1,EL_IN)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)
               LZ_YY_IN = LZ_YY_IN&
              + dg%LZ(K,2,2,EL_IN)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)

#ifdef TRACE
               TZ_X_IN = TZ_X_IN&
              + dg%TZ(K,1,1,EL_IN)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)
               TZ_Y_IN = TZ_Y_IN&
              + dg%TZ(K,2,2,EL_IN)*dg%PHI_EDGE(K,GP_IN,LED_IN,dg%pa)
#endif

            ENDDO

!.....Compute the velocity in the normal and tangental direction

            Q_N_INT = dg%QX_IN*dg%NX + dg%QY_IN*dg%NY
            Q_T_INT = dg%QX_IN*TX + dg%QY_IN*TY

!.....Reflect the velocity in the normal direction

            Q_N_EXT = -Q_N_INT
            Q_T_EXT =  Q_T_INT

!.....Compute the x and y components of the external state flow

            DEN = 1.D0/(dg%NX*TY - dg%NY*TX)
            dg%QX_EX = ( TY*Q_N_EXT - dg%NY*Q_T_EXT)*DEN
            dg%QY_EX = (-TX*Q_N_EXT + dg%NX*Q_T_EXT)*DEN

            dg%ZE_EX = dg%ZE_IN
            dg%HB_EX = dg%HB_IN
            dg%SFAC_EX = dg%SFAC_IN

#ifdef TRACE
            dg%iota_EX = dg%iota_IN
#endif

#ifdef CHEM
            dg%iota_EX = dg%iota_IN
            dg%iota2_EX = dg%iota2_IN
#endif

#ifdef dynp
            dynP_EX = dynP_IN
#endif

#ifdef SED_LAY
            dg%bed_EX(:) = dg%bed_IN(:)
#endif

!.....Compute the numerical flux

!     CALL NUMERICAL_FLUX(IT)
            IF(dg%WDFLG(EL_IN).EQ.0) THEN
               NLEQG_TMP = NLEQG
               NLEQG = 0.D0
               G_TMP = G
               G = 0.D0
               CALL NUMERICAL_FLUX(s,IT,dg%test_el)
               NLEQG = NLEQG_TMP
               G = G_TMP
            ELSE
               CALL NUMERICAL_FLUX(s,IT,dg%test_el)
            ENDIF
            F_HAT = 0.D0        ! Ensure zero mass flux

!.....Add LDG terms for eddy viscosity

#ifdef WAVE_DIF
            F_HAT = F_HAT + HZ_X_IN*dg%NX*dg%SFAC_IN&
           + HZ_Y_IN*dg%NY
#endif
            G_HAT = G_HAT + LZ_XX_IN*dg%NX*dg%SFAC_IN&
           + LZ_XY_IN*dg%NY
            H_HAT = H_HAT + LZ_YX_IN*dg%NX*dg%SFAC_IN&
           + LZ_YY_IN*dg%NY

#ifdef TRACE

            I_HAT = I_HAT + TZ_X_IN*dg%NX*dg%SFAC_IN&
           + TZ_Y_IN*dg%NY

#endif

!.....Add LDG terms for sediment diffusion

#ifdef SED_LAY
            do ll=1,s%layers
               dg%bed_HAT(ll) = dg%bed_HAT(ll) + MZ_X_IN(ll)*dg%SFAC_IN*dg%NX &
              + MZ_Y_IN(ll)*dg%NY
            enddo
#endif

!.....Compute the edge integral
            DO K = 1,dg%DOFS(dg%EL)

               W_IN = EDFAC_IN*dg%EDGEQ(K,GP_IN,LED_IN,dg%pa)

               dg%RHS_ZE(K,EL_IN,dg%IRK) = dg%RHS_ZE(K,EL_IN,dg%IRK) - W_IN*F_HAT
               dg%RHS_QX(K,EL_IN,dg%IRK) = dg%RHS_QX(K,EL_IN,dg%IRK) - W_IN*G_HAT
               dg%RHS_QY(K,EL_IN,dg%IRK) = dg%RHS_QY(K,EL_IN,dg%IRK) - W_IN*H_HAT

#ifdef TRACE
               dg%RHS_iota(K,EL_IN,dg%IRK) = dg%RHS_iota(K,EL_IN,dg%IRK) - W_IN*I_HAT
#endif

#ifdef CHEM
               dg%RHS_iota(K,EL_IN,dg%IRK)  = dg%RHS_iota(K,EL_IN,dg%IRK)  - W_IN*I_HAT
               dg%RHS_iota2(K,EL_IN,dg%IRK) = dg%RHS_iota2(K,EL_IN,dg%IRK) - W_IN*J_HAT
#endif

#ifdef dynp
               dg%RHS_dynP(K,EL_IN,dg%IRK) = dg%RHS_dynP(K,EL_IN,dg%IRK) - W_IN*K_HAT
#endif

#ifdef SED_LAY
               do ll = 1,s%layers
                  dg%RHS_bed(K,EL_IN,dg%IRK,ll) = dg%RHS_bed(K,EL_IN,dg%IRK,ll) - W_IN*dg%bed_HAT(ll) 
               enddo
#endif

            ENDDO

!     ----------------------------- EL_EX ------------------------------

            dg%ZE_EX = dg%ZE(1,EL_EX,dg%IRK)
            dg%QX_EX = dg%QX(1,EL_EX,dg%IRK)
            dg%QY_EX = dg%QY(1,EL_EX,dg%IRK)
            dg%HB_EX = dg%BATHED(GP_EX,LED_EX,EL_EX,dg%pa)

            dg%SFAC_EX = dg%SFACED(GP_EX,LED_EX,EL_EX,dg%pa)

#ifdef TRACE
            dg%iota_EX = dg%iota(1,EL_EX,dg%IRK)
#endif

#ifdef CHEM
            dg%iota_EX = dg%iota(1,EL_EX,dg%IRK)
            dg%iota2_EX = dg%iota2(1,EL_EX,dg%IRK)
#endif

#ifdef dynp
            dynP_EX = dg%dynP(1,EL_EX,dg%IRK)
#endif

            !Must add up the layers as they shift
#ifdef SED_LAY
            dg%HB(1,EL_EX,dg%irk) = 0.D0
            do ll=1,s%layers
               dg%bed_EX(ll) = dg%bed(1,EL_EX,dg%irk,ll)
               dg%HB(1,EL_EX,dg%irk) = dg%HB(1,EL_EX,dg%irk) + dg%bed(1,EL_EX,dg%irk,ll)
               !sediment diffusion
               MZ_X_EX(ll) = dg%MZ(1,1,ll,EL_EX)
               MZ_Y_EX(ll) = dg%MZ(1,2,ll,EL_EX)  
            enddo
            dg%HB_EX = dg%HB(1,EL_EX,dg%irk)
#endif

            !LDG
#ifdef WAVE_DIF
            HZ_X_EX = dg%HZ(1,1,1,EL_EX)
            HZ_Y_EX = dg%HZ(1,2,2,EL_EX)
#endif
            LZ_XX_EX = dg%LZ(1,1,1,EL_EX)
            LZ_XY_EX = dg%LZ(1,1,2,EL_EX)
            LZ_YX_EX = dg%LZ(1,2,1,EL_EX)
            LZ_YY_EX = dg%LZ(1,2,2,EL_EX)

#ifdef TRACE

            TZ_X_EX = dg%TZ(1,1,1,EL_EX)
            TZ_Y_EX = dg%TZ(1,2,2,EL_EX)

#endif

            DO K = 2,dg%DOFS(dg%EL)

               dg%ZE_EX = dg%ZE_EX + dg%ZE(K,EL_EX,dg%IRK)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)
               dg%QX_EX = dg%QX_EX + dg%QX(K,EL_EX,dg%IRK)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)
               dg%QY_EX = dg%QY_EX + dg%QY(K,EL_EX,dg%IRK)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)

#ifdef TRACE
               dg%iota_EX = dg%iota_EX + dg%iota(K,EL_EX,dg%IRK)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)
#endif

#ifdef CHEM
               dg%iota_EX  = dg%iota_EX  + dg%iota(K,EL_EX,dg%IRK)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)
               dg%iota2_EX = dg%iota2_EX + dg%iota2(K,EL_EX,dg%IRK)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)
#endif

#ifdef dynp
               dynP_EX = dynP_EX + dg%dynP(K,EL_EX,dg%IRK)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)
#endif

#ifdef SED_LAY
               do ll=1,s%layers
                  dg%bed_EX(ll) = dg%bed_EX(ll) + dg%bed(k,EL_EX,dg%irk,ll)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)
                  dg%HB_EX = dg%HB_EX + dg%bed(k,EL_EX,dg%irk,ll)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)
                  
                  MZ_X_EX(ll) = MZ_X_EX(ll) + dg%MZ(K,1,ll,EL_EX)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)
                  MZ_Y_EX(ll) = MZ_Y_EX(ll) + dg%MZ(K,2,ll,EL_EX)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa) 
               enddo
#endif

                                ! LDG terms
#ifdef WAVE_DIF
               HZ_X_EX = HZ_X_EX&
              + dg%HZ(K,1,1,EL_EX)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)
               HZ_Y_EX = HZ_Y_EX&
              + dg%HZ(K,2,2,EL_EX)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)
#endif
               LZ_XX_EX = LZ_XX_EX&
              + dg%LZ(K,1,1,EL_EX)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)
               LZ_XY_EX = LZ_XY_EX&
              + dg%LZ(K,1,2,EL_EX)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)
               LZ_YX_EX = LZ_YX_EX&
              + dg%LZ(K,2,1,EL_EX)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)
               LZ_YY_EX = LZ_YY_EX&
              + dg%LZ(K,2,2,EL_EX)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)

#ifdef TRACE

               TZ_X_EX = TZ_X_EX&
              + dg%TZ(K,1,1,EL_EX)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)
               TZ_Y_EX = TZ_Y_EX&
              + dg%TZ(K,2,2,EL_EX)*dg%PHI_EDGE(K,GP_EX,LED_EX,dg%pa)

#endif

                                ! LDG terms for sediment diffusion               

            ENDDO

!.....Compute the velocity in the normal and tangental direction

            Q_N_EXT = dg%QX_EX*dg%NX + dg%QY_EX*dg%NY
            Q_T_EXT = dg%QX_EX*TX + dg%QY_EX*TY

!.....Reflect the velocity in the normal direction

            Q_N_INT = -Q_N_EXT
            Q_T_INT =  Q_T_EXT

!.....Compute the x and y components of the external state flow

            DEN = 1.D0/(dg%NX*TY - dg%NY*TX)
            dg%QX_IN = ( TY*Q_N_INT - dg%NY*Q_T_INT)*DEN
            dg%QY_IN = (-TX*Q_N_INT + dg%NX*Q_T_INT)*DEN

            dg%ZE_IN = dg%ZE_EX
            dg%HB_IN = dg%HB_EX

            dg%SFAC_IN = dg%SFAC_EX

#ifdef TRACE
            dg%iota_IN = dg%iota_EX
#endif

#ifdef CHEM
            dg%iota_IN = dg%iota_EX
            dg%iota2_IN = dg%iota2_EX
#endif

#ifdef dynp
            dynP_IN = dynP_EX
#endif

#ifdef SED_LAY
            dg%bed_IN(:) = dg%bed_EX(:)
#endif

!.....Compute the numerical flux

!     CALL NUMERICAL_FLUX(IT,L)
            IF(dg%WDFLG(EL_EX).EQ.0) THEN
               NLEQG_TMP = NLEQG
               NLEQG = 0.D0
               G_TMP = G
               G = 0.D0
               CALL NUMERICAL_FLUX(s,IT,dg%test_el)
               NLEQG = NLEQG_TMP
               G = G_TMP
            ELSE
               CALL NUMERICAL_FLUX(s,IT,dg%test_el)
            ENDIF
            F_HAT = 0.D0        ! Ensure zero mass flux

!.....Add LDG terms for eddy viscosity
#ifdef WAVE_DIF
            F_HAT = F_HAT + HZ_X_EX*dg%NX*dg%SFAC_EX&
           + HZ_Y_EX*dg%NY
#endif
            G_HAT = G_HAT + LZ_XX_EX*dg%NX*dg%SFAC_EX&
           + LZ_XY_EX*dg%NY
            H_HAT = H_HAT + LZ_YX_EX*dg%NX*dg%SFAC_EX&
           + LZ_YY_EX*dg%NY

#ifdef TRACE

            I_HAT = I_HAT + TZ_X_EX*dg%NX*dg%SFAC_EX&
           + TZ_Y_EX*dg%NY

#endif

!.....Add LDG terms for sediment diffusion

#ifdef SED_LAY
            do ll=1,s%layers
               dg%bed_HAT(ll) = dg%bed_HAT(ll) + MZ_X_EX(ll)*dg%SFAC_EX*dg%NX &
              + MZ_Y_EX(ll)*dg%NY
            enddo
#endif

!.....Compute the edge integral
            DO K = 1,dg%DOFS(dg%EL)

               W_EX = EDFAC_EX*dg%EDGEQ(K,GP_EX,LED_EX,dg%pa)

               dg%RHS_ZE(K,EL_EX,dg%IRK) = dg%RHS_ZE(K,EL_EX,dg%IRK) + W_EX*F_HAT
               dg%RHS_QX(K,EL_EX,dg%IRK) = dg%RHS_QX(K,EL_EX,dg%IRK) + W_EX*G_HAT
               dg%RHS_QY(K,EL_EX,dg%IRK) = dg%RHS_QY(K,EL_EX,dg%IRK) + W_EX*H_HAT

#ifdef TRACE
               dg%RHS_iota(K,EL_EX,dg%IRK) = dg%RHS_iota(K,EL_EX,dg%IRK) + W_EX*I_HAT
#endif

#ifdef CHEM
               dg%RHS_iota(K,EL_EX,dg%IRK)  = dg%RHS_iota(K,EL_EX,dg%IRK)  + W_EX*I_HAT
               dg%RHS_iota2(K,EL_EX,dg%IRK) = dg%RHS_iota2(K,EL_EX,dg%IRK) + W_EX*J_HAT
#endif

#ifdef dynp
               dg%RHS_dynP(K,EL_EX,dg%IRK) = dg%RHS_dynP(K,EL_EX,dg%IRK) + W_EX*K_HAT
#endif

#ifdef SED_LAY
               do ll = 1,s%layers
                  dg%RHS_bed(K,EL_EX,dg%IRK,ll) = dg%RHS_bed(K,EL_EX,dg%IRK,ll) + W_EX*dg%bed_HAT(ll) 
               enddo
#endif

            ENDDO

         ENDDO
         
 1000 CONTINUE

      RETURN
      END SUBROUTINE
