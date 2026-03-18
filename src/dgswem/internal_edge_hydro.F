C***********************************************************************
C
C     SUBROUTINE INTERNAL_EDGE_HYDRO( )
C
C     This subroutine does the following:
C
C     1.  Calculates the values of the necessary variables at the edge
C     gauss points for INTERNAL edges
C     2.  Calls the appropriate subroutine to compute the flux at
C     these points.
C     3.  Calls the appropriate subroutine to compute the boundary
C     integrals.
C
C     Written by Ethan Kubatko (06-11-2004)
C
C     01-10-2011 - cem - adapted for p_enrichment and multicomponent
C     11-11-2011 - cem - adapted for layered sediment
C
C-----------------------------------------------------------------------
C
C     01-02-2007, sb, Modified for LDG
C     C***********************************************************************

      module mod_internal_edge

      USE GLOBAL, only : areas, pdg_el, neled, ne, G
      USE DG, only : rhs_ze, rhs_qx, rhs_qy, cosnx, sinnx, ze, wdflg, qx, qy, xlen,
     &     hb, bathed, dofs, edgeq, negp, max_boa_dt, phi_edge, sfaced, lz, nedel, niedn,
     &     nedsd, elem_edge, ze_edge, qx_edge, qy_edge,
     &     irk, nieds, dofh
      USE NodalAttributes, ONLY :  ESLM
      use sizes, only:  sz, myproc
      use mod_numerical_flux

      contains

      SUBROUTINE INTERNAL_EDGE_HYDRO(IT)


      IMPLICIT NONE

      integer, value :: it

C.....Declare local variables
      real(sz) :: tmp, fluxes(3)
      logical :: is_large_flux
      integer :: el_in, el_ex, pa, el, test_el, ind_ex, ind_in
      real(sz) :: nx, ny, ze_in, qx_in, qy_in, qx_ex, qy_ex, hb_ex, hb_in
      real(sz) :: f_hat, g_hat, h_hat, q_n_ext, q_t_ext, q_n_int, q_t_int, ze_ex
      real(sz) :: sfac_in, sfac_ex
      INTEGER L, LED_IN, LED_EX, GED, GP_IN, GP_EX,k,i,ll
      REAL(SZ), PARAMETER :: ZERO = 1.D-12
      REAL(SZ) TX, TY, W_IN, W_EX
      REAL(SZ) LZ_XX_IN, LZ_XY_IN, LZ_YX_IN, LZ_YY_IN
      REAL(SZ) LZ_XX_EX, LZ_XY_EX, LZ_YX_EX, LZ_YY_EX
      REAL(SZ) HZ_X_EX, HZ_Y_EX, HZ_X_IN, HZ_Y_IN
      REAL(SZ) TZ_X_EX, TZ_Y_EX, TZ_X_IN, TZ_Y_IN
      REAL(SZ) EDFAC_IN, EDFAC_EX, DEN
      real(sz) :: ze_ins(3), ze_exs(3), qx_ins(3), qx_exs(3), qy_ins(3), qy_exs(3)
      REAL(SZ) XLEN_EL_IN, XLEN_EL_EX
      REAL(SZ) MASS_EL_IN, MASS_EL_EX
      REAL(SZ)  ::
     &     RHS_ZE_IN(3), RHS_QX_IN(3), RHS_QY_IN(3),
     &     RHS_ZE_EX(3), RHS_QX_EX(3), RHS_QY_EX(3)
#ifdef TRACE
     &     ,RHS_iota_IN(:), RHS_iota_EX(:)
#endif
#ifdef CHEM
     &     ,RHS_iota_IN(:), RHS_iota_EX(:)
     &     ,RHS_iota2_IN(:), RHS_iota2_EX(:)
#endif
#ifdef DYNP
     &     ,RHS_dynP_IN(:), RHS_dynP_EX(:)
#endif

      REAL(SZ) ARK, BRK
      REAL(SZ) MAX_BOA          ! Maximum of beta_il/alpha_il for all l
      REAL(SZ) F_HAT_O, G_HAT_O, H_HAT_O,i_hat_o,j_hat_o
      REAL(SZ) G_HAT_IN, H_HAT_IN
      REAL(SZ) G_HAT_EX, H_HAT_EX
      REAL(SZ) K_HAT_O

      test_el = 0
!$acc parallel loop gang vector
!$acc+private(rhs_ze_in,rhs_qx_in,rhs_qy_in,rhs_ze_ex,rhs_qx_ex,rhs_qy_ex,
!$acc+ze_ins,qx_ins,qy_ins,ze_exs,qx_exs,qy_exs,fluxes) async(1) default(present)
      do L = 1, NIEDS
!$acc loop seq
        do k = 1,DOFH
          rhs_ze_in(k) = 0.0_sz
          rhs_qx_in(k) = 0.0_sz
          rhs_qy_in(k) = 0.0_sz
          rhs_ze_ex(k) = 0.0_sz
          rhs_qx_ex(k) = 0.0_sz
          rhs_qy_ex(k) = 0.0_sz
        enddo

C.......Retrieve the global and local edge number

         GED = NIEDN(L)
         LED_IN = NEDSD(1,GED)
         LED_EX = NEDSD(2,GED)

C.......Retrieve the elements which share the edge

         EL_IN = NEDEL(1,GED)
         EL_EX = NEDEL(2,GED)

         EL = EL_EX

         if (DOFS(EL_EX).LT.DOFS(EL_IN)) then
            EL = EL_IN
         endif

         pa = PDG_EL(EL)


C.......If both elements on either side of edge are dry then skip

         IF((WDFLG(EL_IN).EQ.0).AND.(WDFLG(EL_EX).EQ.0)) goto 100

         ze_ins = ze(:,el_in,irk)
         qx_ins = qx(:,el_in,irk)
         qy_ins = qy(:,el_in,irk)
         ze_exs = ze(:,el_ex,irk)
         qx_exs = qx(:,el_ex,irk)
         qy_exs = qy(:,el_ex,irk)



C.....Compute the sum of the lengths of three edges

         XLEN_EL_IN = XLEN(NELED(1,EL_IN))
         XLEN_EL_IN = XLEN_EL_IN + XLEN(NELED(2,EL_IN))
         XLEN_EL_IN = XLEN_EL_IN + XLEN(NELED(3,EL_IN))

         XLEN_EL_EX = XLEN(NELED(1,EL_EX))
         XLEN_EL_EX = XLEN_EL_EX + XLEN(NELED(2,EL_EX))
         XLEN_EL_EX = XLEN_EL_EX + XLEN(NELED(3,EL_EX))

C.....Compute the total mass in the elements

         MASS_EL_IN = (ZE(1,EL_IN,IRK)+HB(1,EL_IN,1))*AREAS(EL_IN)*0.50
         MASS_EL_EX = (ZE(1,EL_EX,IRK)+HB(1,EL_EX,1))*AREAS(EL_EX)*0.50



C.....Retrieve the components of the normal vector to the edge

         NX = COSNX(GED)
         NY = SINNX(GED)

C.....Set the components for the tangential vector to the edge

         TX = -NY
         TY =  NX

         EDFAC_IN = XLEN(GED)/AREAS(EL_IN)
         EDFAC_EX = XLEN(GED)/AREAS(EL_EX)
C.....Compute ZE, QX, QY, and HB at each edge Gauss quadrature point

!$acc loop seq
         gauss_point_loop: DO I = 1,NEGP(pa)

            GP_IN = I
            GP_EX = NEGP(pa) - I + 1

            ZE_IN = 0.0_sz
            QX_IN = 0.0_sz
            QY_IN = 0.0_sz
            HB_IN = BATHED(GP_IN,LED_IN,EL_IN,pa)
            SFAC_IN = SFACED(GP_IN,LED_IN,EL_IN,pa)

            ZE_EX = 0.0_sz
            QX_EX = 0.0_sz
            QY_EX = 0.0_sz
            HB_EX = BATHED(GP_EX,LED_EX,EL_EX,pa)
            SFAC_EX = SFACED(GP_EX,LED_EX,EL_EX,pa)

            DO K = 1,DOFH

               ZE_IN = ZE_IN + ZE_ins(K)*PHI_EDGE(K,GP_IN,LED_IN,pa)
               QX_IN = QX_IN + QX_ins(K)*PHI_EDGE(K,GP_IN,LED_IN,pa)
               QY_IN = QY_IN + QY_ins(K)*PHI_EDGE(K,GP_IN,LED_IN,pa)

               ZE_EX = ZE_EX + ZE_exs(K)*PHI_EDGE(K,GP_EX,LED_EX,pa)
               QX_EX = QX_EX + QX_exs(K)*PHI_EDGE(K,GP_EX,LED_EX,pa)
               QY_EX = QY_EX + QY_exs(K)*PHI_EDGE(K,GP_EX,LED_EX,pa)


            ENDDO

C.....Compute the numerical flux

            !CALL NUMERICAL_FLUX(IT,test_el)
            fluxes = numerical_flux(ze_in,ze_ex,qx_in,qx_ex,qy_in,qy_ex,
     &  hb_in,hb_ex,sfac_in,sfac_ex,nx,ny,G)
            f_hat = fluxes(1)
            g_hat = fluxes(2)
            h_hat = fluxes(3)

C......... dummy
            F_HAT_O  = F_HAT

            G_HAT_IN = G_HAT
            H_HAT_IN = H_HAT
            G_HAT_EX = G_HAT
            H_HAT_EX = H_HAT

C........Check to make sure mass flux is not coming from a dry element

            dry_flux: IF (abs(F_HAT).gt.1.d-8) THEN
               !
               IF (WDFLG(EL_IN).EQ.0) THEN
                  ! EL_IN is dry !
                  IF (F_HAT.GT.0.0_sz) THEN
                     ! Flux going from the dry element (in)
                     ! On the wet side (ex): reflect boundary
                     Q_N_EXT = QX_EX*NX + QY_EX*NY
                     Q_T_EXT = QX_EX*TX + QY_EX*TY
                     Q_N_INT = -Q_N_EXT
                     Q_T_INT =  Q_T_EXT
                     DEN = (NX*TY - NY*TX)
                     QX_IN = ( TY*Q_N_INT - NY*Q_T_INT)/DEN
                     QY_IN = (-TX*Q_N_INT + NX*Q_T_INT)/DEN
                     ZE_IN = ZE_EX


                     HB_IN = HB_EX
                     SFAC_IN = SFAC_EX

            fluxes = numerical_flux(ze_in,ze_ex,qx_in,qx_ex,qy_in,qy_ex,
     &  hb_in,hb_ex,sfac_in,sfac_ex,nx,ny,G)
            f_hat = fluxes(1)
            g_hat = fluxes(2)
            h_hat = fluxes(3)
                     F_HAT_O  = F_HAT
                     G_HAT_EX = G_HAT
                     H_HAT_EX = H_HAT

                     ! on the dry side (in): do nothing

                     G_HAT_IN = 0.0_sz
                     H_HAT_IN = 0.0_sz


                  ELSE
                     ! Flux coming from wet side (ex)
                     ! Leave fluxes on ex side intact and use zero gravity flux for dry side
            fluxes = numerical_flux(ze_in,ze_ex,qx_in,qx_ex,qy_in,qy_ex,
     &  hb_in,hb_ex,sfac_in,sfac_ex,nx,ny,0.0_sz)
            f_hat = fluxes(1)
            g_hat = fluxes(2)
            h_hat = fluxes(3)
                     G_HAT_IN = G_HAT
                     H_HAT_IN = H_HAT
                  ENDIF

               ELSEIF (WDFLG(EL_EX).EQ.0) THEN

                  ! EL_EX is dry
                  IF (F_HAT.LT.0.0_sz) THEN
                     ! Flux comming from dry size (ex)
                     ! On the wet side (in): reflect boundary
                     Q_N_INT = QX_IN*NX + QY_IN*NY
                     Q_T_INT = QX_IN*TX + QY_IN*TY
                     Q_N_EXT = -Q_N_INT
                     Q_T_EXT =  Q_T_INT
                     DEN = (NX*TY - NY*TX)
                     QX_EX = ( TY*Q_N_EXT - NY*Q_T_EXT)/DEN
                     QY_EX = (-TX*Q_N_EXT + NX*Q_T_EXT)/DEN
                     ZE_EX = ZE_IN
                     HB_EX = HB_IN
                     SFAC_EX = SFAC_IN

            fluxes = numerical_flux(ze_in,ze_ex,qx_in,qx_ex,qy_in,qy_ex,
     &  hb_in,hb_ex,sfac_in,sfac_ex,nx,ny,G)
            f_hat = fluxes(1)
            g_hat = fluxes(2)
            h_hat = fluxes(3)
                     F_HAT_O  = F_HAT

                     G_HAT_IN = G_HAT
                     H_HAT_IN = H_HAT

                     ! zero out momentum flux on the dry side (ex)
                     G_HAT_EX = 0.0_sz
                     H_HAT_EX = 0.0_sz
                  ELSE
                     ! Flux comming from the wet side (in)
                     ! Leave fluxes on in side intact and use zero gravity for flux for dry (ex) side
                     !CALL NUMERICAL_FLUX(IT,test_el)
            fluxes = numerical_flux(ze_in,ze_ex,qx_in,qx_ex,qy_in,qy_ex,
     &  hb_in,hb_ex,sfac_in,sfac_ex,nx,ny,0.0_sz)
            f_hat = fluxes(1)
            g_hat = fluxes(2)
            h_hat = fluxes(3)
                     G_HAT_EX = G_HAT
                     H_HAT_EX = H_HAT
                  ENDIF
               ENDIF
            ENDIF dry_flux

C.....Compute the edge integral

            DO K = 1,3

               W_IN = EDFAC_IN*EDGEQ(K,GP_IN,LED_IN,pa)
               W_EX = EDFAC_EX*EDGEQ(K,GP_EX,LED_EX,pa)

               rhs_ze_in(k) = rhs_ze_in(k) - w_in * f_hat_o
               rhs_qx_in(k) = rhs_qx_in(k) - w_in * g_hat_in
               rhs_qy_in(k) = rhs_qy_in(k) - w_in * h_hat_in

               rhs_ze_ex(k) = rhs_ze_ex(k) + w_ex * f_hat_o
               rhs_qx_ex(k) = rhs_qx_ex(k) + w_ex * g_hat_ex
               rhs_qy_ex(k) = rhs_qy_ex(k) + w_ex * h_hat_ex
            ENDDO


      ENDDO gauss_point_loop



      ind_ex = -1
      ind_in = -1
      do k = 1,3
         if (elem_edge(el_ex, k) == ged) then
            ind_ex = k
         endif
         if (elem_edge(el_in, k) == ged) then
            ind_in = k
         endif
      enddo


      do k = 1,3
         ze_edge(el_ex,ind_ex,k) = rhs_ze_ex(k)
         qx_edge(el_ex,ind_ex,k) = rhs_qx_ex(k)
         qy_edge(el_ex,ind_ex,k) = rhs_qy_ex(k)

         ze_edge(el_in,ind_in,k) = rhs_ze_in(k)
         qx_edge(el_in,ind_in,k) = rhs_qx_in(k)
         qy_edge(el_in,ind_in,k) = rhs_qy_in(k)
      enddo

100   continue
      enddo

      END SUBROUTINE

      subroutine add_internal_flux()
      implicit none


      real(sz) :: ze_sum, qx_sum, qy_sum
      integer :: l, k, i
!$acc parallel loop gang vector async(1) default(present)
      do l = 1,NE
         do k = 1,3
            ze_sum = 0.
            qx_sum = 0.
            qy_sum = 0.
            do i = 1,3
               ze_sum = ze_sum + ze_edge(L,i,k)
               qx_sum = qx_sum + qx_edge(L,i,k)
               qy_sum = qy_sum + qy_edge(L,i,k)
            enddo

            rhs_ze(k,L, IRK) = rhs_ze(k,L, IRK) + ze_sum
            rhs_qx(k,L, IRK) = rhs_qx(k,L, IRK) + qx_sum
            rhs_qy(k,L, IRK) = rhs_qy(k,L, IRK) + qy_sum
         enddo

         ze_edge(L,:,:) = 0.
         qx_edge(L,:,:) = 0.
         qy_edge(L,:,:) = 0.

      enddo



      end subroutine

      end module mod_internal_edge
