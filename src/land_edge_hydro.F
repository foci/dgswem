C***********************************************************************
C     
C     SUBROUTINE LAND_EDGE_HYDRO( )
C     
C     This subroutine does the following:
C     
C     1.  Calculates the values of the necessary variables at the edge
C     gauss points for NO-NORMAL FLOW edges
C     2.  Calls the appropriate subroutine to compute the flux at
C     these points.
C     3.  Calls the appropriate subroutine to compute the boundary
C     integrals.
C     
C     Written by Ethan Kubatko (06-11-2004)
C     
C-----------------------------------------------------------------------
C     
C     01-02-2007, sb, Modified for LDG
C     08-xx-2005, sb, Modifications for wetting/drying
C     
C     01-10-2011 - cem - adapted for p_enrichment and multicomponent
C     06-02-2012 - cem - adapted for sediment
C     
C***********************************************************************

      SUBROUTINE LAND_EDGE_HYDRO(IT)

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      Use SIZES, only: layers
      use mod_numerical_flux

      IMPLICIT NONE

C.....Declare local variables

      real(sz) :: ze_in,ze_ex,qx_in,qx_ex,qy_in,qy_ex,
     &  hb_in,hb_ex,sfac_in,sfac_ex, nx, ny
      real(sz) :: f_hat, g_hat, h_hat
      INTEGER L, LED, GED, i,k,jj,II,ll,IT,mm
      Real(SZ) DEN2,U_AVG,V_AVG,VEL_NORMAL,q_RoeX, q_RoeY, q_Roe
      REAL(SZ) TX, TY, DEN,ell_1,ell_2,ell_3,HZ_X_IN,HZ_Y_IN
      REAL(SZ) LZ_XX_IN, LZ_XY_IN, LZ_YX_IN, LZ_YY_IN,W_IN
      Real(SZ) MZ_X_IN(layers),MZ_Y_IN(layers),TZ_X_IN,TZ_Y_IN

!$acc parallel loop gang vector async(1) default(present)
      DO 1000 L = 1,NLEDS
         
C.....Retrieve the global and local edge number

         GED = NLEDN(L)
         LED = NEDSD(1,GED)

C.....Retrieve the elements which share the edge

         EL_IN = NEDEL(1,GED)

         pa = PDG_EL(EL_IN)

#ifdef P0
         if (pa.eq.0) then
            pa = 1
         endif
#endif
         
C.....If the element is dry then skip the edge calculation

         IF (WDFLG(EL_IN).EQ.0) GOTO 1000


C.....Retrieve the components of the normal vector to the edge
         
         NX = COSNX(GED)
         NY = SINNX(GED)
         
C.....Set the components for the tangential vector to the edge
         
         TX = -NY
         TY =  NX

C.....Compute ZE, QX, QY, and HB at each Gauss point

!$acc loop seq
         DO I = 1,NEGP(pa)
            
            ZE_IN = ZE(1,EL_IN,IRK)
            QX_IN = QX(1,EL_IN,IRK)
            QY_IN = QY(1,EL_IN,IRK)
            HB_IN = BATHED(I,LED,EL_IN,pa)

            SFAC_IN = SFACED(I,LED,EL_IN,pa)

#ifdef TRACE
            iota_IN = iota(1,EL_IN,IRK)
#endif

#ifdef CHEM
            iota_IN = 0.0 !iota(1,EL_IN,IRK)
            iota2_IN = 0.0 !iota2(1,EL_IN,IRK)
#endif

#ifdef DYNP
            dynP_IN = dynP(1,EL_IN,IRK)
#endif

#ifdef SED_LAY                  !When layered, these change
            HB(1,EL_IN,irk) = 0.0
            do ll=1,layers
               HB(1,EL_IN,irk) = HB(1,EL_IN,irk) + bed(1,EL_IN,irk,ll)

               MZ_X_IN(ll) =  MZ(1,1,ll,EL_IN)
               MZ_Y_IN(ll) =  MZ(1,2,ll,EL_IN)
            enddo
            bed_IN(:) = bed(1,EL_IN,irk,:)
            HB_IN = HB(1,EL_IN,irk)

#endif

            !...do it and do it again ... 
#ifdef WAVE_DIF
            HZ_X_IN = HZ(1,1,1,EL_IN)
            HZ_Y_IN = HZ(1,2,2,EL_IN)
#endif
            LZ_XX_IN = LZ(1,1,1,EL_IN)
            LZ_XY_IN = LZ(1,1,2,EL_IN)
            LZ_YX_IN = LZ(1,2,1,EL_IN)
            LZ_YY_IN = LZ(1,2,2,EL_IN)

#ifdef TRACE

            TZ_X_IN = TZ(1,1,1,EL_IN)
            TZ_Y_IN = TZ(1,2,2,EL_IN)

#endif

            !Compute sediment diffusion contribution

C.....Compute the solution at the interior state
C.....(modified for wetting and drying)

!$acc loop seq
            DO K = 2,DOFS(EL_IN)
               
               ZE_IN = ZE_IN + ZE(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa)
               QX_IN = QX_IN + QX(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa)
               QY_IN = QY_IN + QY(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa)

                                !LDG terms
#ifdef WAVE_DIF
               HZ_X_IN = HZ_X_IN + HZ(K,1,1,EL_IN)*PHI_EDGE(K,I,LED,pa)
               HZ_Y_IN = HZ_Y_IN + HZ(K,2,2,EL_IN)*PHI_EDGE(K,I,LED,pa)
#endif
               LZ_XX_IN = LZ_XX_IN + LZ(K,1,1,EL_IN)*PHI_EDGE(K,I,LED,pa)
               LZ_XY_IN = LZ_XY_IN + LZ(K,1,2,EL_IN)*PHI_EDGE(K,I,LED,pa)
               LZ_YX_IN = LZ_YX_IN + LZ(K,2,1,EL_IN)*PHI_EDGE(K,I,LED,pa)
               LZ_YY_IN = LZ_YY_IN + LZ(K,2,2,EL_IN)*PHI_EDGE(K,I,LED,pa)

#ifdef TRACE

               TZ_X_IN = TZ_X_IN + TZ(K,1,1,EL_IN)*PHI_EDGE(K,I,LED,pa)
               TZ_Y_IN = TZ_Y_IN + TZ(K,2,2,EL_IN)*PHI_EDGE(K,I,LED,pa)

#endif

                                !LDG terms for sediment diffusion

#ifdef TRACE
               iota_IN = iota_IN + iota(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa)
#endif

#ifdef CHEM
               iota_IN = iota_IN + iota(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa)
               iota2_IN = iota2_IN + iota2(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa)
#endif

#ifdef DYNP
               dynP_IN = dynP_IN + dynP(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa)
#endif

#ifdef SED_LAY
               do ll = 1,layers
                  bed_IN(ll) = bed_IN(ll) +
     &                  bed(K,EL_IN,IRK,ll)*PHI_EDGE(K,I,LED,pa)
                  HB_IN=HB_IN + bed(k,EL_IN,irk,ll)*PHI_EDGE(K,I,LED,pa)

                  MZ_X_IN(ll) = MZ_X_IN(ll) +
     &                  MZ(K,1,ll,EL_IN)*PHI_EDGE(K,I,LED,pa)
                  MZ_Y_IN(ll) = MZ_Y_IN(ll) +
     &                  MZ(K,2,ll,EL_IN)*PHI_EDGE(K,I,LED,pa)

               enddo
#endif

            ENDDO


C.....Compute the velocity in the normal and tangental direction
      
            Q_N_INT = QX_IN*NX + QY_IN*NY
            Q_T_INT = QX_IN*TX + QY_IN*TY

C.....Reflect the velocity in the normal direction

            Q_N_EXT = -Q_N_INT
            Q_T_EXT =  Q_T_INT
            
C.....Compute the x and y components of the external state flow

 
            DEN = 1.0/(NX*TY - NY*TX)
            QX_EX = ( TY*Q_N_EXT - NY*Q_T_EXT)*DEN
            QY_EX = (-TX*Q_N_EXT + NX*Q_T_EXT)*DEN

            ZE_EX = ZE_IN
            HB_EX = HB_IN
            SFAC_EX = SFAC_IN

#ifdef TRACE
            iota_EX = 0.0 !iota_IN
            !print*, 'test'
#endif

#ifdef CHEM
            iota_EX = iota_IN
            iota2_EX = iota2_IN
#endif

#ifdef DYNP
            dynP_EX = dynP_IN
#endif

#ifdef SED_LAY
            bed_EX(:) = bed_IN(:)
#endif

C.....Compute the Roe flux

            call numerical_flux(ze_in,ze_ex,qx_in,qx_ex,qy_in,qy_ex,
     &  hb_in,hb_ex,sfac_in,sfac_ex,nx,ny,f_hat,g_hat,h_hat)

C.....Add LDG terms for viscosity

#ifdef WAVE_DIF
            F_HAT = F_HAT + HZ_X_IN*NX*SFAC_IN + HZ_Y_IN*NY
#endif
            G_HAT = G_HAT + LZ_XX_IN*NX*SFAC_IN + LZ_XY_IN*NY
            H_HAT = H_HAT + LZ_YX_IN*NX*SFAC_IN + LZ_YY_IN*NY

#ifdef TRACE
            I_HAT = I_HAT + TZ_X_IN*NX*SFAC_IN + TZ_Y_IN*NY
#endif

C.....Add LDG terms for sediment

#ifdef SED_LAY
            do ll=1,layers
               bed_HAT(ll) = bed_HAT(ll) +
     &                       MZ_X_IN(ll)*NX*SFAC_IN + MZ_Y_IN(ll)*NY
            enddo
#endif
            
C.....Compute the edge integral
C.....(modified for wetting and drying) 


C.....Compute the edge integral
!$acc loop seq
            DO K = 1,DOFS(EL_IN)

               W_IN = 2.0*M_INV(K,pa)/AREAS(EL_IN)*XLEN(GED)*
     &              PHI_EDGE(K,I,LED,pa)*WEGP(I,pa)
!$acc atomic
               RHS_ZE(K,EL_IN,IRK) = RHS_ZE(K,EL_IN,IRK) - W_IN*F_HAT
!$acc atomic
               RHS_QX(K,EL_IN,IRK) = RHS_QX(K,EL_IN,IRK) - W_IN*G_HAT
!$acc atomic
               RHS_QY(K,EL_IN,IRK) = RHS_QY(K,EL_IN,IRK) - W_IN*H_HAT

#ifdef TRACE
               RHS_iota(K,EL_IN,IRK) = RHS_iota(K,EL_IN,IRK) - W_IN*I_HAT
#endif

#ifdef CHEM
               RHS_iota(K,EL_IN,IRK) = RHS_iota(K,EL_IN,IRK) - W_IN*I_HAT
               RHS_iota2(K,EL_IN,IRK) = RHS_iota2(K,EL_IN,IRK) - W_IN*J_HAT
#endif

#ifdef DYNP
               RHS_dynP(K,EL_IN,IRK) = RHS_dynP(K,EL_IN,IRK) - W_IN*K_HAT
#endif

#ifdef SED_LAY
               do ll = 1,layers
                  RHS_bed(K,EL_IN,IRK,ll) = RHS_bed(K,EL_IN,IRK,ll) -
     &                                      W_IN*bed_HAT(ll)
               enddo
#endif

            ENDDO

         ENDDO

 1000 CONTINUE
      
      RETURN
      END SUBROUTINE
