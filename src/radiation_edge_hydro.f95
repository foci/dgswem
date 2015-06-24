!***********************************************************************
!     
!     SUBROUTINE RADIATION_EDGE_HYDRO( )
!     
!     This subroutine does the following:
!     
!     1.  Calculates the values of the necessary variables at the edge
!     gauss points for RADIATION edges
!     2.  Calls the appropriate subroutine to compute the flux at
!     these points.
!     3.  Compute the boundary integrals.
!     
!     Written by Ethan Kubatko (06-11-2004)
!     01-10-2011 - cem - adapted for p_enrichment and multicomponent
!     You know, this file is not really up to date.  
!     If you need to use it, you might fix it first -- cem
!     
!***********************************************************************

      SUBROUTINE RADIATION_EDGE_HYDRO(s,dg,IT)

!.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE sizes
      IMPLICIT NONE
      
      type (sizes_type) :: s
      type (dg_type) :: dg

!.....Declare local variables

      INTEGER L, LED, GED,i,k,ll,IT
      REAL(SZ) TX, TY, W_IN, DEN
      REAL(SZ) LZ_XX_IN, LZ_XY_IN, LZ_YX_IN, LZ_YY_IN
      Real(sz) HZ_X_IN,HZ_Y_IN,TZ_X_IN,TZ_Y_IN
      Real(SZ) MZ_X_IN(s%layers),MZ_Y_IN(s%layers)

      dg%test_el = 0
      DO 1000 L = 1,dg%NREDS
         
!.....Retrieve the global and local edge number

         GED = dg%NREDN(L)
         LED = dg%NEDSD(1,GED)

!.....Retrieve the elements which share the edge

         EL_IN = dg%NEDEL(1,GED)

         dg%pa = PDG_EL(EL_IN)

#ifdef P0
         if (dg%pa.eq.0) then
            dg%pa = 1
         endif
#endif
         
!.....If the element is dry then skip the edge calculation

         IF (dg%WDFLG(EL_IN).EQ.0) GOTO 1000
         dg%test_el = dg%test_el+1

!.....Retrieve the components of the normal vector to the edge
         
         dg%NX = dg%COSNX(GED)
         dg%NY = dg%SINNX(GED)
!.....Set the components for the tangential vector to the edge
         
         TX = -dg%NY
         TY =  dg%NX

!.....Compute dg%ZE, dg%QX, dg%QY, and dg%HB at each Gauss point

         DO I = 1,dg%NEGP(dg%pa)

            dg%ZE_IN = dg%ZE(1,EL_IN,dg%IRK)
            dg%QX_IN = dg%QX(1,EL_IN,dg%IRK)
            dg%QY_IN = dg%QY(1,EL_IN,dg%IRK)
            dg%HB_IN = dg%BATHED(I,LED,EL_IN,dg%pa)

            dg%SFAC_IN = dg%SFACED(I,LED,EL_IN,dg%pa)
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

            dg%iota_IN = dg%iota(1,EL_IN,dg%IRK)
#endif

#ifdef CHEM
            dg%iota_IN = dg%iota(1,EL_IN,dg%IRK)
            dg%iota2_IN = dg%iota2(1,EL_IN,dg%IRK)
#endif

#ifdef DYNP
            dynP_IN = dg%dynP(1,EL_IN,dg%IRK)
#endif

!When layered, these change
#ifdef SED_LAY                  
            dg%HB(1,EL_IN,dg%irk) = 0.D0
            do ll=1,s%layers
               dg%HB(1,EL_IN,dg%irk) = dg%HB(1,EL_IN,dg%irk) + dg%bed(1,EL_IN,dg%irk,ll)

               MZ_X_IN(ll) =  dg%MZ(K,1,ll,EL_IN)
               MZ_Y_IN(ll) =  dg%MZ(K,2,ll,EL_IN)
            enddo
            dg%bed_IN(:) = dg%bed(1,EL_IN,dg%irk,:)
            dg%HB_IN = dg%HB(1,EL_IN,dg%irk)
#endif

!.....Compute the solution at the interior state

            DO K = 2,dg%DOFS(EL_IN)
               
               dg%ZE_IN = dg%ZE_IN + dg%ZE(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
               dg%QX_IN = dg%QX_IN + dg%QX(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
               dg%QY_IN = dg%QY_IN + dg%QY(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)

                                ! LDG terms
#ifdef WAVE_DIF 
               HZ_X_IN = HZ_X_IN + dg%HZ(K,1,1,EL_IN)*dg%PHI_EDGE(K,I,LED,dg%pa)
               HZ_Y_IN = HZ_Y_IN + dg%HZ(K,2,2,EL_IN)*dg%PHI_EDGE(K,I,LED,dg%pa)
#endif

               LZ_XX_IN = LZ_XX_IN + dg%LZ(K,1,1,EL_IN)*dg%PHI_EDGE(K,I,LED,dg%pa)
               LZ_XY_IN = LZ_XY_IN + dg%LZ(K,1,2,EL_IN)*dg%PHI_EDGE(K,I,LED,dg%pa)
               LZ_YX_IN = LZ_YX_IN + dg%LZ(K,2,1,EL_IN)*dg%PHI_EDGE(K,I,LED,dg%pa)
               LZ_YY_IN = LZ_YY_IN + dg%LZ(K,2,2,EL_IN)*dg%PHI_EDGE(K,I,LED,dg%pa)

#ifdef TRACE
               TZ_X_IN = TZ_X_IN + dg%TZ(K,1,1,EL_IN)*dg%PHI_EDGE(K,I,LED,dg%pa)
               TZ_Y_IN = TZ_Y_IN + dg%TZ(K,2,2,EL_IN)*dg%PHI_EDGE(K,I,LED,dg%pa)

               dg%iota_IN = dg%iota_IN + dg%iota(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
#endif

#ifdef CHEM
               dg%iota_IN = dg%iota_IN + dg%iota(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
               dg%iota2_IN = dg%iota2_IN + dg%iota2(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
#endif

#ifdef DYNP
               dynP_IN = dynP_IN + dg%dynP(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
#endif

#ifdef SED_LAY
               do ll = 1,s%layers
                  dg%bed_IN(ll) = dg%bed_IN(ll) + dg%bed(K,EL_IN,dg%IRK,ll)*dg%PHI_EDGE(K,I,LED,dg%pa)
                  dg%HB_IN = dg%HB_IN + dg%bed(k,EL_IN,dg%irk,ll)*dg%PHI_EDGE(K,I,LED,dg%pa)

                  MZ_X_IN(ll) = MZ_X_IN(ll) + dg%MZ(K,1,ll,EL_IN)*dg%PHI_EDGE(K,I,LED,dg%pa)
                  MZ_Y_IN(ll) = MZ_Y_IN(ll) + dg%MZ(K,2,ll,EL_IN)*dg%PHI_EDGE(K,I,LED,dg%pa)
               enddo
#endif

            ENDDO

!.....Compute the velocity in the normal and tangental direction
            
            Q_N_INT = dg%QX_IN*dg%NX + dg%QY_IN*dg%NY
            Q_T_INT = dg%QX_IN*TX + dg%QY_IN*TY

!.....Reflect the velocity in the normal direction

!     Q_N_EXT = -Q_N_INT
            Q_N_EXT = -Q_N_INT+2*SQRT(G/(dg%HB_IN+dg%ZE_IN))*dg%ZE_IN
            Q_T_EXT =  Q_T_INT
            
!.....Compute the x and y components of the external state flow

            DEN = 1.D0/(dg%NX*TY - dg%NY*TX)
            dg%QX_EX = ( TY*Q_N_EXT - dg%NY*Q_T_EXT)*DEN
            dg%QY_EX = (-TX*Q_N_EXT + dg%NX*Q_T_EXT)*DEN

            dg%ZE_EX = dg%ZE_IN
            dg%HB_EX = dg%HB_IN

            dg%SFAC_EX = dg%SFAC_IN

#ifdef SED_LAY
               do ll = 1,s%layers
                  dg%bed_EX(ll) = dg%bed_IN(ll)
               enddo
#endif

#ifdef TRACE
            dg%iota_EX = dg%iota_IN
#endif

#ifdef CHEM
            dg%iota_EX = dg%iota_IN
            dg%iota2_EX = dg%iota2_IN
#endif

#ifdef DYNP
            dynP_EX = dynP_IN
#endif

!.....Compute the flux

            CALL NUMERICAL_FLUX(s,dg,IT,dg%test_el)
            
!.....Compute the edge integral

            
            DO K = 1,dg%DOFS(EL_IN)
               W_IN = 2.0*dg%M_INV(K,dg%pa)/AREAS(EL_IN)*dg%XLEN(GED)*&
                   dg%PHI_EDGE(K,I,LED,dg%pa)*dg%WEGP(I,dg%pa)
               dg%RHS_ZE(K,EL_IN,dg%IRK) = dg%RHS_ZE(K,EL_IN,dg%IRK) - W_IN*F_HAT
               dg%RHS_QX(K,EL_IN,dg%IRK) = dg%RHS_QX(K,EL_IN,dg%IRK) - W_IN*G_HAT
               dg%RHS_QY(K,EL_IN,dg%IRK) = dg%RHS_QY(K,EL_IN,dg%IRK) - W_IN*H_HAT

#ifdef TRACE
               dg%RHS_iota(K,EL_IN,dg%IRK) = dg%RHS_iota(K,EL_IN,dg%IRK) - W_IN*I_HAT
#endif

#ifdef CHEM
               dg%RHS_iota(K,EL_IN,dg%IRK) = dg%RHS_iota(K,EL_IN,dg%IRK) - W_IN*I_HAT
               dg%RHS_iota2(K,EL_IN,dg%IRK) = dg%RHS_iota2(K,EL_IN,dg%IRK) - W_IN*J_HAT
#endif

#ifdef DYNP
               dg%RHS_dynP(K,EL_IN,dg%IRK) = dg%RHS_dynP(K,EL_IN,dg%IRK) - W_IN*K_HAT
#endif

#ifdef SED_LAY
               do ll = 1,s%layers
                  dg%RHS_bed(K,EL_IN,dg%IRK,ll) = dg%RHS_bed(K,EL_IN,dg%IRK,ll) - W_IN*dg%bed_HAT(ll)
               enddo
#endif

               
!     !CALL EDGE_INT_HYDRO(EL_IN, LED, GED, I, F_HAT, G_HAT, H_HAT,k)

            ENDDO

         ENDDO

 1000 CONTINUE      
      RETURN
      END SUBROUTINE
