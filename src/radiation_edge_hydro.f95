!***********************************************************************
!     
!     SUBROUTINE RADIATION_EDGE_HYDRO( )
!     
!     This subroutine does the following:
!     
!     1.  Calculates the values of the necessary variables at the global_here%edge
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

      SUBROUTINE RADIATION_EDGE_HYDRO(s,dg_here,global_here,IT)

!.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE sizes
      IMPLICIT NONE
      
      type (sizes_type) :: s
      type (dg_type) :: dg_here
      type (global_type) :: global_here

!.....Declare local variables

      INTEGER L, LED, GED,i,k,ll,IT
      REAL(SZ) TX, TY, W_IN, DEN
      REAL(SZ) LZ_XX_IN, LZ_XY_IN, LZ_YX_IN, LZ_YY_IN
      Real(sz) HZ_X_IN,HZ_Y_IN,TZ_X_IN,TZ_Y_IN
      Real(SZ) MZ_X_IN(s%layers),MZ_Y_IN(s%layers)

      dg_here%test_el = 0
      DO 1000 L = 1,dg_here%NREDS
         
!.....Retrieve the global and local global_here%edge number

         GED = dg_here%NREDN(L)
         LED = dg_here%NEDSD(1,GED)

!.....Retrieve the elements which share the global_here%edge

         global_here%EL_IN = dg_here%NEDEL(1,GED)

         dg_here%pa = global_here%PDG_EL(global_here%EL_IN)

#ifdef P0
         if (dg_here%pa.eq.0) then
            dg_here%pa = 1
         endif
#endif
         
!.....If the element is dry then skip the global_here%edge calculation

         IF (dg_here%WDFLG(global_here%EL_IN).EQ.0) GOTO 1000
         dg_here%test_el = dg_here%test_el+1

!.....Retrieve the components of the normal vector to the global_here%edge
         
         dg_here%NX = dg_here%COSNX(GED)
         dg_here%NY = dg_here%SINNX(GED)
!.....Set the components for the tangential vector to the global_here%edge
         
         TX = -dg_here%NY
         TY =  dg_here%NX

!.....Compute dg_here%ZE, dg_here%QX, dg_here%QY, and dg_here%HB at each Gauss point

         DO I = 1,dg_here%NEGP(dg_here%pa)

            dg_here%ZE_IN = dg_here%ZE(1,global_here%EL_IN,dg_here%IRK)
            dg_here%QX_IN = dg_here%QX(1,global_here%EL_IN,dg_here%IRK)
            dg_here%QY_IN = dg_here%QY(1,global_here%EL_IN,dg_here%IRK)
            dg_here%HB_IN = dg_here%BATHED(I,LED,global_here%EL_IN,dg_here%pa)

            dg_here%SFAC_IN = dg_here%SFACED(I,LED,global_here%EL_IN,dg_here%pa)
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

            dg_here%iota_IN = dg_here%iota(1,global_here%EL_IN,dg_here%IRK)
#endif

#ifdef CHEM
            dg_here%iota_IN = dg_here%iota(1,global_here%EL_IN,dg_here%IRK)
            dg_here%iota2_IN = dg_here%iota2(1,global_here%EL_IN,dg_here%IRK)
#endif

#ifdef DYNP
            dynP_IN = dg_here%dynP(1,global_here%EL_IN,dg_here%IRK)
#endif

!When layered, these change
#ifdef SED_LAY                  
            dg_here%HB(1,global_here%EL_IN,dg_here%irk) = 0.D0
            do ll=1,s%layers
               dg_here%HB(1,global_here%EL_IN,dg_here%irk) = dg_here%HB(1,global_here%EL_IN,dg_here%irk) + dg_here%bed(1,global_here%EL_IN,dg_here%irk,ll)

               MZ_X_IN(ll) =  dg_here%MZ(K,1,ll,global_here%EL_IN)
               MZ_Y_IN(ll) =  dg_here%MZ(K,2,ll,global_here%EL_IN)
            enddo
            dg_here%bed_IN(:) = dg_here%bed(1,global_here%EL_IN,dg_here%irk,:)
            dg_here%HB_IN = dg_here%HB(1,global_here%EL_IN,dg_here%irk)
#endif

!.....Compute the solution at the interior state

            DO K = 2,dg_here%DOFS(global_here%EL_IN)
               
               dg_here%ZE_IN = dg_here%ZE_IN + dg_here%ZE(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               dg_here%QX_IN = dg_here%QX_IN + dg_here%QX(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               dg_here%QY_IN = dg_here%QY_IN + dg_here%QY(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)

                                ! LDG terms
#ifdef WAVE_DIF 
               HZ_X_IN = HZ_X_IN + dg_here%HZ(K,1,1,global_here%EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               HZ_Y_IN = HZ_Y_IN + dg_here%HZ(K,2,2,global_here%EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
#endif

               LZ_XX_IN = LZ_XX_IN + dg_here%LZ(K,1,1,global_here%EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               LZ_XY_IN = LZ_XY_IN + dg_here%LZ(K,1,2,global_here%EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               LZ_YX_IN = LZ_YX_IN + dg_here%LZ(K,2,1,global_here%EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               LZ_YY_IN = LZ_YY_IN + dg_here%LZ(K,2,2,global_here%EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)

#ifdef TRACE
               TZ_X_IN = TZ_X_IN + dg_here%TZ(K,1,1,global_here%EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               TZ_Y_IN = TZ_Y_IN + dg_here%TZ(K,2,2,global_here%EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)

               dg_here%iota_IN = dg_here%iota_IN + dg_here%iota(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
#endif

#ifdef CHEM
               dg_here%iota_IN = dg_here%iota_IN + dg_here%iota(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               dg_here%iota2_IN = dg_here%iota2_IN + dg_here%iota2(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
#endif

#ifdef DYNP
               dynP_IN = dynP_IN + dg_here%dynP(K,global_here%EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
#endif

#ifdef SED_LAY
               do ll = 1,s%layers
                  dg_here%bed_IN(ll) = dg_here%bed_IN(ll) + dg_here%bed(K,global_here%EL_IN,dg_here%IRK,ll)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
                  dg_here%HB_IN = dg_here%HB_IN + dg_here%bed(k,global_here%EL_IN,dg_here%irk,ll)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)

                  MZ_X_IN(ll) = MZ_X_IN(ll) + dg_here%MZ(K,1,ll,global_here%EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
                  MZ_Y_IN(ll) = MZ_Y_IN(ll) + dg_here%MZ(K,2,ll,global_here%EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               enddo
#endif

            ENDDO

!.....Compute the velocity in the normal and tangental direction
            
            global_here%Q_N_INT = dg_here%QX_IN*dg_here%NX + dg_here%QY_IN*dg_here%NY
            global_here%Q_T_INT = dg_here%QX_IN*TX + dg_here%QY_IN*TY

!.....Reflect the velocity in the normal direction

!     global_here%Q_N_EXT = -global_here%Q_N_INT
            global_here%Q_N_EXT = -global_here%Q_N_INT+2*SQRT(global_here%G/(dg_here%HB_IN+dg_here%ZE_IN))*dg_here%ZE_IN
            global_here%Q_T_EXT =  global_here%Q_T_INT
            
!.....Compute the global_here%x and global_here%y components of the external state flow

            DEN = 1.D0/(dg_here%NX*TY - dg_here%NY*TX)
            dg_here%QX_EX = ( TY*global_here%Q_N_EXT - dg_here%NY*global_here%Q_T_EXT)*DEN
            dg_here%QY_EX = (-TX*global_here%Q_N_EXT + dg_here%NX*global_here%Q_T_EXT)*DEN

            dg_here%ZE_EX = dg_here%ZE_IN
            dg_here%HB_EX = dg_here%HB_IN

            dg_here%SFAC_EX = dg_here%SFAC_IN

#ifdef SED_LAY
               do ll = 1,s%layers
                  dg_here%bed_EX(ll) = dg_here%bed_IN(ll)
               enddo
#endif

#ifdef TRACE
            dg_here%iota_EX = dg_here%iota_IN
#endif

#ifdef CHEM
            dg_here%iota_EX = dg_here%iota_IN
            dg_here%iota2_EX = dg_here%iota2_IN
#endif

#ifdef DYNP
            dynP_EX = dynP_IN
#endif

!.....Compute the flux

            CALL NUMERICAL_FLUX(s,dg_here,global_here,IT,dg_here%test_el)
            
!.....Compute the global_here%edge integral

            
            DO K = 1,dg_here%DOFS(global_here%EL_IN)
               W_IN = 2.0*dg_here%M_INV(K,dg_here%pa)/global_here%AREAS(global_here%EL_IN)*dg_here%XLEN(GED)*&
                   dg_here%PHI_EDGE(K,I,LED,dg_here%pa)*dg_here%WEGP(I,dg_here%pa)
               dg_here%RHS_ZE(K,global_here%EL_IN,dg_here%IRK) = dg_here%RHS_ZE(K,global_here%EL_IN,dg_here%IRK) - W_IN*global_here%F_HAT
               dg_here%RHS_QX(K,global_here%EL_IN,dg_here%IRK) = dg_here%RHS_QX(K,global_here%EL_IN,dg_here%IRK) - W_IN*global_here%G_HAT
               dg_here%RHS_QY(K,global_here%EL_IN,dg_here%IRK) = dg_here%RHS_QY(K,global_here%EL_IN,dg_here%IRK) - W_IN*global_here%H_HAT

#ifdef TRACE
               dg_here%RHS_iota(K,global_here%EL_IN,dg_here%IRK) = dg_here%RHS_iota(K,global_here%EL_IN,dg_here%IRK) - W_IN*global_here%I_HAT
#endif

#ifdef CHEM
               dg_here%RHS_iota(K,global_here%EL_IN,dg_here%IRK) = dg_here%RHS_iota(K,global_here%EL_IN,dg_here%IRK) - W_IN*global_here%I_HAT
               dg_here%RHS_iota2(K,global_here%EL_IN,dg_here%IRK) = dg_here%RHS_iota2(K,global_here%EL_IN,dg_here%IRK) - W_IN*global_here%J_HAT
#endif

#ifdef DYNP
               dg_here%RHS_dynP(K,global_here%EL_IN,dg_here%IRK) = dg_here%RHS_dynP(K,global_here%EL_IN,dg_here%IRK) - W_IN*global_here%K_HAT
#endif

#ifdef SED_LAY
               do ll = 1,s%layers
                  dg_here%RHS_bed(K,global_here%EL_IN,dg_here%IRK,ll) = dg_here%RHS_bed(K,global_here%EL_IN,dg_here%IRK,ll) - W_IN*dg_here%bed_HAT(ll)
               enddo
#endif

               
!     !CALL EDGE_INT_HYDRO(global_here%EL_IN, LED, GED, I, global_here%F_HAT, global_here%G_HAT, global_here%H_HAT,k)

            ENDDO

         ENDDO

 1000 CONTINUE      
      RETURN
      END SUBROUTINE
