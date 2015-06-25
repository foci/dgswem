!***********************************************************************
!     
!     SUBROUTINE OCEAN_EDGE_HYDRO( )
!     
!     This subroutine does the following:
!     
!     1.  Calculates the values of the necessary variables at the edge
!     gauss points for ELEVATION SPECIFIED edges
!     2.  Calls the appropriate subroutine to compute the flux at
!     these points.
!     3.  Calls the appropriate subroutine to compute the boundary
!     integrals.
!     
!     Written by Ethan Kubatko (06-11-2004)
!     
!     01-10-2011 - cem - adapted for p_enrichment and multicomponent  
!     
!***********************************************************************
      SUBROUTINE OCEAN_EDGE_HYDRO(s,dg_here,IT)

!.....Use appropriate modules

      USE GLOBAL
      USE DG
      Use SIZES
!     
      USE NodalAttributes, ONLY: GeoidOffset, LoadGeoidOffset
      use fparser
      use fparser2

      IMPLICIT NONE

      type (sizes_type) :: s
      type (dg_type) :: dg_here

!.....Declare local variables
      
      INTEGER L, LED, GED, i,k,jj,II,ll,IT,w
      Real(SZ) DEN2,U_AVG,V_AVG,VEL_NORMAL,q_RoeX, q_RoeY, q_Roe
      REAL(SZ) TX, TY, HUU, HUV, GH2,FH_NL_IN,F1_NL,FX1_IN,FY1_IN
      REAL(SZ) LZ_XX_IN, LZ_XY_IN, LZ_YX_IN, LZ_YY_IN,W_IN,FX3_IN
      Real(SZ) chi_pref,MZ_X_IN(s%layers),MZ_Y_IN(s%layers)
      Real(SZ) HZ_X_IN,HZ_Y_IN,TZ_X_IN,TZ_Y_IN

      dg_here%test_el = 0
      DO 1000 L = 1, dg_here%needs
         
!.....Retrieve the global and local edge number

         GED = dg_here%NEEDN(L)
         LED = dg_here%NEDSD(1,GED)

!.....Retrieve the elements which share the edge

         EL_IN = dg_here%NEDEL(1,GED)
         
         dg_here%pa = PDG_EL(EL_IN)

#ifdef P0
         if (dg_here%pa.eq.0) then
            dg_here%pa = 1
         endif
#endif
         
!.....If the element is dry then skip the edge calculation

         IF (dg_here%WDFLG(EL_IN).EQ.0) GOTO 1000
         
         dg_here%test_el = dg_here%test_el+1
         
!.....Retrieve the components of the normal vector to the edge
         
         dg_here%NX = dg_here%COSNX(GED)
         dg_here%NY = dg_here%SINNX(GED)
         
!.....Retrieve the nodes of the edge
         
         N1 = dg_here%NEDNO(1,GED)
         N2 = dg_here%NEDNO(2,GED)
         
!.....Compute dg_here%ZE, dg_here%QX, dg_here%QY, and dg_here%HB at each edge Gauss quadrature point

         DO I = 1,dg_here%NEGP(dg_here%pa)

            dg_here%ZE_IN = dg_here%ZE(1,EL_IN,dg_here%IRK)
            dg_here%QX_IN = dg_here%QX(1,EL_IN,dg_here%IRK)
            dg_here%QY_IN = dg_here%QY(1,EL_IN,dg_here%IRK)

            dg_here%ZE_EX = 0.D0
            dg_here%QX_EX = 0.D0
            dg_here%QY_EX = 0.D0
            U_EX  = 0.D0
            V_EX  = 0.D0

#ifdef TRACE
            dg_here%iota_IN = dg_here%iota(1,EL_IN,dg_here%IRK)
            dg_here%iota_EX = 0.D0
#endif

#ifdef CHEM
            dg_here%iota_IN = dg_here%iota(1,EL_IN,dg_here%IRK)
            dg_here%iota2_IN = dg_here%iota2(1,EL_IN,dg_here%IRK)
            dg_here%iota_EX = 0.D0
            dg_here%iota2_EX = 0.D0
#endif

#ifdef dynp
            dynP_IN = dg_here%dynP(1,EL_IN,dg_here%IRK)
            dynP_EX = 0.D0
#endif
            
            dg_here%HB_IN = dg_here%BATHED(I,LED,EL_IN,dg_here%pa)
            dg_here%SFAC_IN = dg_here%SFACED(I,LED,EL_IN,dg_here%pa)
            
            !When layered, these change
#ifdef SED_LAY 
            dg_here%HB(:,EL_IN,dg_here%irk) = 0.D0
            do ll=1,s%layers
               dg_here%HB(1,EL_IN,dg_here%irk) = dg_here%HB(1,EL_IN,dg_here%irk) + dg_here%bed(1,EL_IN,dg_here%irk,ll)

               MZ_X_IN(ll) =  dg_here%MZ(1,1,ll,EL_IN)
               MZ_Y_IN(ll) =  dg_here%MZ(1,2,ll,EL_IN)
            enddo
            dg_here%bed_IN(:) = dg_here%bed(1,EL_IN,dg_here%irk,:)
            dg_here%HB_IN = dg_here%HB(1,EL_IN,dg_here%irk)
#endif
#ifdef WAVE_DIF 
            HZ_X_IN = dg_here%HZ(1,1,1,EL_IN)
            HZ_Y_IN = dg_here%HZ(1,2,2,EL_IN)
#endif

            LZ_XX_IN = dg_here%LZ(1,1,1,EL_IN)
            LZ_XY_IN = dg_here%LZ(1,1,2,EL_IN)
            LZ_YX_IN = dg_here%LZ(1,2,1,EL_IN)
            LZ_YY_IN = dg_here%LZ(1,2,2,EL_IN)

#ifdef TRACE
            TZ_X_IN = dg_here%TZ(1,1,1,EL_IN)
            TZ_Y_IN = dg_here%TZ(1,2,2,EL_IN)
#endif

!.....Compute the specified open ocean elevation
            
            DO JJ=1,NBFR
               
               IF (PER(JJ).EQ.0.D0) THEN
                  NCYC = 0.D0
               ELSE
                  NCYC = INT(dg_here%TIMEDG/PER(JJ))
               ENDIF
               
!...........Surface Elevation

               ARGJ = AMIG(JJ)*(dg_here%TIMEDG - NCYC*PER(JJ)) + FACE(JJ)
               RFF = FF(JJ)*dg_here%RAMPDG
               
               dg_here%EFA_GP = 0.5D0*(dg_here%EFA_DG(JJ,L,1) + dg_here%EFA_DG(JJ,L,2))&
                + 0.5D0*(dg_here%EFA_DG(JJ,L,2) - dg_here%EFA_DG(JJ,L,1))*dg_here%XEGP(I,dg_here%pa)
               dg_here%EMO_GP = 0.5D0*(dg_here%EMO_DG(JJ,L,1) + dg_here%EMO_DG(JJ,L,2))&
                + 0.5D0*(dg_here%EMO_DG(JJ,L,2) - dg_here%EMO_DG(JJ,L,1))*dg_here%XEGP(I,dg_here%pa)

               ARG = ARGJ - dg_here%EFA_GP
               
               dg_here%ZE_EX = dg_here%ZE_EX + dg_here%EMO_GP*RFF*COS(ARG) 
      
            ENDDO


!.....Compute the solution at the interior state

            DO K = 2,dg_here%DOFS(EL_IN)

               dg_here%ZE_IN = dg_here%ZE_IN + dg_here%ZE(K,EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               dg_here%QX_IN = dg_here%QX_IN + dg_here%QX(K,EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               dg_here%QY_IN = dg_here%QY_IN + dg_here%QY(K,EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)

#ifdef TRACE
               dg_here%iota_IN = dg_here%iota_IN + dg_here%iota(K,EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
#endif

#ifdef CHEM
               dg_here%iota_IN = dg_here%iota_IN + dg_here%iota(K,EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               dg_here%iota2_IN = dg_here%iota2_IN + dg_here%iota2(K,EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
#endif

#ifdef dynp
               dynP_IN = dynP_IN + dg_here%dynP(K,EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
#endif

#ifdef SED_LAY
               do ll = 1,s%layers
                  dg_here%bed_IN(ll) = dg_here%bed_IN(ll) + dg_here%bed(K,EL_IN,dg_here%IRK,ll)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
                  dg_here%HB_IN = dg_here%HB_IN + dg_here%bed(k,EL_IN,dg_here%irk,ll)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)

                  MZ_X_IN(ll) = MZ_X_IN(ll) + dg_here%MZ(K,1,ll,EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
                  MZ_Y_IN(ll) = MZ_Y_IN(ll) + dg_here%MZ(K,2,ll,EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               enddo
#endif

                                ! LDG terms
#ifdef WAVE_DIF 
               HZ_X_IN = HZ_X_IN + dg_here%HZ(K,1,1,EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               HZ_Y_IN = HZ_Y_IN + dg_here%HZ(K,2,2,EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
#endif

               LZ_XX_IN = LZ_XX_IN + dg_here%LZ(K,1,1,EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               LZ_XY_IN = LZ_XY_IN + dg_here%LZ(K,1,2,EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               LZ_YX_IN = LZ_YX_IN + dg_here%LZ(K,2,1,EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               LZ_YY_IN = LZ_YY_IN + dg_here%LZ(K,2,2,EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)

#ifdef TRACE
               TZ_X_IN = TZ_X_IN + dg_here%TZ(K,1,1,EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               TZ_Y_IN = TZ_Y_IN + dg_here%TZ(K,2,2,EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
#endif

            ENDDO

!.....Set the exterior value of the bathymetry equal to the interior

            dg_here%HB_EX = dg_here%HB_IN

#ifdef SED_LAY
            dg_here%bed_EX(:) = dg_here%bed_IN(:)
#endif
            dg_here%SFAC_EX = dg_here%SFAC_IN

            IF (LoadGeoidOffset) then
               dg_here%ZE_EX = dg_here%ZE_EX + .5*(GeoidOffset(N1)+GeoidOffset(N2))
            endif


            !dg_here%ZE_EX = dg_here%ZE_IN !this was added

#ifdef TRACE
            dg_here%iota_EX = 0.D0 !dg_here%iota_IN !might need 0.D0
#endif

#ifdef CHEM
            dg_here%iota_EX = 0.D0 !dg_here%iota_IN
            dg_here%iota2_EX = 0.D0 !dg_here%iota2_IN
#endif

#ifdef dynp
            dynP_EX = 0.D0 !What's the right setting here?
#endif

            !When layered, these change
#ifdef SED_LAY
            dg_here%bed_EX(:) = 10.D0 !dg_here%bed_IN(:) 
#endif


!.....Set the exterior state flows equal to the interior state flows

            dg_here%QX_EX = dg_here%QX_IN
            dg_here%QY_EX = dg_here%QY_IN

!.....Compute the Roe flux

            CALL NUMERICAL_FLUX(s,dg_here,IT,dg_here%test_el)

!.....Add LDG terms
            
#ifdef WAVE_DIF 
            F_HAT = F_HAT + HZ_X_IN*dg_here%NX*dg_here%SFAC_IN + HZ_Y_IN*dg_here%NY
#endif
            G_HAT = G_HAT + LZ_XX_IN*dg_here%NX*dg_here%SFAC_IN + LZ_XY_IN*dg_here%NY
            H_HAT = H_HAT + LZ_YX_IN*dg_here%NX*dg_here%SFAC_IN + LZ_YY_IN*dg_here%NY
#ifdef TRACE
            I_HAT = I_HAT + TZ_X_IN*dg_here%NX*dg_here%SFAC_IN + TZ_Y_IN*dg_here%NY
#endif

!.....Add LDG terms for sediment

#ifdef SED_LAY
            do ll=1,s%layers
               dg_here%bed_HAT(ll) = dg_here%bed_HAT(ll) + MZ_X_IN(ll)*dg_here%NX*dg_here%SFAC_IN + MZ_Y_IN(ll)*dg_here%NY
            enddo
#endif

!.....Compute the edge integral
            DO K = 1,dg_here%DOFS(EL_IN)

               W_IN = 2.0*dg_here%M_INV(K,dg_here%pa)/AREAS(EL_IN)*dg_here%XLEN(GED)*&
              dg_here%PHI_EDGE(K,I,LED,dg_here%pa)*dg_here%WEGP(I,dg_here%pa)

               dg_here%RHS_ZE(K,EL_IN,dg_here%IRK) = dg_here%RHS_ZE(K,EL_IN,dg_here%IRK) - W_IN*F_HAT
               dg_here%RHS_QX(K,EL_IN,dg_here%IRK) = dg_here%RHS_QX(K,EL_IN,dg_here%IRK) - W_IN*G_HAT
               dg_here%RHS_QY(K,EL_IN,dg_here%IRK) = dg_here%RHS_QY(K,EL_IN,dg_here%IRK) - W_IN*H_HAT

#ifdef TRACE
               dg_here%RHS_iota(K,EL_IN,dg_here%IRK) = dg_here%RHS_iota(K,EL_IN,dg_here%IRK) - W_IN*I_HAT
#endif

#ifdef CHEM
               dg_here%RHS_iota(K,EL_IN,dg_here%IRK) = dg_here%RHS_iota(K,EL_IN,dg_here%IRK) - W_IN*I_HAT
               dg_here%RHS_iota2(K,EL_IN,dg_here%IRK) = dg_here%RHS_iota2(K,EL_IN,dg_here%IRK) - W_IN*J_HAT
#endif

#ifdef dynp
               dg_here%RHS_dynP(K,EL_IN,dg_here%IRK) = dg_here%RHS_dynP(K,EL_IN,dg_here%IRK) - W_IN*K_HAT
#endif

#ifdef SED_LAY
               
               do ll = 1,s%layers
                  dg_here%RHS_bed(K,EL_IN,dg_here%IRK,ll) = dg_here%RHS_bed(K,EL_IN,dg_here%IRK,ll) - W_IN*dg_here%bed_HAT(ll)
               enddo
#endif

            ENDDO
         ENDDO

 1000 CONTINUE
      
      RETURN
      END SUBROUTINE
