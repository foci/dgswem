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
      SUBROUTINE OCEAN_EDGE_HYDRO(s,IT)

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

!.....Declare local variables
      
      INTEGER L, LED, GED, i,k,jj,II,ll,IT,w
      Real(SZ) DEN2,U_AVG,V_AVG,VEL_NORMAL,q_RoeX, q_RoeY, q_Roe
      REAL(SZ) TX, TY, HUU, HUV, GH2,FH_NL_IN,F1_NL,FX1_IN,FY1_IN
      REAL(SZ) LZ_XX_IN, LZ_XY_IN, LZ_YX_IN, LZ_YY_IN,W_IN,FX3_IN
      Real(SZ) chi_pref,MZ_X_IN(s%layers),MZ_Y_IN(s%layers)
      Real(SZ) HZ_X_IN,HZ_Y_IN,TZ_X_IN,TZ_Y_IN

      dg%test_el = 0
      DO 1000 L = 1, dg%needs
         
!.....Retrieve the global and local edge number

         GED = dg%NEEDN(L)
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
         
!.....Retrieve the nodes of the edge
         
         N1 = dg%NEDNO(1,GED)
         N2 = dg%NEDNO(2,GED)
         
!.....Compute dg%ZE, dg%QX, dg%QY, and dg%HB at each edge Gauss quadrature point

         DO I = 1,dg%NEGP(dg%pa)

            dg%ZE_IN = dg%ZE(1,EL_IN,dg%IRK)
            dg%QX_IN = dg%QX(1,EL_IN,dg%IRK)
            dg%QY_IN = dg%QY(1,EL_IN,dg%IRK)

            dg%ZE_EX = 0.D0
            dg%QX_EX = 0.D0
            dg%QY_EX = 0.D0
            U_EX  = 0.D0
            V_EX  = 0.D0

#ifdef TRACE
            dg%iota_IN = dg%iota(1,EL_IN,dg%IRK)
            dg%iota_EX = 0.D0
#endif

#ifdef CHEM
            dg%iota_IN = dg%iota(1,EL_IN,dg%IRK)
            dg%iota2_IN = dg%iota2(1,EL_IN,dg%IRK)
            dg%iota_EX = 0.D0
            dg%iota2_EX = 0.D0
#endif

#ifdef dynp
            dynP_IN = dg%dynP(1,EL_IN,dg%IRK)
            dynP_EX = 0.D0
#endif
            
            dg%HB_IN = dg%BATHED(I,LED,EL_IN,dg%pa)
            dg%SFAC_IN = dg%SFACED(I,LED,EL_IN,dg%pa)
            
            !When layered, these change
#ifdef SED_LAY 
            dg%HB(:,EL_IN,dg%irk) = 0.D0
            do ll=1,s%layers
               dg%HB(1,EL_IN,dg%irk) = dg%HB(1,EL_IN,dg%irk) + dg%bed(1,EL_IN,dg%irk,ll)

               MZ_X_IN(ll) =  dg%MZ(1,1,ll,EL_IN)
               MZ_Y_IN(ll) =  dg%MZ(1,2,ll,EL_IN)
            enddo
            dg%bed_IN(:) = dg%bed(1,EL_IN,dg%irk,:)
            dg%HB_IN = dg%HB(1,EL_IN,dg%irk)
#endif
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

!.....Compute the specified open ocean elevation
            
            DO JJ=1,NBFR
               
               IF (PER(JJ).EQ.0.D0) THEN
                  NCYC = 0.D0
               ELSE
                  NCYC = INT(dg%TIMEDG/PER(JJ))
               ENDIF
               
!...........Surface Elevation

               ARGJ = AMIG(JJ)*(dg%TIMEDG - NCYC*PER(JJ)) + FACE(JJ)
               RFF = FF(JJ)*dg%RAMPDG
               
               dg%EFA_GP = 0.5D0*(dg%EFA_DG(JJ,L,1) + dg%EFA_DG(JJ,L,2))&
                + 0.5D0*(dg%EFA_DG(JJ,L,2) - dg%EFA_DG(JJ,L,1))*dg%XEGP(I,dg%pa)
               dg%EMO_GP = 0.5D0*(dg%EMO_DG(JJ,L,1) + dg%EMO_DG(JJ,L,2))&
                + 0.5D0*(dg%EMO_DG(JJ,L,2) - dg%EMO_DG(JJ,L,1))*dg%XEGP(I,dg%pa)

               ARG = ARGJ - dg%EFA_GP
               
               dg%ZE_EX = dg%ZE_EX + dg%EMO_GP*RFF*COS(ARG) 
      
            ENDDO


!.....Compute the solution at the interior state

            DO K = 2,dg%DOFS(EL_IN)

               dg%ZE_IN = dg%ZE_IN + dg%ZE(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
               dg%QX_IN = dg%QX_IN + dg%QX(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
               dg%QY_IN = dg%QY_IN + dg%QY(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)

#ifdef TRACE
               dg%iota_IN = dg%iota_IN + dg%iota(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
#endif

#ifdef CHEM
               dg%iota_IN = dg%iota_IN + dg%iota(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
               dg%iota2_IN = dg%iota2_IN + dg%iota2(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
#endif

#ifdef dynp
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
#endif

            ENDDO

!.....Set the exterior value of the bathymetry equal to the interior

            dg%HB_EX = dg%HB_IN

#ifdef SED_LAY
            dg%bed_EX(:) = dg%bed_IN(:)
#endif
            dg%SFAC_EX = dg%SFAC_IN

            IF (LoadGeoidOffset) then
               dg%ZE_EX = dg%ZE_EX + .5*(GeoidOffset(N1)+GeoidOffset(N2))
            endif


            !dg%ZE_EX = dg%ZE_IN !this was added

#ifdef TRACE
            dg%iota_EX = 0.D0 !dg%iota_IN !might need 0.D0
#endif

#ifdef CHEM
            dg%iota_EX = 0.D0 !dg%iota_IN
            dg%iota2_EX = 0.D0 !dg%iota2_IN
#endif

#ifdef dynp
            dynP_EX = 0.D0 !What's the right setting here?
#endif

            !When layered, these change
#ifdef SED_LAY
            dg%bed_EX(:) = 10.D0 !dg%bed_IN(:) 
#endif


!.....Set the exterior state flows equal to the interior state flows

            dg%QX_EX = dg%QX_IN
            dg%QY_EX = dg%QY_IN

!.....Compute the Roe flux

            CALL NUMERICAL_FLUX(s,IT,dg%test_el)

!.....Add LDG terms
            
#ifdef WAVE_DIF 
            F_HAT = F_HAT + HZ_X_IN*dg%NX*dg%SFAC_IN + HZ_Y_IN*dg%NY
#endif
            G_HAT = G_HAT + LZ_XX_IN*dg%NX*dg%SFAC_IN + LZ_XY_IN*dg%NY
            H_HAT = H_HAT + LZ_YX_IN*dg%NX*dg%SFAC_IN + LZ_YY_IN*dg%NY
#ifdef TRACE
            I_HAT = I_HAT + TZ_X_IN*dg%NX*dg%SFAC_IN + TZ_Y_IN*dg%NY
#endif

!.....Add LDG terms for sediment

#ifdef SED_LAY
            do ll=1,s%layers
               dg%bed_HAT(ll) = dg%bed_HAT(ll) + MZ_X_IN(ll)*dg%NX*dg%SFAC_IN + MZ_Y_IN(ll)*dg%NY
            enddo
#endif

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

#ifdef dynp
               dg%RHS_dynP(K,EL_IN,dg%IRK) = dg%RHS_dynP(K,EL_IN,dg%IRK) - W_IN*K_HAT
#endif

#ifdef SED_LAY
               
               do ll = 1,s%layers
                  dg%RHS_bed(K,EL_IN,dg%IRK,ll) = dg%RHS_bed(K,EL_IN,dg%IRK,ll) - W_IN*dg%bed_HAT(ll)
               enddo
#endif

            ENDDO
         ENDDO

 1000 CONTINUE
      
      RETURN
      END SUBROUTINE
