!***********************************************************************
!     
!     SUBROUTINE FLOW_EDGE_HYDRO( )
!     
!     This subroutine does the following:
!     
!     1.  Calculates the values of the necessary variables at the edge
!     gauss points for NON-ZERO FLUX edges
!     2.  Calls the appropriate subroutine to compute the flux at
!     these points.
!     3.  Calls the appropriate subroutine to compute the boundary
!     integrals.
!     
!     Written by Ethan Kubatko (06-11-2004)

!     01-10-2011 - cem - adapted for p_enrichment and multicomponent
!     
!***********************************************************************

      SUBROUTINE FLOW_EDGE_HYDRO(s,dg_here,IT)

!.....Use appropriate modules

      USE GLOBAL
      USE DG
      use SIZES
      use fparser
      use fparser2

      IMPLICIT NONE

      type (sizes_type) :: s
      type (dg_type) :: dg_here

!.....Declare local variables

      INTEGER L, LED, GED, i,k,jj,II,ll,IT,w
      Real(SZ) DEN2,U_AVG,V_AVG,VEL_NORMAL,q_RoeX, q_RoeY, q_Roe
      REAL(SZ) TX, TY,HZ_X_IN,HZ_Y_IN,TZ_X_IN,TZ_Y_IN
      REAL(SZ) LZ_XX_IN, LZ_XY_IN, LZ_YX_IN, LZ_YY_IN,W_IN
      Real(SZ) VEI(3,3),chi_pref,C_char,MZ_X_IN(s%layers),MZ_Y_IN(s%layers)

      dg_here%test_el = 0
      DO 1000 L = 1,dg_here%NFEDS
         
!.....Retrieve the global and local edge number

         GED = dg_here%NFEDN(L)
         LED = dg_here%NEDSD(1,GED)

!.....Retrieve the element to which the edge belongs

         EL_IN = dg_here%NEDEL(1,GED)

!.....Organize the dg_here%dofs for the p_adaptivity

         dg_here%PA = PDG_EL(EL_IN)

#ifdef P0
         if (dg_here%pa.eq.0) then
            dg_here%pa = 1
         endif
#endif

!.....If the element is dry then skip the edge calculation

         IF (dg_here%WDFLG(EL_IN).EQ.0) GOTO 1000
         
!.....Retrieve the components of the normal vector to the edge
         
         dg_here%NX = dg_here%COSNX(GED)
         dg_here%NY = dg_here%SINNX(GED)
         
!.....Set the components for the tangential vector to the edge

         TX = -dg_here%NY
         TY =  dg_here%NX
         
!.....Compute dg_here%ZE, dg_here%QX, dg_here%QY, and dg_here%HB at each edge Gauss quadrature point

         DO I = 1,dg_here%NEGP(dg_here%pa)

            dg_here%ZE_IN = 0.D0
            dg_here%QX_IN = 0.D0
            dg_here%QY_IN = 0.D0

            dg_here%ZE_EX = 0.D0
            dg_here%QX_EX = 0.D0
            dg_here%QY_EX = 0.D0
            
            dg_here%HB_IN = 0.D0
            dg_here%SFAC_IN = dg_here%SFACED(I,LED,EL_IN,dg_here%pa)

#ifdef TRACE
            dg_here%iota_IN = 0.D0
            dg_here%iota_EX = 0.D0
#endif

#ifdef CHEM
            dg_here%iota_IN = 0.D0
            dg_here%iota_EX = 0.D0

            dg_here%iota2_IN = 0.D0
            dg_here%iota2_EX = 0.D0
#endif     


#ifdef DYNP 
            dynP_IN = 0.D0
            dynP_EX = 0.D
#endif

#ifdef TRACE
            TZ_X_IN = 0.D0
            TZ_Y_IN = 0.D0
#endif


#ifdef SED_LAY
            dg_here%bed_IN(:) = 0.D0 
            dg_here%bed_EX(:) = 0.D0
            
            MZ_X_IN(:) = 0.D0
            MZ_Y_IN(:) = 0.D0
#endif

#ifdef WAVE_DIF
            HZ_X_IN = 0.D0
            HZ_Y_IN = 0.D0
#endif

            LZ_XX_IN = 0.D0
            LZ_XY_IN = 0.D0
            LZ_YX_IN = 0.D0
            LZ_YY_IN = 0.D0

            
!.....Compute the specified flow boundaries for the exterior state

            Q_N_EXT = 0.D0
            DO JJ = 1,NFFR
               
               IF(FPER(JJ).EQ.0.D0) THEN
                  NCYC = 0.D0
               ELSE
                  NCYC = INT(dg_here%TIMEDG/FPER(JJ))
               ENDIF

               ARGJ = FAMIG(JJ)*(dg_here%TIMEDG - NCYC*FPER(JJ)) + FFACE(JJ)
               RFF  = FFF(JJ)*dg_here%RAMPDG
               
               dg_here%QNAM_GP = 0.5D0*(dg_here%QNAM_DG(JJ,L,1) + dg_here%QNAM_DG(JJ,L,2))&
              + 0.5D0*(dg_here%QNAM_DG(JJ,L,2) - dg_here%QNAM_DG(JJ,L,1))*dg_here%XEGP(I,dg_here%PA)
               dg_here%QNPH_GP = 0.5D0*(dg_here%QNPH_DG(JJ,L,1) + dg_here%QNPH_DG(JJ,L,2))&
              + 0.5D0*(dg_here%QNPH_DG(JJ,L,2) - dg_here%QNPH_DG(JJ,L,1))*dg_here%XEGP(I,dg_here%PA)
               
               ARG = ARGJ - dg_here%QNPH_GP
               
               Q_N_EXT = Q_N_EXT + dg_here%QNAM_GP*RFF*COS(ARG)
                                !Q_N_EXT = Q_N_EXT + dg_here%QNAM_GP*RFF*COS(ARG)+rff*RAMPExtFlux*0.05D0 
               Q_T_EXT =  0.D0

               dg_here%QX_EX = -( TY*Q_N_EXT - dg_here%NY*Q_T_EXT)/(dg_here%NX*TY - dg_here%NY*TX)
               dg_here%QY_EX = -(-TX*Q_N_EXT + dg_here%NX*Q_T_EXT)/(dg_here%NX*TY - dg_here%NY*TX)

            ENDDO
            
!.....Compute the solution at the interior state

                                !dg_here%HB_IN = dg_here%BATHED(I,LED,EL_IN,dg_here%pa)

            DO K = 1,dg_here%DOFS(EL_IN)

               dg_here%ZE_IN = dg_here%ZE_IN + dg_here%ZE(K,EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               dg_here%QX_IN = dg_here%QX_IN + dg_here%QX(K,EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               dg_here%QY_IN = dg_here%QY_IN + dg_here%QY(K,EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               
#ifdef TRACE    
               dg_here%iota_IN = dg_here%iota_IN + dg_here%iota(K,EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
#endif

#ifdef CHEM
               dg_here%iota_IN =  dg_here%iota_IN  + dg_here%iota(K,EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               dg_here%iota2_IN = dg_here%iota2_IN + dg_here%iota2(K,EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
#endif

#ifdef DYNP    
               dynP_IN = dynP_IN + dg_here%dynP(K,EL_IN,dg_here%IRK)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
#endif

#ifdef SED_LAY
               do ll = 1,s%layers
                  dg_here%bed_IN(ll) = dg_here%bed_IN(ll) + dg_here%bed(K,EL_IN,dg_here%IRK,ll)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
                  dg_here%HB_IN = dg_here%HB_IN + dg_here%bed(k,EL_IN,dg_here%irk,ll)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               
                  !MZ_X_IN(ll) = MZ_X_IN(ll) + dg_here%MZ(K,1,ll,EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
                  !MZ_Y_IN(ll) = MZ_Y_IN(ll) + dg_here%MZ(K,2,ll,EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)     
               enddo
          

#else
               dg_here%HB_IN = dg_here%HB_IN + dg_here%HB(K,EL_IN,1)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
#endif


! Diffusion on BCs?  Induce boundary layers?
                                ! LDG terms

               !HZ_X_IN = HZ_X_IN + dg_here%HZ(K,1,1,EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               !HZ_Y_IN = HZ_Y_IN + dg_here%HZ(K,2,2,EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               !LZ_XX_IN = LZ_XX_IN + dg_here%LZ(K,1,1,EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               !LZ_XY_IN = LZ_XY_IN + dg_here%LZ(K,1,2,EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               !LZ_YX_IN = LZ_YX_IN + dg_here%LZ(K,2,1,EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               !LZ_YY_IN = LZ_YY_IN + dg_here%LZ(K,2,2,EL_IN)*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)

            ENDDO

            !dg_here%QX_EX = -( TY*Q_N_EXT - dg_here%NY*Q_T_EXT)/(dg_here%NX*TY - dg_here%NY*TX)
            !dg_here%QY_EX = -(-TX*Q_N_EXT + dg_here%NX*Q_T_EXT)/(dg_here%NX*TY - dg_here%NY*TX)


#ifdef P0
            IF (dg_here%DOFS(EL_IN).EQ.1) THEN
               DO K = 2,3
                  dg_here%HB_IN = dg_here%HB_IN + dg_here%HB(K,EL_IN,1  )*dg_here%PHI_EDGE(K,I,LED,dg_here%pa)
               ENDDO
            ENDIF
#endif
            
!.....Set the exterior dg_here%bed and surface elevation equal to the interior

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

#ifdef DYNP
            dynP_EX = dg_here%iota_IN
#endif

#ifdef SED_LAY

            dg_here%bed_EX(:) = dg_here%bed_IN(:)

#endif

            CALL NUMERICAL_FLUX(s,dg_here,IT,L)

!.....Add LDG terms
#ifdef WAVE_DIF
            F_HAT = F_HAT + HZ_X_IN*dg_here%NX*dg_here%SFAC_IN + HZ_Y_IN*dg_here%NY
#endif
            G_HAT = G_HAT + LZ_XX_IN*dg_here%NX*dg_here%SFAC_IN + LZ_XY_IN*dg_here%NY
            H_HAT = H_HAT + LZ_YX_IN*dg_here%NX*dg_here%SFAC_IN + LZ_YY_IN*dg_here%NY

#ifdef TRACE
            I_HAT = I_HAT + TZ_X_IN*dg_here%NX*dg_here%SFAC_IN + TZ_Y_IN*dg_here%NY
#endif

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

#ifdef DYNP
               dg_here%RHS_dynP(K,EL_IN,dg_here%IRK) = dg_here%RHS_dynP(K,EL_IN,dg_here%IRK) - W_IN*K_HAT
#endif

#ifdef SED_LAY
               do ll=1,s%layers
                  dg_here%RHS_bed(K,EL_IN,dg_here%IRK,ll) =  dg_here%RHS_bed(K,EL_IN,dg_here%IRK,ll) - W_IN*dg_here%bed_HAT(ll)
               enddo
#endif

            ENDDO

!$$$  DO K = 1,dg_here%DOFS(EL_IN)
!$$$  CALL EDGE_INT_HYDRO(EL_IN, LED, GED, I, F_HAT, G_HAT, H_HAT,i_hat,j_hat,k,dg_here%pa)
!$$$  ENDDO
            
         ENDDO

 1000 CONTINUE
      
      RETURN
      END SUBROUTINE
