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

      SUBROUTINE FLOW_EDGE_HYDRO(s,IT)

!.....Use appropriate modules

      USE GLOBAL
      USE DG
      use SIZES
      use fparser
      use fparser2

      IMPLICIT NONE

      type (sizes_type) :: s

!.....Declare local variables

      INTEGER L, LED, GED, i,k,jj,II,ll,IT,w
      Real(SZ) DEN2,U_AVG,V_AVG,VEL_NORMAL,q_RoeX, q_RoeY, q_Roe
      REAL(SZ) TX, TY,HZ_X_IN,HZ_Y_IN,TZ_X_IN,TZ_Y_IN
      REAL(SZ) LZ_XX_IN, LZ_XY_IN, LZ_YX_IN, LZ_YY_IN,W_IN
      Real(SZ) VEI(3,3),chi_pref,C_char,MZ_X_IN(s%layers),MZ_Y_IN(s%layers)

      dg%test_el = 0
      DO 1000 L = 1,dg%NFEDS
         
!.....Retrieve the global and local edge number

         GED = dg%NFEDN(L)
         LED = dg%NEDSD(1,GED)

!.....Retrieve the element to which the edge belongs

         EL_IN = dg%NEDEL(1,GED)

!.....Organize the dg%dofs for the p_adaptivity

         dg%PA = PDG_EL(EL_IN)

#ifdef P0
         if (dg%pa.eq.0) then
            dg%pa = 1
         endif
#endif

!.....If the element is dry then skip the edge calculation

         IF (dg%WDFLG(EL_IN).EQ.0) GOTO 1000
         
!.....Retrieve the components of the normal vector to the edge
         
         dg%NX = dg%COSNX(GED)
         dg%NY = dg%SINNX(GED)
         
!.....Set the components for the tangential vector to the edge

         TX = -dg%NY
         TY =  dg%NX
         
!.....Compute dg%ZE, dg%QX, dg%QY, and dg%HB at each edge Gauss quadrature point

         DO I = 1,dg%NEGP(dg%pa)

            dg%ZE_IN = 0.D0
            dg%QX_IN = 0.D0
            dg%QY_IN = 0.D0

            dg%ZE_EX = 0.D0
            dg%QX_EX = 0.D0
            dg%QY_EX = 0.D0
            
            dg%HB_IN = 0.D0
            dg%SFAC_IN = dg%SFACED(I,LED,EL_IN,dg%pa)

#ifdef TRACE
            dg%iota_IN = 0.D0
            dg%iota_EX = 0.D0
#endif

#ifdef CHEM
            dg%iota_IN = 0.D0
            dg%iota_EX = 0.D0

            dg%iota2_IN = 0.D0
            dg%iota2_EX = 0.D0
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
            dg%bed_IN(:) = 0.D0 
            dg%bed_EX(:) = 0.D0
            
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
                  NCYC = INT(dg%TIMEDG/FPER(JJ))
               ENDIF

               ARGJ = FAMIG(JJ)*(dg%TIMEDG - NCYC*FPER(JJ)) + FFACE(JJ)
               RFF  = FFF(JJ)*dg%RAMPDG
               
               dg%QNAM_GP = 0.5D0*(dg%QNAM_DG(JJ,L,1) + dg%QNAM_DG(JJ,L,2))&
              + 0.5D0*(dg%QNAM_DG(JJ,L,2) - dg%QNAM_DG(JJ,L,1))*dg%XEGP(I,dg%PA)
               dg%QNPH_GP = 0.5D0*(dg%QNPH_DG(JJ,L,1) + dg%QNPH_DG(JJ,L,2))&
              + 0.5D0*(dg%QNPH_DG(JJ,L,2) - dg%QNPH_DG(JJ,L,1))*dg%XEGP(I,dg%PA)
               
               ARG = ARGJ - dg%QNPH_GP
               
               Q_N_EXT = Q_N_EXT + dg%QNAM_GP*RFF*COS(ARG)
                                !Q_N_EXT = Q_N_EXT + dg%QNAM_GP*RFF*COS(ARG)+rff*RAMPExtFlux*0.05D0 
               Q_T_EXT =  0.D0

               dg%QX_EX = -( TY*Q_N_EXT - dg%NY*Q_T_EXT)/(dg%NX*TY - dg%NY*TX)
               dg%QY_EX = -(-TX*Q_N_EXT + dg%NX*Q_T_EXT)/(dg%NX*TY - dg%NY*TX)

            ENDDO
            
!.....Compute the solution at the interior state

                                !dg%HB_IN = dg%BATHED(I,LED,EL_IN,dg%pa)

            DO K = 1,dg%DOFS(EL_IN)

               dg%ZE_IN = dg%ZE_IN + dg%ZE(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
               dg%QX_IN = dg%QX_IN + dg%QX(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
               dg%QY_IN = dg%QY_IN + dg%QY(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
               
#ifdef TRACE    
               dg%iota_IN = dg%iota_IN + dg%iota(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
#endif

#ifdef CHEM
               dg%iota_IN =  dg%iota_IN  + dg%iota(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
               dg%iota2_IN = dg%iota2_IN + dg%iota2(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
#endif

#ifdef DYNP    
               dynP_IN = dynP_IN + dg%dynP(K,EL_IN,dg%IRK)*dg%PHI_EDGE(K,I,LED,dg%pa)
#endif

#ifdef SED_LAY
               do ll = 1,s%layers
                  dg%bed_IN(ll) = dg%bed_IN(ll) + dg%bed(K,EL_IN,dg%IRK,ll)*dg%PHI_EDGE(K,I,LED,dg%pa)
                  dg%HB_IN = dg%HB_IN + dg%bed(k,EL_IN,dg%irk,ll)*dg%PHI_EDGE(K,I,LED,dg%pa)
               
                  !MZ_X_IN(ll) = MZ_X_IN(ll) + dg%MZ(K,1,ll,EL_IN)*dg%PHI_EDGE(K,I,LED,dg%pa)
                  !MZ_Y_IN(ll) = MZ_Y_IN(ll) + dg%MZ(K,2,ll,EL_IN)*dg%PHI_EDGE(K,I,LED,dg%pa)     
               enddo
          

#else
               dg%HB_IN = dg%HB_IN + dg%HB(K,EL_IN,1)*dg%PHI_EDGE(K,I,LED,dg%pa)
#endif


! Diffusion on BCs?  Induce boundary layers?
                                ! LDG terms

               !HZ_X_IN = HZ_X_IN + dg%HZ(K,1,1,EL_IN)*dg%PHI_EDGE(K,I,LED,dg%pa)
               !HZ_Y_IN = HZ_Y_IN + dg%HZ(K,2,2,EL_IN)*dg%PHI_EDGE(K,I,LED,dg%pa)
               !LZ_XX_IN = LZ_XX_IN + dg%LZ(K,1,1,EL_IN)*dg%PHI_EDGE(K,I,LED,dg%pa)
               !LZ_XY_IN = LZ_XY_IN + dg%LZ(K,1,2,EL_IN)*dg%PHI_EDGE(K,I,LED,dg%pa)
               !LZ_YX_IN = LZ_YX_IN + dg%LZ(K,2,1,EL_IN)*dg%PHI_EDGE(K,I,LED,dg%pa)
               !LZ_YY_IN = LZ_YY_IN + dg%LZ(K,2,2,EL_IN)*dg%PHI_EDGE(K,I,LED,dg%pa)

            ENDDO

            !dg%QX_EX = -( TY*Q_N_EXT - dg%NY*Q_T_EXT)/(dg%NX*TY - dg%NY*TX)
            !dg%QY_EX = -(-TX*Q_N_EXT + dg%NX*Q_T_EXT)/(dg%NX*TY - dg%NY*TX)


#ifdef P0
            IF (dg%DOFS(EL_IN).EQ.1) THEN
               DO K = 2,3
                  dg%HB_IN = dg%HB_IN + dg%HB(K,EL_IN,1  )*dg%PHI_EDGE(K,I,LED,dg%pa)
               ENDDO
            ENDIF
#endif
            
!.....Set the exterior dg%bed and surface elevation equal to the interior

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

#ifdef DYNP
            dynP_EX = dg%iota_IN
#endif

#ifdef SED_LAY

            dg%bed_EX(:) = dg%bed_IN(:)

#endif

            CALL NUMERICAL_FLUX(s,IT,L)

!.....Add LDG terms
#ifdef WAVE_DIF
            F_HAT = F_HAT + HZ_X_IN*dg%NX*dg%SFAC_IN + HZ_Y_IN*dg%NY
#endif
            G_HAT = G_HAT + LZ_XX_IN*dg%NX*dg%SFAC_IN + LZ_XY_IN*dg%NY
            H_HAT = H_HAT + LZ_YX_IN*dg%NX*dg%SFAC_IN + LZ_YY_IN*dg%NY

#ifdef TRACE
            I_HAT = I_HAT + TZ_X_IN*dg%NX*dg%SFAC_IN + TZ_Y_IN*dg%NY
#endif

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

#ifdef DYNP
               dg%RHS_dynP(K,EL_IN,dg%IRK) = dg%RHS_dynP(K,EL_IN,dg%IRK) - W_IN*K_HAT
#endif

#ifdef SED_LAY
               do ll=1,s%layers
                  dg%RHS_bed(K,EL_IN,dg%IRK,ll) =  dg%RHS_bed(K,EL_IN,dg%IRK,ll) - W_IN*dg%bed_HAT(ll)
               enddo
#endif

            ENDDO

!$$$  DO K = 1,dg%DOFS(EL_IN)
!$$$  CALL EDGE_INT_HYDRO(EL_IN, LED, GED, I, F_HAT, G_HAT, H_HAT,i_hat,j_hat,k,dg%pa)
!$$$  ENDDO
            
         ENDDO

 1000 CONTINUE
      
      RETURN
      END SUBROUTINE
