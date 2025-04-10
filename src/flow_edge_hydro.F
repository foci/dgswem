C***********************************************************************
C     
C     SUBROUTINE FLOW_EDGE_HYDRO( )
C     
C     This subroutine does the following:
C     
C     1.  Calculates the values of the necessary variables at the edge
C     gauss points for NON-ZERO FLUX edges
C     2.  Calls the appropriate subroutine to compute the flux at
C     these points.
C     3.  Calls the appropriate subroutine to compute the boundary
C     integrals.
C     
C     Written by Ethan Kubatko (06-11-2004)

C     01-10-2011 - cem - adapted for p_enrichment and multicomponent
C     
C***********************************************************************

      SUBROUTINE FLOW_EDGE_HYDRO(IT)

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      use SIZES, only: layers
      use fparser
      use fparser2
      use mod_numerical_flux

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LED, GED, i,k,jj,II,ll,IT,w
      Real(SZ) DEN2,U_AVG,V_AVG,VEL_NORMAL,q_RoeX, q_RoeY, q_Roe
      REAL(SZ) TX, TY,HZ_X_IN,HZ_Y_IN,TZ_X_IN,TZ_Y_IN
      REAL(SZ) LZ_XX_IN, LZ_XY_IN, LZ_YX_IN, LZ_YY_IN,W_IN
      Real(SZ) chi_pref,C_char,MZ_X_IN(layers),MZ_Y_IN(layers)
      real(sz) :: ze_in,ze_ex,qx_in,qx_ex,qy_in,qy_ex,
     &  hb_in,hb_ex,sfac_in,sfac_ex, nx, ny
      real(sz) :: f_hat, g_hat, h_hat

      test_el = 0
!$acc parallel loop gang vector async(1) default(present)
      DO 1000 L = 1,NFEDS
         
C.....Retrieve the global and local edge number

         GED = NFEDN(L)
         LED = NEDSD(1,GED)

C.....Retrieve the element to which the edge belongs

         EL_IN = NEDEL(1,GED)

C.....Organize the dofs for the p_adaptivity

         PA = PDG_EL(EL_IN)

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
         
C.....Compute ZE, QX, QY, and HB at each edge Gauss quadrature point

!$acc loop seq
         DO I = 1,NEGP(pa)

            ZE_IN = 0.0
            QX_IN = 0.0
            QY_IN = 0.0

            ZE_EX = 0.0
            QX_EX = 0.0
            QY_EX = 0.0
            
            HB_IN = 0.0
            SFAC_IN = SFACED(I,LED,EL_IN,pa)

#ifdef TRACE
            iota_IN = 0.0
            iota_EX = 0.0
#endif

#ifdef CHEM
            iota_IN = 0.0
            iota_EX = 0.0

            iota2_IN = 0.0
            iota2_EX = 0.0
#endif     


#ifdef DYNP 
            dynP_IN = 0.0
            dynP_EX = 0.D
#endif

#ifdef TRACE
            TZ_X_IN = 0.0
            TZ_Y_IN = 0.0
#endif


#ifdef SED_LAY
            bed_IN(:) = 0.0 
            bed_EX(:) = 0.0
            
            MZ_X_IN(:) = 0.0
            MZ_Y_IN(:) = 0.0
#endif

#ifdef WAVE_DIF
            HZ_X_IN = 0.0
            HZ_Y_IN = 0.0
#endif

            LZ_XX_IN = 0.0
            LZ_XY_IN = 0.0
            LZ_YX_IN = 0.0
            LZ_YY_IN = 0.0

            
C.....Compute the specified flow boundaries for the exterior state

            Q_N_EXT = 0.0
!$acc loop seq
            DO JJ = 1,MNFFR
               IF (NFFR.EQ.0) THEN
                  ARGJ = 0.0
                  RFF = RAMPDG
               ELSEIF (FPER(JJ).EQ.0.0) THEN
                  NCYC = 0.0
                  ARGJ = FAMIG(JJ)*(TIMEDG - NCYC*FPER(JJ)) + FFACE(JJ)
                  RFF  = FFF(JJ)*RAMPDG
               ELSE
                  NCYC = INT(TIMEDG/FPER(JJ))
                  ARGJ = FAMIG(JJ)*(TIMEDG - NCYC*FPER(JJ)) + FFACE(JJ)
                  RFF  = FFF(JJ)*RAMPDG
               ENDIF
               
               QNAM_GP = 0.5D0*(QNAM_DG(JJ,L,1) + QNAM_DG(JJ,L,2))
     &            + 0.5D0*(QNAM_DG(JJ,L,2) - QNAM_DG(JJ,L,1))*XEGP(I,PA)
               QNPH_GP = 0.5D0*(QNPH_DG(JJ,L,1) + QNPH_DG(JJ,L,2))
     &            + 0.5D0*(QNPH_DG(JJ,L,2) - QNPH_DG(JJ,L,1))*XEGP(I,PA)
               
               ARG = ARGJ - QNPH_GP
               
               Q_N_EXT = Q_N_EXT + QNAM_GP*RFF*COS(ARG)
        !Q_N_EXT = Q_N_EXT + QNAM_GP*RFF*COS(ARG)+rff*RAMPExtFlux*0.05D0 
               Q_T_EXT =  0.0
               
               QX_EX = -( TY*Q_N_EXT - NY*Q_T_EXT)/(NX*TY - NY*TX)
               QY_EX = -(-TX*Q_N_EXT + NX*Q_T_EXT)/(NX*TY - NY*TX)
                 
            ENDDO
            
C.....Compute the solution at the interior state

                                !HB_IN = BATHED(I,LED,EL_IN,pa)

!$acc loop seq
            DO K = 1,DOFS(EL_IN)

               ZE_IN = ZE_IN + ZE(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa)
               QX_IN = QX_IN + QX(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa)
               QY_IN = QY_IN + QY(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa)
               
#ifdef TRACE    
             iota_IN = iota_IN + iota(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa)
#endif

#ifdef CHEM
            iota_IN =  iota_IN  + iota(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa)
           iota2_IN = iota2_IN + iota2(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa)
#endif

#ifdef DYNP    
               dynP_IN = dynP_IN + dynP(K,EL_IN,IRK)
     &                       *PHI_EDGE(K,I,LED,pa)
#endif

#ifdef SED_LAY
               do ll = 1,layers
      bed_IN(ll) = bed_IN(ll) + bed(K,EL_IN,IRK,ll)*PHI_EDGE(K,I,LED,pa)
               HB_IN = HB_IN + bed(k,EL_IN,irk,ll)*PHI_EDGE(K,I,LED,pa)
               
      !MZ_X_IN(ll) = MZ_X_IN(ll) + MZ(K,1,ll,EL_IN)*PHI_EDGE(K,I,LED,pa)
      !MZ_Y_IN(ll) = MZ_Y_IN(ll) + MZ(K,2,ll,EL_IN)*PHI_EDGE(K,I,LED,pa)     
               enddo
          

#else
               HB_IN = HB_IN + HB(K,EL_IN,1)*PHI_EDGE(K,I,LED,pa)
#endif


! Diffusion on BCs?  Induce boundary layers?
                                ! LDG terms

             !HZ_X_IN = HZ_X_IN + HZ(K,1,1,EL_IN)*PHI_EDGE(K,I,LED,pa)
             !HZ_Y_IN = HZ_Y_IN + HZ(K,2,2,EL_IN)*PHI_EDGE(K,I,LED,pa)
             !LZ_XX_IN = LZ_XX_IN + LZ(K,1,1,EL_IN)*PHI_EDGE(K,I,LED,pa)
             !LZ_XY_IN = LZ_XY_IN + LZ(K,1,2,EL_IN)*PHI_EDGE(K,I,LED,pa)
             !LZ_YX_IN = LZ_YX_IN + LZ(K,2,1,EL_IN)*PHI_EDGE(K,I,LED,pa)
             !LZ_YY_IN = LZ_YY_IN + LZ(K,2,2,EL_IN)*PHI_EDGE(K,I,LED,pa)

            ENDDO

            !QX_EX = -( TY*Q_N_EXT - NY*Q_T_EXT)/(NX*TY - NY*TX)
            !QY_EX = -(-TX*Q_N_EXT + NX*Q_T_EXT)/(NX*TY - NY*TX)


#ifdef P0
            IF (DOFS(EL_IN).EQ.1) THEN
               DO K = 2,3
                  HB_IN = HB_IN + HB(K,EL_IN,1  )*PHI_EDGE(K,I,LED,pa)
               ENDDO
            ENDIF
#endif
            
C.....Set the exterior bed and surface elevation equal to the interior

            ZE_EX = ZE_IN
            HB_EX = HB_IN

            SFAC_EX = SFAC_IN

#ifdef TRACE
            iota_EX = iota_IN
#endif

#ifdef CHEM
            iota_EX = iota_IN
            iota2_EX = iota2_IN 
#endif

#ifdef DYNP
            dynP_EX = iota_IN
#endif

#ifdef SED_LAY

            bed_EX(:) = bed_IN(:)

#endif

            call numerical_flux(ze_in,ze_ex,qx_in,qx_ex,qy_in,qy_ex,
     &  hb_in,hb_ex,sfac_in,sfac_ex,nx,ny,f_hat,g_hat,h_hat)

C.....Add LDG terms
#ifdef WAVE_DIF
            F_HAT = F_HAT + HZ_X_IN*NX*SFAC_IN + HZ_Y_IN*NY
#endif
            G_HAT = G_HAT + LZ_XX_IN*NX*SFAC_IN + LZ_XY_IN*NY
            H_HAT = H_HAT + LZ_YX_IN*NX*SFAC_IN + LZ_YY_IN*NY

#ifdef TRACE
            I_HAT = I_HAT + TZ_X_IN*NX*SFAC_IN + TZ_Y_IN*NY
#endif

#ifdef SED_LAY
            do ll=1,layers
      bed_HAT(ll)= bed_HAT(ll) + MZ_X_IN(ll)*NX*SFAC_IN + MZ_Y_IN(ll)*NY
            enddo
#endif
            
C.....Compute the edge integral

!$acc loop seq
            DO K = 1,DOFS(EL_IN)

               W_IN = 2.0*M_INV(K,pa)/AREAS(EL_IN)*XLEN(GED)*
     &              PHI_EDGE(K,I,LED,pa)*WEGP(I,pa)

               RHS_ZE(K,EL_IN,IRK) = RHS_ZE(K,EL_IN,IRK) - W_IN*F_HAT
               RHS_QX(K,EL_IN,IRK) = RHS_QX(K,EL_IN,IRK) - W_IN*G_HAT
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
               do ll=1,layers
                  RHS_bed(K,EL_IN,IRK,ll) =  RHS_bed(K,EL_IN,IRK,ll) -
     & W_IN*bed_HAT(ll)
               enddo
#endif

            ENDDO

c$$$  DO K = 1,DOFS(EL_IN)
c$$$  CALL EDGE_INT_HYDRO(EL_IN, LED, GED, I, F_HAT, G_HAT, H_HAT,i_hat,
c$$$     &                                                   j_hat,k,pa)
c$$$  ENDDO
            
         ENDDO

 1000 CONTINUE
      
      RETURN
      END SUBROUTINE
