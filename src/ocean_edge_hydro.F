C***********************************************************************
C     
C     SUBROUTINE OCEAN_EDGE_HYDRO( )
C     
C     This subroutine does the following:
C     
C     1.  Calculates the values of the necessary variables at the edge
C     gauss points for ELEVATION SPECIFIED edges
C     2.  Calls the appropriate subroutine to compute the flux at
C     these points.
C     3.  Calls the appropriate subroutine to compute the boundary
C     integrals.
C     
C     Written by Ethan Kubatko (06-11-2004)
C     
C     01-10-2011 - cem - adapted for p_enrichment and multicomponent  
C     
C***********************************************************************
      SUBROUTINE OCEAN_EDGE_HYDRO(IT)

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      Use SIZES, only: layers
      use mod_numerical_flux
c     
      USE NodalAttributes, ONLY: GeoidOffset, LoadGeoidOffset
      use fparser
      use fparser2

      IMPLICIT NONE

C.....Declare local variables

      real(sz) :: ze_in,ze_ex,qx_in,qx_ex,qy_in,qy_ex,
     &  hb_in,hb_ex,sfac_in,sfac_ex, nx, ny
      real(sz) :: f_hat, g_hat, h_hat
      INTEGER L, LED, GED, i,k,jj,II,ll,IT,w
      Real(SZ) DEN2,U_AVG,V_AVG,VEL_NORMAL,q_RoeX, q_RoeY, q_Roe
      REAL(SZ) TX, TY, HUU, HUV, GH2,FH_NL_IN,F1_NL,FX1_IN,FY1_IN
      REAL(SZ) LZ_XX_IN, LZ_XY_IN, LZ_YX_IN, LZ_YY_IN,W_IN,FX3_IN
      Real(SZ) chi_pref,MZ_X_IN(layers),MZ_Y_IN(layers)
      Real(SZ) HZ_X_IN,HZ_Y_IN,TZ_X_IN,TZ_Y_IN

      INTEGER:: IPT
      REAL(SZ):: ZEFREQ      
      REAL(SZ), dimension(2):: EFA_GPT, EMO_GPT, ARG_GPT, ZE_GPT

!$acc parallel loop gang vector private(efa_gpt, emo_gpt, arg_gpt, ze_gpt) async(1) default(present)
      DO 1000 L = 1, needs
         
C.....Retrieve the global and local edge number

         GED = NEEDN(L)
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
         
C.....Retrieve the nodes of the edge
         
         N1 = NEDNO(1,GED)
         N2 = NEDNO(2,GED)
         
C.....Compute ZE, QX, QY, and HB at each edge Gauss quadrature point
!$acc loop seq
         DO I = 1,NEGP(pa)

            ZE_IN = ZE(1,EL_IN,IRK)
            QX_IN = QX(1,EL_IN,IRK)
            QY_IN = QY(1,EL_IN,IRK)

            ZE_EX = 0.0
            QX_EX = 0.0
            QY_EX = 0.0
            U_EX  = 0.0
            V_EX  = 0.0

#ifdef TRACE
            iota_IN = iota(1,EL_IN,IRK)
            iota_EX = 0.0
#endif

#ifdef CHEM
            iota_IN = iota(1,EL_IN,IRK)
            iota2_IN = iota2(1,EL_IN,IRK)
            iota_EX = 0.0
            iota2_EX = 0.0
#endif

#ifdef DYNP
            dynP_IN = dynP(1,EL_IN,IRK)
            dynP_EX = 0.0
#endif
            
            HB_IN = BATHED(I,LED,EL_IN,pa)
            SFAC_IN = SFACED(I,LED,EL_IN,pa)
            
#ifdef SED_LAY                  !When layered, these change
            HB(:,EL_IN,irk) = 0.0
            do ll=1,layers
               HB(1,EL_IN,irk) = HB(1,EL_IN,irk) + bed(1,EL_IN,irk,ll)

               MZ_X_IN(ll) =  MZ(1,1,ll,EL_IN)
               MZ_Y_IN(ll) =  MZ(1,2,ll,EL_IN)
            enddo
            bed_IN(:) = bed(1,EL_IN,irk,:)
            HB_IN = HB(1,EL_IN,irk)
#endif
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

C.....Compute the specified open ocean elevation

!$acc loop seq
            DO JJ=1,NBFR
               
               IF (PER(JJ).EQ.0.0) THEN
                  NCYC = 0.0
               ELSE
                  NCYC = INT(TIMEDG/PER(JJ))
               ENDIF
               
C...........Surface Elevation
               
               ARGJ = AMIG(JJ)*(TIMEDG - NCYC*PER(JJ)) + FACE(JJ)
               RFF = FF(JJ)*RAMPDG

C..............linearly interpolate from the time-series of the harmonic forcing
!$acc loop seq
               DO IPT = 1, 2
                  EFA_GPT(IPT) = EFA_DG(JJ,L,IPT) ;
                  EMO_GPT(IPT) = EMO_DG(JJ,L,IPT) ; 

                  ARG_GPT(IPT) = ARGJ - EFA_GPT(IPT) ;

                  ZE_GPT(IPT) = EMO_GPT(IPT)*RFF*COS(ARG_GPT(IPT)) ;
               END DO

               ZEFREQ = 0.5D0*( ZE_GPT(1) + ZE_GPT(2) ) +
     &              0.5D0*(ZE_GPT(2) - ZE_GPT(1))*XEGP(I,pa) ;


               ZE_EX = ZE_EX + ZEFREQ ;
            ENDDO


C.....For non-periodic elevation bcs 

            IF((NBFR.EQ.0).AND.(NOPE.GT.0)) THEN
              ZE_EX = 0.5D0*(ETA2(N1) + ETA2(N2))
     &              + 0.5D0*(ETA2(N2) - ETA2(N1))*XEGP(I,pa)
            ENDIF


C.....Compute the solution at the interior state

!$acc loop seq
            DO K = 2,DOFS(EL_IN)

               ZE_IN = ZE_IN + ZE(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa)
               QX_IN = QX_IN + QX(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa)
               QY_IN = QY_IN + QY(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa)

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
     & bed(K,EL_IN,IRK,ll)*PHI_EDGE(K,I,LED,pa)
                  HB_IN = HB_IN + bed(k,EL_IN,irk,ll)*
     & PHI_EDGE(K,I,LED,pa)

                  MZ_X_IN(ll) = MZ_X_IN(ll) +
     & MZ(K,1,ll,EL_IN)*PHI_EDGE(K,I,LED,pa)
                  MZ_Y_IN(ll) = MZ_Y_IN(ll) +
     & MZ(K,2,ll,EL_IN)*PHI_EDGE(K,I,LED,pa)
               enddo
#endif

                                ! LDG terms
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

            ENDDO

C.....Set the exterior value of the bathymetry equal to the interior

            HB_EX = HB_IN

#ifdef SED_LAY
            bed_EX(:) = bed_IN(:)
#endif
            SFAC_EX = SFAC_IN

            IF (LoadGeoidOffset) then
               ZE_EX = ZE_EX + .5*(GeoidOffset(N1)+GeoidOffset(N2))
            endif

            ! Eirik's fix
            if ((ZE_EX*IFNLFA + HB_EX).LE.0) then
              ZE_EX =  ABS(HB_EX) + H0
            endif
            ! End fix

            !ZE_EX = ZE_IN !this was added

#ifdef TRACE
            iota_EX = 0.0 !iota_IN !might need 0.0
#endif

#ifdef CHEM
            iota_EX = 0.0 !iota_IN
            iota2_EX = 0.0 !iota2_IN
#endif

#ifdef DYNP
            dynP_EX = 0.0 !What's the right setting here?
#endif

#ifdef SED_LAY                  !When layered, these change
            bed_EX(:) = 10.0 !bed_IN(:) 
#endif


C.....Set the exterior state flows equal to the interior state flows

            QX_EX = QX_IN
            QY_EX = QY_IN

C.....Compute the Roe flux

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

C.....Add LDG terms for sediment

#ifdef SED_LAY
            do ll=1,layers
               bed_HAT(ll) = bed_HAT(ll) +
     & MZ_X_IN(ll)*NX*SFAC_IN + MZ_Y_IN(ll)*NY
            enddo
#endif

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
     & W_IN*bed_HAT(ll)
               enddo
#endif

            ENDDO
         ENDDO

 1000 CONTINUE
      
      RETURN
      END SUBROUTINE
