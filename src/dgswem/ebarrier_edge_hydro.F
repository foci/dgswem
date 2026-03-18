      SUBROUTINE EBARRIER_EDGE_HYDRO(IT)
      
C.....Use appropriate modules

      USE GLOBAL
      USE DG
      use mod_numerical_flux

      IMPLICIT NONE

C.....Declare local variables

      INTEGER GED, LED, L,i,k,ll,IT
      REAL(SZ) ABOVE, SUBSUP, W_IN, TX, TY
      real(sz) :: ze_in,ze_ex,qx_in,qx_ex,qy_in,qy_ex,
     &  hb_in,hb_ex,sfac_in,sfac_ex, nx, ny
      real(sz) :: f_hat, g_hat, h_hat, fluxes(3)

C.....Loop over the external barrier segments

!$acc parallel loop gang vector private(fluxes) async(1) default(present)
      DO 1000 L = 1,NEBSEG
         
C.......Obtain the global and local edge numbers

         GED = NEBSEGN(L)
         LED = NEDSD(1,GED)

C.......Obtain the element

         EL = NEDEL(1,GED)
         PA = PDG_EL(EL)

#ifdef P0
         if (pa.eq.0) then
            pa = 1
         endif
#endif
         
         IF (WDFLG(EL).EQ.0) GOTO 1000

C.....Retrieve the components of the normal vector to the edge

         NX = COSNX(GED)
         NY = SINNX(GED)
         
C.....Set the components for the tangential vector to the edge

         TX = -NY
         TY =  NX
         
C.......Compute the variables at the quadrature points

!$acc loop seq
         DO I = 1,NEGP(pa)
            
C.........Obtain the height of the barrier at the quadrature point

            ZE_IN = ZE(1,EL,IRK)
            QX_IN = QX(1,EL,IRK)
            QY_IN = QY(1,EL,IRK)
            HB_IN = BATHED(I,LED,EL,pa)

#ifdef TRACE
            iota_IN = iota(1,EL,IRK)
#endif

#ifdef CHEM
            iota_IN = iota(1,EL,IRK)
            iota2_IN = iota2(1,EL,IRK)
#endif

#ifdef DYNP
            dynP_IN = dynP(1,EL,IRK)
#endif

#ifdef SED_LAY
            HB(:,EL,irk) = 0.0
            do ll=1,layers
               HB(1,EL,irk) = HB(1,EL,irk) + bed(1,EL,irk,ll)
            enddo
            bed_IN(:) = bed(1,EL,irk,:)
            HB_IN = HB(1,EL,irk)
#endif


!$acc loop seq
            DO K = 2,DOFS(EL)
               ZE_IN = ZE_IN + ZE(K,EL,IRK)*PHI_EDGE(K,I,LED,pa)
               QX_IN = QX_IN + QX(K,EL,IRK)*PHI_EDGE(K,I,LED,pa)
               QY_IN = QY_IN + QY(K,EL,IRK)*PHI_EDGE(K,I,LED,pa)
               !HB_IN = HB_IN + HB(K,EL,IRK)*PHI_EDGE(K,I,LED,pa)

#ifdef TRACE
               iota_IN = iota_IN + iota(K,EL,IRK)*PHI_EDGE(K,I,LED,pa)
#endif

#ifdef CHEM
               iota_IN = iota_IN + iota(K,EL,IRK)*PHI_EDGE(K,I,LED,pa)
               iota2_IN = iota2_IN +iota2(K,EL,IRK)*PHI_EDGE(K,I,LED,pa)
#endif

#ifdef DYNP
               dynP_IN = dynP_IN + dynP(K,EL,IRK)*PHI_EDGE(K,I,LED,pa)
#endif

#ifdef SED_LAY

            do ll = 1,layers
               
               bed_IN(ll) = bed_IN(ll) + 
     &            bed(K,EL,IRK,ll)*PHI_EDGE(K,I,LED,pa)
               HB_IN = HB_IN + bed(k,EL,irk,ll)*PHI_EDGE(K,I,LED,pa)
               
            enddo
#endif

            ENDDO

            SFAC_IN = SFACED(I,LED,EL,pa)

            ABOVE  = ZE_IN - EBHT(L)
            SUBSUP = 2.0*ABOVE/3.0
            
C.........Case 1:  Water is below barrier
C     ---------------------------------------------
            
            IF (ABOVE.LT.BARMIN) THEN
               Q_N_EXT = -(QX_IN*NX + QY_IN*NY)
               Q_T_EXT =   QX_IN*TX + QY_IN*TY
               
C.........Case 2:  Water is above barrier
C     ------------------------------------------------------------
               
            ELSE
               Q_N_EXT = RAMPDG*EBCFSP(L)*SUBSUP*SQRT(SUBSUP*G)
               Q_T_EXT = 0.0
            ENDIF
            
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
            dynP_EX = dynP_IN
#endif

#ifdef SED_LAY

            bed_EX(:) = bed_IN(:)

#endif

C.........Compute the x and y components of the external state flow

            QX_EX = ( TY*Q_N_EXT - NY*Q_T_EXT)/(NX*TY - NY*TX)
            QY_EX = (-TX*Q_N_EXT + NX*Q_T_EXT)/(NX*TY - NY*TX)
            
C.........Compute numerical flux
            
            fluxes = numerical_flux(ze_in,ze_ex,qx_in,qx_ex,qy_in,qy_ex,
     &  hb_in,hb_ex,sfac_in,sfac_ex,nx,ny,G)
            f_hat = fluxes(1)
            g_hat = fluxes(2)
            h_hat = fluxes(3)

!$acc loop seq
            DO K = 1,DOFS(EL)

               W_IN = 2.0*M_INV(K,pa)/AREAS(EL)*XLEN(GED)
     &              *PHI_EDGE(K,I,LED,pa)*WEGP(I,pa)

!$acc atomic
               RHS_ZE(K,EL,IRK) = RHS_ZE(K,EL,IRK) - W_IN*F_HAT
!$acc atomic
               RHS_QX(K,EL,IRK) = RHS_QX(K,EL,IRK) - W_IN*G_HAT
!$acc atomic
               RHS_QY(K,EL,IRK) = RHS_QY(K,EL,IRK) - W_IN*H_HAT

#ifdef TRACE
               RHS_iota(K,EL,IRK) = RHS_iota(K,EL,IRK) - W_IN*I_HAT
#endif

#ifdef CHEM
               RHS_iota(K,EL,IRK) = RHS_iota(K,EL,IRK) - W_IN*I_HAT
               RHS_iota2(K,EL,IRK) = RHS_iota2(K,EL,IRK) - W_IN*J_HAT
#endif

#ifdef DYNP
               RHS_dynP(K,EL,IRK) = RHS_dynP(K,EL,IRK) - W_IN*I_HAT
#endif

#ifdef SED_LAY

               do ll=1,layers

                  RHS_bed(K,EL,IRK,ll) = RHS_bed(K,EL,IRK,ll) 
     &               - W_IN*Bed_HAT(ll)

               enddo

#endif

            ENDDO
         ENDDO
 1000 CONTINUE

      RETURN
      END SUBROUTINE EBARRIER_EDGE_HYDRO
