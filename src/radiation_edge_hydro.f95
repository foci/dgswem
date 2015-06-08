C***********************************************************************
C     
C     SUBROUTINE RADIATION_EDGE_HYDRO( )
C     
C     This subroutine does the following:
C     
C     1.  Calculates the values of the necessary variables at the edge
C     gauss points for RADIATION edges
C     2.  Calls the appropriate subroutine to compute the flux at
C     these points.
C     3.  Compute the boundary integrals.
C     
C     Written by Ethan Kubatko (06-11-2004)
C     01-10-2011 - cem - adapted for p_enrichment and multicomponent
C     You know, this file is not really up to date.  
C     If you need to use it, you might fix it first -- cem
C     
C***********************************************************************

      SUBROUTINE RADIATION_EDGE_HYDRO(IT)

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE sizes, ONLY: myproc,layers
      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LED, GED,i,k,ll,IT
      REAL(SZ) TX, TY, W_IN, DEN
      REAL(SZ) LZ_XX_IN, LZ_XY_IN, LZ_YX_IN, LZ_YY_IN
      Real(sz) HZ_X_IN,HZ_Y_IN,TZ_X_IN,TZ_Y_IN
      Real(SZ) MZ_X_IN(layers),MZ_Y_IN(layers)

      test_el = 0
      DO 1000 L = 1,NREDS
         
C.....Retrieve the global and local edge number

         GED = NREDN(L)
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
         test_el = test_el+1

C.....Retrieve the components of the normal vector to the edge
         
         NX = COSNX(GED)
         NY = SINNX(GED)
C.....Set the components for the tangential vector to the edge
         
         TX = -NY
         TY =  NX

C.....Compute ZE, QX, QY, and HB at each Gauss point

         DO I = 1,NEGP(pa)

            ZE_IN = ZE(1,EL_IN,IRK)
            QX_IN = QX(1,EL_IN,IRK)
            QY_IN = QY(1,EL_IN,IRK)
            HB_IN = BATHED(I,LED,EL_IN,pa)

            SFAC_IN = SFACED(I,LED,EL_IN,pa)
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

            iota_IN = iota(1,EL_IN,IRK)
#endif

#ifdef CHEM
            iota_IN = iota(1,EL_IN,IRK)
            iota2_IN = iota2(1,EL_IN,IRK)
#endif

#ifdef DYNP
            dynP_IN = dynP(1,EL_IN,IRK)
#endif

#ifdef SED_LAY                  !When layered, these change
            HB(1,EL_IN,irk) = 0.D0
            do ll=1,layers
               HB(1,EL_IN,irk) = HB(1,EL_IN,irk) + bed(1,EL_IN,irk,ll)

               MZ_X_IN(ll) =  MZ(K,1,ll,EL_IN)
               MZ_Y_IN(ll) =  MZ(K,2,ll,EL_IN)
            enddo
            bed_IN(:) = bed(1,EL_IN,irk,:)
            HB_IN = HB(1,EL_IN,irk)
#endif

C.....Compute the solution at the interior state

            DO K = 2,DOFS(EL_IN)
               
               ZE_IN = ZE_IN + ZE(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa)
               QX_IN = QX_IN + QX(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa)
               QY_IN = QY_IN + QY(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa)

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
                  bed_IN(ll) = bed_IN(ll) + bed(K,EL_IN,IRK,ll)*PHI_EDGE(K,I,LED,pa)
                  HB_IN = HB_IN + bed(k,EL_IN,irk,ll)*PHI_EDGE(K,I,LED,pa)

                  MZ_X_IN(ll) = MZ_X_IN(ll) + MZ(K,1,ll,EL_IN)*PHI_EDGE(K,I,LED,pa)
                  MZ_Y_IN(ll) = MZ_Y_IN(ll) + MZ(K,2,ll,EL_IN)*PHI_EDGE(K,I,LED,pa)
               enddo
#endif

            ENDDO

C.....Compute the velocity in the normal and tangental direction
            
            Q_N_INT = QX_IN*NX + QY_IN*NY
            Q_T_INT = QX_IN*TX + QY_IN*TY

C.....Reflect the velocity in the normal direction

c     Q_N_EXT = -Q_N_INT
            Q_N_EXT = -Q_N_INT+2*SQRT(G/(HB_IN+ZE_IN))*ZE_IN
            Q_T_EXT =  Q_T_INT
            
C.....Compute the x and y components of the external state flow

            DEN = 1.D0/(NX*TY - NY*TX)
            QX_EX = ( TY*Q_N_EXT - NY*Q_T_EXT)*DEN
            QY_EX = (-TX*Q_N_EXT + NX*Q_T_EXT)*DEN

            ZE_EX = ZE_IN
            HB_EX = HB_IN

            SFAC_EX = SFAC_IN

#ifdef SED_LAY
               do ll = 1,layers
                  bed_EX(ll) = bed_IN(ll)
               enddo
#endif

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

C.....Compute the flux

            CALL NUMERICAL_FLUX(IT,test_el)
            
C.....Compute the edge integral

            
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
               do ll = 1,layers
                  RHS_bed(K,EL_IN,IRK,ll) = RHS_bed(K,EL_IN,IRK,ll) - W_IN*bed_HAT(ll)
               enddo
#endif

               
c     !CALL EDGE_INT_HYDRO(EL_IN, LED, GED, I, F_HAT, G_HAT, H_HAT,k)

            ENDDO

         ENDDO

 1000 CONTINUE      
      RETURN
      END SUBROUTINE
