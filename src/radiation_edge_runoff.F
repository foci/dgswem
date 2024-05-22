C***********************************************************************
C     
C     SUBROUTINE RADIATION_EDGE_RUNOFF( )
C        
C     This is KWE version of RADIATION_EDGE_HYDRO.
C     
C     Written by Younghun Kang (01-28-2020)
C     
C***********************************************************************

      SUBROUTINE RADIATION_EDGE_RUNOFF(IT)

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      Use SIZES, only: layers

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LED, GED, i,k,jj,II,ll,IT,mm
      Real(SZ) DEN2,U_AVG,V_AVG,VEL_NORMAL,q_RoeX, q_RoeY, q_Roe
      REAL(SZ) TX, TY, DEN,ell_1,ell_2,ell_3,HZ_X_IN,HZ_Y_IN
      REAL(SZ) LZ_XX_IN, LZ_XY_IN, LZ_YX_IN, LZ_YY_IN,W_IN
      Real(SZ) MZ_X_IN(layers),MZ_Y_IN(layers),TZ_X_IN,TZ_Y_IN
      
      REAL(SZ) H_IN, HQX_IN, HQY_IN, HQ_HAT
      REAL(SZ) mX,mY,INFLAG,ONFLAG
      REAL(SZ) CS

      test_el = 0
      DO L = 1,NREDS
         
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
         N1 = NEDNO(1,GED)
         N2 = NEDNO(2,GED)
         
         mX = (X(N1) + X(N2))/2
         mY = (Y(N1) + Y(N2))/2
         
C.....Retrieve the components of the normal vector to the edge
         
            NX = COSNX(GED)
            NY = SINNX(GED)
         
C.....Set the components for the tangential vector to the edge
         
            TX = -NY
            TY =  NX

C.....Compute ZE, QX, QY, and HB at each Gauss point

            DO I = 1,NEGP(pa)
               H_IN = H(1,EL_IN,IRK)
           
               DO K = 2,DOFS(EL_IN)
                  
                  H_IN = H_IN + H(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa)
               
               ENDDO
               
               IF (H_IN .lt. 0.D0) THEN
                  H_IN = 0.D0
               ENDIF
               
               ! Manning_n = 0.015
               HQX_IN = S0X(EL_IN)*(H_IN**(5.0d0/3.0d0)) /
     &          (manningKWE2(EL_IN)*sqrt(S0mag(EL_IN)))
               HQY_IN = S0Y(EL_IN)*(H_IN**(5.0d0/3.0d0)) /
     &          (manningKWE2(EL_IN)*sqrt(S0mag(EL_IN)))

C.....Compute the numerical flux (upwind scheme)
                            
               HQ_HAT = HQX_IN*NX + HQY_IN*NY
               
C.....Compute the edge integral
               DO K = 1,DOFS(EL_IN)
                  W_IN = 2.0*M_INV(K,pa)/AREAS(EL_IN)*XLEN(GED)*
     &                 PHI_EDGE(K,I,LED,pa)*WEGP(I,pa)
                  RHS_KW2(K,EL_IN,IRK) = RHS_KW2(K,EL_IN,IRK) -
     &             W_IN*HQ_HAT
               ENDDO
            
            ENDDO
      ENDDO
      RETURN
      END SUBROUTINE

