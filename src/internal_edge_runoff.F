C***********************************************************************
C     
C     SUBROUTINE INTERNAL_EDGE_RUNOFF( )
C     
C     This is KWE version of INTERNAL_EDGE_HYDRO using simple upwind
C     numerical flux.
C     
C     Written by Younghun Kang (01-28-2020)
C
C***********************************************************************

      SUBROUTINE INTERNAL_EDGE_RUNOFF(IT)

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      USE NodalAttributes, ONLY :  ESLM
      use sizes, only:  sz, myproc, layers

      IMPLICIT NONE

C.....Declare local variables

      INTEGER L, LED_IN, LED_EX, GED, GP_IN, GP_EX,k,i,ll,IT
      REAL(SZ), PARAMETER :: ZERO = 1.D-12
      REAL(SZ) W_IN, W_EX
      REAL(SZ) EDFAC_IN, EDFAC_EX

      REAL(SZ) ARK, BRK
      REAL(SZ) MAX_BOA          ! Maximum of beta_il/alpha_il for all l
      REAL(SZ) H_IN, H_EX, HQX_IN, HQY_IN, HQX_EX, HQY_EX, HQ_HAT_IN, HQ_HAT_EX
c     

      test_el = 0
      
      DO L = 1,NIEDS

C.......Retrieve the global and local edge number

         GED = NIEDN(L)

C.......Retrieve the elements which share the edge
         ! Check direction of edge to make sure water flows from EL_IN to EL_EX
         IF (flipNEDEL(GED).eq.0) THEN
            LED_IN = NEDSD(1,GED)
            LED_EX = NEDSD(2,GED)
            
            EL_IN = NEDEL(1,GED)
            EL_EX = NEDEL(2,GED)
            
            EL = EL_EX

            if (DOFS(EL_EX).LT.DOFS(EL_IN)) then
               EL = EL_IN
            endif
         
            pa = PDG_EL(EL)

#ifdef P0
            if (pa.eq.0) then
               pa = 1
            endif
#endif

C.....Retrieve the components of the normal vector to the edge
        
            NX = COSNX(GED)
            NY = SINNX(GED)

C.....Set the components for the tangential vector to the edge

            EDFAC_IN = XLEN(GED)/AREAS(EL_IN)
            EDFAC_EX = XLEN(GED)/AREAS(EL_EX)

C.....Younghun added:
C.....Compute H (water depth), Q in 2D KWE at each edge Gauss quadrature point
            
            DO I = 1,NEGP(pa)
            
               GP_IN = I
               GP_EX = NEGP(pa) - I + 1

               H_IN = H(1,EL_IN,IRK)
               H_EX = H(1,EL_EX,IRK)
               
               DO K = 2,DOFS(EL_IN)
               
                  H_IN = H_IN + H(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN,pa)
                  H_EX = H_EX + H(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX,pa)
                  
               ENDDO
               
               ! Check if water depth is below 0
               IF (H_IN .lt. 0.D0) THEN
                  H_IN = 0.D0
               ENDIF
               
               IF (H_EX .lt. 0.D0) THEN
                  H_EX = 0.D0
               ENDIF
               
               ! Compute flux
               HQX_IN = S0X(EL_IN)*(H_IN**(5.D0/3.D0)) / (manningKWE2(EL_IN)*sqrt(S0mag(EL_IN)))
               HQY_IN = S0Y(EL_IN)*(H_IN**(5.D0/3.D0)) / (manningKWE2(EL_IN)*sqrt(S0mag(EL_IN)))

               HQX_EX = S0X(EL_EX)*(H_EX**(5.D0/3.D0)) / (manningKWE2(EL_EX)*sqrt(S0mag(EL_EX)))
               HQY_EX = S0Y(EL_EX)*(H_EX**(5.D0/3.D0)) / (manningKWE2(EL_EX)*sqrt(S0mag(EL_EX)))

C.....Compute the numerical flux (upwind scheme)
                            
               HQ_HAT_IN = HQX_IN*NX + HQY_IN*NY
               HQ_HAT_EX = HQX_EX*NX + HQY_EX*NY
               
C.....Compute the edge integral
               
               DO K = 1,DOFS(EL)
                  
                  W_IN = EDFAC_IN*EDGEQ(K,GP_IN,LED_IN,pa)
                  W_EX = EDFAC_EX*EDGEQ(K,GP_EX,LED_EX,pa)
                  
                  ! Check edge type
                  IF (ET(GED).eq.-1) THEN ! flowing edge (from EL_IN to EL_EX)
                     RHS_KW2(K,EL_IN,IRK) = RHS_KW2(K,EL_IN,IRK) - W_IN*HQ_HAT_IN
                     RHS_KW2(K,EL_EX,IRK) = RHS_KW2(K,EL_EX,IRK) + W_EX*HQ_HAT_IN
                     
                  ELSEIF (ET(GED).eq.0) THEN ! diverging edge
                     ! No flow
                  ELSEIF (ET(GED).eq.-2) THEN ! converging edge (it should not happen after introducing channel_edge_runoff subroutine)
                     RHS_KW2(K,EL_IN,IRK) = RHS_KW2(K,EL_IN,IRK) - W_IN*HQ_HAT_IN
                     RHS_KW2(K,EL_EX,IRK) = RHS_KW2(K,EL_EX,IRK) + W_EX*HQ_HAT_EX
                  ENDIF

               ENDDO
        
            ENDDO
            
         ELSEIF (flipNEDEL(GED).eq.1) THEN
            LED_IN = NEDSD(2,GED)
            LED_EX = NEDSD(1,GED)
            
            EL_IN = NEDEL(2,GED)
            EL_EX = NEDEL(1,GED)
           
            EL = EL_EX
           
            if (DOFS(EL_EX).LT.DOFS(EL_IN)) then
               EL = EL_IN
            endif
            
            pa = PDG_EL(EL)

#ifdef P0
            if (pa.eq.0) then
               pa = 1
            endif
#endif

C.....Retrieve the components of the normal vector to the edge
        
            NX = -COSNX(GED)
            NY = -SINNX(GED)
         
C.....Set the components for the tangential vector to the edge

            EDFAC_IN = XLEN(GED)/AREAS(EL_IN)
            EDFAC_EX = XLEN(GED)/AREAS(EL_EX)

C.....Younghun added:
C.....Compute H (water depth), Q in 2D KWE at each edge Gauss quadrature point

            DO I = 1,NEGP(pa)
            
               GP_IN = NEGP(pa) - I + 1
               GP_EX = I
               
               H_IN = H(1,EL_IN,IRK)
               H_EX = H(1,EL_EX,IRK)
               
               DO K = 2,DOFS(EL_IN)
               
                  H_IN = H_IN + H(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN,pa)
                  H_EX = H_EX + H(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX,pa)
                  
               ENDDO
               
               ! Check if water depth is below 0
               IF (H_IN .lt. 0.D0) THEN
                  H_IN = 0.D0
               ENDIF
               
               IF (H_EX .lt. 0.D0) THEN
                  H_EX = 0.D0
               ENDIF
               
               ! Compute flux
               HQX_IN = S0X(EL_IN)*(H_IN**(5.D0/3.D0)) / (manningKWE2(EL_IN)*sqrt(S0mag(EL_IN)))
               HQY_IN = S0Y(EL_IN)*(H_IN**(5.D0/3.D0)) / (manningKWE2(EL_IN)*sqrt(S0mag(EL_IN)))
                                                  
               HQX_EX = S0X(EL_EX)*(H_EX**(5.D0/3.D0)) / (manningKWE2(EL_EX)*sqrt(S0mag(EL_EX)))
               HQY_EX = S0Y(EL_EX)*(H_EX**(5.D0/3.D0)) / (manningKWE2(EL_EX)*sqrt(S0mag(EL_EX)))

C.....Compute the numerical flux (upwind scheme)
                            
               HQ_HAT_IN = HQX_IN*NX + HQY_IN*NY
               HQ_HAT_EX = HQX_EX*NX + HQY_EX*NY
               
C.....Compute the edge integral

               DO K = 1,DOFS(EL)
                  
                  W_IN = EDFAC_IN*EDGEQ(K,GP_IN,LED_IN,pa)
                  W_EX = EDFAC_EX*EDGEQ(K,GP_EX,LED_EX,pa)
                  
                  IF (ET(GED).eq.-1) THEN
     
                     RHS_KW2(K,EL_IN,IRK) = RHS_KW2(K,EL_IN,IRK) - W_IN*HQ_HAT_IN
                     RHS_KW2(K,EL_EX,IRK) = RHS_KW2(K,EL_EX,IRK) + W_EX*HQ_HAT_IN
                  ELSEIF (ET(GED).eq.0) THEN ! diverging edge
                     ! No flow
                  ELSEIF (ET(GED).eq.-2) THEN ! converging edge (it should not happen after introducing channel_edge_runoff subroutine)
                     RHS_KW2(K,EL_IN,IRK) = RHS_KW2(K,EL_IN,IRK) - W_IN*HQ_HAT_IN
                     RHS_KW2(K,EL_EX,IRK) = RHS_KW2(K,EL_EX,IRK) + W_EX*HQ_HAT_EX
                  ENDIF
                  
               ENDDO
            
            ENDDO
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE
