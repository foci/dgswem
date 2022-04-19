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
      REAL(SZ) HWD_IN, HWD_EX, HQX_IN, HQY_IN, HQX_EX, HQY_EX, HQ_HAT_IN, HQ_HAT_EX
c     

      test_el = 0
      
      DO L = 1,NIEDS

C.......Retrieve the global and local edge number

         GED = NIEDN(L)

C.......Retrieve the elements which share the edge
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

               HWD_IN = HWD(1,EL_IN,IRK)
               HWD_EX = HWD(1,EL_EX,IRK)
               
               DO K = 2,DOFS(EL_IN)
               
                  HWD_IN = HWD_IN + HWD(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN,pa)
                  HWD_EX = HWD_EX + HWD(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX,pa)
                  
               ENDDO
               
               IF (HWD_IN .lt. 0.D0) THEN
                  HWD_IN = 0.D0
                  ! DO K = 1,DOFS(EL_IN)
                     ! HWD(K,EL_IN,IRK) = 0.D0
                  ! ENDDO
               ENDIF
               
               IF (HWD_EX .lt. 0.D0) THEN
                  HWD_EX = 0.D0
                  ! DO K = 1,DOFS(EL_IN)
                     ! HWD(K,EL_IN,IRK) = 0.D0
                  ! ENDDO
               ENDIF
               
               ! Manning_n = 0.015
               HQX_IN = S0X(EL_IN)*(HWD_IN**(5.D0/3.D0)) / (manningKWE2(EL_IN)*sqrt(S0mag(EL_IN)))
               HQY_IN = S0Y(EL_IN)*(HWD_IN**(5.D0/3.D0)) / (manningKWE2(EL_IN)*sqrt(S0mag(EL_IN)))

               HQX_EX = S0X(EL_EX)*(HWD_EX**(5.D0/3.D0)) / (manningKWE2(EL_EX)*sqrt(S0mag(EL_EX)))
               HQY_EX = S0Y(EL_EX)*(HWD_EX**(5.D0/3.D0)) / (manningKWE2(EL_EX)*sqrt(S0mag(EL_EX)))
               
               ! print*,S0X(EL_EX) / (manningKWE2(EL_EX)*sqrt(S0mag(EL_EX)))

C.....Compute the numerical flux (upwind scheme)
                            
               HQ_HAT_IN = HQX_IN*NX + HQY_IN*NY
               HQ_HAT_EX = HQX_EX*NX + HQY_EX*NY
               
               ! WRITE(*,*) HQ_HAT_IN,IT*DTDP
C.....Compute the edge integral
               
               DO K = 1,DOFS(EL)
                  
                  W_IN = EDFAC_IN*EDGEQ(K,GP_IN,LED_IN,pa)
                  W_EX = EDFAC_EX*EDGEQ(K,GP_EX,LED_EX,pa)
                  
C                  W_IN = 2.0d0*XLEN(GED)/AREAS(EL_IN)*M_INV(K,pa)*
C     &                     PHI_EDGE(K,GP_IN,LED_IN,pa)*WEGP(GP_IN,pa)
C                  W_EX = 2.0d0*XLEN(GED)/AREAS(EL_EX)*M_INV(K,pa)*
C     &                     PHI_EDGE(K,GP_EX,LED_EX,pa)*WEGP(GP_EX,pa)
     
                  IF (ET(GED).eq.-1) THEN
                     W_IN = 2.0d0*XLEN(GED)/AREAS(EL_IN)*M_INV(K,pa)*
     &                        PHI_EDGE(K,GP_IN,LED_IN,pa)*WEGP(GP_IN,pa)
                     W_EX = 2.0d0*XLEN(GED)/AREAS(EL_IN)*M_INV(K,pa)*
     &                        PHI_EDGE(K,GP_EX,LED_EX,pa)*WEGP(GP_EX,pa)
                  
                     RHS_HWD(K,EL_IN,IRK) = RHS_HWD(K,EL_IN,IRK) - W_IN*HQ_HAT_IN
                     ! RHS_HWD(K,EL_EX,IRK) = RHS_HWD(K,EL_EX,IRK) + W_EX*HQ_HAT_IN
                     IF (NEDCEDN(GED) .gt. 0) THEN
                        W_EX = EDFAC_EX*EDGEQ(K,GP_EX,LED_EX,pa)
                        ! RHS_HWD(K,EL_EX,IRK) = RHS_HWD(K,EL_EX,IRK) + W_EX*HQ_HAT_EX
                        ! stop
                     else
                        RHS_HWD(K,EL_EX,IRK) = RHS_HWD(K,EL_EX,IRK) + W_EX*HQ_HAT_IN
                     ENDIF
                     
                  ELSEIF (ET(GED).eq.0) THEN
                     ! No flow
                  ELSEIF (ET(GED).eq.-2) THEN
                     RHS_HWD(K,EL_IN,IRK) = RHS_HWD(K,EL_IN,IRK) - W_IN*HQ_HAT_IN
                     RHS_HWD(K,EL_EX,IRK) = RHS_HWD(K,EL_EX,IRK) + W_EX*HQ_HAT_EX
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

C.......If both elements on either side of edge are dry then skip

         ! IF((WDFLG(EL_IN).EQ.0).AND.(WDFLG(EL_EX).EQ.0)) GOTO 1000

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
               
               HWD_IN = HWD(1,EL_IN,IRK)
               HWD_EX = HWD(1,EL_EX,IRK)
               
               DO K = 2,DOFS(EL_IN)
               
                  HWD_IN = HWD_IN + HWD(K,EL_IN,IRK)*PHI_EDGE(K,GP_IN,LED_IN,pa)
                  HWD_EX = HWD_EX + HWD(K,EL_EX,IRK)*PHI_EDGE(K,GP_EX,LED_EX,pa)
                  
               ENDDO
               
               IF (HWD_IN .lt. 0.D0) THEN
                  HWD_IN = 0.D0
                  ! DO K = 1,DOFS(EL_IN)
                     ! HWD(K,EL_IN,IRK) = 0.D0
                  ! ENDDO
               ENDIF
               
               IF (HWD_EX .lt. 0.D0) THEN
                  HWD_EX = 0.D0
                  ! DO K = 1,DOFS(EL_IN)
                     ! HWD(K,EL_IN,IRK) = 0.D0
                  ! ENDDO
               ENDIF
               
               
               ! Manning_n = 0.015
               HQX_IN = S0X(EL_IN)*(HWD_IN**(5.D0/3.D0)) / (manningKWE2(EL_IN)*sqrt(S0mag(EL_IN)))
               HQY_IN = S0Y(EL_IN)*(HWD_IN**(5.D0/3.D0)) / (manningKWE2(EL_IN)*sqrt(S0mag(EL_IN)))
                                                  
               HQX_EX = S0X(EL_EX)*(HWD_EX**(5.D0/3.D0)) / (manningKWE2(EL_EX)*sqrt(S0mag(EL_EX)))
               HQY_EX = S0Y(EL_EX)*(HWD_EX**(5.D0/3.D0)) / (manningKWE2(EL_EX)*sqrt(S0mag(EL_EX)))


C.....Compute the numerical flux (upwind scheme)
                            
               HQ_HAT_IN = HQX_IN*NX + HQY_IN*NY
               HQ_HAT_EX = HQX_EX*NX + HQY_EX*NY
               
C.....Compute the edge integral

               DO K = 1,DOFS(EL)
                  
                  W_IN = EDFAC_IN*EDGEQ(K,GP_IN,LED_IN,pa)
                  W_EX = EDFAC_EX*EDGEQ(K,GP_EX,LED_EX,pa)
                  
                  IF (ET(GED).eq.-1) THEN
                     W_IN = 2.0d0*XLEN(GED)/AREAS(EL_IN)*M_INV(K,pa)*
     &                        PHI_EDGE(K,GP_IN,LED_IN,pa)*WEGP(GP_IN,pa)
                     W_EX = 2.0d0*XLEN(GED)/AREAS(EL_IN)*M_INV(K,pa)*
     &                        PHI_EDGE(K,GP_EX,LED_EX,pa)*WEGP(GP_EX,pa)
     
                     RHS_HWD(K,EL_IN,IRK) = RHS_HWD(K,EL_IN,IRK) - W_IN*HQ_HAT_IN
                     IF (NEDCEDN(GED) .gt. 0) THEN
                        W_EX = EDFAC_EX*EDGEQ(K,GP_EX,LED_EX,pa)
                        ! RHS_HWD(K,EL_EX,IRK) = RHS_HWD(K,EL_EX,IRK) + W_EX*HQ_HAT_EX
                        ! stop
                     else
                        RHS_HWD(K,EL_EX,IRK) = RHS_HWD(K,EL_EX,IRK) + W_EX*HQ_HAT_IN
                     ENDIF
                  ELSEIF (ET(GED).eq.0) THEN
                     ! No flow
                  ELSEIF (ET(GED).eq.-2) THEN
                     RHS_HWD(K,EL_IN,IRK) = RHS_HWD(K,EL_IN,IRK) - W_IN*HQ_HAT_IN
                     RHS_HWD(K,EL_EX,IRK) = RHS_HWD(K,EL_EX,IRK) + W_EX*HQ_HAT_EX
                  ENDIF
                  
               ENDDO
            
            ENDDO
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE
