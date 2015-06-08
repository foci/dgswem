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
!     
!***********************************************************************

                                !Is this even used? -- CM

      SUBROUTINE OCEAN_FLOW_EDGE_HYDRO(L,QBC)

!.....Use appropriate modules

      USE GLOBAL
      USE DG

      IMPLICIT NONE

!.....Declare local variables

      INTEGER L, LED, GED,k,i
      REAL(SZ) TX, TY, QBC
      
!.....Retrieve the global and local edge number

      GED = NEEDN(L)
      LED = NEDSD(1,GED)

!.....Retrieve the element to which the edge belongs

      EL_IN = NEDEL(1,GED)
      PA = PDG_EL(EL_IN)

#ifdef P0
      if (pa.eq.0) then
         pa = 1
      endif
#endif

!.....Retrieve the components of the normal vector to the edge
      
      NX = COSNX(GED)
      NY = SINNX(GED)
      
!.....Set the components for the tangential vector to the edge

      TX = -NY
      TY =  NX
      
!.....Compute ZE, QX, QY, and HB at each edge Gauss quadrature point

      DO I = 1,NEGP(pa)

         ZE_IN = 0.D0
         QX_IN = 0.D0
         QY_IN = 0.D0

         ZE_EX = 0.D0
         QX_EX = 0.D0
         QY_EX = 0.D0
         
         HB_IN = 0.D0

#ifdef TRACE
         iota_IN = 0.D0
         iota_EX = 0.D0
#endif

#ifdef CHEM
         iota_IN = 0.D0
         iota_EX = 0.D0
         iota2_IN = 0.D0
         iota2_EX = 0.D0
#endif

#ifdef DYNP
         dynP_IN = 0.D0
         dynP_EX = 0.D0
#endif
         
!.....Compute the specified flow boundaries for the exterior state

         Q_N_EXT = QBC
         Q_T_EXT = 0.D0

         QX_EX = -( TY*Q_N_EXT - NY*Q_T_EXT)/(NX*TY - NY*TX)
         QY_EX = -(-TX*Q_N_EXT + NX*Q_T_EXT)/(NX*TY - NY*TX)

!.....Compute the solution at the interior state

         DO K = 1,DOFS(EL_IN)

            ZE_IN = ZE_IN + ZE(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa)
            QX_IN = QX_IN + QX(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa)
            QY_IN = QY_IN + QY(K,EL_IN,IRK)*PHI_EDGE(K,I,LED,pa)
            
            HB_IN = HB_IN + HB(K,EL_IN,1  )*PHI_EDGE(K,I,LED,pa)

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

         ENDDO
         
!.....Set the exterior bed and surface elevation equal to the interior

         ZE_EX = ZE_IN
         HB_EX = HB_IN
         
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
         
         
!.....Compute the Roe flux

         CALL NUMERICAL_FLUX()
         
!.....Compute the edge integral

         DO K = 1,DOFS(EL_IN)
            CALL EDGE_INT_HYDRO(EL_IN, LED, GED, I, F_HAT, G_HAT, H_HAT, i_hat, j_hat, k)
         ENDDO
         
      ENDDO
      
      RETURN
      END SUBROUTINE
