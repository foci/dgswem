C***********************************************************************
C     
C     SUBROUTINE EDGE_INT_HYDRO()
C     
C     This subroutine computes the edge integrals for the DG hydro using
C     Gauss quadrature and adds them to the RHS.
C     
C     Written by Ethan Kubatko (03-05-2005)
C     
C     02-23-2007 - sb - rewritten for better performance
C     01-10-2011 - cem - adapted for p_enrichment and multicomponent
C     
C***********************************************************************

#if 0
      SUBROUTINE EDGE_INT_HYDRO(EL_I,LED,GED,GP,FLUX1,FLUX2,FLUX3,flux4,flux5,flux6,k,pah)
                                ! <ezpp-noinst>
      
C.....Use appropriate modules

      USE GLOBAL
      USE DG

      IMPLICIT NONE
      
C.....Declare local variables

      INTEGER EL_I, LED, GED, GP, K, pah
      REAL(SZ) FLUX1, FLUX2, FLUX3, flux4, flux5, flux6, W

C.....Retrieve the elements which share the edge

#ifdef P0
      if (pah.eq.0) then
         pah = 1
      endif
#endif
      
C.....Comput the edge integral

      W = 2.0*M_INV(K,pah)/AREAS(EL_I)*XLEN(GED)*PHI_EDGE(K,GP,LED,pah)
     &     *WEGP(GP,pah)

      RHS_ZE(K,EL_I,IRK) = RHS_ZE(K,EL_I,IRK) - W*FLUX1
      RHS_QX(K,EL_I,IRK) = RHS_QX(K,EL_I,IRK) - W*FLUX2
      RHS_QY(K,EL_I,IRK) = RHS_QY(K,EL_I,IRK) - W*FLUX3

#ifdef TRACE
      RHS_iota(k,EL_I,IRK) = RHS_iota(k,EL_I,IRK) - W*FLUX4
#endif

#ifdef CHEM
      RHS_iota(k,EL_I,IRK) = RHS_iota(k,EL_I,IRK) - W*FLUX4
      RHS_iota2(k,EL_I,IRK) = RHS_iota2(k,EL_I,IRK) - W*FLUX5
#endif

#ifdef DYNP
      RHS_dynP(k,EL_I,IRK) = RHS_dynP(k,EL_I,IRK) - W*FLUX6
#endif
      
      RETURN
      END SUBROUTINE
#else
      SUBROUTINE EDGE_INT_HYDRO(EL_I,LED,GED,GP,FLUX1,FLUX2,FLUX3,flux4,flux5,flux6,k,pah)
                                ! <ezpp-noinst>
      
C.....Use appropriate modules

      USE GLOBAL
      USE DG

      IMPLICIT NONE
      
C.....Declare local variables

      INTEGER EL_I, LED, GED, GP,k,pah
      REAL(SZ) AREA, IMASS, FLUX1, FLUX2, FLUX3,flux4,flux5,flux6, W

C.....Retrieve the elements which share the edge

#ifdef P0
      if (pah.eq.0) then
         pah = 1
      endif
#endif

C.....Retrieve the element area
      
      AREA = 0.5D0*AREAS(EL_I)

C.....Comput the edge integral

      IMASS = M_INV(K,pah)/(0.5D0*AREA)

      RHS_ZE(K,EL_I,IRK) = RHS_ZE(K,EL_I,IRK)
     &     - IMASS*XLEN(GED)/2.D0*FLUX1*PHI_EDGE(K,GP,LED,pah)*WEGP(GP,pah)
      RHS_QX(K,EL_I,IRK) = RHS_QX(K,EL_I,IRK)
     &     - IMASS*XLEN(GED)/2.D0*FLUX2*PHI_EDGE(K,GP,LED,pah)*WEGP(GP,pah)
      RHS_QY(K,EL_I,IRK) = RHS_QY(K,EL_I,IRK)
     &     - IMASS*XLEN(GED)/2.D0*FLUX3*PHI_EDGE(K,GP,LED,pah)*WEGP(GP,pah)

#ifdef TRACE
      RHS_iota(k,EL_I,IRK) = RHS_iota(k,EL_I,IRK)
     &     - IMASS*XLEN(GED)/2.D0*FLUX4*PHI_EDGE(k,GP,LED,pah)*WEGP(GP,pah)
#endif

#ifdef CHEM
      RHS_iota(k,EL_I,IRK) = RHS_iota(k,EL_I,IRK)
     &     - IMASS*XLEN(GED)/2.D0*FLUX4*PHI_EDGE(k,GP,LED,pah)*WEGP(GP,pah)
      RHS_iota2(k,EL_I,IRK) = RHS_iota2(k,EL_I,IRK)
     &     - IMASS*XLEN(GED)/2.D0*FLUX5*PHI_EDGE(k,GP,LED,pah)*WEGP(GP,pah)
#endif

#ifdef DYNP
      RHS_dynP(k,EL_I,IRK) = RHS_dynP(k,EL_I,IRK)
     &     - IMASS*XLEN(GED)/2.D0*FLUX6*PHI_EDGE(k,GP,LED,pah)*WEGP(GP,pah)
#endif
      
      RETURN
      END SUBROUTINE
#endif
