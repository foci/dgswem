C***********************************************************************
C
C     SUBROUTINE CALC_TRISLOPE( )
C
C     This subroutine calculates the slope vector for a given element
C
C     Written by Younghun Kang (01-28-2020)
C
C***********************************************************************

      SUBROUTINE CALC_TRISLOPE()

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      
      IMPLICIT NONE
      
C.....Declare local variables
      INTEGER I
      REAL(SZ) Z1, Z2, Z3
      REAL(SZ) X21, X31, Y21, Y31, Z21, Z31
      REAL(SZ) II1, II2, II3, JJ1, JJ2, JJ3

      ! INTEGER IEL, IED,i

C.....Loop over the elements

      DO I = 1,NE
      
C.....Retrieve the node coefficients for the given element
        
        ! X1 = X(NM(I,1))
        ! X2 = X(NM(I,2))
        ! X3 = X(NM(I,3))
        
        ! Y1 = Y(NM(I,1))
        ! Y2 = Y(NM(I,2))
        ! Y3 = Y(NM(I,3))
        
        Z1 = DP(NM(I,1))
        Z2 = DP(NM(I,2))
        Z3 = DP(NM(I,3))
        
C.....Compute the slope
        
        X21 = X(NM(I,2)) - X(NM(I,1))
        X31 = X(NM(I,3)) - X(NM(I,1))

        Y21 = Y(NM(I,2)) - Y(NM(I,1))
        Y31 = Y(NM(I,3)) - Y(NM(I,1))
        
        ! minus sign comes from bathymetry < 0 on dry nodes
        Z21 = -(DP(NM(I,2)) - DP(NM(I,1)))
        Z31 = -(DP(NM(I,3)) - DP(NM(I,1)))
        
        S0X(I) = (Z31*Y21 - Z21*Y31)/(X21*Y31 - X31*Y21)
        S0Y(I) = (Z21*X31 - Z31*X21)/(X21*Y31 - X31*Y21)
        
        S0mag(I) = sqrt(S0X(I)**2.d0 + S0Y(I)**2.d0)
        
      ENDDO
        
      RETURN
      END


