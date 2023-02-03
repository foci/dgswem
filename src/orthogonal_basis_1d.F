C-----------------------------------------------------------------------
C
C     Subroutine ORTHOGONAL_BASIS_1D
C
C     Written by Younghun Kang (01-28-2020)
C
C     This subroutine evaluates the orthogonal basis functions on 1D.
C     The recurrence relations are from [1].
C
C-----------------------------------------------------------------------
C
C     Input:
C     ------
C       PT:      Coordinate point (of dimension 1)
C       D:       Polynomial degree of basis
C
C     Output:
C     -------
C       BASIS:  Array of length D containing the values of the basis
C               functions evaluated at PT.
C       DBASIS: Array of length D containing the values of the gradients
C               of the basis functions evaluated at PT.
C
C-----------------------------------------------------------------------
C
C     Reference: [1] Robert C. Kirby, "Singularity-free evaluation of
C                    collapsed-coordinate orthogonal polynomials", ACM
C                    Transactions on Mathematical Software, Volume 37
C                    Issue 1, January 2010, Article No. 5, Pages 5:1--
C                    5:16.
C
C                    http://doi.acm.org/10.1145/1644001.1644006
C
C-----------------------------------------------------------------------


      SUBROUTINE ORTHOGONAL_BASIS_1D(PT, D, BASIS, DBASIS )
      
      IMPLICIT NONE

C.....Declare subroutine input and output

      INTEGER, INTENT(IN)  :: D
      REAL(8), INTENT(IN)  :: PT
      
      REAL(8), INTENT(OUT) :: BASIS(D+1)
      REAL(8), INTENT(OUT) :: DBASIS(D+1)

      
C.....Explicitly declare remaining variables
      
      INTEGER :: P
      REAL(8) :: X
      REAL(8), DIMENSION(0:D) :: P00, P11, DP00

C.....Check to make sure input is correct

      IF (D.LT.0) THEN
        PRINT*,'  ****** ERROR!!! D must be a positve integer *******  '
        PRINT*,'  Execution terminated in subroutine orthogonal_basis  '
        STOP
      ! ELSEIF ((DIM.GT.1).OR.(DIM.LT.1)) THEN
        ! PRINT*,'  ********* ERROR!!! DIM must be 1 *********  '
        ! PRINT*,'  Execution terminated in subroutine orthogonal_basis  '
        ! STOP
      ! ELSEIF ((ELEM.NE.1).AND.(ELEM.NE.2)) THEN
       ! PRINT*,'  **** ERROR!!! ELEM must be TRIA = 1 or QUAD = 2 ****  '
       ! PRINT*,'  Execution terminated in subroutine orthogonal_basis  '
       ! STOP
      ENDIF

C.....Assign X, Y, and Z coordinates

      X = PT
      ! IF (DIM.EQ.3) THEN
        ! Z = PT(3)
      ! ELSE
        ! Z = 0.D0
      ! ENDIF
      
C.....Case D = 0 for 2D or 3D

      IF (D.EQ.0) THEN
        BASIS  = 1.D0
        DBASIS = 0.D0
        RETURN
      ENDIF
      
C-----------------------------------------------------------------------
C
C     1D Bases
C
C-----------------------------------------------------------------------      
      
      P00(0) = 1.D0
      P00(1) = X
      P11(0) = 1.D0
      P11(1) = 2.d0*X
      DP00(0) = 0.D0
      DP00(1) = 1.d0
      
      DO P = 1, D-1
         
         ! write(*,*) p, (P+1.D0)/2.D0
         P00(P+1) = (2.d0*P+1.d0)/(P+1.d0) * X * P00(P) 
     &                           - (P*1.d0)/(P+1.d0)*P00(P-1)
      ENDDO
      
      DO P = 1,D-1
         DP00(P+1) = (2.d0*P+1.d0)/(P+1.d0) * P00(P) 
     &             + (2.d0*P+1.d0)/(P+1.d0) * X * DP00(P) 
     &             - (P*1.d0)/(P+1.d0)*DP00(P-1)
      ENDDO
      
      BASIS = P00
      DBASIS = DP00
      
      

      RETURN
      END SUBROUTINE ORTHOGONAL_BASIS_1D
