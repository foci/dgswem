C-----------------------------------------------------------------------
C
C     Subroutine ORTHOGONAL_BASIS
C
C     Written by Ethan Kubatko (11-12-2009)
C
C     This subroutine evaluates the orthogonal basis functions and its
C     derivatives at a given point for a 2D triangle, a 3D
C     triangular prism, a 2D quadrilateral or a 3D
C     hexahedral.  The recurrence relations used in the 2D triangle
C     are from [1].
C
C-----------------------------------------------------------------------
C
C     Input:
C     ------
C       ELEM:    Element type (1 for Triangle 2 for Quad)
C       PT:      Coordinate point (of dimension 2 or 3)
C       D:       Polynomial degree of basis
C       DIM:     Dimension (2 or 3)
C
C     Output:
C     -------
C       BASIS:  Array of length (D+2)*(D+1)**(DIM-1)/2 containing the
C               values of the basis functions evaluated at PT.
C       DBASIS: Array of size (D+2)*(D+1)**(DIM-1)/2 by DIM containing
C               the values of the gradients of the basis functions
C               evaluated at PT.
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
C
C     Updates:
C
C     Quadrilateral and hexahedral element bases added
C       - Ashley Maggi, 07-13-2010
C
C-----------------------------------------------------------------------


      SUBROUTINE ORTHOGONAL_BASIS(ELEM, PT, D, DIM, BASIS, DBASIS )
      
      IMPLICIT NONE

C.....Declare subroutine input and output

      INTEGER, INTENT(IN)  :: D, DIM, ELEM
      REAL(8), INTENT(IN)  :: PT(DIM)
      REAL(8), INTENT(OUT) :: BASIS((D+2)*(D+1)**(DIM-1)/2)
      REAL(8), INTENT(OUT) :: DBASIS((D+2)*(D+1)**(DIM-1)/2,DIM)
      
C.....Explicitly declare remaining variables
      
      INTEGER :: P, Q, R
      REAL(8) :: AQ, BQ, CQ, X, Y, Z
      REAL(8), DIMENSION(0:D) :: LEGENDRE
      REAL(8), DIMENSION(0:D) :: LEGENDRE_X, LEGENDRE_Y, LEGENDRE_Z
      REAL(8), DIMENSION(0:D) :: DLEGENDRE
      REAL(8), DIMENSION(0:D) :: DLEGENDRE_X, DLEGENDRE_Y, DLEGENDRE_Z
      REAL(8), DIMENSION(0:D,0:D) :: PHI
      REAL(8), DIMENSION(0:D,0:D,2) :: DPHI
      REAL(8), DIMENSION((D+2)*(D+1)/2) :: PHI_2D
      REAL(8), DIMENSION((D+2)*(D+1)/2,2) :: DPHI_2D

C.....Check to make sure input is correct

      IF (D.LT.0) THEN
        PRINT*,'  ****** ERROR!!! D must be a positve integer *******  '
        PRINT*,'  Execution terminated in subroutine orthogonal_basis  '
        STOP
      ELSEIF ((DIM.GT.3).OR.(DIM.LT.2)) THEN
        PRINT*,'  *********** ERROR!!! DIM must be 2 or 3 ***********  '
        PRINT*,'  Execution terminated in subroutine orthogonal_basis  '
        STOP
      ELSEIF ((ELEM.NE.1).AND.(ELEM.NE.2)) THEN
       PRINT*,'  **** ERROR!!! ELEM must be TRIA = 1 or QUAD = 2 ****  '
       PRINT*,'  Execution terminated in subroutine orthogonal_basis  '
       STOP
      ENDIF

C.....Assign X, Y, and Z coordinates

      X = PT(1)
      Y = PT(2)
      IF (DIM.EQ.3) THEN
        Z = PT(3)
      ELSE
        Z = 0.
      ENDIF
      
C.....Case D = 0 for 2D or 3D

      IF (D.EQ.0) THEN
        BASIS  = 1.
        DBASIS = 0.
        RETURN
      ENDIF

C-----------------------------------------------------------------------
C
C     2D Triangular basis
C
C-----------------------------------------------------------------------
      IF (ELEM.EQ.1) THEN
      
C.......Check to make sure point PT is contained within the reference
C.......element.  Issue warning if not but continue execution.

        IF ((ABS(X).LE.1).AND.(ABS(Y).LE.1).AND.
     &                    ((X+Y).LE.0).AND.(ABS(Z).LE.1)) THEN
          CONTINUE
        ELSE
         PRINT*,' ***** WARNING!! in subroutine orthogonal_basis ***** '
         PRINT*,'       Point is outside of reference element!!        '
         PRINT*,'               Execution will continue                '
        ENDIF
      
C.......First basis function and derivatives

        PHI(0,0)    = 1.
        DPHI(0,0,1) = 0.
        DPHI(0,0,2) = 0.
      
C.......Second basis function and derivatives
      
        PHI(1,0)    = (1+2*X+Y)/2.
        DPHI(1,0,1) = 1.
        DPHI(1,0,2) = 1./2
      
C.......Recurrence relations
      
        DO P = 1,D-1
          PHI(P+1,0)    = (2.*P+1.)/(P+1.)*PHI(1,0)*PHI(P,0)
     &                  - P/(P+1.)*((1.-Y)/2.)**2.*PHI(P-1,0)
          DPHI(P+1,0,1) = (2.*P+1.)/(P+1.)*(DPHI(1,0,1)*PHI(P,0)
     &                  + PHI(1,0)*DPHI(P,0,1))
     &                  - P/(P+1.)*((1.-Y)/2.)**2.*DPHI(P-1,0,1)
          DPHI(P+1,0,2) = (2.*P+1.)/(P+1.)*(DPHI(1,0,2)*PHI(P,0)
     &                  + PHI(1,0)*DPHI(P,0,2)) - P/(P+1.)*((Y-1)/2
     &                  * PHI(P-1,0) + ((1.-Y)/2.)**2.*DPHI(P-1,0,2))
        ENDDO
        DO P = 0,D-1
          PHI(P,1)    =  PHI(P,0)*(1.+2.*P+(3.+2.*P)*Y)/2.
          DPHI(P,1,1) = DPHI(P,0,1)*(1.+2.*P+(3.+2.*P)*Y)/2.
          DPHI(P,1,2) = DPHI(P,0,2)*(1.+2.*P+(3.+2.*P)*Y)/2.
     &                + PHI(P,0)*(3.+2.*P)/2.
        ENDDO

        DO Q = 1,D-1
          DO P = 0,D-Q-1
            AQ = (2.*Q+2+2.*P)*(2.*Q+3+2.*P)/(2.*Q+2)/(Q+2+2.*P)
            BQ = (2.*P+1)**2.*(2.*Q+2+2.*P)/
     &                                 (2.*Q+2)/(2.*Q+2.*P+1)/(Q+2+2.*P)
            CQ = (Q+2.*P+1)*Q*(2.*Q+3+2.*P)/
     &                                   (Q+1.)/(Q+2+2.*P)/(2.*Q+2.*P+1)
            PHI(P,Q+1)    = (AQ*Y+BQ)*PHI(P,Q) - CQ*PHI(P,Q-1)
            DPHI(P,Q+1,1) = (AQ*Y+BQ)*DPHI(P,Q,1) - CQ*DPHI(P,Q-1,1)
            DPHI(P,Q+1,2) = (AQ*Y+BQ)*DPHI(P,Q,2) + AQ*PHI(P,Q)
     &                                            - CQ*DPHI(P,Q-1,2)
          ENDDO
        ENDDO
        
C-----------------------------------------------------------------------
C
C     2D Quadrilateral basis
C
C-----------------------------------------------------------------------
        
      ELSEIF (ELEM.EQ.2) THEN
      
C.......Check to make sure point PT is contained within the reference
C.......element.  Issue warning if not but continue execution.
      
        IF ((ABS(X).LE.1).AND.(ABS(Y).LE.1).AND.(ABS(Z).LE.1)) THEN
          CONTINUE
        ELSE
          PRINT*,'  ***** WARNING in suroutine orthogonal_basis *****  '
          PRINT*,'       Point is outside of reference element!!       '
          PRINT*,'               Execution will continue               '
        ENDIF

C.......Legendre Polynomials in X

        LEGENDRE_X(0)  = 1.
        DLEGENDRE_X(0) = 0.
        LEGENDRE_X(1)  = X
        DLEGENDRE_X(1) = 1.
        DO P = 1,D-1
         LEGENDRE_X(P+1) = ((2.*P+1.)*X*LEGENDRE_X(P)-P*LEGENDRE_X(P-1))
     &                                                           /(P+1.)
         DLEGENDRE_X(P+1) = ((2.*P+1.)*(X*DLEGENDRE_X(P) +
     &                       LEGENDRE_X(P)) - P*DLEGENDRE_X(P-1))/(P+1.)
        ENDDO

C.......Legendre Polynomials in Y

        LEGENDRE_Y(0)  = 1.
        DLEGENDRE_Y(0) = 0.
        LEGENDRE_Y(1)  = Y
        DLEGENDRE_Y(1) = 1.
        DO P = 1,D-1
         LEGENDRE_Y(P+1) = ((2.*P+1.)*Y*LEGENDRE_Y(P)-P*LEGENDRE_Y(P-1))
     &                                                           /(P+1.)
         DLEGENDRE_Y(P+1) = ((2.*P+1.)*(Y*DLEGENDRE_Y(P) +
     &                       LEGENDRE_Y(P)) - P*DLEGENDRE_Y(P-1))/(P+1.)
        ENDDO

C.......Forming the basis functions
        DO P = 0,D
          DO Q = 0,D-P
            PHI(Q,P) =  LEGENDRE_X(Q)*LEGENDRE_Y(P)
C.......DPHI with respect to X
            DPHI(Q,P,1) = DLEGENDRE_X(Q)*LEGENDRE_Y(P)
C.......DPHI with respect to Y
            DPHI(Q,P,2) = LEGENDRE_X(Q)*DLEGENDRE_Y(P)
          ENDDO
        ENDDO
      ENDIF
      
C.....Re-order basis functions hierarchically

      R = 1
      DO Q = 0,D
        DO P = 0,Q
          PHI_2D(R)    = PHI(P,Q-P)
          DPHI_2D(R,1) = DPHI(P,Q-P,1)
          DPHI_2D(R,2) = DPHI(P,Q-P,2)
          IF (ABS(PHI_2D(R)).LT.1.0E-15) PHI_2D(R) = 0.
          IF (ABS(DPHI_2D(R,1)).LT.1.0E-15) DPHI_2D(R,1) = 0.
          IF (ABS(DPHI_2D(R,2)).LT.1.0E-15) DPHI_2D(R,2) = 0.
          R = R + 1
        ENDDO
      ENDDO
      

      IF (DIM.EQ.2) THEN
         BASIS  = PHI_2D
         DBASIS = DPHI_2D
      ELSEIF (DIM.EQ.3) THEN

C-----------------------------------------------------------------------
C
C     3D Bases
C
C-----------------------------------------------------------------------

         LEGENDRE(0)  = 1.
         DLEGENDRE(0) = 0.
         LEGENDRE(1)  = Z
         DLEGENDRE(1) = 1.
         DO P = 1,D-1
            LEGENDRE(P+1) =((2*P+1)*Z*LEGENDRE(P)-P*LEGENDRE(P-1))/(P+1)
            DLEGENDRE(P+1) = ((2*P+1)*(Z*DLEGENDRE(P) + LEGENDRE(P))
     &                                         - P*DLEGENDRE(P-1))/(P+1)
         ENDDO
         R = 1
         DO Q = 0,D
            DO P = 1,(D+2)*(D+1)/2
               BASIS(R)    = PHI_2D(P)*LEGENDRE(Q)
               DBASIS(R,1) = DPHI_2D(P,1)*LEGENDRE(Q)
               DBASIS(R,2) = DPHI_2D(P,2)*LEGENDRE(Q)
               DBASIS(R,3) = PHI_2D(P)*DLEGENDRE(Q)
               IF (ABS(BASIS(R)).LT.1.0E-15) BASIS(R) = 0.
               IF (ABS(DBASIS(R,1)).LT.1.0E-15) DBASIS(R,1) = 0.
               IF (ABS(DBASIS(R,2)).LT.1.0E-15) DBASIS(R,2) = 0.
               IF (ABS(DBASIS(R,3)).LT.1.0E-15) DBASIS(R,3) = 0.
               R = R + 1
            ENDDO
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE ORTHOGONAL_BASIS
