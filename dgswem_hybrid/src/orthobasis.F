C***********************************************************************
C     
C     Subroutine to compute the values of the orthogonal basis at the
C     given point
C     
C***********************************************************************

      SUBROUTINE ORTHOBASIS( X, Y, P, DOF, PHI, DPHIDZ1, DPHIDZ2 )
      
      IMPLICIT NONE

      INTEGER A, B, I, J, K, P, DOF
      
      REAL(8) JACOBI(P+1,2*P+3,2)
      REAL(8) PHI2(P+1,P+1), DPHIDZ12(P+1,P+1), DPHIDZ22(P+1,P+1)
      REAL(8) PHI(DOF), DPHIDZ1(DOF), DPHIDZ2(DOF)
      REAL(8) X, Y, Z
      REAL(8) NUM1, DEN, DEN1, COEFF1, COEFF2, POLY1, POLY2

C.....Check to see if the point is one of the vertices of the element

      IF ( (X.EQ.-1).AND.(Y.EQ.-1) ) THEN
         DO I = 0,P
            DO J = 0,P-I
               PHI2(I+1,J+1) = (-1.D0)**(I+J+2)
            ENDDO
         ENDDO
         GOTO 111
      ENDIF

      IF ( (X.EQ.1).AND.(Y.EQ.-1) ) THEN
         DO I = 0,P
            DO J = 0,P-I
               PHI2(I+1,J+1) = (-1.D0)**(J+2)
            ENDDO
         ENDDO
         GOTO 111
      ENDIF

      IF ( (X.EQ.-1).AND.(Y.EQ.1) ) THEN
         DO I = 0,P
            DO J = 0,P-I
               IF ((I+1).EQ.1) THEN
                  PHI2(I+1,J+1) = J+1.D0
               ELSE
                  PHI2(I+1,J+1) = 0.D0
               ENDIF
            ENDDO
         ENDDO
         GOTO 111
      ENDIF

C.....If not construct and evaluate the necessary Jacobi polynomials at
C.....the given point

      DO A = 0,2*P+2
         DO B = 0,1
            DO J = 0,P

               IF ((A.EQ.0).OR.((A.EQ.1).AND.(B.EQ.1))) THEN
                  Z = 2.D0*(1.D0 + X)/(1.D0 - Y) - 1.D0
               ELSE
                  Z = Y
               ENDIF

               IF (J.EQ.0) THEN
                  JACOBI(J+1,A+1,B+1) = 1.D0
               ELSEIF (J.EQ.1) THEN
                  JACOBI(J+1,A+1,B+1) =(2*(A+1) + (A+B+2)*(Z-1.D0))/2.D0
               ELSEIF (J.GE.2) THEN

                  NUM1 = 1.D0
                  DO I = 1,2*(J-1)+A+B+2
                     NUM1 = NUM1*I
                  ENDDO

                  DEN1 = 1.D0
                  DO I=1,2*(J-1)+A+B-1
                     DEN1 = DEN1*I
                  ENDDO

                  COEFF1 = (2*(J-1)+A+B+1)*(A**2-B**2) + Z*NUM1/DEN1
                  COEFF2 = -2*(J-1+A)*(J-1+B)*(2*(J-1)+A+B+2)
                  POLY1  = JACOBI(J,A+1,B+1)
                  POLY2  = JACOBI(J-1,A+1,B+1)
                  DEN = 2*(J-1+1)*(J-1+A+B+1)*(2*(J-1)+A+B)

                  JACOBI(J+1,A+1,B+1) = (COEFF1*POLY1 + COEFF2*POLY2)/DEN

               ENDIF

            ENDDO
         ENDDO
      ENDDO

C.....Construct and evaluate the orthogonal basis and its derivatives at
C.....the given point

      DO I = 0,P
         DO J = 0,P-I
            PHI2(I+1,J+1) = JACOBI(I+1,1,1)*((1.D0-Y)/2.D0)**I
     &           *JACOBI(J+1,2*I+2,1)
            
            IF (I.EQ.0) THEN
               DPHIDZ12(I+1,J+1) = 0.D0
            ELSE
               DPHIDZ12(I+1,J+1) = 2.D0/(1.D0-Y)*(I+1.D0)/2.D0*
     &              JACOBI(I,2,2)*((1-Y)/2.D0)**I*JACOBI(J+1,2*I+2,1)
            ENDIF

            IF ((I.EQ.0).AND.(J.EQ.0)) THEN
               DPHIDZ22(I+1,J+1) = 0.D0
            ELSEIF ((I.EQ.0).AND.(J.GE.1)) THEN
               DPHIDZ22(I+1,J+1) = (J+2)/2.D0*JACOBI(J,3,2)
            ELSEIF ((J.EQ.0).AND.(I.GE.1)) THEN
               DPHIDZ22(I+1,J+1) = 2.D0*(1.D0+X)/(1.D0-Y)**2*(I+1)/2.D0*
     &              JACOBI(I,2,2)*((1.D0-Y)/2.D0)**I*
     &              JACOBI(J+1,2*I+2,1) + JACOBI(I+1,1,1)*(-I/2.D0*
     &              ((1.D0-Y)/2.D0)**(I-1))*JACOBI(J+1,2*I+2,1)
            ELSE
               DPHIDZ22(I+1,J+1) = 2.D0*(1.D0+X)/(1.D0-Y)**2*(I+1)/2.D0*
     &              JACOBI(I,2,2)*((1.D0-Y)/2.D0)**I*
     &              JACOBI(J+1,2*I+2,1) + JACOBI(I+1,1,1)*(-I/2.D0*
     &              ((1.D0-Y)/2.D0)**(I-1))*JACOBI(J+1,2*I+2,1) +
     &              JACOBI(I+1,1,1)*((1.D0-Y)/2.D0)**I*
     &              (J+2*I+2)/2.D0*JACOBI(J,2*I+3,2)
            ENDIF
         ENDDO
      ENDDO

C.....Re-order the basis functions into hierarchical order

 111  K = 1
      DO J = 0,P
         DO I = 0,J
            PHI(K) = PHI2(I+1,J-I+1)
            IF (ABS(PHI(K)).LT.1.0E-15) PHI(K) = 0.D0
            DPHIDZ1(K) = DPHIDZ12(I+1,J-I+1)
            DPHIDZ2(K) = DPHIDZ22(I+1,J-I+1)
            K = K+1
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE ORTHOBASIS
