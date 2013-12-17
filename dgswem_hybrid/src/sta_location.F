C***********************************************************************
C
C     SUBROUTINE STA_LOCATION()
C
C     This subroutine locates which element a station is located within
C     and computes the basis functions at that location.
C
C     Written by Ethan Kubatko (02-22-2007)
C
C***********************************************************************

      SUBROUTINE STA_LOCATION( XSTA, YSTA, PHI_STA, ELSTA )
      
      USE GLOBAL, ONLY : X, Y, NE, NM, AREAS
      USE DG, ONLY : DOF, PDG_EL
      USE SIZES

      IMPLICIT NONE

      INTEGER XL(NE), XR(NE), YB(NE), YT(NE)
      INTEGER J, ELSTA, NSTA, N1, N2, N3
      REAL(SZ) PHI_STA(DOF), DPHIDZ1(DOF), DPHIDZ2(DOF)
      REAL(SZ) XSTA, YSTA
      REAL(SZ) Z1, Z2, XP, YP, X1, X2, X3, Y1, Y2, Y3
      
C.....Determine the boundary of the rectangle which enclose the element

      DO J = 1,NE
      
        N1 = NM(J,1)
        N2 = NM(J,2)
        N3 = NM(J,3)

        YT(J) = MAX( Y(N1), Y(N2), Y(N3) )
        YB(J) = MIN( Y(N1), Y(N2), Y(N3) )

        XR(J) = MAX( X(N1), X(N2), X(N3) )
        XL(J) = MIN( X(N1), X(N2), X(N3) )
        
      ENDDO
      
C.....Loop over stations

      XP = XSTA
      YP = YSTA
      ELSTA = 0.D0
      

      DO 200 J = 1,NE

C.....Deterime if the points falls within the bounding box of element

        IF ( (XP.LE.XR(J)).AND.(XP.GE.XL(J)).AND.
     &       (YP.LE.YT(J)).AND.(YP.GE.YB(J)) ) THEN

C.....Retrieve the vertices of the element

          X1 = X(NM(J,1))
          X2 = X(NM(J,2))
          X3 = X(NM(J,3))

          Y1 = Y(NM(J,1))
          Y2 = Y(NM(J,2))
          Y3 = Y(NM(J,3))

C.....Transform to local Z element coordinates

          Z1 = -1.D0/(AREAS(J))*( XP*(2.D0*Y1 - 2.D0*Y3) +
     &                            YP*(2.D0*X3 - 2.D0*X1) +
     &                   X1*Y2 + X1*Y3 - X2*Y1 + X2*Y3 - X3*Y1 - X3*Y2 )

          Z2 =  1.D0/(AREAS(J))*( XP*(2.D0*Y1 - 2.D0*Y2) +
     &                            YP*(2.D0*X2 - 2.D0*X1) +
     &                   X1*Y2 + X1*Y3 - X2*Y1 - X2*Y3 - X3*Y1 + X3*Y2 )

C.....Determine if the node falls within the given element

          IF ( (Z1.GE.(-1)).AND.(Z2.GE.(-1)).AND.
     &        ((Z1 + Z2).LE.0) ) THEN
     
C.............Store the element in the appropriate array

            ELSTA = J
              
C.............Compute the basis functions at that point and store
              
            CALL ORTHOBASIS( Z1, Z2, pdg_el(j), DOF, PHI_STA, DPHIDZ1, DPHIDZ2 )
            GOTO 100
          ENDIF
        ENDIF
200   CONTINUE
100   CONTINUE

      RETURN
      END SUBROUTINE STA_LOCATION
      
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
     &                                            *JACOBI(J+1,2*I+2,1)

          IF (I.EQ.0) THEN
            DPHIDZ12(I+1,J+1) = 0.D0
          ELSE
            DPHIDZ12(I+1,J+1) = 2.D0/(1.D0-Y)*(I+1.D0)/2.D0*
     &                JACOBI(I,2,2)*((1-Y)/2.D0)**I*JACOBI(J+1,2*I+2,1)
          ENDIF

          IF ((I.EQ.0).AND.(J.EQ.0)) THEN
            DPHIDZ22(I+1,J+1) = 0.D0
          ELSEIF ((I.EQ.0).AND.(J.GE.1)) THEN
            DPHIDZ22(I+1,J+1) = (J+2)/2.D0*JACOBI(J,3,2)
          ELSEIF ((J.EQ.0).AND.(I.GE.1)) THEN
            DPHIDZ22(I+1,J+1) = 2.D0*(1.D0+X)/(1.D0-Y)**2*(I+1)/2.D0*
     &                   JACOBI(I,2,2)*((1.D0-Y)/2.D0)**I*
     &                   JACOBI(J+1,2*I+2,1) + JACOBI(I+1,1,1)*(-I/2.D0*
     &                   ((1.D0-Y)/2.D0)**(I-1))*JACOBI(J+1,2*I+2,1)
          ELSE
            DPHIDZ22(I+1,J+1) = 2.D0*(1.D0+X)/(1.D0-Y)**2*(I+1)/2.D0*
     &                   JACOBI(I,2,2)*((1.D0-Y)/2.D0)**I*
     &                   JACOBI(J+1,2*I+2,1) + JACOBI(I+1,1,1)*(-I/2.D0*
     &                   ((1.D0-Y)/2.D0)**(I-1))*JACOBI(J+1,2*I+2,1) +
     &                   JACOBI(I+1,1,1)*((1.D0-Y)/2.D0)**I*
     &                   (J+2*I+2)/2.D0*JACOBI(J,2*I+3,2)
          ENDIF
        ENDDO
      ENDDO

C.....Re-order the basis functions into hierarchical order

111   K = 1
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
