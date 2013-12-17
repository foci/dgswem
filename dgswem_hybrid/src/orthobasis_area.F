C***********************************************************************
C     
C     SUBROUTINE ORTHOBASIS_AREA()
C     
C     Evaluates the polynomials for the P degree orthogonal basis
C     for a triangle at the area integral gauss points.
C     
C     Reference for orthogonal basis:
C     
C     "A Discontinuous Galerkin Method for the Viscous MHD Equations",
C     Journal of Computational Physics 152, 608-641 (1999)
C     
C     Written by Ethan Kubatko (06-08-2004)
C     01-10-2011 - cem - adapted for p_enrichment
C     
C***********************************************************************

      SUBROUTINE ORTHOBASIS_AREA(P,L)

C.....Use appropriate modules

      USE GLOBAL
      USE DG

      IMPLICIT NONE

C.....Declare local variables

      INTEGER Q,A,B,M,k,i,jj,j,P,L
      REAL(SZ) COEFF1,COEFF2,POLY1,POLY2,NUM,DEN,NUM1,DEN1,DEN2,DEN3,Z
      REAL(SZ) XP,YP

C.....Allocate JACOBI

      if (L.eq.1.and.p.eq.1) then

         CALL ALLOC_JACOBI()

      endif

C.....Construct and evaluate the necessary Jacobi polynomials at the
C.....area integral gauss points and the barycenter of the element
      
      DO Q=1,NAGP(L)+1
         DO A=0,2*P+2
            DO B=0,1
               DO J=0,P
                  
                  IF (Q.EQ.NAGP(L)+1) THEN
                     XP = -1.D0/3.D0
                     YP = -1.D0/3.D0
                  ELSE
                     XP = XAGP(Q,L)
                     YP = YAGP(Q,L)
                  ENDIF
                  
                  IF ((A.EQ.0).OR.((A.EQ.1).AND.(B.EQ.1))) THEN
                     Z = 2.D0*(1.D0 + XP)/(1.D0 - YP) - 1.D0
                  ELSE
                     Z = YP
                  ENDIF

                  IF (J.EQ.0) THEN
                     JACOBI(J+1,A+1,B+1,Q) = 1.D0
                  ELSEIF (J.EQ.1) THEN
                     JACOBI(J+1,A+1,B+1,Q) = 0.5D0*(2*(A+1)+(A+B+2)*(Z-1.D0))
                  ELSEIF (J.GE.2) THEN
                     
                     NUM1 = 1.D0
                     DO I=1,2*(J-1)+A+B+2
                        NUM1 = NUM1*I
                     ENDDO

                     DEN1 = 1.D0
                     DO I=1,2*(J-1)+A+B-1
                        DEN1 = DEN1*I
                     ENDDO

                     COEFF1 = (2*(J-1)+A+B+1)*(A**2-B**2) + Z*NUM1/DEN1
                     COEFF2 = -2*(J-1+A)*(J-1+B)*(2*(J-1)+A+B+2)
                     POLY1  = JACOBI(J,A+1,B+1,Q)
                     POLY2  = JACOBI(J-1,A+1,B+1,Q)
                     DEN = 2*(J-1+1)*(J-1+A+B+1)*(2*(J-1)+A+B)
                     
                     JACOBI(J+1,A+1,B+1,Q) = (COEFF1*POLY1+COEFF2*POLY2)/DEN
                     
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      
C.....Construct and evaluate the orthogonal basis and its derivatives at
C.....the area integral gauss points

      DO Q=1,NAGP(L)+1
         DO I=0,P
            DO J=0,P-I
               
               IF (Q.EQ.NAGP(L)+1) THEN
                  XP = -1.D0/3.D0
                  YP = -1.D0/3.D0
               ELSE
                  XP = XAGP(Q,L)
                  YP = YAGP(Q,L)
               ENDIF

               PHI2(I+1,J+1,Q) = JACOBI(I+1,1,1,Q)*((1.D0-YP)/2.D0)**I
     &              *JACOBI(J+1,2*I+2,1,Q)

               IF (I.EQ.0) THEN
                  DXPHI2(I+1,J+1,Q) = 0.D0
               ELSE
                  DXPHI2(I+1,J+1,Q) = 2.D0/(1.D0-YP)*(I+1.D0)/2.D0*
     &                 JACOBI(I,2,2,Q)*((1-YP)/2.D0)**I*JACOBI(J+1,2*I+2,1,Q)
               ENDIF
               
               IF ((I.EQ.0).AND.(J.EQ.0)) THEN
                  DYPHI2(I+1,J+1,Q) = 0.D0
               ELSEIF ((I.EQ.0).AND.(J.GE.1)) THEN
                  DYPHI2(I+1,J+1,Q) = (J+2)/2.D0*JACOBI(J,3,2,Q)
               ELSEIF ((J.EQ.0).AND.(I.GE.1)) THEN
                  DYPHI2(I+1,J+1,Q) = 2.D0*(1.D0+XP)/(1.D0-YP)**2*
     &                 (I+1)/2.D0*JACOBI(I,2,2,Q)*((1.D0-YP)/2.D0)**I
     &                 *JACOBI(J+1,2*I+2,1,Q) + JACOBI(I+1,1,1,Q)*(-I/2.D0*
     &                 ((1.D0-YP)/2.D0)**(I-1))*JACOBI(J+1,2*I+2,1,Q)
               ELSE

                  DYPHI2(I+1,J+1,Q) = 2.D0*(1.D0+XP)/(1.D0-YP)**2*
     &                 (I+1)/2.D0*JACOBI(I,2,2,Q)*((1.D0-YP)/2.D0)**I
     &                 *JACOBI(J+1,2*I+2,1,Q) + JACOBI(I+1,1,1,Q)*(-I/2.D0*
     &                 ((1.D0-YP)/2.D0)**(I-1))*JACOBI(J+1,2*I+2,1,Q)
     &                 + JACOBI(I+1,1,1,Q)*((1.D0-YP)/2.D0)**I*
     &                 (J+2*I+2)/2.D0*JACOBI(J,2*I+3,2,Q)

               ENDIF
               
               IF (Q.EQ.1) THEN
                  PHI_CORNER1(I+1,J+1,1,P) = (-1.D0)**(I+J+2)
                  PHI_CORNER1(I+1,J+1,2,P) = (-1.D0)**(J+2)
                  IF ((I+1).EQ.1) THEN
                     PHI_CORNER1(I+1,J+1,3,P) = J+1.D0 
                  ELSE
                     PHI_CORNER1(I+1,J+1,3,P) = 0.D0
                  ENDIF
               ENDIF
               
            ENDDO
         ENDDO
      ENDDO
      
C.....Re-order the basis functions into hierarchical order

      DO Q=1,NAGP(L)+1
         M = 1
         DO J=0,P
            DO I=0,J
               JJ = J - I
               IF (Q.LE.NAGP(L)) THEN
                  PHI_AREA(M,Q,P) = PHI2(I+1,JJ+1,Q)
                  DRPHI(M,Q,P) = DXPHI2(I+1,JJ+1,Q)
                  DSPHI(M,Q,P) = DYPHI2(I+1,JJ+1,Q)
                  DO K=1,3
                     PHI_CORNER(M,K,P) = PHI_CORNER1(I+1,JJ+1,K,P)
                  ENDDO
               ELSE
                  PHI_CENTER(M,P) = PHI2(I+1,JJ+1,Q)
               ENDIF
               M = M+1            
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE
