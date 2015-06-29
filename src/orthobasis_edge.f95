!***********************************************************************
!     
!     SUBROUTINE ORTHOBASIS_EDGE()
!     
!     Evaluates the polynomials for the P degree orthogonal basis
!     for a triangle at the global_here%edge integral gauss points.
!     
!     Reference for orthogonal basis:
!     
!     "A Discontinuous Galerkin Method for the Viscous MHD Equations",
!     Journal of Computational Physics 152, 608-641 (1999)
!     
!     Written by Ethan Kubatko (06-08-2004)
!     01-10-2011 - cem - adapted for p_enrichment
!     
!***********************************************************************

      SUBROUTINE ORTHOBASIS_EDGE(dg_here,P,L)

!.....Use appropriate modules

      USE GLOBAL
      USE DG

      IMPLICIT NONE

      type (dg_type) :: dg_here

      INTEGER Q, A, B, M, II, III,jj,i,j,P,L
      REAL(SZ) COEFF1,COEFF2,POLY1,POLY2,NUM,DEN,NUM1,DEN1,DEN2,DEN3,Z
      REAL(SZ) XGP,YGP
      
!.....Loop over the 3 edges

      III = 4
      DO 100 II = 1,3

!.....Construct and evaluate the necessary dg_here%Jacobi polynomials at the
!.....global_here%edge integral gauss points and the global_here%edge mid point
         
         DO Q=1,dg_here%NEGP(L)+1
            IF (II.EQ.1) THEN
               IF (Q.LT.(dg_here%NEGP(L)+1)) THEN
                  XGP = -dg_here%XEGP(Q,L)
                  YGP =  dg_here%XEGP(Q,L)
               ELSE
                  XGP = 0.D0
                  YGP = 0.D0
               ENDIF
            ELSEIF (II.EQ.2) THEN
               IF (Q.LT.(dg_here%NEGP(L)+1)) THEN
                  XGP = -1.D0
                  YGP = -dg_here%XEGP(Q,L)
               ELSE
                  XGP = -1.D0
                  YGP =  0.D0
               ENDIF
            ELSEIF (II.EQ.3) THEN
               IF (Q.LT.(dg_here%NEGP(L)+1)) THEN
                  XGP = dg_here%XEGP(Q,L)
                  YGP = -1.D0
               ELSE
                  XGP =  0.D0
                  YGP = -1.D0
               ENDIF
            ENDIF
            DO A=0,2*P+2
               DO B=0,1
                  DO J=0,P
                     
                     IF ((A.EQ.0).OR.((A.EQ.1).AND.(B.EQ.1))) THEN
                        Z = 2.D0*(1.D0 + XGP)/(1.D0 - YGP) - 1.D0
                     ELSE
                        Z = YGP
                     ENDIF

                     IF (J.EQ.0) THEN
                        dg_here%JACOBI(J+1,A+1,B+1,Q) = 1.D0
                     ELSEIF (J.EQ.1) THEN
                        dg_here%JACOBI(J+1,A+1,B+1,Q) =(2*(A+1)+(A+B+2)*(Z-1.D0))/2.D0
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
                        POLY1  = dg_here%JACOBI(J,A+1,B+1,Q)
                        POLY2  = dg_here%JACOBI(J-1,A+1,B+1,Q)
                        DEN = 2*(J-1+1)*(J-1+A+B+1)*(2*(J-1)+A+B)
                        
                        dg_here%JACOBI(J+1,A+1,B+1,Q) = (COEFF1*POLY1+COEFF2*POLY2)&
                       /DEN
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         
!.....Construct and evaluate the orthogonal basis at the global_here%edge quadrature
!.....points and the global_here%edge midpoint

         DO Q=1,dg_here%NEGP(L)+1
            
            IF (Q.LT.(dg_here%NEGP(L)+1)) THEN
               IF (II.EQ.1) YGP =  dg_here%XEGP(Q,L)
               IF (II.EQ.2) YGP = -dg_here%XEGP(Q,L)
               IF (II.EQ.3) YGP = -1.D0
            ELSE
               IF (II.EQ.1) YGP =  0.D0
               IF (II.EQ.2) YGP =  0.D0
               IF (II.EQ.3) YGP = -1.D0
            ENDIF

            DO I=0,P
               DO J=0,P-I
                  dg_here%PHI2(I+1,J+1,Q) = dg_here%JACOBI(I+1,1,1,Q)*((1.D0-YGP)/2.D0)**I&
                 *dg_here%JACOBI(J+1,2*I+2,1,Q)
               ENDDO
            ENDDO
         ENDDO

!.....Re-order the basis functions into hierarchical order and compute
!.....the inverse of the mass matrix (diagonal)

         DO Q = 1,dg_here%NEGP(L)+1
            M = 1
            DO J = 0,P
               DO I = 0,J
                  JJ = J - I
                  IF (Q.LT.(dg_here%NEGP(L)+1)) THEN
                     dg_here%PHI_EDGE(M,Q,II,P) = dg_here%PHI2(I+1,JJ+1,Q)
                  ELSE
                     dg_here%PHI_MID(M,II,P) = dg_here%PHI2(I+1,JJ+1,Q)
                  ENDIF
                  dg_here%M_INV(M,P) = ((2.D0*I+1.D0)*(2.D0*JJ+2.D0*I+2.D0)/4.D0)
                  M = M+1
               ENDDO
            ENDDO
            IF (P.GT.1) THEN
               IF (Q.LT.(dg_here%NEGP(L)+1)) THEN
                  dg_here%PSI_CHECK(1,III) = -1.D0/2.D0*(XGP + YGP )
                  dg_here%PSI_CHECK(2,III) =  1.D0/2.D0*(XGP + 1.D0)
                  dg_here%PSI_CHECK(3,III) =  1.D0/2.D0*(YGP + 1.D0)
                  III = III + 1
               ENDIF
            ENDIF
         ENDDO

 100  CONTINUE

      RETURN
      END SUBROUTINE
