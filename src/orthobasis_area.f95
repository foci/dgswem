!***********************************************************************
!     
!     SUBROUTINE ORTHOBASIS_AREA()
!     
!     Evaluates the polynomials for the P degree orthogonal basis
!     for a triangle at the area integral gauss points.
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

      SUBROUTINE ORTHOBASIS_AREA(dg,P,L)

!.....Use appropriate modules

      USE GLOBAL
      USE DG

      IMPLICIT NONE

      type (dg_type) :: dg

!.....Declare local variables

      INTEGER Q,A,B,M,k,i,jj,j,P,L
      REAL(SZ) COEFF1,COEFF2,POLY1,POLY2,NUM,DEN,NUM1,DEN1,DEN2,DEN3,Z
      REAL(SZ) XP,YP

!.....Allocate dg%JACOBI

      if (L.eq.1.and.p.eq.1) then

         CALL ALLOC_JACOBI()

      endif

!.....Construct and evaluate the necessary dg%Jacobi polynomials at the
!.....area integral gauss points and the barycenter of the element
      
      DO Q=1,dg%NAGP(L)+1
         DO A=0,2*P+2
            DO B=0,1
               DO J=0,P
                  
                  IF (Q.EQ.dg%NAGP(L)+1) THEN
                     XP = -1.D0/3.D0
                     YP = -1.D0/3.D0
                  ELSE
                     XP = dg%XAGP(Q,L)
                     YP = dg%YAGP(Q,L)
                  ENDIF
                  
                  IF ((A.EQ.0).OR.((A.EQ.1).AND.(B.EQ.1))) THEN
                     Z = 2.D0*(1.D0 + XP)/(1.D0 - YP) - 1.D0
                  ELSE
                     Z = YP
                  ENDIF

                  IF (J.EQ.0) THEN
                     dg%JACOBI(J+1,A+1,B+1,Q) = 1.D0
                  ELSEIF (J.EQ.1) THEN
                     dg%JACOBI(J+1,A+1,B+1,Q) = 0.5D0*(2*(A+1)+(A+B+2)*(Z-1.D0))
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
                     POLY1  = dg%JACOBI(J,A+1,B+1,Q)
                     POLY2  = dg%JACOBI(J-1,A+1,B+1,Q)
                     DEN = 2*(J-1+1)*(J-1+A+B+1)*(2*(J-1)+A+B)
                     
                     dg%JACOBI(J+1,A+1,B+1,Q) = (COEFF1*POLY1+COEFF2*POLY2)/DEN
                     
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      
!.....Construct and evaluate the orthogonal basis and its derivatives at
!.....the area integral gauss points

      DO Q=1,dg%NAGP(L)+1
         DO I=0,P
            DO J=0,P-I
               
               IF (Q.EQ.dg%NAGP(L)+1) THEN
                  XP = -1.D0/3.D0
                  YP = -1.D0/3.D0
               ELSE
                  XP = dg%XAGP(Q,L)
                  YP = dg%YAGP(Q,L)
               ENDIF

               dg%PHI2(I+1,J+1,Q) = dg%JACOBI(I+1,1,1,Q)*((1.D0-YP)/2.D0)**I&
              *dg%JACOBI(J+1,2*I+2,1,Q)

               IF (I.EQ.0) THEN
                  dg%DXPHI2(I+1,J+1,Q) = 0.D0
               ELSE
                  dg%DXPHI2(I+1,J+1,Q) = 2.D0/(1.D0-YP)*(I+1.D0)/2.D0*&
                 dg%JACOBI(I,2,2,Q)*((1-YP)/2.D0)**I*dg%JACOBI(J+1,2*I+2,1,Q)
               ENDIF
               
               IF ((I.EQ.0).AND.(J.EQ.0)) THEN
                  dg%DYPHI2(I+1,J+1,Q) = 0.D0
               ELSEIF ((I.EQ.0).AND.(J.GE.1)) THEN
                  dg%DYPHI2(I+1,J+1,Q) = (J+2)/2.D0*dg%JACOBI(J,3,2,Q)
               ELSEIF ((J.EQ.0).AND.(I.GE.1)) THEN
                  dg%DYPHI2(I+1,J+1,Q) = 2.D0*(1.D0+XP)/(1.D0-YP)**2*&
                 (I+1)/2.D0*dg%JACOBI(I,2,2,Q)*((1.D0-YP)/2.D0)**I&
                 *dg%JACOBI(J+1,2*I+2,1,Q) + dg%JACOBI(I+1,1,1,Q)*(-I/2.D0*&
                 ((1.D0-YP)/2.D0)**(I-1))*dg%JACOBI(J+1,2*I+2,1,Q)
               ELSE

                  dg%DYPHI2(I+1,J+1,Q) = 2.D0*(1.D0+XP)/(1.D0-YP)**2*&
                 (I+1)/2.D0*dg%JACOBI(I,2,2,Q)*((1.D0-YP)/2.D0)**I&
                 *dg%JACOBI(J+1,2*I+2,1,Q) + dg%JACOBI(I+1,1,1,Q)*(-I/2.D0*&
                 ((1.D0-YP)/2.D0)**(I-1))*dg%JACOBI(J+1,2*I+2,1,Q)&
                 + dg%JACOBI(I+1,1,1,Q)*((1.D0-YP)/2.D0)**I*&
                 (J+2*I+2)/2.D0*dg%JACOBI(J,2*I+3,2,Q)

               ENDIF
               
               IF (Q.EQ.1) THEN
                  dg%PHI_CORNER1(I+1,J+1,1,P) = (-1.D0)**(I+J+2)
                  dg%PHI_CORNER1(I+1,J+1,2,P) = (-1.D0)**(J+2)
                  IF ((I+1).EQ.1) THEN
                     dg%PHI_CORNER1(I+1,J+1,3,P) = J+1.D0 
                  ELSE
                     dg%PHI_CORNER1(I+1,J+1,3,P) = 0.D0
                  ENDIF
               ENDIF
               
            ENDDO
         ENDDO
      ENDDO
      
!.....Re-order the basis functions into hierarchical order

      DO Q=1,dg%NAGP(L)+1
         M = 1
         DO J=0,P
            DO I=0,J
               JJ = J - I
               IF (Q.LE.dg%NAGP(L)) THEN
                  dg%PHI_AREA(M,Q,P) = dg%PHI2(I+1,JJ+1,Q)
                  dg%DRPHI(M,Q,P) = dg%DXPHI2(I+1,JJ+1,Q)
                  dg%DSPHI(M,Q,P) = dg%DYPHI2(I+1,JJ+1,Q)
                  DO K=1,3
                     dg%PHI_CORNER(M,K,P) = dg%PHI_CORNER1(I+1,JJ+1,K,P)
                  ENDDO
               ELSE
                  dg%PHI_CENTER(M,P) = dg%PHI2(I+1,JJ+1,Q)
               ENDIF
               M = M+1            
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE
