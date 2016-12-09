C***********************************************************************
C     
C     SUBROUTINE ORTHOBASIS_PLOT()
C     
C     Evaluates the polynomials for the P degree orthogonal basis
C     for a triangle at the nodal points for plotting.
C     
C     Written by Ethan Kubatko (06-08-2004)
C
C***********************************************************************

      SUBROUTINE ORTHOBASIS_PLOT(XPNP, YPNP, PHI_PLOT,P,NPNP)
      use sizes, only : sz
      IMPLICIT NONE
C.....Input variables
      REAL(SZ), allocatable, intent(in) :: XPNP(:), YPNP(:)
      REAL(SZ), allocatable, intent(inout) :: PHI_PLOT(:,:)
      INTEGER, intent(in) :: P
      INTEGER, intent(inout) :: NPNP
      
C.....Declare local variables

      INTEGER Q,A,B,M,k,i,jj,j
      REAL(SZ) COEFF1,COEFF2,POLY1,POLY2,NUM,DEN,NUM1,DEN1,DEN2,DEN3,Z
      REAL(SZ) XP,YP
      REAL(SZ), allocatable ::  JACOBI(:,:,:,:)
      REAL(SZ), allocatable ::  PHI2(:,:,:) 

C.....Construct and evaluate the necessary Jacobi polynomials at the
C.....area integral gauss points and the barycenter of the element
      
      NPNP = SIZE(XP)
      ALLOCATE(JACOBI(p+1,2*p+3,2,NAGP+1), PHI2(p+1,p+1,NAGP+1) )
      ALLOCATE(PHI_PLOT((p+1)*(p+2)/2,NPNP))

      DO Q=1,NPNP
         DO A=0,2*P+2
            DO B=0,1
               DO J=0,P
                  
                  XP = XAGP(Q)
                  YP = YAGP(Q)  
                  
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

      DO Q=1,NPNP
         DO I=0,P
            DO J=0,P-I
               
               XP = XAGP(Q)
               YP = YAGP(Q)
               
               PHI2(I+1,J+1,Q) = JACOBI(I+1,1,1,Q)*((1.D0-YP)/2.D0)**I
     &              *JACOBI(J+1,2*I+2,1,Q)
               
            ENDDO
         ENDDO
      ENDDO
      
C.....Re-order the basis functions into hierarchical order

      DO Q=1,NPNP
         M = 1
         DO J=0,P
            DO I=0,J
               JJ = J - I

               PHI_PLOT(M,Q) = PHI2(I+1,JJ+1,Q)
               
               M = M+1            
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE
