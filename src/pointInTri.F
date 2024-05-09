C***********************************************************************
C
C     SUBROUTINE pointInTri( )
C
C     This subroutine determines whether point (p1,p2) is inside a I's element
C     
C     The output FLAG indicates three cases:
C        FLAG = 1  if the given point is inside of element
C        FLAG = 0  if the given point is on the edge of element
C        FLAG = -1 if the given point is outside of element
C
C     Written by Dustin West in MATLAB and copied by Younghun Kang
C
C***********************************************************************

      SUBROUTINE pointInTri(p1,p2,I,FLAG)

C.....Use appropriate modules

      USE GLOBAL
      USE DG
      IMPLICIT NONE
      
C.....Declare local variables
      REAL(SZ), INTENT(IN) :: p1, p2
      INTEGER, INTENT(IN) :: I
      INTEGER, INTENT(OUT) :: FLAG
      
      INTEGER II
      REAL(SZ) D, W(3), q1, q2
      
C.....Retrieve the node coefficients for the given element
      X1 = X(NM(I,1))
      X2 = X(NM(I,2)) - X1
      X3 = X(NM(I,3)) - X1
      
      Y1 = Y(NM(I,1))
      Y2 = Y(NM(I,2)) - Y1
      Y3 = Y(NM(I,3)) - Y1
      
      q1 = p1 - X1
      q2 = p2 - Y1

C.....Compute denominator for calculating barycentric weights
      D = X2*Y3 - X3*Y2
        
C.....Calculate barycentric weights
      W(1) = (q1*(Y2 - Y3) + q2*(X3 - X2) + X2*Y3 - X3*Y2)/D
      W(2) = (q1*Y3 - q2*X3)/D
      W(3) = (q2*X2 - q1*Y2)/D
      
      if ( abs(W(1)+W(2)+W(3) - 1.D0) .gt. 1e-8) then
         print*,'Huh?',abs(W(1)+W(2)+W(3) - 1.D0)
      endif
      
C.....Cap to zero
      DO II = 1, 3
         IF (abs(W(II)).lt.1e-8) THEN
           W(II) = 0.D0
         ENDIF
      ENDDO

C.....Set default flags = 1
      FLAG = 1
      
C.....Check if the point is out of element
      DO II = 1, 3
         IF (W(II).lt.0.D0 .OR. W(II).gt.1.D0) THEN
            FLAG = -1
            EXIT
         ENDIF
      ENDDO
      
C.....Check if the point is on the edge of element
      DO II = 1, 3
         IF (W(II).EQ.0.D0) THEN
            FLAG = 0
            EXIT
         ENDIF
      ENDDO
      
      ! WRITE(999999,*) W(1), W(2), W(3)
      ! WRITE(999999,*) I, p1, p2
      ! WRITE(999999,*) X(NM(I,1)), X(NM(I,2)), X(NM(I,3))
      ! WRITE(999999,*) Y(NM(I,1)), Y(NM(I,2)), Y(NM(I,3))
      ! WRITE(999999,*) INFLAG, ONFLAG
        
      RETURN
      END


