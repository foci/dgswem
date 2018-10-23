C***********************************************************************
C
C     SUBROUTINE STA_BASIS()
C
C     This subroutine computes the basis functions at a given station
C
C     Written by Ethan Kubatko
C     
C***********************************************************************

      SUBROUTINE STA_BASIS( ELEM, DIM, XSTA, YSTA, ELSTA, PHI_STA )

      USE GLOBAL
      USE DG
      USE SIZES

      IMPLICIT NONE

      INTEGER ELSTA, NSTA, ELEM, SZ2
	  INTEGER, intent(in) :: DIM
      REAL(SZ) PHI_STA(DOF)
      REAL(SZ) XSTA, YSTA, AREA
      REAL(SZ) Z1, Z2, TOL
      REAL(8) PT(DIM)
      REAL(8),Allocatable  :: BASIS(:), DBASIS(:,:)
      INTEGER :: i

C.....Retrieve element coordinates and area

      X1 = X(NM(ELSTA,1))
      X2 = X(NM(ELSTA,2))
      X3 = X(NM(ELSTA,3))

      Y1 = Y(NM(ELSTA,1))
      Y2 = Y(NM(ELSTA,2))
      Y3 = Y(NM(ELSTA,3))

      AREA = 0.5D0*AREAS(ELSTA)
      
C.....Transform to local Z element coordinates

      Z1 = -1.D0/(2.D0*AREA)*( XSTA*(2.D0*Y1 - 2.D0*Y3) +
     &                         YSTA*(2.D0*X3 - 2.D0*X1) +
     &                   X1*Y2 + X1*Y3 - X2*Y1 + X2*Y3 - X3*Y1 - X3*Y2 )

      Z2 =  1.D0/(2.D0*AREA)*( XSTA*(2.D0*Y1 - 2.D0*Y2) +
     &                         YSTA*(2.D0*X2 - 2.D0*X1) +
     &                   X1*Y2 + X1*Y3 - X2*Y1 - X2*Y3 - X3*Y1 + X3*Y2 )

C.....Compute the basis functions at that point and store
      PT(1) = Z1
      PT(2) = Z2

      SZ2 = (PDG_EL(ELSTA)+2)*(PDG_EL(ELSTA)+1)**(DIM-1)/2
      Allocate ( BASIS(SZ2),DBASIS(SZ2,DIM) )
      
      do i = 1,DIM
         if (abs(NINT(PT(i)) - PT(i)) .le. 1.0d-12) then
            PT(i) = REAL(NINT(PT(i)),8)
         endif
      enddo
   
      CALL ORTHOGONAL_BASIS( ELEM,PT,PDG_EL(ELSTA),DIM,BASIS,DBASIS )

      PHI_STA(:) = BASIS

      RETURN
      END SUBROUTINE STA_BASIS
      
