C***********************************************************************
C
C     SUBROUTINE STA_BASIS()
C
C     This subroutine computes the basis functions at a given station
C
C     Written by Ethan Kubatko
C     
C***********************************************************************

      SUBROUTINE STA_BASIS( D, XSTA, YSTA, ELSTA, PHI_STA )

      USE GLOBAL
      USE DG
      USE SIZES

      IMPLICIT NONE

      INTEGER ELSTA, NSTA, SZ2, ELT
      INTEGER, intent(in) :: D
      REAL(SZ) PHI_STA(DOF)
      REAL(SZ) XSTA, YSTA, AREA
      REAL(SZ) Z1, Z2, TOL
      REAL(8) PT(D)
      REAL(8),Allocatable  :: BASIS(:), DBASIS(:,:)
      INTEGER :: i

C.....Retrieve element coordinates and area

      ELT = EDN(ELSTA) - 2
      if (EDN(ELSTA) .eq. 1) then
         X1 = X(NM(ELSTA,1))
         X2 = X(NM(ELSTA,2))
         X3 = X(NM(ELSTA,3))
   
         Y1 = Y(NM(ELSTA,1))
         Y2 = Y(NM(ELSTA,2))
         Y3 = Y(NM(ELSTA,3))
   
         AREA = 0.5D0*AREAS(ELSTA)
      
C.....Transform to local Z element coordinates

         Z1 = -1.D0/(2.D0*AREA)*( XSTA*(2.D0*Y1 - 2.D0*Y3) +
     &                                        YSTA*(2.D0*X3 - 2.D0*X1) +
     &                   X1*Y2 + X1*Y3 - X2*Y1 + X2*Y3 - X3*Y1 - X3*Y2 )

         Z2 =  1.D0/(2.D0*AREA)*( XSTA*(2.D0*Y1 - 2.D0*Y2) +
     &                                        YSTA*(2.D0*X2 - 2.D0*X1) +
     &                   X1*Y2 + X1*Y3 - X2*Y1 - X2*Y3 - X3*Y1 + X3*Y2 )
      elseif (EDN(ELSTA) .eq. 2) then
         X1 = X(NM(ELSTA,1))
         X2 = X(NM(ELSTA,2))
         X3 = X(NM(ELSTA,3))
         X4 = X(NM(ELSTA,4))
   
         Y1 = Y(NM(ELSTA,1))
         Y2 = Y(NM(ELSTA,2))
         Y3 = Y(NM(ELSTA,3))
         Y4 = Y(NM(ELSTA,4))
   
         AREA = 0.5D0*AREAS(ELSTA)
      
C.....Transform to local Z element coordinates
      ! FIND APPROPRIATE EXPRESSIONS FOR QUADS!
         Z1 = -1.D0/(2.D0*AREA)*( XSTA*(2.D0*Y1 - 2.D0*Y3) +
     &                                        YSTA*(2.D0*X3 - 2.D0*X1) +
     &                   X1*Y2 + X1*Y3 - X2*Y1 + X2*Y3 - X3*Y1 - X3*Y2 )

         Z2 =  1.D0/(2.D0*AREA)*( XSTA*(2.D0*Y1 - 2.D0*Y2) +
     &                                        YSTA*(2.D0*X2 - 2.D0*X1) +
     &                   X1*Y2 + X1*Y3 - X2*Y1 - X2*Y3 - X3*Y1 + X3*Y2 )
      endif
C.....Compute the basis functions at that point and store
      PT(1) = Z1
      PT(2) = Z2

      SZ2 = (PDG_EL(ELSTA)+2)*(PDG_EL(ELSTA)+1)**(D-1)/2
      Allocate ( BASIS(SZ2),DBASIS(SZ2,D) )
      
      do i = 1,D
         if (abs(NINT(PT(i)) - PT(i)) .le. 1.0d-12) then
            PT(i) = REAL(NINT(PT(i)),8)
         endif
      enddo
   
      CALL ORTHOGONAL_BASIS( ELT,PT,PDG_EL(ELSTA),D,BASIS,DBASIS )

      PHI_STA(:) = BASIS

      RETURN
      END SUBROUTINE STA_BASIS
      
