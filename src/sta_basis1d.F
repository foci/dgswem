C***********************************************************************
C
C     SUBROUTINE STA_BASIS1D()
C
C     This subroutine computes the 1D basis functions at a given station
C
C     Written by Younghun Kang (01-28-2020)
C     
C***********************************************************************

      SUBROUTINE STA_BASIS1D( XSTA, YSTA, EDSTA, PHI_STA1D )

      USE GLOBAL
      USE DG
      USE SIZES

      IMPLICIT NONE

      INTEGER J, SZ2, EDSTA
      REAL(SZ) PHI_STA1D(ph+1)
      REAL(SZ) XSTA, YSTA
      REAL(SZ) D1, D2, D3, DMIN
      REAL(8) PT
      REAL(8),Allocatable  :: BASIS(:), DBASIS1(:)
      INTEGER :: i, JMIN
      
      EDSTA = 0
      DMIN = 1.0E+25
      DO J = 1,NCEDS
         IF (CBWIDTH(1,J).eq.0.D0 .or. CBWIDTH(2,J).eq.0.D0) THEN ! skip if routing edge
            CYCLE
         ENDIF
         X1 = X(NCEDNO(1,J))
         X2 = X(NCEDNO(2,J))
         
         Y1 = Y(NCEDNO(1,J))
         Y2 = Y(NCEDNO(2,J))
         
         D1 = sqrt((XSTA-X1)**2 + (YSTA-Y1)**2)
         D2 = sqrt((XSTA-X2)**2 + (YSTA-Y2)**2)
         D3 = sqrt((X1-X2)**2 + (Y1-Y2)**2)
         
         if (abs(D1+D2-D3)/D3.lt.1e-6) then
            ! EDSTA = NCEDN(J)
            EDSTA = J
            DMIN = abs(D1+D2-D3)/D3
         endif 
         
         if (abs(D1+D2-D3)/D3.lt.DMIN) then
            DMIN = abs(D1+D2-D3)/D3
            JMIN = J
         endif
         
      ENDDO
      
      if (EDSTA.eq.0) then
         ! EDSTA = NCEDN(JMIN)
         EDSTA = JMIN
         print *, 'PROGRAM ESTIMATED NEAREST EDGE',NCEDN(JMIN)
         print *, 'PROXIMITY INDEX FOR THIS EDGE EQUALS',DMIN
      endif
      
      ! if (DMIN.lt.1e-6) then
         ! EDSTA = JMIN
      ! endif
      
      X1 = X(NCEDNO(1,EDSTA))
      X2 = X(NCEDNO(2,EDSTA))
      
      Y1 = Y(NCEDNO(1,EDSTA))
      Y2 = Y(NCEDNO(2,EDSTA))
      
      D1 = sqrt((XSTA-X1)**2 + (YSTA-Y1)**2)
      D2 = sqrt((XSTA-X2)**2 + (YSTA-Y2)**2)
      D3 = sqrt((X1-X2)**2 + (Y1-Y2)**2)
      
      
      ! if (NCEDNO(3,EDSTA).eq.0) then
         ! PT = D1/D3*2.D0 - 1
      ! elseif (NCEDNO(3,EDSTA).eq.1) then
         ! PT = 1 - D1/D3*2.D0
      ! endif
      
      PT = D1/(D1+D2)*2.D0 - 1
      ! if (I.eq.3) print*,D1,D2,D3,PT, NCEDNO(1,EDSTA), NCEDNO(2,EDSTA),XEL(I) ,YEL(I)
                                     
      ALLOCATE (BASIS(ph+1),DBASIS1(ph+1))
      CALL ORTHOGONAL_BASIS_1D(PT,ph,BASIS,DBASIS1)
      PHI_STA1D(:) = BASIS

      RETURN
      END SUBROUTINE STA_BASIS1D
      
