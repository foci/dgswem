      MODULE spherical_mod
      
      USE globals, ONLY: rp,deg2rad      
      
      CONTAINS

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

      
      SUBROUTINE cpp_transformation(coord_sys,r_earth,slam0,sphi0,nn,xy)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: coord_sys
      REAL(rp), INTENT(IN) :: r_earth
      REAL(rp), INTENT(IN) :: slam0
      REAL(rp), INTENT(IN) :: sphi0        
      INTEGER, INTENT(IN) :: nn
      REAL(rp), DIMENSION(:,:), INTENT(INOUT) :: xy     
      
      INTEGER :: i
      
      ! Transform from lon,lat (in degrees) to x,y
      
      IF (coord_sys /= 1) THEN
      
        DO i = 1,nn
          xy(1,i) = xy(1,i)*deg2rad
          xy(2,i) = xy(2,i)*deg2rad
        
          xy(1,i) = r_earth*(xy(1,i)-slam0)*cos(sphi0)
          xy(2,i) = r_earth*xy(2,i)
        ENDDO      
      
      ENDIF
      
      RETURN
      END SUBROUTINE cpp_transformation
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      


      SUBROUTINE cpp_factor(coord_sys,r_earth,slam0,sphi0,ypt,Sp)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: coord_sys
      REAL(rp), INTENT(IN) :: r_earth
      REAL(rp), INTENT(IN) :: slam0
      REAL(rp), INTENT(IN) :: sphi0        
      REAL(rp), INTENT(IN) :: ypt
      REAL(rp), INTENT(OUT) :: Sp

      IF (coord_sys == 1) THEN
        Sp = 1d0
      ELSE
        Sp = cos(sphi0)/cos(ypt/r_earth)        
      ENDIF      

      RETURN
      END SUBROUTINE cpp_factor
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      SUBROUTINE cpp_inverse(coord_sys,r_earth,slam0,sphi0,nn,xy)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: coord_sys
      REAL(rp), INTENT(IN) :: r_earth
      REAL(rp), INTENT(IN) :: slam0
      REAL(rp), INTENT(IN) :: sphi0        
      INTEGER, INTENT(IN) :: nn
      REAL(rp), DIMENSION(:,:), INTENT(INOUT) :: xy  
      
      INTEGER :: i
      
      
      IF (coord_sys /= 1) THEN
        DO i = 1,nn      
          xy(1,i) = xy(1,i)/(r_earth*cos(sphi0)) + slam0
          xy(2,i) = xy(2,i)/r_earth
          
          xy(1,i) = xy(1,i)/deg2rad
          xy(2,i) = xy(2,i)/deg2rad
        ENDDO      
      ENDIF
      
      
      END SUBROUTINE cpp_inverse


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      END MODULE spherical_mod