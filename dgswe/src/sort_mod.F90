      MODULE sort_mod
      
      USE globals, ONLY: rp
      
      IMPLICIT NONE
      
      CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      
      SUBROUTINE insertion_sort(n,val,lst)
      
      INTEGER, INTENT(IN) :: n
      REAL(rp), DIMENSION(:), INTENT(INOUT) :: val
      INTEGER , DIMENSION(:), INTENT(INOUT) :: lst
      
      INTEGER :: i,j
      REAL(rp) :: tmpv
      INTEGER :: tmpl
      
      DO i = 2,n
        tmpv = val(i)
        tmpl = lst(i)
        j = i-1

 shift: DO                       ! DO WHILE (j>=1 .and. val(j)>tmpv)
          IF (j < 1) THEN
            EXIT shift
          ENDIF
          IF (val(j) > tmpv) THEN 
             val(j+1) = val(j)
             lst(j+1) = lst(j)
             j = j-1
           ELSE   
             EXIT shift
           ENDIF
        ENDDO shift
        
        val(j+1) = tmpv
        lst(j+1) = tmpl
      ENDDO


      RETURN
      END SUBROUTINE insertion_sort
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       

      END MODULE sort_mod
