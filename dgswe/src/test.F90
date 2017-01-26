      PROGRAM test
      IMPLICIT NONE
      INTEGER :: i,j
      REAL :: array(1000,10)

      OPEN(unit=10,file='test.out')

      DO i = 1,10
        DO j = 1,1000
          array(j,i) = i*j/j
        ENDDO
      ENDDO

      DO i = 1,10
        WRITE(10,*) (array(j,i), j = 1,1000)
      ENDDO
      
       DO i = 1,10
        DO j = 1,1000
          WRITE(10,*) array(j,i)
        ENDDO
      ENDDO

      DO i = 1,10
        WRITE(10,*) array(:,i)
      ENDDO

      END PROGRAM test