      PROGRAM pointer_test
      IMPLICIT NONE
      INTEGER, PARAMETER :: rp = kind(1d0)

      INTEGER :: m,i,j

      TYPE :: ptr_array
        REAL(rp), POINTER :: ptr
      END TYPE ptr_array

      TYPE(ptr_array), ALLOCATABLE, DIMENSION(:) :: test

      REAL(rp), TARGET, DIMENSION(5,2) :: H
      REAL(rp), DIMENSION(10) :: H2

!       TYPE(ptr_array), POINTER :: ptr2type

      ALLOCATE(test(10))

      m = 0
      DO i = 1,2
        DO j = 1,5
          m = m+1
          H(j,i) = m
          test(m)%ptr => H(j,i)
        ENDDO
      ENDDO

      DO i = 1,5
        PRINT*, (H(i,j),j=1,2)
      ENDDO

      PRINT*, ' ' 

      DO i = 1,m
        PRINT*, test(i)%ptr
      ENDDO

      PRINT*, ' '


      m = 10
      DO i = 1,2
        DO j = 1,5
          m = m+1
          H(j,i) = m
        ENDDO
      ENDDO

      DO i = 1,5
        PRINT*, (H(i,j),j=1,2)
      ENDDO

      PRINT*, ' ' 

      DO i = 1,10
        PRINT*, test(i)%ptr
      ENDDO

      PRINT*, ' '

      DO i = 1,10
        H2(i) = test(i)%ptr*test(i)%ptr
      ENDDO

      DO i = 1,10
        PRINT*, H2(i)
      ENDDO
      
      PRINT*, ' '      
        
!       ptr2type => test
! 
!       DO i = 1,10
!         PRINT*, ptr2type(i)
!       ENDDO
      
      CALL pass_test(test,H)

      END PROGRAM pointer_test
      
      SUBROUTINE pass_test(test_pass,H_pass)
       
      IMPLICIT NONE
      
      INTEGER, PARAMETER :: rp = kind(1d0)      
      INTEGER :: i,m,j
      
      REAL(rp) :: H_pass(5,2)

      TYPE :: ptr_array
        REAL(rp), POINTER :: ptr
      END TYPE ptr_array
      
      TYPE(ptr_array) :: test_pass(10)      

      
      m = 0
      DO i = 1,2
        DO j = 1,5
          m = m+1
          H_pass(j,i) = m
        ENDDO
      ENDDO

      DO i = 1,5
        PRINT*, (H_pass(i,j),j=1,2)
      ENDDO
      
      DO i = 1,10
        PRINT*, test_pass(i)%ptr
!         PRINT*, test_pass(i)
      ENDDO
      
      
      END SUBROUTINE pass_test