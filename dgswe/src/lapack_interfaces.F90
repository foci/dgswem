MODULE lapack_interfaces

      INTERFACE 

        SUBROUTINE DGESV(n,nrhs,a,lda,ipiv,b,ldb,info)

          ! DGESV computes the solution to a real system of linear equations 
          !    A * X = B,
          ! where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
          !
          ! The LU decomposition with partial pivoting and row interchanges is
          ! used to factor A as
          !     A = P * L * U,
          ! where P is a permutation matrix, L is unit lower triangular, and U is
          ! upper triangular.  The factored form of A is then used to solve the
          ! system of equations A * X = B.


          INTEGER, INTENT(IN) :: n                                  ! The number of linear equations, i.e., the order of the
                                                                    ! matrix A.  N >= 0.

          INTEGER, INTENT(IN) :: nrhs                               ! The number of right hand sides, i.e., the number of columns
                                                                    ! of the matrix B.  NRHS >= 0.

          INTEGER, INTENT(IN) :: lda                                ! The leading dimension of the array A.  LDA >= max(1,N).

          DOUBLE PRECISION, DIMENSION(lda,*), INTENT(INOUT) :: a    ! dimension (LDA,N)
                                                                    ! On entry, the N-by-N coefficient matrix A.
                                                                    ! On exit, the factors L and U from the factorization
                                                                    ! A = P*L*U; the unit diagonal elements of L are not stored.

          INTEGER, DIMENSION(*), INTENT(OUT) :: ipiv                ! dimension (N)
                                                                    ! The pivot indices that define the permutation matrix P;
                                                                    ! row i of the matrix was interchanged with row IPIV(i).

          INTEGER, INTENT(IN) :: ldb                                ! The leading dimension of the array B.  LDB >= max(1,N).

          DOUBLE PRECISION, DIMENSION(ldb,*), INTENT(INOUT) :: b    ! dimension (LDB,NRHS)
                                                                    ! On entry, the N-by-NRHS matrix of right hand side matrix B.
                                                                    ! On exit, if INFO = 0, the N-by-NRHS solution matrix X.

          INTEGER, INTENT(OUT) :: info                              ! = 0:  successful exit
                                                                    ! < 0:  if INFO = -i, the i-th argument had an illegal value
                                                                    ! > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
                                                                    !       has been completed, but the factor U is exactly
                                                                    !       singular, so the solution could not be computed.
        END SUBROUTINE DGESV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE DGETRF(m,n,a,lda,ipiv,info)  

          ! DGETRF computes an LU factorization of a general M-by-N matrix A
          ! using partial pivoting with row interchanges.
          !
          ! The factorization has the form
          !    A = P * L * U
          ! where P is a permutation matrix, L is lower triangular with unit
          ! diagonal elements (lower trapezoidal if m > n), and U is upper
          ! triangular (upper trapezoidal if m < n).
          
          INTEGER, INTENT(IN) :: m                                 ! The number of rows of the matrix A.  M >= 0.
          
          INTEGER, INTENT(IN) :: n                                 ! The number of columns of the matrix A.  N >= 0.

          INTEGER, INTENT(IN) :: lda                               ! The leading dimension of the array A.  LDA >= max(1,M).

          DOUBLE PRECISION, DIMENSION(lda,*), INTENT(INOUT) :: a   ! dimension (LDA,N)
                                                                   ! On entry, the M-by-N matrix to be factored.
                                                                   ! On exit, the factors L and U from the factorization
                                                                   ! A = P*L*U; the unit diagonal elements of L are not stored.

          INTEGER, DIMENSION(*), INTENT(OUT) :: ipiv               ! dimension (min(M,N))
                                                                   ! The pivot indices; for 1 <= i <= min(M,N), row i of the
                                                                   ! matrix was interchanged with row IPIV(i).

          INTEGER, INTENT(OUT) :: info                             !  = 0:  successful exit
                                                                   ! < 0:  if INFO = -i, the i-th argument had an illegal value
                                                                   ! > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
                                                                   !       has been completed, but the factor U is exactly
                                                                   !       singular, and division by zero will occur if it is used
                                                                   !       to solve a system of equations.        
        END SUBROUTINE DGETRF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE DGETRS(trans,n,nrhs,a,lda,ipiv,b,ldb,info)

          ! DGETRS solves a system of linear equations
          !    A * X = B  or  A**T * X = B
          ! with a general N-by-N matrix A using the LU factorization computed
          ! by DGETRF.

          CHARACTER, INTENT(IN) :: trans                           ! Specifies the form of the system of equations:
                                                                   ! = 'N':  A * X = B  (No transpose)
                                                                   ! = 'T':  A**T* X = B  (Transpose)
                                                                   ! = 'C':  A**T* X = B  (Conjugate transpose = Transpose)

          INTEGER, INTENT(IN) :: n                                 ! The order of the matrix A.  N >= 0. 

          INTEGER, INTENT(IN) :: nrhs                              ! The number of right hand sides, i.e., the number of columns
                                                                   ! of the matrix B.  NRHS >= 0.

          INTEGER, INTENT(IN) :: lda                               ! The leading dimension of the array A.  LDA >= max(1,N).

          DOUBLE PRECISION, DIMENSION(lda,*), INTENT(IN) :: a      ! dimension (LDA,N)
                                                                   ! The factors L and U from the factorization A = P*L*U
                                                                   ! as computed by DGETRF.

          INTEGER, DIMENSION(*), INTENT(IN) :: ipiv               ! dimension (N)
                                                                   ! The pivot indices from DGETRF; for 1<=i<=N, row i of the
                                                                   ! matrix was interchanged with row IPIV(i).

          INTEGER, INTENT(IN) :: ldb                               ! The leading dimension of the array B.  LDB >= max(1,N).

          DOUBLE PRECISION, DIMENSION(ldb,*), INTENT(INOUT) :: b   ! dimension (LDB,NRHS)
                                                                   ! On entry, the right hand side matrix B.
                                                                   ! On exit, the solution matrix X.

          INTEGER, INTENT(OUT) :: info                             ! The leading dimension of the array B.  LDB >= max(1,N).

        END SUBROUTINE DGETRS
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
        
        SUBROUTINE DGETRI(n, a, lda, ipiv, work, lwork, info)
        
        !  DGETRI computes the inverse of a matrix using the LU factorization
        !  computed by DGETRF.
        !
        !  This method inverts U and then computes inv(A) by solving the system
        !  inv(A)*L = inv(U) for inv(A).   
        
        INTEGER, INTENT(IN) :: n                                   ! The order of the matrix A.  N >= 0.
        
        DOUBLE PRECISION, DIMENSION(lda,*), INTENT(INOUT) :: a     ! dimension (LDA,N)
                                                                   ! On entry, the factors L and U from the factorization
                                                                   ! A = P*L*U as computed by DGETRF.
                                                                   ! On exit, if INFO = 0, the inverse of the original matrix A.
                                                                   
        INTEGER, INTENT(IN) :: lda                                 ! The leading dimension of the array A.  LDA >= max(1,N).
        
        INTEGER, DIMENSION(*), INTENT(IN) :: ipiv                  ! dimension (N)
                                                                   ! The pivot indices from DGETRF; for 1<=i<=N, row i of the
                                                                   ! matrix was interchanged with row IPIV(i).
                                                                   
        DOUBLE PRECISION, DIMENSION(*), INTENT(INOUT) :: work      ! dimension (MAX(1,LWORK))
                                                                   ! On exit, if INFO=0, then WORK(1) returns the optimal LWORK.                                                                   
                                                                   
        INTEGER, INTENT(IN) :: lwork                                ! The dimension of the array WORK.  LWORK >= max(1,N).
                                                                   ! For optimal performance LWORK >= N*NB, where NB is
                                                                   ! the optimal blocksize returned by ILAENV.
                                                                   !
                                                                   ! If LWORK = -1, then a workspace query is assumed; the routine
                                                                   ! only calculates the optimal size of the WORK array, returns
                                                                   ! this value as the first entry of the WORK array, and no error
                                                                   ! message related to LWORK is issued by XERBLA.
                                                                   
        INTEGER, INTENT(OUT) :: info                               ! = 0:  successful exit
                                                                   ! < 0:  if INFO = -i, the i-th argument had an illegal value
                                                                   ! > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
                                                                   !       singular and its inverse could not be computed.
                                                                   
        
        END SUBROUTINE DGETRI
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                      
        
        SUBROUTINE DGTSV(n, nrhs, dl, d, du, b, ldb, info)
        
          !  DGTSV  solves the equation
          !
          !     A*X = B,
          !
          !  where A is an n by n tridiagonal matrix, by Gaussian elimination with
          !  partial pivoting.
          !
          !  Note that the equation  A**T*X = B  may be solved by interchanging the
          !  order of the arguments DU and DL.        
        
          INTEGER, INTENT(IN) :: n                                   ! The order of the matrix A.  N >= 0.
        
          INTEGER, INTENT(IN) :: nrhs                                ! The number of right hand sides, i.e., the number of columns
                                                                     ! of the matrix B.  NRHS >= 0.
                                                                   
          DOUBLE PRECISION, DIMENSION(*), INTENT(INOUT) :: dl        ! dimension (N-1)
                                                                     ! On entry, DL must contain the (n-1) sub-diagonal elements of A.
                                                                     ! On exit, DL is overwritten by the (n-2) elements of the
                                                                     ! second super-diagonal of the upper triangular matrix U from
                                                                     ! the LU factorization of A, in DL(1), ..., DL(n-2).
                                                                   
          DOUBLE PRECISION, DIMENSION(*), INTENT(INOUT) :: d         ! dimension (N)
                                                                     ! On entry, D must contain the diagonal elements of A.
                                                                     ! On exit, D is overwritten by the n diagonal elements of U.
                                                                   
          DOUBLE PRECISION, DIMENSION(*), INTENT(INOUT) :: du        ! dimension (N-1)
                                                                     ! On entry, DU must contain the (n-1) super-diagonal elements of A.
                                                                     ! On exit, DU is overwritten by the (n-1) elements of the first
                                                                     ! super-diagonal of U.
                                                                   
          DOUBLE PRECISION, DIMENSION(ldb,*), INTENT(INOUT) :: b     ! dimension (LDB,NRHS)
                                                                     ! On entry, the N by NRHS matrix of right hand side matrix B.
                                                                     ! On exit, if INFO = 0, the N by NRHS solution matrix X.
                                                                   
          INTEGER, INTENT(IN) :: ldb                                 ! The leading dimension of the array B.  LDB >= max(1,N).
        
          INTEGER, INTENT(OUT) :: info                               ! = 0: successful exit
                                                                     ! < 0: if INFO = -i, the i-th argument had an illegal value
                                                                     ! > 0: if INFO = i, U(i,i) is exactly zero, and the solution
                                                                     !      has not been computed.  The factorization has not been
                                                                     !      completed unless i = N.
          
        END SUBROUTINE DGTSV            

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

        SUBROUTINE DSTEQR(compz, n, d, e, z, ldz, work, info)
        
          !  DSTEQR computes all eigenvalues and, optionally, eigenvectors of a
          !  symmetric tridiagonal matrix using the implicit QL or QR method.
          !  The eigenvectors of a full or band symmetric matrix can also be found
          !  if DSYTRD or DSPTRD or DSBTRD has been used to reduce this matrix to
          !  tridiagonal form.          
          
          CHARACTER, INTENT(IN) :: compz                           ! CHARACTER*1
                                                                   ! = 'N':  Compute eigenvalues only.
                                                                   ! = 'V':  Compute eigenvalues and eigenvectors of the original 
                                                                   !         symmetric matrix.  On entry, Z must contain the
                                                                   !         orthogonal matrix used to reduce the original matrix
                                                                   !         to tridiagonal form.
                                                                   ! = 'I':  Compute eigenvalues and eigenvectors of the
                                                                   !         tridiagonal matrix.  Z is initialized to the identity
                                                                   !         matrix.
                                                                   
          INTEGER, INTENT(IN) :: n                                 ! The order of the matrix.  N >= 0.
          
          DOUBLE PRECISION, DIMENSION(*), INTENT(INOUT) :: d       ! dimension (N)
                                                                   ! On entry, the diagonal elements of the tridiagonal matrix.
                                                                   ! On exit, if INFO = 0, the eigenvalues in ascending order.                                                                   
                                                                   
          DOUBLE PRECISION, DIMENSION(*), INTENT(INOUT) :: e       ! dimension (N-1)
                                                                   ! On entry, the (n-1) subdiagonal elements of the tridiagonal matrix.
                                                                   ! On exit, E has been destroyed.
                                                                   
          DOUBLE PRECISION, DIMENSION(ldz,*), INTENT(INOUT) :: z   ! dimension (LDZ, N)
                                                                   ! On entry, if  COMPZ = 'V', then Z contains the orthogonal
                                                                   ! matrix used in the reduction to tridiagonal form.
                                                                   ! On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the
                                                                   ! orthonormal eigenvectors of the original symmetric matrix,
                                                                   ! and if COMPZ = 'I', Z contains the orthonormal eigenvectors
                                                                   ! of the symmetric tridiagonal matrix.
                                                                   ! If COMPZ = 'N', then Z is not referenced.
                                                                   
          INTEGER, INTENT(IN) :: ldz                               ! The leading dimension of the array Z.  LDZ >= 1, and if
                                                                   ! eigenvectors are desired, then  LDZ >= max(1,N).
                                                                   
          DOUBLE PRECISION, DIMENSION(*), INTENT(INOUT) :: work    ! dimension (max(1,2*N-2))
                                                                   ! If COMPZ = 'N', then WORK is not referenced.
                                                                   
          INTEGER, INTENT(OUT) :: info                             ! = 0:  successful exit
                                                                   ! < 0:  if INFO = -i, the i-th argument had an illegal value
                                                                   ! > 0:  the algorithm has failed to find all the eigenvalues in
                                                                   ! a total of 30*N iterations; if INFO = i, then i
                                                                   ! elements of E have not converged to zero; on exit, D
                                                                   ! and E contain the elements of a symmetric tridiagonal
                                                                   ! matrix which is orthogonally similar to the original
                                                                   ! matrix.
                                                                   
          
          
          
        
        END SUBROUTINE  
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE DGELS(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info )
        
          ! DGELS solves overdetermined or underdetermined real linear systems
          ! involving an M-by-N matrix A, or its transpose, using a QR or LQ
          ! factorization of A.  It is assumed that A has full rank.       
          

          CHARACTER, INTENT(IN) :: trans                           ! = 'N': the linear system involves A;
                                                                   ! = 'T': the linear system involves A**T.
                                                                   
          INTEGER, INTENT(IN) :: m                                 ! The number of rows of the matrix A.  M >= 0.
          
          INTEGER, INTENT(IN) :: n                                 ! The number of columns of the matrix A.  N >= 0.
          
          INTEGER, INTENT(IN) :: nrhs                              ! The number of right hand sides, i.e., the number of
                                                                   ! columns of the matrices B and X. NRHS >=0
                                                                   
          DOUBLE PRECISION, DIMENSION(lda,*), INTENT(INOUT) :: a   ! dimension (LDA,N)
                                                                   ! On entry, the M-by-N matrix A.
                                                                   ! On exit,
                                                                   ! if M >= N, A is overwritten by details of its QR
                                                                   !            factorization as returned by DGEQRF;
                                                                   ! if M <  N, A is overwritten by details of its LQ
                                                                   !            factorization as returned by DGELQF.
                                                                   
          INTEGER, INTENT(IN) :: lda                               ! The leading dimension of the array A.  LDA >= max(1,M).
          
          DOUBLE PRECISION, DIMENSION(ldb,*), INTENT(INOUT) :: b   ! dimension (LDB,NRHS)
                                                                   ! On entry, the matrix B of right hand side vectors, stored
                                                                   ! columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS
                                                                   ! if TRANS = 'T'.                                                                  
                                                                   ! On exit, if INFO = 0, B is overwritten by the solution
                                                                   ! vectors, stored columnwise:                                                                   
                                                                   ! if TRANS = 'N' and m >= n, rows 1 to n of B contain the least                                                                   
                                                                   ! squares solution vectors; the residual sum of squares for the                                                                   
                                                                   ! solution in each column is given by the sum of squares of                                                                   
                                                                   ! elements N+1 to M in that column;                                                                   
                                                                   ! if TRANS = 'N' and m < n, rows 1 to N of B contain the                                                                   
                                                                   ! minimum norm solution vectors;                                                                   
                                                                   ! if TRANS = 'T' and m >= n, rows 1 to M of B contain the
                                                                   ! minimum norm solution vectors;
                                                                   ! if TRANS = 'T' and m < n, rows 1 to M of B contain the
                                                                   ! least squares solution vectors; the residual sum of squares
                                                                   ! for the solution in each column is given by the sum of
                                                                   ! squares of elements M+1 to N in that column.                                                             
                                                                                                                                                                                                        
          INTEGER, INTENT(IN) :: ldb                               ! The leading dimension of the array B. LDB >= MAX(1,M,N).
          
          DOUBLE PRECISION, DIMENSION(*), INTENT(INOUT) :: work    ! dimension (MAX(1,LWORK))
                                                                   ! On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
                                                                   
          INTEGER, INTENT(IN) :: lwork                             ! The dimension of the array WORK.
                                                                   ! LWORK >= max( 1, MN + max( MN, NRHS ) ).
                                                                   ! For optimal performance,
                                                                   ! LWORK >= max( 1, MN + max( MN, NRHS )*NB ).
                                                                   ! where MN = min(M,N) and NB is the optimum block size.
                                                                   !
                                                                   ! If LWORK = -1, then a workspace query is assumed; the routine
                                                                   ! only calculates the optimal size of the WORK array, returns
                                                                   ! this value as the first entry of the WORK array, and no error
                                                                   ! message related to LWORK is issued by XERBLA.                                                                             

          INTEGER, INTENT(OUT) :: info                             ! = 0:  successful exit
                                                                   ! < 0:  if INFO = -i, the i-th argument had an illegal value
                                                                   ! > 0:  if INFO =  i, the i-th diagonal element of the
                                                                   !       triangular factor of A is zero, so that A does not have
                                                                   !       full rank; the least squares solution could not be
                                                                   !       computed.
        
        
        END SUBROUTINE DGELS        

      END INTERFACE

END MODULE lapack_interfaces