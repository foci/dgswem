
      SUBROUTINE MESSAGE_INIT (MNPROC,MYRANK)
      IMPLICIT NONE
#ifdef CRAY
      INCLUDE './mpif.h.t3e'
      INTEGER*8 IERR, MNPROC, MYRANK
#else
      INCLUDE 'mpif.h'
      INTEGER IERR, MNPROC, MYRANK
#endif
C
C--------------------------------------------------------------------------
C  Routine performs following steps:
C   (1)  initialize MPI, 
C   (2)  get number of processors,
C   (3)  get MPI rank of processor, 
C  vjp  8/06/1999
C--------------------------------------------------------------------------
C
      CALL MPI_INIT(IERR)                               ! Initialize MPI
c     print *, "status from mpi_init: ",ierr
C
      CALL MPI_COMM_SIZE (MPI_COMM_WORLD,MNPROC,IERR)   ! Get number of procs
C
      CALL MPI_COMM_RANK (MPI_COMM_WORLD,MYRANK,IERR)   ! Get MPI rank
c     print *, "status from mpi_comm_rank: ",ierr
c     print *, "myrank = ",MYRANK
C
      RETURN
      END SUBROUTINE



      SUBROUTINE SYNCHRONIZE()
      IMPLICIT NONE
#ifdef CRAY
      INCLUDE './mpif.h.t3e'
      INTEGER*8 IERR
#else
      INCLUDE 'mpif.h'
      INTEGER IERR    
#endif
C
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
C
      RETURN
      END SUBROUTINE



      SUBROUTINE MESSAGE_FINI()
      IMPLICIT NONE
#ifdef CRAY
      INCLUDE './mpif.h.t3e'
      INTEGER*8 IERR
#else
      INCLUDE 'mpif.h'
      INTEGER IERR    
#endif
C
      CALL MPI_FINALIZE(IERR)
C
      RETURN
      END SUBROUTINE