      MODULE PDG_DEBUG
      
      USE SIZES
      USE GLOBAL
      USE MESSENGER_EDGE
      USE COMM_EDGE

#ifdef HAVE_MPI_MOD
      use mpi  
      IMPLICIT NONE
#else
      IMPLICIT NONE
#endif

      CONTAINS

      SUBROUTINE COMPARE_EDGELENGTH()
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif

      INTEGER I,J,K,IEL,IED,IERR

      DO J=1,NEDG_SHARD
        IEL = IEDG_SHARD(1,J)
        IED = IEDG_SHARD(2,J)
        NM1 = NM(IEL,MOD(IED+0,3)+1)
        NM2 = NM(IEL,MOD(IED+1,3)+1)
        ZE_ES(1,J) = SQRT((X(NM2)-X(NM1))**2+(Y(NM2)-Y(NM1))**2)
      ENDDO

      CALL UPDATER_EDG(ZE_ES,QX_ES,QY_ES,ZE_EG,QX_EG,QY_EG,1)

      WRITE(*,*)
      WRITE(*,*) '           (CHECK OF EDGE LENGTH)              '
      WRITE(*,*) '   PROC   EDG     LOCAL    RECEIVED     DIFF   '
      WRITE(*,*) '   ----- ----- ---------- ---------- ----------'
      DO J=1,NEDG_GHOST
        IEL = IEDG_GHOST(1,J)
        IED = IEDG_GHOST(2,J)
        NM1 = NM(IEL,MOD(IED+0,3)+1)
        NM2 = NM(IEL,MOD(IED+1,3)+1)
        WRITE(*,1000) MYPROC, J,
     &       ZE_EG(1,J), SQRT((X(NM2)-X(NM1))**2+(Y(NM2)-Y(NM1))**2),
     &       ZE_EG(1,J) - SQRT((X(NM2)-X(NM1))**2+(Y(NM2)-Y(NM1))**2),
     &       IEL,IED
      ENDDO

 1000 FORMAT('  ',2I6,3E11.3,2I6)
      END SUBROUTINE

      END MODULE PDG_DEBUG
