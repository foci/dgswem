!******************************************************************************
!  last changes in this file VERSION 12.sb01                                  *
!  S.Bunya changed this file a bit. 07/13/2005                                *
!  S.Bunya changed this file a bit. 01/04/2007                                *
!                                                                             *
!     01-10-2011 - cem - adapted for p_enrichment and multicomponent          * C****************************************************************************** 
! 
      MODULE MESSENGER_ELEM
!
      USE SIZES
      USE GLOBAL
      USE DG
      USE DIFF45_41,ONLY : MNODES


#ifdef HAVE_MPI_MOD
      use mpi  
      IMPLICIT NONE
#else
      IMPLICIT NONE
#endif
!
!--------------------------------------------------------------------------
!  This module supplies the MPI Message-Passing Interface
!
!  Uses asynchronous and persistent communication with buffer packing
!  as performance enhancements for "cluster" architectures.
!
!  vjp  8/29/1999
!--------------------------------------------------------------------------
!
!

!
!  Message-Passing Array space
!
!sb-PDG1
      PUBLIC

      INTEGER, SAVE :: MPI_COMM_ADCIRC
      INTEGER, SAVE ::  COMM_COMP     
      INTEGER, SAVE ::  GROUP_WORLD, GROUP_COMP

      INTEGER,SAVE,PRIVATE :: REALTYPE, DBLETYPE, COMM   
      INTEGER,SAVE,PRIVATE ::  NEIGHPROC_R, NEIGHPROC_S, RDIM, IERR
      INTEGER,SAVE,PRIVATE ::  TAG = 101
      LOGICAL,SAVE,ALLOCATABLE :: RESELEM(:)
!
      INTEGER, PRIVATE, ALLOCATABLE ::IPROC_R(:),IPROC_S(:),NELEMLOC(:)
      INTEGER, PRIVATE, ALLOCATABLE ::NELEMSEND(:), NELEMRECV(:),ISENDLOC(:,:), IBELONGTO(:)
      INTEGER, PRIVATE, ALLOCATABLE ::IRECVLOC(:,:), ISENDBUF(:,:), IRECVBUF(:,:)
!
      INTEGER, PRIVATE, ALLOCATABLE :: REQ_I1(:), REQ_I2(:)
      INTEGER, PRIVATE, ALLOCATABLE :: STAT_I1(:,:), STAT_I2(:,:)
      INTEGER, PRIVATE, ALLOCATABLE :: REQ_R1(:), REQ_R2(:), REQ_R3(:),&
                                 REQ_LZ(:) 
      INTEGER, PRIVATE, ALLOCATABLE :: STAT_R1(:,:), STAT_R2(:,:), &
                                 STAT_R3(:,:), STAT_LZ(:,:)
      INTEGER, PRIVATE, ALLOCATABLE :: REQ_R3D(:), STAT_R3D(:,:)
      INTEGER, PRIVATE, ALLOCATABLE :: REQ_C3D(:), STAT_C3D(:,:)
      INTEGER, PRIVATE, ALLOCATABLE :: INDEX(:)
      REAL(SZ), PRIVATE,ALLOCATABLE :: SENDBUF(:,:), RECVBUF(:,:)
!--
!

!---------------------end of data declarations--------------------------------C


      CONTAINS


      SUBROUTINE MSG_TYPES_ELEM()
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
!
#ifdef  REAL4
      REALTYPE = MPI_REAL
      DBLETYPE = MPI_DOUBLE_PRECISION
#else
      REALTYPE = MPI_DOUBLE_PRECISION
      DBLETYPE = MPI_DOUBLE_PRECISION
#endif
!
      RETURN
      END  SUBROUTINE


      SUBROUTINE MESSAGE_INIT (s)
!--------------------------------------------------------------------------
!  Routine performs following steps:
!   (1)  initialize MPI, 
!   (2)  get number of processors,
!   (3)  get MPI rank of processor 
!  vjp  8/06/1999
!--------------------------------------------------------------------------
        use sizes
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
!ND
      type (sizes_type) :: s
      Integer I
      INTEGER, ALLOCATABLE :: RANKS(:)

      CALL MPI_INIT(IERR)                               ! Initialize MPI
      CALL MPI_COMM_SIZE (MPI_COMM_WORLD,s%MNPROC,IERR)   ! Get number of procs
      CALL MPI_COMM_RANK (MPI_COMM_WORLD,s%MYPROC,IERR)   ! Get MPI rank
!ND
      ALLOCATE(RANKS(s%MNPROC+1))
      DO I=1,s%MNPROC
         RANKS(I)=I-1
      ENDDO
      CALL MPI_COMM_GROUP(MPI_COMM_WORLD,GROUP_WORLD,IERR)
      CALL MPI_GROUP_INCL(GROUP_WORLD,s%MNPROC,RANKS,GROUP_COMP,IERR)
      CALL MPI_COMM_CREATE(MPI_COMM_WORLD,GROUP_COMP,COMM_COMP,IERR)
      DEALLOCATE(RANKS)
      COMM = COMM_COMP
      RETURN
      END SUBROUTINE

!$$$      IMPLICIT NONE
!$$$#ifndef HAVE_MPI_MOD
!$$$      include 'mpif.h'
!$$$#endif
!$$$      CALL MPI_INIT(IERR)                               ! Initialize MPI
!$$$      CALL MPI_COMM_SIZE (MPI_COMM_WORLD,MNPROC,IERR)   ! Get number of procs
!$$$      CALL MPI_COMM_RANK (MPI_COMM_WORLD,MYPROC,IERR)   ! Get MPI rank
!$$$      RETURN
!$$$      END SUBROUTINE

      SUBROUTINE ErrorElevSum( ErrorElevExceeded )
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      INTEGER ErrorElevExceeded !=1 if this subdomain has exceeded warning elev
      INTEGER SumErrorElevExceeded !sum total of all flags from all subdomains
      INTEGER kount             ! to avoid compiler bug on certain platforms
      
      SumErrorElevExceeded = 0
      kount=1
      call MPI_ALLREDUCE( ErrorElevExceeded, SumErrorElevExceeded, kount,&
     MPI_INTEGER, MPI_SUM, MPI_COMM_world, ierr)
      ErrorElevExceeded = SumErrorElevExceeded
      END SUBROUTINE ErrorElevSum

 
      SUBROUTINE MESSAGE_FINI (s)
!--------------------------------------------------------------------------
!  Shutdown MPI library.
!--------------------------------------------------------------------------
        use sizes
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      type (sizes_type) :: s
      INTEGER IERR,I
!
      CALL MPI_FINALIZE(IERR)
      IF (s%MYPROC.EQ.0)  &
  print *, "MPI terminated with Status = ",IERR      
      RETURN
      END SUBROUTINE

 
      SUBROUTINE MESSAGE_ABORT (s)
!--------------------------------------------------------------------------
!  Shutdown MPI library.
!--------------------------------------------------------------------------
        use sizes
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      type (sizes_type) :: s
      INTEGER IERR,I
!
      CALL MPI_ABORT(MPI_COMM_WORLD,IERR)
      IF (s%MYPROC.EQ.0)  &
  print *, "MPI Aborted with Status = ",IERR      
      RETURN
      END SUBROUTINE



      SUBROUTINE MSG_TABLE_ELEM (s) 
!
!--------------------------------------------------------------------------
!  Routine preforms following steps:
!
!   (1) Read Message-Passing Information from file "DG.18"
!   (2) Determine resident nodes: RESNODE(I) is true  if I is resident node
!   (3) Determine ghost nodes:    RESNODE(I) is false if I is ghost node    
!   (4) Determine number of neighbor subdomains
!   (5) MPI rank of each neighbor and number of ghosts nodes to receive
!   (6) Read Message-Passing Receive List
!   (7) MPI rank of each neighbor and number of ghosts nodes to send
!   (8) Read Message-Passing Send List
!
!  vjp  8/06/1999
!--------------------------------------------------------------------------
!
      use sizes
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      type (sizes_type) :: s
      INTEGER :: IDPROC,NLOCAL,I,J,JJ,NEIGH
!
      OPEN(18,FILE=s%DIRNAME(1:s%LNAME)//'/'//'DG.18')
!
      READ(18,3010) IDPROC,NLOCAL
!
      ALLOCATE ( NELEMLOC(NLOCAL) )
!
      READ(18,1130) (NELEMLOC(I), I=1,NLOCAL)
!
      ALLOCATE ( IBELONGTO(S%MNE),RESELEM(S%MNE) )
!
      DO I=1,S%MNE
         IBELONGTO(I) = 0
      ENDDO
      DO I=1,NLOCAL
         IBELONGTO(NELEMLOC(I)) = IDPROC + 1
      ENDDO
      DO I=1, S%MNE
         IF (IBELONGTO(I)-1.EQ.s%MYPROC) THEN
           RESELEM(I) = .TRUE.
         ELSE 
           RESELEM(I) = .FALSE.
         ENDIF
      ENDDO
!
      READ(18,3010) NEIGHPROC_R,NEIGHPROC_S
!
      RDIM = NEIGHPROC_R + NEIGHPROC_S
      ALLOCATE( INDEX(RDIM) )
!
      ALLOCATE( IPROC_R(NEIGHPROC_R),NELEMRECV(NEIGHPROC_R) )
      ALLOCATE( IRECVLOC(S%MNE,NEIGHPROC_R) )
!
      DO JJ=1,NEIGHPROC_R
         J = MOD(JJ-1+s%MYPROC,NEIGHPROC_R)+1
         READ(18,3010) IPROC_R(J),NELEMRECV(J)
         READ(18,1130) (IRECVLOC(I,J), I=1,NELEMRECV(J))
      ENDDO
!
      ALLOCATE( IPROC_S(NEIGHPROC_S),NELEMSEND(NEIGHPROC_S) )
      ALLOCATE( ISENDLOC(S%MNE,NEIGHPROC_S) )
!
      DO JJ=1,NEIGHPROC_S
         J = MOD(JJ-1+s%MYPROC,NEIGHPROC_S)+1
         READ(18,3010) IPROC_S(J),NELEMSEND(J)
         READ(18,1130) (ISENDLOC(I,J), I=1,NELEMSEND(J))
      ENDDO
!
      CLOSE(18)
      RETURN
!
1130  FORMAT(8X,9I8)
3010  FORMAT(8X,2I8)
      END SUBROUTINE


      SUBROUTINE MESSAGE_START_ELEM (dg_here,s)
!
!--------------------------------------------------------------------------
!  Routine preforms following steps:
!   (1)  allocate message-passing space
!   (2)  setup MPI data structures for "persistent" message-passing.
!
!  vjp  8/06/1999
!--------------------------------------------------------------------------
!
        use dg
        use sizes
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      type (sizes_type) :: s
      type (dg_type) :: dg_here
      INTEGER :: J,NCOMMELEM_S,NCOMMELEM_R
!
      NCOMMELEM_R = 0
      NCOMMELEM_S = 0
      DO J=1,NEIGHPROC_R
        NCOMMELEM_R = NCOMMELEM_R + NELEMRECV(J)
      ENDDO
      DO J=1,NEIGHPROC_S
        NCOMMELEM_S = NCOMMELEM_S + NELEMSEND(J)
      ENDDO
!
      ALLOCATE ( ISENDBUF(NCOMMELEM_S*2,NEIGHPROC_S) )
      ALLOCATE ( IRECVBUF(NCOMMELEM_R*2,NEIGHPROC_R) )
!
      IF (s%C3D) THEN
!         ALLOCATE ( SENDBUF(2*MNP*MNODES,NEIGHPROC) )
!         ALLOCATE ( RECVBUF(2*MNP*MNODES,NEIGHPROC) )
        STOP 'NOT UPDATED'
      ELSE
         ALLOCATE ( SENDBUF(NCOMMELEM_S*DG_HERE%DOFH*4,NEIGHPROC_S) )
         ALLOCATE ( RECVBUF(NCOMMELEM_R*DG_HERE%DOFH*4,NEIGHPROC_R) )
      ENDIF
!
      ALLOCATE ( REQ_I1(RDIM),REQ_I2(RDIM) )
      ALLOCATE ( REQ_R1(RDIM),REQ_R2(RDIM),REQ_R3(RDIM),REQ_LZ(RDIM) )
!
      ALLOCATE ( STAT_I1(MPI_STATUS_SIZE,RDIM),       &
           STAT_I2(MPI_STATUS_SIZE,RDIM) )

      ALLOCATE ( STAT_R1(MPI_STATUS_SIZE,RDIM),       &
           STAT_R2(MPI_STATUS_SIZE,RDIM),&
           STAT_R3(MPI_STATUS_SIZE,RDIM),&
           STAT_LZ(MPI_STATUS_SIZE,RDIM) )
!
      IF (s%C3D) THEN
!         ALLOCATE ( REQ_R3D(RDIM) )
!         ALLOCATE ( STAT_R3D(MPI_STATUS_SIZE,RDIM) )
!         ALLOCATE ( REQ_C3D(RDIM) )
!         ALLOCATE ( STAT_C3D(MPI_STATUS_SIZE,RDIM) )
        STOP 'NOT UPDATED'
      ENDIF
!
             !  Setup persistent structures for integer arrays
!
      DO J=1,NEIGHPROC_R
         CALL MPI_RECV_INIT ( IRECVBUF(1,J), NELEMRECV(J), &
     MPI_INTEGER,IPROC_R(J), TAG, MPI_COMM_WORLD,&
     REQ_I1(J),IERR)
      ENDDO
      DO J=1,NEIGHPROC_S
         CALL MPI_SEND_INIT ( ISENDBUF(1,J), NELEMSEND(J), &
    MPI_INTEGER,IPROC_S(J), TAG, MPI_COMM_WORLD,&
    REQ_I1(J+NEIGHPROC_R),IERR )
      ENDDO
!
!
      DO J=1,NEIGHPROC_R
         CALL MPI_RECV_INIT ( IRECVBUF(1,J), 2*NELEMRECV(J), &
     MPI_INTEGER,IPROC_R(J), TAG, MPI_COMM_WORLD,&
     REQ_I2(J),IERR)
      ENDDO
      DO J=1,NEIGHPROC_S
         CALL MPI_SEND_INIT ( ISENDBUF(1,J), 2*NELEMSEND(J), &
    MPI_INTEGER,IPROC_S(J), TAG, MPI_COMM_WORLD,&
    REQ_I2(J+NEIGHPROC_R),IERR )
      ENDDO
!
            !  Setup persistent structures for real arrays
!
      DO J=1,NEIGHPROC_R
         CALL MPI_RECV_INIT ( RECVBUF(1,J), DG_HERE%DOFH*NELEMRECV(J), &
     REALTYPE,IPROC_R(J), TAG, MPI_COMM_WORLD,&
     REQ_R1(J),IERR)
      ENDDO
      DO J=1,NEIGHPROC_S
         CALL MPI_SEND_INIT ( SENDBUF(1,J), DG_HERE%DOFH*NELEMSEND(J), &
     REALTYPE,IPROC_S(J), TAG, MPI_COMM_WORLD,&
     REQ_R1(J+NEIGHPROC_R),IERR)
      ENDDO
!
      DO J=1,NEIGHPROC_R
         CALL MPI_RECV_INIT ( RECVBUF(1,J), 2*DG_HERE%DOFH*NELEMRECV(J), &
     REALTYPE,IPROC_R(J), TAG, MPI_COMM_WORLD,&
     REQ_R2(J),IERR)
      ENDDO
      DO J=1,NEIGHPROC_S
         CALL MPI_SEND_INIT ( SENDBUF(1,J), 2*DG_HERE%DOFH*NELEMSEND(J), &
     REALTYPE,IPROC_S(J), TAG, MPI_COMM_WORLD,&
     REQ_R2(J+NEIGHPROC_R),IERR)
      ENDDO
!
      DO J=1,NEIGHPROC_R
         CALL MPI_RECV_INIT ( RECVBUF(1,J), 3*DG_HERE%DOFH*NELEMRECV(J), &
     REALTYPE,IPROC_R(J), TAG, MPI_COMM_WORLD,&
     REQ_R3(J),IERR)
      ENDDO
      DO J=1,NEIGHPROC_S
         CALL MPI_SEND_INIT ( SENDBUF(1,J), 3*DG_HERE%DOFH*NELEMSEND(J), &
     REALTYPE,IPROC_S(J), TAG, MPI_COMM_WORLD,&
     REQ_R3(J+NEIGHPROC_R),IERR)
      ENDDO
!
      DO J=1,NEIGHPROC_R
         CALL MPI_RECV_INIT ( RECVBUF(1,J), 4*DG_HERE%DOFH*NELEMRECV(J), &
     REALTYPE,IPROC_R(J), TAG, MPI_COMM_WORLD,&
     REQ_LZ(J),IERR)
      ENDDO
      DO J=1,NEIGHPROC_S
         CALL MPI_SEND_INIT ( SENDBUF(1,J), 4*DG_HERE%DOFH*NELEMSEND(J), &
     REALTYPE,IPROC_S(J), TAG, MPI_COMM_WORLD,&
     REQ_LZ(J+NEIGHPROC_R),IERR)
      ENDDO
!
      IF (s%C3D) THEN
!         DO J=1,NEIGHPROC  
!            CALL MPI_RECV_INIT ( RECVBUF(1,J), MNODES*NNODRECV(J), 
!     &        REALTYPE,IPROC(J), TAG, MPI_COMM_WORLD,
!     &        REQ_R3D(J),IERR)
!         ENDDO
!         DO J=1,NEIGHPROC  
!            CALL MPI_SEND_INIT ( SENDBUF(1,J), MNODES*NNODSEND(J), 
!     &        REALTYPE,IPROC(J), TAG, MPI_COMM_WORLD,
!     &        REQ_R3D(J+NEIGHPROC),IERR)
!         ENDDO
!         DO J=1,NEIGHPROC  
!            CALL MPI_RECV_INIT ( RECVBUF(1,J), 2*MNODES*NNODRECV(J), 
!     &        REALTYPE,IPROC(J), TAG, MPI_COMM_WORLD,
!     &        REQ_C3D(J),IERR)
!         ENDDO
!         DO J=1,NEIGHPROC  
!            CALL MPI_SEND_INIT ( SENDBUF(1,J), 2*MNODES*NNODSEND(J), 
!     &        REALTYPE,IPROC(J), TAG, MPI_COMM_WORLD,
!     &        REQ_C3D(J+NEIGHPROC),IERR)
!         ENDDO
        STOP 'NOT UPDATED'
      ENDIF
!
      RETURN
      END SUBROUTINE


      SUBROUTINE UPDATEI_ELEM( IVEC1, IVEC2, NMSG )
!
!--------------------------------------------------------------------------
!  Update 1 or 2 Integer Arrays's Ghost Cells using asynchronous
!  and persistent message-passing.
!
!  vjp  8/06/1999
!--------------------------------------------------------------------------
! 
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
!
      INTEGER,   INTENT(IN) :: NMSG
      Real(sz),   INTENT(INOUT) :: IVEC1(:),IVEC2(:)
      INTEGER :: N,I,J,NCOUNT,NFINI,TOT
!
                             !..Pack 1 or 2 Messages
      DO J=1,NEIGHPROC_S
         NCOUNT = 0
         DO I=1,NELEMSEND(J)
            NCOUNT = NCOUNT+1
            SENDBUF(NCOUNT,J)=IVEC1(ISENDLOC(I,J))
         ENDDO
         IF (NMSG.GT.1) THEN
           DO I=1,NELEMSEND(J)
              NCOUNT = NCOUNT+1
              SENDBUF(NCOUNT,J)=IVEC2(ISENDLOC(I,J))
           ENDDO
         ENDIF
      ENDDO
!                     
                          ! Send/receive messages to/from all neighbors
      IF (NMSG.EQ.1) THEN
        CALL MPI_STARTALL ( RDIM, REQ_R1, IERR )
      ELSE
        CALL MPI_STARTALL ( RDIM, REQ_R2, IERR )
      ENDIF
!
                          !..Unpack Received messages as they arrive  

      IF (NMSG.EQ.1) THEN   
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDEX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R1,NFINI,INDEX,STAT_R1,IERR )
           TOT = TOT + NFINI
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                     NCOUNT = NCOUNT+1
                     IVEC1(IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
      ELSE
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDEX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R2,NFINI,INDEX,STAT_R2,IERR )
           TOT = TOT + NFINI
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                     NCOUNT = NCOUNT+1
                     IVEC1(IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
                  ENDDO
                  DO I=1,NELEMRECV(J)
                     NCOUNT = NCOUNT+1
                     IVEC2(IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
      ENDIF
! 
 999  continue
      RETURN
      END SUBROUTINE


      SUBROUTINE UPDATER_ELEM(dg_here, VEC1, VEC2, VEC3, IRK, NMSG )
!
!--------------------------------------------------------------------------
!  Update 1, 2, or 3 Real Arrays's Ghost Cells using asynchronous
!  and persistent message-passing.
!
!  vjp  8/06/1999
!--------------------------------------------------------------------------
! 
        use dg
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
!
      type (dg_type) :: dg_here
      INTEGER,  INTENT(IN) ::  IRK, NMSG
      REAL(SZ), INTENT(INOUT) ::  VEC1(:,:,:),VEC2(:,:,:),VEC3(:,:,:)
      INTEGER :: N,I,J,K,NCOUNT,NFINI,TOT
!
                             !..Pack 1, 2, or 3 Messages
      DO J=1,NEIGHPROC_S
         NCOUNT = 0
         DO I=1,NELEMSEND(J)
            DO K=1,DG_HERE%DOFH
              NCOUNT = NCOUNT+1
              SENDBUF(NCOUNT,J)=VEC1(K,ISENDLOC(I,J),IRK)
            ENDDO
         ENDDO
         IF (NMSG.GT.1) THEN
           DO I=1,NELEMSEND(J)
             DO K=1,DG_HERE%DOFH
               NCOUNT = NCOUNT+1
               SENDBUF(NCOUNT,J)=VEC2(K,ISENDLOC(I,J),IRK)
             ENDDO
           ENDDO
         ENDIF
         IF (NMSG.GT.2) THEN
           DO I=1,NELEMSEND(J)
             DO K=1,DG_HERE%DOFH
               NCOUNT = NCOUNT+1
               SENDBUF(NCOUNT,J)=VEC3(K,ISENDLOC(I,J),IRK)
             ENDDO
           ENDDO
         ENDIF
      ENDDO
!                    
              ! Send/receive messages to/from all neighbors
!
      IF (NMSG.EQ.1) THEN
        CALL MPI_STARTALL ( RDIM, REQ_R1, IERR )
      ELSEIF (NMSG.EQ.2) THEN
        CALL MPI_STARTALL ( RDIM, REQ_R2, IERR )
      ELSE
        CALL MPI_STARTALL ( RDIM, REQ_R3, IERR )
      ENDIF
              !..Unpack Received messages as they arrive     
!
      IF (NMSG.EQ.1) THEN   
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDEX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R1,NFINI,INDEX,STAT_R1,IERR )
           TOT = TOT + NFINI
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                    DO K=1,DG_HERE%DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(K,IRECVLOC(I,J),IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
        GOTO 999
      ELSEIF (NMSG.EQ.2) THEN
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDEX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R2,NFINI,INDEX,STAT_R2,IERR )
           TOT = TOT + NFINI
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                    DO K=1,DG_HERE%DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(K,IRECVLOC(I,J),IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DG_HERE%DOFH
                      NCOUNT = NCOUNT+1
                      VEC2(K,IRECVLOC(I,J),IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
        GOTO 999
      ELSE
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDEX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R3,NFINI,INDEX,STAT_R3,IERR )
           TOT = TOT + NFINI
!debug     print *, myproc, tot,nfini,index(1),index(2)
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                    DO K=1,DG_HERE%DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(K,IRECVLOC(I,J),IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DG_HERE%DOFH
                      NCOUNT = NCOUNT+1
                      VEC2(K,IRECVLOC(I,J),IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DG_HERE%DOFH
                      NCOUNT = NCOUNT+1
                      VEC3(K,IRECVLOC(I,J),IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
        GOTO 999
      ENDIF
! 
 999  continue
      RETURN
      END SUBROUTINE


      SUBROUTINE UPDATELZ_ELEM(dg_here, LZ )
!
!--------------------------------------------------------------------------
!  Update LZ Real Arrays's Ghost Cells using asynchronous
!  and persistent message-passing.
!
!  sb  1/04/2007
!--------------------------------------------------------------------------
! 
        use dg
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
!
      type (dg_type) :: dg_here
      REAL(SZ), INTENT(INOUT) ::  LZ(:,:,:,:)
      INTEGER :: N,I,J,K,NCOUNT,NFINI,TOT
!
                             !..Pack 1, 2, or 3 Messages
      DO J=1,NEIGHPROC_S
         NCOUNT = 0
         DO I=1,NELEMSEND(J)
            DO K=1,DG_HERE%DOFH
              SENDBUF(NCOUNT+1,J)=LZ(K,1,1,ISENDLOC(I,J))
              SENDBUF(NCOUNT+2,J)=LZ(K,1,2,ISENDLOC(I,J))
              SENDBUF(NCOUNT+3,J)=LZ(K,2,1,ISENDLOC(I,J))
              SENDBUF(NCOUNT+4,J)=LZ(K,2,2,ISENDLOC(I,J))
              NCOUNT = NCOUNT+4
            ENDDO
         ENDDO
      ENDDO
!                    
              ! Send/receive messages to/from all neighbors
!
      CALL MPI_STARTALL ( RDIM, REQ_LZ, IERR )

              !..Unpack Received messages as they arrive     
!
      TOT = 0
      DO WHILE (TOT.LT.RDIM)
        DO N=1, RDIM
          INDEX(N) = 0
        ENDDO
        CALL MPI_WAITSOME( RDIM,REQ_LZ,NFINI,INDEX,STAT_LZ,IERR )
        TOT = TOT + NFINI
        DO N=1, NFINI
          IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
            IF (INDEX(N).LE.NEIGHPROC_R) THEN
              J = INDEX(N)
              NCOUNT = 0
              DO I=1,NELEMRECV(J)
                DO K=1,DG_HERE%DOFH
                  LZ(K,1,1,IRECVLOC(I,J)) = RECVBUF(NCOUNT+1,J)
                  LZ(K,1,2,IRECVLOC(I,J)) = RECVBUF(NCOUNT+2,J)
                  LZ(K,2,1,IRECVLOC(I,J)) = RECVBUF(NCOUNT+3,J)
                  LZ(K,2,2,IRECVLOC(I,J)) = RECVBUF(NCOUNT+4,J)
                  NCOUNT = NCOUNT+4
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE

      SUBROUTINE UPDATEMZ_ELEM(dg_here, MZ )
!
!--------------------------------------------------------------------------
!  Update MZ Real Arrays's Ghost Cells using asynchronous
!  and persistent message-passing.
!
!  sb  1/04/2007
!--------------------------------------------------------------------------
! 
        use dg
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
!
      type (dg_type) :: dg_here
      REAL(SZ), INTENT(INOUT) ::  MZ(:,:,:,:)
      INTEGER :: N,I,J,K,NCOUNT,NFINI,TOT
!
                             !..Pack 1, 2, or 3 Messages
      DO J=1,NEIGHPROC_S
         NCOUNT = 0
         DO I=1,NELEMSEND(J)
            DO K=1,dg_here%DOFH
              SENDBUF(NCOUNT+1,J)=MZ(K,1,1,ISENDLOC(I,J))
              SENDBUF(NCOUNT+2,J)=MZ(K,2,1,ISENDLOC(I,J))
              NCOUNT = NCOUNT+2
            ENDDO
         ENDDO
      ENDDO
!                    
              ! Send/receive messages to/from all neighbors
!
      CALL MPI_STARTALL ( RDIM, REQ_LZ, IERR )

              !..Unpack Received messages as they arrive     
!
      TOT = 0
      DO WHILE (TOT.LT.RDIM)
        DO N=1, RDIM
          INDEX(N) = 0
        ENDDO
        CALL MPI_WAITSOME( RDIM,REQ_LZ,NFINI,INDEX,STAT_LZ,IERR )
        TOT = TOT + NFINI
        DO N=1, NFINI
          IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
            IF (INDEX(N).LE.NEIGHPROC_R) THEN
              J = INDEX(N)
              NCOUNT = 0
              DO I=1,NELEMRECV(J)
                DO K=1,dg_here%DOFH
                  MZ(K,1,1,IRECVLOC(I,J)) = RECVBUF(NCOUNT+1,J)
                  MZ(K,2,1,IRECVLOC(I,J)) = RECVBUF(NCOUNT+2,J)
                  NCOUNT = NCOUNT+2
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE


      SUBROUTINE UPDATER_ELEM_MOD(dg_here, VEC1, VEC2, VEC3, IRK, NMSG )
!
!--------------------------------------------------------------------------
!  Update 1, 2, or 3 Real Arrays's Ghost Cells using asynchronous
!  and persistent message-passing.
!
!  vjp  8/06/1999
!--------------------------------------------------------------------------
! 
        use dg
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
!
      type (dg_type) :: dg_here
      INTEGER,  INTENT(IN) ::  IRK, NMSG
      REAL(SZ), INTENT(INOUT) ::  VEC1(:,:,:),VEC2(:,:,:),VEC3(:,:,:)
      INTEGER :: N,I,J,JJ,K,NCOUNT,NFINI,TOT
!
                             !..Pack 1, 2, or 3 Messages
!      DO JJ=1,NEIGHPROC_S
!         J = MOD(JJ-1+MYPROC,NEIGHPROC_S)+1
      DO J=1,NEIGHPROC_S
         NCOUNT = 0
         DO I=1,NELEMSEND(J)
            DO K=1,DG_HERE%DOFH
              NCOUNT = NCOUNT+1
              SENDBUF(NCOUNT,J)=VEC1(K,ISENDLOC(I,J),IRK)
            ENDDO
         ENDDO
         IF (NMSG.GT.1) THEN
           DO I=1,NELEMSEND(J)
             DO K=1,DG_HERE%DOFH
               NCOUNT = NCOUNT+1
               SENDBUF(NCOUNT,J)=VEC2(K,ISENDLOC(I,J),IRK)
             ENDDO
           ENDDO
         ENDIF
         IF (NMSG.GT.2) THEN
           DO I=1,NELEMSEND(J)
             DO K=1,DG_HERE%DOFH
               NCOUNT = NCOUNT+1
               SENDBUF(NCOUNT,J)=VEC3(K,ISENDLOC(I,J),IRK)
             ENDDO
           ENDDO
         ENDIF

         ! Start sending a message
         IF (NMSG.EQ.1) THEN
            CALL MPI_START ( REQ_R1(J+NEIGHPROC_R), IERR )
         ELSEIF (NMSG.EQ.2) THEN
            CALL MPI_START ( REQ_R2(J+NEIGHPROC_R), IERR )
         ELSE
            CALL MPI_START ( REQ_R3(J+NEIGHPROC_R), IERR )
         ENDIF
      ENDDO
!                    
              ! Send/receive messages to/from all neighbors
!
      IF (NMSG.EQ.1) THEN
        CALL MPI_STARTALL ( RDIM-NEIGHPROC_S, REQ_R1, IERR )
      ELSEIF (NMSG.EQ.2) THEN
        CALL MPI_STARTALL ( RDIM-NEIGHPROC_S, REQ_R2, IERR )
      ELSE
        CALL MPI_STARTALL ( RDIM-NEIGHPROC_S, REQ_R3, IERR )
      ENDIF
!
              !..Unpack Received messages as they arrive     
!
      IF (NMSG.EQ.1) THEN   
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDEX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R1,NFINI,INDEX,STAT_R1,IERR )
           TOT = TOT + NFINI
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                    DO K=1,DG_HERE%DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(K,IRECVLOC(I,J),IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
        GOTO 999
      ELSEIF (NMSG.EQ.2) THEN
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDEX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R2,NFINI,INDEX,STAT_R2,IERR )
           TOT = TOT + NFINI
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                    DO K=1,DG_HERE%DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(K,IRECVLOC(I,J),IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DG_HERE%DOFH
                      NCOUNT = NCOUNT+1
                      VEC2(K,IRECVLOC(I,J),IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
        GOTO 999
      ELSE
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDEX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R3,NFINI,INDEX,STAT_R3,IERR )
           TOT = TOT + NFINI
!debug     print *, myproc, tot,nfini,index(1),index(2)
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                    DO K=1,DG_HERE%DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(K,IRECVLOC(I,J),IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DG_HERE%DOFH
                      NCOUNT = NCOUNT+1
                      VEC2(K,IRECVLOC(I,J),IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DG_HERE%DOFH
                      NCOUNT = NCOUNT+1
                      VEC3(K,IRECVLOC(I,J),IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
        GOTO 999
      ENDIF
! 
 999  continue
      RETURN
      END SUBROUTINE

      SUBROUTINE UPDATER_ELEM_MOD2(dg_here, VEC1, VEC2, VEC3, IRK, NMSG )
!
!--------------------------------------------------------------------------
!  Update 1, 2, or 3 Real Arrays's Ghost Cells using asynchronous
!  and persistent message-passing.
!
!  vjp  8/06/1999
!--------------------------------------------------------------------------
! 
        use dg
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      type (dg_type) :: dg_here
!
      INTEGER,  INTENT(IN) ::  IRK, NMSG
      REAL(SZ), INTENT(INOUT) ::  VEC1(:,:,:),VEC2(:,:,:),VEC3(:,:,:)
      INTEGER :: N,I,J,JJ,K,NCOUNT,NFINI,TOT
!
                             !..Pack 1, 2, or 3 Messages
!      DO JJ=1,NEIGHPROC_S
!         J = MOD(JJ-1+MYPROC,NEIGHPROC_S)+1
      DO J=1,NEIGHPROC_S
         NCOUNT = 0
         DO I=1,NELEMSEND(J)
            DO K=1,DG_HERE%DOFH
              NCOUNT = NCOUNT+1
              SENDBUF(NCOUNT,J)=VEC1(ISENDLOC(I,J),K,IRK)
            ENDDO
         ENDDO
         IF (NMSG.GT.1) THEN
           DO I=1,NELEMSEND(J)
             DO K=1,DG_HERE%DOFH
               NCOUNT = NCOUNT+1
               SENDBUF(NCOUNT,J)=VEC2(ISENDLOC(I,J),K,IRK)
             ENDDO
           ENDDO
         ENDIF
         IF (NMSG.GT.2) THEN
           DO I=1,NELEMSEND(J)
             DO K=1,DG_HERE%DOFH
               NCOUNT = NCOUNT+1
               SENDBUF(NCOUNT,J)=VEC3(ISENDLOC(I,J),K,IRK)
             ENDDO
           ENDDO
         ENDIF

         ! Start sending a message
         IF (NMSG.EQ.1) THEN
            CALL MPI_START ( REQ_R1(J+NEIGHPROC_R), IERR )
         ELSEIF (NMSG.EQ.2) THEN
            CALL MPI_START ( REQ_R2(J+NEIGHPROC_R), IERR )
         ELSE
            CALL MPI_START ( REQ_R3(J+NEIGHPROC_R), IERR )
         ENDIF
      ENDDO
!                    
              ! Send/receive messages to/from all neighbors
!
      IF (NMSG.EQ.1) THEN
        CALL MPI_STARTALL ( RDIM-NEIGHPROC_S, REQ_R1, IERR )
      ELSEIF (NMSG.EQ.2) THEN
        CALL MPI_STARTALL ( RDIM-NEIGHPROC_S, REQ_R2, IERR )
      ELSE
        CALL MPI_STARTALL ( RDIM-NEIGHPROC_S, REQ_R3, IERR )
      ENDIF
!
              !..Unpack Received messages as they arrive     
!
      IF (NMSG.EQ.1) THEN   
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDEX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R1,NFINI,INDEX,STAT_R1,IERR )
           TOT = TOT + NFINI
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                    DO K=1,DG_HERE%DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(IRECVLOC(I,J),K,IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
        GOTO 9998
      ELSEIF (NMSG.EQ.2) THEN
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDEX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R2,NFINI,INDEX,STAT_R2,IERR )
           TOT = TOT + NFINI
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                    DO K=1,DG_HERE%DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(IRECVLOC(I,J),K,IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DG_HERE%DOFH
                      NCOUNT = NCOUNT+1
                      VEC2(IRECVLOC(I,J),K,IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
        GOTO 9998
      ELSE
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDEX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R3,NFINI,INDEX,STAT_R3,IERR )
           TOT = TOT + NFINI
!debug     print *, myproc, tot,nfini,index(1),index(2)
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                    DO K=1,DG_HERE%DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(IRECVLOC(I,J),K,IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DG_HERE%DOFH
                      NCOUNT = NCOUNT+1
                      VEC2(IRECVLOC(I,J),K,IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,dg_here%DOFH
                      NCOUNT = NCOUNT+1
                      VEC3(IRECVLOC(I,J),K,IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
        GOTO 9998
      ENDIF
! 
 9998 continue
      RETURN
      END SUBROUTINE

      SUBROUTINE UPDATER_ELEM_MOD3(dg_here, VEC1, VEC2, VEC3, NMSG )
!
!--------------------------------------------------------------------------
!  Update 1, 2, or 3 Real Arrays's Ghost Cells using asynchronous
!  and persistent message-passing.
!
!  vjp  8/06/1999
!--------------------------------------------------------------------------
! 
        use dg
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
!
      type (dg_type) :: dg_here
      INTEGER,  INTENT(IN) ::  NMSG
      REAL(SZ), INTENT(INOUT) ::  VEC1(:,:),VEC2(:,:),VEC3(:,:)
      INTEGER :: N,I,J,JJ,K,NCOUNT,NFINI,TOT
!
                             !..Pack 1, 2, or 3 Messages
!      DO JJ=1,NEIGHPROC_S
!         J = MOD(JJ-1+MYPROC,NEIGHPROC_S)+1
      DO J=1,NEIGHPROC_S
         NCOUNT = 0
         DO I=1,NELEMSEND(J)
            DO K=1,DG_HERE%DOFH
              NCOUNT = NCOUNT+1
              SENDBUF(NCOUNT,J)=VEC1(K,ISENDLOC(I,J))
            ENDDO
         ENDDO
         IF (NMSG.GT.1) THEN
           DO I=1,NELEMSEND(J)
             DO K=1,DG_HERE%DOFH
               NCOUNT = NCOUNT+1
               SENDBUF(NCOUNT,J)=VEC2(K,ISENDLOC(I,J))
             ENDDO
           ENDDO
         ENDIF
         IF (NMSG.GT.2) THEN
           DO I=1,NELEMSEND(J)
             DO K=1,DG_HERE%DOFH
               NCOUNT = NCOUNT+1
               SENDBUF(NCOUNT,J)=VEC3(K,ISENDLOC(I,J))
             ENDDO
           ENDDO
         ENDIF

         ! Start sending a message
         IF (NMSG.EQ.1) THEN
            CALL MPI_START ( REQ_R1(J+NEIGHPROC_R), IERR )
         ELSEIF (NMSG.EQ.2) THEN
            CALL MPI_START ( REQ_R2(J+NEIGHPROC_R), IERR )
         ELSE
            CALL MPI_START ( REQ_R3(J+NEIGHPROC_R), IERR )
         ENDIF
      ENDDO
!                    
              ! Send/receive messages to/from all neighbors
!
      IF (NMSG.EQ.1) THEN
        CALL MPI_STARTALL ( RDIM-NEIGHPROC_S, REQ_R1, IERR )
      ELSEIF (NMSG.EQ.2) THEN
        CALL MPI_STARTALL ( RDIM-NEIGHPROC_S, REQ_R2, IERR )
      ELSE
        CALL MPI_STARTALL ( RDIM-NEIGHPROC_S, REQ_R3, IERR )
      ENDIF
!
              !..Unpack Received messages as they arrive     
!
      IF (NMSG.EQ.1) THEN   
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDEX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R1,NFINI,INDEX,STAT_R1,IERR )
           TOT = TOT + NFINI
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                    DO K=1,DG_HERE%DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(K,IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
        GOTO 9998
      ELSEIF (NMSG.EQ.2) THEN
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDEX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R2,NFINI,INDEX,STAT_R2,IERR )
           TOT = TOT + NFINI
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                    DO K=1,DG_HERE%DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(K,IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DG_HERE%DOFH
                      NCOUNT = NCOUNT+1
                      VEC2(K,IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
        GOTO 9998
      ELSE
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDEX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_R3,NFINI,INDEX,STAT_R3,IERR )
           TOT = TOT + NFINI
!debug     print *, myproc, tot,nfini,index(1),index(2)
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                    DO K=1,DG_HERE%DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(K,IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DG_HERE%DOFH
                      NCOUNT = NCOUNT+1
                      VEC2(K,IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DG_HERE%DOFH
                      NCOUNT = NCOUNT+1
                      VEC3(K,IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
        GOTO 9998
      ENDIF
! 
 9998 continue
      RETURN
      END SUBROUTINE


#if 0

      SUBROUTINE UPDATER3D( VEC )
!
!--------------------------------------------------------------------------
!  Update 1 Three-dimensional Real Arrays's Ghost Cells using asynchronous
!  and persistent message-passing.
!
!  tjc  6/24/2002
!--------------------------------------------------------------------------
! 
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
!
      REAL(SZ), INTENT(INOUT) ::  VEC(MNP,MNODES)
      INTEGER :: N,I,J,K,NCOUNT,NFINI,TOT
!$$$C
!$$$                             !..Pack Messages
!$$$      DO J=1,NEIGHPROC
!$$$         NCOUNT = 0
!$$$         DO I=1,NELEMSEND(J)
!$$$            DO K=1,MNODES
!$$$               NCOUNT = NCOUNT+1
!$$$               SENDBUF(NCOUNT,J)=VEC(ISENDLOC(I,J),K)
!$$$            ENDDO
!$$$         ENDDO
!$$$      ENDDO
!$$$C                    
!$$$              ! Send/receive messages to/from all neighbors
!$$$C
!$$$      CALL MPI_STARTALL ( RDIM, REQ_R3D, IERR )
!$$$C
!$$$              !..Unpack Received messages as they arrive     
!$$$C
!$$$      TOT = 0
!$$$      DO WHILE (TOT.LT.RDIM)
!$$$         DO N=1, RDIM
!$$$            INDEX(N) = 0
!$$$         ENDDO
!$$$         CALL MPI_WAITSOME( RDIM,REQ_R3D,NFINI,INDEX,STAT_R3D,IERR )
!$$$         TOT = TOT + NFINI
!$$$         DO N=1, NFINI
!$$$            IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
!$$$              IF (INDEX(N).LE.NEIGHPROC) THEN
!$$$                J = INDEX(N)
!$$$                NCOUNT = 0
!$$$                DO I=1,NELEMRECV(J)
!$$$                   DO K=1,MNODES
!$$$                      NCOUNT = NCOUNT+1
!$$$                      VEC(IRECVLOC(I,J),K) = RECVBUF(NCOUNT,J)
!$$$                   ENDDO
!$$$                ENDDO
!$$$              ENDIF
!$$$            ENDIF
!$$$         ENDDO
!$$$      ENDDO
!$$$C 
      STOP 'NOT UPDATED'
      RETURN
      END SUBROUTINE

      SUBROUTINE UPDATEC3D( VEC )
!
!--------------------------------------------------------------------------
!  Update 1 Three-dimensional Complex Arrays's Ghost Cells using asynchronous
!  and persistent message-passing.
!
!  tjc  6/24/2002
!--------------------------------------------------------------------------
! 
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
!
      COMPLEX, INTENT(INOUT) ::  VEC(MNP,MNODES)
      INTEGER :: N,I,J,K,NCOUNT,NFINI,TOT
!$$$C
!$$$                             !..Pack Messages
!$$$      DO J=1,NEIGHPROC
!$$$         NCOUNT = 0
!$$$         DO I=1,NELEMSEND(J)
!$$$            DO K=1,MNODES
!$$$               NCOUNT = NCOUNT+1
!$$$               SENDBUF(NCOUNT,J)=REAL(VEC(ISENDLOC(I,J),K))
!$$$               NCOUNT = NCOUNT+1
!$$$               SENDBUF(NCOUNT,J)=AIMAG(VEC(ISENDLOC(I,J),K))
!$$$            ENDDO
!$$$         ENDDO
!$$$      ENDDO
!$$$C                    
!$$$              ! Send/receive messages to/from all neighbors
!$$$C
!$$$      CALL MPI_STARTALL ( RDIM, REQ_C3D, IERR )
!$$$C
!$$$              !..Unpack Received messages as they arrive     
!$$$C
!$$$      TOT = 0
!$$$      DO WHILE (TOT.LT.RDIM)
!$$$         DO N=1, RDIM
!$$$            INDEX(N) = 0
!$$$         ENDDO
!$$$         CALL MPI_WAITSOME( RDIM,REQ_C3D,NFINI,INDEX,STAT_C3D,IERR )
!$$$         TOT = TOT + NFINI
!$$$         DO N=1, NFINI
!$$$            IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
!$$$              IF (INDEX(N).LE.NEIGHPROC) THEN
!$$$                J = INDEX(N)
!$$$                NCOUNT = 0
!$$$                DO I=1,NELEMRECV(J)
!$$$                   DO K=1,MNODES
!$$$                      VEC(IRECVLOC(I,J),K) = 
!$$$     &                   CMPLX(RECVBUF(NCOUNT+1,J),RECVBUF(NCOUNT+2,J))
!$$$                      NCOUNT = NCOUNT+2
!$$$                   ENDDO
!$$$                ENDDO
!$$$              ENDIF
!$$$            ENDIF
!$$$         ENDDO
!$$$      ENDDO
!$$$C 
      STOP 'NOT UPDATED'
      RETURN
      END SUBROUTINE
#endif

      SUBROUTINE MSG_BLOCKSYNC_START(s)
        use sizes
        IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      
      type (sizes_type) :: s

      INTEGER,PARAMETER :: BLOCK = 10
      INTEGER :: SRC,DUM,STAT(MPI_STATUS_SIZE)

      SRC = s%MYPROC - BLOCK
      IF(SRC.LT.0) RETURN

      CALL MPI_RECV(DUM,1,MPI_INTEGER,SRC,TAG+1,MPI_COMM_WORLD,STAT,IERR)

      END SUBROUTINE 

      SUBROUTINE MSG_BLOCKSYNC_FINISH(s)
      use sizes
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      type (sizes_type) :: s

      INTEGER,PARAMETER :: BLOCK = 10
      INTEGER :: DST,DUM,STAT(MPI_STATUS_SIZE)

      DST = S%MYPROC + BLOCK
      IF(DST.GE.S%MNPROC) RETURN

      DUM = 1

      CALL MPI_SEND(DUM,1,MPI_INTEGER,DST,TAG+1,MPI_COMM_WORLD,STAT,IERR)

      END SUBROUTINE 

      SUBROUTINE para_sum( sum_this )
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif


      Real(SZ), INTENT(INOUT) ::sum_this
      Real(SZ) sum_of_this 
      Integer kount       

      sum_of_this = 0.D0

      kount=1

      call MPI_ALLREDUCE( sum_this, sum_of_this, kount,&
     REALTYPE, MPI_SUM, COMM, ierr)

      sum_this = sum_of_this 

      END SUBROUTINE para_sum

      SUBROUTINE para_max( max_this )
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif

      Real(SZ), INTENT(INOUT) :: max_this
      Real(SZ) max_of_this     
      Integer kount            

      max_of_this = -100.D0

      kount=1

      call MPI_ALLREDUCE( max_this, max_of_this, kount,&
     REALTYPE, MPI_MAX, COMM, ierr)

      max_this = max_of_this

      END SUBROUTINE para_max

      SUBROUTINE para_min( min_this )
      
      IMPLICIT NONE

#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif


      Real(SZ), INTENT(INOUT) :: min_this
      Real(SZ) min_of_this    
      Integer kount           

      min_of_this = 100.D0

      kount=1

      call MPI_ALLREDUCE( min_this, min_of_this, kount,&
     REALTYPE, MPI_MIN, COMM, ierr)

      min_this = min_of_this

      END SUBROUTINE para_min

      END MODULE MESSENGER_ELEM
       
