C******************************************************************************
C  last changes in this file VERSION 12.sb01                                  *
C  S.Bunya changed this file a bit. 07/13/2005                                *
C  S.Bunya changed this file a bit. 01/04/2007                                *
C                                                                             *
C     01-10-2011 - cem - adapted for p_enrichment and multicomponent          * C****************************************************************************** 
C 
      MODULE MESSENGER_ELEM
C
      USE SIZES
      USE GLOBAL!, ONLY: C3D
      USE DG, ONLY: DOFH
      USE DIFF45_41,ONLY : MNODES


#ifdef HAVE_MPI_MOD
      use mpi  
      IMPLICIT NONE
#else
      IMPLICIT NONE
#endif
C
C--------------------------------------------------------------------------
C  This module supplies the MPI Message-Passing Interface
C
C  Uses asynchronous and persistent communication with buffer packing
C  as performance enhancements for "cluster" architectures.
C
C  vjp  8/29/1999
C--------------------------------------------------------------------------
C
C

C
C  Message-Passing Array space
C
Csb-PDG1
      PUBLIC

      INTEGER, SAVE :: MPI_COMM_ADCIRC
      INTEGER, SAVE ::  COMM_COMP     
      INTEGER, SAVE ::  GROUP_WORLD, GROUP_COMP

      INTEGER,SAVE,PRIVATE :: REALTYPE, DBLETYPE, COMM   
      INTEGER,SAVE,PRIVATE ::  NEIGHPROC_R, NEIGHPROC_S, RDIM, IERR
      INTEGER,SAVE,PRIVATE ::  TAG = 101
      LOGICAL,SAVE,ALLOCATABLE :: RESELEM(:)
C
      INTEGER, PRIVATE, ALLOCATABLE ::IPROC_R(:),IPROC_S(:),NELEMLOC(:),
     &    NELEMSEND(:), NELEMRECV(:),ISENDLOC(:,:), IBELONGTO(:),
     &    IRECVLOC(:,:), ISENDBUF(:,:), IRECVBUF(:,:)
C
      INTEGER, PRIVATE, ALLOCATABLE :: REQ_I1(:), REQ_I2(:)
      INTEGER, PRIVATE, ALLOCATABLE :: STAT_I1(:,:), STAT_I2(:,:)
      INTEGER, PRIVATE, ALLOCATABLE :: REQ_R1(:), REQ_R2(:), REQ_R3(:),
     &                                 REQ_LZ(:) 
      INTEGER, PRIVATE, ALLOCATABLE :: STAT_R1(:,:), STAT_R2(:,:), 
     &                                 STAT_R3(:,:), STAT_LZ(:,:)
      INTEGER, PRIVATE, ALLOCATABLE :: REQ_R3D(:), STAT_R3D(:,:)
      INTEGER, PRIVATE, ALLOCATABLE :: REQ_C3D(:), STAT_C3D(:,:)
      INTEGER, PRIVATE, ALLOCATABLE :: INDEX(:)
      REAL(SZ), PRIVATE,ALLOCATABLE :: SENDBUF(:,:), RECVBUF(:,:)
C--
C

C---------------------end of data declarations--------------------------------C


      CONTAINS


      SUBROUTINE MSG_TYPES_ELEM()
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
c
#ifdef  REAL4
      REALTYPE = MPI_REAL
      DBLETYPE = MPI_DOUBLE_PRECISION
#else
      REALTYPE = MPI_DOUBLE_PRECISION
      DBLETYPE = MPI_DOUBLE_PRECISION
#endif
c
      RETURN
      END  SUBROUTINE


      SUBROUTINE MESSAGE_INIT ()
C--------------------------------------------------------------------------
C  Routine performs following steps:
C   (1)  initialize MPI, 
C   (2)  get number of processors,
C   (3)  get MPI rank of processor 
C  vjp  8/06/1999
C--------------------------------------------------------------------------
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
CND
      Integer I
      INTEGER, ALLOCATABLE :: RANKS(:)

      CALL MPI_INIT(IERR)                               ! Initialize MPI
      CALL MPI_COMM_SIZE (MPI_COMM_WORLD,MNPROC,IERR)   ! Get number of procs
      CALL MPI_COMM_RANK (MPI_COMM_WORLD,MYPROC,IERR)   ! Get MPI rank
CND
      ALLOCATE(RANKS(MNPROC+1))
      DO I=1,MNPROC
         RANKS(I)=I-1
      ENDDO
      CALL MPI_COMM_GROUP(MPI_COMM_WORLD,GROUP_WORLD,IERR)
      CALL MPI_GROUP_INCL(GROUP_WORLD,MNPROC,RANKS,GROUP_COMP,IERR)
      CALL MPI_COMM_CREATE(MPI_COMM_WORLD,GROUP_COMP,COMM_COMP,IERR)
      DEALLOCATE(RANKS)
      COMM = COMM_COMP
      RETURN
      END SUBROUTINE

c$$$      IMPLICIT NONE
c$$$#ifndef HAVE_MPI_MOD
c$$$      include 'mpif.h'
c$$$#endif
c$$$      CALL MPI_INIT(IERR)                               ! Initialize MPI
c$$$      CALL MPI_COMM_SIZE (MPI_COMM_WORLD,MNPROC,IERR)   ! Get number of procs
c$$$      CALL MPI_COMM_RANK (MPI_COMM_WORLD,MYPROC,IERR)   ! Get MPI rank
c$$$      RETURN
c$$$      END SUBROUTINE

      SUBROUTINE ErrorElevSum( ErrorElevExceeded )
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      INTEGER ErrorElevExceeded !=1 if this subdomain has exceeded warning elev
      INTEGER SumErrorElevExceeded !sum total of all flags from all subdomains
      INTEGER kount             ! to avoid compiler bug on certain platforms
      
      SumErrorElevExceeded = 0
      kount=1
      call MPI_ALLREDUCE( ErrorElevExceeded, SumErrorElevExceeded, kount,
     &     MPI_INTEGER, MPI_SUM, MPI_COMM_world, ierr)
      ErrorElevExceeded = SumErrorElevExceeded
      END SUBROUTINE ErrorElevSum

 
      SUBROUTINE MESSAGE_FINI ()
C--------------------------------------------------------------------------
C  Shutdown MPI library.
C--------------------------------------------------------------------------
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      INTEGER IERR,I
C
      CALL MPI_FINALIZE(IERR)
      IF (MYPROC.EQ.0)  
     &  print *, "MPI terminated with Status = ",IERR      
      RETURN
      END SUBROUTINE

 
      SUBROUTINE MESSAGE_ABORT ()
C--------------------------------------------------------------------------
C  Shutdown MPI library.
C--------------------------------------------------------------------------
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      INTEGER IERR,I
C
      CALL MPI_ABORT(MPI_COMM_WORLD,IERR)
      IF (MYPROC.EQ.0)  
     &  print *, "MPI Aborted with Status = ",IERR      
      RETURN
      END SUBROUTINE



      SUBROUTINE MSG_TABLE_ELEM () 
C
C--------------------------------------------------------------------------
C  Routine preforms following steps:
C
C   (1) Read Message-Passing Information from file "DG.18"
C   (2) Determine resident nodes: RESNODE(I) is true  if I is resident node
C   (3) Determine ghost nodes:    RESNODE(I) is false if I is ghost node    
C   (4) Determine number of neighbor subdomains
C   (5) MPI rank of each neighbor and number of ghosts nodes to receive
C   (6) Read Message-Passing Receive List
C   (7) MPI rank of each neighbor and number of ghosts nodes to send
C   (8) Read Message-Passing Send List
C
C  vjp  8/06/1999
C--------------------------------------------------------------------------
C
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      INTEGER :: IDPROC,NLOCAL,I,J,JJ,NEIGH
C
      OPEN(18,FILE=DIRNAME(1:LNAME)//'/'//'DG.18')
C
      READ(18,3010) IDPROC,NLOCAL
C
      ALLOCATE ( NELEMLOC(NLOCAL) )
C
      READ(18,1130) (NELEMLOC(I), I=1,NLOCAL)
C
      ALLOCATE ( IBELONGTO(MNE),RESELEM(MNE) )
C
      DO I=1,MNE
         IBELONGTO(I) = 0
      ENDDO
      DO I=1,NLOCAL
         IBELONGTO(NELEMLOC(I)) = IDPROC + 1
      ENDDO
      DO I=1, MNE
         IF (IBELONGTO(I)-1.EQ.MYPROC) THEN
           RESELEM(I) = .TRUE.
         ELSE 
           RESELEM(I) = .FALSE.
         ENDIF
      ENDDO
C
      READ(18,3010) NEIGHPROC_R,NEIGHPROC_S
C
      RDIM = NEIGHPROC_R + NEIGHPROC_S
      ALLOCATE( INDEX(RDIM) )
C
      ALLOCATE( IPROC_R(NEIGHPROC_R),NELEMRECV(NEIGHPROC_R) )
      ALLOCATE( IRECVLOC(MNE,NEIGHPROC_R) )
C
      DO JJ=1,NEIGHPROC_R
         J = MOD(JJ-1+MYPROC,NEIGHPROC_R)+1
         READ(18,3010) IPROC_R(J),NELEMRECV(J)
         READ(18,1130) (IRECVLOC(I,J), I=1,NELEMRECV(J))
      ENDDO
C
      ALLOCATE( IPROC_S(NEIGHPROC_S),NELEMSEND(NEIGHPROC_S) )
      ALLOCATE( ISENDLOC(MNE,NEIGHPROC_S) )
C
      DO JJ=1,NEIGHPROC_S
         J = MOD(JJ-1+MYPROC,NEIGHPROC_S)+1
         READ(18,3010) IPROC_S(J),NELEMSEND(J)
         READ(18,1130) (ISENDLOC(I,J), I=1,NELEMSEND(J))
      ENDDO
C
      CLOSE(18)
      RETURN
C
1130  FORMAT(8X,9I8)
3010  FORMAT(8X,2I8)
      END SUBROUTINE


      SUBROUTINE MESSAGE_START_ELEM ()
C
C--------------------------------------------------------------------------
C  Routine preforms following steps:
C   (1)  allocate message-passing space
C   (2)  setup MPI data structures for "persistent" message-passing.
C
C  vjp  8/06/1999
C--------------------------------------------------------------------------
C
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      INTEGER :: J,NCOMMELEM_S,NCOMMELEM_R
C
      NCOMMELEM_R = 0
      NCOMMELEM_S = 0
      DO J=1,NEIGHPROC_R
        NCOMMELEM_R = NCOMMELEM_R + NELEMRECV(J)
      ENDDO
      DO J=1,NEIGHPROC_S
        NCOMMELEM_S = NCOMMELEM_S + NELEMSEND(J)
      ENDDO
C
      ALLOCATE ( ISENDBUF(NCOMMELEM_S*2,NEIGHPROC_S) )
      ALLOCATE ( IRECVBUF(NCOMMELEM_R*2,NEIGHPROC_R) )
C
      IF (C3D) THEN
C         ALLOCATE ( SENDBUF(2*MNP*MNODES,NEIGHPROC) )
C         ALLOCATE ( RECVBUF(2*MNP*MNODES,NEIGHPROC) )
        STOP 'NOT UPDATED'
      ELSE
         ALLOCATE ( SENDBUF(NCOMMELEM_S*DOFH*4,NEIGHPROC_S) )
         ALLOCATE ( RECVBUF(NCOMMELEM_R*DOFH*4,NEIGHPROC_R) )
      ENDIF
C
      ALLOCATE ( REQ_I1(RDIM),REQ_I2(RDIM) )
      ALLOCATE ( REQ_R1(RDIM),REQ_R2(RDIM),REQ_R3(RDIM),REQ_LZ(RDIM) )
C
      ALLOCATE ( STAT_I1(MPI_STATUS_SIZE,RDIM),       
     &           STAT_I2(MPI_STATUS_SIZE,RDIM) )

      ALLOCATE ( STAT_R1(MPI_STATUS_SIZE,RDIM),       
     &           STAT_R2(MPI_STATUS_SIZE,RDIM),
     &           STAT_R3(MPI_STATUS_SIZE,RDIM),
     &           STAT_LZ(MPI_STATUS_SIZE,RDIM) )
C
      IF (C3D) THEN
C         ALLOCATE ( REQ_R3D(RDIM) )
C         ALLOCATE ( STAT_R3D(MPI_STATUS_SIZE,RDIM) )
C         ALLOCATE ( REQ_C3D(RDIM) )
C         ALLOCATE ( STAT_C3D(MPI_STATUS_SIZE,RDIM) )
        STOP 'NOT UPDATED'
      ENDIF
C
             !  Setup persistent structures for integer arrays
C
      DO J=1,NEIGHPROC_R
         CALL MPI_RECV_INIT ( IRECVBUF(1,J), NELEMRECV(J), 
     &     MPI_INTEGER,IPROC_R(J), TAG, MPI_COMM_WORLD,
     &     REQ_I1(J),IERR)
      ENDDO
      DO J=1,NEIGHPROC_S
         CALL MPI_SEND_INIT ( ISENDBUF(1,J), NELEMSEND(J), 
     &    MPI_INTEGER,IPROC_S(J), TAG, MPI_COMM_WORLD,
     &    REQ_I1(J+NEIGHPROC_R),IERR )
      ENDDO
C
C
      DO J=1,NEIGHPROC_R
         CALL MPI_RECV_INIT ( IRECVBUF(1,J), 2*NELEMRECV(J), 
     &     MPI_INTEGER,IPROC_R(J), TAG, MPI_COMM_WORLD,
     &     REQ_I2(J),IERR)
      ENDDO
      DO J=1,NEIGHPROC_S
         CALL MPI_SEND_INIT ( ISENDBUF(1,J), 2*NELEMSEND(J), 
     &    MPI_INTEGER,IPROC_S(J), TAG, MPI_COMM_WORLD,
     &    REQ_I2(J+NEIGHPROC_R),IERR )
      ENDDO
C
            !  Setup persistent structures for real arrays
C
      DO J=1,NEIGHPROC_R
         CALL MPI_RECV_INIT ( RECVBUF(1,J), DOFH*NELEMRECV(J), 
     &     REALTYPE,IPROC_R(J), TAG, MPI_COMM_WORLD,
     &     REQ_R1(J),IERR)
      ENDDO
      DO J=1,NEIGHPROC_S
         CALL MPI_SEND_INIT ( SENDBUF(1,J), DOFH*NELEMSEND(J), 
     &     REALTYPE,IPROC_S(J), TAG, MPI_COMM_WORLD,
     &     REQ_R1(J+NEIGHPROC_R),IERR)
      ENDDO
C
      DO J=1,NEIGHPROC_R
         CALL MPI_RECV_INIT ( RECVBUF(1,J), 2*DOFH*NELEMRECV(J), 
     &     REALTYPE,IPROC_R(J), TAG, MPI_COMM_WORLD,
     &     REQ_R2(J),IERR)
      ENDDO
      DO J=1,NEIGHPROC_S
         CALL MPI_SEND_INIT ( SENDBUF(1,J), 2*DOFH*NELEMSEND(J), 
     &     REALTYPE,IPROC_S(J), TAG, MPI_COMM_WORLD,
     &     REQ_R2(J+NEIGHPROC_R),IERR)
      ENDDO
C
      DO J=1,NEIGHPROC_R
         CALL MPI_RECV_INIT ( RECVBUF(1,J), 3*DOFH*NELEMRECV(J), 
     &     REALTYPE,IPROC_R(J), TAG, MPI_COMM_WORLD,
     &     REQ_R3(J),IERR)
      ENDDO
      DO J=1,NEIGHPROC_S
         CALL MPI_SEND_INIT ( SENDBUF(1,J), 3*DOFH*NELEMSEND(J), 
     &     REALTYPE,IPROC_S(J), TAG, MPI_COMM_WORLD,
     &     REQ_R3(J+NEIGHPROC_R),IERR)
      ENDDO
C
      DO J=1,NEIGHPROC_R
         CALL MPI_RECV_INIT ( RECVBUF(1,J), 4*DOFH*NELEMRECV(J), 
     &     REALTYPE,IPROC_R(J), TAG, MPI_COMM_WORLD,
     &     REQ_LZ(J),IERR)
      ENDDO
      DO J=1,NEIGHPROC_S
         CALL MPI_SEND_INIT ( SENDBUF(1,J), 4*DOFH*NELEMSEND(J), 
     &     REALTYPE,IPROC_S(J), TAG, MPI_COMM_WORLD,
     &     REQ_LZ(J+NEIGHPROC_R),IERR)
      ENDDO
C
      IF (C3D) THEN
C         DO J=1,NEIGHPROC  
C            CALL MPI_RECV_INIT ( RECVBUF(1,J), MNODES*NNODRECV(J), 
C     &        REALTYPE,IPROC(J), TAG, MPI_COMM_WORLD,
C     &        REQ_R3D(J),IERR)
C         ENDDO
C         DO J=1,NEIGHPROC  
C            CALL MPI_SEND_INIT ( SENDBUF(1,J), MNODES*NNODSEND(J), 
C     &        REALTYPE,IPROC(J), TAG, MPI_COMM_WORLD,
C     &        REQ_R3D(J+NEIGHPROC),IERR)
C         ENDDO
C         DO J=1,NEIGHPROC  
C            CALL MPI_RECV_INIT ( RECVBUF(1,J), 2*MNODES*NNODRECV(J), 
C     &        REALTYPE,IPROC(J), TAG, MPI_COMM_WORLD,
C     &        REQ_C3D(J),IERR)
C         ENDDO
C         DO J=1,NEIGHPROC  
C            CALL MPI_SEND_INIT ( SENDBUF(1,J), 2*MNODES*NNODSEND(J), 
C     &        REALTYPE,IPROC(J), TAG, MPI_COMM_WORLD,
C     &        REQ_C3D(J+NEIGHPROC),IERR)
C         ENDDO
        STOP 'NOT UPDATED'
      ENDIF
C
      RETURN
      END SUBROUTINE


      SUBROUTINE UPDATEI_ELEM( IVEC1, IVEC2, NMSG )
C
C--------------------------------------------------------------------------
C  Update 1 or 2 Integer Arrays's Ghost Cells using asynchronous
C  and persistent message-passing.
C
C  vjp  8/06/1999
C--------------------------------------------------------------------------
C 
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
C
      INTEGER,   INTENT(IN) :: NMSG
      Real(sz),   INTENT(INOUT) :: IVEC1(:),IVEC2(:)
      INTEGER :: N,I,J,NCOUNT,NFINI,TOT
C
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
C                     
                          ! Send/receive messages to/from all neighbors
      IF (NMSG.EQ.1) THEN
        CALL MPI_STARTALL ( RDIM, REQ_R1, IERR )
      ELSE
        CALL MPI_STARTALL ( RDIM, REQ_R2, IERR )
      ENDIF
C
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
C 
 999  continue
      RETURN
      END SUBROUTINE


      SUBROUTINE UPDATER_ELEM( VEC1, VEC2, VEC3, IRK, NMSG )
C
C--------------------------------------------------------------------------
C  Update 1, 2, or 3 Real Arrays's Ghost Cells using asynchronous
C  and persistent message-passing.
C
C  vjp  8/06/1999
C--------------------------------------------------------------------------
C 
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
C
      INTEGER,  INTENT(IN) ::  IRK, NMSG
      REAL(SZ), INTENT(INOUT) ::  VEC1(:,:,:),VEC2(:,:,:),VEC3(:,:,:)
      INTEGER :: N,I,J,K,NCOUNT,NFINI,TOT
C
                             !..Pack 1, 2, or 3 Messages
      DO J=1,NEIGHPROC_S
         NCOUNT = 0
         DO I=1,NELEMSEND(J)
            DO K=1,DOFH
              NCOUNT = NCOUNT+1
              SENDBUF(NCOUNT,J)=VEC1(K,ISENDLOC(I,J),IRK)
            ENDDO
         ENDDO
         IF (NMSG.GT.1) THEN
           DO I=1,NELEMSEND(J)
             DO K=1,DOFH
               NCOUNT = NCOUNT+1
               SENDBUF(NCOUNT,J)=VEC2(K,ISENDLOC(I,J),IRK)
             ENDDO
           ENDDO
         ENDIF
         IF (NMSG.GT.2) THEN
           DO I=1,NELEMSEND(J)
             DO K=1,DOFH
               NCOUNT = NCOUNT+1
               SENDBUF(NCOUNT,J)=VEC3(K,ISENDLOC(I,J),IRK)
             ENDDO
           ENDDO
         ENDIF
      ENDDO
C                    
              ! Send/receive messages to/from all neighbors
C
      IF (NMSG.EQ.1) THEN
        CALL MPI_STARTALL ( RDIM, REQ_R1, IERR )
      ELSEIF (NMSG.EQ.2) THEN
        CALL MPI_STARTALL ( RDIM, REQ_R2, IERR )
      ELSE
        CALL MPI_STARTALL ( RDIM, REQ_R3, IERR )
      ENDIF
              !..Unpack Received messages as they arrive     
C
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
                    DO K=1,DOFH
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
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(K,IRECVLOC(I,J),IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
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
cdebug     print *, myproc, tot,nfini,index(1),index(2)
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(K,IRECVLOC(I,J),IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC2(K,IRECVLOC(I,J),IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
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
C 
 999  continue
      RETURN
      END SUBROUTINE


      SUBROUTINE UPDATELZ_ELEM( LZ )
C
C--------------------------------------------------------------------------
C  Update LZ Real Arrays's Ghost Cells using asynchronous
C  and persistent message-passing.
C
C  sb  1/04/2007
C--------------------------------------------------------------------------
C 
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
C
      REAL(SZ), INTENT(INOUT) ::  LZ(:,:,:,:)
      INTEGER :: N,I,J,K,NCOUNT,NFINI,TOT
C
                             !..Pack 1, 2, or 3 Messages
      DO J=1,NEIGHPROC_S
         NCOUNT = 0
         DO I=1,NELEMSEND(J)
            DO K=1,DOFH
              SENDBUF(NCOUNT+1,J)=LZ(K,1,1,ISENDLOC(I,J))
              SENDBUF(NCOUNT+2,J)=LZ(K,1,2,ISENDLOC(I,J))
              SENDBUF(NCOUNT+3,J)=LZ(K,2,1,ISENDLOC(I,J))
              SENDBUF(NCOUNT+4,J)=LZ(K,2,2,ISENDLOC(I,J))
              NCOUNT = NCOUNT+4
            ENDDO
         ENDDO
      ENDDO
C                    
              ! Send/receive messages to/from all neighbors
C
      CALL MPI_STARTALL ( RDIM, REQ_LZ, IERR )

              !..Unpack Received messages as they arrive     
C
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
                DO K=1,DOFH
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

      SUBROUTINE UPDATEMZ_ELEM( MZ )
C
C--------------------------------------------------------------------------
C  Update MZ Real Arrays's Ghost Cells using asynchronous
C  and persistent message-passing.
C
C  sb  1/04/2007
C--------------------------------------------------------------------------
C 
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
C
      REAL(SZ), INTENT(INOUT) ::  MZ(:,:,:,:)
      INTEGER :: N,I,J,K,NCOUNT,NFINI,TOT
C
                             !..Pack 1, 2, or 3 Messages
      DO J=1,NEIGHPROC_S
         NCOUNT = 0
         DO I=1,NELEMSEND(J)
            DO K=1,DOFH
              SENDBUF(NCOUNT+1,J)=MZ(K,1,1,ISENDLOC(I,J))
              SENDBUF(NCOUNT+2,J)=MZ(K,2,1,ISENDLOC(I,J))
              NCOUNT = NCOUNT+2
            ENDDO
         ENDDO
      ENDDO
C                    
              ! Send/receive messages to/from all neighbors
C
      CALL MPI_STARTALL ( RDIM, REQ_LZ, IERR )

              !..Unpack Received messages as they arrive     
C
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
                DO K=1,DOFH
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


      SUBROUTINE UPDATER_ELEM_MOD( VEC1, VEC2, VEC3, IRK, NMSG )
C
C--------------------------------------------------------------------------
C  Update 1, 2, or 3 Real Arrays's Ghost Cells using asynchronous
C  and persistent message-passing.
C
C  vjp  8/06/1999
C--------------------------------------------------------------------------
C 
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
C
      INTEGER,  INTENT(IN) ::  IRK, NMSG
      REAL(SZ), INTENT(INOUT) ::  VEC1(:,:,:),VEC2(:,:,:),VEC3(:,:,:)
      INTEGER :: N,I,J,JJ,K,NCOUNT,NFINI,TOT
C
                             !..Pack 1, 2, or 3 Messages
C      DO JJ=1,NEIGHPROC_S
C         J = MOD(JJ-1+MYPROC,NEIGHPROC_S)+1
      DO J=1,NEIGHPROC_S
         NCOUNT = 0
         DO I=1,NELEMSEND(J)
            DO K=1,DOFH
              NCOUNT = NCOUNT+1
              SENDBUF(NCOUNT,J)=VEC1(K,ISENDLOC(I,J),IRK)
            ENDDO
         ENDDO
         IF (NMSG.GT.1) THEN
           DO I=1,NELEMSEND(J)
             DO K=1,DOFH
               NCOUNT = NCOUNT+1
               SENDBUF(NCOUNT,J)=VEC2(K,ISENDLOC(I,J),IRK)
             ENDDO
           ENDDO
         ENDIF
         IF (NMSG.GT.2) THEN
           DO I=1,NELEMSEND(J)
             DO K=1,DOFH
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
C                    
              ! Send/receive messages to/from all neighbors
C
      IF (NMSG.EQ.1) THEN
        CALL MPI_STARTALL ( RDIM-NEIGHPROC_S, REQ_R1, IERR )
      ELSEIF (NMSG.EQ.2) THEN
        CALL MPI_STARTALL ( RDIM-NEIGHPROC_S, REQ_R2, IERR )
      ELSE
        CALL MPI_STARTALL ( RDIM-NEIGHPROC_S, REQ_R3, IERR )
      ENDIF
C
              !..Unpack Received messages as they arrive     
C
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
                    DO K=1,DOFH
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
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(K,IRECVLOC(I,J),IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
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
cdebug     print *, myproc, tot,nfini,index(1),index(2)
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(K,IRECVLOC(I,J),IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC2(K,IRECVLOC(I,J),IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
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
C 
 999  continue
      RETURN
      END SUBROUTINE

      SUBROUTINE UPDATER_ELEM_MOD2( VEC1, VEC2, VEC3, IRK, NMSG )
C
C--------------------------------------------------------------------------
C  Update 1, 2, or 3 Real Arrays's Ghost Cells using asynchronous
C  and persistent message-passing.
C
C  vjp  8/06/1999
C--------------------------------------------------------------------------
C 
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
C
      INTEGER,  INTENT(IN) ::  IRK, NMSG
      REAL(SZ), INTENT(INOUT) ::  VEC1(:,:,:),VEC2(:,:,:),VEC3(:,:,:)
      INTEGER :: N,I,J,JJ,K,NCOUNT,NFINI,TOT
C
                             !..Pack 1, 2, or 3 Messages
C      DO JJ=1,NEIGHPROC_S
C         J = MOD(JJ-1+MYPROC,NEIGHPROC_S)+1
      DO J=1,NEIGHPROC_S
         NCOUNT = 0
         DO I=1,NELEMSEND(J)
            DO K=1,DOFH
              NCOUNT = NCOUNT+1
              SENDBUF(NCOUNT,J)=VEC1(ISENDLOC(I,J),K,IRK)
            ENDDO
         ENDDO
         IF (NMSG.GT.1) THEN
           DO I=1,NELEMSEND(J)
             DO K=1,DOFH
               NCOUNT = NCOUNT+1
               SENDBUF(NCOUNT,J)=VEC2(ISENDLOC(I,J),K,IRK)
             ENDDO
           ENDDO
         ENDIF
         IF (NMSG.GT.2) THEN
           DO I=1,NELEMSEND(J)
             DO K=1,DOFH
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
C                    
              ! Send/receive messages to/from all neighbors
C
      IF (NMSG.EQ.1) THEN
        CALL MPI_STARTALL ( RDIM-NEIGHPROC_S, REQ_R1, IERR )
      ELSEIF (NMSG.EQ.2) THEN
        CALL MPI_STARTALL ( RDIM-NEIGHPROC_S, REQ_R2, IERR )
      ELSE
        CALL MPI_STARTALL ( RDIM-NEIGHPROC_S, REQ_R3, IERR )
      ENDIF
C
              !..Unpack Received messages as they arrive     
C
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
                    DO K=1,DOFH
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
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(IRECVLOC(I,J),K,IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
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
cdebug     print *, myproc, tot,nfini,index(1),index(2)
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(IRECVLOC(I,J),K,IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC2(IRECVLOC(I,J),K,IRK) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
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
C 
 9998 continue
      RETURN
      END SUBROUTINE

      SUBROUTINE UPDATER_ELEM_MOD3( VEC1, VEC2, VEC3, NMSG )
C
C--------------------------------------------------------------------------
C  Update 1, 2, or 3 Real Arrays's Ghost Cells using asynchronous
C  and persistent message-passing.
C
C  vjp  8/06/1999
C--------------------------------------------------------------------------
C 
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
C
      INTEGER,  INTENT(IN) ::  NMSG
      REAL(SZ), INTENT(INOUT) ::  VEC1(:,:),VEC2(:,:),VEC3(:,:)
      INTEGER :: N,I,J,JJ,K,NCOUNT,NFINI,TOT
C
                             !..Pack 1, 2, or 3 Messages
C      DO JJ=1,NEIGHPROC_S
C         J = MOD(JJ-1+MYPROC,NEIGHPROC_S)+1
      DO J=1,NEIGHPROC_S
         NCOUNT = 0
         DO I=1,NELEMSEND(J)
            DO K=1,DOFH
              NCOUNT = NCOUNT+1
              SENDBUF(NCOUNT,J)=VEC1(K,ISENDLOC(I,J))
            ENDDO
         ENDDO
         IF (NMSG.GT.1) THEN
           DO I=1,NELEMSEND(J)
             DO K=1,DOFH
               NCOUNT = NCOUNT+1
               SENDBUF(NCOUNT,J)=VEC2(K,ISENDLOC(I,J))
             ENDDO
           ENDDO
         ENDIF
         IF (NMSG.GT.2) THEN
           DO I=1,NELEMSEND(J)
             DO K=1,DOFH
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
C                    
              ! Send/receive messages to/from all neighbors
C
      IF (NMSG.EQ.1) THEN
        CALL MPI_STARTALL ( RDIM-NEIGHPROC_S, REQ_R1, IERR )
      ELSEIF (NMSG.EQ.2) THEN
        CALL MPI_STARTALL ( RDIM-NEIGHPROC_S, REQ_R2, IERR )
      ELSE
        CALL MPI_STARTALL ( RDIM-NEIGHPROC_S, REQ_R3, IERR )
      ENDIF
C
              !..Unpack Received messages as they arrive     
C
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
                    DO K=1,DOFH
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
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(K,IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
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
cdebug     print *, myproc, tot,nfini,index(1),index(2)
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC_R) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(K,IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC2(K,IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NELEMRECV(J)
                    DO K=1,DOFH
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
C 
 9998 continue
      RETURN
      END SUBROUTINE



      SUBROUTINE UPDATER3D( VEC )
C
C--------------------------------------------------------------------------
C  Update 1 Three-dimensional Real Arrays's Ghost Cells using asynchronous
C  and persistent message-passing.
C
C  tjc  6/24/2002
C--------------------------------------------------------------------------
C 
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
C
      REAL(SZ), INTENT(INOUT) ::  VEC(MNP,MNODES)
      INTEGER :: N,I,J,K,NCOUNT,NFINI,TOT
c$$$C
c$$$                             !..Pack Messages
c$$$      DO J=1,NEIGHPROC
c$$$         NCOUNT = 0
c$$$         DO I=1,NELEMSEND(J)
c$$$            DO K=1,MNODES
c$$$               NCOUNT = NCOUNT+1
c$$$               SENDBUF(NCOUNT,J)=VEC(ISENDLOC(I,J),K)
c$$$            ENDDO
c$$$         ENDDO
c$$$      ENDDO
c$$$C                    
c$$$              ! Send/receive messages to/from all neighbors
c$$$C
c$$$      CALL MPI_STARTALL ( RDIM, REQ_R3D, IERR )
c$$$C
c$$$              !..Unpack Received messages as they arrive     
c$$$C
c$$$      TOT = 0
c$$$      DO WHILE (TOT.LT.RDIM)
c$$$         DO N=1, RDIM
c$$$            INDEX(N) = 0
c$$$         ENDDO
c$$$         CALL MPI_WAITSOME( RDIM,REQ_R3D,NFINI,INDEX,STAT_R3D,IERR )
c$$$         TOT = TOT + NFINI
c$$$         DO N=1, NFINI
c$$$            IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
c$$$              IF (INDEX(N).LE.NEIGHPROC) THEN
c$$$                J = INDEX(N)
c$$$                NCOUNT = 0
c$$$                DO I=1,NELEMRECV(J)
c$$$                   DO K=1,MNODES
c$$$                      NCOUNT = NCOUNT+1
c$$$                      VEC(IRECVLOC(I,J),K) = RECVBUF(NCOUNT,J)
c$$$                   ENDDO
c$$$                ENDDO
c$$$              ENDIF
c$$$            ENDIF
c$$$         ENDDO
c$$$      ENDDO
c$$$C 
      STOP 'NOT UPDATED'
      RETURN
      END SUBROUTINE


      SUBROUTINE UPDATEC3D( VEC )
C
C--------------------------------------------------------------------------
C  Update 1 Three-dimensional Complex Arrays's Ghost Cells using asynchronous
C  and persistent message-passing.
C
C  tjc  6/24/2002
C--------------------------------------------------------------------------
C 
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
C
      COMPLEX, INTENT(INOUT) ::  VEC(MNP,MNODES)
      INTEGER :: N,I,J,K,NCOUNT,NFINI,TOT
c$$$C
c$$$                             !..Pack Messages
c$$$      DO J=1,NEIGHPROC
c$$$         NCOUNT = 0
c$$$         DO I=1,NELEMSEND(J)
c$$$            DO K=1,MNODES
c$$$               NCOUNT = NCOUNT+1
c$$$               SENDBUF(NCOUNT,J)=REAL(VEC(ISENDLOC(I,J),K))
c$$$               NCOUNT = NCOUNT+1
c$$$               SENDBUF(NCOUNT,J)=AIMAG(VEC(ISENDLOC(I,J),K))
c$$$            ENDDO
c$$$         ENDDO
c$$$      ENDDO
c$$$C                    
c$$$              ! Send/receive messages to/from all neighbors
c$$$C
c$$$      CALL MPI_STARTALL ( RDIM, REQ_C3D, IERR )
c$$$C
c$$$              !..Unpack Received messages as they arrive     
c$$$C
c$$$      TOT = 0
c$$$      DO WHILE (TOT.LT.RDIM)
c$$$         DO N=1, RDIM
c$$$            INDEX(N) = 0
c$$$         ENDDO
c$$$         CALL MPI_WAITSOME( RDIM,REQ_C3D,NFINI,INDEX,STAT_C3D,IERR )
c$$$         TOT = TOT + NFINI
c$$$         DO N=1, NFINI
c$$$            IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
c$$$              IF (INDEX(N).LE.NEIGHPROC) THEN
c$$$                J = INDEX(N)
c$$$                NCOUNT = 0
c$$$                DO I=1,NELEMRECV(J)
c$$$                   DO K=1,MNODES
c$$$                      VEC(IRECVLOC(I,J),K) = 
c$$$     &                   CMPLX(RECVBUF(NCOUNT+1,J),RECVBUF(NCOUNT+2,J))
c$$$                      NCOUNT = NCOUNT+2
c$$$                   ENDDO
c$$$                ENDDO
c$$$              ENDIF
c$$$            ENDIF
c$$$         ENDDO
c$$$      ENDDO
c$$$C 
      STOP 'NOT UPDATED'
      RETURN
      END SUBROUTINE

      SUBROUTINE MSG_BLOCKSYNC_START()
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif

      INTEGER,PARAMETER :: BLOCK = 10
      INTEGER :: SRC,DUM,STAT(MPI_STATUS_SIZE)

      SRC = MYPROC - BLOCK
      IF(SRC.LT.0) RETURN

      CALL MPI_RECV(DUM,1,MPI_INTEGER,SRC,TAG+1,MPI_COMM_WORLD,STAT,IERR)

      END SUBROUTINE 

      SUBROUTINE MSG_BLOCKSYNC_FINISH()
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif

      INTEGER,PARAMETER :: BLOCK = 10
      INTEGER :: DST,DUM,STAT(MPI_STATUS_SIZE)

      DST = MYPROC + BLOCK
      IF(DST.GE.MNPROC) RETURN

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

      call MPI_ALLREDUCE( sum_this, sum_of_this, kount,
     &     REALTYPE, MPI_SUM, COMM, ierr)

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

      call MPI_ALLREDUCE( max_this, max_of_this, kount,
     &     REALTYPE, MPI_MAX, COMM, ierr)

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

      call MPI_ALLREDUCE( min_this, min_of_this, kount,
     &     REALTYPE, MPI_MIN, COMM, ierr)

      min_this = min_of_this

      END SUBROUTINE para_min

      END MODULE MESSENGER_ELEM
       
