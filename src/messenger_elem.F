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
C.....2021-05-13 Younghun: Declare variables for 1D KWE... should it be private?      
      INTEGER, ALLOCATABLE :: IRECVBUF_CH(:),ISENDBUF_CH(:)
      INTEGER :: NUM_PE_RECV_CH,NUM_PE_SEND_CH,NUM_PE_COMM_CH
      INTEGER,ALLOCATABLE :: IPE_RECV_CH(:),IPE_SEND_CH(:)
      INTEGER,ALLOCATABLE :: NUM_CH_RECV(:),NUM_CH_SEND(:)
      INTEGER,ALLOCATABLE :: IRECVLOC_CH(:,:),ISENDLOC_CH(:,:)
      
      INTEGER, PRIVATE, ALLOCATABLE :: REQ_R1_CH(:), REQ_R2_CH(:)
      INTEGER, PRIVATE, ALLOCATABLE :: STAT_R1_CH(:,:), STAT_R2_CH(:,:)
      INTEGER, PRIVATE, ALLOCATABLE :: INDEX_CH(:)
      REAL(SZ), PRIVATE,ALLOCATABLE :: SENDBUF_CH(:,:), RECVBUF_CH(:,:)
C      
C.....2021-05-27 Younghun: Declare variables for passing lateral inflow... should it be private?      
      INTEGER, ALLOCATABLE :: IRECVBUF_LI(:),ISENDBUF_LI(:)
      INTEGER :: NUM_PE_RECV_LI,NUM_PE_SEND_LI,NUM_PE_COMM_LI
      INTEGER,ALLOCATABLE :: IPE_RECV_LI(:),IPE_SEND_LI(:)
      INTEGER,ALLOCATABLE :: NUM_LI_RECV(:),NUM_LI_SEND(:)
      INTEGER,ALLOCATABLE :: IRECVLOC_LI(:,:),ISENDLOC_LI(:,:)
      
      INTEGER, PRIVATE, ALLOCATABLE :: REQ_R1_LI(:), REQ_R2_LI(:)
      INTEGER, PRIVATE, ALLOCATABLE :: STAT_R1_LI(:,:), STAT_R2_LI(:,:)
      INTEGER, PRIVATE, ALLOCATABLE :: INDEX_LI(:)
      REAL(SZ), PRIVATE,ALLOCATABLE :: SENDBUF_LI(:,:), RECVBUF_LI(:,:)
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

C.....2021-05-12 Younghun:
      SUBROUTINE MSG_TABLE_CHANNEL()
      USE DG
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      INTEGER :: IDPROC,NLOCAL,I,J,JJ,NEIGH,K

      INTEGER,ALLOCATABLE :: IRECV_CH(:,:,:),ISEND_CH(:,:,:),ISEG_SEND(:,:)
      INTEGER,ALLOCATABLE :: IRECV_LI(:,:,:),ISEND_LI(:,:,:)
      
      ! INTEGER,ALLOCATABLE :: IPROC_RECV(:,:),IPROC_SEND(:,:,:)
      
      integer :: irecv,ged,isend
      ! allocate (IRECV_CHANNEL(nedges),ISEND_CHANNEL(nedges))
C
      OPEN(182,FILE=DIRNAME(1:LNAME)//'/'//'DG.19')
      READ(182,3010) NUM_PE_RECV_CH,NUM_PE_SEND_CH
      
      NUM_PE_COMM_CH = NUM_PE_RECV_CH + NUM_PE_SEND_CH
      ALLOCATE (INDEX_CH(NUM_PE_COMM_CH))
      
      ALLOCATE (IPE_RECV_CH(NUM_PE_RECV_CH),IPE_SEND_CH(NUM_PE_SEND_CH))
      ALLOCATE (NUM_CH_RECV(NUM_PE_RECV_CH),NUM_CH_SEND(NUM_PE_SEND_CH))
      ALLOCATE (IRECVLOC_CH(NEDGES,NUM_PE_RECV_CH),ISENDLOC_CH(NEDGES,NUM_PE_SEND_CH))

      ALLOCATE (IRECV_CH(NUM_PE_RECV_CH,NBOU,2))
      ALLOCATE (ISEND_CH(NUM_PE_SEND_CH,NBOU,2),ISEG_SEND(NUM_PE_SEND_CH,NBOU))

      DO K = 1,NUM_PE_RECV_CH
         READ(182,3010) IPE_RECV_CH(K),NUM_CH_RECV(K)
         READ(182,1130) (IRECV_CH(K,I,1),IRECV_CH(K,I,2),I=1,NUM_CH_RECV(K))
      ENDDO
      
      DO K = 1,NUM_PE_SEND_CH
         READ(182,3010) IPE_SEND_CH(K),NUM_CH_SEND(K)
         READ(182,1130) (ISEND_CH(K,I,1),ISEND_CH(K,I,2),I=1,NUM_CH_SEND(K))
      ENDDO
      
      DO K = 1,NUM_PE_RECV_CH
         DO I = 1,NUM_CH_RECV(K)
            DO J = 1, NCEDS
               IF ((NCEDNO(1,J).EQ.IRECV_CH(K,I,1) .AND. NCEDNO(2,J).EQ.IRECV_CH(K,I,2))
     &           .OR. (NCEDNO(2,J).EQ.IRECV_CH(K,I,1) .AND. NCEDNO(1,J).EQ.IRECV_CH(K,I,2))) THEN
                  ged = ncedn(J)
                  IRECVLOC_CH(I,K) = ged
                  EXIT
               ENDIF
               IF (J.EQ.NCEDS) THEN
                  print*,myproc,'recv ch null'
                  stop
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      
      DO K = 1,NUM_PE_SEND_CH
         DO I = 1,NUM_CH_SEND(K)
            DO J = 1, NCEDS
               IF ((NCEDNO(1,J).EQ.ISEND_CH(K,I,1) .AND. NCEDNO(2,J).EQ.ISEND_CH(K,I,2))
     &           .OR. (NCEDNO(2,J).EQ.ISEND_CH(K,I,1) .AND. NCEDNO(1,J).EQ.ISEND_CH(K,I,2))) THEN
                  ged = ncedn(J)
                  ISENDLOC_CH(I,K) = ged
                  EXIT
               ENDIF
               IF (J.EQ.NCEDS) THEN
                  print*,myproc,'send ch null'
                  stop
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      
      ! Do same thing for lateral inflow
      READ(182,3010) NUM_PE_RECV_LI,NUM_PE_SEND_LI
      
      NUM_PE_COMM_LI = NUM_PE_RECV_LI + NUM_PE_SEND_LI
      ALLOCATE (INDEX_LI(NUM_PE_COMM_LI))
      
      ALLOCATE (IPE_RECV_LI(NUM_PE_RECV_LI),IPE_SEND_LI(NUM_PE_SEND_LI))
      ALLOCATE (NUM_LI_RECV(NUM_PE_RECV_LI),NUM_LI_SEND(NUM_PE_SEND_LI))
      ALLOCATE (IRECVLOC_LI(NEDGES,NUM_PE_RECV_LI),ISENDLOC_LI(NEDGES,NUM_PE_SEND_LI))

      ALLOCATE (IRECV_LI(NUM_PE_RECV_LI,NBOU,2),ISEND_LI(NUM_PE_SEND_LI,NBOU,2))

      DO K = 1,NUM_PE_RECV_LI
         READ(182,3010) IPE_RECV_LI(K),NUM_LI_RECV(K)
         READ(182,1130) (IRECV_LI(K,I,1),IRECV_LI(K,I,2),I=1,NUM_LI_RECV(K))
      ENDDO
      
      DO K = 1,NUM_PE_SEND_LI
         READ(182,3010) IPE_SEND_LI(K),NUM_LI_SEND(K)
         READ(182,1130) (ISEND_LI(K,I,1),ISEND_LI(K,I,2),I=1,NUM_LI_SEND(K))
      ENDDO
      
      DO K = 1,NUM_PE_RECV_LI
         DO I = 1,NUM_LI_RECV(K)
            DO J = 1, NCEDS
               IF ((NCEDNO(1,J).EQ.IRECV_LI(K,I,1) .AND. NCEDNO(2,J).EQ.IRECV_LI(K,I,2))
     &           .OR. (NCEDNO(2,J).EQ.IRECV_LI(K,I,1) .AND. NCEDNO(1,J).EQ.IRECV_LI(K,I,2))) THEN
                  ged = ncedn(J)
                  IRECVLOC_LI(I,K) = J ! array QL in rhs_dg_runoff.F uses local channel numbering
                  EXIT
               ENDIF
               IF (J.EQ.NCEDS) THEN
                  print*,myproc,'RECV LI null'
                  stop
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      
      DO K = 1,NUM_PE_SEND_LI
         DO I = 1,NUM_LI_SEND(K)
            DO J = 1, NCEDS
               IF ((NCEDNO(1,J).EQ.ISEND_LI(K,I,1) .AND. NCEDNO(2,J).EQ.ISEND_LI(K,I,2))
     &           .OR. (NCEDNO(2,J).EQ.ISEND_LI(K,I,1) .AND. NCEDNO(1,J).EQ.ISEND_LI(K,I,2))) THEN
                  ged = ncedn(J)
                  ISENDLOC_LI(I,K) = J ! array QL in rhs_dg_runoff.F uses local channel numbering
                  EXIT
               ENDIF
               IF (J.EQ.NCEDS) THEN
                  print*,myproc,'send LI null'
                  stop
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      
      ! DO K = 1,NUM_PE_RECV_CH
         ! WRITE(183,1130) (IRECVLOC_CH(I,K),I=1,NUM_CH_RECV(K))
      ! ENDDO
      ! DO K = 1,NUM_PE_SEND_CH
         ! WRITE(184,1130) (ISENDLOC_CH(I,K),I=1,NUM_CH_SEND(K))
      ! ENDDO
      
      ! DO K = 1,NUM_PE_SEND_CH
         ! DO I = 1,NUM_CH_SEND(K)
            ! iSEND = 0
            ! DO J = 1,ISEND_CH(K,I)
               ! if (segtype(J).eq.6) iSEND = iSEND + nvell(J)-1
            ! END
            ! isend = isend + ISEG_SEND(K,I)
            ! ged = ncedn(iSEND)
            ! ISENDLOC_CH(I,K) = ged
         ! ENDDO
      ! ENDDO      
      
      
      
      ! DO J = 1,NBOU
         ! DO K = 1,NUM_PE_RECV_CH
            ! DO I = 1, NUM_CH_RECV(K)
               ! IF (IRECV_CH(K,I).EQ.J) THEN
                  ! IPROC_RECV(J,1) = IPROC_RECV(J,1) + 1
                  ! JJ = IPROC_RECV(J,1)+1
                  ! IPROC_RECV(J,JJ) = IPE_RECV_CH(K)
               ! ENDIF
            ! ENDDO
         ! ENDDO
      ! ENDDO
      
      ! DO J = 1,NBOU
         ! DO K = 1,NUM_PE_SEND_CH
            ! DO I = 1, NUM_CH_SEND(K)
               ! IF (ISEND_CH(K,I).EQ.J) THEN
                  ! IPROC_SEND(J,1,1) = IPROC_SEND(J,1,1) + 1
                  ! JJ = IPROC_SEND(J,1,1)+1
                  ! IPROC_SEND(J,JJ,1) = IPE_SEND_CH(K)
                  ! IPROC_SEND(J,JJ,2) = ISEG_SEND(K,I)
               ! ENDIF
            ! ENDDO
         ! ENDDO
      ! ENDDO
      
      
      ! OPEN(183,FILE=DIRNAME(1:LNAME)//'/'//'fort.183')
      ! OPEN(184,FILE=DIRNAME(1:LNAME)//'/'//'fort.184')
      ! DO J = 1,NBOU
         ! ! WRITE(183,3010) IPE_RECV_CH(K),NUM_CH_RECV(K)
         ! WRITE(183,1130) IPROC_RECV(J,1),(IPROC_RECV(J,I+1),I=1,NUM_PE_RECV_CH)         
         ! ! WRITE(183,3010) IPE_RECV_CH(K),NUM_CH_RECV(K)
         ! WRITE(184,1130) IPROC_SEND(J,1,1),(IPROC_SEND(J,I+1,1),IPROC_SEND(J,I+1,2),I=1,NUM_PE_SEND_CH)
      ! ENDDO
      
      ! irecv = 1
      ! OPEN(185,FILE=DIRNAME(1:LNAME)//'/'//'fort.185')
      ! DO K = 1,NBOU
         ! if (iproc_recv(K,1).ne.0) then
            ! ! WRITE(185,*) irecv
            ! ged = ncedn(irecv)
            ! IRECV_CHANNEL(ged) = iproc_recv(K,2)
            ! ! WRITE(185,*) ged,IRECV_CHANNEL(ged)
            ! WRITE(185,*)IRECV_CHANNEL(ged),NCEDNO(1,irecv),NCEDNO(2,irecv)
         ! endif
         ! if (segtype(k).eq.6) irecv = irecv + nvell(K)-1
      ! enddo

      ! isend = 0
      ! OPEN(186,FILE=DIRNAME(1:LNAME)//'/'//'fort.186')
      ! DO K = 1,NBOU
         ! if (iproc_send(K,1,1).ne.0) then
            ! DO I = 1,iproc_send(K,1,1)
               ! J = isend + IPROC_SEND(K,I+1,2)
               ! ged = ncedn(J)
               ! ISEND_CHANNEL(ged) = iproc_send(K,I+1,1)
               ! ! WRITE(186,*) J
               ! ! WRITE(186,*) ged,ISEND_CHANNEL(ged)
               ! WRITE(186,*)ISEND_CHANNEL(ged),NCEDNO(1,J),NCEDNO(2,J)
            ! ENDDO
         ! endif
         ! if (segtype(k).eq.6) isend = isend + nvell(K)-1
      ! enddo
      
1130  FORMAT(8X,9I8)
3010  FORMAT(8X,2I8)

      END SUBROUTINE MSG_TABLE_CHANNEL
      
C.....2021-05-13 Younghun
      SUBROUTINE MESSAGE_START_CHANNEL()
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      INTEGER :: J,NCOMM_CH_S,NCOMM_CH_R
      INTEGER :: NCOMM_LI_S,NCOMM_LI_R
C
      NCOMM_CH_R = 0
      NCOMM_CH_S = 0
      DO J=1,NUM_PE_RECV_CH
        NCOMM_CH_R = NCOMM_CH_R + NUM_CH_RECV(J)
      ENDDO
      DO J=1,NUM_PE_SEND_CH
        NCOMM_CH_S = NCOMM_CH_S + NUM_CH_SEND(J)
      ENDDO
C
      ALLOCATE ( SENDBUF_CH(NCOMM_CH_S*DOFH*4,NUM_PE_SEND_CH) )
      ALLOCATE ( RECVBUF_CH(NCOMM_CH_R*DOFH*4,NUM_PE_RECV_CH) )
C
      ALLOCATE ( REQ_R1_CH(NUM_PE_COMM_CH),REQ_R2_CH(NUM_PE_COMM_CH))
C
      ALLOCATE ( STAT_R1_CH(MPI_STATUS_SIZE,NUM_PE_COMM_CH),       
     &           STAT_R2_CH(MPI_STATUS_SIZE,NUM_PE_COMM_CH))
C
            !  Setup persistent structures for real arrays
C
      DO J=1,NUM_PE_RECV_CH
         ! print*,'recv',NUM_CH_RECV(J),IPE_RECV_CH(J)
         CALL MPI_RECV_INIT ( RECVBUF_CH(1,J), DOFH*NUM_CH_RECV(J), 
     &     REALTYPE,IPE_RECV_CH(J), TAG, MPI_COMM_WORLD,
     &     REQ_R1_CH(J),IERR)
      ENDDO
      DO J=1,NUM_PE_SEND_CH
         ! print*,'send',NUM_CH_SEND(J),IPE_SEND_CH(J)
         CALL MPI_SEND_INIT ( SENDBUF_CH(1,J), DOFH*NUM_CH_SEND(J), 
     &     REALTYPE,IPE_SEND_CH(J), TAG, MPI_COMM_WORLD,
     &     REQ_R1_CH(J+NUM_PE_RECV_CH),IERR)
      ENDDO
C      
C.....Lateral inflow
      NCOMM_LI_R = 0
      NCOMM_LI_S = 0
      DO J=1,NUM_PE_RECV_LI
        NCOMM_LI_R = NCOMM_LI_R + NUM_LI_RECV(J)
      ENDDO
      DO J=1,NUM_PE_SEND_LI
        NCOMM_LI_S = NCOMM_LI_S + NUM_LI_SEND(J)
      ENDDO
C
      ALLOCATE ( SENDBUF_LI(NCOMM_LI_S*DOFH*4,NUM_PE_SEND_LI) )
      ALLOCATE ( RECVBUF_LI(NCOMM_LI_R*DOFH*4,NUM_PE_RECV_LI) )
C
      ALLOCATE ( REQ_R1_LI(NUM_PE_COMM_LI),REQ_R2_LI(NUM_PE_COMM_LI))
C
      ALLOCATE ( STAT_R1_LI(MPI_STATUS_SIZE,NUM_PE_COMM_LI),       
     &           STAT_R2_LI(MPI_STATUS_SIZE,NUM_PE_COMM_LI))
C
            !  Setup persistent structures for real arrays
C
      DO J=1,NUM_PE_RECV_LI
         ! print*,'recv',NUM_LI_RECV(J),IPE_RECV_LI(J)
         CALL MPI_RECV_INIT ( RECVBUF_LI(1,J), DOFH*NUM_LI_RECV(J), 
     &     REALTYPE,IPE_RECV_LI(J), TAG, MPI_COMM_WORLD,
     &     REQ_R1_LI(J),IERR)
      ENDDO
      DO J=1,NUM_PE_SEND_LI
         ! print*,'send',NUM_LI_SEND(J),IPE_SEND_LI(J)
         CALL MPI_SEND_INIT ( SENDBUF_LI(1,J), DOFH*NUM_LI_SEND(J), 
     &     REALTYPE,IPE_SEND_LI(J), TAG, MPI_COMM_WORLD,
     &     REQ_R1_LI(J+NUM_PE_RECV_LI),IERR)
      ENDDO

      RETURN
      END SUBROUTINE
      
      SUBROUTINE UPDATER_CHANNEL( VEC1, VEC2, IRK, NMSG )
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
C
      INTEGER,  INTENT(IN) ::  IRK, NMSG
      REAL(SZ), INTENT(INOUT) ::  VEC1(:,:,:),VEC2(:,:,:)
      INTEGER :: N,I,J,JJ,K,NCOUNT,NFINI,TOT
C
                             !..Pack 1, 2, or 3 Messages
      DO J=1,NUM_PE_SEND_CH
         NCOUNT = 0
         DO I=1,NUM_CH_SEND(J)
            DO K=1,DOFH
              NCOUNT = NCOUNT+1
              SENDBUF_CH(NCOUNT,J)=VEC1(K,ISENDLOC_CH(I,J),IRK)
            ENDDO
         ENDDO
         IF (NMSG.GT.1) THEN
           DO I=1,NUM_CH_SEND(J)
             DO K=1,DOFH
               NCOUNT = NCOUNT+1
               SENDBUF_CH(NCOUNT,J)=VEC2(K,ISENDLOC_CH(I,J),IRK)
             ENDDO
           ENDDO
         ENDIF

         ! Start sending a message
         IF (NMSG.EQ.1) THEN
            CALL MPI_START ( REQ_R1_CH(J+NUM_PE_RECV_CH), IERR )
         ELSE
            CALL MPI_START ( REQ_R2_CH(J+NUM_PE_RECV_CH), IERR )
         ENDIF
      ENDDO
C                    
              ! Send/receive messages to/from all neighbors
C
      IF (NMSG.EQ.1) THEN
        CALL MPI_STARTALL ( NUM_PE_COMM_CH-NUM_PE_SEND_CH, REQ_R1_CH, IERR )
      ELSE
        CALL MPI_STARTALL ( NUM_PE_COMM_CH-NUM_PE_SEND_CH, REQ_R2_CH, IERR )
      ENDIF
C
              !..Unpack Received messages as they arrive     
C
      IF (NMSG.EQ.1) THEN   
        TOT = 0
        DO WHILE (TOT.LT.NUM_PE_COMM_CH)
           DO N=1, NUM_PE_COMM_CH
              INDEX_CH(N) = 0
           ENDDO
           CALL MPI_WAITSOME( NUM_PE_COMM_CH,REQ_R1_CH,NFINI,INDEX_CH,STAT_R1_CH,IERR )
           TOT = TOT + NFINI
           DO N=1, NFINI
              IF (INDEX_CH(N).GT.0.AND.INDEX_CH(N).LE.NUM_PE_COMM_CH)  THEN
                IF (INDEX_CH(N).LE.NUM_PE_RECV_CH) THEN
                  J = INDEX_CH(N)
                  NCOUNT = 0
                  DO I=1,NUM_CH_RECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(K,IRECVLOC_CH(I,J),IRK) = RECVBUF_CH(NCOUNT,J)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
        GOTO 999
      ELSEIF (NMSG.EQ.2) THEN
        TOT = 0
        DO WHILE (TOT.LT.NUM_PE_COMM_CH)
           DO N=1, NUM_PE_COMM_CH
              INDEX_CH(N) = 0
           ENDDO
           CALL MPI_WAITSOME( NUM_PE_COMM_CH,REQ_R2_CH,NFINI,INDEX_CH,STAT_R2_CH,IERR )
           TOT = TOT + NFINI
           DO N=1, NFINI
              IF (INDEX_CH(N).GT.0.AND.INDEX_CH(N).LE.NUM_PE_COMM_CH)  THEN
                IF (INDEX_CH(N).LE.NUM_PE_RECV_CH) THEN
                  J = INDEX_CH(N)
                  NCOUNT = 0
                  DO I=1,NUM_CH_RECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(K,IRECVLOC_CH(I,J),IRK) = RECVBUF_CH(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NUM_CH_RECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC2(K,IRECVLOC_CH(I,J),IRK) = RECVBUF_CH(NCOUNT,J)
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
      
      SUBROUTINE UPDATER_LATERALINFLOW( VEC1, VEC2, IRK, NMSG )
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
C
      INTEGER,  INTENT(IN) ::  IRK, NMSG
      REAL(SZ), INTENT(INOUT) ::  VEC1(:,:),VEC2(:,:,:)
      INTEGER :: N,I,J,JJ,K,NCOUNT,NFINI,TOT
C
                             !..Pack 1, 2, or 3 Messages
      DO J=1,NUM_PE_SEND_LI
         NCOUNT = 0
         DO I=1,NUM_LI_SEND(J)
            DO K=1,DOFH
              NCOUNT = NCOUNT+1
              SENDBUF_LI(NCOUNT,J)=VEC1(K,ISENDLOC_LI(I,J))
            ENDDO
         ENDDO
         IF (NMSG.GT.1) THEN
           DO I=1,NUM_LI_SEND(J)
             DO K=1,DOFH
               NCOUNT = NCOUNT+1
               SENDBUF_LI(NCOUNT,J)=VEC2(K,ISENDLOC_LI(I,J),IRK)
             ENDDO
           ENDDO
         ENDIF

         ! Start sending a message
         IF (NMSG.EQ.1) THEN
            CALL MPI_START ( REQ_R1_LI(J+NUM_PE_RECV_LI), IERR )
         ELSE
            CALL MPI_START ( REQ_R2_LI(J+NUM_PE_RECV_LI), IERR )
         ENDIF
      ENDDO
C                    
              ! Send/receive messages to/from all neighbors
C
      IF (NMSG.EQ.1) THEN
        CALL MPI_STARTALL ( NUM_PE_COMM_LI-NUM_PE_SEND_LI, REQ_R1_LI, IERR )
      ELSE
        CALL MPI_STARTALL ( NUM_PE_COMM_LI-NUM_PE_SEND_LI, REQ_R2_LI, IERR )
      ENDIF
C
              !..Unpack Received messages as they arrive     
C
      IF (NMSG.EQ.1) THEN   
        TOT = 0
        DO WHILE (TOT.LT.NUM_PE_COMM_LI)
           DO N=1, NUM_PE_COMM_LI
              INDEX_LI(N) = 0
           ENDDO
           CALL MPI_WAITSOME( NUM_PE_COMM_LI,REQ_R1_LI,NFINI,INDEX_LI,STAT_R1_LI,IERR )
           TOT = TOT + NFINI
           DO N=1, NFINI
              IF (INDEX_LI(N).GT.0.AND.INDEX_LI(N).LE.NUM_PE_COMM_LI)  THEN
                IF (INDEX_LI(N).LE.NUM_PE_RECV_LI) THEN
                  J = INDEX_LI(N)
                  NCOUNT = 0
                  DO I=1,NUM_LI_RECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(K,IRECVLOC_LI(I,J)) = RECVBUF_LI(NCOUNT,J)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
           ENDDO
        ENDDO
        GOTO 999
      ELSEIF (NMSG.EQ.2) THEN
        TOT = 0
        DO WHILE (TOT.LT.NUM_PE_COMM_LI)
           DO N=1, NUM_PE_COMM_LI
              INDEX_LI(N) = 0
           ENDDO
           CALL MPI_WAITSOME( NUM_PE_COMM_LI,REQ_R2_LI,NFINI,INDEX_LI,STAT_R2_LI,IERR )
           TOT = TOT + NFINI
           DO N=1, NFINI
              IF (INDEX_LI(N).GT.0.AND.INDEX_LI(N).LE.NUM_PE_COMM_LI)  THEN
                IF (INDEX_LI(N).LE.NUM_PE_RECV_LI) THEN
                  J = INDEX_LI(N)
                  NCOUNT = 0
                  DO I=1,NUM_LI_RECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC1(K,IRECVLOC_LI(I,J)) = RECVBUF_LI(NCOUNT,J)
                    ENDDO
                  ENDDO
                  DO I=1,NUM_LI_RECV(J)
                    DO K=1,DOFH
                      NCOUNT = NCOUNT+1
                      VEC2(K,IRECVLOC_LI(I,J),IRK) = RECVBUF_LI(NCOUNT,J)
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
      
      END MODULE MESSENGER_ELEM
       
