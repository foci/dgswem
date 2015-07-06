C******************************************************************************
C  last changes in this file VERSION 44.18                                    *
C  S.Bunya changed this file a bit. 07/13/2005                                *
C****************************************************************************** 
C 
      MODULE MESSENGER
C
      USE SIZES
      USE GLOBAL, ONLY: C3D
      USE DG, ONLY: dofh
Csb-PDG1
      USE DIFF45_41,ONLY : MNODES
C--

#ifdef HAVE_MPI_MOD
      use mpi  
      IMPLICIT NONE
#else
      IMPLICIT NONE
#endif
C
C--------------------------------------------------------------------------
C  This module supplies the MPI Message-Passing Interface.
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

      INTEGER,SAVE,PRIVATE :: REALTYPE, DBLETYPE   
Casey 121126: Changed some variables to public.
      INTEGER,SAVE :: NEIGHPROC
      INTEGER,SAVE,PRIVATE ::  RDIM, IERR
      INTEGER,SAVE,PRIVATE ::  TAG = 100
      LOGICAL,SAVE,PRIVATE,ALLOCATABLE :: RESNODE(:)
C
Casey 121126: Changed some variables to public.
      INTEGER, ALLOCATABLE :: NNODRECV(:),IRECVLOC(:,:)
      INTEGER, PRIVATE, ALLOCATABLE :: IPROC(:), NNODELOC(:),
     &    NNODSEND(:), IBELONGTO(:),ISENDLOC(:,:), 
     &    ISENDBUF(:,:), IRECVBUF(:,:)
C
      INTEGER, PRIVATE, ALLOCATABLE :: REQ_I1(:), REQ_I2(:)
      INTEGER, PRIVATE, ALLOCATABLE :: STAT_I1(:,:), STAT_I2(:,:)
      INTEGER, PRIVATE, ALLOCATABLE :: REQ_R1(:), REQ_R2(:), REQ_R3(:)
      INTEGER, PRIVATE, ALLOCATABLE :: STAT_R1(:,:), STAT_R2(:,:), 
     &                                 STAT_R3(:,:)
      INTEGER, PRIVATE, ALLOCATABLE :: REQ_R3D(:), STAT_R3D(:,:)
      INTEGER, PRIVATE, ALLOCATABLE :: REQ_C3D(:), STAT_C3D(:,:)
      INTEGER, PRIVATE, ALLOCATABLE :: INDEX(:)
      REAL(SZ), PRIVATE,ALLOCATABLE :: SENDBUF(:,:), RECVBUF(:,:)
C--
C

C---------------------end of data declarations--------------------------------C


      CONTAINS


      SUBROUTINE MSG_TYPES()
#ifdef SWAN
Casey 121128: Added for global I/O.
      USE GLOBAL, ONLY: float_type
#endif
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
c
#ifdef CRAY
#ifdef REAL4
      REALTYPE = MPI_REAL4
      DBLETYPE = MPI_REAL8
#else
      REALTYPE = MPI_REAL8
      DBLETYPE = MPI_REAL8
#endif
#else
#ifdef REAL4
      REALTYPE = MPI_REAL
      DBLETYPE = MPI_DOUBLE_PRECISION
#else
      REALTYPE = MPI_DOUBLE_PRECISION
      DBLETYPE = MPI_DOUBLE_PRECISION
#endif
#endif
c
#ifdef SWAN
Casey 121128: Added for global I/O.
      float_type = REALTYPE
#endif
      RETURN
      END  SUBROUTINE


 
      SUBROUTINE MSG_TABLE () 
C
C--------------------------------------------------------------------------
C  Routine preforms following steps:
C
C   (1) Read Message-Passing Information from file "fort.18"
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
#ifdef SWAN
Casey 121128: Need the following variables.
      USE GLOBAL, ONLY: NODES_LG,NP_G
#endif
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      INTEGER :: IDPROC,NLOCAL,I,J
C
      OPEN(18,FILE=DIRNAME(1:LNAME)//'/'//'fort.18')
C
#ifdef SWAN
Casey 121128: Need the local-to-global vertex numbering
C             for the output of global files from SWAN.
      READ(18,'(8X,3I8)') NP_G
      ALLOCATE ( NODES_LG(MNP) )
      DO I=1,MNP
         READ(18,*) NODES_LG(I)
      ENDDO
#endif
C
      READ(18,3010) IDPROC,NLOCAL    
C
      ALLOCATE ( NNODELOC(NLOCAL) )
C
      READ(18,1130) (NNODELOC(I), I=1,NLOCAL)
C
      ALLOCATE ( IBELONGTO(MNP),RESNODE(MNP) )
C
      DO I=1,MNP
         IBELONGTO(I) = 0
      ENDDO
      DO I=1,NLOCAL
         IBELONGTO(NNODELOC(I)) = IDPROC + 1
      ENDDO
      DO I=1, MNP
         IF (IBELONGTO(I)-1.EQ.MYPROC) THEN
           RESNODE(I) = .TRUE.
         ELSE 
           RESNODE(I) = .FALSE.
         ENDIF
      ENDDO
C
      READ(18,3015) NEIGHPROC
C
      RDIM = 2*NEIGHPROC
      ALLOCATE( INDEX(RDIM) )
C
      ALLOCATE( IPROC(NEIGHPROC),NNODRECV(NEIGHPROC) )
      ALLOCATE( IRECVLOC(MNP,NEIGHPROC) )
C
      DO J=1,NEIGHPROC
         READ(18,3010) IPROC(J),NNODRECV(J)
         READ(18,1130) (IRECVLOC(I,J), I=1,NNODRECV(J))
      ENDDO
C
      ALLOCATE( NNODSEND(NEIGHPROC) )
      ALLOCATE( ISENDLOC(MNP,NEIGHPROC) )
C
      DO J=1,NEIGHPROC
         READ(18,3010) IPROC(J),NNODSEND(J)
         READ(18,1130) (ISENDLOC(I,J), I=1,NNODSEND(J))
      ENDDO
C
      CLOSE(18)
      RETURN
C
1130  FORMAT(8X,9I8)
3010  FORMAT(8X,2I8)
3015  FORMAT(8X,I8)
      END SUBROUTINE


      SUBROUTINE MESSAGE_START ()
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
      INTEGER :: J
C
      ALLOCATE ( ISENDBUF(MNP,NEIGHPROC), IRECVBUF(MNP,NEIGHPROC) )
C
      IF (C3D) THEN
         ALLOCATE ( SENDBUF(2*MNP*MNODES,NEIGHPROC) )
         ALLOCATE ( RECVBUF(2*MNP*MNODES,NEIGHPROC) )
      ELSE
         ALLOCATE ( SENDBUF(MNP,NEIGHPROC) )
         ALLOCATE ( RECVBUF(MNP,NEIGHPROC) )
      ENDIF
C
      ALLOCATE ( REQ_I1(RDIM),REQ_I2(RDIM) )
      ALLOCATE ( REQ_R1(RDIM),REQ_R2(RDIM),REQ_R3(RDIM) )
C
      ALLOCATE ( STAT_I1(MPI_STATUS_SIZE,RDIM),       
     &           STAT_I2(MPI_STATUS_SIZE,RDIM) )

      ALLOCATE ( STAT_R1(MPI_STATUS_SIZE,RDIM),       
     &           STAT_R2(MPI_STATUS_SIZE,RDIM),
     &           STAT_R3(MPI_STATUS_SIZE,RDIM) )
C
      IF (C3D) THEN
         ALLOCATE ( REQ_R3D(RDIM) )
         ALLOCATE ( STAT_R3D(MPI_STATUS_SIZE,RDIM) )
         ALLOCATE ( REQ_C3D(RDIM) )
         ALLOCATE ( STAT_C3D(MPI_STATUS_SIZE,RDIM) )
      ENDIF
C
             !  Setup persistent structures for integer arrays
C
      DO J=1,NEIGHPROC   
         CALL MPI_RECV_INIT ( IRECVBUF(1,J), NNODRECV(J), 
     &     MPI_INTEGER,IPROC(J), TAG, MPI_COMM_WORLD,
     &     REQ_I1(J),IERR)
      ENDDO
      DO J=1,NEIGHPROC   
         CALL MPI_SEND_INIT ( ISENDBUF(1,J), NNODSEND(J), 
     &    MPI_INTEGER,IPROC(J), TAG, MPI_COMM_WORLD,
     &    REQ_I1(J+NEIGHPROC),IERR )
      ENDDO
C
C
      DO J=1,NEIGHPROC   
         CALL MPI_RECV_INIT ( IRECVBUF(1,J), 2*NNODRECV(J), 
     &     MPI_INTEGER,IPROC(J), TAG, MPI_COMM_WORLD,
     &     REQ_I2(J),IERR)
      ENDDO
      DO J=1,NEIGHPROC   
         CALL MPI_SEND_INIT ( ISENDBUF(1,J), 2*NNODSEND(J), 
     &    MPI_INTEGER,IPROC(J), TAG, MPI_COMM_WORLD,
     &    REQ_I2(J+NEIGHPROC),IERR )
      ENDDO
C
            !  Setup persistent structures for real arrays
C
      DO J=1,NEIGHPROC  
         CALL MPI_RECV_INIT ( RECVBUF(1,J), NNODRECV(J), 
     &     REALTYPE,IPROC(J), TAG, MPI_COMM_WORLD,
     &     REQ_R1(J),IERR)
      ENDDO
      DO J=1,NEIGHPROC  
         CALL MPI_SEND_INIT ( SENDBUF(1,J), NNODSEND(J), 
     &     REALTYPE,IPROC(J), TAG, MPI_COMM_WORLD,
     &     REQ_R1(J+NEIGHPROC),IERR)
      ENDDO
C
      DO J=1,NEIGHPROC  
         CALL MPI_RECV_INIT ( RECVBUF(1,J), 2*NNODRECV(J), 
     &     REALTYPE,IPROC(J), TAG, MPI_COMM_WORLD,
     &     REQ_R2(J),IERR)
      ENDDO
      DO J=1,NEIGHPROC  
         CALL MPI_SEND_INIT ( SENDBUF(1,J), 2*NNODSEND(J), 
     &     REALTYPE,IPROC(J), TAG, MPI_COMM_WORLD,
     &     REQ_R2(J+NEIGHPROC),IERR)
      ENDDO
C
      DO J=1,NEIGHPROC  
         CALL MPI_RECV_INIT ( RECVBUF(1,J), 3*NNODRECV(J), 
     &     REALTYPE,IPROC(J), TAG, MPI_COMM_WORLD,
     &     REQ_R3(J),IERR)
      ENDDO
      DO J=1,NEIGHPROC  
         CALL MPI_SEND_INIT ( SENDBUF(1,J), 3*NNODSEND(J), 
     &     REALTYPE,IPROC(J), TAG, MPI_COMM_WORLD,
     &     REQ_R3(J+NEIGHPROC),IERR)
      ENDDO
C
      IF (C3D) THEN
         DO J=1,NEIGHPROC  
            CALL MPI_RECV_INIT ( RECVBUF(1,J), MNODES*NNODRECV(J), 
     &        REALTYPE,IPROC(J), TAG, MPI_COMM_WORLD,
     &        REQ_R3D(J),IERR)
         ENDDO
         DO J=1,NEIGHPROC  
            CALL MPI_SEND_INIT ( SENDBUF(1,J), MNODES*NNODSEND(J), 
     &        REALTYPE,IPROC(J), TAG, MPI_COMM_WORLD,
     &        REQ_R3D(J+NEIGHPROC),IERR)
         ENDDO
         DO J=1,NEIGHPROC  
            CALL MPI_RECV_INIT ( RECVBUF(1,J), 2*MNODES*NNODRECV(J), 
     &        REALTYPE,IPROC(J), TAG, MPI_COMM_WORLD,
     &        REQ_C3D(J),IERR)
         ENDDO
         DO J=1,NEIGHPROC  
            CALL MPI_SEND_INIT ( SENDBUF(1,J), 2*MNODES*NNODSEND(J), 
     &        REALTYPE,IPROC(J), TAG, MPI_COMM_WORLD,
     &        REQ_C3D(J+NEIGHPROC),IERR)
         ENDDO
      ENDIF
C
      RETURN
      END SUBROUTINE


      SUBROUTINE UPDATEI( IVEC1, IVEC2, NMSG )
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
      INTEGER,   INTENT(INOUT) :: IVEC1(*),IVEC2(*)
      INTEGER :: N,I,J,NCOUNT,NFINI,TOT
C
                             !..Pack 1 or 2 Messages
      DO J=1,NEIGHPROC
         NCOUNT = 0
         DO I=1,NNODSEND(J)
            NCOUNT = NCOUNT+1
            ISENDBUF(NCOUNT,J)=IVEC1(ISENDLOC(I,J))
         ENDDO
         IF (NMSG.GT.1) THEN
           DO I=1,NNODSEND(J)
              NCOUNT = NCOUNT+1
              ISENDBUF(NCOUNT,J)=IVEC2(ISENDLOC(I,J))
           ENDDO
         ENDIF
      ENDDO
C                     
                          ! Send/receive messages to/from all neighbors
      IF (NMSG.EQ.1) THEN
        CALL MPI_STARTALL ( RDIM, REQ_I1, IERR )
      ELSE
        CALL MPI_STARTALL ( RDIM, REQ_I2, IERR )
      ENDIF
C
                          !..Unpack Received messages as they arrive  

      IF (NMSG.EQ.1) THEN   
        TOT = 0
        DO WHILE (TOT.LT.RDIM)
           DO N=1, RDIM
              INDEX(N) = 0
           ENDDO
           CALL MPI_WAITSOME( RDIM,REQ_I1,NFINI,INDEX,STAT_I1,IERR )
           TOT = TOT + NFINI
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NNODRECV(J)
                     NCOUNT = NCOUNT+1
                     IVEC1(IRECVLOC(I,J)) = IRECVBUF(NCOUNT,J)
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
           CALL MPI_WAITSOME( RDIM,REQ_I2,NFINI,INDEX,STAT_I2,IERR )
           TOT = TOT + NFINI
           DO N=1, NFINI
              IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
                IF (INDEX(N).LE.NEIGHPROC) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NNODRECV(J)
                     NCOUNT = NCOUNT+1
                     IVEC1(IRECVLOC(I,J)) = IRECVBUF(NCOUNT,J)
                  ENDDO
                  DO I=1,NNODRECV(J)
                     NCOUNT = NCOUNT+1
                     IVEC2(IRECVLOC(I,J)) = IRECVBUF(NCOUNT,J)
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


      SUBROUTINE UPDATER( VEC1, VEC2, VEC3, NMSG )
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
      REAL(SZ), INTENT(INOUT) ::  VEC1(*),VEC2(*),VEC3(*)
      INTEGER :: N,I,J,NCOUNT,NFINI,TOT
C
                             !..Pack 1, 2, or 3 Messages
      DO J=1,NEIGHPROC
         NCOUNT = 0
         DO I=1,NNODSEND(J)
            NCOUNT = NCOUNT+1
            SENDBUF(NCOUNT,J)=VEC1(ISENDLOC(I,J))
         ENDDO
         IF (NMSG.GT.1) THEN
           DO I=1,NNODSEND(J)
              NCOUNT = NCOUNT+1
              SENDBUF(NCOUNT,J)=VEC2(ISENDLOC(I,J))
           ENDDO
         ENDIF
         IF (NMSG.GT.2) THEN
           DO I=1,NNODSEND(J)
              NCOUNT = NCOUNT+1
              SENDBUF(NCOUNT,J)=VEC3(ISENDLOC(I,J))
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
                IF (INDEX(N).LE.NEIGHPROC) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NNODRECV(J)
                     NCOUNT = NCOUNT+1
                     VEC1(IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
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
                IF (INDEX(N).LE.NEIGHPROC) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NNODRECV(J)
                     NCOUNT = NCOUNT+1
                     VEC1(IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
                  ENDDO
                  DO I=1,NNODRECV(J)
                     NCOUNT = NCOUNT+1
                     VEC2(IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
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
                IF (INDEX(N).LE.NEIGHPROC) THEN
                  J = INDEX(N)
                  NCOUNT = 0
                  DO I=1,NNODRECV(J)
                     NCOUNT = NCOUNT+1
                     VEC1(IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
                  ENDDO
                  DO I=1,NNODRECV(J)
                     NCOUNT = NCOUNT+1
                     VEC2(IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
                  ENDDO
                  DO I=1,NNODRECV(J)
                     NCOUNT = NCOUNT+1
                     VEC3(IRECVLOC(I,J)) = RECVBUF(NCOUNT,J)
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
C
                             !..Pack Messages
      DO J=1,NEIGHPROC
         NCOUNT = 0
         DO I=1,NNODSEND(J)
            DO K=1,MNODES
               NCOUNT = NCOUNT+1
               SENDBUF(NCOUNT,J)=VEC(ISENDLOC(I,J),K)
            ENDDO
         ENDDO
      ENDDO
C                    
              ! Send/receive messages to/from all neighbors
C
      CALL MPI_STARTALL ( RDIM, REQ_R3D, IERR )
C
              !..Unpack Received messages as they arrive     
C
      TOT = 0
      DO WHILE (TOT.LT.RDIM)
         DO N=1, RDIM
            INDEX(N) = 0
         ENDDO
         CALL MPI_WAITSOME( RDIM,REQ_R3D,NFINI,INDEX,STAT_R3D,IERR )
         TOT = TOT + NFINI
         DO N=1, NFINI
            IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
              IF (INDEX(N).LE.NEIGHPROC) THEN
                J = INDEX(N)
                NCOUNT = 0
                DO I=1,NNODRECV(J)
                   DO K=1,MNODES
                      NCOUNT = NCOUNT+1
                      VEC(IRECVLOC(I,J),K) = RECVBUF(NCOUNT,J)
                   ENDDO
                ENDDO
              ENDIF
            ENDIF
         ENDDO
      ENDDO
C 
      RETURN
      END SUBROUTINE

      SUBROUTINE UPDATERV( VEC )
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
      REAL(SZ), INTENT(INOUT) ::  VEC(MNP,dofh)
      INTEGER :: N,I,J,K,NCOUNT,NFINI,TOT
C
                             !..Pack Messages
      DO J=1,NEIGHPROC
         NCOUNT = 0
         DO I=1,NNODSEND(J)
            DO K=1,dofh
               NCOUNT = NCOUNT+1
               SENDBUF(NCOUNT,J)=VEC(ISENDLOC(I,J),K)
            ENDDO
         ENDDO
      ENDDO
C                    
              ! Send/receive messages to/from all neighbors
C
      CALL MPI_STARTALL ( RDIM, REQ_R3D, IERR )
C
              !..Unpack Received messages as they arrive     
C
      TOT = 0
      DO WHILE (TOT.LT.RDIM)
         DO N=1, RDIM
            INDEX(N) = 0
         ENDDO
         CALL MPI_WAITSOME( RDIM,REQ_R3D,NFINI,INDEX,STAT_R3D,IERR )
         TOT = TOT + NFINI
         DO N=1, NFINI
            IF (INDEX(N).GT.0.AND.INDEX(N).LE.RDIM)  THEN
              IF (INDEX(N).LE.NEIGHPROC) THEN
                J = INDEX(N)
                NCOUNT = 0
                DO I=1,NNODRECV(J)
                   DO K=1,dofh
                      NCOUNT = NCOUNT+1
                      VEC(IRECVLOC(I,J),K) = RECVBUF(NCOUNT,J)
                   ENDDO
                ENDDO
              ENDIF
            ENDIF
         ENDDO
      ENDDO
C 
      RETURN
      END SUBROUTINE

      END MODULE MESSENGER
       


c$$$C------------------------------------------------------------------------------
c$$$C               S U B R O U T I N E   Error  E L E V  S U M 
c$$$C------------------------------------------------------------------------------
c$$$C
c$$$      SUBROUTINE ErrorElevSum( ErrorElevExceeded )
c$$$#ifndef HAVE_MPI_MOD
c$$$      include 'mpif.h'
c$$$#endif
c$$$      INTEGER ErrorElevExceeded !=1 if this subdomain has exceeded warning elev
c$$$      INTEGER SumErrorElevExceeded !sum total of all flags from all subdomains
c$$$      INTEGER kount       ! to avoid compiler bug on certain platforms
c$$$
c$$$      SumErrorElevExceeded = 0
c$$$      kount=1
c$$$      call MPI_ALLREDUCE( ErrorElevExceeded, SumErrorElevExceeded, kount,
c$$$     &     MPI_INTEGER, MPI_SUM, COMM, ierr)
c$$$      ErrorElevExceeded = SumErrorElevExceeded
c$$$      END SUBROUTINE ErrorElevSum
c$$$
c$$$
c$$$
c$$$      END MODULE MESSENGER
       
