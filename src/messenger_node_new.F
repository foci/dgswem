C******************************************************************************
C PADCIRC VERSION 45.01 09/30/2004                                            *
C  last changes in this file VERSION 44.18                                    *
C  S.Bunya changed this file a bit. 07/13/2005                                *
C****************************************************************************** 
C 
      MODULE MESSENGER_NODE_NEW
C
      USE SIZES
      USE GLOBAL, ONLY: C3D,COMM
      USE DG, ONLY:  Ndom, tsgroupnode
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
C  This module supplies the MPI Message-Passing Interface for PADCIRC.
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
      INTEGER,SAVE,PRIVATE ::  NEIGHPROC, RDIM, IERR
      INTEGER,SAVE,PRIVATE ::  TAG = 100
      LOGICAL,SAVE,PRIVATE,ALLOCATABLE :: RESNODE(:)
C
      INTEGER, PRIVATE, ALLOCATABLE :: IPROC(:), NNODELOC(:),
     &    NNODSEND(:), NNODRECV(:), IBELONGTO(:),ISENDLOC(:,:), 
     &    IRECVLOC(:,:), ISENDBUF(:,:), IRECVBUF(:,:),
     $     IRECVLOCLVL(:,:,:),ISENDLOCLVL(:,:,:),
     $     NNODERECVLVL(:,:),NNODESENDLVL(:,:)
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
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
c
#ifdef CRAY
#ifdef  REAL4
      REALTYPE = MPI_REAL4
      DBLETYPE = MPI_REAL8
#else
      REALTYPE = MPI_REAL8
      DBLETYPE = MPI_REAL8
#endif
#else
#ifdef  REAL4
      REALTYPE = MPI_REAL
      DBLETYPE = MPI_DOUBLE_PRECISION
#else
      REALTYPE = MPI_DOUBLE_PRECISION
      DBLETYPE = MPI_DOUBLE_PRECISION
#endif
#endif
c
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
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      INTEGER :: IDPROC,NLOCAL,I,J,LVL,KK
C
      OPEN(18,FILE=DIRNAME(1:LNAME)//'/'//'fort.18')
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
      ALLOCATE (NNODERECVLVL(NEIGHPROC,NDOM))
      ALLOCATE (NNODESENDLVL(NEIGHPROC,NDOM))
      ALLOCATE( IRECVLOCLVL(MNP,NEIGHPROC,NDOM))
      ALLOCATE( ISENDLOCLVL(MNP,NEIGHPROC,NDOM))
C
      DO J=1,NEIGHPROC
         READ(18,3010) IPROC(J),NNODRECV(J)
         READ(18,1130) (IRECVLOC(I,J), I=1,NNODRECV(J))
C
         DO LVL=1,NDOM
            READ(18,3010) KK,NNODERECVLVL(J,LVL)
            IF (NNODERECVLVL(J,LVL).NE.0) THEN
               READ(18,1130) 
     $              (IRECVLOCLVL(I,J,LVL),I=1,NNODERECVLVL(J,LVL))
            ENDIF
         ENDDO
      ENDDO
C
      ALLOCATE( NNODSEND(NEIGHPROC) )
      ALLOCATE( ISENDLOC(MNP,NEIGHPROC) )
C
      DO J=1,NEIGHPROC
         READ(18,3010) IPROC(J),NNODSEND(J)
         READ(18,1130) (ISENDLOC(I,J), I=1,NNODSEND(J))
         DO LVL=1,NDOM
            READ(18,3010) KK,NNODESENDLVL(J,LVL)
            IF (NNODESENDLVL(J,LVL).NE.0) THEN
               READ(18,1130) 
     $              (ISENDLOCLVL(I,J,LVL),I=1,NNODESENDLVL(J,LVL))
            ENDIF
         ENDDO
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

      SUBROUTINE UPDATER_NODE_MOD( VEC1, lvl )
C
C 
      IMPLICIT NONE
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
C
      INTEGER,  INTENT(IN) ::  lvl
      REAL(SZ), INTENT(INOUT) ::  VEC1(:)
      INTEGER :: N,I,J,JJ,K,NCOUNT,NFINI,TOT,IERR
c
      DO J=1,NEIGHPROC
         NCOUNT = 0
         DO I=1,NNODESENDLVL(J,lvl)
            NCOUNT = NCOUNT+1
            SENDBUF(NCOUNT,J)=VEC1(ISENDLOCLVL(I,J,lvl))
         ENDDO
         if (nnodesendlvl(j,lvl).ne.0) then
            call mpi_send(sendbuf(1,j),ncount,realtype,iproc(j),tag,
     $           MPI_COMM_WORLD,ierr)
         endif
      ENDDO
c
      DO J=1,NEIGHPROC
         if (nnoderecvlvl(j,lvl).ne.0) then
c            write(*,*) myproc,'receiving message from ',iproc_r(j)
            call mpi_recv(recvbuf(1,j),nnoderecvlvl(j,lvl),realtype,
     $           iproc(j),tag,mpi_comm_world,stat_r1(1,j),ierr)
            NCOUNT = 0
            DO I=1,NNODERECVLVL(J,lvl)
               NCOUNT = NCOUNT+1
               VEC1(IRECVLOCLVL(I,J,lvl)) = RECVBUF(NCOUNT,J)
            ENDDO
         ENDIF
      ENDDO
c
      RETURN
      END SUBROUTINE


C
C               S U B R O U T I N E   Error  E L E V  S U M 
C------------------------------------------------------------------------------
C
      SUBROUTINE ErrorElevSum( ErrorElevExceeded )
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      INTEGER ErrorElevExceeded !=1 if this subdomain has exceeded warning elev
      INTEGER SumErrorElevExceeded !sum total of all flags from all subdomains
      INTEGER kount       ! to avoid compiler bug on certain platforms

      SumErrorElevExceeded = 0
      kount=1
      call MPI_ALLREDUCE( ErrorElevExceeded, SumErrorElevExceeded, kount,
     &     MPI_INTEGER, MPI_SUM, MPI_COMM_world, ierr)
      ErrorElevExceeded = SumErrorElevExceeded
      END SUBROUTINE ErrorElevSum

      SUBROUTINE blockit
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      END SUBROUTINE blockit

      SUBROUTINE FindMax( dtemp, ddummy )
#ifndef HAVE_MPI_MOD
      include 'mpif.h'
#endif
      REAL*8 dtemp,ddummy
      INTEGER kount       ! to avoid compiler bug on certain platforms

      kount=1
      call MPI_ALLREDUCE( dtemp, ddummy, kount,
     &     MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_world, ierr)
      
      END SUBROUTINE FindMax

      subroutine update_msgpass_node
c
      implicit none
      integer lvl,j,k
      do j=1,neighproc
         do lvl=1,Ndom
            nnoderecvlvl(j,lvl)=0
            do k=1,nnodrecv(j)
               if (tsgroupnode(irecvloc(k,j),lvl).eq.lvl) then
                  nnoderecvlvl(j,lvl)=nnoderecvlvl(j,lvl)+1
                  irecvloclvl(nnoderecvlvl(j,lvl),j,lvl)=irecvloc(k,j)
               endif
            enddo
c            write(888+myproc,*) 'recv ',j,nnoderecv(j),lvl,nnoderecvlvl(j,lvl),
c     $           (irecvloclvl(k,j,lvl),k=1,nnoderecvlvl(j,lvl))
            nnodesendlvl(j,lvl)=0
            do k=1,nnodsend(j)
               if (tsgroupnode(isendloc(k,j),lvl).eq.lvl) then
                  nnodesendlvl(j,lvl)=nnodesendlvl(j,lvl)+1
                  isendloclvl(nnodesendlvl(j,lvl),j,lvl)=isendloc(k,j)
               endif
            enddo
c            write(888+myproc,*) 'send ',j,lvl,nnodesendlvl(j,lvl),
c     $           (isendloclvl(k,j,lvl),k=1,nnodesendlvl(j,lvl))
         enddo
      enddo
      return
      end subroutine


      END MODULE MESSENGER_NODE_NEW
       
