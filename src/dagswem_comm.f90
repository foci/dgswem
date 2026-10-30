MODULE DAGSWEM_COMM
   USE SIZES
   USE GLOBAL
   USE DG, ONLY: DOFH, ZE, QX, QY, ZE_MIN1, ZE_MAX1, QX_MIN1, QX_MAX1, QY_MIN1, QY_MAX1
   USE FSTARPU_MOD
   USE ISO_C_BINDING
   IMPLICIT NONE

   TYPE SubdomainComm
      !! Topological information for a subdomain (elements)
      INTEGER :: IDPROC                                       !! Subdomain ID
      INTEGER :: NLOCAL                                       !! Number of elements in the subdomain (including ghost elements)
      INTEGER, ALLOCATABLE :: NELEMLOC(:)                     !! Local element IDs
      INTEGER :: NEIGHPROC_R                                  !! Number of receiving neighbor subdomains
      INTEGER :: NEIGHPROC_S                                  !! Number of sending neighbor subdomains
      INTEGER, ALLOCATABLE :: IPROC_R(:)                      !! IDs of receiving neighbor subdomains
      INTEGER, ALLOCATABLE :: NELEMRECV(:)                    !! Number of elements to receive from neighbors
      INTEGER, ALLOCATABLE :: IRECVLOC(:,:)                   !! Elements received from neighbors (list/jagged array)
      INTEGER, ALLOCATABLE :: IPROC_S(:)                      !! IDs of sending neighbor subdomains
      INTEGER, ALLOCATABLE :: NELEMSEND(:)                    !! Number of elements to send to each neighbor
      INTEGER, ALLOCATABLE :: ISENDLOC(:,:)                   !! Elements sent to neighbors (list/jagged array)
   END TYPE SubdomainComm

   TYPE SubdomainNodeComm
      !! Topological information for a subdomain (nodes)
      INTEGER :: IDPROC                                       !! Subdomain ID
      INTEGER :: NLOCAL                                       !! Number of resident nodes
      INTEGER, ALLOCATABLE :: NNODELOC(:)                     !! Local resident node IDs
      INTEGER :: NEIGHPROC_R                                  !! Number of receiving neighbor subdomains
      INTEGER :: NEIGHPROC_S                                  !! Number of sending neighbor subdomains
      INTEGER, ALLOCATABLE :: IPROC_R(:)                      !! IDs of receiving neighbor subdomains
      INTEGER, ALLOCATABLE :: NNODRECV(:)                     !! Number of nodes to receive from neighbors
      INTEGER, ALLOCATABLE :: IRECVLOC(:,:)                   !! Nodes received from neighbors (list/jagged array)
      INTEGER, ALLOCATABLE :: IPROC_S(:)                      !! IDs of sending neighbor subdomains
      INTEGER, ALLOCATABLE :: NNODSEND(:)                     !! Number of nodes to send to each neighbor
      INTEGER, ALLOCATABLE :: ISENDLOC(:,:)                   !! Nodes sent to neighbors (list/jagged array)
   END TYPE SubdomainNodeComm

   TYPE CommBuffer
      !! A buffer to store ghost data
      REAL(SZ), ALLOCATABLE :: buf(:)
   END TYPE CommBuffer

   TYPE(SubdomainComm), ALLOCATABLE, TARGET :: sub_comm(:)
   TYPE(SubdomainNodeComm), ALLOCATABLE, TARGET :: sub_node_comm(:)
   TYPE(C_PTR), ALLOCATABLE, TARGET :: comm_handles(:,:)      !! StarPU handles for final/postprocessed data
   TYPE(CommBuffer), ALLOCATABLE, TARGET :: comm_bufs(:,:)    !! Write buffer for final/postprocessed data
   TYPE(C_PTR), ALLOCATABLE, TARGET :: comm_handles_inter(:,:)!! StarPU handles for intermediate data
   TYPE(CommBuffer), ALLOCATABLE, TARGET :: comm_bufs_inter(:,:)!! Write buffer for intermediate data
   TYPE(C_PTR), ALLOCATABLE, TARGET :: comm_handles_node(:,:) !! StarPU handles for nodal extrema data
   TYPE(CommBuffer), ALLOCATABLE, TARGET :: comm_bufs_node(:,:)  !! Write buffer for nodal extrema data

CONTAINS

   SUBROUTINE READ_SUBDOMAIN_COMM(sub_id, comm_data)
      !! Read the topological information of a subdomain
      INTEGER, INTENT(IN) :: sub_id
      TYPE(SubdomainComm), INTENT(OUT) :: comm_data
      INTEGER :: I, J, JJ, IDPROC, NLOCAL, NEIGHPROC_R, NEIGHPROC_S
      CHARACTER(256) :: dg18_path
      INTEGER :: iu
      LOGICAL :: exists

      WRITE(dg18_path, '(A,I4.4,A)') 'PE', sub_id, '/DG.18'
      INQUIRE(FILE=TRIM(dg18_path), EXIST=exists)
      IF (.NOT. exists) THEN
         PRINT *, "Error: DG.18 not found at ", TRIM(dg18_path)
         STOP 1
      END IF

      iu = 180 + sub_id
      OPEN(iu, FILE=TRIM(dg18_path), STATUS='OLD')

      READ(iu, 3010) IDPROC, NLOCAL
      comm_data%IDPROC = IDPROC
      comm_data%NLOCAL = NLOCAL

      ALLOCATE(comm_data%NELEMLOC(NLOCAL))
      READ(iu, 1130) (comm_data%NELEMLOC(I), I=1,NLOCAL)

      READ(iu, 3010) NEIGHPROC_R, NEIGHPROC_S
      comm_data%NEIGHPROC_R = NEIGHPROC_R
      comm_data%NEIGHPROC_S = NEIGHPROC_S

      ALLOCATE(comm_data%IPROC_R(NEIGHPROC_R), comm_data%NELEMRECV(NEIGHPROC_R))
      ALLOCATE(comm_data%IRECVLOC(MNE, NEIGHPROC_R))

      DO JJ = 1, NEIGHPROC_R
         J = MOD(JJ - 1 + sub_id, NEIGHPROC_R) + 1
         READ(iu, 3010) comm_data%IPROC_R(J), comm_data%NELEMRECV(J)
         READ(iu, 1130) (comm_data%IRECVLOC(I, J), I=1,comm_data%NELEMRECV(J))
      END DO

      ALLOCATE(comm_data%IPROC_S(NEIGHPROC_S), comm_data%NELEMSEND(NEIGHPROC_S))
      ALLOCATE(comm_data%ISENDLOC(MNE, NEIGHPROC_S))

      DO JJ = 1, NEIGHPROC_S
         J = MOD(JJ - 1 + sub_id, NEIGHPROC_S) + 1
         READ(iu, 3010) comm_data%IPROC_S(J), comm_data%NELEMSEND(J)
         READ(iu, 1130) (comm_data%ISENDLOC(I, J), I=1,comm_data%NELEMSEND(J))
      END DO

      CLOSE(iu)

1130  FORMAT(8X,9I8)
3010  FORMAT(8X,2I8)
   END SUBROUTINE READ_SUBDOMAIN_COMM


   SUBROUTINE UNPACK_COMM_BUFFER(comm_vec, n_elems, recv_loc, stage)
      !! Unpack state data from ghost elements for a RK stage
      REAL(SZ), INTENT(IN) :: comm_vec(:)           !! Buffer containing data from neighbor subdomains
      INTEGER, INTENT(IN) :: n_elems                !! Number of elements to unpack
      INTEGER, INTENT(IN) :: recv_loc(:)            !! Local indices of the elements to unpack
      INTEGER, INTENT(IN) :: stage                  !! RK stage
      INTEGER :: I, K, ncount

      ncount = 0
      ! Unpack ZE
      DO I = 1, n_elems
         DO K = 1, DOFH
            ncount = ncount + 1
            ZE(K, recv_loc(I), stage) = comm_vec(ncount)
         END DO
      END DO

      ! Unpack QX
      DO I = 1, n_elems
         DO K = 1, DOFH
            ncount = ncount + 1
            QX(K, recv_loc(I), stage) = comm_vec(ncount)
         END DO
      END DO

      ! Unpack QY
      DO I = 1, n_elems
         DO K = 1, DOFH
            ncount = ncount + 1
            QY(K, recv_loc(I), stage) = comm_vec(ncount)
         END DO
      END DO
   END SUBROUTINE UNPACK_COMM_BUFFER


   SUBROUTINE PACK_COMM_BUFFER(comm_vec, n_elems, send_loc, stage)
      !! Pack state data for a neighbor subdomain for a given RK stage
      REAL(SZ), INTENT(OUT) :: comm_vec(:)
      INTEGER, INTENT(IN) :: n_elems                !! Number of elements to pack
      INTEGER, INTENT(IN) :: send_loc(:)            !! Local indices of the elements to pack
      INTEGER, INTENT(IN) :: stage                  !! RK stage
      INTEGER :: I, K, ncount

      ncount = 0
      ! Pack ZE
      DO I = 1, n_elems
         DO K = 1, DOFH
            ncount = ncount + 1
            comm_vec(ncount) = ZE(K, send_loc(I), stage)
         END DO
      END DO

      ! Pack QX
      DO I = 1, n_elems
         DO K = 1, DOFH
            ncount = ncount + 1
            comm_vec(ncount) = QX(K, send_loc(I), stage)
         END DO
      END DO

      ! Pack QY
      DO I = 1, n_elems
         DO K = 1, DOFH
            ncount = ncount + 1
            comm_vec(ncount) = QY(K, send_loc(I), stage)
         END DO
      END DO
   END SUBROUTINE PACK_COMM_BUFFER


   SUBROUTINE READ_SUBDOMAIN_NODE_COMM(sub_id, comm_data)
      !! Read nodal message-passing information from fort.18
      INTEGER, INTENT(IN) :: sub_id
      TYPE(SubdomainNodeComm), INTENT(OUT) :: comm_data
      INTEGER :: I, J, JJ, IDPROC, NLOCAL, NEIGHPROC
      CHARACTER(256) :: f18_path
      INTEGER :: iu
      LOGICAL :: exists

      WRITE(f18_path, '(A,I4.4,A)') 'PE', sub_id, '/fort.18'
      INQUIRE(FILE=TRIM(f18_path), EXIST=exists)
      IF (.NOT. exists) THEN
         PRINT *, "Error: fort.18 not found at ", TRIM(f18_path)
         STOP 1
      END IF

      iu = 190 + sub_id
      OPEN(iu, FILE=TRIM(f18_path), STATUS='OLD')

      READ(iu, 3010) IDPROC, NLOCAL
      comm_data%IDPROC = IDPROC
      comm_data%NLOCAL = NLOCAL

      ALLOCATE(comm_data%NNODELOC(NLOCAL))
      READ(iu, 1130) (comm_data%NNODELOC(I), I=1,NLOCAL)

      READ(iu, 3015) NEIGHPROC
      comm_data%NEIGHPROC_R = NEIGHPROC
      comm_data%NEIGHPROC_S = NEIGHPROC

      ALLOCATE(comm_data%IPROC_R(NEIGHPROC), comm_data%NNODRECV(NEIGHPROC))
      ALLOCATE(comm_data%IRECVLOC(MNP, NEIGHPROC))

      DO JJ = 1, NEIGHPROC
         J = MOD(JJ - 1 + sub_id, NEIGHPROC) + 1
         READ(iu, 3010) comm_data%IPROC_R(J), comm_data%NNODRECV(J)
         READ(iu, 1130) (comm_data%IRECVLOC(I, J), I=1,comm_data%NNODRECV(J))
      END DO

      ALLOCATE(comm_data%IPROC_S(NEIGHPROC), comm_data%NNODSEND(NEIGHPROC))
      ALLOCATE(comm_data%ISENDLOC(MNP, NEIGHPROC))

      DO JJ = 1, NEIGHPROC
         J = MOD(JJ - 1 + sub_id, NEIGHPROC) + 1
         READ(iu, 3010) comm_data%IPROC_S(J), comm_data%NNODSEND(J)
         READ(iu, 1130) (comm_data%ISENDLOC(I, J), I=1,comm_data%NNODSEND(J))
      END DO

      CLOSE(iu)

1130  FORMAT(8X,9I8)
3010  FORMAT(8X,2I8)
3015  FORMAT(8X,I8)
   END SUBROUTINE READ_SUBDOMAIN_NODE_COMM


   SUBROUTINE PACK_NODE_COMM_BUFFER(comm_vec, n_nodes, send_loc)
      !! Pack nodal extrema for a neighbor subdomain
      REAL(SZ), INTENT(OUT) :: comm_vec(:)
      INTEGER, INTENT(IN) :: n_nodes
      INTEGER, INTENT(IN) :: send_loc(:)
      INTEGER :: I, ncount

      ncount = 0
      ! 1. ZE_MIN1
      DO I = 1, n_nodes
         ncount = ncount + 1
         comm_vec(ncount) = ZE_MIN1(send_loc(I))
      END DO
      ! 2. ZE_MAX1
      DO I = 1, n_nodes
         ncount = ncount + 1
         comm_vec(ncount) = ZE_MAX1(send_loc(I))
      END DO
      ! 3. QX_MIN1
      DO I = 1, n_nodes
         ncount = ncount + 1
         comm_vec(ncount) = QX_MIN1(send_loc(I))
      END DO
      ! 4. QX_MAX1
      DO I = 1, n_nodes
         ncount = ncount + 1
         comm_vec(ncount) = QX_MAX1(send_loc(I))
      END DO
      ! 5. QY_MIN1
      DO I = 1, n_nodes
         ncount = ncount + 1
         comm_vec(ncount) = QY_MIN1(send_loc(I))
      END DO
      ! 6. QY_MAX1
      DO I = 1, n_nodes
         ncount = ncount + 1
         comm_vec(ncount) = QY_MAX1(send_loc(I))
      END DO
   END SUBROUTINE PACK_NODE_COMM_BUFFER


   SUBROUTINE UNPACK_NODE_COMM_BUFFER(comm_vec, n_nodes, recv_loc)
      !! Unpack nodal extrema received from a neighbor subdomain
      REAL(SZ), INTENT(IN) :: comm_vec(:)
      INTEGER, INTENT(IN) :: n_nodes
      INTEGER, INTENT(IN) :: recv_loc(:)
      INTEGER :: I, ncount

      ncount = 0
      ! 1. ZE_MIN1
      DO I = 1, n_nodes
         ncount = ncount + 1
         ZE_MIN1(recv_loc(I)) = comm_vec(ncount)
      END DO
      ! 2. ZE_MAX1
      DO I = 1, n_nodes
         ncount = ncount + 1
         ZE_MAX1(recv_loc(I)) = comm_vec(ncount)
      END DO
      ! 3. QX_MIN1
      DO I = 1, n_nodes
         ncount = ncount + 1
         QX_MIN1(recv_loc(I)) = comm_vec(ncount)
      END DO
      ! 4. QX_MAX1
      DO I = 1, n_nodes
         ncount = ncount + 1
         QX_MAX1(recv_loc(I)) = comm_vec(ncount)
      END DO
      ! 5. QY_MIN1
      DO I = 1, n_nodes
         ncount = ncount + 1
         QY_MIN1(recv_loc(I)) = comm_vec(ncount)
      END DO
      ! 6. QY_MAX1
      DO I = 1, n_nodes
         ncount = ncount + 1
         QY_MAX1(recv_loc(I)) = comm_vec(ncount)
      END DO
   END SUBROUTINE UNPACK_NODE_COMM_BUFFER

END MODULE DAGSWEM_COMM
