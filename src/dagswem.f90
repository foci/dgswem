PROGRAM DAGSWEM
   USE SIZES
   USE GLOBAL
   USE DG
   USE DAGSWEM_STATE
   USE DAGSWEM_COMM
   USE DAGSWEM_CL
   USE FSTARPU_MOD
   USE ISO_C_BINDING, ONLY: C_LOC, C_SIZEOF, C_PTR, C_NULL_PTR, C_CHAR, C_NULL_CHAR, C_INT, C_FUNLOC, C_ASSOCIATED

   IMPLICIT NONE

   ! Orchestrator variables
   INTEGER :: I, J, k, ret, IT, stage_idx
   INTEGER :: target_id, sender_id, receiver_id, n_send, buf_len, nbufs
   INTEGER(C_INT), TARGET :: it_c, irk_c, sub_id_c, nbufs_c, nbufs_fin_c
   TYPE SubdomainData
      TYPE(C_PTR) :: handles(NUM_STATE_HANDLES)
   END TYPE
   TYPE(SubdomainData), ALLOCATABLE, TARGET :: subdomains(:)
   LOGICAL :: fileFound

   TYPE(C_PTR) :: cl_advance
   TYPE(C_PTR) :: cl_finalize
   TYPE(C_PTR), ALLOCATABLE, TARGET :: descrs_advance(:)            !! StarPU handle/access info for subdomains
   TYPE(C_PTR), ALLOCATABLE, TARGET :: descrs_finalize(:)
   INTEGER(C_INT), ALLOCATABLE, TARGET :: nbuffers_advance(:)

   ! Find total number of subdomains (MNPROC)
   MNPROC = 0
   DO
      WRITE(DIRNAME(3:6),'(I4.4)') MNPROC
      INQUIRE(FILE=TRIM(DIRNAME)//'/fort.14', EXIST=fileFound)
      IF (.NOT. fileFound) EXIT
      MNPROC = MNPROC + 1
   END DO

   IF (MNPROC == 0) THEN
      PRINT *, "No PE**** directories found."
      STOP
   END IF

   PRINT *, "Found ", MNPROC, " subdomains."

   ALLOCATE(subdomains(0:MNPROC-1))
   ALLOCATE(sub_comm(0:MNPROC-1))

   ! Initialize each subdomain
   DO I = 0, MNPROC - 1
      MYPROC = I
      CALL MAKE_DIRNAME()

      PRINT *, "Initializing Subdomain ", I
      CALL READ_SUBDOMAIN_COMM(I, sub_comm(I))
      CALL COLDSTART()

      ! Save context (pointers are saved in StarPU, and nullified)
      CALL DGSWEM_STATE_REGISTER(subdomains(I)%handles)
   END DO

   ! StarPU initialization
   ret = fstarpu_init(C_NULL_PTR)
   IF (ret /= 0) THEN
      PRINT *, "Error: fstarpu_init failed with code ", ret
      STOP 1
   END IF

   ! Allocate and register send-directed inter-subdomain data handles and buffers
   ALLOCATE(comm_handles(0:MNPROC-1, 0:MNPROC-1))
   ALLOCATE(comm_bufs(0:MNPROC-1, 0:MNPROC-1))
   comm_handles = C_NULL_PTR

   DO I = 0, MNPROC - 1
      DO J = 1, sub_comm(I)%NEIGHPROC_S
         target_id = sub_comm(I)%IPROC_S(J)
         n_send = sub_comm(I)%NELEMSEND(J)
         buf_len = 3 * n_send * DOFH
         ALLOCATE(comm_bufs(I, target_id)%buf(buf_len))
         comm_bufs(I, target_id)%buf = 0.0_SZ
         CALL fstarpu_vector_data_register(comm_handles(I, target_id), 0, &
            C_LOC(comm_bufs(I, target_id)%buf(1)), &
            buf_len, C_SIZEOF(comm_bufs(I, target_id)%buf(1)))
      END DO
   END DO

   ! Setup Codelets
   cl_advance = fstarpu_codelet_allocate()
   CALL fstarpu_codelet_set_name(cl_advance, C_CHAR_"ADVANCE_STAGE_DISTRIBUTED"//C_NULL_CHAR)
   CALL fstarpu_codelet_add_cpu_func(cl_advance, C_FUNLOC(ADVANCE_STAGE_DISTRIBUTED))
   CALL fstarpu_codelet_set_variable_nbuffers(cl_advance)

   cl_finalize = fstarpu_codelet_allocate()
   CALL fstarpu_codelet_set_name(cl_finalize, C_CHAR_"FINALIZE_AND_IO"//C_NULL_CHAR)
   CALL fstarpu_codelet_add_cpu_func(cl_finalize, C_FUNLOC(FINALIZE_AND_IO))
   CALL fstarpu_codelet_set_variable_nbuffers(cl_finalize)

   ! Pre-allocate task buffer descriptor arrays
   ALLOCATE(descrs_advance(0:MNPROC-1))
   ALLOCATE(descrs_finalize(0:MNPROC-1))
   ALLOCATE(nbuffers_advance(0:MNPROC-1))

   DO I = 0, MNPROC - 1
      nbufs = NUM_STATE_HANDLES + sub_comm(I)%NEIGHPROC_R + sub_comm(I)%NEIGHPROC_S
      nbuffers_advance(I) = nbufs
      descrs_advance(I) = fstarpu_data_descr_array_alloc(nbufs)

      ! 1. State handles (RW)
      DO k = 0, NUM_STATE_HANDLES - 1
         CALL fstarpu_data_descr_array_set(descrs_advance(I), k, subdomains(I)%handles(k + 1), FSTARPU_RW)
      END DO

      ! 2. Incoming comm handles (R)
      DO J = 1, sub_comm(I)%NEIGHPROC_R
         sender_id = sub_comm(I)%IPROC_R(J)
         CALL fstarpu_data_descr_array_set(descrs_advance(I), NUM_STATE_HANDLES + J - 1, comm_handles(sender_id, I), FSTARPU_R)
      END DO

      ! 3. Outgoing comm handles (W)
      DO J = 1, sub_comm(I)%NEIGHPROC_S
         receiver_id = sub_comm(I)%IPROC_S(J)
         CALL fstarpu_data_descr_array_set(descrs_advance(I), NUM_STATE_HANDLES + sub_comm(I)%NEIGHPROC_R + J - 1, comm_handles(I, receiver_id), FSTARPU_W)
      END DO

      ! Finalize task handles (RW)
      descrs_finalize(I) = fstarpu_data_descr_array_alloc(NUM_STATE_HANDLES)
      DO k = 0, NUM_STATE_HANDLES - 1
         CALL fstarpu_data_descr_array_set(descrs_finalize(I), k, subdomains(I)%handles(k + 1), FSTARPU_RW)
      END DO
   END DO

   PRINT *, "Starting Task-Based Timestepping: NT = ", NT, ", NRK = ", NRK

   ! Outer Timestep Loop
   DO IT = ITHS + 1, NT
      ! Loop over Runge-Kutta stages
      DO stage_idx = 1, NRK
         DO I = 0, MNPROC - 1
            it_c = IT
            irk_c = stage_idx
            sub_id_c = I
            nbufs_c = nbuffers_advance(I)

            CALL fstarpu_task_insert((/ cl_advance, &
               FSTARPU_VALUE, C_LOC(it_c), FSTARPU_SZ_C_INT, &
               FSTARPU_VALUE, C_LOC(irk_c), FSTARPU_SZ_C_INT, &
               FSTARPU_VALUE, C_LOC(sub_id_c), FSTARPU_SZ_C_INT, &
               FSTARPU_DATA_MODE_ARRAY, descrs_advance(I), C_LOC(nbufs_c), &
               C_NULL_PTR /))
         END DO
      END DO

      ! Finalize and I/O tasks for each subdomain
      DO I = 0, MNPROC - 1
         it_c = IT
         sub_id_c = I
         nbufs_fin_c = NUM_STATE_HANDLES

         CALL fstarpu_task_insert((/ cl_finalize, &
            FSTARPU_VALUE, C_LOC(it_c), FSTARPU_SZ_C_INT, &
            FSTARPU_VALUE, C_LOC(sub_id_c), FSTARPU_SZ_C_INT, &
            FSTARPU_DATA_MODE_ARRAY, descrs_finalize(I), C_LOC(nbufs_fin_c), &
            C_NULL_PTR /))
      END DO
   END DO

   ! Wait for all submitted tasks to complete
   CALL fstarpu_task_wait_for_all()
   PRINT *, "Simulation completed successfully."

   ! Cleanup
   DO I = 0, MNPROC - 1
      CALL fstarpu_data_descr_array_free(descrs_advance(I))
      CALL fstarpu_data_descr_array_free(descrs_finalize(I))
   END DO

   DO I = 0, MNPROC - 1
      DO J = 0, MNPROC - 1
         IF (C_ASSOCIATED(comm_handles(I, J))) THEN
            CALL fstarpu_data_unregister(comm_handles(I, J))
         END IF
      END DO
   END DO

   DO I = 0, MNPROC - 1
      DO k = 1, NUM_STATE_HANDLES
         IF (C_ASSOCIATED(subdomains(I)%handles(k))) THEN
            CALL fstarpu_data_unregister(subdomains(I)%handles(k))
         END IF
      END DO
   END DO

   CALL fstarpu_codelet_free(cl_advance)
   CALL fstarpu_codelet_free(cl_finalize)

   CALL fstarpu_shutdown()

END PROGRAM DAGSWEM
