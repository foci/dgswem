MODULE DAGSWEM_CL
   !! Codelets defining tasks for submission

   use sizes
   use dg
   use global
   use fstarpu_mod
   use dagswem_state
   use dagswem_comm
   use iso_c_binding, only: c_ptr, c_f_pointer, c_loc, c_int, c_double

   implicit none

contains

   RECURSIVE SUBROUTINE ADVANCE_STAGE_INTEGRATE(buffers, cl_args) bind(C)
      !! Unpacks ghost elements (level irk), executes stage integration (IRK -> IRK+1),
      !! and packs resident elements (level irk+1) into intermediate communication buffers.
      TYPE(C_PTR), VALUE, INTENT(IN) :: BUFFERS
      TYPE(C_PTR), VALUE, INTENT(IN) :: CL_ARGS
      INTEGER(C_INT), TARGET :: it, irk, sub_id
      REAL(C_DOUBLE), TARGET :: time_a_in
      INTEGER :: j, buf_idx
      REAL(SZ), POINTER :: in_vec(:), out_vec(:)

      CALL fstarpu_unpack_arg(cl_args, (/ C_LOC(it), C_LOC(time_a_in), C_LOC(irk), C_LOC(sub_id) /))

      CALL DGSWEM_STATE_ACTIVATE(buffers)
      MYPROC = sub_id
      TIME_A = time_a_in

      ! Unpack incoming ghost elements from neighbor subdomains into level irk
      DO j = 1, sub_comm(sub_id)%NEIGHPROC_R
         buf_idx = NUM_STATE_HANDLES + j - 1
         CALL c_f_pointer(fstarpu_vector_get_ptr(buffers, buf_idx), in_vec, &
            shape=[fstarpu_vector_get_nx(buffers, buf_idx)])
         CALL UNPACK_COMM_BUFFER(in_vec, sub_comm(sub_id)%NELEMRECV(j), &
            sub_comm(sub_id)%IRECVLOC(:, j), irk)
      END DO

      ! Execute spatial and stage time integration
      CALL DG_HYDRO_TIMESTEP_STAGE_INTEGRATE(it, irk)

      ! Pack newly updated resident elements (level irk + 1) into outgoing intermediate buffers
      DO j = 1, sub_comm(sub_id)%NEIGHPROC_S
         buf_idx = NUM_STATE_HANDLES + sub_comm(sub_id)%NEIGHPROC_R + j - 1
         CALL c_f_pointer(fstarpu_vector_get_ptr(buffers, buf_idx), out_vec, &
            shape=[fstarpu_vector_get_nx(buffers, buf_idx)])
         CALL PACK_COMM_BUFFER(out_vec, sub_comm(sub_id)%NELEMSEND(j), &
            sub_comm(sub_id)%ISENDLOC(:, j), irk + 1)
      END DO

      ! Pack scalars back into StarPU buffers
      CALL DGSWEM_STATE_DEACTIVATE(buffers)

   end subroutine ADVANCE_STAGE_INTEGRATE


   RECURSIVE SUBROUTINE ADVANCE_STAGE_POSTPROCESS(buffers, cl_args) bind(C)
      !! Unpacks intermediate ghost elements (level irk+1), executes slope limiter
      !! and wetting/drying treatment, packs postprocessed resident elements (level irk+1),
      !! and if final stage (irk == NRK), updates state variables for next timestep.
      TYPE(C_PTR), VALUE, INTENT(IN) :: BUFFERS
      TYPE(C_PTR), VALUE, INTENT(IN) :: CL_ARGS
      INTEGER(C_INT), TARGET :: it, irk, sub_id
      REAL(C_DOUBLE), TARGET :: time_a_in
      INTEGER :: j, buf_idx
      REAL(SZ), POINTER :: in_vec(:), out_vec(:)

      CALL fstarpu_unpack_arg(cl_args, (/ C_LOC(it), C_LOC(time_a_in), C_LOC(irk), C_LOC(sub_id) /))

      CALL DGSWEM_STATE_ACTIVATE(buffers)
      MYPROC = sub_id
      TIME_A = time_a_in

      ! Unpack incoming intermediate ghost elements into level irk + 1
      DO j = 1, sub_comm(sub_id)%NEIGHPROC_R
         buf_idx = NUM_STATE_HANDLES + j - 1
         CALL c_f_pointer(fstarpu_vector_get_ptr(buffers, buf_idx), in_vec, &
            shape=[fstarpu_vector_get_nx(buffers, buf_idx)])
         CALL UNPACK_COMM_BUFFER(in_vec, sub_comm(sub_id)%NELEMRECV(j), &
            sub_comm(sub_id)%IRECVLOC(:, j), irk + 1)
      END DO

      ! Apply post-processing treatments (slope limiter and wet/dry treatment)
      CALL DG_HYDRO_TIMESTEP_STAGE_POSTPROCESS(it, irk)

      ! Pack post-processed resident elements into outgoing buffers for next stage/step
      DO j = 1, sub_comm(sub_id)%NEIGHPROC_S
         buf_idx = NUM_STATE_HANDLES + sub_comm(sub_id)%NEIGHPROC_R + j - 1
         CALL c_f_pointer(fstarpu_vector_get_ptr(buffers, buf_idx), out_vec, &
            shape=[fstarpu_vector_get_nx(buffers, buf_idx)])
         CALL PACK_COMM_BUFFER(out_vec, sub_comm(sub_id)%NELEMSEND(j), &
            sub_comm(sub_id)%ISENDLOC(:, j), irk + 1)
      END DO

      ! At the end of the final Runge-Kutta stage, reset and finalize state variables for next timestep
      IF (irk == NRK) THEN
         CALL DG_HYDRO_TIMESTEP_FINALIZE(it)
      END IF

      ! Pack scalars back into StarPU buffers
      CALL DGSWEM_STATE_DEACTIVATE(buffers)

   end subroutine ADVANCE_STAGE_POSTPROCESS


   RECURSIVE SUBROUTINE WRITE_OUTPUT_DISTRIBUTED(buffers, cl_args) bind(C)
      !! Pure I/O task: unpacks ghost elements into slot 1, converts modal to nodal representation,
      !! and writes fort.63
      TYPE(C_PTR), VALUE, INTENT(IN) :: BUFFERS
      TYPE(C_PTR), VALUE, INTENT(IN) :: CL_ARGS
      INTEGER(C_INT), TARGET :: it, sub_id
      REAL(C_DOUBLE), TARGET :: time_a_in
      INTEGER :: j, buf_idx
      REAL(SZ), POINTER :: in_vec(:)

      CALL fstarpu_unpack_arg(cl_args, (/ C_LOC(it), C_LOC(time_a_in), C_LOC(sub_id) /))

      CALL DGSWEM_STATE_ACTIVATE(buffers)
      MYPROC = sub_id
      TIME_A = time_a_in

      ! Resident elements in slot 1 are already updated by finalize.
      ! Unpack incoming ghost elements into slot 1 from buffer vectors.
      DO j = 1, sub_comm(sub_id)%NEIGHPROC_R
         buf_idx = NUM_STATE_HANDLES + j - 1
         CALL c_f_pointer(fstarpu_vector_get_ptr(buffers, buf_idx), in_vec, &
            shape=[fstarpu_vector_get_nx(buffers, buf_idx)])
         CALL UNPACK_COMM_BUFFER(in_vec, sub_comm(sub_id)%NELEMRECV(j), &
            sub_comm(sub_id)%IRECVLOC(:, j), 1)
      END DO

      CALL modal2nodal()
      CALL MAKE_DIRNAME()
      CALL WRITE_GLOBAL_ELEVATION(it, time_a_in)

      ! Pack scalars back into StarPU buffers
      CALL DGSWEM_STATE_DEACTIVATE(buffers)

   end subroutine WRITE_OUTPUT_DISTRIBUTED

END MODULE DAGSWEM_CL
