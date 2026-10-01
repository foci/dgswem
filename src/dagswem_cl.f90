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

        RECURSIVE SUBROUTINE ADVANCE_STAGE_DISTRIBUTED(buffers, cl_args) bind(C)
            !! Advance a node-level subdomain through a single RK stage. 
            TYPE(C_PTR), VALUE, INTENT(IN) :: BUFFERS
            TYPE(C_PTR), VALUE, INTENT(IN) :: CL_ARGS
            INTEGER(C_INT), TARGET :: it, irk, sub_id
            REAL(C_DOUBLE), TARGET :: time_a_in
            INTEGER :: j, buf_idx
            REAL(SZ), POINTER :: in_vec(:), out_vec(:)

            CALL fstarpu_unpack_arg(cl_args, (/ C_LOC(it), C_LOC(time_a_in), C_LOC(irk), C_LOC(sub_id) /))

            CALL DGSWEM_STATE_ACTIVATE(buffers)
            MYPROC = sub_id

            ! Unpack ghost cell communication buffers from neighbors (if not first stage of initial timestep)
            IF (irk > 1 .OR. it > ITHS + 1) THEN
                DO j = 1, sub_comm(sub_id)%NEIGHPROC_R
                    buf_idx = NUM_STATE_HANDLES + j - 1
                    CALL c_f_pointer(fstarpu_vector_get_ptr(buffers, buf_idx), in_vec, &
                                     shape=[fstarpu_vector_get_nx(buffers, buf_idx)])
                    CALL UNPACK_COMM_BUFFER(in_vec, sub_comm(sub_id)%NELEMRECV(j), &
                                            sub_comm(sub_id)%IRECVLOC(:, j), irk)
                END DO
            END IF

            ! Execute spatial and stage time advancement
            TIME_A = time_a_in
            CALL DG_HYDRO_TIMESTEP_STAGE(it, irk)

            ! Pack resident element communication buffers for outgoing neighbors
            DO j = 1, sub_comm(sub_id)%NEIGHPROC_S
                buf_idx = NUM_STATE_HANDLES + sub_comm(sub_id)%NEIGHPROC_R + j - 1
                CALL c_f_pointer(fstarpu_vector_get_ptr(buffers, buf_idx), out_vec, &
                                 shape=[fstarpu_vector_get_nx(buffers, buf_idx)])
                CALL PACK_COMM_BUFFER(out_vec, sub_comm(sub_id)%NELEMSEND(j), &
                                      sub_comm(sub_id)%ISENDLOC(:, j), irk + 1)
            END DO

        end subroutine ADVANCE_STAGE_DISTRIBUTED


        RECURSIVE SUBROUTINE FINALIZE_AND_IO(buffers, cl_args) bind(C)
            !! Rotate time levels, convert modal to nodal basis, and write output
            TYPE(C_PTR), VALUE, INTENT(IN) :: BUFFERS
            TYPE(C_PTR), VALUE, INTENT(IN) :: CL_ARGS
            INTEGER(C_INT), TARGET :: it, sub_id
            REAL(C_DOUBLE), TARGET :: time_a_in

            CALL fstarpu_unpack_arg(cl_args, (/ C_LOC(it), C_LOC(time_a_in), C_LOC(sub_id) /))

            CALL DGSWEM_STATE_ACTIVATE(buffers)
            MYPROC = sub_id
            TIME_A = time_a_in

            CALL DG_HYDRO_TIMESTEP_FINALIZE(it)
            CALL modal2nodal()

        end subroutine FINALIZE_AND_IO


        RECURSIVE SUBROUTINE WRITE_OUTPUT_DISTRIBUTED(buffers, cl_args) bind(C)
            !! Write fort.63 independently
            TYPE(C_PTR), VALUE, INTENT(IN) :: BUFFERS
            TYPE(C_PTR), VALUE, INTENT(IN) :: CL_ARGS
            INTEGER(C_INT), TARGET :: it, sub_id
            REAL(C_DOUBLE), TARGET :: time_a_in

            CALL fstarpu_unpack_arg(cl_args, (/ C_LOC(it), C_LOC(time_a_in), C_LOC(sub_id) /))

            CALL DGSWEM_STATE_ACTIVATE(buffers)
            MYPROC = sub_id
            CALL MAKE_DIRNAME()

            CALL WRITE_GLOBAL_ELEVATION(it, time_a_in)
            
        end subroutine WRITE_OUTPUT_DISTRIBUTED

END MODULE DAGSWEM_CL