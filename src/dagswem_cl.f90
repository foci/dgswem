MODULE DAGSWEM_CL
    !! Codelets defining tasks for submission
    
    use dg 
    use global 
    use fstarpu_mod
    use dagswem_state
    use iso_c_binding, only: c_ptr, c_f_pointer

    implicit none 

    contains 

        RECURSIVE SUBROUTINE ADVANCE_STAGE_DISTRIBUTED(buffers, cl_args) bind(C)
            !! Advance a node-level subdomain through a single RK stage. 
            TYPE(C_PTR), VALUE, INTENT(IN) :: BUFFERS
            TYPE(C_PTR), VALUE, INTENT(IN) :: CL_ARGS

            call DGSWEM_STATE_ACTIVATE(buffers)
            
            ! Code to unpack buffers and execute RK stage will go here
            
        end subroutine ADVANCE_STAGE_DISTRIBUTED
END MODULE DAGSWEM_CL