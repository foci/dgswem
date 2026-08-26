MODULE DAGSWEM_CL
    !! Codelets defining tasks for submission
    
    use dg 
    use global 
    use fstarpu_mod, only: &
    fstarpu_vector_get_ptr, fstarpu_matrix_get_ptr, fstarpu_block_get_ptr, fstarpu_tensor_get_ptr, &
    fstarpu_vector_get_nx, &
    fstarpu_matrix_get_nx, fstarpu_matrix_get_ny, &
    fstarpu_block_get_nx, fstarpu_block_get_ny, fstarpu_block_get_nz, & 
    fstarpu_tensor_get_nx, fstarpu_tensor_get_ny, fstarpu_tensor_get_nz, fstarpu_tensor_get_nw 
    use iso_c_binding, only: c_ptr, c_f_pointer

    implicit none 

    contains 

        RECURSIVE SUBROUTINE ADVANCE_STAGE_DISTRIBUTED(buffers, cl_args) bind(C)
            !! Advance a node-level subdomain through a single RK stage. 
            TYPE(C_PTR), VALUE, INTENT(IN) :: BUFFERS
            TYPE(C_PTR), VALUE, INTENT(IN) :: CL_ARGS

            call dgswem_state_activate(buffers) 
        end subroutine ADVANCE_STAGE_DISTRIBUTED

        SUBROUTINE DGSWEM_STATE_ACTIVATE(state)
            !! Repoint all header-level declarations to the provided
            !! StarPU data handle. 
            !!
            !! Scoped currently to v0.1.0's hydro-only, non-p-adaptive, SLOPE5 (not SLOPEALL)
            !! build: names whose declaration lives inside a currently-disabled
            !! #ifdef block (SLOPEALL, TRACE, CHEM, DYNP, SED_LAY, P_AD, SWAN, CMPI)
            !! are intentionally excluded -- those symbols do not exist in this build.
            !!
            !! Unfortunately, StarPU Fortran does not currently provide a way to 
            !! query the interface of a data handle, so all of these have to be hardcoded and 
            !! must be passed in the exact order specified here. 
            
            TYPE(C_PTR), VALUE, INTENT(IN) :: STATE                                     
            INTEGER, POINTER ::     INT_VEC_STATE_VARS(:),                              &  
                                    INT_MAT_STATE_VARS(:),                              &  
                                    INT_BLOCK_STATE_VARS(:),                            &  
                                    INT_TENSOR_STATE_VARS(:)                                
            REAL(SZ), POINTER ::    REAL_VEC_STATE_VARS(:),                             &  
                                    REAL_MAT_STATE_VARS(:),                             &
                                    REAL_BLOCK_STATE_VARS(:),                           &  
                                    REAL_TENSOR_STATE_VARS(:)                              
            
            integer :: IDX
            
            int_vec_state_vars = [wdflg]
            int_mat_state_vars = []
            int_block_state_vars = []
            int_tensor_state_vars = []
            
            real_vec_state_vars = []
            real_mat_state_vars = []
            real_block_state_vars = []
            real_tensor_state_vars = []

            do idx = 0, size(int_vec_state_vars)-1
                call c_f_pointer(   c_ptr = fstarpu_vector_get_ptr(state, idx),             &            
                                    f_ptr = int_vec_state_vars(idx),                        &           
                                    shape = [fstarpu_vector_get_nx(state, idx)])    
            end do  
            do idx = 0, size(int_mat_state_vars)-1  
                call c_f_pointer(   c_ptr = fstarpu_matrix_get_ptr(state, idx),             &           
                                    f_ptr = int_mat_state_vars(idx),                        &            
                                    shape = [fstarpu_vector_get_nx(state, idx),             &
                                             fstarpu_vector_get_ny(state, idx)])    
            end do  
            do idx = 0, size(int_block_state_vars)-1    
                call c_f_pointer(   c_ptr = fstarpu_block_get_ptr(state, idx),              &          
                                    f_ptr = int_block_state_vars(idx),                      &          
                                    shape = [fstarpu_vector_get_nx(state, idx),             &
                                             fstarpu_vector_get_ny(state, idx),             &
                                             fstarpu_vector_get_nz(state, idx)])    
            end do  
            do idx = 0, size(int_tensor_state_vars)-1   
                call c_f_pointer(   c_ptr = fstarpu_tensor_get_ptr(state, idx),             &                                    
                                    f_ptr = int_tensor_state_vars(idx),                     &                                   
                                    shape = [fstarpu_vector_get_nx(state, idx),             &
                                             fstarpu_vector_get_ny(state, idx),             &
                                             fstarpu_vector_get_nz(state, idx),             &
                                             fstarpu_vector_get_nw(state, idx)])    
            end do
            do idx = 0, size(int_vec_state_vars)-1
                call c_f_pointer(   c_ptr = fstarpu_vector_get_ptr(state, idx),             &       
                                    f_ptr = int_vec_state_vars(idx),                        &      
                                    shape = [fstarpu_vector_get_nx(state, idx)])    
            end do  
            do idx = 0, size(int_mat_state_vars)-1  
                call c_f_pointer(   c_ptr = fstarpu_matrix_get_ptr(state, idx),             &     
                                    f_ptr = int_mat_state_vars(idx),                        &    
                                    shape = [fstarpu_vector_get_nx(state, idx),&
                                             fstarpu_vector_get_ny(state, idx)])    
            end do  
            do idx = 0, size(int_block_state_vars)-1    
                call c_f_pointer(   c_ptr = fstarpu_block_get_ptr(state, idx),              &       
                                    f_ptr = int_block_state_vars(idx),                      &      
                                    shape = [fstarpu_vector_get_nx(state, idx),             &
                                             fstarpu_vector_get_ny(state, idx),             &
                                             fstarpu_vector_get_nz(state, idx)])    
            end do  
            do idx = 0, size(int_tensor_state_vars)-1   
                call c_f_pointer(   c_ptr = fstarpu_tensor_get_ptr(state, idx),             &       
                                    f_ptr = int_tensor_state_vars(idx),                     &      
                                    shape = [fstarpu_vector_get_nx(state, idx),             &
                                             fstarpu_vector_get_ny(state, idx),             &
                                             fstarpu_vector_get_nz(state, idx),             &
                                             fstarpu_vector_get_nw(state, idx)])    
            end do
            do idx = 0, size(real_vec_state_vars)-1
                call c_f_pointer(   c_ptr = fstarpu_vector_get_ptr(state, idx),             &       
                                    f_ptr = real_vec_state_vars(idx),                       &       
                                    shape = [fstarpu_vector_get_nx(state, idx)])    
            end do  
            do idx = 0, size(real_mat_state_vars)-1  
                call c_f_pointer(   c_ptr = fstarpu_matrix_get_ptr(state, idx),             &     
                                    f_ptr = real_mat_state_vars(idx),                       &     
                                    shape = [fstarpu_vector_get_nx(state, idx),             &
                                             fstarpu_vector_get_ny(state, idx)])    
            end do  
            do idx = 0, size(real_block_state_vars)-1    
                call c_f_pointer(   c_ptr = fstarpu_block_get_ptr(state, idx),              &       
                                    f_ptr = real_block_state_vars(idx),                     &       
                                    shape = [fstarpu_vector_get_nx(state, idx),             &
                                             fstarpu_vector_get_ny(state, idx),             &
                                             fstarpu_vector_get_nz(state, idx)])    
            end do  
            do idx = 0, size(real_tensor_state_vars)-1   
                call c_f_pointer(   c_ptr = fstarpu_tensor_get_ptr(state, idx),             &       
                                    f_ptr = real_tensor_state_vars(idx),                    &       
                                    shape = [fstarpu_vector_get_nx(state, idx),             &
                                             fstarpu_vector_get_ny(state, idx),             &
                                             fstarpu_vector_get_nz(state, idx),             &
                                             fstarpu_vector_get_nw(state, idx)])    
            end do              
        end subroutine dgswem_state_activate