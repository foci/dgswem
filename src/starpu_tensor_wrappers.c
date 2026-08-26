#include <starpu.h>
#include <starpu_data_interfaces.h>

/*
 * StarPU's Fortran support library (src/util/fstarpu.c) implements
 * fstarpu_variable_*, fstarpu_vector_*, fstarpu_matrix_*, and fstarpu_block_*,
 * but omitted the implementation for fstarpu_tensor_* despite declaring them
 * in fstarpu_mod.f90.
 *
 * These wrappers provide the exact missing functions matching StarPU's implementation.
 */

void * fstarpu_tensor_get_ptr(void *buffers[], int i)
{
    return (void *)STARPU_TENSOR_GET_PTR(buffers[i]);
}

int fstarpu_tensor_get_nx(void *buffers[], int i)
{
    return STARPU_TENSOR_GET_NX(buffers[i]);
}

int fstarpu_tensor_get_ny(void *buffers[], int i)
{
    return STARPU_TENSOR_GET_NY(buffers[i]);
}

int fstarpu_tensor_get_nz(void *buffers[], int i)
{
    return STARPU_TENSOR_GET_NZ(buffers[i]);
}

int fstarpu_tensor_get_nt(void *buffers[], int i)
{
    return STARPU_TENSOR_GET_NT(buffers[i]);
}
