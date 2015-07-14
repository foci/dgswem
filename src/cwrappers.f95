#ifdef HPX
subroutine dgswem_init_fort(sizes_c_ptr,dg_c_ptr,global_c_ptr,nodalattr_c_ptr,n_timesteps,n_domains,id)
  use, intrinsic :: iso_c_binding
  use sizes
  use dg
  use global
  use NodalAttributes
  implicit none

  type (C_PTR) :: sizes_c_ptr
  type (C_PTR) :: dg_c_ptr
  type (C_PTR) :: global_c_ptr
  type (C_PTR) :: nodalattr_c_ptr
  integer :: n_timesteps
  integer :: id
  integer :: n_domains

  type (sizes_type), pointer :: s
  type (dg_type), pointer :: dg_here
  type (global_type), pointer :: global_here
  type (nodalattr_type), pointer :: nodalattr_here

  allocate(s)
  allocate(dg_here)
  allocate(global_here)
  allocate(nodalattr_here)

  s%myproc = id;
  s%mnproc = n_domains;

  call dgswem_init(s,dg_here,global_here,nodalattr_here,n_timesteps)
  
  sizes_c_ptr = C_LOC(s)
  dg_c_ptr = C_LOC(dg_here)
  global_c_ptr = C_LOC(global_here)
  nodalattr_c_ptr = C_LOC(nodalattr_here)

end subroutine dgswem_init_fort

subroutine dg_timestep_fort(sizes_c_ptr,dg_c_ptr,global_c_ptr,nodalattr_c_ptr,timestep)
  use, intrinsic :: iso_c_binding
  use sizes
  use dg
  use global
  use NodalAttributes
  implicit none

  type (C_PTR) :: sizes_c_ptr
  type (C_PTR) :: dg_c_ptr
  type (C_PTR) :: global_c_ptr
  type (C_PTR) :: nodalattr_c_ptr
  integer :: timestep

  type (sizes_type), pointer :: s
  type (dg_type), pointer :: dg_here
  type (global_type), pointer :: global_here
  type (nodalattr_type), pointer :: nodalattr_here

  call C_F_POINTER(sizes_c_ptr,s)
  call C_F_POINTER(dg_c_ptr,dg_here)
  call C_F_POINTER(global_c_ptr,global_here)
  call C_F_POINTER(nodalattr_c_ptr,nodalattr_here)

  call dg_timestep(s,dg_here,global_here,nodalattr_here,timestep)

end subroutine dg_timestep_fort

subroutine get_neighbors_fort(sizes_c_ptr,dg_c_ptr,global_c_ptr,neighbors,num_neighbors)
  use, intrinsic :: iso_c_binding
  use sizes
  use dg
  use global
  implicit none

  type (C_PTR) :: sizes_c_ptr
  type (C_PTR) :: dg_c_ptr
  type (C_PTR) :: global_c_ptr
  integer :: neighbors(MAX_DOMAIN_NEIGHBORS)
  integer :: num_neighbors

  integer :: i

  type (sizes_type), pointer :: s
  type (dg_type), pointer :: dg_here
  type (global_type), pointer :: global_here

  call C_F_POINTER(sizes_c_ptr,s)
  call C_F_POINTER(dg_c_ptr,dg_here)
  call C_F_POINTER(global_c_ptr,global_here)

  call get_neighbors(s,dg_here,global_here,neighbors,num_neighbors)

end subroutine get_neighbors_fort


subroutine hpx_get_elems_fort(dg_c_ptr,neighbor,volume,sendbuf)
  use, intrinsic :: iso_c_binding
  use dg
  use sizes
  implicit none

  type (C_PTR) :: dg_c_ptr
  real(sz) :: sendbuf(MAX_BUFFER_SIZE)
  integer :: volume
  integer :: neighbor

  integer :: i

  type (dg_type), pointer :: dg_here

  call C_F_POINTER(dg_c_ptr,dg_here)

  call hpx_get_elems(dg_here,neighbor,volume,sendbuf)

end subroutine hpx_get_elems_fort

subroutine hpx_put_elems_fort(dg_c_ptr,neighbor,volume,recvbuf)
  use, intrinsic :: iso_c_binding
  use dg
  use sizes
  implicit none

  type (C_PTR) :: dg_c_ptr
  real(sz) :: recvbuf(MAX_BUFFER_SIZE)
  integer :: volume
  integer :: neighbor

  integer :: i

  type (dg_type), pointer :: dg_here

  call C_F_POINTER(dg_c_ptr,dg_here)

  call hpx_put_elems(dg_here,neighbor,volume,recvbuf)

end subroutine hpx_put_elems_fort

#endif
