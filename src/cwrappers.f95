#ifdef HPX
subroutine dgswem_init_fort(sizes_c_ptr,dg_c_ptr,global_c_ptr,nodalattr_c_ptr,id,n_timesteps,n_rksteps)
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
  integer :: n_rksteps

  type (sizes_type), pointer :: s
  type (dg_type), pointer :: dg_here
  type (global_type), pointer :: global_here
  type (nodalattr_type), pointer :: nodalattr_here

  allocate(s)
  allocate(dg_here)
  allocate(global_here)
  allocate(nodalattr_here)

  s%myproc = id ! I think this should be before the call to dgswem_init

  call dgswem_init(s,dg_here,global_here,nodalattr_here)

  ! Pass these variables to the c++ side
  print*, "FORTRAN: s%mnproc = ", s%mnproc
  print*, "FORTRAN: dg_here%nrk = ", dg_here%nrk
  print*, "FORTRAN: global_here%NT = ", global_here%NT

  n_domains = s%mnproc
  n_rksteps = dg_here%nrk
  n_timesteps = global_here%NT

  sizes_c_ptr = C_LOC(s)
  dg_c_ptr = C_LOC(dg_here)
  global_c_ptr = C_LOC(global_here)
  nodalattr_c_ptr = C_LOC(nodalattr_here)

end subroutine dgswem_init_fort

subroutine dg_hydro_timestep_fort(sizes_c_ptr,dg_c_ptr,global_c_ptr,nodalattr_c_ptr,timestep,rkstep)
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
  integer :: rkstep

  type (sizes_type), pointer :: s
  type (dg_type), pointer :: dg_here
  type (global_type), pointer :: global_here
  type (nodalattr_type), pointer :: nodalattr_here

  call C_F_POINTER(sizes_c_ptr,s)
  call C_F_POINTER(dg_c_ptr,dg_here)
  call C_F_POINTER(global_c_ptr,global_here)
  call C_F_POINTER(nodalattr_c_ptr,nodalattr_here)

  call dg_hydro_timestep(s,dg_here,global_here,nodalattr_here,timestep,rkstep)

end subroutine dg_hydro_timestep_fort

SUBROUTINE DG_TIMESTEP_ADVANCE_fort(sizes_c_ptr,dg_c_ptr,global_c_ptr,nodalattr_c_ptr,timestep)
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
  integer :: it

  it = timestep

  call C_F_POINTER(sizes_c_ptr,s)
  call C_F_POINTER(dg_c_ptr,dg_here)
  call C_F_POINTER(global_c_ptr,global_here)
  call C_F_POINTER(nodalattr_c_ptr,nodalattr_here)

  call DG_TIMESTEP_ADVANCE(s,dg_here,global_here,nodalattr_here,IT)

end SUBROUTINE DG_TIMESTEP_ADVANCE_fort

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


subroutine hpx_swap_elems_fort(dg_domain_c_ptr,dg_neighbor_c_ptr)
  use, intrinsic :: iso_c_binding
  use dg
  use sizes
  implicit none

  type (C_PTR) :: dg_domain_c_ptr
  type (C_PTR) :: dg_neighbor_c_ptr
  real(sz) :: recvbuf(MAX_BUFFER_SIZE)
  integer :: volume
  integer :: neighbor

  integer :: i

  type (dg_type), pointer :: dg_here_domain
  type (dg_type), pointer :: dg_here_neighbor

  call C_F_POINTER(dg_domain_c_ptr,dg_here_domain)
  call C_F_POINTER(dg_neighbor_c_ptr,dg_here_neighbor)

  call hpx_swap_elems(dg_here_domain,dg_here_neighbor)

end subroutine hpx_swap_elems_fort

#endif
