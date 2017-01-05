#ifdef HPX
subroutine dgswem_init_fort(sizes_c_ptr,dg_c_ptr,global_c_ptr,nodalattr_c_ptr,id,single_domain)
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
  integer,intent(in) :: id
  logical*1,intent(in) :: single_domain

!  integer :: n_timesteps
!  integer,intent(out) :: n_domains
!  integer,intent(out) :: n_rksteps

  type (sizes_type), pointer :: s
  type (dg_type), pointer :: dg_here
  type (global_type), pointer :: global_here
  type (nodalattr_type), pointer :: nodalattr_here

  write(99,*) "Entering dgswem_init_fort, id =", id
  flush(99)

#ifdef VERBOSE
!  print*, "FORTRAN: Entering dgswem_init_fort"
!  print*, "FORTRAN: id =", id
#endif

  allocate(s)
  allocate(dg_here)
  allocate(global_here)
  allocate(nodalattr_here)

  s%myproc = id ! I think this should be before the call to dgswem_init
  s%cpp_timestep = 0
  s%cpp_rkstep = 1
  
   s%cpp_single_domain = single_domain
!  s%cpp_single_domain = .false.
  print*, "cpp_single_domain (from cwrappers)= ", s%cpp_single_domain

  call dgswem_init(s,dg_here,global_here,nodalattr_here)

  ! Pass these variables to the c++ side
#ifdef VERBOSE
!  print*, "FORTRAN: s%myproc = ", s%myproc
!  print*, "FORTRAN: dg_here%nrk = ", dg_here%nrk
!  print*, "FORTRAN: global_here%NT = ", global_here%NT
#endif

!  n_domains = s%mnproc
!  n_rksteps = dg_here%nrk
!  n_timesteps = global_here%NT

  sizes_c_ptr = C_LOC(s)
  dg_c_ptr = C_LOC(dg_here)
  global_c_ptr = C_LOC(global_here)
  nodalattr_c_ptr = C_LOC(nodalattr_here)

!  print*, "FORTRAN dgswem_init, cpp_rkstep = ", s%cpp_rkstep
!  print*, "FORTRAN dgswem_init, cpp_timestep = ", s%cpp_timestep

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
  
!  print*, "FORTRAN: sizes_c_ptr = ", sizes_c_ptr
  
#ifdef VERBOSE
  write(99,*) "Entering dg_hydro_timestep_fort, id =", s%myproc, " timestep = ", timestep, " rkstep = ", rkstep
!  print*, "FORTRAN: Entering dg_hydro_timestep_fort"
!  print*, "FORTRAN: myproc =", s%myproc
!  print*, "FORTRAN: timestep =", timestep
!  print*, "FORTRAN: rkstep =", rkstep
#endif

  call dg_hydro_timestep(s,dg_here,global_here,nodalattr_here,timestep,rkstep)

end subroutine dg_hydro_timestep_fort

subroutine slopelimiter_partA_fort(sizes_c_ptr,dg_c_ptr,global_c_ptr)
  use, intrinsic :: iso_c_binding
  use sizes
  use dg
  use global
  implicit none

  type (C_PTR) :: sizes_c_ptr
  type (C_PTR) :: dg_c_ptr
  type (C_PTR) :: global_c_ptr

  type (sizes_type), pointer :: s
  type (dg_type), pointer :: dg_here
  type (global_type), pointer :: global_here

  call C_F_POINTER(sizes_c_ptr,s)
  call C_F_POINTER(dg_c_ptr,dg_here)
  call C_F_POINTER(global_c_ptr,global_here)
  
!  print*, "FORTRAN: sizes_c_ptr = ", sizes_c_ptr
  
#ifdef VERBOSE
  write(99,*) "Entering slopelimiter_fort"
!  print*, "FORTRAN: Entering dg_hydro_timestep_fort"
!  print*, "FORTRAN: myproc =", s%myproc
!  print*, "FORTRAN: timestep =", timestep
!  print*, "FORTRAN: rkstep =", rkstep
#endif

#ifdef SLOPE5
  call slopelimiter5_partA(s,dg_here,global_here)
#endif
end subroutine slopelimiter_partA_fort

subroutine slopelimiter_partB_fort(sizes_c_ptr,dg_c_ptr,global_c_ptr)
  use, intrinsic :: iso_c_binding
  use sizes
  use dg
  use global
  implicit none

  type (C_PTR) :: sizes_c_ptr
  type (C_PTR) :: dg_c_ptr
  type (C_PTR) :: global_c_ptr

  type (sizes_type), pointer :: s
  type (dg_type), pointer :: dg_here
  type (global_type), pointer :: global_here

  call C_F_POINTER(sizes_c_ptr,s)
  call C_F_POINTER(dg_c_ptr,dg_here)
  call C_F_POINTER(global_c_ptr,global_here)
  
!  print*, "FORTRAN: sizes_c_ptr = ", sizes_c_ptr
  
#ifdef VERBOSE
  write(99,*) "Entering slopelimiter_fort"
!  print*, "FORTRAN: Entering dg_hydro_timestep_fort"
!  print*, "FORTRAN: myproc =", s%myproc
!  print*, "FORTRAN: timestep =", timestep
!  print*, "FORTRAN: rkstep =", rkstep
#endif

#ifdef SLOPE5
  call slopelimiter5_partB(s,dg_here,global_here)
#endif
  IF (global_here%NOLIFA .GE. 2) THEN
     call wetdry(dg_here,global_here)
  ENDIF

end subroutine slopelimiter_partB_fort


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

#ifdef VERBOSE
  write(99,*) "Entering dgswem_timestep_advance_fort, id =", s%myproc, " timestep = ", timestep
!  print*, "FORTRAN: Entering dgswem_timestep_advance_fort"
#endif

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

#ifdef VERBOSE
  write(99,*) "FORTRAN: Entering get_neighbors_fort, id = ", s%myproc
#endif

  print*, 's%myproc = ', s%myproc

  call get_neighbors(s,dg_here,global_here,neighbors,num_neighbors)

end subroutine get_neighbors_fort


subroutine hpx_get_elems_fort(dg_c_ptr,neighbor,volume,sendbuf,rkindex)
  use, intrinsic :: iso_c_binding
  use dg
  use sizes
  implicit none

  type (C_PTR) :: dg_c_ptr
  real(sz) :: sendbuf(MAX_BUFFER_SIZE)
  integer :: volume
  integer :: neighbor
  integer :: rkindex

  integer :: i

  type (dg_type), pointer :: dg_here

  call C_F_POINTER(dg_c_ptr,dg_here)

!  print*,'rkindex=',rkindex

#ifdef VERBOSE
  write(99,*) "Entering hpx_get_elems_fort"
#endif

  call hpx_get_elems(dg_here,neighbor,volume,sendbuf,rkindex)

end subroutine hpx_get_elems_fort

subroutine hpx_put_elems_fort(dg_c_ptr,neighbor,volume,recvbuf,rkindex)
  use, intrinsic :: iso_c_binding
  use dg
  use sizes
  implicit none

  type (C_PTR) :: dg_c_ptr
  real(sz) :: recvbuf(MAX_BUFFER_SIZE)
  integer :: volume
  integer :: neighbor
  integer :: rkindex

  integer :: i

  type (dg_type), pointer :: dg_here

  call C_F_POINTER(dg_c_ptr,dg_here)

!  print*,'rkindex=',rkindex

#ifdef VERBOSE
!  print*, "FORTRAN: Entering hpx_put_elems_fort"
  write(99,*) "Entering hpx_put_elems_fort"
#endif

  call hpx_put_elems(dg_here,neighbor,volume,recvbuf,rkindex)

end subroutine hpx_put_elems_fort



subroutine hpx_get_nodes_fort(dg_c_ptr,neighbor,volume,sendbuf)
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


#ifdef VERBOSE
  write(99,*) "Entering hpx_get_nodes_fort"
#endif

  call hpx_get_nodes(dg_here,neighbor,volume,sendbuf)

end subroutine hpx_get_nodes_fort

subroutine hpx_put_nodes_fort(dg_c_ptr,neighbor,volume,recvbuf)
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

#ifdef VERBOSE
!  print*, "FORTRAN: Entering hpx_put_nodes_fort"
  write(99,*) "Entering hpx_put_nodes_fort"
#endif

  call hpx_put_nodes(dg_here,neighbor,volume,recvbuf)

end subroutine hpx_put_nodes_fort




subroutine term_fort(sizes_c_ptr,dg_c_ptr,global_c_ptr,nodalattr_c_ptr)
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

  type (sizes_type), pointer :: s
  type (dg_type), pointer :: dg_here
  type (global_type), pointer :: global_here
  type (nodalattr_type), pointer :: nodalattr_here


  call C_F_POINTER(sizes_c_ptr,s)
  call C_F_POINTER(dg_c_ptr,dg_here)
  call C_F_POINTER(global_c_ptr,global_here)
  call C_F_POINTER(nodalattr_c_ptr,nodalattr_here)

#ifdef VERBOSE
  write(99,*) "FORTRAN: Entering term_fort, id =", s%myproc
#endif

!  deallocate(s)
!  deallocate(dg_here)
!  deallocate(global_here)
!  deallocate(nodalattr_here)

end subroutine term_fort


subroutine cpp_vars_to_fort(sizes_c_ptr,cpp_rkstep,cpp_timestep)

  use, intrinsic :: iso_c_binding
  use sizes

  implicit none
  
  type (C_PTR) :: sizes_c_ptr
  integer :: cpp_rkstep, cpp_timestep
  type (sizes_type), pointer :: s

  call C_F_POINTER(sizes_c_ptr,s)

!  print*, "FORTRAN: cpp_vars_to_fort: rkstep = ", cpp_rkstep, " timestep = ", cpp_timestep
  
  s%cpp_rkstep = cpp_rkstep
  s%cpp_timestep = cpp_timestep  

end subroutine cpp_vars_to_fort

subroutine cpp_vars_from_fort(sizes_c_ptr,cpp_rkstep,cpp_timestep)

  use, intrinsic :: iso_c_binding
  use sizes

  implicit none
  
  type (C_PTR) :: sizes_c_ptr
  integer :: cpp_rkstep, cpp_timestep
  type (sizes_type), pointer :: s

  call C_F_POINTER(sizes_c_ptr,s)
  
  cpp_rkstep = s%cpp_rkstep
  cpp_timestep = s%cpp_timestep

!  print*, "FORTRAN: cpp_vars_from_fort: rkstep = ", cpp_rkstep, " timestep = ", cpp_timestep

end subroutine cpp_vars_from_fort

subroutine check_c_ptr(cptr,ret_val)
 
  use, intrinsic :: iso_c_binding
  use sizes

  implicit none

  type (C_PTR) :: cptr
  integer :: ret_val

!  print*, "FORTRAN: check_c_ptr: cptr = ", cptr

!  if (s.eq.0) then
!     ret_val = 0
!  else
!     ret_val = 1
!  end if

end subroutine check_c_ptr

#endif


