program dgswem_hpx
  use global_data
  implicit none
  
  type (global_data_type) :: g

  call init_hpx(g)
  
  do itime_a = iths+1,nt
     call dg_timestep(itime_a,g)
  end do

