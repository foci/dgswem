#ifdef HPX
subroutine get_neighbors(s,dg_here,global_here,neighbors,num_neighbors)
use sizes
use global
use dg
implicit none

type (sizes_type) :: s
type (dg_type) :: dg_here
type (global_type) :: global_here

integer :: num_neighbors
integer :: neighbors(MAX_DOMAIN_NEIGHBORS)

integer :: i

num_neighbors = dg_here%NEIGHPROC_S

do i=1,num_neighbors
   neighbors(i) = dg_here%IPROC_S(i)
end do

end subroutine get_neighbors
#endif
