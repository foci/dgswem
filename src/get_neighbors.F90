#ifdef HPX
subroutine get_neighbors(s,dg_here,global_here,neighbors,num_neighbors)
Use sizes
use global
use dg
implicit none

type (sizes_type) :: s
type (dg_type) :: dg_here
type (global_type) :: global_here

integer :: num_send
integer :: num_recv
integer :: num_neighbors
integer :: neighbors(MAX_DOMAIN_NEIGHBORS)

integer :: i
integer :: j
integer :: flag

num_send = dg_here%NEIGHPROC_S

do i=1,num_send
   neighbors(i) = dg_here%IPROC_S(i)
end do

num_neighbors = num_send

num_recv = dg_here%NEIGHPROC_R

do i=1,num_recv

   flag = 0
   do j = 1,num_send
     if (neighbors(j) == dg_here%IPROC_R(i)) then
       flag = 1
       exit
     end if
   end do
   
   if (flag == 0) then
     num_neighbors = num_neighbors + 1
     neighbors(num_neighbors) = dg_here%IPROC_R(i)
   end if
     
end do

end subroutine get_neighbors
#endif
