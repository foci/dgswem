module c_interface

  interface
     subroutine create_adj_list(np, nedges, nedno1, nedno2) bind(C, name='create_adj_list')
       use, intrinsic :: iso_c_binding, only : c_int

       implicit none

       integer(c_int), value :: np, nedges
       integer(c_int) :: nedno1(nedges), nedno2(nedges)

     end subroutine create_adj_list

     integer function get_edge_no(n1, n2) bind(C, name='get_edge_no')
       use, intrinsic :: iso_c_binding, only : c_int
       implicit none

       integer(c_int), value :: n1, n2

     end function get_edge_no

  end interface
end module c_interface
