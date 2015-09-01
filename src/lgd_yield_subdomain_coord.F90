! Helper function for LibGeoDecomp which returns a representative
! coordinate of the subdomain. LibGeoDecomp needs this coordinate for
! its own internal geometric decomposition (i.e. which subdomain get
! mapped to which rank/locality?)
subroutine lgd_yield_subdomain_coord(g, coord)
  use global
  implicit none

  type(global_type), pointer :: g
  real(8), intent(out) :: coord(2)

!  if (g%NE .ne. 0) then
     coord(1) = g%X(1)
     coord(2) = g%Y(1)
!  endif

end subroutine lgd_yield_subdomain_coord
