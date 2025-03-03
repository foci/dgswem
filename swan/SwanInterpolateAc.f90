subroutine SwanInterpolateAc ( acintp, x, y, ac2, excpt )
!
!   --|-----------------------------------------------------------|--
!     | Delft University of Technology                            |
!     | Faculty of Civil Engineering and Geosciences              |
!     | Environmental Fluid Mechanics Section                     |
!     | P.O. Box 5048, 2600 GA  Delft, The Netherlands            |
!     |                                                           |
!     | Programmer: Marcel Zijlema                                |
!   --|-----------------------------------------------------------|--
!
!
!     SWAN (Simulating WAves Nearshore); a third generation wave model
!     Copyright (C) 2008  Delft University of Technology
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!     GNU General Public License for more details.
!
!     A copy of the GNU General Public License is available at
!     http://www.gnu.org/copyleft/gpl.html#SEC3
!     or by writing to the Free Software Foundation, Inc.,
!     59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
!   Authors
!
!   40.80: Marcel Zijlema
!   40.90: Nico Booij
!
!   Updates
!
!   40.80, August 2007: New subroutine
!   40.90,   June 2008: improved interpolation near obstacles
!
!   Purpose
!
!   Interpolates action density to given point
!
!   Method
!
!   Look for closest vertex and determine triangle in which given point is located
!   Determine weighting coefficients for the corresponding vertices
!   Set weighting coeff to zero if there is an obstacle between given point and vertex
!   Interpolate action density using the resulting weighting coefficients
!
!   Modules used
!
    use ocpcomm4
    use swcomm2
    use swcomm3
    use m_obsta
    use SwanGriddata
    use SwanGridobjects
!
    implicit none
!
!   Argument variables
!
    real, dimension(MDC,MSC,nverts), intent(in) :: ac2    ! action density at current time
    real, dimension(MDC,MSC), intent(out)       :: acintp ! interpolated action density
    real, intent(in)                            :: x      ! x-coordinate of given point
    real, intent(in)                            :: y      ! y-coordinate of given point
    logical, intent(out)                        :: excpt  ! if true, value is undefined
!
!   Local variables
!
    integer                               :: icell     ! cell index
    integer, save                         :: ient = 0  ! number of entries in this subroutine
    integer, dimension(3)                 :: ivc       ! vertex indices in cyclic order
    integer                               :: ivert     ! vertex index
    integer                               :: jc        ! loop counter
    integer                               :: k         ! loop counter
    integer                               :: l         ! loop counter
    integer                               :: numcor    ! number of corner points in an obstacle
    integer, dimension(3)                 :: v         ! vertices in present cell
    !
    real                                  :: dxp       ! distance between given point and present vertex in x-direction
    real, dimension (3)                   :: dxv       ! difference of vertices of opposite side in x-coordinate
    real                                  :: dyp       ! distance between given point and present vertex in y-direction
    real, dimension (3)                   :: dyv       ! difference of vertices of opposite side in y-coordinate
    real                                  :: eps       ! a small number
    real                                  :: sumww     ! sum of the interpolation weights
    real                                  :: th        ! direction of given point to present vertex
    real                                  :: th1       ! direction of one face pointing to present vertex
    real                                  :: th2       ! direction of another face pointing to present vertex
    real                                  :: thdiff    ! difference between th and th2
    real                                  :: xb        ! user x-coordinate of begin of obstacle side
    real                                  :: xe        ! user x-coordinate of end of obstacle side
    real, dimension (3)                   :: xv        ! x-coordinate of the vertex
    real                                  :: yb        ! user y-coordinate of begin of obstacle side
    real                                  :: ye        ! user y-coordinate of end of obstacle side
    real, dimension (3)                   :: yv        ! y-coordinate of the vertex
    real, dimension (3)                   :: ww        ! weight of each vertex in the interpolation
    !
    character(80)                         :: msgstr    ! string to pass message
    !
    logical                               :: cellfound ! indicate whether cell containing given point is found or not
    logical, dimension (3)                :: cross     ! if true there is an obstacle between given point and vertex
    logical                               :: EQREAL    ! indicate whether two reals are equal or not
    logical                               :: obstcell  ! if true there is an obstacle in cell
    logical                               :: TCROSS    ! determines whether two line segments cross
    logical                               :: xonobst   ! not used
    !
    type(OBSTDAT), pointer                :: COBST     ! pointer to obstacle data
    !
    type(celltype), dimension(:), pointer :: cell      ! datastructure for cells with their attributes
    type(verttype), dimension(:), pointer :: vert      ! datastructure for vertices with their attributes
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanInterpolateAc')
    !
    ! point to vertex and cell objects
    !
    vert => gridobject%vert_grid
    cell => gridobject%cell_grid
    !
    ! initialize array for interpolated action density
    !
    acintp = 0.
    !
    excpt = .true.
    !
    ! find closest vertex for given point
    !
    call SwanFindPoint ( x, y, ivert )
    !
    ! if point not found, return
    !
    if ( ivert < 0 ) return
    !
    ! determine direction of given point to closest vertex
    !
    dxp = xcugrd(ivert) - x
    dyp = ycugrd(ivert) - y
    th = atan2(dyp,dxp)
    !
    ! if given point equals closest vertex, determine determine action density and return
    !
    if ( EQREAL(dxp,0.) .and. EQREAL(dyp,0.) ) then
       excpt  = .false.
       acintp = ac2(:,:,ivert)
       return
    endif
    !
    cellfound = .false.
    !
    ! loop over cells around closest vertex
    !
    celloop: do jc = 1, vert(ivert)%noc
       !
       ! get cell and its vertices
       !
       icell = vert(ivert)%cell(jc)%atti(CELLID)
       !
       v(1) = cell(icell)%atti(CELLV1)
       v(2) = cell(icell)%atti(CELLV2)
       v(3) = cell(icell)%atti(CELLV3)
       !
       ! get directions of faces to closest vertex
       !
       do k = 1, 3
          if ( v(k) == ivert ) then
             th1 = cell(icell)%geom(k)%th1
             th2 = cell(icell)%geom(k)%th2
             exit
          endif
       enddo
       !
       thdiff = th - th2
       do
          if ( abs(thdiff) <= PI ) exit
          th = th - sign (2., thdiff) * PI
          thdiff = th - th2
       enddo
       !
       ! is given point inside considered cell?
       !
       if ( vert(ivert)%atti(VMARKER) == 1 ) then   ! boundary vertex
          eps = PI/360.
       else
          eps = 0.
       endif
       !
       if ( th > th1-eps .and. th <= th2+eps ) then
          cellfound = .true.
          exit celloop
       endif
       !
    enddo celloop
    !
    ! if cell containing given point not found, give warning and return
    !
    if ( .not.cellfound ) then
       write (msgstr, '(a,f12.4,a,f12.4,a)') ' No triangle containing point (',x+XOFFS,',',y+YOFFS,') is found'
       call msgerr( 1, trim(msgstr) )
       return
    endif
    !
    ! 2D linear interpolation on considered triangle is carried out now
    !
    excpt = .false.
    !
    !  get coordinates of the vertices
    !
    do k = 1, 3
       xv(k) = xcugrd(v(k))
       yv(k) = ycugrd(v(k))
       cross(k) = .false.
    enddo
    !
    ! determine difference in x and y of opposite side
    !
    do k = 1, 3
       ivc(2) = mod(k  ,3)+1
       ivc(3) = mod(k+1,3)+1
       dxv(k) = xv(ivc(3)) - xv(ivc(2))
       dyv(k) = yv(ivc(3)) - yv(ivc(2))
    enddo
    !
    ! determine whether there is an obstacle between given point and vertices
    !
    if ( NUMOBS > 0 ) then
       !
       COBST => FOBSTAC
       !
       do jc = 1, NUMOBS
          !
          numcor = COBST%NCRPTS
          if ( ITEST >= 120 ) write (PRINTF,10) jc, numcor
          !
          xb = COBST%XCRP(1)
          yb = COBST%YCRP(1)
          if ( ITEST >= 120 ) write (PRINTF,20) 1, xb+XOFFS, yb+YOFFS
          !
          do l = 2, numcor
             !
             xe = COBST%XCRP(l)
             ye = COBST%YCRP(l)
             if ( ITEST >= 120 ) write (PRINTF,20) l, xe+XOFFS, ye+YOFFS
             !
             ! loop over vertices
             !
             do k = 1, 3
                if ( TCROSS(x, xv(k), xb, xe, y, yv(k), yb, ye, xonobst) ) cross(k) = .true.
             enddo
             !
             xb = xe
             yb = ye
             !
          enddo
          !
          if (.not.associated(COBST%NEXTOBST)) exit
          COBST => COBST%NEXTOBST
          !
       enddo
       !
    endif
    !
    ! determine weighting coefficients
    !
    obstcell = .false.
    do k = 1, 3
       if (cross(k)) then
          ww(k) = 0.
          obstcell = .true.
       else
          ivc(1) = k
          ivc(2) = mod(k  ,3)+1
          ivc(3) = mod(k+1,3)+1
          ww(k) = ((x - xv(ivc(3))) * dyv(ivc(1)) - (y - yv(ivc(3))) * dxv(ivc(1))) / ( dxv(ivc(2)) * dyv(ivc(1)) - dyv(ivc(2)) * dxv(ivc(1)) )
       endif
    enddo
    if (obstcell) sumww = sum(ww)
    !
    ! use weighting coefficients to determine interpolated action density
    !
    do k = 1, 3
       if ( ww(k) > 1.e-10 ) then
          if (obstcell) ww(k) = ww(k) / sumww
          acintp(:,:) = acintp(:,:) + ww(k) * ac2(:,:,v(k))
       endif
    enddo
    !
 10 format (' Obstacle number : ', i4,'  has ', i4, ' corners')
 20 format (' Corner number:', i4,'    Xp: ', e10.4, ' Yp: ', e11.4)
    !
end subroutine SwanInterpolateAc
