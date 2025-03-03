subroutine SwanDispParm ( kwave, cgo, dep2, spcsig )
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
!
!   Updates
!
!   40.80, July 2007: New subroutine
!
!   Purpose
!
!   computes dispersion parameters, wave number and group velocity,
!   in vertices of computational stencil
!
!   Modules used
!
    use ocpcomm4
    use swcomm3
    use SwanGriddata
    use SwanCompdata
!
    implicit none
!
!   Argument variables
!
    real, dimension(MSC,ICMAX), intent(out) :: cgo    ! group velocity
    real, dimension(nverts), intent(in)     :: dep2   ! water depth at current time level
    real, dimension(MSC,ICMAX), intent(out) :: kwave  ! wave number
    real, dimension(MSC), intent(in)        :: spcsig ! relative frequency bins
!
!   Local variables
!
    integer              :: ic       ! loop counter over stencil
    integer, save        :: ient = 0 ! number of entries in this subroutine
    integer              :: is       ! loop counter over frequency bins
    integer              :: ivert    ! vertex index
    !
    real                 :: deploc   ! local depth
    real, dimension(MSC) :: n        ! ratio of group and phase velocity
    real, dimension(MSC) :: nd       ! derivative of N with respect to depth
!
!   Structure
!
!   Description of the pseudo code
!
!   Source text
!
    if (ltrace) call strace (ient,'SwanDispParm')

    do ic = 1, ICMAX
       !
       ivert  = vs(ic)      ! points in computational stencil
       deploc = dep2(ivert)
       !
       if ( deploc > DEPMIN ) then
          !
          call KSCIP1 (MSC, spcsig, deploc, kwave(1,ic), cgo(1,ic), n, nd)
          !
       else
          !
          do is = 1, MSC
             kwave(is,ic) = -1.
             cgo  (is,ic) =  0.
          enddo
          !
       endif
    enddo
    !
end subroutine SwanDispParm
