module nctest_util

    use ncdg, only: ncdg_63
    
    implicit none

    contains

    subroutine print_metadata(file)
        type(ncdg_63), intent(in) :: file

        print *, "Printing metadata:"

        print *, "DGSWEM Dimension IDs"
        print *, file%time_dimname, file%time_dimid
        
        print *, "Nodal Dimension IDs"
        print *, file%node_dimname, file%node_dimid
        print *, file%nele_dimname, file%nele_dimid
        print *, file%nvertex_dimname, file%nvertex_dimid
        print *, file%nope_dimname, file%nope_dimid
        print *, file%neta_dimname, file%neta_dimid
        print *, file%max_nvdll_dimname, file%max_nvdll_dimid
        print *, file%nbou_dimname, file%nbou_dimid
        print *, file%nvel_dimname, file%nvel_dimid
        print *, file%max_nvell_dimname, file%max_nvell_dimid
        print *, file%mesh_dimname, file%mesh_dimid
        
        print *, "DGSWEM Variable IDs"
        print *, file%time_varname, file%time_varid
        
        print *, "Nodal Variable IDs"
        print *, file%x_varname, file%x_varid
        print *, file%y_varname, file%y_varid
        print *, file%element_varname, file%element_varid
        print *, file%adcirc_mesh_varname, file%adcirc_mesh_varid
        print *, file%neta_varname, file%neta_varid
        print *, file%nvdll_varname, file%nvdll_varid
        print *, file%max_nvdll_varname, file%max_nvdll_varid
        print *, file%ibtypee_varname, file%ibtypee_varid
        print *, file%nbdv_varname, file%nbdv_varid
        print *, file%nvel_varname, file%nvel_varid
        print *, file%nvell_varname, file%nvell_varid
        print *, file%max_nvell_varname, file%max_nvell_varid
        print *, file%ibtype_varname, file%ibtype_varid
        print *, file%nbvv_varname, file%nbvv_varid
        print *, file%depth_varname, file%depth_varid

        print *, "fort.63.nc Variable IDs"
        print *, file%zeta_varname, file%zeta_varid
    end subroutine print_metadata
    
end module nctest_util



program nctest

    use ncdg, only : ncdg_63
    use nctest_util
    use netcdf
    
    implicit none
    
    type(ncdg_63) :: my_63

    ! Bare-bones mesh consisting of a single triangle
    integer, parameter :: ics = 1 ! Cartesian
    integer, parameter :: np = 3 ! Number of nodes
    integer, parameter :: ne = 1 ! Number of elements
    integer, parameter :: nhy = 3 ! Number of nodes per element
    integer :: nm(ne, nhy) ! Node numbers for each element

    ! Coordinates and bathymetry
    real :: x(np) = [0.0, 1.0, 1.0]
    real :: y(np) = [0.0, 0.0, 1.0]
    real :: dp(np) = [0.0, 1.0, 2.0]

    ! Elevation boundary
    integer, parameter :: nope = 1 ! Number of elevation segments
    integer, parameter :: neta = 2 ! Number of elevation nodes
    integer :: ibtypee(nope) = [0] ! Segment types
    integer :: nvdll(nope) = [2] ! Nodes per segment
    integer :: nbdv(nope, neta) ! Node numbers for each segment

    ! Normal flow (discharge) boundary
    integer, parameter :: nbou = 1 ! Number of normal flow segments
    integer, parameter :: nvel = 3 ! Number of normal flow nodes
    integer :: ibtype(nbou) = [1] ! Segment types
    integer :: nvell(nbou) = [3] ! Nodes per segment
    integer :: nbvv(nbou, nvel) ! Node numbers for each segment

    ! Fill in arrays
    nm(1, :) = [1, 2, 3]
    nbdv(1, :) = [1, 2]
    nbvv(1, :) = [2, 3, 1]

    ! Test the NetCDF tools
    print *, "Creating fort.63.nc"
    ! call print_metadata(my_63)
    call my_63%init("fort.63.nc")
    print *, "Setting metadata"
    call my_63%ncdg_63_set_metadata( &
        time_dimsize = nf90_unlimited, &
        node_dimsize = np, &
        nele_dimsize = ne, &
        nvertex_dimsize = nhy, &
        nope_dimsize = nope, &
        neta_dimsize = neta, &
        max_nvdll_dimsize = maxval(nvdll), &
        nbou_dimsize = nbou, &
        nvel_dimsize = nvel, &
        max_nvell_dimsize = maxval(nvell), &
        mesh_dimsize = 1, &
        nope = nope, &
        nbou = nbou, &
        ics = ics &
    )
    ! call print_metadata(my_63)
    print *, "Closing fort.63.nc"
    call my_63%close()
    ! call print_metadata(my_63)
    print *, "Opening fort.63.nc"
    call my_63%open(nf90_write)
    ! call print_metadata(my_63)
    print *, "Writing mesh data"
    call my_63%ncdg_nodal_write_mesh( &
        x = x, &
        y = y, &
        dp = dp, &
        nm = nm, &
        nhy = nhy, &
        ne = ne, &
        nope = nope, &
        neta = neta, &
        ibtypee = ibtypee, &
        nvdll = nvdll, &
        nbdv = nbdv, &
        nbou = nbou, &
        nvel = nvel, &
        ibtype = ibtype, &
        nvell = nvell, &
        nbvv = nbvv &
    )
    print *, "Closing fort.63.nc"
    call my_63%close()
    print *, "Program complete"
    
end program nctest
