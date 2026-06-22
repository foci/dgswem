module nctest_util

    use ncdg, only: ncdg_63, ncdg_64

    implicit none

    contains

    subroutine print_metadata_63(file)
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
    end subroutine print_metadata_63

    subroutine print_metadata_64(file)
        type(ncdg_64), intent(in) :: file

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

        print *, "fort.64.nc Variable IDs"
        print *, file%u_vel_varname, file%u_vel_varid
        print *, file%v_vel_varname, file%v_vel_varid
    end subroutine print_metadata_64

end module nctest_util



program nctest

    use ncdg, only : ncdg_63, ncdg_64, ncdg_maxele
    use nctest_util
    use netcdf

    implicit none

    type(ncdg_63) :: my_63
    type(ncdg_64) :: my_64
    type(ncdg_maxele) :: my_maxele
    !type(ncdg_maxvel) :: my_maxvel

    logical, parameter :: write_metadata = .false.

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
    call my_63%init()
    print *, "Creating fort.64.nc"
    call my_64%init()
    print *, "Creating maxele.63.nc"
    call my_maxele%init()
    !print *, "Creating maxvel.64.nc"
    !call my_maxvel%init()
    if (write_metadata) then
        call print_metadata_63(my_63)
        call print_metadata_64(my_64)
    endif
    print *, "Setting metadata"
    call my_63%ncdg_63_set_metadata( &
        nt=nf90_unlimited, &
        np=np, &
        ne=ne, &
        nhy=nhy, &
        nope=nope, &
        neta=neta, &
        max_nvdll=maxval(nvdll), &
        nbou=nbou, &
        nvel=nvel, &
        max_nvell=maxval(nvell), &
        ics=ics &
    )
    call my_64%ncdg_64_set_metadata( &
        nt=nf90_unlimited, &
        np=np, &
        ne=ne, &
        nhy=nhy, &
        nope=nope, &
        neta=neta, &
        max_nvdll=maxval(nvdll), &
        nbou=nbou, &
        nvel=nvel, &
        max_nvell=maxval(nvell), &
        ics=ics &
    )
    call my_maxele%ncdg_maxele_set_metadata( &
        nt=nf90_unlimited, &
        np=np, &
        ne=ne, &
        nhy=nhy, &
        nope=nope, &
        neta=neta, &
        max_nvdll=maxval(nvdll), &
        nbou=nbou, &
        nvel=nvel, &
        max_nvell=maxval(nvell), &
        ics=ics &
    )
    if (write_metadata) then
        call print_metadata_63(my_63)
        call print_metadata_64(my_64)
    endif
    print *, "Closing fort.63.nc"
    call my_63%close()
    print *, "Closing fort.64.nc"
    call my_64%close()
    print *, "Closing maxele.63.nc"
    call my_maxele%close()
    if (write_metadata) then
        call print_metadata_63(my_63)
        call print_metadata_64(my_64)
    endif
    print *, "Opening fort.63.nc"
    call my_63%open(mode=nf90_write)
    print *, "Opening fort.64.nc"
    call my_64%open(mode=nf90_write)
    print *, "Opening maxele.63.nc"
    call my_maxele%open(mode=nf90_write)
    if (write_metadata) then
        call print_metadata_63(my_63)
        call print_metadata_64(my_64)
    endif
    print *, "Writing mesh data"
    call my_63%write_mesh( &
        x=x, &
        y=y, &
        dp=dp, &
        nm=nm, &
        nhy=nhy, &
        ne=ne, &
        nope=nope, &
        neta=neta, &
        ibtypee=ibtypee, &
        nvdll=nvdll, &
        nbdv=nbdv, &
        nbou=nbou, &
        nvel=nvel, &
        ibtype=ibtype, &
        nvell=nvell, &
        nbvv=nbvv &
    )
    call my_64%write_mesh( &
        x=x, &
        y=y, &
        dp=dp, &
        nm=nm, &
        nhy=nhy, &
        ne=ne, &
        nope=nope, &
        neta=neta, &
        ibtypee=ibtypee, &
        nvdll=nvdll, &
        nbdv=nbdv, &
        nbou=nbou, &
        nvel=nvel, &
        ibtype=ibtype, &
        nvell=nvell, &
        nbvv=nbvv &
    )
    call my_maxele%write_mesh( &
        x=x, &
        y=y, &
        dp=dp, &
        nm=nm, &
        nhy=nhy, &
        ne=ne, &
        nope=nope, &
        neta=neta, &
        ibtypee=ibtypee, &
        nvdll=nvdll, &
        nbdv=nbdv, &
        nbou=nbou, &
        nvel=nvel, &
        ibtype=ibtype, &
        nvell=nvell, &
        nbvv=nbvv &
    )
    print *, "Writing step"
    call my_63%ncdg_63_write_step( &
        t=0.0, &
        zeta=[1.0, 2.0, 3.0] &
    )
    call my_64%ncdg_64_write_step( &
        t=0.0, &
        u_vel=[1.0, 2.0, 3.0], &
        v_vel=[4.0, 5.0, 6.0] &
    )
    call my_maxele%ncdg_maxele_write_step( &
        t=0.0, &
        zeta_max=[1.0, 2.0, 3.0] &
    )
    print *, "Writing step"
    call my_63%ncdg_63_write_step( &
        t=1.0, &
        zeta=[4.0, 5.0, 6.0] &
    )
    call my_64%ncdg_64_write_step( &
        t=1.0, &
        u_vel=[7.0, 8.0, 9.0], &
        v_vel=[10.0, 11.0, 12.0] &
    )
    call my_maxele%ncdg_maxele_write_step( &
        t=1.0, &
        zeta_max=[2.0, 2.0, 3.0] &
    )
    print *, "Writing step"
    call my_63%ncdg_63_write_step( &
        t=2.0, &
        zeta=[7.0, 8.0, 9.0] &
    )
    call my_64%ncdg_64_write_step( &
        t=2.0, &
        u_vel=[13.0, 14.0, 15.0], &
        v_vel=[16.0, 17.0, 18.0] &
    )
    call my_maxele%ncdg_maxele_write_step( &
        t=2.0, &
        zeta_max=[2.0, 1.0, 4.0] & ! Bad zeta_max update at position 2
    )
    print *, "Closing fort.63.nc"
    call my_63%close()
    print *, "Closing fort.64.nc"
    call my_64%close()
    print *, "Closing maxele.63.nc"
    call my_maxele%close()
    print *, "Program complete"

end program nctest
