module netcdf_io
!! Module for interacting with NetCDF files.

    !use dg, only : none
    use netcdf, only : nf90_unlimited
    use ncdg, only : ncdg_63, ncdg_64
    
    implicit none

    ! Committing the sin of global variables, for now
    type(ncdg_63) :: my_63
    type(ncdg_64) :: my_64
    
    contains
    
    subroutine netcdf_create_files_serial()
        !! Initialize all NetCDF output files for the simulation.
        !!
        !! @warning "Serial Use Only"
        !! This subroutine is not designed for use with MPI. If used in parallel,
        !! it should only be called from a single processor.
        !! @endwarning
        
        use global, only : dp, ics, mnbou, mnope, nbdv, nbou, nbvv, ne, neta, &
            nhy, nm, nope, np, nvel, nvdll, nvell, segtype, x, y
        
        implicit none

        integer, allocatable :: ibtypee(:)
        integer, allocatable :: ibtype(:)

        allocate(ibtypee(nope))
        allocate(ibtype(nbou))

        ibtypee = 0 ! This is valid syntax for array of 0s
        ibtype = segtype

        ! fort.63.nc
        call my_63%init("fort.63.nc")
        call my_63%ncdg_63_set_metadata( &
            time_dimsize = nf90_unlimited, &
            node_dimsize = np, &
            nele_dimsize = ne, &
            nvertex_dimsize = 3, &
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

        ! fort.64.nc
        call my_64%init("fort.64.nc")
        call my_64%ncdg_64_set_metadata( &
            time_dimsize = nf90_unlimited, &
            node_dimsize = np, &
            nele_dimsize = ne, &
            nvertex_dimsize = 3, &
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
        call my_64%ncdg_nodal_write_mesh( &
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

        deallocate(ibtypee)
        deallocate(ibtype)
    end subroutine netcdf_create_files_serial

    subroutine netcdf_write_step_serial()
        !! Write current state to all NetCDF output files for the simulation.

        use global, only : eta2, time_a, uu2, vv2

        implicit none
        
        call my_63%ncdg_63_write_step(t=time_a, zeta=eta2)
        call my_64%ncdg_64_write_step(t=time_a, u_vel=uu2, v_vel=vv2)
    end subroutine netcdf_write_step_serial

    subroutine netcdf_close_files_serial()
        !! Close all NetCDF output files
        
        call my_63%close()
        call my_64%close()
    end subroutine netcdf_close_files_serial

end module netcdf_io
