module netcdf_create_files
!! Module for initializing NetCDF files.

    use global, only : dp, mnbou, mnope, nbdv, nbou, nbvv, ne, neta, nhy, nm, nope, np, nvel, nvdll, nvell, segtype, x, y
    !use dg, only : none
    use netcdf, only : nf90_unlimited
    use ncdg, only : ncdg_63
    
    implicit none
    
    contains
    
    subroutine netcdf_create_files_serial()
    !! Initialize all NetCDF output files for the simulation.
    !!
    !! @warning "Serial Use Only"
    !! This subroutine is not designed for use with MPI. If used in parallel,
    !! it should only be called from a single processor.
    !! @endwarning
        
        implicit none

        type(ncdg_63) :: my_63

        integer, allocatable :: ibtypee(:)
        integer, allocatable :: ibtype(:)

        allocate(ibtypee(nope))
        allocate(ibtype(nbou))

        ibtypee = 0 ! This is valid syntax
        ibtype = segtype
        
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
            mesh_dimsize = 1 &
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
        call my_63%close()

        deallocate(ibtypee)
        deallocate(ibtype)
    end subroutine netcdf_create_files_serial

end module netcdf_create_files
