module netcdf_create_files
!! Module for initializing NetCDF files.

    use global, only : nbou, ne, neta, nope, np, nvel, nvell, nvdll
    !use dg, only : none
    use netcdf, only : nf90_unlimited
    use netcdf_fort, only : fort_63_nc
    
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

        type(fort_63_nc) :: my_fort_63_nc
        
        call my_fort_63_nc%init("fort.63.nc")
        call my_fort_63_nc%set_metadata( &
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
        call my_fort_63_nc%close()
    end subroutine netcdf_create_files_serial

end module netcdf_create_files
