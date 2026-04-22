module netcdf_create_files
!! Module for initializing NetCDF files.

    !use global, only : none
    !use dg, only : none
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
            time_dimsize = 1, &
            node_dimsize = 1, &
            nele_dimsize = 1, &
            nvertex_dimsize = 1, &
            nope_dimsize = 1, &
            neta_dimsize = 1, &
            max_nvdll_dimsize = 1, &
            nbou_dimsize = 1, &
            nvel_dimsize = 1, &
            max_nvell_dimsize = 1, &
            mesh_dimsize = 1 &
        )
        call my_fort_63_nc%close()
    end subroutine netcdf_create_files_serial

end module netcdf_create_files
