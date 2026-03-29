module create_files_netcdf
!! Module for initializing NetCDF files.

    use netcdf
    
    implicit none
    
    contains
    
    subroutine create_netcdf_files_serial()
    !! Initialize all NetCDF output files for the simulation.
    !!
    !! @warning "Serial Use Only"
    !! This subroutine is not designed for use with MPI. If used in parallel,
    !! it should only be called from a single processor.
    !! @endwarning
    
        ! use/import statements go here
        
        implicit none
        
        !! 
        
    end subroutine create_netcdf_files_serial

end module create_files_netcdf
