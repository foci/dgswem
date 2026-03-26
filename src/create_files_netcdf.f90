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
        
        call create_fort_63()
        
    end subroutine create_netcdf_files_serial

    subroutine create_fort_63()
    
        implicit none
        
        integer :: stat ! Status of each NetCDF function call
        integer :: ncid ! ID of current NetCDF file being operated on
        integer :: dimid ! ID of current dimension being operated on
        integer :: varid ! ID of current variable being operated on
        
        ! fort.63.nc: elevation time series
        stat = nf90_create('fort.63.nc', nf90_clobber, ncid)
        if (stat.neq.nf90_noerr) call handle_err(stat)
        stat = nf90_def_dim(ncid, "time", nf90_unlimited, dimid)
        if (stat.neq.nf90_noerr) call handle_err(stat)
        stat = nf90_def_dim(ncid, "node", ???, dimid)
        if (stat.neq.nf90_noerr) call handle_err(stat)
        stat = nf90_def_dim(ncid, "nele", ???, dimid)
        if (stat.neq.nf90_noerr) call handle_err(stat)
        stat = nf90_def_dim(ncid, "nvertex", ???, dimid)
        if (stat.neq.nf90_noerr) call handle_err(stat)
        stat = nf90_def_dim(ncid, "nope", ???, dimid)
        if (stat.neq.nf90_noerr) call handle_err(stat)
        stat = nf90_def_dim(ncid, "neta", ???, dimid)
        if (stat.neq.nf90_noerr) call handle_err(stat)
        stat = nf90_def_dim(ncid, "max_nvdll", ???, dimid)
        if (stat.neq.nf90_noerr) call handle_err(stat)
        stat = nf90_def_dim(ncid, "nbou", ???, dimid)
        if (stat.neq.nf90_noerr) call handle_err(stat)
        stat = nf90_def_dim(ncid, "nvel", ???, dimid)
        if (stat.neq.nf90_noerr) call handle_err(stat)
        stat = nf90_def_dim(ncid, "max_nvell", ???, dimid)
        if (stat.neq.nf90_noerr) call handle_err(stat)
        stat = nf90_def_dim(ncid, "mesh", ???, dimid)
        if (stat.neq.nf90_noerr) call handle_err(stat)
        stat = nf90_def_var(ncid, "time", nf90_double, &
            (/  /), varid)
        if (stat.neq.nf90_noerr) call handle_err(stat)
        stat = nf90_put_att
        if (stat.neq.nf90_noerr) call handle_err(stat)
        stat = nf90_enddef(ncid)
        if (stat.neq.nf90_noerr) call handle_err(stat)
        stat = nf90_close(ncid)
        if (stat.neq.nf90_noerr) call handle_err(stat)
    end subroutine create_fort_63()

end module create_files_netcdf
