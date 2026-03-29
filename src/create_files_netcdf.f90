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
        integer :: fort_63_ncid ! ID of fort.63.nc NetCDF file.
        integer :: time_dimid ! ID of time dimension
        integer :: node_dimid ! ID of node dimension
        integer :: nele_dimid ! ID of nele dimension
        integer :: nvertex_dimid ! ID of nvertex dimension
        integer :: nope_dimid ! ID of nope dimension
        integer :: neta_dimid ! ID of neta dimension
        integer :: max_nvdll_dimid ! ID of max_nvdll dimension
        integer :: nbou_dimid ! ID of nbou dimension
        integer :: nvel_dimid ! ID of nvel dimension
        integer :: max_nvell_dimid ! ID of max_nvell dimension
        integer :: mesh_dimid ! ID of mesh dimension
        integer :: time_varid ! ID of time variable
        integer :: x_varid ! ID of x variable
        integer :: y_varid ! ID of y variable
        integer :: element_varid ! ID of element variable
        integer :: adcirc_mesh_varid ! ID of adcirc_mesh variable
        integer :: _varid ! ID of  variable
        integer :: _varid ! ID of  variable
        integer :: _varid ! ID of  variable
        integer :: _varid ! ID of  variable
        integer :: _varid ! ID of  variable
        integer :: _varid ! ID of  variable
        integer :: _varid ! ID of  variable
        integer :: _varid ! ID of  variable
        integer :: _varid ! ID of  variable
        integer :: _varid ! ID of  variable
        integer :: _varid ! ID of  variable
        integer :: _varid ! ID of  variable
        
        ! fort.63.nc: elevation time series
        stat = nf90_def_dim(fort_63_ncid, "time", nf90_unlimited, dimid)
        if (stat.neq.nf90_noerr) call handle_err(stat)
        stat = nf90_def_dim(fort_63_ncid, "node", ???, dimid)
        if (stat.neq.nf90_noerr) call handle_err(stat)
        stat = nf90_def_dim(fort_63_ncid, "nele", ???, dimid)
        if (stat.neq.nf90_noerr) call handle_err(stat)
        stat = nf90_def_dim(fort_63_ncid, "nvertex", ???, dimid)
        if (stat.neq.nf90_noerr) call handle_err(stat)
        stat = nf90_def_dim(fort_63_ncid, "nope", ???, dimid)
        if (stat.neq.nf90_noerr) call handle_err(stat)
        stat = nf90_def_dim(fort_63_ncid, "neta", ???, dimid)
        if (stat.neq.nf90_noerr) call handle_err(stat)
        stat = nf90_def_dim(fort_63_ncid, "max_nvdll", ???, dimid)
        if (stat.neq.nf90_noerr) call handle_err(stat)
        stat = nf90_def_dim(fort_63_ncid, "nbou", ???, dimid)
        if (stat.neq.nf90_noerr) call handle_err(stat)
        stat = nf90_def_dim(fort_63_ncid, "nvel", ???, dimid)
        if (stat.neq.nf90_noerr) call handle_err(stat)
        stat = nf90_def_dim(fort_63_ncid, "max_nvell", ???, dimid)
        if (stat.neq.nf90_noerr) call handle_err(stat)
        stat = nf90_def_dim(fort_63_ncid, "mesh", ???, dimid)
        if (stat.neq.nf90_noerr) call handle_err(stat)
        stat = nf90_def_var(fort_63_ncid, "time", nf90_double, &
            (/  /), varid)
        if (stat.neq.nf90_noerr) call handle_err(stat)
        stat = nf90_put_att
        if (stat.neq.nf90_noerr) call handle_err(stat)
        stat = nf90_enddef(fort_63_ncid)
        if (stat.neq.nf90_noerr) call handle_err(stat)
        stat = nf90_close(fort_63_ncid)
        if (stat.neq.nf90_noerr) call handle_err(stat)
    end subroutine create_fort_63

end module create_files_netcdf
