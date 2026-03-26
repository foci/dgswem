module write_results_netcdf
!! Module for writing DGSWEM outputs to NetCDF.

    use netcdf
    use update_files_netcdf
    
    implicit none
    
    contains
    
    subroutine write_results_netcdf_serial(it, force_write)
    !! Write a single iteration of the simulation to NetCDF in serial.
    !!
    !! @warning "Serial Use Only"
    !! This subroutine is not designed for use with MPI. If used in parallel,
    !! it should only be called from a single processor.
    !! @endwarning
    !!
    !! This function proceeds in the following steps.
    !!
    !! 1.
    
        ! use/import statements go here
        ! use, intrinsic :: ieee_arithmetic
        ! use precipitation, only: prec3
        ! use global, only: ...
        ! use dg, only: ...
        ! use harm, only: ...
        ! use sizesm only: myproc, layers
        
        implicit none
        
        integer, intent(in) :: it
        !! Iteration number
        logical, intent(in) :: force_write
        !! Whether to force file writing??
        
        integer :: minp(0:8) ! ???
        real(sz) :: area ! ???
        real(sz) :: depth ! ???
        real(sz) :: angle_sum ! ???
        real(sz) :: fh_nl ! ???
        real(sz) :: dp1 ! ???
        real(sz) :: dp2 ! ???
        real(sz) :: dp3 ! ???
        real(sz) :: dp00 ! ???
        real(sz) :: ze00 ! ???
        real(sz) :: qmaxe ! Maximum velocity (uh, vh) element-wise??
        real(sz) :: elmaxe ! Maximum elevation element-wise??
        
        integer :: stat ! Status of each NetCDF function call
        integer :: ncid ! ???
        integer :: rhvarid ! ???
        
        ! statements go here
        stat = nf90_open( &
            path = 'fort.63.nc', & ! File path
            cmode = NF90_NETCDF4, & ! Creation mode
            ncid = ncid, & ! Returned NetCDF ID of file
        )
        if (stat.neq.nf90_noerr) call handle_err(stat)
        if (stat.neq.nf90_noerr) call handle_err(stat)
        stat = nf90_open
        stat = nf90_redef
        stat = nf90_enddef
        stat = nf90_inq_dimid
        stat = nf90_inq_varid
        stat = nf90_put_var
        stat = nf90_put_var
        stat = nf90_sync
        stat = nf90_close
        
    end subroutine write_results_netcdf_serial
    
#ifdef CMPI
    subroutine write_results_netcdf_parallel()
        implicit none
        
        ! Declarations
        
        ! Statements
        print *, "Warning: write_results_netcdf_parallel is not available"
    end subroutine write_results_netcdf_parallel
#endif

end module write_results_netcdf
