module write_results_netcdf
    implicit none
    
    ! Type definitions
    ! N/A
    
    ! Constants
    ! N/A
    
    ! Procedures
    contains
    
    subroutine write_results_netcdf_serial(it, force_write)
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
        
        ! statements go here
        
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
