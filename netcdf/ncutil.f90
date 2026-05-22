module ncutil
!! Module for manipulating a generic NetCDF file.

    use netcdf
    
    implicit none
    
    type ncfile
        !! A lightweight structure to hold information necesary for manipulating
        !! a NetCDF file.
        
        character(len=:), allocatable :: path
        !! Path of the NetCDF file
        integer :: ncid = -1
        !! ID of the NetCDF file
        !!
        !! @warning "Transient ID"
        !! NetCDF IDs are not guaranteed to be consistent across open/close
        !! operations.
        !! @endwarning
        
        contains
        
        procedure, public :: init => ncfile_init
        procedure, public :: open => ncfile_open
        procedure, public :: close => ncfile_close
        
    end type ncfile
    
    contains
    
    subroutine ncfile_init(self, path, cmode)
        !! Initialization function for
        !! [[netcdf_file(module):netcdf_file(type)]]. Left in define mode.
        !!
        !! @warning "File Not Closed"
        !! After initialization, the underlying file object is in define mode.
        !! The user must remember to close it.
        !! @endwarning
        
        implicit none
        
        class(ncfile), intent(inout) :: self
        !! The wrapper object being initialized
        character(len=*), optional, intent(in) :: path
        !! Path and name for the NetCDF file. Default is untitled.nc.
        integer, optional, intent(in) :: cmode
        !! NetCDF creation mode. Default is nf90_clobber.
        
        integer :: ncstat ! Status of most recent operation

        ! Set path
        if (present(path)) then
            self%path = path
        else
            self%path = "untitled.nc"
        endif
        
        ! Create file
        if (present(cmode)) then
            ncstat = nf90_create(self%path, cmode, self%ncid)
        else
            ncstat = nf90_create(self%path, nf90_clobber, self%ncid)
        endif
        call ncfile_check_error(ncstat)
    end subroutine ncfile_init
    
    subroutine ncfile_open(self, mode)
        !! Open the NetCDF file for manipulation.
        
        implicit none
        
        class(ncfile), intent(inout) :: self
        !! The wrapper object of the file
        integer, optional, intent(in) :: mode
        !! NetCDF open mode. Default is nf90_nowrite, i.e. read-only.
        
        integer :: ncstat ! Status of most recent operation
        
        if (present(mode)) then
            ncstat = nf90_open(self%path, mode, self%ncid)
        else
            ncstat = nf90_open(self%path, nf90_nowrite, self%ncid)
        endif
        call ncfile_check_error(ncstat)
    end subroutine ncfile_open
    
    subroutine ncfile_close(self)
        !! Close the NetCDF file.
        
        implicit none
        
        class(ncfile), intent(inout) :: self
        !! The wrapper object of the file
        
        integer :: ncstat ! Status of most recent operation
        
        ncstat = nf90_close(self%ncid)
        call ncfile_check_error(ncstat)
        
        self%ncid = -1
    end subroutine ncfile_close
    
    subroutine ncfile_check_error(ncstat)
        !! Look for and handle an error code from the NetCDF library.
        
        implicit none
        
        integer, intent(in) :: ncstat
        !! Status output by NetCDF library function call
        
        if (ncstat /= nf90_noerr) call ncfile_handle_error(ncstat)
    end subroutine ncfile_check_error
    
    subroutine ncfile_handle_error(ncstat)
        !! Handle an error code from the NetCDF library.
        
        implicit none
        
        integer, intent(in) :: ncstat
        !! Status of most recent operation
        
        print *, "NetCDF Error: ", nf90_strerror(ncstat)
        stop
    end subroutine ncfile_handle_error

end module ncutil
