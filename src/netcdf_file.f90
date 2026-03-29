module netcdf_file
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
        !! Create a new NetCDF file
        procedure, public :: open => ncfile_open
        !! Open the file
        procedure, public :: close => ncfile_close
        !! Close the file
        
    end type ncfile
    
    contains
    
    subroutine ncfile_init(self, path, cmode)
        !! Initialization function for
        !! [[netcdf_file(module):netcdf_file(type)]].
        !!
        !! @warning "File Not Closed"
        !! After initialization, the underlying file object is in define mode.
        !! The user must remember to close it.
        !! @endwarning
        
        implicit none
        
        class(ncfile), intent(inout) :: self
        !! The wrapper object being initialized
        character(len=*), intent(in) :: path
        !! Path and name for the NetCDF file
        integer, optional, intent(in) :: cmode
        !! NetCDF creation mode. Default is nf90_clobber.
        
        integer :: stat ! Status of most recent operation
        
        self%path = path
        
        ! Create file
        if (present(cmode)) then
            stat = nf90_create(self%path, cmode, self%ncid)
        else
            stat = nf90_create(self%path, nf90_clobber, self%ncid)
        endif
        if (stat.neq.nf90_noerr) call ncfile_handle_error(stat)
    end subroutine ncfile_init
    
    subroutine ncfile_open(self, mode)
        !! Open the NetCDF file for manipulation.
        
        implicit none
        
        class(ncfile), intent(inout) :: self
        !! The wrapper object of the file
        integer, optional, intent(in) :: mode
        !! NetCDF open mode. Default is nf90_nowrite, i.e. read-only.
        
        integer :: stat ! Status of most recent operation
        
        if (present(mode)) then
            stat = nf90_open(self%path, mode, self%ncid)
        else
            stat = nf90_open(self%path, nf90_nowrite, self%ncid)
        endif
        if (stat.neq.nf90_noerr) call ncfile_handle_error(stat)
    end subroutine ncfile_open
    
    subroutine ncfile_close(self)
        !! Close the NetCDF file.
        
        implicit none
        
        class(ncfile), intent(inout) :: self
        !! The wrapper object of the file
        
        integer :: stat ! Status of most recent operation
        
        stat = nf90_close(self%ncid)
        if (stat.neq.nf90_noerr) call ncfile_handle_error(stat)
        
        self%ncid = -1
    end subroutine ncfile_close
    
    subroutine ncfile_handle_error(stat)
        !! Handle an error code from the NetCDF library.
        
        implicit none
        
        integer, intent(in) :: stat
        !! Status of most recent operation
        
        print *, "NetCDF Error: ", nf90_strerror(stat)
        stop
    end subroutine ncfile_handle_error

end module netcdf_file