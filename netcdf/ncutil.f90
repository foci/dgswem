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
    
    subroutine ncfile_init(self, path, cmode, comm, info)
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
        !! NetCDF creation mode. Default is ior(nf90_clobber, nf90_netcdf4).
        integer, optional, intent(in) :: comm
        !! MPI communicator variable for parallel I/O
        integer, optional, intent(in) :: info
        !! MPI info variable for parallel I/O
        
        integer :: ncstat ! Status of most recent operation
        integer :: the_cmode ! For defaulting

        ! Set path
        if (present(path)) then
            self%path = path
        else
            self%path = "untitled.nc"
        end if

        ! Set cmode
        if (present(cmode)) then
            the_cmode = cmode
        else
            the_cmode = ior(nf90_clobber, nf90_netcdf4)
        end if
        
        ! Create file
        if (present(comm) .and. present(info)) then
            ncstat = nf90_create(self%path, the_cmode, self%ncid, &
                comm=comm, info=info)
        else
            ncstat = nf90_create(self%path, the_cmode, self%ncid)
        end if
        call ncfile_check_error(ncstat)
    end subroutine ncfile_init
    
    subroutine ncfile_open(self, mode, comm, info)
        !! Open the NetCDF file for manipulation.
        
        implicit none
        
        class(ncfile), intent(inout) :: self
        !! The wrapper object of the file
        integer, optional, intent(in) :: mode
        !! NetCDF open mode. Default is nf90_nowrite, i.e. read-only.
        integer, optional, intent(in) :: comm
        !! MPI communicator variable for parallel I/O
        integer, optional, intent(in) :: info
        !! MPI info variable for parallel I/O
        
        integer :: ncstat ! Status of most recent operation
        integer :: the_mode ! For defaulting

        ! Set mode
        if (present(mode)) then
            the_mode = mode
        else
            the_mode = nf90_nowrite
        end if

        ! Open file
        if (present(comm) .and. present(info)) then
            ncstat = nf90_open(self%path, mode, self%ncid, &
                comm=comm, info=info)
        else
            ncstat = nf90_open(self%path, mode, self%ncid)
        end if
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
