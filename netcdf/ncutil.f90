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

        procedure, public :: create => ncfile_create
        procedure, public :: open => ncfile_open
        procedure, public :: close => ncfile_close

    end type ncfile

    contains

    subroutine ncfile_create(self, path, cmode)
        !! Initialization function for
        !! [[ncutil(module):ncfile(type)]] which creates a new file. Left in
        !! define mode.
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
        !! NetCDF creation mode. Default is ior(nf90_noclobber, nf90_netcdf4).

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
            the_cmode = ior(nf90_noclobber, nf90_netcdf4)
        end if

        ! Create file
        ncstat = nf90_create(self%path, the_cmode, self%ncid)
        call ncfile_check_error(ncstat)
    end subroutine ncfile_create

    subroutine ncfile_open(self, path, mode)
        !! Initialization function for
        !! [[ncutil(module):ncfile(type)]] which opens an existing file. Left in
        !! define mode.
        !!
        !! @warning "File Not Closed"
        !! After initialization, the underlying file object is in define mode.
        !! The user must remember to close it.
        !! @endwarning

        implicit none

        class(ncfile), intent(inout) :: self
        !! The wrapper object of the file
        character(len=*), optional, intent(in) :: path
        !! Path and name for the NetCDF file. Default is whatever path is
        !! already set. If the path is not already set and nothing is provided,
        !! the program will stop with an error.
        integer, optional, intent(in) :: mode
        !! NetCDF open mode. Default is nf90_nowrite, i.e. read-only.

        integer :: ncstat ! Status of most recent operation
        integer :: the_mode ! For defaulting

        ! Set path
        ! Default is whatever is already there
        if (present(path)) then
            self%path = path
        else if (.not. allocated(self%path)) then
            error stop "NC util error: Tried to open file with no path."
        end if

        ! Set mode
        if (present(mode)) then
            the_mode = mode
        else
            the_mode = nf90_nowrite
        end if

        ! Open file
        ncstat = nf90_open(self%path, the_mode, self%ncid)
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

        if (ncstat /= nf90_noerr) then
            error stop "NetCDF Error: " // nf90_strerror(ncstat)
        end if
    end subroutine ncfile_check_error

end module ncutil
