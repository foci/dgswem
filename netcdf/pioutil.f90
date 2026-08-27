module pioutil
    !! Module for manipulating a generic NetCDF file in PIO.

    use pio

    implicit none

    type :: piofile
        !! A lightweight structure to hold information necessary for
        !! manipulating a NetCDF file with PIO.

        character(len=:), allocatable :: path
        !! Path of the NetCDF file
        type(file_desc_t) :: piofiledesc
        !! Reference to the NetCDF file

        contains

        procedure, public :: create => piofile_create
        procedure, public :: open => piofile_open
        procedure, public :: close => piofile_close
    end type piofile

    contains

    subroutine piofile_create(self, piosystem, path, piotype, cmode)
        !! Parallel-access initialization function for
        !! [[pioutil(module):piofile(type)]] which creates a new file. Left in
        !! define mode.
        !!
        !! @warning "File Not Closed"
        !! After initialization, the underlying file object is in define mode.
        !! The user must remember to close it.
        !! @endwarning

        implicit none

        class(piofile), intent(inout) :: self
        !! The wrapper object being initialized
        type(iosystem_desc_t), intent(inout), target :: piosystem
        !! The PIO system object
        character(len=*), optional, intent(in) :: path
        !! Path and name for the NetCDF file. Default is untitled.nc.
        integer, optional, intent(in) :: piotype
        !! PIO file and I/O type. Default is pio_iotype_netcdf4p.
        integer, optional, intent(in) :: cmode
        !! PIO creation mode. Default is pio_noclobber.

        integer :: piostat ! Status of most recent operation
        integer :: the_piotype ! For defaulting
        integer :: the_cmode ! For defaulting

        ! Set path
        if (present(path)) then
            self%path = path
        else
            self%path = "untitled.nc"
        end if

        ! Set piotype
        if (present(piotype)) then
            the_piotype = piotype
        else
            the_piotype = pio_iotype_netcdf4p
        end if

        ! Set cmode
        if (present(cmode)) then
            the_cmode = cmode
        else
            the_cmode = pio_noclobber
        end if

        ! Create file
        piostat = pio_createfile(piosystem, self%piofiledesc, &
            the_piotype, self%path, the_cmode)
        call piofile_check_error(piostat)
    end subroutine piofile_create

    subroutine piofile_open(self, piosystem, path, piotype, omode)
        !! Parallel-access initialization function for
        !! [[pioutil(module):piofile(type)]] which opens an existing file. Left
        !! in define mode.
        !!
        !! @warning "File Not Closed"
        !! After initialization, the underlying file object is in define mode.
        !! The user must remember to close it.
        !! @endwarning

        implicit none

        class(piofile), intent(inout) :: self
        !! The wrapper object of the file
        type(iosystem_desc_t), intent(inout), target :: piosystem
        !! The PIO system object
        character(len=*), optional, intent(in) :: path
        !! Path and name for the NetCDF file. Default is whatever path is
        !! already set. If the path is not already set and nothing is provided,
        !! the program will stop with an error.
        integer, optional, intent(in) :: piotype
        !! PIO file and I/O type. Default is pio_iotype_netcdf4p.
        integer, optional, intent(in) :: omode
        !! PIO open mode. Default is pio_nowrite.

        integer :: piostat ! Status of most recent operation
        integer :: the_piotype ! For defaulting
        integer :: the_omode ! For defaulting

        ! Set path
        ! Default is whatever is already there
        if (present(path)) then
            self%path = path
        else if (.not. allocated(self%path)) then
            error stop "PIO util error: Tried to open file with no path."
        end if

        ! Set piotype
        if (present(piotype)) then
            the_piotype = piotype
        else
            the_piotype = pio_iotype_netcdf4p
        end if

        ! Set cmode
        if (present(omode)) then
            the_omode = omode
        else
            the_omode = pio_nowrite
        end if

        ! Open file
        piostat = pio_openfile(piosystem, self%piofiledesc, &
            the_piotype, self%path, the_omode)
        call piofile_check_error(piostat)
    end subroutine piofile_open

    subroutine piofile_close(self)
        !! Close the NetCDF file.

        implicit none

        class(piofile), intent(inout) :: self
        !! The wrapper object of the file

        ! integer :: piostat ! Status of most recent operation
        type(file_desc_t) :: piofiledesc_unset ! Blank struct

        call pio_closefile(self%piofiledesc)

        ! Flush file struct
        self%piofiledesc = piofiledesc_unset
    end subroutine piofile_close

    subroutine piofile_check_error(piostat)
        !! Look for and handle an error code from the PIO library.

        implicit none

        integer, intent(in) :: piostat
        !! Status output by the PIO library function call

        integer :: ierr ! Ironically, an error code for the error function
        character(len=256) :: msg ! For retrieving the PIO message

        if (piostat /= pio_noerr) then
            ierr = pio_strerror(piostat, msg)
            error stop "PIO Error: " // trim(msg)
        end if
    end subroutine piofile_check_error

end module pioutil
