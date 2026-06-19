module pioutil
    !! Module for manipulating a generic NetCDF file in PIO.

    use pio

    implicit none

    type :: piofile
        !! A lightweight structure to hold information necessary for
        !! manipulating a NetCDF file with PIO.

        character(len=:), allocatable :: path
        !! Path of the NetCDF file
        type(file_desc_t) :: piofile
        !! Reference to the NetCDF file

        contains

        procedure, public :: create => piofile_create
        procedure, public :: open => piofile_open
        procedure, public :: pclose => piofile_close
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
        type(iosystem_desc_t), intent(inout) :: piosystem
        !! The PIO system object
        character(len=*), optional, intent(in) :: path
        !! Path and name for the NetCDF file. Default is untitled.nc.
        integer, optional, intent(in) :: piotype
        !! PIO file and I/O type. Default is pio_iotype_netcdf4p.
        integer, optional, intent(in) :: cmode
        !! PIO creation mode. Default is pio_write.

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
            the_cmode = pio_write
        end if

        ! Create file
        piostat = pio_createfile(self%piosystem, self%piofile, &
            the_piotype, self%path, the_cmode)
        call piofile_check_error(piostat)
    end subroutine piofile_init

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
        type(iosystem_desc_t), intent(inout) :: piosystem
        !! The PIO system object
        character(len=*), optional, intent(in) :: path
        !! Path and name for the NetCDF file. Default is whatever path is
        !! already set or untitled.nc.
        integer, optional, intent(in) :: piotype
        !! PIO file and I/O type. Default is pio_iotype_netcdf4p.
        integer, optional, intent(in) :: omode
        !! PIO open mode. Default is pio_write.

        integer :: piostat ! Status of most recent operation
        integer :: the_piotype ! For defaulting
        integer :: the_omode ! For defaulting

        ! Set path
        ! Default is either untitled.nc
        ! or whatever is already there
        if (present(path)) then
            self%path = path
        else if (.not. allocated(self%path))
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
            the_cmode = pio_write
        end if

        ! Open file
    end subroutine piofile_open

    subroutine piofile_check_error(piostat)
        !! Look for and handle an error code from the PIO library.

        implicit none

        integer, intent(in) :: piostat
        !! Status output by the PIO library function call

        integer :: ierr ! Ironically, an error code for the error function
        character(len=:), allocatable :: msg ! For retrieving the PIO message

        if (piostat /= pio_noerr) then
            ierr = pio_strerror(piostat, msg)
            error stop "PIO Error: " // msg
        end if
    end subroutine piofile_check_error

end module pioutil
