submodule (piodg) piodg_file_sub
!! Implementation of functions common to all DGSWEM output files.

    implicit none

    contains

    module procedure piodg_file_open
        ! See interface for arguments and documentation

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

        ! Set omode
        if (present(omode)) then
            the_omode = omode
        else
            the_omode = pio_nowrite
        end if

        ! Call parent function
        call self%piofile%open(piosystem, self%path, the_piotype, the_omode)

        ! Set dimension IDs
        piostat = pio_inq_dimid(self%piofiledesc, self%time_dimname, &
            self%time_dimid)
        call piofile_check_error(piostat)

        ! Set variable IDs
        piostat = pio_inq_varid(self%piofiledesc, self%time_varname, &
            self%time_vardesc)
        call piofile_check_error(piostat)
    end procedure piodg_file_open

    module procedure piodg_file_close
        ! See interface for arguments and documentation

        type(var_desc_t) :: time_vardesc_unset ! Blank struct

        ! Call parent function
        call self%piofile%close()

        ! Flush dimension IDs
        self%time_dimid = -1

        ! Flush variable IDs
        self%time_vardesc = time_vardesc_unset
    end procedure piodg_file_close

end submodule piodg_file_sub
