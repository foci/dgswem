submodule (piodg:piodg_nodal_sub) piodg_63_sub
    !! Implementation of fort.63.nc

    implicit none

    contains

    module procedure piodg_63_open
        ! See interface for arguments and documentation

        integer :: piostat ! Status of most recent operation
        integer :: the_piotype ! For defaulting
        integer :: the_omode ! For defaulting

        ! Set path
        ! Default is whatever is already there
        if (present(path)) then
            self%path = path
        else if (.not. allocated(self%path)) then
            self%path = "fort.63.nc"
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
        call self%piodg_nodal%open(piosystem, self%path, the_piotype, the_omode)

        ! Set dimension IDs
        ! N/A

        ! Set variable IDs
        piostat = pio_inq_varid(self%piofiledesc, self%zeta_varname, &
            self%zeta_vardesc)
        call piofile_check_error(piostat)
    end procedure piodg_63_open

    module procedure piodg_63_close
        ! See interface for arguments and documentation

        type(var_desc_t) :: zeta_vardesc_unset ! Blank struct

        ! Call parent function
        call self%piodg_nodal%close()

        ! Flush dimension IDs
        ! N/A

        ! Flush variable IDs
        self%zeta_vardesc = zeta_vardesc_unset
    end procedure piodg_63_close

end submodule piodg_63_sub
