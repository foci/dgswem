submodule (piodg:piodg_nodal_sub) piodg_maxele_sub
    !! Implementation of maxele.63.nc

    implicit none

    contains

    module procedure piodg_maxele_open
        ! See interface for arguments and documentation

        integer :: piostat ! Status of most recent operation
        integer :: the_piotype ! For defaulting
        integer :: the_omode ! For defaulting

        ! Set path
        ! Default is whatever is already there
        if (present(path)) then
            self%path = path
        else if (.not. allocated(self%path)) then
            self%path = "maxele.63.nc"
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
        piostat = pio_inq_varid(self%piofiledesc, self%zeta_max_varname, &
            self%zeta_max_vardesc)
        call piofile_check_error(piostat)
        piostat = pio_inq_varid(self%piofiledesc, self%time_of_zeta_max_varname, &
            self%time_of_zeta_max_vardesc)
        call piofile_check_error(piostat)
    end procedure piodg_maxele_open

    module procedure piodg_maxele_close
        ! See interface for arguments and documentation

        type(var_desc_t) :: zeta_max_vardesc_unset ! Blank struct
        type(var_desc_t) :: time_of_zeta_max_vardesc_unset ! Blank struct

        ! Call parent function
        call self%piodg_nodal%close()

        ! Flush dimension IDs
        ! N/A

        ! Flush variable IDs
        self%zeta_max_vardesc = zeta_max_vardesc_unset
        self%time_of_zeta_max_vardesc = time_of_zeta_max_vardesc_unset
    end procedure piodg_maxele_close

end submodule piodg_maxele_sub
