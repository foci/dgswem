submodule (piodg:piodg_nodal_scalar_max_sub) piodg_maxele_sub
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

        ! Set scalar and time_of_scalar_max names
        self%scalar_max_varname = "zeta_max"
        self%time_of_scalar_max_varname = "time_of_zeta_max"

        ! Call parent function
        call self%piodg_nodal_scalar_max%open(piosystem, self%path, &
            the_piotype, the_omode)
    end procedure piodg_maxele_open

end submodule piodg_maxele_sub
