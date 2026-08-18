submodule (piodg:piodg_nodal_sub) piodg_64_sub
    !! Implementation of fort.64.nc

    implicit none

    contains

    module procedure piodg_64_open
        ! See interface for arguments and documentation

        integer :: piostat ! Status of most recent operation
        integer :: the_piotype ! For defaulting
        integer :: the_omode ! For defaulting

        ! Set path
        ! Default is whatever is already there
        if (present(path)) then
            self%path = path
        else if (.not. allocated(self%path)) then
            self%path = "fort.64.nc"
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

        ! Set vector names
        self%vector_u_varname = "u_vel"
        self%vector_v_varname = "v_vel"

        ! Call parent function
        call self%piodg_nodal_vector%open(piosystem, self%path, the_piotype, &
            the_omode)
    end procedure piodg_64_open

end submodule piodg_64_sub
