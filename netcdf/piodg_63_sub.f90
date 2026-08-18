submodule (piodg:piodg_nodal_scalar_sub) piodg_63_sub
    !! Implementation of fort.63.nc

    use netcdf, only : nf90_enotindefine

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

        ! Set scalar name
        self%scalar_varname = "zeta"

        ! Call parent function
        call self%piodg_nodal_scalar%open(piosystem, self%path, the_piotype, &
            the_omode)
    end procedure piodg_63_open

end submodule piodg_63_sub
