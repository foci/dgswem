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

    module procedure piodg_maxele_write_step
        ! See interface for arguments and documentation

        integer :: piostat ! Status of most recent operation
        logical :: the_sync ! For defaulting
        real, allocatable :: zeta_max_old(:)
        real, allocatable :: time_of_zeta_max(:)

        ! Set sync
        if (present(sync)) then
            the_sync = sync
        else
            the_sync = .false.
        end if

        ! File is expected in data mode

        ! Get existing data for comparison
        allocate(zeta_max_old(np))
        allocate(time_of_zeta_max(np))
        call pio_read_darray(self%piofiledesc, self%zeta_max_vardesc, &
            piodesc, zeta_max_old, piostat)
        call piofile_check_error(piostat)
        call pio_read_darray(self%piofiledesc, self%time_of_zeta_max_vardesc, &
            piodesc, time_of_zeta_max, piostat)
        call piofile_check_error(piostat)

        ! zeta_max
        call pio_write_darray(self%piofiledesc, self%zeta_max_vardesc, &
            piodesc, zeta_max, piostat)
        call piofile_check_error(piostat)

        ! time_of_zeta_max
        time_of_zeta_max = merge( &
            spread(t, 1, np), & ! Used if true: current time
            time_of_zeta_max, & ! Used if false: old time
            zeta_max > zeta_max_old & ! Condition: new maximum is strictly greater
        )
        call pio_write_darray(self%piofiledesc, self%time_of_zeta_max_vardesc, &
            piodesc, time_of_zeta_max, piostat)
        call piofile_check_error(piostat)
        deallocate(zeta_max_old)
        deallocate(time_of_zeta_max)

        ! Call sync to write to file
        if (the_sync) call pio_syncfile(self%piofiledesc)
    end procedure piodg_maxele_write_step

end submodule piodg_maxele_sub
