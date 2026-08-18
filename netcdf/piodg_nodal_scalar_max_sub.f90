submodule (piodg:piodg_nodal_sub) piodg_nodal_scalar_max_sub
    !! Implementation of maxele.63.nc

    implicit none

    contains

    module procedure piodg_nodal_scalar_max_open
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
        piostat = pio_inq_varid(self%piofiledesc, self%scalar_max_varname, &
            self%scalar_max_vardesc)
        call piofile_check_error(piostat)
        piostat = pio_inq_varid(self%piofiledesc, self%time_of_scalar_max_varname, &
            self%time_of_scalar_max_vardesc)
        call piofile_check_error(piostat)
    end procedure piodg_nodal_scalar_max_open

    module procedure piodg_nodal_scalar_max_close
        ! See interface for arguments and documentation

        type(var_desc_t) :: scalar_max_vardesc_unset ! Blank struct
        type(var_desc_t) :: time_of_scalar_max_vardesc_unset ! Blank struct

        ! Call parent function
        call self%piodg_nodal%close()

        ! Flush dimension IDs
        ! N/A

        ! Flush variable IDs
        self%scalar_max_vardesc = scalar_max_vardesc_unset
        self%time_of_scalar_max_vardesc = time_of_scalar_max_vardesc_unset
    end procedure piodg_nodal_scalar_max_close

    module procedure piodg_nodal_scalar_max_write_step
        ! See interface for arguments and documentation

        integer :: piostat ! Status of most recent operation
        logical :: the_compute_max ! For defaulting
        logical :: the_sync ! For defaulting
        real, allocatable :: scalar_max(:)
        real, allocatable :: scalar_max_old(:)
        real, allocatable :: time_of_scalar_max(:)

        ! Set compute_max
        if (present(compute_max)) then
            the_compute_max = compute_max
        else
            the_compute_max = .false.
        end if

        ! Set sync
        if (present(sync)) then
            the_sync = sync
        else
            the_sync = .false.
        end if

        ! File is expected in data mode

        ! Get existing data for comparison
        allocate(scalar_max_old(np))
        allocate(time_of_scalar_max(np))
        call pio_read_darray(self%piofiledesc, self%scalar_max_vardesc, &
            piodesc, scalar_max_old, piostat)
        call piofile_check_error(piostat)
        call pio_read_darray(self%piofiledesc, self%time_of_scalar_max_vardesc, &
            piodesc, time_of_scalar_max, piostat)
        call piofile_check_error(piostat)

        ! scalar_max
        if (the_compute_max) then
            scalar_max = merge( &
                scalar, & ! Used if true: new scalar
                scalar_max_old, & ! Used if false: old scalar max
                scalar > scalar_max_old & ! Condition: new scalar is strictly greater
            )
        else
            scalar_max = scalar
        end if
        call pio_write_darray(self%piofiledesc, self%scalar_max_vardesc, &
            piodesc, scalar_max, piostat)
        call piofile_check_error(piostat)

        ! time_of_scalar_max
        time_of_scalar_max = merge( &
            spread(t, 1, np), & ! Used if true: current time
            time_of_scalar_max, & ! Used if false: old time
            scalar_max > scalar_max_old & ! Condition: new scalar is strictly greater
        )
        call pio_write_darray(self%piofiledesc, self%time_of_scalar_max_vardesc, &
            piodesc, time_of_scalar_max, piostat)
        call piofile_check_error(piostat)

        ! Deallocate arrays
        deallocate(scalar_max)
        deallocate(scalar_max_old)
        deallocate(time_of_scalar_max)

        ! Call sync to write to file
        if (the_sync) call pio_syncfile(self%piofiledesc)
    end procedure piodg_nodal_scalar_max_write_step

end submodule piodg_nodal_scalar_max_sub
