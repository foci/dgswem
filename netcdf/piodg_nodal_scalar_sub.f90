submodule (piodg:piodg_nodal_sub) piodg_nodal_scalar_sub
    !! Implementation of generic nodal scalar

    use netcdf, only : nf90_enotindefine

    implicit none

    contains

    module procedure piodg_nodal_scalar_open
        ! See interface for arguments and documentation

        integer :: piostat ! Status of most recent operation
        integer :: the_piotype ! For defaulting
        integer :: the_omode ! For defaulting

        ! Set path
        ! Default is whatever is already there
        if (present(path)) then
            self%path = path
        else if (.not. allocated(self%path)) then
            self%path = "nodalscalar.nc"
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
        piostat = pio_inq_varid(self%piofiledesc, self%scalar_varname, &
            self%scalar_vardesc)
        call piofile_check_error(piostat)
    end procedure piodg_nodal_scalar_open

    module procedure piodg_nodal_scalar_close
        ! See interface for arguments and documentation

        type(var_desc_t) :: scalar_vardesc_unset ! Blank struct

        ! Call parent function
        call self%piodg_nodal%close()

        ! Flush dimension IDs
        ! N/A

        ! Flush variable IDs
        self%scalar_vardesc = scalar_vardesc_unset
    end procedure piodg_nodal_scalar_close

    module procedure piodg_nodal_scalar_write_step
        ! See interface for arguments and documentation

        integer :: piostat ! Status of most recent operation
        logical :: the_sync ! For defaulting
        integer :: record_count ! For indexing time step

        ! Set sync
        if (present(sync)) then
            the_sync = sync
        else
            the_sync = .false.
        end if

        ! File is expected in data mode

        ! Get time slice from time dimension length
        piostat = pio_inquire_dimension(self%piofiledesc, self%time_dimid, &
            len=record_count)
        call piofile_check_error(piostat)

        ! time
        piostat = pio_put_var(self%piofiledesc, self%time_vardesc, &
            [record_count + 1], t)
        call piofile_check_error(piostat)

        ! scalar
        call pio_setframe(self%piofiledesc, self%scalar_vardesc, &
            int(record_count + 1, pio_offset_kind))
        call pio_write_darray(self%piofiledesc, self%scalar_vardesc, &
            piodesc, scalar, piostat)
        call piofile_check_error(piostat)

        ! Call sync to write to file
        if (the_sync) call pio_syncfile(self%piofiledesc)
    end procedure piodg_nodal_scalar_write_step

end submodule piodg_nodal_scalar_sub
