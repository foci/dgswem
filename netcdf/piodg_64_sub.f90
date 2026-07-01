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

        ! Call parent function
        call self%piodg_nodal%open(piosystem, self%path, the_piotype, the_omode)

        ! Set dimension IDs
        ! N/A

        ! Set variable IDs
        piostat = pio_inq_varid(self%piofiledesc, self%u_vel_varname, &
            self%u_vel_vardesc)
        call piofile_check_error(piostat)
        piostat = pio_inq_varid(self%piofiledesc, self%v_vel_varname, &
            self%v_vel_vardesc)
        call piofile_check_error(piostat)
    end procedure piodg_64_open

    module procedure piodg_64_close
        ! See interface for arguments and documentation

        type(var_desc_t) :: u_vel_vardesc_unset ! Blank struct
        type(var_desc_t) :: v_vel_vardesc_unset ! Blank struct

        ! Call parent function
        call self%piodg_nodal%close()

        ! Flush dimension IDs
        ! N/A

        ! Flush variable IDs
        self%u_vel_vardesc = u_vel_vardesc_unset
        self%v_vel_vardesc = v_vel_vardesc_unset
    end procedure piodg_64_close

    module procedure piodg_64_write_step
        ! See interface for arguments and documentation

        integer :: piostat ! Status of most recent operation
        integer :: record_count ! For indexing time step

        ! File is expected in data mode

        ! Get time slice from time dimension length
        piostat = pio_inquire_dimension(self%piofiledesc, self%time_dimid, &
            len=record_count)
        call piofile_check_error(piostat)

        ! time
        piostat = pio_put_var(self%piofiledesc, self%time_vardesc, &
            [record_count + 1], t)
        call piofile_check_error(piostat)

        ! u_vel
        call pio_setframe(self%piofiledesc, self%u_vel_vardesc, &
            int(record_count + 1, pio_offset_kind))
        call pio_write_darray(self%piofiledesc, self%u_vel_vardesc, &
            piodesc, u_vel, piostat)
        call piofile_check_error(piostat)

        ! v_vel
        call pio_setframe(self%piofiledesc, self%v_vel_vardesc, &
            int(record_count + 1, pio_offset_kind))
        call pio_write_darray(self%piofiledesc, self%v_vel_vardesc, &
            piodesc, v_vel, piostat)
        call piofile_check_error(piostat)

        ! Call sync to write to file
        call pio_syncfile(self%piofiledesc)
    end procedure piodg_64_write_step

end submodule piodg_64_sub
