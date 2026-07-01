submodule (piodg:piodg_nodal_sub) piodg_63_sub
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

    module procedure piodg_63_write_step
        ! See interface for arguments and documentation

        integer :: piostat ! Status of most recent operation
        ! real :: t_arr(1) ! Vector structure for writing t
        ! integer :: time_start(1) ! Starting position for time data
        ! integer :: time_count(1) ! time data block size
        integer :: record_count ! For indexing time step

        ! File is expected in data mode

        ! Get time slice from time dimension length
        piostat = pio_inquire_dimension(self%piofiledesc, self%time_dimid, &
            len=record_count)
        call piofile_check_error(piostat)

        ! time
        ! t_arr(1) = t
        ! time_start(1) = record_count + 1 ! Position to write at
        ! time_count(1) = 1 ! Number of steps to write
        piostat = pio_put_var(self%piofiledesc, self%time_vardesc, &
            [record_count + 1], t)
        call piofile_check_error(piostat)

        ! zeta
        call pio_setframe(self%piofiledesc, self%zeta_vardesc, &
            int(record_count + 1, pio_offset_kind))
        call pio_write_darray(self%piofiledesc, self%zeta_vardesc, &
            piodesc, zeta, piostat)
        call piofile_check_error(piostat)
        ! call pio_advanceframe(self%piofiledesc, self%zeta_vardesc) ! Broken

        ! Call sync to write to file
        call pio_syncfile(self%piofiledesc)
    end procedure piodg_63_write_step

end submodule piodg_63_sub
