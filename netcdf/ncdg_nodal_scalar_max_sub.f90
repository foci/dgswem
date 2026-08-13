submodule (ncdg:ncdg_nodal_scalar_sub) ncdg_nodal_scalar_max_sub
!! Implementation of generic nodal scalar max

    implicit none

    contains

    module procedure ncdg_nodal_scalar_max_open
        ! See interface for arguments and documentation

        integer :: ncstat ! Status of most recent operation
        integer :: the_mode ! For defaulting

        ! Set path
        ! Default is whatever is already there
        if (present(path)) then
            self%path = path
        else if (.not. allocated(self%path)) then
            self%path = "maxele.63.nc"
        end if

        ! Set mode
        if (present(mode)) then
            the_mode = mode
        else
            the_mode = nf90_nowrite
        end if

        ! Call parent function
        call self%ncdg_nodal%open(path=self%path, mode=the_mode)

        ! Set dimension IDs
        ! N/A

        ! Set variable IDs
        ncstat = nf90_inq_varid(self%ncid, self%scalar_max_varname, &
            self%scalar_max_varid)
        ncstat = nf90_inq_varid(self%ncid, self%time_of_scalar_max_varname, &
            self%time_of_scalar_max_varid)
        call ncfile_check_error(ncstat)
    end procedure ncdg_nodal_scalar_max_open

    module procedure ncdg_nodal_scalar_max_close
        ! See interface for arguments and documentation

        ! integer :: ncstat ! Status of most recent operation

        ! Call parent function
        call self%ncdg_nodal%close()

        ! Flush dimension IDs
        ! N/A

        ! Flush variable IDs
        self%scalar_max_varid = -1
        self%time_of_scalar_max_varid = -1
    end procedure ncdg_nodal_scalar_max_close

    module procedure ncdg_nodal_scalar_max_write_step
         ! See interface for arguments and documentation

        integer :: ncstat ! Status of most recent operation
        logical :: the_compute_max ! For defaulting
        logical :: the_sync ! For defaulting
        integer :: np ! Number of nodes
        real, allocatable :: scalar_max(:) ! Current maximums,
                                               ! either assumed or computed
        real, allocatable :: scalar_max_old(:) ! Previous maximums,
                                               ! if present
        real, allocatable :: time_of_scalar_max(:) ! Nodewise times when
                                                   ! scalar_max is achieved

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

        ! Exit define mode
        ncstat = nf90_enddef(self%ncid)
        if (ncstat /= nf90_enotindefine) then
            call ncfile_check_error(ncstat)
        endif

        ! Get number of nodes from node dimension length
        ncstat = nf90_inquire_dimension(self%ncid, self%node_dimid, &
            len=np)
        call ncfile_check_error(ncstat)

        ! Get existing data for comparison
        allocate(scalar_max_old(np))
        allocate(time_of_scalar_max(np))
        ncstat = nf90_get_var(self%ncid, self%scalar_max_varid, &
            scalar_max_old)
        call ncfile_check_error(ncstat)
        ncstat = nf90_get_var(self%ncid, self%time_of_scalar_max_varid, &
            time_of_scalar_max)
        call ncfile_check_error(ncstat)

        ! scalar_max
        if (the_compute_max) then
            scalar_max = merge( &
                scalar, &! Used if true: new scalar
                scalar_max_old, & ! Used if false: old scalar max
                scalar > scalar_max_old & ! Condition: new scalar is strictly greater
            )
        else
            scalar_max = scalar
        end if
        ncstat = nf90_put_var(self%ncid, self%scalar_max_varid, &
            scalar_max)
        call ncfile_check_error(ncstat)

        ! time_of_scalar_max
        time_of_scalar_max = merge( &
            spread(t, 1, np), & ! Used if true: current time
            time_of_scalar_max, & ! Used if false: old time
            scalar_max > scalar_max_old & ! Condition: new maximum is strictly greater
        )
        ncstat = nf90_put_var(self%ncid, self%time_of_scalar_max_varid, &
            time_of_scalar_max)
        call ncfile_check_error(ncstat)
        deallocate(scalar_max_old)
        deallocate(time_of_scalar_max)

        ! After each write, sync
        if (the_sync) then
            ncstat = nf90_sync(self%ncid)
            call ncfile_check_error(ncstat)
        end if
    end procedure ncdg_nodal_scalar_max_write_step

end submodule ncdg_nodal_scalar_max_sub
