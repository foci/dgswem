submodule (ncdg:ncdg_nodal_sub) ncdg_nodal_scalar_sub
!! Implementation of fort.63.nc

    implicit none

    contains

    module procedure ncdg_nodal_scalar_create
        ! See interface for arguments and documentation

        character(len=:), allocatable :: the_path ! For defaulting
        integer :: the_cmode ! For defaulting

        ! Set path
        if (present(path)) then
            the_path = path
        else
            the_path = "nodalscalar.nc"
        end if

        ! Set cmode
        if (present(cmode)) then
            the_cmode = cmode
        else
            the_cmode = ior(nf90_noclobber, nf90_netcdf4)
        end if

        ! Call parent function
        call self%ncfile%create(path=the_path, cmode=the_cmode)

        deallocate(the_path)
    end procedure ncdg_nodal_scalar_create

    module procedure ncdg_nodal_scalar_open
        ! See interface for arguments and documentation

        integer :: ncstat ! Status of most recent operation
        integer :: the_mode ! For defaulting

        ! Set path
        ! Default is whatever is already there
        if (present(path)) then
            self%path = path
        else if (.not. allocated(self%path)) then
            self%path = "nodalscalar.nc"
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
        ncstat = nf90_inq_varid(self%ncid, self%scalar_varname, &
            self%scalar_varid)
        call ncfile_check_error(ncstat)
    end procedure ncdg_nodal_scalar_open

    module procedure ncdg_nodal_scalar_close
        ! See interface for arguments and documentation

        ! integer :: ncstat ! Status of most recent operation

        ! Call parent function
        call self%ncdg_nodal%close()

        ! Flush dimension IDs
        ! N/A

        ! Flush variable IDs
        self%scalar_varid = -1
    end procedure ncdg_nodal_scalar_close

    module procedure ncdg_nodal_scalar_write_step
        ! See interface for arguments and documentation

        integer :: ncstat ! Status of most recent operation
        logical :: the_sync ! For defaulting
        integer :: np ! Number of nodes
        real :: t_arr(1) ! Vector structure for writing t
        integer :: time_start(1) ! Starting position for time data
        integer :: time_count(1) ! time data block size
        integer :: scalar_start(2) ! Starting position for zeta data
        integer :: scalar_count(2) ! zeta data block size
        integer :: record_count ! For indexing time step

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

        ! Get time slice from time dimension length
        ncstat = nf90_inquire_dimension(self%ncid, self%time_dimid, &
            len=record_count)
        call ncfile_check_error(ncstat)

        ! time
        t_arr(1) = t
        time_start(1) = record_count + 1 ! Position to write at
        time_count(1) = 1 ! Number of steps to write
        ncstat = nf90_put_var(self%ncid, self%time_varid, &
            t_arr, time_start, time_count)
        call ncfile_check_error(ncstat)

        ! zeta
        scalar_start(1) = 1 ! node dimension start
        scalar_start(2) = record_count + 1 ! time dimension start
        scalar_count(1) = np ! node dimension span
        scalar_count(2) = 1 ! time dimension span
        ncstat = nf90_put_var(self%ncid, self%scalar_varid, &
            scalar, scalar_start, scalar_count)
        call ncfile_check_error(ncstat)

        ! After each write, sync
        if (the_sync) then
            ncstat = nf90_sync(self%ncid)
            call ncfile_check_error(ncstat)
        end if
    end procedure ncdg_nodal_scalar_write_step

end submodule ncdg_nodal_scalar_sub
