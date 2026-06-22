submodule (ncdg:ncdg_nodal_sub) ncdg_64_sub
!! Implementation of fort.64.nc

    implicit none

    contains

    module procedure ncdg_64_create
        ! See interface for arguments and documentation

        character(len=:), allocatable :: the_path ! For defaulting
        integer :: the_cmode ! For defaulting

        ! Set path
        if (present(path)) then
            the_path = path
        else
            the_path = "fort.64.nc"
        end if

        ! Set cmode
        if (present(cmode)) then
            the_cmode = cmode
        else
            the_cmode = ior(nf90_clobber, nf90_netcdf4)
        end if

        ! Call parent function
        call self%ncfile%create(path=the_path, cmode=the_cmode)

        deallocate(the_path)
    end procedure ncdg_64_create

    module procedure ncdg_64_open
        ! See interface for arguments and documentation

        integer :: ncstat ! Status of most recent operation
        integer :: the_mode ! For defaulting

        ! Set path
        ! Default is whatever is already there
        if (present(path)) then
            self%path = path
        else if (.not. allocated(self%path)) then
            self%path = "fort.64.nc"
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
        ncstat = nf90_inq_varid(self%ncid, self%u_vel_varname, &
            self%u_vel_varid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_inq_varid(self%ncid, self%v_vel_varname, &
            self%v_vel_varid)
        call ncfile_check_error(ncstat)
    end procedure ncdg_64_open

    module procedure ncdg_64_close
        ! See interface for arguments and documentation

        ! integer :: ncstat ! Status of most recent operation

        ! Call parent function
        call self%ncdg_nodal%close()

        ! Flush dimension IDs
        ! N/A

        ! Flush variable IDs
        self%u_vel_varid = -1
        self%v_vel_varid = -1
    end procedure ncdg_64_close

    module procedure ncdg_64_set_metadata
        ! See interface for arguments and documentation

        integer :: ncstat ! Status of most recent operation

        ! Put in define mode
        ncstat = nf90_redef(self%ncid)
        if (ncstat /= nf90_eindefine) then
            call ncfile_check_error(ncstat)
        endif

        ! Call parent function
        call self%ncdg_nodal%ncdg_nodal_set_metadata( &
            nt=nt, &
            np=np, &
            ne=ne, &
            nhy=nhy, &
            nope=nope, &
            neta=neta, &
            max_nvdll=max_nvdll, &
            nbou=nbou, &
            nvel=nvel, &
            max_nvell=max_nvell, &
            ics=ics &
        )

        ! Define dimensions
        ! N/A

        ! Define variables and their attributes

        ! u_vel
        if (ics == 1) then ! Cartesian
            ncstat = nf90_def_var(self%ncid, self%u_vel_varname, &
                nf90_double, (/ self%node_dimid, self%time_dimid /), &
                varid=self%u_vel_varid)
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%u_vel_varid, &
                "long_name", "x component of depth-averaged velocity")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%u_vel_varid, &
                "standard_name", "sea_water_x_velocity")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%u_vel_varid, &
                "positive", "right")
            call ncfile_check_error(ncstat)
        else if (ics == 2) then ! Spherical
            ncstat = nf90_def_var(self%ncid, self%u_vel_varname, &
                nf90_double, (/ self%node_dimid, self%time_dimid /), &
                varid=self%u_vel_varid)
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%u_vel_varid, &
                "long_name", "eastward component of depth-averaged " &
                // "velocity")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%u_vel_varid, &
                "standard_name", "eastward_water_velocity")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%u_vel_varid, &
                "positive", "east")
            call ncfile_check_error(ncstat)
        end if
        if (ics == 1 .or. ics == 2) then ! For all cases
            ncstat = nf90_put_att(self%ncid, self%u_vel_varid, &
                "coordinates", "time y x")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%u_vel_varid, &
                "location", "node")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%u_vel_varid, &
                "mesh", "adcirc_mesh")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%u_vel_varid, &
                "units", "m s-1")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%u_vel_varid, &
                "_FillValue", real_fill_value)
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%u_vel_varid, &
                "dry_Value", real_fill_value)
            call ncfile_check_error(ncstat)
        endif

        ! v_vel
        if (ics == 1) then ! Cartesian
            ncstat = nf90_def_var(self%ncid, self%v_vel_varname, &
                nf90_double, (/ self%node_dimid, self%time_dimid /), &
                varid=self%v_vel_varid)
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%v_vel_varid, &
                "long_name", "y component of depth-averaged velocity")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%v_vel_varid, &
                "standard_name", "sea_water_y_velocity")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%v_vel_varid, &
                "positive", "90 degrees counterclockwise from x " &
                // "water velocity")
            call ncfile_check_error(ncstat)
        else if (ics == 2) then ! Spherical
            ncstat = nf90_def_var(self%ncid, self%v_vel_varname, &
                nf90_double, (/ self%node_dimid, self%time_dimid /), &
                varid=self%v_vel_varid)
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%v_vel_varid, &
                "long_name", "northward component of depth-averaged " &
                // "velocity")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%v_vel_varid, &
                "standard_name", "northward_water_velocity")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%v_vel_varid, &
                "positive", "north")
            call ncfile_check_error(ncstat)
        end if
        if (ics == 1 .or. ics == 2) then ! For all cases
            ncstat = nf90_put_att(self%ncid, self%v_vel_varid, &
                "coordinates", "time y x")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%v_vel_varid, &
                "location", "node")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%v_vel_varid, &
                "mesh", "adcirc_mesh")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%v_vel_varid, &
                "units", "m s-1")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%v_vel_varid, &
                "_FillValue", real_fill_value)
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%v_vel_varid, &
                "dry_Value", real_fill_value)
            call ncfile_check_error(ncstat)
        endif

        ! Define global attributes
        ! N/A, for now
    end procedure ncdg_64_set_metadata

    module procedure ncdg_64_write_step
        ! See interface for arguments and documentation

        integer :: ncstat ! Status of most recent operation
        integer :: np ! Number of nodes
        real :: t_arr(1) ! Vector structure for writing t
        integer :: time_start(1) ! Starting position for time data
        integer :: time_count(1) ! time data block size
        integer :: vel_start(2) ! Starting position for vel data
        integer :: vel_count(2) ! vel data block size
        integer :: record_count ! For indexing time step

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

        ! u_vel and v_vel
        vel_start(1) = 1 ! node dimension start
        vel_start(2) = record_count + 1 ! time dimension start
        vel_count(1) = np ! node dimension span
        vel_count(2) = 1 ! time dimension span

        ! u_vel
        ncstat = nf90_put_var(self%ncid, self%u_vel_varid, &
            u_vel, vel_start, vel_count)
        call ncfile_check_error(ncstat)

        ! v_vel
        ncstat = nf90_put_var(self%ncid, self%v_vel_varid, &
            v_vel, vel_start, vel_count)
        call ncfile_check_error(ncstat)

        ! After each write, sync
        ncstat = nf90_sync(self%ncid)
        call ncfile_check_error(ncstat)
    end procedure ncdg_64_write_step

end submodule ncdg_64_sub
