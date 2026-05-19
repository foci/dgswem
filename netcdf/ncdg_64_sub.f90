submodule (ncdg:ncdg_nodal_sub) ncdg_64_sub
!! Implementation of fort.64.nc

    implicit none

    contains

    module procedure ncdg_64_open
        ! See interface for arguments and documentation
        
        integer :: ncstat ! Status of most recent operation

        ! Call parent function
        if (present(mode)) then
            call self%ncdg_nodal%open(mode)
        else
            call self%ncdg_nodal%open()
        endif

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
            time_dimsize = time_dimsize, &
            node_dimsize = node_dimsize, &
            nele_dimsize = nele_dimsize, &
            nvertex_dimsize = nvertex_dimsize, &
            nope_dimsize = nope_dimsize, &
            neta_dimsize = neta_dimsize, &
            max_nvdll_dimsize = max_nvdll_dimsize, &
            nbou_dimsize = nbou_dimsize, &
            nvel_dimsize = nvel_dimsize, &
            max_nvell_dimsize = max_nvell_dimsize, &
            mesh_dimsize = mesh_dimsize, &
            nope = nope, &
            nbou = nbou, &
            ics = ics &
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
                "long_name", "water column vertically averaged " &
                // "in x-direction")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%u_vel_varid, &
                "standard_name", "x_water_velocity_depth_averaged")
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
                "long_name", "water column vertically averaged " &
                // "east/west velocity")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%u_vel_varid, &
                "standard_name", "eastward_water_velocity")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%u_vel_varid, &
                "positive", "east")
            call ncfile_check_error(ncstat)
        end if
        if (ics == 0 .or. ics == 1) then ! For all cases
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
                "long_name", "water column vertically averaged " &
                // "in y-direction")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%v_vel_varid, &
                "standard_name", "y_water_velocity_depth_averaged")
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
                "long_name", "water column vertically averaged " &
                // "north/south velocity")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%v_vel_varid, &
                "standard_name", "northward_water_velocity")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%v_vel_varid, &
                "positive", "north")
            call ncfile_check_error(ncstat)
        end if
        if (ics == 0 .or. ics == 1) then ! For all cases
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