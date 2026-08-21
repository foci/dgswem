submodule (ncdg:ncdg_nodal_vector_sub) ncdg_vel_sub
!! Implementation of fort.64.nc

    implicit none

    contains

    module procedure ncdg_vel_create
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
            the_cmode = ior(nf90_noclobber, nf90_netcdf4)
        end if

        ! Set vector names
        self%vector_u_varname = "u_vel"
        self%vector_v_varname = "v_vel"

        ! Call parent function
        call self%ncfile%create(path=the_path, cmode=the_cmode)

        deallocate(the_path)
    end procedure ncdg_vel_create

    module procedure ncdg_vel_open
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

        ! Set vector names
        self%vector_u_varname = "u_vel"
        self%vector_v_varname = "v_vel"

        ! Call parent function
        call self%ncdg_nodal_vector%open(path=self%path, mode=the_mode)
    end procedure ncdg_vel_open

    module procedure ncdg_vel_set_metadata
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
            ncstat = nf90_def_var(self%ncid, self%vector_u_varname, &
                nf90_double, (/ self%node_dimid, self%time_dimid /), &
                varid=self%vector_u_varid)
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%vector_u_varid, &
                "long_name", "x component of depth-averaged velocity")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%vector_u_varid, &
                "standard_name", "sea_water_x_velocity")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%vector_u_varid, &
                "positive", "right")
            call ncfile_check_error(ncstat)
        else if (ics == 2) then ! Spherical
            ncstat = nf90_def_var(self%ncid, self%vector_u_varname, &
                nf90_double, (/ self%node_dimid, self%time_dimid /), &
                varid=self%vector_u_varid)
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%vector_u_varid, &
                "long_name", "eastward component of depth-averaged " &
                // "velocity")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%vector_u_varid, &
                "standard_name", "eastward_water_velocity")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%vector_u_varid, &
                "positive", "east")
            call ncfile_check_error(ncstat)
        end if
        if (ics == 1 .or. ics == 2) then ! For all cases
            ncstat = nf90_put_att(self%ncid, self%vector_u_varid, &
                "coordinates", "time y x")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%vector_u_varid, &
                "location", "node")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%vector_u_varid, &
                "mesh", "adcirc_mesh")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%vector_u_varid, &
                "units", "m s-1")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%vector_u_varid, &
                "_FillValue", real_fill_value)
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%vector_u_varid, &
                "dry_Value", real_fill_value)
            call ncfile_check_error(ncstat)
        end if

        ! v_vel
        if (ics == 1) then ! Cartesian
            ncstat = nf90_def_var(self%ncid, self%vector_v_varname, &
                nf90_double, (/ self%node_dimid, self%time_dimid /), &
                varid=self%vector_v_varid)
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%vector_v_varid, &
                "long_name", "y component of depth-averaged velocity")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%vector_v_varid, &
                "standard_name", "sea_water_y_velocity")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%vector_v_varid, &
                "positive", "90 degrees counterclockwise from x " &
                // "water velocity")
            call ncfile_check_error(ncstat)
        else if (ics == 2) then ! Spherical
            ncstat = nf90_def_var(self%ncid, self%vector_v_varname, &
                nf90_double, (/ self%node_dimid, self%time_dimid /), &
                varid=self%vector_v_varid)
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%vector_v_varid, &
                "long_name", "northward component of depth-averaged " &
                // "velocity")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%vector_v_varid, &
                "standard_name", "northward_water_velocity")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%vector_v_varid, &
                "positive", "north")
            call ncfile_check_error(ncstat)
        end if
        if (ics == 1 .or. ics == 2) then ! For all cases
            ncstat = nf90_put_att(self%ncid, self%vector_v_varid, &
                "coordinates", "time y x")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%vector_v_varid, &
                "location", "node")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%vector_v_varid, &
                "mesh", "adcirc_mesh")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%vector_v_varid, &
                "units", "m s-1")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%vector_v_varid, &
                "_FillValue", real_fill_value)
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%vector_v_varid, &
                "dry_Value", real_fill_value)
            call ncfile_check_error(ncstat)
        end if

        ! Define global attributes
        ! N/A, for now
    end procedure ncdg_vel_set_metadata

end submodule ncdg_vel_sub
