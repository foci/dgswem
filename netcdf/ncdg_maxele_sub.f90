submodule (ncdg:ncdg_nodal_scalar_max_sub) ncdg_maxele_sub
!! Implementation of maxele.63.nc

    implicit none

    contains

    module procedure ncdg_maxele_create
        ! See interface for arguments and documentation

        character(len=:), allocatable :: the_path ! For defaulting
        integer :: the_cmode ! For defaulting

        ! Set path
        if (present(path)) then
            the_path = path
        else
            the_path = "maxele.63.nc"
        end if

        ! Set cmode
        if (present(cmode)) then
            the_cmode = cmode
        else
            the_cmode = ior(nf90_noclobber, nf90_netcdf4)
        end if

        ! Set scalar and time_of_scalar_max names
        self%scalar_max_varname = "zeta_max"
        self%time_of_scalar_max_varname = "time_of_zeta_max"

        ! Call parent function
        call self%ncfile%create(path=the_path, cmode=the_cmode)

        deallocate(the_path)
    end procedure ncdg_maxele_create

    module procedure ncdg_maxele_open
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

        ! Set scalar and time_of_scalar_max names
        self%scalar_max_varname = "zeta_max"
        self%time_of_scalar_max_varname = "time_of_zeta_max"

        ! Call parent function
        call self%ncdg_nodal_scalar_max%open(path=self%path, mode=the_mode)
    end procedure ncdg_maxele_open

    module procedure ncdg_maxele_set_metadata
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

        ! zeta_max
        ncstat = nf90_def_var(self%ncid, self%scalar_max_varname, &
            nf90_double, self%node_dimid, &
            varid=self%scalar_max_varid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%scalar_max_varid, &
            "long_name", "maximum water surface elevation above geoid")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%scalar_max_varid, &
            "standard_name", "maximum_sea_surface_height_above_geoid")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%scalar_max_varid, &
            "coordinates", "y x")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%scalar_max_varid, &
            "location", "node")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%scalar_max_varid, &
            "mesh", "adcirc_mesh")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%scalar_max_varid, &
            "units", "m")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%scalar_max_varid, &
            "_FillValue", real_fill_value)
        call ncfile_check_error(ncstat)

        ! time_of_zeta_max
        ncstat = nf90_def_var(self%ncid, self%time_of_scalar_max_varname, &
            nf90_double, self%node_dimid, &
            varid=self%time_of_scalar_max_varid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%time_of_scalar_max_varid, &
            "long_name", "time of maximum water surface elevation above geoid")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%time_of_scalar_max_varid, &
            "standard_name", "time_of_maximum_sea_surface_height_above_geoid")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%time_of_scalar_max_varid, &
            "coordinates", "y x")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%time_of_scalar_max_varid, &
            "location", "node")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%time_of_scalar_max_varid, &
            "mesh", "adcirc_mesh")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%time_of_scalar_max_varid, &
            "units", "s")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%time_of_scalar_max_varid, &
            "_FillValue", real_fill_value)
        call ncfile_check_error(ncstat)

        ! Define global attributes
        ! N/A, for now
    end procedure ncdg_maxele_set_metadata

end submodule ncdg_maxele_sub
