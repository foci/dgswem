submodule (ncdg:ncdg_nodal_scalar_sub) ncdg_tracer_sub
!! Implementation of fort.63.nc

    implicit none

    contains

    module procedure ncdg_tracer_create
        ! See interface for arguments and documentation

        character(len=:), allocatable :: the_path ! For defaulting
        integer :: the_cmode ! For defaulting

        ! Set path
        if (present(path)) then
            the_path = path
        else
            the_path = "fort.88.nc"
        end if

        ! Set cmode
        if (present(cmode)) then
            the_cmode = cmode
        else
            the_cmode = ior(nf90_noclobber, nf90_netcdf4)
        end if

        ! Set scalar name
        self%scalar_varname = "concentration"

        ! Call parent function
        call self%ncfile%create(path=the_path, cmode=the_cmode)

        deallocate(the_path)
    end procedure ncdg_tracer_create

    module procedure ncdg_tracer_open
        ! See interface for arguments and documentation

        integer :: ncstat ! Status of most recent operation
        integer :: the_mode ! For defaulting

        ! Set path
        ! Default is whatever is already there
        if (present(path)) then
            self%path = path
        else if (.not. allocated(self%path)) then
            self%path = "fort.88.nc"
        end if

        ! Set mode
        if (present(mode)) then
            the_mode = mode
        else
            the_mode = nf90_nowrite
        end if

        ! Set scalar name
        self%scalar_varname = "concentration"

        ! Call parent function
        call self%ncdg_nodal_scalar%open(path=self%path, mode=the_mode)
    end procedure ncdg_tracer_open

    module procedure ncdg_tracer_set_metadata
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

        ! zeta
        ncstat = nf90_def_var(self%ncid, self%scalar_varname, &
            nf90_double, (/ self%node_dimid, self%time_dimid /), &
            varid=self%scalar_varid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%scalar_varid, &
            "long_name", "tracer concentration")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%scalar_varid, &
            "standard_name", "mass_concentration_of_tracer_in_sea_water")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%scalar_varid, &
            "coordinates", "time y x")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%scalar_varid, &
            "location", "node")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%scalar_varid, &
            "mesh", "adcirc_mesh")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%scalar_varid, &
            "units", "m")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%scalar_varid, &
            "_FillValue", real_fill_value)
        call ncfile_check_error(ncstat)

        ! Define global attributes
        ! N/A, for now
    end procedure ncdg_tracer_set_metadata

end submodule ncdg_tracer_sub
