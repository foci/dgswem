submodule (ncdg:ncdg_nodal_sub) ncdg_maxele_sub
!! Implementation of fort.63.nc

    implicit none

    contains

    module procedure ncdg_maxele_init
        ! See interface for arguments and documentation

        character(len=:), allocatable :: the_path ! For defaulting

        ! Set path
        if (present(path)) then
            the_path = path
        else
            the_path = "maxele.63.nc"
        end if

        ! Call parent function
        if (present(cmode)) then
            call self%ncfile%init(path=the_path, cmode=cmode)
        else
            call self%ncfile%init(path=the_path)
        end if

        deallocate(the_path)
    end procedure ncdg_maxele_init

    module procedure ncdg_maxele_open
        ! See interface for arguments and documentation

        integer :: ncstat ! Status of most recent operation

        ! Call parent function
        if (present(mode)) then
            call self%ncdg_nodal%open(mode=mode)
        else
            call self%ncdg_nodal%open()
        end if

        ! Set dimension IDs
        ! N/A

        ! Set variable IDs
        ncstat = nf90_inq_varid(self%ncid, self%zeta_max_varname, &
            self%zeta_max_varid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_inq_varid(self%ncid, self%time_of_zeta_max_varname, &
            self%time_of_zeta_max_varid)
        call ncfile_check_error(ncstat)
    end procedure ncdg_maxele_open

    module procedure ncdg_maxele_close
        ! See interface for arguments and documentation

        ! integer :: ncstat ! Status of most recent operation

        ! Call parent function
        call self%ncdg_nodal%close()

        ! Flush dimension IDs
        ! N/A

        ! Flush variable IDs
        self%zeta_max_varid = -1
        self%time_of_zeta_max_varid = -1
    end procedure ncdg_maxele_close

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
        ncstat = nf90_def_var(self%ncid, self%zeta_max_varname, &
            nf90_double, self%node_dimid, &
            varid=self%zeta_max_varid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%zeta_max_varid, &
            "long_name", "maximum water surface elevation above geoid")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%zeta_max_varid, &
            "standard_name", "maximum_sea_surface_height_above_geoid")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%zeta_max_varid, &
            "coordinates", "y x")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%zeta_max_varid, &
            "location", "node")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%zeta_max_varid, &
            "mesh", "adcirc_mesh")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%zeta_max_varid, &
            "units", "m")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%zeta_max_varid, &
            "_FillValue", real_fill_value)
        call ncfile_check_error(ncstat)

        ! time_of_zeta_max
        ncstat = nf90_def_var(self%ncid, self%time_of_zeta_max_varname, &
            nf90_double, self%node_dimid, &
            varid=self%time_of_zeta_max_varid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%time_of_zeta_max_varid, &
            "long_name", "time of maximum water surface elevation above geoid")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%time_of_zeta_max_varid, &
            "standard_name", "time_of_maximum_sea_surface_height_above_geoid")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%time_of_zeta_max_varid, &
            "coordinates", "y x")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%time_of_zeta_max_varid, &
            "location", "node")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%time_of_zeta_max_varid, &
            "mesh", "adcirc_mesh")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%time_of_zeta_max_varid, &
            "units", "s")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%time_of_zeta_max_varid, &
            "_FillValue", real_fill_value)
        call ncfile_check_error(ncstat)

        ! Define global attributes
        ! N/A, for now
    end procedure ncdg_maxele_set_metadata

    module procedure ncdg_maxele_write_step
         ! See interface for arguments and documentation

        integer :: ncstat ! Status of most recent operation
        integer :: np ! Number of nodes
        real, allocatable :: zeta_max_old(:) ! Previous maximums,
                                             ! if present
        real, allocatable :: time_of_zeta_max(:) ! Nodewise times when
                                                 ! zeta_max is achieved

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
        allocate(zeta_max_old(np))
        allocate(time_of_zeta_max(np))
        ncstat = nf90_get_var(self%ncid, self%zeta_max_varid, &
            zeta_max_old)
        call ncfile_check_error(ncstat)
        ncstat = nf90_get_var(self%ncid, self%time_of_zeta_max_varid, &
            time_of_zeta_max)
        call ncfile_check_error(ncstat)

        ! zeta_max
        ncstat = nf90_put_var(self%ncid, self%zeta_max_varid, &
            zeta_max)
        call ncfile_check_error(ncstat)

        ! time_of_zeta_max
        time_of_zeta_max = merge( &
            spread(t, 1, np), & ! Used if true: current time
            time_of_zeta_max, & ! Used if false: old time
            zeta_max > zeta_max_old & ! Condition: new maximum is strictly greater
        )
        ncstat = nf90_put_var(self%ncid, self%time_of_zeta_max_varid, &
            time_of_zeta_max)
        call ncfile_check_error(ncstat)
        deallocate(zeta_max_old)
        deallocate(time_of_zeta_max)

        ! After each write, sync
        ncstat = nf90_sync(self%ncid)
        call ncfile_check_error(ncstat)
    end procedure ncdg_maxele_write_step

end submodule ncdg_maxele_sub
