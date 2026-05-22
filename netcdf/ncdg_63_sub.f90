submodule (ncdg:ncdg_nodal_sub) ncdg_63_sub
!! Implementation of fort.63.nc

    implicit none

    contains

    module procedure ncdg_63_open
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
        ncstat = nf90_inq_varid(self%ncid, self%zeta_varname, &
            self%zeta_varid)
        call ncfile_check_error(ncstat)
    end procedure ncdg_63_open

    module procedure ncdg_63_close
        ! See interface for arguments and documentation
        
        ! integer :: ncstat ! Status of most recent operation

        ! Call parent function
        call self%ncdg_nodal%close()

        ! Flush dimension IDs
        ! N/A

        ! Flush variable IDs
        self%zeta_varid = -1
    end procedure ncdg_63_close

    module procedure ncdg_63_set_metadata
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
        ncstat = nf90_def_var(self%ncid, self%zeta_varname, &
            nf90_double, (/ self%node_dimid, self%time_dimid /), &
            varid=self%zeta_varid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%zeta_varid, &
            "long_name", "water surface elevation above geoid")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%zeta_varid, &
            "standard_name", "sea_surface_height_above_geoid")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%zeta_varid, &
            "coordinates", "time y x")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%zeta_varid, &
            "location", "node")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%zeta_varid, &
            "mesh", "adcirc_mesh")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%zeta_varid, &
            "units", "m")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%zeta_varid, &
            "_FillValue", real_fill_value)
        call ncfile_check_error(ncstat)
        
        ! Define global attributes
        ! N/A, for now
    end procedure ncdg_63_set_metadata

    module procedure ncdg_63_write_step
        ! See interface for arguments and documentation
        
        integer :: ncstat ! Status of most recent operation
        integer :: np ! Number of nodes
        real :: t_arr(1) ! Vector structure for writing t
        integer :: time_start(1) ! Starting position for time data
        integer :: time_count(1) ! time data block size
        integer :: zeta_start(2) ! Starting position for zeta data
        integer :: zeta_count(2) ! zeta data block size
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

        ! zeta
        zeta_start(1) = 1 ! node dimension start
        zeta_start(2) = record_count + 1 ! time dimension start
        zeta_count(1) = np ! node dimension span
        zeta_count(2) = 1 ! time dimension span
        ncstat = nf90_put_var(self%ncid, self%zeta_varid, &
            zeta, zeta_start, zeta_count)
        call ncfile_check_error(ncstat)

        ! After each write, sync
        ncstat = nf90_sync(self%ncid)
        call ncfile_check_error(ncstat)
    end procedure ncdg_63_write_step
    
end submodule ncdg_63_sub