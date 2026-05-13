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
        
        ! integer :: ncstat ! Status of most recent operation

        ! N/A
    end procedure ncdg_63_write_step
    
end submodule ncdg_63_sub