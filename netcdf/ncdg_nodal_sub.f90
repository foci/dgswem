submodule (ncdg : ncdg_file_sub) ncdg_nodal_sub
!! Implementation of functions common to all DGSWEM nodal output files.

    implicit none

    contains

    module procedure ncdg_nodal_open
        ! See interface for arguments and documentation
        
        ! integer :: ncstat ! Status of most recent operation

        ! Call parent function
        if (present(mode)) then
            call self%ncdg_file%open(mode)
        else
            call self%ncdg_file%open()
        endif

        ! Set dimension IDs
        ! N/A, for now

        ! Set variable IDs
        ! N/A, for now
    end procedure ncdg_nodal_open

    module procedure ncdg_nodal_close
        ! See interface for arguments and documentation
        
        ! integer :: ncstat ! Status of most recent operation

        ! Call parent function
        call self%ncdg_file%close()

        ! Flush dimension IDs
        ! N/A, for now

        ! Flush variable IDs
        ! N/A, for now
    end procedure ncdg_nodal_close

    module procedure ncdg_nodal_set_metadata
        ! See interface for arguments and documentation
        
        integer :: ncstat ! Status of most recent operation

        ! Call parent function
        call self%ncdg_file%ncdg_file_set_metadata( &
            time_dimsize = time_dimsize &
        )

        ! Define dimensions and their attributes
        ncstat = nf90_def_dim(self%ncid, self%node_dimname, &
            node_dimsize, self%node_dimid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_def_dim(self%ncid, self%nele_dimname, &
            nele_dimsize, self%nele_dimid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_def_dim(self%ncid, self%nvertex_dimname, &
            nvertex_dimsize, self%nvertex_dimid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_def_dim(self%ncid, self%nope_dimname, &
            nope_dimsize, self%nope_dimid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_def_dim(self%ncid, self%neta_dimname, &
            neta_dimsize, self%neta_dimid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_def_dim(self%ncid, self%max_nvdll_dimname, &
            max_nvdll_dimsize, self%max_nvdll_dimid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_def_dim(self%ncid, self%nbou_dimname, &
            nbou_dimsize, self%nbou_dimid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_def_dim(self%ncid, self%nvel_dimname, &
            nvel_dimsize, self%nvel_dimid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_def_dim(self%ncid, self%max_nvell_dimname, &
            max_nvell_dimsize, self%max_nvell_dimid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_def_dim(self%ncid, self%mesh_dimname, &
            mesh_dimsize, self%mesh_dimid)
        call ncfile_check_error(ncstat)

        ! Define variables and their attributes
        ncstat = nf90_def_var(self%ncid, self%x_varname, &
            nf90_double, self%node_dimid, varid=self%x_varid)
        call ncfile_check_error(ncstat)
        
        ncstat = nf90_def_var(self%ncid, self%y_varname, &
            nf90_double, self%node_dimid, varid=self%y_varid)
        call ncfile_check_error(ncstat)
        
        ncstat = nf90_def_var(self%ncid, self%element_varname, &
            nf90_int, (/ self%nvertex_dimid, self%nele_dimid /), &
            varid=self%element_varid)
        call ncfile_check_error(ncstat)
        
        ncstat = nf90_def_var(self%ncid, self%adcirc_mesh_varname, &
            nf90_int, self%mesh_dimid, varid=self%adcirc_mesh_varid)
        call ncfile_check_error(ncstat)
        
        ncstat = nf90_def_var(self%ncid, self%neta_varname, &
            nf90_int, varid=self%neta_varid)
        call ncfile_check_error(ncstat)
        
        ncstat = nf90_def_var(self%ncid, self%nvdll_varname, &
            nf90_int, self%nope_dimid, varid=self%nvdll_varid)
        call ncfile_check_error(ncstat)
        
        ncstat = nf90_def_var(self%ncid, self%max_nvdll_varname, &
            nf90_int, varid=self%max_nvdll_varid)
        call ncfile_check_error(ncstat)
        
        ncstat = nf90_def_var(self%ncid, self%ibtypee_varname, &
            nf90_int, self%nope_dimid, varid=self%ibtypee_varid)
        call ncfile_check_error(ncstat)
        
        ncstat = nf90_def_var(self%ncid, self%nbdv_varname, &
            nf90_int, self%neta_dimid, varid=self%nbdv_varid)
        call ncfile_check_error(ncstat)
        
        ncstat = nf90_def_var(self%ncid, self%nvel_varname, &
            nf90_int, varid=self%nvel_varid)
        call ncfile_check_error(ncstat)
        
        ncstat = nf90_def_var(self%ncid, self%nvell_varname, &
            nf90_int, self%nbou_dimid, varid=self%nvell_varid)
        call ncfile_check_error(ncstat)
        
        ncstat = nf90_def_var(self%ncid, self%max_nvell_varname, &
            nf90_int, varid=self%max_nvell_varid)
        call ncfile_check_error(ncstat)
        
        ncstat = nf90_def_var(self%ncid, self%ibtype_varname, &
            nf90_int, self%nbou_dimid, varid=self%ibtype_varid)
        call ncfile_check_error(ncstat)
        
        ncstat = nf90_def_var(self%ncid, self%nbvv_varname, &
            nf90_int, self%nvel_dimid, varid=self%nbvv_varid)
        call ncfile_check_error(ncstat)
        
        ncstat = nf90_def_var(self%ncid, self%depth_varname, &
            nf90_double, self%node_dimid, varid=self%depth_varid)
        call ncfile_check_error(ncstat)
        
        ! Define attributes
        ! N/A, for now
    end procedure ncdg_nodal_set_metadata

end submodule ncdg_nodal_sub