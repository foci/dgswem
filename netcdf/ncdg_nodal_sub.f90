submodule (ncdg:ncdg_file_sub) ncdg_nodal_sub
!! Implementation of functions common to all DGSWEM nodal output files.

    implicit none

    contains

    module procedure ncdg_nodal_open
        ! See interface for arguments and documentation
        
        integer :: ncstat ! Status of most recent operation

        ! Call parent function
        if (present(mode)) then
            call self%ncdg_file%open(mode)
        else
            call self%ncdg_file%open()
        endif

        ! Set dimension IDs
        ncstat = nf90_inq_dimid(self%ncid, self%node_dimname, &
            self%node_dimid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_inq_dimid(self%ncid, self%nele_dimname, &
            self%nele_dimid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_inq_dimid(self%ncid, self%nvertex_dimname, &
            self%nvertex_dimid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_inq_dimid(self%ncid, self%nope_dimname, &
            self%nope_dimid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_inq_dimid(self%ncid, self%neta_dimname, &
            self%neta_dimid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_inq_dimid(self%ncid, self%max_nvdll_dimname, &
            self%max_nvdll_dimid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_inq_dimid(self%ncid, self%nbou_dimname, &
            self%nbou_dimid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_inq_dimid(self%ncid, self%nvel_dimname, &
            self%nvel_dimid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_inq_dimid(self%ncid, self%max_nvell_dimname, &
            self%max_nvell_dimid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_inq_dimid(self%ncid, self%mesh_dimname, &
            self%mesh_dimid)
        call ncfile_check_error(ncstat)

        ! Set variable IDs
        ncstat = nf90_inq_varid(self%ncid, self%x_varname, &
            self%x_varid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_inq_varid(self%ncid, self%y_varname, &
            self%y_varid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_inq_varid(self%ncid, self%element_varname, &
            self%element_varid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_inq_varid(self%ncid, self%adcirc_mesh_varname, &
            self%adcirc_mesh_varid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_inq_varid(self%ncid, self%neta_varname, &
            self%neta_varid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_inq_varid(self%ncid, self%nvdll_varname, &
            self%nvdll_varid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_inq_varid(self%ncid, self%max_nvdll_varname, &
            self%max_nvdll_varid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_inq_varid(self%ncid, self%ibtypee_varname, &
            self%ibtypee_varid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_inq_varid(self%ncid, self%nbdv_varname, &
            self%nbdv_varid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_inq_varid(self%ncid, self%nvel_varname, &
            self%nvel_varid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_inq_varid(self%ncid, self%nvell_varname, &
            self%nvell_varid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_inq_varid(self%ncid, self%max_nvell_varname, &
            self%max_nvell_varid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_inq_varid(self%ncid, self%ibtype_varname, &
            self%ibtype_varid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_inq_varid(self%ncid, self%nbvv_varname, &
            self%nbvv_varid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_inq_varid(self%ncid, self%depth_varname, &
            self%depth_varid)
        call ncfile_check_error(ncstat)
    end procedure ncdg_nodal_open

    module procedure ncdg_nodal_close
        ! See interface for arguments and documentation
        
        ! integer :: ncstat ! Status of most recent operation

        ! Call parent function
        call self%ncdg_file%close()

        ! Flush dimension IDs
        self%node_dimid = -1
        self%nele_dimid = -1
        self%nvertex_dimid = -1
        self%nope_dimid = -1
        self%neta_dimid = -1
        self%max_nvdll_dimid = -1
        self%nbou_dimid = -1
        self%nvel_dimid = -1
        self%max_nvell_dimid = -1
        self%mesh_dimid = -1

        ! Flush variable IDs
        self%x_varid = -1
        self%y_varid = -1
        self%element_varid = -1
        self%adcirc_mesh_varid = -1
        self%neta_varid = -1
        self%nvdll_varid = -1
        self%max_nvdll_varid = -1
        self%ibtypee_varid = -1
        self%nbdv_varid = -1
        self%nvel_varid = -1
        self%nvell_varid = -1
        self%max_nvell_varid = -1
        self%ibtype_varid = -1
        self%nbvv_varid = -1
        self%depth_varid = -1
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