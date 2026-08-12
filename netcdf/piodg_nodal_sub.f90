submodule (piodg:piodg_file_sub) piodg_nodal_sub
    !! Implementation of functions common to all DGSWEM nodal output files.

    implicit none

    contains

    module procedure piodg_nodal_open
        ! See interface for arguments and documentation

        integer :: piostat ! Status of most recent operation
        integer :: the_piotype ! For defaulting
        integer :: the_omode ! For defaulting

        ! Set path
        ! Default is whatever is already there
        if (present(path)) then
            self%path = path
        else if (.not. allocated(self%path)) then
            error stop "PIO util error: Tried to open file with no path."
        end if

        ! Set piotype
        if (present(piotype)) then
            the_piotype = piotype
        else
            the_piotype = pio_iotype_netcdf4p
        end if

        ! Set omode
        if (present(omode)) then
            the_omode = omode
        else
            the_omode = pio_nowrite
        end if

        ! Call parent function
        call self%piodg_file%open(piosystem, self%path, the_piotype, the_omode)

        ! Set dimension IDs
        piostat = pio_inq_dimid(self%piofiledesc, self%node_dimname, &
            self%node_dimid)
        call piofile_check_error(piostat)
        piostat = pio_inq_dimid(self%piofiledesc, self%nele_dimname, &
            self%nele_dimid)
        call piofile_check_error(piostat)
        piostat = pio_inq_dimid(self%piofiledesc, self%nvertex_dimname, &
            self%nvertex_dimid)
        call piofile_check_error(piostat)
        piostat = pio_inq_dimid(self%piofiledesc, self%nope_dimname, &
            self%nope_dimid)
        call piofile_check_error(piostat)
        piostat = pio_inq_dimid(self%piofiledesc, self%neta_dimname, &
            self%neta_dimid)
        call piofile_check_error(piostat)
        piostat = pio_inq_dimid(self%piofiledesc, self%max_nvdll_dimname, &
            self%max_nvdll_dimid)
        call piofile_check_error(piostat)
        piostat = pio_inq_dimid(self%piofiledesc, self%nbou_dimname, &
            self%nbou_dimid)
        call piofile_check_error(piostat)
        piostat = pio_inq_dimid(self%piofiledesc, self%nvel_dimname, &
            self%nvel_dimid)
        call piofile_check_error(piostat)
        piostat = pio_inq_dimid(self%piofiledesc, self%max_nvell_dimname, &
            self%max_nvell_dimid)
        call piofile_check_error(piostat)
        piostat = pio_inq_dimid(self%piofiledesc, self%mesh_dimname, &
            self%mesh_dimid)
        call piofile_check_error(piostat)

        ! Set variable IDs
        piostat = pio_inq_varid(self%piofiledesc, self%x_varname, &
            self%x_vardesc)
        call piofile_check_error(piostat)
        piostat = pio_inq_varid(self%piofiledesc, self%y_varname, &
            self%y_vardesc)
        call piofile_check_error(piostat)
        piostat = pio_inq_varid(self%piofiledesc, self%element_varname, &
            self%element_vardesc)
        call piofile_check_error(piostat)
        piostat = pio_inq_varid(self%piofiledesc, self%adcirc_mesh_varname, &
            self%adcirc_mesh_vardesc)
        call piofile_check_error(piostat)
        piostat = pio_inq_varid(self%piofiledesc, self%neta_varname, &
            self%neta_vardesc)
        call piofile_check_error(piostat)
        piostat = pio_inq_varid(self%piofiledesc, self%nvdll_varname, &
            self%nvdll_vardesc)
        call piofile_check_error(piostat)
        piostat = pio_inq_varid(self%piofiledesc, self%max_nvdll_varname, &
            self%max_nvdll_vardesc)
        call piofile_check_error(piostat)
        piostat = pio_inq_varid(self%piofiledesc, self%ibtypee_varname, &
            self%ibtypee_vardesc)
        call piofile_check_error(piostat)
        piostat = pio_inq_varid(self%piofiledesc, self%nbdv_varname, &
            self%nbdv_vardesc)
        call piofile_check_error(piostat)
        piostat = pio_inq_varid(self%piofiledesc, self%nvel_varname, &
            self%nvel_vardesc)
        call piofile_check_error(piostat)
        piostat = pio_inq_varid(self%piofiledesc, self%nvell_varname, &
            self%nvell_vardesc)
        call piofile_check_error(piostat)
        piostat = pio_inq_varid(self%piofiledesc, self%max_nvell_varname, &
            self%max_nvell_vardesc)
        call piofile_check_error(piostat)
        piostat = pio_inq_varid(self%piofiledesc, self%ibtype_varname, &
            self%ibtype_vardesc)
        call piofile_check_error(piostat)
        piostat = pio_inq_varid(self%piofiledesc, self%nbvv_varname, &
            self%nbvv_vardesc)
        call piofile_check_error(piostat)
        piostat = pio_inq_varid(self%piofiledesc, self%depth_varname, &
            self%depth_vardesc)
        call piofile_check_error(piostat)
    end procedure piodg_nodal_open

    module procedure piodg_nodal_close
        ! See interface for arguments and documentation

        type(var_desc_t) :: x_vardesc_unset ! Blank struct
        type(var_desc_t) :: y_vardesc_unset ! Blank struct
        type(var_desc_t) :: element_vardesc_unset ! Blank struct
        type(var_desc_t) :: adcirc_mesh_vardesc_unset ! Blank struct
        type(var_desc_t) :: neta_vardesc_unset ! Blank struct
        type(var_desc_t) :: nvdll_vardesc_unset ! Blank struct
        type(var_desc_t) :: max_nvdll_vardesc_unset ! Blank struct
        type(var_desc_t) :: ibtypee_vardesc_unset ! Blank struct
        type(var_desc_t) :: nbdv_vardesc_unset ! Blank struct
        type(var_desc_t) :: nvel_vardesc_unset ! Blank struct
        type(var_desc_t) :: nvell_vardesc_unset ! Blank struct
        type(var_desc_t) :: max_nvell_vardesc_unset ! Blank struct
        type(var_desc_t) :: ibtype_vardesc_unset ! Blank struct
        type(var_desc_t) :: nbvv_vardesc_unset ! Blank struct
        type(var_desc_t) :: depth_vardesc_unset ! Blank struct

        ! Call parent function
        call self%piodg_file%close()

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
        self%x_vardesc = x_vardesc_unset
        self%y_vardesc = y_vardesc_unset
        self%element_vardesc = element_vardesc_unset
        self%adcirc_mesh_vardesc = adcirc_mesh_vardesc_unset
        self%neta_vardesc = neta_vardesc_unset
        self%nvdll_vardesc = nvdll_vardesc_unset
        self%max_nvdll_vardesc = max_nvdll_vardesc_unset
        self%ibtypee_vardesc = ibtypee_vardesc_unset
        self%nbdv_vardesc = nbdv_vardesc_unset
        self%nvel_vardesc = nvel_vardesc_unset
        self%nvell_vardesc = nvell_vardesc_unset
        self%max_nvell_vardesc = max_nvell_vardesc_unset
        self%ibtype_vardesc = ibtype_vardesc_unset
        self%nbvv_vardesc = nbvv_vardesc_unset
        self%depth_vardesc = depth_vardesc_unset
    end procedure piodg_nodal_close

end submodule piodg_nodal_sub
