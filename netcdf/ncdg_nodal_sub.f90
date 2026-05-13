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

        ! Define dimensions
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
        
        ! x
        ncstat = nf90_def_var(self%ncid, self%x_varname, &
            nf90_double, self%node_dimid, varid=self%x_varid)
        call ncfile_check_error(ncstat)
        if (ics == 1) then ! Cartesian
            ncstat = nf90_put_att(self%ncid, self%x_varid, &
                "long_name", "Cartesian coordinate x")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%x_varid, &
                "standard_name", "projection_x_coordinate")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%x_varid, &
                "units", "m")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%x_varid, &
                "positive", "right")
            call ncfile_check_error(ncstat)
        else if (ics == 2) then ! Spherical
            ncstat = nf90_put_att(self%ncid, self%x_varid, &
                "long_name", "longitude")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%x_varid, &
                "standard_name", "longitude")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%x_varid, &
                "units", "degrees_east")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%x_varid, &
                "positive", "east")
            call ncfile_check_error(ncstat)
        end if

        ! y
        ncstat = nf90_def_var(self%ncid, self%y_varname, &
            nf90_double, self%node_dimid, varid=self%y_varid)
        call ncfile_check_error(ncstat)
        if (ics == 1) then ! Cartesian
            ncstat = nf90_put_att(self%ncid, self%y_varid, &
                "long_name", "Cartesian coordinate y")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%y_varid, &
                "standard_name", "projection_y_coordinate")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%y_varid, &
                "units", "m")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%y_varid, &
                "positive", "90 degrees counterclockwise from x")
            call ncfile_check_error(ncstat)
        else if (ics == 2) then ! Spherical
            ncstat = nf90_put_att(self%ncid, self%y_varid, &
                "long_name", "latitude")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%y_varid, &
                "standard_name", "latitude")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%y_varid, &
                "units", "degrees_north")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%y_varid, &
                "positive", "north")
            call ncfile_check_error(ncstat)
        end if

        ! element
        ncstat = nf90_def_var(self%ncid, self%element_varname, &
            nf90_int, (/ self%nvertex_dimid, self%nele_dimid /), &
            varid=self%element_varid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%element_varid, &
            "long_name", "element")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%element_varid, &
            "cf_role", "face_node_connectivity")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%element_varid, &
            "start_index", 1)
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%element_varid, &
            "units", "nondimensional")
        call ncfile_check_error(ncstat)

        ! adcirc_mesh
        ncstat = nf90_def_var(self%ncid, self%adcirc_mesh_varname, &
            nf90_int, self%mesh_dimid, varid=self%adcirc_mesh_varid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%adcirc_mesh_varid, &
            "long_name", "mesh_topology")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%adcirc_mesh_varid, &
            "cf_role", "mesh_topology")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%adcirc_mesh_varid, &
            "topology_dimension", "2")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%adcirc_mesh_varid, &
            "node_coordinates", "x y")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%adcirc_mesh_varid, &
            "face_node_connectivity", "element")
        call ncfile_check_error(ncstat)

        ! Additional 
        if (ics == 1) then
            print *, "Cartesian coordinate system detected. GIS " &
                // "rendering may fail."
        end if

        ! Elevation boundary variables
        if (nope /= 0) then
            ! neta
            ncstat = nf90_def_var(self%ncid, self%neta_varname, &
                nf90_int, varid=self%neta_varid)
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%neta_varid, &
                "long_name", "total numbar of elevation specified " &
                // "boundary nodes")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%neta_varid, &
                "units", "nondimensional")
            call ncfile_check_error(ncstat)
    
            ! nvdll
            ncstat = nf90_def_var(self%ncid, self%nvdll_varname, &
                nf90_int, self%nope_dimid, varid=self%nvdll_varid)
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%nvdll_varid, &
                "long_name", "number of nodes in each elevation " &
                // "specified boundary segment")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%nvdll_varid, &
                "units", "nondimensional")
            call ncfile_check_error(ncstat)
    
            ! max_nvdll
            ncstat = nf90_def_var(self%ncid, self%max_nvdll_varname, &
                nf90_int, varid=self%max_nvdll_varid)
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%max_nvdll_varid, &
                "long_name", "maximum number of nodes in any " &
                // "elevation specified boundary segment")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%max_nvdll_varid, &
                "units", "nondimensional")
            call ncfile_check_error(ncstat)
    
            ! ibtypee
            ncstat = nf90_def_var(self%ncid, self%ibtypee_varname, &
                nf90_int, self%nope_dimid, varid=self%ibtypee_varid)
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%ibtypee_varid, &
                "long_name", "type of each elevation speficied boundary " &
                // "segment")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%ibtypee_varid, &
                "units", "nondimensional")
            call ncfile_check_error(ncstat)
    
            ! nbdv
            ncstat = nf90_def_var(self%ncid, self%nbdv_varname, &
                nf90_int, self%neta_dimid, varid=self%nbdv_varid)
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%nbdv_varid, &
                "long_name", "node numers on each elevation " &
                // "specified boundary segment")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%nbdv_varid, &
                "units", "nondimensional")
            call ncfile_check_error(ncstat)
        end if

        ! Normal flow (discharge) boundary variables
        if (nbou /= 0) then
            ! nvel
            ncstat = nf90_def_var(self%ncid, self%nvel_varname, &
                nf90_int, varid=self%nvel_varid)
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%nvel_varid, &
                "long_name", "total number of normal flow " &
                // "(discharge) specified boundary nodes")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%nvel_varid, &
                "units", "nondimensional")
            call ncfile_check_error(ncstat)
    
            ! nvell
            ncstat = nf90_def_var(self%ncid, self%nvell_varname, &
                nf90_int, self%nbou_dimid, varid=self%nvell_varid)
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%nvell_varid, &
                "long_name", "number of nodes in each normal flow " &
                // "(discharge) specified boundary segment")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%nvell_varid, &
                "units", "nondimensional")
            call ncfile_check_error(ncstat)
    
            ! max_nvell
            ncstat = nf90_def_var(self%ncid, self%max_nvell_varname, &
                nf90_int, varid=self%max_nvell_varid)
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%max_nvell_varid, &
                "long_name", "maximum number of nodes in any normal " &
                // "flow (discharge) specified boundary segment")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%max_nvell_varid, &
                "units", "nondimensional")
            call ncfile_check_error(ncstat)
    
            ! ibtype
            ncstat = nf90_def_var(self%ncid, self%ibtype_varname, &
                nf90_int, self%nbou_dimid, varid=self%ibtype_varid)
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%ibtype_varid, &
                "long_name", "type of each normal flow (discharge) " &
                // " specified boundary segment")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%ibtype_varid, &
                "units", "nondimensional")
            call ncfile_check_error(ncstat)
    
            ! nbvv
            ncstat = nf90_def_var(self%ncid, self%nbvv_varname, &
                nf90_int, self%nvel_dimid, varid=self%nbvv_varid)
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%nbvv_varid, &
                "long_name", "node numbers on each normal flow " &
                // "(discharge) specified boundary segment")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, self%nbvv_varid, &
                "units", "nondimensional")
            call ncfile_check_error(ncstat)
        end if

        ! depth
        ncstat = nf90_def_var(self%ncid, self%depth_varname, &
            nf90_double, self%node_dimid, varid=self%depth_varid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%depth_varid, &
            "long_name", "distance below geoid")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%depth_varid, &
            "standard_name", "depth below geoid")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%depth_varid, &
            "coordinates", "time y x")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%depth_varid, &
            "location", "node")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%depth_varid, &
            "mesh", "adcirc_mesh")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%depth_varid, &
            "units", "m")
        call ncfile_check_error(ncstat)
        
        ! Define global attributes
        ! N/A
    end procedure ncdg_nodal_set_metadata

    module procedure ncdg_nodal_write_mesh
        ! See interface for arguments and documentation
        
        integer :: ncstat ! Status of most recent operation
        real :: element(nhy, ne) ! Transpose of nm
        integer :: element_start(2) ! Starting position for element data
        integer :: element_count(2) ! Element data block size
        real :: nbdv2(neta) ! Rearrangement of nbdv for compactness
        real :: nbvv2(nvel) ! Rearrangement of nbvv for compactness
        integer :: i, j, k, kj ! Indices

        ! Exit define mode
        ncstat = nf90_enddef(self%ncid)
        if (ncstat /= nf90_enotindefine) then
            call ncfile_check_error(ncstat)
        end if

        print *, "re-arranging data arrays"

        ! Re-arrange data arrays
        ! nm -> element
        do i = 1, ne
            do j = 1, nhy
                element(j, i) = nm(i, j)
            end do
        end do
        ! nbdv -> nbdv2
        kj = 1
        do k = 1, nope
            do j = 1, nvdll(k)
                nbdv2(kj) = nbdv(k, j)
                kj = kj + 1
            end do
        end do
        ! nbvv -> nbvv2
        kj = 1
        do k = 1, nbou
            do j = 1, nvell(k)
                nbvv2(kj) = nbvv(k, j)
                kj = kj + 1
            end do
        end do

        ! Coordinates
        ncstat = nf90_put_var(self%ncid, self%x_varid, x)
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_var(self%ncid, self%y_varid, y)
        call ncfile_check_error(ncstat)

        ! Bathymetry
        ncstat = nf90_put_var(self%ncid, self%depth_varid, dp)
        call ncfile_check_error(ncstat)

        ! Elements
        element_start(1) = 1
        element_start(2) = 1
        element_count(1) = nhy
        element_count(2) = ne
        ncstat = nf90_put_var(self%ncid, self%element_varid, element, &
            element_start, element_count)
        call ncfile_check_error(ncstat)
        
        ! Elevation boundaries
        if (nope /= 0) then
            ! Boundary segment types
            ncstat = nf90_put_var(self%ncid, self%ibtypee_varid, ibtypee)
            call ncfile_check_error(ncstat)
            ! Number of nodes per segment
            ncstat = nf90_put_var(self%ncid, self%nvdll_varid, nvdll)
            call ncfile_check_error(ncstat)
            ! Node numbers for each segment
            ncstat = nf90_put_var(self%ncid, self%nbdv_varid, nbdv2)
            call ncfile_check_error(ncstat)
        endif
        
        ! Normal flow (discharge) boundaries
        if (nbou /= 0) then
            ! Boundary segment types
            ncstat = nf90_put_var(self%ncid, self%ibtype_varid, ibtype)
            call ncfile_check_error(ncstat)
            ! Number of nodes per segment
            ncstat = nf90_put_var(self%ncid, self%nvell_varid, nvell)
            call ncfile_check_error(ncstat)
            ! Node numbers for each segment
            ncstat = nf90_put_var(self%ncid, self%nbvv_varid, nbvv2)
            call ncfile_check_error(ncstat)
            ! Weirs and pipes would be implemented here
        endif
    end procedure ncdg_nodal_write_mesh

end submodule ncdg_nodal_sub