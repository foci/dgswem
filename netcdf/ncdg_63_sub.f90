submodule (ncdg : ncdg_nodal_sub) ncdg_63_sub
!! Implementation of fort.63.nc

    implicit none

    contains

    module procedure ncdg_63_set_metadata
        ! See interface for arguments and documentation
        
        integer :: stat ! Status of most recent operation
        
        ! Dimension IDs
        integer :: time_dimid ! ID of time dimension
        integer :: node_dimid ! ID of node dimension
        integer :: nele_dimid ! ID of nele dimension
        integer :: nvertex_dimid ! ID of nvertex dimension
        integer :: nope_dimid ! ID of nope dimension
        integer :: neta_dimid ! ID of neta dimension
        integer :: max_nvdll_dimid ! ID of max_nvdll dimension
        integer :: nbou_dimid ! ID of nbou dimension
        integer :: nvel_dimid ! ID of nvel dimension
        integer :: max_nvell_dimid ! ID of max_nvell dimension
        integer :: mesh_dimid ! ID of mesh dimension
        
        ! Variable IDs
        integer :: time_varid ! ID of time variable
        integer :: x_varid ! ID of x variable
        integer :: y_varid ! ID of y variable
        integer :: element_varid ! ID of element variable
        integer :: adcirc_mesh_varid ! ID of adcirc_mesh variable
        integer :: neta_var_varid ! ID of neta variable, renamed
        integer :: nvdll_varid ! ID of nvdll variable
        integer :: max_nvdll_varid ! ID of max_nvdll variable
        integer :: ibtypee_varid ! ID of ibtypee variable
        integer :: nbdv_varid ! ID of nbdv variable
        integer :: nvel_var_varid ! ID of nvel variable, renamed
        integer :: nvell_varid ! ID of nvell variable
        integer :: max_nvell_varid ! ID of max_nvell variable
        integer :: ibtype_varid ! ID of ibtype variable
        integer :: nbvv_varid ! ID of nbvv variable
        integer :: depth_varid ! ID of depth variable
        integer :: zeta_varid ! ID of zeta variable
        
        ! Put in define mode
        stat = nf90_redef(self%ncid)
        if (stat .ne. nf90_eindefine) then
            call ncfile_check_error(stat)
        endif
        
        ! Define dimensions
        stat = nf90_def_dim(self%ncid, "time", time_dimsize, time_dimid)
        call ncfile_check_error(stat)
        stat = nf90_def_dim(self%ncid, "node", node_dimsize, node_dimid)
        call ncfile_check_error(stat)
        stat = nf90_def_dim(self%ncid, "nele", nele_dimsize, nele_dimid)
        call ncfile_check_error(stat)
        stat = nf90_def_dim(self%ncid, "nvertex", nvertex_dimsize, &
            nvertex_dimid)
        call ncfile_check_error(stat)
        stat = nf90_def_dim(self%ncid, "nope", nope_dimsize, nope_dimid)
        call ncfile_check_error(stat)
        stat = nf90_def_dim(self%ncid, "neta", neta_dimsize, neta_dimid)
        call ncfile_check_error(stat)
        stat = nf90_def_dim(self%ncid, "max_nvdll", 1, &
            max_nvdll_dimid)
        call ncfile_check_error(stat)
        stat = nf90_def_dim(self%ncid, "nbou", nbou_dimsize, nbou_dimid)
        call ncfile_check_error(stat)
        stat = nf90_def_dim(self%ncid, "nvel", nvel_dimsize, nvel_dimid)
        call ncfile_check_error(stat)
        stat = nf90_def_dim(self%ncid, "max_nvell", max_nvell_dimsize, &
            max_nvell_dimid)
        call ncfile_check_error(stat)
        stat = nf90_def_dim(self%ncid, "mesh", mesh_dimsize, mesh_dimid)
        call ncfile_check_error(stat)
        
        ! Define variables and their attributes
        stat = nf90_def_var(self%ncid, "time", nf90_double, &
            time_dimid, varid=time_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "x", nf90_double, &
            node_dimid, varid=x_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "y", nf90_double, &
            node_dimid, varid=y_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "element", nf90_int, &
            (/ nvertex_dimid, nele_dimid /), varid=element_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "adcirc_mesh", nf90_int, &
            mesh_dimid, varid=adcirc_mesh_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "neta_var", nf90_int, &
            varid=neta_var_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "nvdll", nf90_int, &
            nope_dimid, varid=nvdll_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "max_nvdll", nf90_int, &
            varid=max_nvdll_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "ibtypee", nf90_int, &
            nope_dimid, varid=ibtypee_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "nbdv", nf90_int, &
            neta_dimid, varid=nbdv_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "nvel_var", nf90_int, &
            varid=nvel_var_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "nvell", nf90_int, &
            nbou_dimid, varid=nvell_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "max_nvell", nf90_int, &
            varid=max_nvell_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "ibtype", nf90_int, &
            nbou_dimid, varid=ibtype_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "nbvv", nf90_int, &
            nvel_dimid, varid=nbvv_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "depth", nf90_double, &
            node_dimid, varid=depth_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "zeta", nf90_double, &
            (/ node_dimid, time_dimid /), varid=zeta_varid)
        call ncfile_check_error(stat)
        
        ! Add attributes to each variable
        ! stat = nf90_put_att()
        ! call ncfile_check_error(stat)
    end procedure ncdg_63_set_metadata

    module procedure ncdg_63_write_step
        ! See interface for arguments and documentation

        ! N/A
    end procedure ncdg_63_write_step
    
end submodule ncdg_63_sub