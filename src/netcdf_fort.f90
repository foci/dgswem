module netcdf_fort
    !! Module for creating and writing to ADCIRC-style fort files.
    
    use netcdf
    use netcdf_file, only : ncfile, ncfile_check_error
    
    implicit none
    
    type, extends(ncfile) :: fort_63_nc
        !! fort.63.nc file editing object
        
        contains
        
        procedure, public :: set_metadata => fort_63_nc_set_metadata
        !! Set variables, dimensions, and attributes
        procedure, public :: write_step => fort_63_nc_write_step
        !! Write data to the file for a single time step
        
    end type fort_63_nc
    
    contains
    
    subroutine fort_63_nc_set_metadata(self, time_dimsize, node_dimsize, &
        nele_dimsize, nvertex_dimsize, nope_dimsize, neta_dimsize, &
        max_nvdll_dimsize, nbou_dimsize, nvel_dimsize, max_nvell_dimsize, &
        mesh_dimsize)
        !! Set variables, dimensions, and metadata for
        !! [[netcdf_fort(module):fort_63_nc(type)]].
        
        implicit none
        
        class(fort_63_nc), intent(inout) :: self
        !! The wrapper object being initialized
        integer :: time_dimsize
        !! Size of time dimension
        integer :: node_dimsize
        !! Size of node dimension
        integer :: nele_dimsize
        !! Size of nele dimension
        integer :: nvertex_dimsize
        !! Size of nvertex dimension
        integer :: nope_dimsize
        !! Size of nope dimension
        integer :: neta_dimsize 
        !! Size of neta dimension
        integer :: max_nvdll_dimsize
        !! Size of max_nvdll dimension
        integer :: nbou_dimsize
        !! Size of nbou dimension
        integer :: nvel_dimsize
        !! Size of nvel dimension
        integer :: max_nvell_dimsize
        !! Size of max_nvell dimension
        integer :: mesh_dimsize
        !! Size of mesh dimension
        
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
            (/ nele_dimid, nvertex_dimid /), varid=element_varid)
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
            (/ time_dimid, node_dimid /), varid=zeta_varid)
        call ncfile_check_error(stat)
        
        ! Add attributes to each variable
        ! stat = nf90_put_att()
        ! call ncfile_check_error(stat)
    end subroutine fort_63_nc_set_metadata
    
    subroutine fort_63_nc_write_step(self)
        !! Timestep writing function for
        !! [[netcdf_fort(module):fort_63_nc(type)]].
        
        implicit none
        
        class(fort_63_nc), intent(inout) :: self
        !! The wrapper object being initialized
    end subroutine fort_63_nc_write_step
    
end module netcdf_fort
