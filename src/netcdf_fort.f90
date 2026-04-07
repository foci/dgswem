module netcdf_fort
    !! Module for creating and writing to ADCIRC-style fort files.
    
    use netcdf
    use netcdf_file, only : ncfile, ncfile_check_error
    
    implicit none
    
    type, extends(ncfile) fort_63_nc
        !! fort.63.nc
        
        contains
            
        procedure, public :: init => fort_63_nc_init
        !! Create a new file
        procedure, public :: write_step => fort_63_nc_write_step
        !! Write data to the file for a single time step
        
    end type
    
    contains
        
    subroutine fort_63_nc_init(self, path)
        !! Initialization function for
        !! [[netcdf_fort(module):fort_63_nc(type)]].
        !!
        !! @warning "File Not Closed"
        !! After initialization, the underlying file object is in define mode.
        !! The user must remember to close it.
        !! @endwarning
        
        implicit none
        
        class(ncfile), intent(inout) :: self
        !! The wrapper object being initialized
        character(len=*), intent(in) :: path
        !! Path and name for the NetCDF file
        
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
        
        ! Call the parent initialization function
        call self%ncfile%init(path)
        
        ! Starts in define mode
        
        ! Define dimensions
        stat = nf90_def_dim(self%ncid, "time", nf90_unlimited, time_dimid)
        call ncfile_check_error(stat)
        stat = nf90_def_dim(self%ncid, "node", nf90_unlimited, node_dimid)
        call ncfile_check_error(stat)
        stat = nf90_def_dim(self%ncid, "nele", nf90_unlimited, nele_dimid)
        call ncfile_check_error(stat)
        stat = nf90_def_dim(self%ncid, "nvertex", nf90_unlimited, nvertex_dimid)
        call ncfile_check_error(stat)
        stat = nf90_def_dim(self%ncid, "nope", nf90_unlimited, nope_dimid)
        call ncfile_check_error(stat)
        stat = nf90_def_dim(self%ncid, "neta", nf90_unlimited, neta_dimid)
        call ncfile_check_error(stat)
        stat = nf90_def_dim(self%ncid, "max_nvdll", nf90_unlimited, &
            max_nvdll_dimid)
        call ncfile_check_error(stat)
        stat = nf90_def_dim(self%ncid, "nbou", nf90_unlimited, nbou_dimid)
        call ncfile_check_error(stat)
        stat = nf90_def_dim(self%ncid, "nvel", nf90_unlimited, nvel_dimid)
        call ncfile_check_error(stat)
        stat = nf90_def_dim(self%ncid, "max_nvell", nf90_unlimited, &
            max_nvell_dimid)
        call ncfile_check_error(stat)
        stat = nf90_def_dim(self%ncid, "mesh", nf90_unlimited, mesh_dimid)
        call ncfile_check_error(stat)
        
        ! Define variables and their attributes
        stat = nf90_def_var(self%ncid, "time", nf90_double, &
            time, time_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "x", nf90_double, &
            node, x_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "y", nf90_double, &
            node, y_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "element", nf90_int, &
            (/ nele_dimid, nvertex_dimid /), element_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "adcirc_mesh", nf90_int, &
            mesh_dimid, adcirc_mesh_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "neta_var", nf90_int, &
            (/ /), neta_var_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "nvdll", nf90_int, &
            nope_dimid, nvdll_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "max_nvdll", nf90_int, &
            (/ /), max_nvdll_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "ibtypee", nf90_int, &
            nope_dimid, ibtypee_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "nbdv", nf90_int, &
            neta_dimid, nbdv_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "nvel_var", nf90_int, &
            (/ /), nvel_var_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "nvell", nf90_int, &
            nbou_dimid, nvell_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "max_nvell", nf90_int, &
            (/ /), max_nvell_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "ibtype", nf90_int, &
            nbou_dimid, ibtype_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "nbvv", nf90_int, &
            nvel_dimid, nbvv_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "depth", nf90_double, &
            node_dimid, depth_varid)
        call ncfile_check_error(stat)
        
        stat = nf90_def_var(self%ncid, "zeta", nf90_double, &
            (/ time_dimid, node_dimid /), zeta_varid)
        call ncfile_check_error(stat)
        
        ! Add attributes to each variable
        ! stat = nf90_put_att()
        ! call ncfile_check_error(stat)
        
        ! Leave define mode?
        ! stat = nf90_enddef(self%ncid)
        ! call ncfile_check_error(stat)
    end subroutine init
    
    subroutine fort_63_nc_write_step(self)
        !! Timestep writing function for
        !! [[netcdf_fort(module):fort_63_nc(type)]].
        
        implicit none
        
        class(ncfile), intent(inout) :: self
        !! The wrapper object being initialized
    end subroutine fort_63_nc_write_step
    
end module netcdf_fort