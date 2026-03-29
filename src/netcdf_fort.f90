module netcdf_fort
    !! Module for creating and writing to ADCIRC-style fort files.
    
    use netcdf
    use netcdf_file
    
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
        !! [[netcdf_file(module):netcdf_file(type)]].
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
        integer :: time_varid ! ID of time variable
        integer :: x_varid ! ID of x variable
        integer :: y_varid ! ID of y variable
        integer :: element_varid ! ID of element variable
        integer :: adcirc_mesh_varid ! ID of adcirc_mesh variable
        
        ! Variable IDs
        integer :: _varid ! ID of  variable
        integer :: _varid ! ID of  variable
        integer :: _varid ! ID of  variable
        integer :: _varid ! ID of  variable
        integer :: _varid ! ID of  variable
        integer :: _varid ! ID of  variable
        integer :: _varid ! ID of  variable
        integer :: _varid ! ID of  variable
        integer :: _varid ! ID of  variable
        integer :: _varid ! ID of  variable
        integer :: _varid ! ID of  variable
        integer :: _varid ! ID of  variable
        
        ! Call the parent initialization function
        call self%ncfile%init(path)
        
        ! Put in define mode
        !!!!!!!!!!!!!!!!
        
        ! Define dimensions
        stat = nf90_def_dim(self%ncid, "time", ???, dimid)
        call ncfile_check_error(stat)
        stat = nf90_def_dim(self%ncid, "node", ???, dimid)
        call ncfile_check_error(stat)
        stat = nf90_def_dim(self%ncid, "nele", ???, dimid)
        call ncfile_check_error(stat)
        stat = nf90_def_dim(self%ncid, "nvertex", ???, dimid)
        call ncfile_check_error(stat)
        stat = nf90_def_dim(self%ncid, "nope", ???, dimid)
        call ncfile_check_error(stat)
        stat = nf90_def_dim(self%ncid, "neta", ???, dimid)
        call ncfile_check_error(stat)
        stat = nf90_def_dim(self%ncid, "max_nvdll", ???, dimid)
        call ncfile_check_error(stat)
        stat = nf90_def_dim(self%ncid, "nbou", ???, dimid)
        call ncfile_check_error(stat)
        stat = nf90_def_dim(self%ncid, "nvel", ???, dimid)
        call ncfile_check_error(stat)
        stat = nf90_def_dim(self%ncid, "max_nvell", ???, dimid)
        call ncfile_check_error(stat)
        stat = nf90_def_dim(self%ncid, "mesh", ???, dimid)
        call ncfile_check_error(stat)
        
        ! Define variables
        stat = nf90_def_var(fort_63_ncid, "time", nf90_double, &
            (/  /), varid)
        call ncfile_check_error(stat)
        stat = nf90_put_att
        call ncfile_check_error(stat)
        
        ! Add attributes to each variable
        
        ! Close file
    end subroutine init
    
end module netcdf_fort