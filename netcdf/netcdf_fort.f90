module netcdf_fort
    !! Module for creating and writing to ADCIRC-style fort files.
    
    use netcdf
    use netcdf_file, only : ncfile, ncfile_check_error
    
    implicit none

    private

    public :: fort_63_nc
    
    type, extends(ncfile) :: fort_63_nc
        !! fort.63.nc file editing object
        
        contains
        
        procedure, public :: set_metadata => fort_63_nc_set_metadata
        procedure, public :: write_step => fort_63_nc_write_step
        
    end type fort_63_nc

    interface
        module subroutine fort_63_nc_set_metadata(self, time_dimsize, &
            node_dimsize, nele_dimsize, nvertex_dimsize, &
            nope_dimsize, neta_dimsize, max_nvdll_dimsize, &
            nbou_dimsize, nvel_dimsize, max_nvell_dimsize, &
            mesh_dimsize)
            !! Set variables, dimensions, and metadata for
            !! [[netcdf_fort(module):fort_63_nc(type)]].

            class(fort_63_nc), intent(inout) :: self
            !! The wrapper object being initialized
            integer, intent(in) :: time_dimsize
            !! Size of time dimension
            integer, intent(in):: node_dimsize
            !! Size of node dimension
            integer, intent(in) :: nele_dimsize
            !! Size of nele dimension
            integer, intent(in) :: nvertex_dimsize
            !! Size of nvertex dimension
            integer, intent(in) :: nope_dimsize
            !! Size of nope dimension
            integer, intent(in) :: neta_dimsize 
            !! Size of neta dimension
            integer, intent(in) :: max_nvdll_dimsize
            !! Size of max_nvdll dimension
            integer, intent(in) :: nbou_dimsize
            !! Size of nbou dimension
            integer, intent(in) :: nvel_dimsize
            !! Size of nvel dimension
            integer, intent(in) :: max_nvell_dimsize
            !! Size of max_nvell dimension
            integer, intent(in) :: mesh_dimsize
            !! Size of mesh dimension
        end subroutine fort_63_nc_set_metadata

        module subroutine fort_63_nc_write_step(self)
        !! Timestep writing function for
        !! [[netcdf_fort(module):fort_63_nc(type)]].
        
        class(fort_63_nc), intent(inout) :: self
        !! The wrapper object being written to
        end subroutine fort_63_nc_write_step
    end interface
    
end module netcdf_fort
