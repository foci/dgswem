module ncdg
    !! Module for creating and writing to ADCIRC-style fort files.
    
    use netcdf
    use ncutil, only : ncfile, ncfile_check_error
    
    implicit none

    private

    ! Only these specific child types are exposed to the public. The
    ! parent types are just for organization, exploiting inheritance.
    public :: ncdg_63

    type, extends(ncfile) :: ncdg_file
        !! Generic DGSWEM NetCDF file parent type. This contains data
        !! common to all DGSWEM output files.

        ! Dimension IDs
        integer :: time_varid = -1
        !! ID of time dimension

        ! Variable IDs
        integer :: time_dimid = -1
        !! ID of time variable

        contains

        procedure, public :: open => ncdg_file_open
        procedure, public :: close => ncdg_file_close
        procedure, public :: ncdg_file_set_metadata
        
    end type ncdg_file

    type, extends(ncdg_file) :: ncdg_nodal
        !! Generic DGSWEM nodal file parent type. This contains data
        !! common to all nodal files.

        ! Dimension IDs

        ! Variable IDs

        contains

        procedure, public :: open => ncdg_nodal_open
        procedure, public :: close => ncdg_nodal_close
        procedure, public :: ncdg_nodal_set_metadata
        procedure, public :: ncdg_nodal_write_mesh

    end type ncdg_nodal

    type, extends(ncdg_file) :: ncdg_station
        !! Generic DGSWEM station file parent type. This contains data
        !! common to all station files.

        ! Dimension IDs

        ! Variable IDs

        contains

        procedure, public :: open => ncdg_station_open
        procedure, public :: close => ncdg_station_close
        procedure, public :: ncdg_station_set_metadata
        procedure, public :: ncdg_station_write_stations

    end type ncdg_station
    
    type, extends(ncdg_nodal) :: ncdg_63
        !! fort.63.nc file editing object

        ! Dimension IDs

        ! Variable IDs

        contains

        procedure, public :: open => ncdg_63_open
        procedure, public :: close => ncdg_63_close
        procedure, public :: ncdg_63_set_metadata
        procedure, public :: ncdg_63_write_step
        
    end type ncdg_63

    interface

        ! ncdg_file subroutines
        ! Common to all DGSWEM output files

        module subroutine ncdg_file_open(self, mode)
            !! Open a DGSWEM output file, and set all IDs.

            implicit none

            class(ncdg_file), intent(inout) :: self
            !! The wrapper object of the file
            integer, optional, intent(in) :: mode
            !! NetCDF open mode. Default is nf90_nowrite, i.e. read-only.
        end subroutine ncdg_file_open
        
        module subroutine ncdg_file_close(self)
            !! Close a DGSWEM output file, and flush all IDs.
            
            implicit none
            
            class(ncdg_file), intent(inout) :: self
            !! The wrapper object of the file
        end subroutine ncdg_file_close

        ! ncdg_nodal subroutines
        ! Common to all nodal files

        module subroutine ncdg_nodal_open(self, mode)
            !! Open a nodal output file, and set all IDs.

            implicit none

            class(ncdg_nodal), intent(inout) :: self
            !! The wrapper object of the file
            integer, optional, intent(in) :: mode
            !! NetCDF open mode. Default is nf90_nowrite, i.e. read-only.
        end subroutine ncdg_nodal_open
        
        module subroutine ncdg_nodal_close(self)
            !! Close a nodal output file, and flush all IDs.
            
            implicit none
            
            class(ncdg_nodal), intent(inout) :: self
            !! The wrapper object of the file
        end subroutine ncdg_nodal_close

        ! ncdg_station subroutines
        ! Common to all station files

        module subroutine ncdg_station_open(self, mode)
            !! Open a station output file, and set all IDs.

            implicit none

            class(ncdg_station), intent(inout) :: self
            !! The wrapper object of the file
            integer, optional, intent(in) :: mode
            !! NetCDF open mode. Default is nf90_nowrite, i.e. read-only.
        end subroutine ncdg_station_open
        
        module subroutine ncdg_station_close(self)
            !! Close a station output file, and flush all IDs.
            
            implicit none
            
            class(ncdg_file), intent(inout) :: self
            !! The wrapper object of the file
        end subroutine ncdg_station_close

        ! ncdg_63 subroutines
        ! For fort.63.nc

        module subroutine ncdg_63_open(self, mode)
            !! Open a fort.63.nc output file, and set all IDs.

            implicit none

            class(ncdg_63), intent(inout) :: self
            !! The wrapper object of the file
            integer, optional, intent(in) :: mode
            !! NetCDF open mode. Default is nf90_nowrite, i.e. read-only.
        end subroutine ncdg_63_open
        
        module subroutine ncdg_63_close(self)
            !! Close a fort.63.nc output file, and flush all IDs.
            
            implicit none
            
            class(ncdg_63), intent(inout) :: self
            !! The wrapper object of the file
        end subroutine ncdg_63_close
        
        module subroutine ncdg_63_set_metadata(self, time_dimsize, &
            node_dimsize, nele_dimsize, nvertex_dimsize, &
            nope_dimsize, neta_dimsize, max_nvdll_dimsize, &
            nbou_dimsize, nvel_dimsize, max_nvell_dimsize, &
            mesh_dimsize)
            !! Set variables, dimensions, and metadata for
            !! [[netcdf_fort(module):fort_63_nc(type)]].

            implicit none

            class(ncdg_63), intent(inout) :: self
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
        end subroutine ncdg_63_set_metadata

        module subroutine ncdg_63_write_step(self)
        !! Timestep writing function for
        !! [[netcdf_fort(module):fort_63_nc(type)]].

        implicit none
        
        class(ncdg_63), intent(inout) :: self
        !! The wrapper object being written to
        end subroutine ncdg_63_write_step
    end interface
    
end module ncdg
