module ncdg
    !! Module for creating and writing to ADCIRC-style fort files.
    
    use netcdf
    use ncutil, only : ncfile, ncfile_check_error
    
    implicit none

    private

    ! Only these specific child types are exposed to the public. The
    ! parent types are just for organization, exploiting inheritance.
    public :: ncdg_63, ncdg_64

    ! Intra-module constants

    integer, parameter :: integer_fill_value = -99999 ! Consistent with ADCIRC
    real(8), parameter :: real_fill_value = -99999.d0 ! Consistent with ADCIRC

    ! Begin type definitions

    type, extends(ncfile) :: ncdg_file
        !! Generic DGSWEM NetCDF file parent type. This contains data
        !! common to all DGSWEM output files.

        ! Dimension IDs
        integer :: time_dimid = -1
        !! ID of time dimension

        ! Dimension names
        character(len=4) :: time_dimname = "time"
        !! Name of time dimension

        ! Variable IDs
        integer :: time_varid = -1
        !! ID of time variable

        ! Variable names
        character(len=4) :: time_varname = "time"
        !! Name of time variable

        contains

        procedure, public :: open => ncdg_file_open
        procedure, public :: close => ncdg_file_close
        procedure, public :: ncdg_file_set_metadata
        
    end type ncdg_file

    type, extends(ncdg_file) :: ncdg_nodal
        !! Generic DGSWEM nodal file parent type. This contains data
        !! common to all nodal files.

        ! Dimension IDs
        integer :: node_dimid = -1
        !! ID of node dimension
        integer :: nele_dimid = -1
        !! ID of nele dimension
        integer :: nvertex_dimid = -1
        !! ID of nvertex dimension
        integer :: nope_dimid = -1
        !! ID of nope dimension
        integer :: neta_dimid = -1
        !! ID of neta dimension
        integer :: max_nvdll_dimid = -1
        !! ID of max_nvdll dimension
        integer :: nbou_dimid = -1
        !! ID of nbou dimension
        integer :: nvel_dimid = -1
        !! ID of nvel dimension
        integer :: max_nvell_dimid = -1
        !! ID of max_nvell dimension
        integer :: mesh_dimid = -1
        !! ID of mesh dimension

        ! Dimension names
        character(len=4) :: node_dimname = "node"
        !! Name of node dimension
        character(len=4) :: nele_dimname = "nele"
        !! Name of nele dimension
        character(len=7) :: nvertex_dimname = "nvertex"
        !! Name of nvertex dimension
        character(len=4) :: nope_dimname = "nope"
        !! Name of nope dimension
        character(len=4) :: neta_dimname = "neta"
        !! Name of neta dimension
        character(len=9) :: max_nvdll_dimname = "max_nvdll"
        !! Name of max_nvdll dimension
        character(len=4) :: nbou_dimname = "nbou"
        !! Name of nbou dimension
        character(len=4) :: nvel_dimname = "nvel"
        !! Name of nvel dimension
        character(len=9) :: max_nvell_dimname = "max_nvell"
        !! Name of max_nvell dimension
        character(len=4) :: mesh_dimname = "mesh"
        !! Name of mesh dimension

        ! Variable IDs
        integer :: x_varid = -1
        !! ID of x variable
        integer :: y_varid = -1
        !! ID of y variable
        integer :: element_varid = -1
        !! ID of element variable
        integer :: adcirc_mesh_varid = -1
        !! ID of adcirc_mesh variable
        integer :: neta_varid = -1
        !! ID of neta variable, renamed
        integer :: nvdll_varid = -1
        !! ID of nvdll variable
        integer :: max_nvdll_varid = -1
        !! ID of max_nvdll variable
        integer :: ibtypee_varid = -1
        !! ID of ibtypee variable
        integer :: nbdv_varid = -1
        !! ID of nbdv variable
        integer :: nvel_varid = -1
        !! ID of nvel variable, renamed
        integer :: nvell_varid = -1
        !! ID of nvell variable
        integer :: max_nvell_varid = -1
        !! ID of max_nvell variable
        integer :: ibtype_varid = -1
        !! ID of ibtype variable
        integer :: nbvv_varid = -1
        !! ID of nbvv variable
        integer :: depth_varid = -1
        !! ID of depth variable

        ! Variable names
        character(len=1) :: x_varname = "x"
        !! Name of x variable
        character(len=1) :: y_varname = "y"
        !! Name of y variable
        character(len=7) :: element_varname = "element"
        !! Name of element variable
        character(len=11) :: adcirc_mesh_varname = "adcirc_mesh"
        !! Name of adcirc_mesh variable
        character(len=8) :: neta_varname = "neta_var"
        !! Name of neta variable
        character(len=5) :: nvdll_varname = "nvdll"
        !! Name of nvdll variable
        character(len=9) :: max_nvdll_varname = "max_nvdll"
        !! Name of max_nvdll variable
        character(len=8) :: ibtypee_varname = "ibtypee"
        !! Name of ibtypee variable
        character(len=4) :: nbdv_varname = "nbdv"
        !! Name of nbdv variable
        character(len=8) :: nvel_varname = "nvel_var"
        !! Name of nvel variable
        character(len=5) :: nvell_varname = "nvell"
        !! Name of nvell variable
        character(len=9) :: max_nvell_varname = "max_nvell"
        !! Name of max_nvell variable
        character(len=6) :: ibtype_varname = "ibtype"
        !! Name of ibtype variable
        character(len=4) :: nbvv_varname = "nbvv"
        !! Name of nbvv variable
        character(len=5) :: depth_varname = "depth"
        !! Name of depth variable

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

        ! Dimension names

        ! Variable IDs

        ! Variable names

        contains

        procedure, public :: open => ncdg_station_open
        procedure, public :: close => ncdg_station_close
        procedure, public :: ncdg_station_set_metadata
        ! procedure, public :: ncdg_station_write_stations

    end type ncdg_station
    
    type, extends(ncdg_nodal) :: ncdg_63
        !! fort.63.nc file editing object

        ! Dimension IDs

        ! Dimension names

        ! Variable IDs
        integer :: zeta_varid = -1
        !! ID of zeta variable

        ! Variable names
        character(len=4) :: zeta_varname = "zeta"
        !! Name of zeta variable

        contains

        procedure, public :: open => ncdg_63_open
        procedure, public :: close => ncdg_63_close
        procedure, public :: ncdg_63_set_metadata
        procedure, public :: ncdg_63_write_step
        
    end type ncdg_63

    type, extends(ncdg_nodal) :: ncdg_64
        !! fort.64.nc file editing object

        ! Dimension IDs

        ! Dimension names

        ! Variable IDs
        integer :: u_vel_varid = -1
        !! ID of u_vel variable
        integer :: v_vel_varid = -1
        !! ID of v_vel variable

        ! Variable names
        character(len=5) :: u_vel_varname = "u_vel"
        !! Name of u_vel variable
        character(len=5) :: v_vel_varname = "v_vel"
        !! Name of v_vel variable

        contains

        procedure, public :: open => ncdg_64_open
        procedure, public :: close => ncdg_64_close
        procedure, public :: ncdg_64_set_metadata
        procedure, public :: ncdg_64_write_step
        
    end type ncdg_64

    ! Begin interface

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

        module subroutine ncdg_file_set_metadata(self, time_dimsize)
            !! Set variables, dimensions, and metadata for a DGSWEM
            !! output file.

            implicit none

            class(ncdg_file), intent(inout) :: self
            !! The wrapper object of the file
            integer, intent(in) :: time_dimsize
            !! Size of time dimension
        end subroutine ncdg_file_set_metadata

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

        module subroutine ncdg_nodal_set_metadata(self, time_dimsize, &
            node_dimsize, nele_dimsize, nvertex_dimsize, &
            nope_dimsize, neta_dimsize, max_nvdll_dimsize, &
            nbou_dimsize, nvel_dimsize, max_nvell_dimsize, &
            mesh_dimsize, nope, nbou, ics)
            !! Set variables, dimensions, and metadata for a nodal
            !! output file.

            implicit none

            class(ncdg_nodal), intent(inout) :: self
            !! The wrapper object of the file
            integer, intent(in):: time_dimsize
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
            ! BEGIN ADCIRC VARIABLES
            integer, intent(in) :: nope
            !! The number of elevation boundary segments. This corresponds
            !! to [[global(module):nope(variable)]].
            integer, intent(in) :: nbou
            !! The number of normal flow (discharge) boundary
            !! segments. This corresponds to
            !! [[global(module):nbou(variable)]].
            integer, intent(in) :: ics
            !! The number of elevation boundary segments. This corresponds
            !! to [[global(module):ics(variable)]].
        end subroutine ncdg_nodal_set_metadata

        module subroutine ncdg_nodal_write_mesh(self, x, y, dp, &
            nm, nhy, ne, nope, neta, ibtypee, nvdll, &
            nbdv, nbou, nvel, ibtype, nvell, nbvv)
            !! Write the mesh data for a nodal file.

            implicit none

            class(ncdg_nodal), intent(inout) :: self
            !! The wrapper object of the file
            real, intent(in) :: x(:)
            !! The x-coordinate of each node. This corresponds to
            !! [[global(module):x(variable)]], of size
            !! [[global(module):mnp(variable)]].
            real, intent(in) :: y(:)
            !! The y-coordinate of each node. This corresponds to
            !! [[global(module):x(variable)]], of size
            !! [[global(module):mnp(variable)]].
            real, intent(in) :: dp(:)
            !! The bathymetry of each node. This corresponds to
            !! [[global(module):dp(variable)]], of size
            !! [[global(module):mnp(variable)]].
            integer, intent(in) :: nm(:, :)
            !! The node numbers for each element. This corresponds to
            !! [[global(module):nm(variable)]], of size
            !! [[global(module):mne(variable)]] by
            !! [[global(module):nhy(variable)]].
            integer, intent(in) :: nhy
            !! The number of vertices per element. This corresponds to
            !! [[global(module):nhy(variable)]].
            integer, intent(in) :: ne
            !! The number of elements in the mesh. This corresponds to
            !! [[global(module):ne(variable)]].
            integer, intent(in) :: nope
            !! The number of elevation boundary segments. This corresponds
            !! to [[global(module):nope(variable)]].
            integer, intent(in) :: neta
            !! The number of elevation boundary nodes. This corresponds
            !! to [[global(module):neta(variable)]].
            integer, intent(in) :: ibtypee(:)
            !! The elevation boundary type for each boundary segment.
            !! This does not yet correspond to a DGSWEM variable. The
            !! only accepted flag currently is 0. For now, just
            !! substitute a zero array of size
            !! [[global(module):nope(variable)]].
            integer, intent(in) :: nvdll(:)
            !! The number of nodes in a given elevation boundary
            !! segment. This corresponds to
            !! [[global(module):nvdll(variable)]], of size
            !! [[global(module):mnope(variable)]].
            integer, intent(in) :: nbdv(:, :)
            !! The node numbers for each elevation boundary segment.
            !! This corresponds to
            !! [[global(module):nbdv(variable)]], of size
            !! [[global(module):mnope(variable)]] by
            !! [[global(module):mneta(variable)]].
            integer, intent(in) :: nbou
            !! The number of normal flow (discharge) boundary
            !! segments. This corresponds to
            !! [[global(module):nbou(variable)]].
            integer, intent(in) :: nvel
            !! The number of normal flow (discharge) boundary nodes.
            !! This corresponds to
            !! [[global(module):nvel(variable)]].
            integer, intent(in) :: ibtype(:)
            !! The normal flow (discharge) boundary type for each
            !! boundary segment. This corresponds to
            !! [[global(module):segtype(variable)]] of size
            !! [[global(module):mnbou(variable)]].
            integer, intent(in) :: nvell(:)
            !! The number of nodes in a given normal flow (discharge)
            !! boundary segment. This corresponds to
            !! [[global(module):nvell(variable)]], of size
            !! [[global(module):mnbou(variable)]].
            integer, intent(in) :: nbvv(:, :)
            !! The node numbers for each normal flow (discharge)
            !! boundary segment. This corresponds to
            !! [[global(module):nbvv(variable)]], of size
            !! [[global(module):mnbou(variable)]] by
            !! 0:[[global(module):mnvel(variable)]].
        end subroutine ncdg_nodal_write_mesh

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
            
            class(ncdg_station), intent(inout) :: self
            !! The wrapper object of the file
        end subroutine ncdg_station_close

        module subroutine ncdg_station_set_metadata(self)
            !! Set variables, dimensions, and metadata for a nodal
            !! output file.

            implicit none

            class(ncdg_station), intent(inout) :: self
            !! The wrapper object of the file
        end subroutine ncdg_station_set_metadata

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
            mesh_dimsize, nope, nbou, ics)
            !! Set variables, dimensions, and metadata for a fort.63.nc
            !! output file.

            implicit none

            class(ncdg_63), intent(inout) :: self
            !! The wrapper object being initialized
            integer, intent(in):: time_dimsize
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
            ! BEGIN ADCIRC VARIABLES
            integer, intent(in) :: nope
            !! The number of elevation boundary segments. This corresponds
            !! to [[global(module):nope(variable)]].
            integer, intent(in) :: nbou
            !! The number of normal flow (discharge) boundary
            !! segments. This corresponds to
            !! [[global(module):nbou(variable)]].
            integer, intent(in) :: ics
            !! The number of elevation boundary segments. This corresponds
            !! to [[global(module):ics(variable)]].
        end subroutine ncdg_63_set_metadata

        module subroutine ncdg_63_write_step(self, t, zeta)
            !! Timestep writing function for a fort.63.nc output file.
    
            implicit none
            
            class(ncdg_63), intent(inout) :: self
            !! The wrapper object being written to
            real, intent(in) :: t
            !! The current time
            real, intent(in) :: zeta(:)
            !! The current zeta values
        end subroutine ncdg_63_write_step

        ! ncdg_64 subroutines
        ! For fort.64.nc

        module subroutine ncdg_64_open(self, mode)
            !! Open a fort.64.nc output file, and set all IDs.

            implicit none

            class(ncdg_64), intent(inout) :: self
            !! The wrapper object of the file
            integer, optional, intent(in) :: mode
            !! NetCDF open mode. Default is nf90_nowrite, i.e. read-only.
        end subroutine ncdg_64_open
        
        module subroutine ncdg_64_close(self)
            !! Close a fort.64.nc output file, and flush all IDs.
            
            implicit none
            
            class(ncdg_64), intent(inout) :: self
            !! The wrapper object of the file
        end subroutine ncdg_64_close

        module subroutine ncdg_64_set_metadata(self, time_dimsize, &
            node_dimsize, nele_dimsize, nvertex_dimsize, &
            nope_dimsize, neta_dimsize, max_nvdll_dimsize, &
            nbou_dimsize, nvel_dimsize, max_nvell_dimsize, &
            mesh_dimsize, nope, nbou, ics)
            !! Set variables, dimensions, and metadata for a fort.64.nc
            !! output file.

            implicit none

            class(ncdg_64), intent(inout) :: self
            !! The wrapper object being initialized
            integer, intent(in):: time_dimsize
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
            ! BEGIN ADCIRC VARIABLES
            integer, intent(in) :: nope
            !! The number of elevation boundary segments. This corresponds
            !! to [[global(module):nope(variable)]].
            integer, intent(in) :: nbou
            !! The number of normal flow (discharge) boundary
            !! segments. This corresponds to
            !! [[global(module):nbou(variable)]].
            integer, intent(in) :: ics
            !! The number of elevation boundary segments. This corresponds
            !! to [[global(module):ics(variable)]].
        end subroutine ncdg_64_set_metadata

        module subroutine ncdg_64_write_step(self, t, u_vel, v_vel)
            !! Timestep writing function for a fort.64.nc output file.
    
            implicit none
            
            class(ncdg_64), intent(inout) :: self
            !! The wrapper object being written to
            real, intent(in) :: t
            !! The current time
            real, intent(in) :: u_vel(:)
            !! The current u_vel values
            real, intent(in) :: v_vel(:)
            !! The current v_vel values
        end subroutine ncdg_64_write_step
        
    end interface
    
end module ncdg
