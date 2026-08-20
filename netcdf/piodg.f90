module piodg
    !! Module for creating and writing ADCIRC-style fort files in parallel.
    !!
    !! The derived type inheritance hierarchy is as follows.
    !! - piofile
    !!     - piodg_file
    !!         - piog_nodal
    !!             - piodg_nodal_scalar
    !!                 - piodg_63
    !!                 - piodg_73
    !!             - piodg_nodal_scalar_max
    !!                 - piodg_maxele
    !!             - piodg_nodal_vector
    !!                 - piodg_64
    !!                 - piodg_74
    !!         - piodg_station

    use netcdf, only : nf90_enotindefine
    use pio
    use pioutil, only: piofile, piofile_check_error

    implicit none

    private

    ! Only these specific child types are exposed to the public. The
    ! parent types are just for organization, exploiting inheritance.
    public :: piodg_63, piodg_64, piodg_73, piodg_74, piodg_maxele

    ! Begin type definitions

    type, extends(piofile) :: piodg_file
        !! Generic DGSWEM PIO parent type. This contains data common to all
        !! DGSWEM output files.

        ! Dimension IDs
        integer :: time_dimid = -1
        !! ID of time dimension

        ! Dimension names
        character(len=4) :: time_dimname = "time"
        !! Name of time dimension

        ! Variable IDs
        type(var_desc_t) :: time_vardesc
        !! ID struct of time variable

        ! Variable names
        character(len=4) :: time_varname = "time"
        !! Name of time variable

        contains

        procedure, public :: open => piodg_file_open
        procedure, public :: close => piodg_file_close

    end type piodg_file

    type, extends(piodg_file) :: piodg_nodal
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
        type(var_desc_t) :: x_vardesc
        !! ID struct of x variable
        type(var_desc_t) :: y_vardesc
        !! ID struct of y variable
        type(var_desc_t) :: element_vardesc
        !! ID struct of element variable
        type(var_desc_t) :: adcirc_mesh_vardesc
        !! ID struct of adcirc_mesh variable
        type(var_desc_t) :: neta_vardesc
        !! ID struct of neta variable, renamed
        type(var_desc_t) :: nvdll_vardesc
        !! ID struct of nvdll variable
        type(var_desc_t) :: max_nvdll_vardesc
        !! ID struct of max_nvdll variable
        type(var_desc_t) :: ibtypee_vardesc
        !! ID struct of ibtypee variable
        type(var_desc_t) :: nbdv_vardesc
        !! ID struct of nbdv variable
        type(var_desc_t) :: nvel_vardesc
        !! ID struct of nvel variable, renamed
        type(var_desc_t) :: nvell_vardesc
        !! ID struct of nvell variable
        type(var_desc_t) :: max_nvell_vardesc
        !! ID struct of max_nvell variable
        type(var_desc_t) :: ibtype_vardesc
        !! ID struct of ibtype variable
        type(var_desc_t) :: nbvv_vardesc
        !! ID struct of nbvv variable
        type(var_desc_t) :: depth_vardesc
        !! ID struct of depth variable

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

        procedure, public :: open => piodg_nodal_open
        procedure, public :: close => piodg_nodal_close

    end type piodg_nodal

    type, extends(piodg_nodal) :: piodg_nodal_scalar
        !! Generic nodal scalar file editing object

        ! Dimension IDs

        ! Dimension names

        ! Variable IDs
        type(var_desc_t) :: scalar_vardesc
        !! ID of scalar variable

        ! Variable names
        character(len=:), allocatable :: scalar_varname
        !! Name of scalar variable

        contains

        procedure, public :: open => piodg_nodal_scalar_open
        procedure, public :: close => piodg_nodal_scalar_close
        procedure, public :: write_step => piodg_nodal_scalar_write_step

    end type piodg_nodal_scalar

    type, extends(piodg_nodal_scalar) :: piodg_63
        !! fort.63.nc file editing object

        contains

        procedure, public :: open => piodg_63_open

    end type piodg_63

    type, extends(piodg_nodal_scalar) :: piodg_73
        !! fort.73.nc file editing object

        contains

        procedure, public :: open => piodg_73_open

    end type piodg_73

    type, extends(piodg_nodal) :: piodg_nodal_scalar_max
        !! Generic nodal scalar max file editing object

        ! Dimension IDs

        ! Dimension names

        ! Variable IDs
        type(var_desc_t) :: scalar_max_vardesc
        !! ID of scalar_max variable
        type(var_desc_t) :: time_of_scalar_max_vardesc
        !! ID of time_of_scalar_max variable

        ! Variable names
        character(len=:), allocatable :: scalar_max_varname
        !! Name of scalar_max variable
        character(len=:), allocatable :: time_of_scalar_max_varname
        !! Name of time_of_scalar_max variable

        contains

        procedure, public :: open => piodg_nodal_scalar_max_open
        procedure, public :: close => piodg_nodal_scalar_max_close
        procedure, public :: write_step => piodg_nodal_scalar_max_write_step

    end type piodg_nodal_scalar_max

    type, extends(piodg_nodal_scalar_max) :: piodg_maxele
        !! maxele.63.nc file editing object

        contains

        procedure, public :: open => piodg_maxele_open

    end type piodg_maxele

    type, extends(piodg_nodal) :: piodg_nodal_vector
        !! Generic nodal vector file editing object

        ! Dimension IDs

        ! Dimension names

        ! Variable IDs
        type(var_desc_t) :: vector_u_vardesc
        !! ID of vector_u variable
        type(var_desc_t) :: vector_v_vardesc
        !! ID of vector_v variable

        ! Variable names
        character(len=:), allocatable :: vector_u_varname
        !! Name of vector_u variable
        character(len=:), allocatable :: vector_v_varname
        !! Name of vector_v variable

        contains

        procedure, public :: open => piodg_nodal_vector_open
        procedure, public :: close => piodg_nodal_vector_close
        procedure, public :: write_step => piodg_nodal_vector_write_step

    end type piodg_nodal_vector

    type, extends(piodg_nodal_vector) :: piodg_64
        !! fort.64.nc file editing object

        contains

        procedure, public :: open => piodg_64_open

    end type piodg_64

    type, extends(piodg_nodal_vector) :: piodg_74
        !! fort.74.nc file editing object

        contains

        procedure, public :: open => piodg_74_open

    end type piodg_74

    ! type, extends(piodg_file) :: piodg_station
    !     !! Generic DGSWEM station file parent type. This contains data
    !     !! common to all station files.

    !     ! Dimension IDs

    !     ! Dimension names

    !     ! Variable IDs

    !     ! Variable names

    !     contains

    !     procedure, public :: open => piodg_station_open
    !     procedure, public :: close => piodg_station_close

    ! end type piodg_station

    ! Begin interface

    interface

        ! piodg_file subroutines
        ! Common to all DGSWEM output files

        module subroutine piodg_file_open(self, piosystem, path, piotype, omode)
            !! Open a DGSWEM output file, and set all IDs.

            implicit none

            class(piodg_file), intent(inout) :: self
            !! The wrapper object of the file
            type(iosystem_desc_t), intent(inout), target :: piosystem
            !! The PIO system object
            character(len=*), optional, intent(in) :: path
            !! Path and name for the NetCDF file. Default is whatever path is
            !! already set. If the path is not already set and nothing is provided,
            !! the program will stop with an error.
            integer, optional, intent(in) :: piotype
            !! PIO file and I/O type. Default is pio_iotype_netcdf4p.
            integer, optional, intent(in) :: omode
            !! PIO open mode. Default is pio_nowrite.
        end subroutine piodg_file_open

        module subroutine piodg_file_close(self)
            !! Close a DGSWEM output file, and flush all IDs.

            implicit none

            class(piodg_file), intent(inout) :: self
            !! The wrapper object of the file
        end subroutine piodg_file_close

        ! piodg_nodal subroutines
        ! Common to all nodal files

        module subroutine piodg_nodal_open(self, piosystem, path, piotype, &
            omode)
            !! Open a nodal output file, and set all IDs.

            implicit none

            class(piodg_nodal), intent(inout) :: self
            !! The wrapper object of the file
            type(iosystem_desc_t), intent(inout), target :: piosystem
            !! The PIO system object
            character(len=*), optional, intent(in) :: path
            !! Path and name for the NetCDF file. Default is whatever path is
            !! already set. If the path is not already set and nothing is provided,
            !! the program will stop with an error.
            integer, optional, intent(in) :: piotype
            !! PIO file and I/O type. Default is pio_iotype_netcdf4p.
            integer, optional, intent(in) :: omode
            !! PIO open mode. Default is pio_nowrite.
        end subroutine piodg_nodal_open

        module subroutine piodg_nodal_close(self)
            !! Close a nodal output file, and flush all IDs.

            implicit none

            class(piodg_nodal), intent(inout) :: self
            !! The wrapper object of the file
        end subroutine piodg_nodal_close

        ! piodg_nodal_scalar subroutines
        ! Common to all nodal scalar files
        ! open, close, write step

        module subroutine piodg_nodal_scalar_open(self, piosystem, path, &
            piotype, omode)
            !! Open a nodal scalar output file, and set all IDs

            implicit none

            class(piodg_nodal_scalar), intent(inout) :: self
            !! The wrapper object of the file
            type(iosystem_desc_t), intent(inout), target :: piosystem
            !! The PIO system object
            character(len=*), optional, intent(in) :: path
            !! Path and name for the NetCDF file. Default is whatever path is
            !! already set. If the path is not already set and nothing is provided,
            !! the program will assume nodalscalar.nc.
            integer, optional, intent(in) :: piotype
            !! PIO file and I/O type. Default is pio_iotype_netcdf4p.
            integer, optional, intent(in) :: omode
            !! PIO open mode. Default is pio_nowrite.
        end subroutine piodg_nodal_scalar_open

        module subroutine piodg_nodal_scalar_close(self)
            !! Close a nodal scalar output file, and flush all IDs.

            implicit none

            class(piodg_nodal_scalar), intent(inout) :: self
            !! The wrapper object of the file
        end subroutine piodg_nodal_scalar_close

        module subroutine piodg_nodal_scalar_write_step(self, t, scalar, &
            piodesc, sync)
            !! Parallel timestep writing function for a nodal scalar output
            !! file.

            implicit none

            class(piodg_nodal_scalar), intent(inout) :: self
            !! The wrapper object being written to
            real, intent(in) :: t
            !! The current time
            real, intent(in) :: scalar(:)
            !! The current zeta values
            type(io_desc_t), intent(inout) :: piodesc
            !! An I/O parameter struct associated with a decomposition
            logical, optional, intent(in) :: sync
            !! Whether to write to disk immediately. Default is .false.
        end subroutine piodg_nodal_scalar_write_step

        ! piodg_63 subroutines
        ! For fort.63.nc

        module subroutine piodg_63_open(self, piosystem, path, piotype, omode)
            !! Open a fort.63.nc output file, and set all IDs

            implicit none

            class(piodg_63), intent(inout) :: self
            !! The wrapper object of the file
            type(iosystem_desc_t), intent(inout), target :: piosystem
            !! The PIO system object
            character(len=*), optional, intent(in) :: path
            !! Path and name for the NetCDF file. Default is whatever path is
            !! already set. If the path is not already set and nothing is provided,
            !! the program will assume fort.63.nc.
            integer, optional, intent(in) :: piotype
            !! PIO file and I/O type. Default is pio_iotype_netcdf4p.
            integer, optional, intent(in) :: omode
            !! PIO open mode. Default is pio_nowrite.
        end subroutine piodg_63_open

        ! piodg_73 subroutines
        ! For fort.73.nc

        module subroutine piodg_73_open(self, piosystem, path, piotype, omode)
            !! Open a fort.73.nc output file, and set all IDs

            implicit none

            class(piodg_73), intent(inout) :: self
            !! The wrapper object of the file
            type(iosystem_desc_t), intent(inout), target :: piosystem
            !! The PIO system object
            character(len=*), optional, intent(in) :: path
            !! Path and name for the NetCDF file. Default is whatever path is
            !! already set. If the path is not already set and nothing is provided,
            !! the program will assume fort.73.nc.
            integer, optional, intent(in) :: piotype
            !! PIO file and I/O type. Default is pio_iotype_netcdf4p.
            integer, optional, intent(in) :: omode
            !! PIO open mode. Default is pio_nowrite.
        end subroutine piodg_73_open

        ! piodg_nodal_scalar_max subroutines
        ! Common to all nodal scalar max files

        module subroutine piodg_nodal_scalar_max_open(self, piosystem, path, &
            piotype, omode)
            !! Open a nodal scalar max output file, and set all IDs

            implicit none

            class(piodg_nodal_scalar_max), intent(inout) :: self
            !! The wrapper object of the file
            type(iosystem_desc_t), intent(inout), target :: piosystem
            !! The PIO system object
            character(len=*), optional, intent(in) :: path
            !! Path and name for the NetCDF file. Default is whatever path is
            !! already set. If the path is not already set and nothing is provided,
            !! the program will assume fort.63.nc.
            integer, optional, intent(in) :: piotype
            !! PIO file and I/O type. Default is pio_iotype_netcdf4p.
            integer, optional, intent(in) :: omode
            !! PIO open mode. Default is pio_nowrite.
        end subroutine piodg_nodal_scalar_max_open

        module subroutine piodg_nodal_scalar_max_close(self)
            !! Close a nodal scalar max output file, and flush all IDs.

            implicit none

            class(piodg_nodal_scalar_max), intent(inout) :: self
            !! The wrapper object of the file
        end subroutine piodg_nodal_scalar_max_close

        module subroutine piodg_nodal_scalar_max_write_step(self, t, scalar, &
            np, piodesc, compute_max, sync)
            !! Timestep writing function for a nodal scalar max output file.

            implicit none

            class(piodg_nodal_scalar_max), intent(inout) :: self
            !! The wrapper object being written to
            real, intent(in) :: t
            !! The current time
            real, intent(in) :: scalar(:)
            !! The current or current maximum scalar values
            integer, intent(in) :: np
            !! The total number of nodes on the current process
            type(io_desc_t), intent(inout) :: piodesc
            !! An I/O parameter struct associated with a decomposition
            logical, optional, intent(in) :: compute_max
            !! Whether to compute a new maximum, or overwrite everything with
            !! the value in scalar. Default is .false.
            logical, optional, intent(in) :: sync
            !! Whether to write to disk immediately. Default is .false.
        end subroutine piodg_nodal_scalar_max_write_step

        ! piodg_maxele subroutines
        ! For maxele.63.nc

        module subroutine piodg_maxele_open(self, piosystem, path, piotype, omode)
            !! Open a maxele.63.nc output file, and set all IDs

            implicit none

            class(piodg_maxele), intent(inout) :: self
            !! The wrapper object of the file
            type(iosystem_desc_t), intent(inout), target :: piosystem
            !! The PIO system object
            character(len=*), optional, intent(in) :: path
            !! Path and name for the NetCDF file. Default is whatever path is
            !! already set. If the path is not already set and nothing is provided,
            !! the program will assume fort.63.nc.
            integer, optional, intent(in) :: piotype
            !! PIO file and I/O type. Default is pio_iotype_netcdf4p.
            integer, optional, intent(in) :: omode
            !! PIO open mode. Default is pio_nowrite.
        end subroutine piodg_maxele_open

        ! piodg_nodal_vector subroutines
        ! Common to all nodal vector files

        module subroutine piodg_nodal_vector_open(self, piosystem, path, &
            piotype, omode)
            !! Open a nodal vector output file, and set all IDs

            implicit none

            class(piodg_nodal_vector), intent(inout) :: self
            !! The wrapper object of the file
            type(iosystem_desc_t), intent(inout), target :: piosystem
            !! The PIO system object
            character(len=*), optional, intent(in) :: path
            !! Path and name for the NetCDF file. Default is whatever path is
            !! already set. If the path is not already set and nothing is provided,
            !! the program will assume nodalvector.nc.
            integer, optional, intent(in) :: piotype
            !! PIO file and I/O type. Default is pio_iotype_netcdf4p.
            integer, optional, intent(in) :: omode
            !! PIO open mode. Default is pio_nowrite.
        end subroutine piodg_nodal_vector_open

        module subroutine piodg_nodal_vector_close(self)
            !! Close a nodal vector output file, and flush all IDs.

            implicit none

            class(piodg_nodal_vector), intent(inout) :: self
            !! The wrapper object of the file
        end subroutine piodg_nodal_vector_close

        module subroutine piodg_nodal_vector_write_step(self, t, vector_u, &
            vector_v, piodesc, sync)
            !! Parallel timestep writing function for a nodal vector output
            !! file.

            implicit none

            class(piodg_nodal_vector), intent(inout) :: self
            !! The wrapper object being written to
            real, intent(in) :: t
            !! The current time
            real, intent(in) :: vector_u(:)
            !! The current vector_u values
            real, intent(in) :: vector_v(:)
            !! The current vector_v values
            type(io_desc_t), intent(inout) :: piodesc
            !! An I/O parameter struct associated with a decomposition
            logical, optional, intent(in) :: sync
            !! Whether to write to disk immediately. Default is .false.
        end subroutine piodg_nodal_vector_write_step

        ! piodg_64 subroutines
        ! For fort.64.nc

        module subroutine piodg_64_open(self, piosystem, path, piotype, omode)
            !! Open a fort.64.nc output file, and set all IDs

            implicit none

            class(piodg_64), intent(inout) :: self
            !! The wrapper object of the file
            type(iosystem_desc_t), intent(inout), target :: piosystem
            !! The PIO system object
            character(len=*), optional, intent(in) :: path
            !! Path and name for the NetCDF file. Default is whatever path is
            !! already set. If the path is not already set and nothing is provided,
            !! the program will assume fort.64.nc.
            integer, optional, intent(in) :: piotype
            !! PIO file and I/O type. Default is pio_iotype_netcdf4p.
            integer, optional, intent(in) :: omode
            !! PIO open mode. Default is pio_nowrite.
        end subroutine piodg_64_open

        ! piodg_74 subroutines
        ! For fort.74.nc

        module subroutine piodg_74_open(self, piosystem, path, piotype, omode)
            !! Open a fort.74.nc output file, and set all IDs

            implicit none

            class(piodg_74), intent(inout) :: self
            !! The wrapper object of the file
            type(iosystem_desc_t), intent(inout), target :: piosystem
            !! The PIO system object
            character(len=*), optional, intent(in) :: path
            !! Path and name for the NetCDF file. Default is whatever path is
            !! already set. If the path is not already set and nothing is provided,
            !! the program will assume fort.74.nc.
            integer, optional, intent(in) :: piotype
            !! PIO file and I/O type. Default is pio_iotype_netcdf4p.
            integer, optional, intent(in) :: omode
            !! PIO open mode. Default is pio_nowrite.
        end subroutine piodg_74_open

        ! piodg_station subroutines
        ! Common to all station files

    end interface

end module piodg
