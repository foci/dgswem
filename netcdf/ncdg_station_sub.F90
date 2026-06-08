submodule (ncdg:ncdg_file_sub) ncdg_station_sub
!! Implementation of functions common to all DGSWEM station output files.

    implicit none

    contains

    module procedure ncdg_station_open
        ! See interface for arguments and documentation

        ! integer :: ncstat ! Status of most recent operation
        integer :: the_mode ! For defaulting

        ! Set mode
        if (present(mode)) then
            the_mode = mode
        else
            the_mode = nf90_nowrite
        end if

        ! Call parent function
        call self%ncdg_file%open(mode=the_mode)

        ! Set dimension IDs
        ! N/A, for now

        ! Set variable IDs
        ! N/A, for now
    end procedure ncdg_station_open

#ifdef CMPI
    module procedure ncdg_station_open_parallel
        ! See interface for arguments and documentation

        ! integer :: ncstat ! Status of most recent operation
        integer :: the_mode ! For defaulting
        integer :: the_comm ! For defaulting
        integer :: the_info ! For defaulting

        ! Set mode
        if (present(mode)) then
            the_mode = mode
        else
            the_mode = nf90_nowrite
        end if

        ! Set comm
        if (present(comm)) then
            the_comm = comm
        else
            the_comm = mpi_comm_world
        end if

        ! Set info
        if (present(info)) then
            the_info = info
        else
            the_info = mpi_info_null
        end if

        ! Call parent function
        call self%ncdg_file%popen(mode=the_mode, comm=the_comm, info=the_info)

        ! Set dimension IDs
        ! N/A, for now

        ! Set variable IDs
        ! N/A, for now
    end procedure ncdg_station_open_parallel
#endif

    module procedure ncdg_station_close
        ! See interface for arguments and documentation

        ! integer :: ncstat ! Status of most recent operation

        ! Call parent function
        call self%ncdg_file%close()

        ! Flush dimension IDs
        ! N/A, for now

        ! Flush variable IDs
        ! N/A, for now
    end procedure ncdg_station_close

    module procedure ncdg_station_set_metadata
        ! See interface for arguments and documentation

        ! integer :: ncstat ! Status of most recent operation

        ! Call parent function

        ! Define dimensions and their attributes
        ! N/A, for now

        ! Define variables and their attributes
        ! N/A, for now

        ! Define attributes
        ! N/A, for now
    end procedure ncdg_station_set_metadata

end submodule ncdg_station_sub
