submodule (ncdg:ncdg_file_sub) ncdg_station_sub
!! Implementation of functions common to all DGSWEM station output files.

    implicit none

    contains

    module procedure ncdg_station_open
        ! See interface for arguments and documentation

        ! integer :: ncstat ! Status of most recent operation
        integer :: the_mode ! For defaulting

        ! Set path
        ! Default is whatever is already there
        if (present(path)) then
            self%path = path
        else if (.not. allocated(self%path)) then
            error stop "NCDG station error: Tried to open file with no path."
        end if

        ! Set mode
        if (present(mode)) then
            the_mode = mode
        else
            the_mode = nf90_nowrite
        end if

        ! Call parent function
        call self%ncdg_file%open(path=self%path, mode=the_mode)

        ! Set dimension IDs
        ! N/A, for now

        ! Set variable IDs
        ! N/A, for now
    end procedure ncdg_station_open

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
