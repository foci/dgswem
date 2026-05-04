submodule (ncdg) ncdg_nodal_sub
!! Implementation of functions common to all DGSWEM nodal output files.

    implicit none

    contains

    module procedure ncdg_nodal_open
        ! See interface for arguments and documentation
        
        integer :: ncstat ! Status of most recent operation

        ! Call parent function
        if (present(mode)) then
            call self%ncdg_file%open(mode)
        else
            call self%ncdg_file%open()
        endif

        ! Set dimension IDs
        ! N/A, for now

        ! Set variable IDs
        ! N/A, for now
    end procedure ncdg_nodal_open

    module procedure ncdg_nodal_close
        ! See interface for arguments and documentation

        ! Call parent function
        call self%ncdg_file%close()

        ! Flush dimension IDs
        ! N/A, for now

        ! Flush variable IDs
        ! N/A, for now
    end procedure ncdg_nodal_close

    module procedure ncdg_nodal_set_metadata
        ! See interface for arguments and documentation
        
        integer :: stat ! Status of most recent operation

        ! Call parent function

        ! Define dimensions and their attributes
        ! N/A, for now

        ! Define variables and their attributes
        ! N/A, for now

        ! Define attributes
        ! N/A, for now
    end procedure ncdg_nodal_set_metadata

end submodule ncdg_nodal_sub