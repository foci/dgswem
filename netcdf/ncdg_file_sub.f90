submodule (ncdg) ncdg_common_sub
!! Implementation of functions common to all DGSWEM output files.

    implicit none

    contains

    module procedure ncdg_file_open
        ! See interface for arguments and documentation
        
        integer :: ncstat ! Status of most recent operation

        ! Call parent function
        if (present(mode)) then
            call self%ncfile%open(mode)
        else
            call self%ncfile%open()
        endif

        ! Set dimension IDs
        ncstat = nf90_inq_dimid(self%ncid, self%time_dimname, &
            self%time_dimid)
        call ncfile_check_error(ncstat)

        ! Set variable IDs
        ncstat = nf90_inq_varid(self%ncid, self%time_varname, &
            self%time_varid)
        call ncfile_check_error(ncstat)
    end procedure ncdg_file_open

    module procedure ncdg_file_close
        ! See interface for arguments and documentation

        ! Call parent function
        call self%ncfile%close()

        ! Flush dimension IDs
        self%time_dimid = -1

        ! Flush variable IDs
        self%time_varid = -1
    end procedure ncdg_file_close

    module procedure ncdg_file_set_metadata
        ! See interface for arguments and documentation

        ! Define dimensions and their attributes
        ncstat = nf90_def_dim(self%ncid, self%time_dimname, &
            time_dimsize, self%time_dimid)
        call ncfile_check_error(ncstat)

        ! Define variables and their attributes
        ncstat = nf90_def_var(self%ncid, self%time_varname, &
            nf90_double, self%time_dimid, varid=self%time_varid)
        call ncfile_check_error(ncstat)

        ! Define attributes
        ! N/A, for now
    end procedure ncdg_file_set_metadata

end submodule netcdf_fort_common