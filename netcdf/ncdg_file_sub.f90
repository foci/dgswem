submodule (ncdg) ncdg_file_sub
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
        
        integer :: ncstat ! Status of most recent operation

        ! Define dimensions
        ncstat = nf90_def_dim(self%ncid, self%time_dimname, &
            time_dimsize, self%time_dimid)
        call ncfile_check_error(ncstat)

        ! Define variables and their attributes
        
        ! time
        ncstat = nf90_def_var(self%ncid, self%time_varname, &
            nf90_double, self%time_dimid, varid=self%time_varid)
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%time_varid, &
            "long_name", "model time")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%time_varid, &
            "standard_name", "time")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%time_varid, &
            "units", "seconds since 1970-01-01 00:00:00")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%time_varid, &
            "base_date", "1970-01-01 00:00:00")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, self%time_varid, &
            "axis", "T")
        call ncfile_check_error(ncstat)

        ! Define global attributes
        ncstat = nf90_put_att(self%ncid, nf90_global, &
            "model", "DGSWEM")
        call ncfile_check_error(ncstat)
        ! Begin excluded block
        ! Implementation deferred
        if (.false.) then
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "version", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "git_hash", "")
            call ncfile_check_error(ncstat)
        end if
        ! End excluded block
        ncstat = nf90_put_att(self%ncid, nf90_global, &
            "grid_type", "Triangular")
        call ncfile_check_error(ncstat)
        ! Begin excluded block
        ! Implementation deferred
        if (.false.) then
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "description", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "agrid", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "title", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "institution", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "source", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "history", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "comments", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "host", "")
            call ncfile_check_error(ncstat)
        end if
        ! End excluded block
        ncstat = nf90_put_att(self%ncid, nf90_global, &
            "convention", "CF")
        call ncfile_check_error(ncstat)
        ncstat = nf90_put_att(self%ncid, nf90_global, &
            "Conventions", "UGRID-1.0")
        call ncfile_check_error(ncstat)
        ! Begin excluded block
        ! Implementation deferred
        if (.false.) then
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "contact", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "creation_date", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "modification_date", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "dt", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "ihot", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "ics", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "nolibf", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "nolifa", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "nolica", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "nolicat", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "nwp", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "ncor", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "ntip", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "nws", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "nramp", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "tau0", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "statim", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "reftim", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "rnday", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "dramp", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "a00", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "b00", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "c00", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "h0", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "slam0", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "sfea0", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "cf", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "eslm", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "cori", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "ntif", "")
            call ncfile_check_error(ncstat)
            ncstat = nf90_put_att(self%ncid, nf90_global, &
                "nbfr", "")
            call ncfile_check_error(ncstat)
        end if
        ! End excluded block
    end procedure ncdg_file_set_metadata

end submodule ncdg_file_sub