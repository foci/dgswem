module netcdf_prep
    !! Module for initializing NetCDF files.

    use netcdf, only : nf90_unlimited
    use ncdg, only : ncdg_63, ncdg_64, ncdg_maxele

    implicit none

    ! Committing the sin of global variables, for now
    type(ncdg_63) :: my_63
    type(ncdg_64) :: my_64
    type(ncdg_maxele) :: my_maxele

    contains

    subroutine netcdf_prep_files_serial()
        !! Initialize and close all NetCDF output files for the simulation.
        !!
        !! @warning "Serial Use Only"
        !! This subroutine is not designed for use with MPI. If used in parallel,
        !! it should only be called from a single processor.
        !! @endwarning

        use pre_global, only : dp, ics, nbdv, nbou, nbvv, nelg, neta, &
            nneg, nope, nnodg, nvel, nvdll, nvell, ibtype, x, y

        implicit none

        integer, allocatable :: ibtypee(:) ! Not included in prep
        integer :: nhy = 3 ! Not included in prep

        allocate(ibtypee(nope))

        ! Only current supported value is 0
        ibtypee = 0 ! This is valid syntax for array of 0s

        ! fort.63.nc
        call my_63%init()
        call my_63%ncdg_63_set_metadata( &
            nt=nf90_unlimited, &
            np=nnodg, &
            ne=nelg, &
            nhy=nhy, &
            nope=nope, &
            neta=neta, &
            max_nvdll=maxval(nvdll), &
            nbou=nbou, &
            nvel=nvel, &
            max_nvell=maxval(nvell), &
            ics=ics &
        )
        call my_63%write_mesh( &
            x=x, &
            y=y, &
            dp=dp, &
            nm=nneg, &
            nhy=nhy, &
            ne=nelg, &
            nope=nope, &
            neta=neta, &
            ibtypee=ibtypee, &
            nvdll=nvdll, &
            nbdv=nbdv, &
            nbou=nbou, &
            nvel=nvel, &
            ibtype=ibtype, &
            nvell=nvell, &
            nbvv=nbvv(:, 1:) &
        )
        call my_63%close()

        ! fort.64.nc
        call my_64%init()
        call my_64%ncdg_64_set_metadata( &
            nt=nf90_unlimited, &
            np=nnodg, &
            ne=nelg, &
            nhy=nhy, &
            nope=nope, &
            neta=neta, &
            max_nvdll=maxval(nvdll), &
            nbou=nbou, &
            nvel=nvel, &
            max_nvell=maxval(nvell), &
            ics=ics &
        )
        call my_64%write_mesh( &
            x=x, &
            y=y, &
            dp=dp, &
            nm=nneg, &
            nhy=nhy, &
            ne=nelg, &
            nope=nope, &
            neta=neta, &
            ibtypee=ibtypee, &
            nvdll=nvdll, &
            nbdv=nbdv, &
            nbou=nbou, &
            nvel=nvel, &
            ibtype=ibtype, &
            nvell=nvell, &
            nbvv=nbvv(:, 1:) &
        )
        call my_64%close()

        ! maxele.63.nc
        call my_maxele%init()
        call my_maxele%ncdg_maxele_set_metadata( &
            nt=nf90_unlimited, &
            np=nnodg, &
            ne=nelg, &
            nhy=nhy, &
            nope=nope, &
            neta=neta, &
            max_nvdll=maxval(nvdll), &
            nbou=nbou, &
            nvel=nvel, &
            max_nvell=maxval(nvell), &
            ics=ics &
        )
        call my_maxele%write_mesh( &
            x=x, &
            y=y, &
            dp=dp, &
            nm=nneg, &
            nhy=nhy, &
            ne=nelg, &
            nope=nope, &
            neta=neta, &
            ibtypee=ibtypee, &
            nvdll=nvdll, &
            nbdv=nbdv, &
            nbou=nbou, &
            nvel=nvel, &
            ibtype=ibtype, &
            nvell=nvell, &
            nbvv=nbvv(:, 1:) &
        )
        call my_maxele%close()

        deallocate(ibtypee)
    end subroutine netcdf_prep_files_serial
end module netcdf_prep
