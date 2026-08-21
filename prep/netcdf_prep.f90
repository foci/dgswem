module netcdf_prep
    !! Module for initializing NetCDF files.

    use netcdf, only : nf90_unlimited
    use ncdg, only : ncdg_ele, ncdg_vel, ncdg_maxele

    implicit none

    ! Committing the sin of global variables, for now
    type(ncdg_ele) :: my_ele
    type(ncdg_vel) :: my_vel
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
            nneg, nope, nnodg, nvel, nvdll, nvell, ibtype, x, y, &
            noute, noutv, noutc, noutm, & ! Whether to write to NetCDF
            noutge, noutgv, noutgc, noutgw

        implicit none

        integer, allocatable :: ibtypee(:) ! Not included in prep
        integer :: nhy ! Not included in prep

        allocate(ibtypee(nope))

        nhy = 3

        ! Only current supported value is 0
        ibtypee = 0 ! This is valid syntax for array of 0s

        ! Create new files

        ! station elevation
        if (abs(noute) == 3) then
        end if

        ! station velocity
        if (abs(noutv) == 3) then
        end if

        ! station concentration
        if (abs(noutc) == 3) then
        end if

        ! station air pressure/wind velocity
        if (abs(noutm) == 3) then
        end if

        ! nodal elevation
        if (abs(noutge) == 3) then
            ! Create fort.63
            call my_ele%create()
            call my_ele%set_metadata( &
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
            call my_ele%write_mesh( &
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
                nbvv=nbvv(:, 1:), &
                transpose_nm=.false. &
            )
            call my_ele%close()
            ! Create maxele.63
            call my_maxele%create()
            call my_maxele%set_metadata( &
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
                nbvv=nbvv(:, 1:), &
                transpose_nm=.false. &
            )
            call my_maxele%close()
        end if

        ! nodal velocity
        if (abs(noutgv) == 3) then
            ! Create fort.64
            call my_vel%create()
            call my_vel%set_metadata( &
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
            call my_vel%write_mesh( &
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
                nbvv=nbvv(:, 1:), &
                transpose_nm=.false. &
            )
            call my_vel%close()
        end if

        ! nodal concentration
        if (abs(noutgc) == 3) then
        end if

        ! nodal air pressure/wind velocity
        if (abs(noutgw) == 3) then
        end if

        deallocate(ibtypee)
    end subroutine netcdf_prep_files_serial
end module netcdf_prep
