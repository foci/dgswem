module netcdf_io
    !! Module for interacting with NetCDF files.

    ! Output flag guide
    ! nout*: whether to write the file
    ! touts*: when to start
    ! toutf*: when to stop
    ! ntcys*: start interval
    ! ntcyf*: stop interval
    ! nspool* : write interval
    ! *e: station elevation
    ! *v: station velocity
    ! *c: station concentration
    ! *m: station air pressure/wind velocity
    ! *ge: nodal elevation
    ! *gv: nodal velocity
    ! *gc: nodal concentration
    ! *gw: nodal air pressure/wind velocity

    use netcdf, only : nf90_unlimited, nf90_write
    use ncdg, only : ncdg_ele, ncdg_vel, ncdg_pr, ncdg_wvel, ncdg_maxele

    implicit none

    ! Committing the sin of global variables, for now
    type(ncdg_ele) :: my_ele
    type(ncdg_vel) :: my_vel
    type(ncdg_pr) :: my_pr
    type(ncdg_wvel) :: my_wvel
    type(ncdg_maxele) :: my_maxele

    ! Output counters
    integer :: nscoue
    !! Local counter for station elevation output
    integer :: nscouv
    !! Local counter for station velocity output
    integer :: nscouc
    !! Local counter for station concentration output
    integer :: nscoum
    !! Local counter for station weather output
    integer :: nscouge
    !! Local counter for nodal elevation output
    integer :: nscougv
    !! Local counter for nodal velocity output
    integer :: nscougc
    !! Local counter for nodal concentration output
    integer :: nscougw
    !! Local counter for nodal weather output

    contains

    subroutine netcdf_init_files_serial()
        !! Create new and/or open existing NetCDF output files for the
        !! simulation.
        !!
        !! @warning "Serial Use Only"
        !! This subroutine is not designed for use with MPI. If used in parallel,
        !! it should only be called from a single processor.
        !! @endwarning

        use global, only : dp, ics, nbdv, nbou, nbvv, ne, neta, &
            nhy, nm, nope, np, nvel, nvdll, nvell, segtype, x, y, &
            noute, noutv, noutc, noutm, & ! Whether to write to NetCDF
            noutge, noutgv, noutgc, noutgw

        implicit none

        integer, allocatable :: ibtypee(:) ! Not included in src
        integer, allocatable :: ibtype(:) ! Not included in src

        logical :: file_exists ! Whether a particular file has been created yet

        allocate(ibtypee(nope))
        allocate(ibtype(nbou))

        ! Only current supported value is 0
        ibtypee = 0 ! This is valid syntax for array of 0s
        ibtype = segtype

        ! Set counters to 0
        nscoue = 0
        nscouv = 0
        nscouc = 0
        nscoum = 0
        nscouge = 0
        nscougv = 0
        nscougc = 0
        nscougw = 0

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
            ! Create or open fort.63
            inquire(file="fort.63.nc", exist=file_exists)
            if (noutge == -3 .or. .not. file_exists) then ! create new
                call my_ele%create()
                call my_ele%set_metadata( &
                    nt=nf90_unlimited, &
                    np=np, &
                    ne=ne, &
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
                    nm=nm, &
                    nhy=nhy, &
                    ne=ne, &
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
            else if (noutge == 3) then ! append to existing
                call my_ele%open(mode=nf90_write)
            end if
            ! Create or open maxele.63
            inquire(file="maxele.63.nc", exist=file_exists)
            if (noutge == -3 .or. .not. file_exists) then ! create new
                call my_maxele%create()
                call my_maxele%set_metadata( &
                    nt=nf90_unlimited, &
                    np=np, &
                    ne=ne, &
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
                    nm=nm, &
                    nhy=nhy, &
                    ne=ne, &
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
            else if (noutge == 3) then ! append to existing
                call my_maxele%open(mode=nf90_write)
            end if
        end if

        ! nodal velocity
        if (abs(noutgv) == 3) then
            ! Create or open fort.64.nc
            inquire(file="fort.64.nc", exist=file_exists)
            if (noutgv == -3 .or. .not. file_exists) then ! create new
                call my_vel%create()
                call my_vel%set_metadata( &
                    nt=nf90_unlimited, &
                    np=np, &
                    ne=ne, &
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
                    nm=nm, &
                    nhy=nhy, &
                    ne=ne, &
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
            else if (noutgv == 3) then ! append to existing
                call my_vel%open(mode=nf90_write)
            end if
        end if

        ! nodal concentration
        if (abs(noutgc) == 3) then
        end if

        ! nodal air pressure/wind velocity
        if (abs(noutgw) == 3) then
            ! Create or open fort.73
            inquire(file="fort.73.nc", exist=file_exists)
            if (noutgw == -3 .or. .not. file_exists) then ! create new
                call my_pr%create()
                call my_pr%set_metadata( &
                    nt=nf90_unlimited, &
                    np=np, &
                    ne=ne, &
                    nhy=nhy, &
                    nope=nope, &
                    neta=neta, &
                    max_nvdll=maxval(nvdll), &
                    nbou=nbou, &
                    nvel=nvel, &
                    max_nvell=maxval(nvell), &
                    ics=ics &
                )
                call my_pr%write_mesh( &
                    x=x, &
                    y=y, &
                    dp=dp, &
                    nm=nm, &
                    nhy=nhy, &
                    ne=ne, &
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
            else if (noutgw == 3) then ! append to existing
                call my_pr%open(mode=nf90_write)
            end if
            ! Create or open fort.74.nc
            inquire(file="fort.74.nc", exist=file_exists)
            if (noutgw == -3 .or. .not. file_exists) then ! create new
                call my_wvel%create()
                call my_wvel%set_metadata( &
                    nt=nf90_unlimited, &
                    np=np, &
                    ne=ne, &
                    nhy=nhy, &
                    nope=nope, &
                    neta=neta, &
                    max_nvdll=maxval(nvdll), &
                    nbou=nbou, &
                    nvel=nvel, &
                    max_nvell=maxval(nvell), &
                    ics=ics &
                )
                call my_wvel%write_mesh( &
                    x=x, &
                    y=y, &
                    dp=dp, &
                    nm=nm, &
                    nhy=nhy, &
                    ne=ne, &
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
            else if (noutgw == 3) then ! append to existing
                call my_wvel%open(mode=nf90_write)
            end if
        end if

        deallocate(ibtypee)
        deallocate(ibtype)
    end subroutine netcdf_init_files_serial

    subroutine netcdf_write_step_serial(it, force_write)
        !! Write current state to all NetCDF output files for the simulation.

        use global, only : &
            eta2, etamax, time_a, uu2, vv2, &
            pr2, wvnxout, wvnyout, &
            noute, noutv, noutc, noutm, & ! Whether to write to NetCDF
            noutge, noutgv, noutgc, noutgw, &
            ntcyse, ntcysv, ntcysc, ntcysm, & ! Start iteration
            ntcysge, ntcysgv, ntcysgc, ntcysgw, &
            ntcyfe, ntcyfv, ntcyfc, ntcyfm, & ! Stop iteration
            ntcyfge, ntcyfgv, ntcyfgc, ntcyfgw, &
            nspoole, nspoolv, nspoolc, nspoolm, & ! Write every
            nspoolge, nspoolgv, nspoolgc, nspoolgw

        implicit none

        integer :: it
        !! DGSWEM iteration index
        logical :: force_write
        !! Whether to force

        ! Write step
        ! The algorithm is
        ! 1. Check if the appropriate flag is +/- 3
        ! 2. Check if the iteration is in the appropriate range
        ! 3. Check if the iteration has reached the write interval
        ! 4. If all true, write!

        ! station elevation
        if (abs(noute) == 3) then
            if ((it > ntcyse) .and. (it <= ntcyfe) .or. force_write) then
                nscoue = nscoue + 1
                if (nscoue == nspoole .or. force_write) then
                    nscoue = 0
                end if
            end if
        end if

        ! station velocity
        if (abs(noutv) == 3) then
            if ((it > ntcysv) .and. (it <= ntcyfv) .or. force_write) then
                nscouv = nscouv + 1
                if (nscouv == nspoolv .or. force_write) then
                    nscouv = 0
                end if
            end if
        end if

        ! station concentration
        if (abs(noutc) == 3) then
            if ((it > ntcysc) .and. (it <= ntcyfc) .or. force_write) then
                nscouc = nscouc + 1
                if (nscouc == nspoolc .or. force_write) then
                    nscouc = 0
                end if
            end if
        end if

        ! station air pressure/wind velocity
        if (abs(noutm) == 3) then
            if ((it > ntcysm) .and. (it <= ntcyfm) .or. force_write) then
                nscoum = nscoum + 1
                if (nscoum == nspoolm .or. force_write) then
                    nscoum = 0
                end if
            end if
        end if

        ! nodal elevation
        if (abs(noutge) == 3) then
            if ((it > ntcysge) .and. (it <= ntcyfge) .or. force_write) then
                nscouge = nscouge + 1
                if (nscouge == nspoolge .or. force_write) then
                    call my_ele%write_step( &
                        t=time_a, &
                        scalar=eta2, &
                        sync=.false. &
                    )
                    call my_maxele%write_step( &
                        t=time_a, &
                        scalar=etamax, &
                        compute_max=.false., &
                        sync=.false. &
                    )
                    nscouge = 0
                end if
            end if
        end if

        ! nodal velocity
        if (abs(noutgv) == 3) then
            if ((it > ntcysgv) .and. (it <= ntcyfgv) .or. force_write) then
                nscougv = nscougv + 1
                if (nscougv == nspoolgv .or. force_write) then
                    call my_vel%write_step( &
                        t=time_a, &
                        vector_u=uu2, &
                        vector_v=vv2, &
                        sync=.false. &
                    )
                    nscougv = 0
                end if
            end if
        end if

        ! nodal concentration
        if (abs(noutgc) == 3) then
            if ((it > ntcysgc) .and. (it <= ntcyfgc) .or. force_write) then
                nscougc = nscougc + 1
                if (nscougc == nspoolgc .or. force_write) then
                    nscougc = 0
                end if
            end if
        end if

        ! nodal air pressure/wind velocity
        if (abs(noutgw) == 3) then
            if ((it > ntcysgw) .and. (it <= ntcyfgw) .or. force_write) then
                nscougw = nscougw + 1
                if (nscougw == nspoolgw .or. force_write) then
                    call my_pr%write_step( &
                        t=time_a, &
                        scalar=pr2, &
                        sync=.false. &
                    )
                    call my_wvel%write_step( &
                        t=time_a, &
                        vector_u=wvnxout, &
                        vector_v=wvnyout, &
                        sync=.false. &
                    )
                    nscougw = 0
                end if
            end if
        end if
    end subroutine netcdf_write_step_serial

    subroutine netcdf_close_files_serial()
        !! Close all NetCDF output files

        use global, only : &
            noute, noutv, noutc, noutm, & ! Whether to write to NetCDF
            noutge, noutgv, noutgc, noutgw

        implicit none

        ! Close files

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
            call my_ele%close()
            call my_maxele%close()
        end if

        ! nodal velocity
        if (abs(noutgv) == 3) then
            call my_vel%close()
        end if

        ! nodal concentration
        if (abs(noutgc) == 3) then
        end if

        ! nodal air pressure/wind velocity
        if (abs(noutgw) == 3) then
            call my_pr%close()
            call my_wvel%close()
        end if
    end subroutine netcdf_close_files_serial

end module netcdf_io
