module parallelio_io
    !! Module for interacting with NetCDF files in parallel.

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

    use mpi, only : mpi_comm_world
    use pio, only : iosystem_desc_t, io_desc_t, pio_init, pio_finalize, &
        pio_iotype_netcdf4p, pio_rearr_box, pio_rearr_subset, pio_write, &
        pio_double, pio_initdecomp
    use piodg, only : piodg_ele, piodg_vel, piodg_pr, piodg_wvel, piodg_maxele
    use pioutil, only : piofile_check_error

    implicit none

    ! Committing the sin of global variables, for now

    ! PIODG variables
    type(piodg_ele) :: my_ele
    type(piodg_vel) :: my_vel
    type(piodg_pr) :: my_pr
    type(piodg_wvel) :: my_wvel
    type(piodg_maxele) :: my_maxele

    ! PIO variables
    type(iosystem_desc_t) :: my_iosystem
    type(iosystem_desc_t) :: my_iosystems(1) ! Only 1 kind of compute task
    type(io_desc_t) :: my_iodesc

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

    subroutine parallelio_init_parallel()
        !! Initialize ParallelIO system. This filters out processes for
        !! computing and writing, as well as handles domain decomposition.
        !! Called immediately before time stepping loop.
        !!
        !! ## Single Communicator
        !! If the number of processes is equal to the number of subdomains,
        !! every process will return from pio_init and proceed with the DGSWEM
        !! program.
        !!
        !! ## Computing and Writing Communicators
        !! If the number of processes is greater than the number of subdomains,
        !! only the first X processes will return and proceed with computing
        !! their subdomains. The surplus will perform the actual writing to
        !! file, and do not return from pio_init until the compute processes
        !! collectively call pio_finalize. Because of this, independent write
        !! processes call mpi_finalize and stop within this subroutine.

        use sizes, only : mnproc, myproc, nsubdom, nnodg, imap_nod_gh0, &
            comm_comp, comm_io
        use messenger_elem, only : message_fini

        implicit none

        integer :: i ! For looping
        integer :: my_comm_comps(1) ! For nprocs > 2, and only 1 kind of program

        ! Check number of processes vs. subdomains
        ! Option 1: compute processes write their own data
        ! Option 2: surplus writer processes aggregate and write data
        if (mnproc < nsubdom) then
            if (myproc == 0) print *, "Error: Fewer processes than domains."
            stop
        else if (mnproc == nsubdom) then
            if (myproc == 0) print *, "No surplus, so every process will write."
            call pio_init( &
                comp_rank=myproc, &
                comp_comm=mpi_comm_world, &
                num_iotasks=mnproc, &
                num_aggregator=mnproc, &
                stride=1, &
                rearr=pio_rearr_subset, &
                iosystem=my_iosystem &
            )
            comm_comp = mpi_comm_world
            comm_io = mpi_comm_world
            call pio_initdecomp( &
                iosystem=my_iosystem, &
                basepiotype=pio_double, &
                dims=[nnodg], &
                compdof=imap_nod_gh0, &
                iodesc=my_iodesc &
            )
        else if (mnproc > nsubdom) then
            if (myproc == 0) print *, "Surplus of ", mnproc - nsubdom, &
                " processes detected for writing."
            call pio_init( &
                iosystem=my_iosystems, & ! One for each compute program type
                incomm=mpi_comm_world, &
                procs_per_component=[nsubdom], & ! Only one kind of computation
                comp_proc_list=reshape([(i, i=0, nsubdom - 1)], &
                    shape=[nsubdom, 1]), & ! One list for each program type
                io_proc_list=[(i, i=nsubdom, mnproc - 1)], &
                rearranger=pio_rearr_subset, &
                comp_comm=my_comm_comps, & ! One for each compute program type
                io_comm=comm_io &
            ) ! I/O processes do not return from this call
            my_iosystem = my_iosystems(1)
            comm_comp = my_comm_comps(1)
            ! comm_io already set
            if (myproc < nsubdom) then
                call pio_initdecomp( &
                    iosystem=my_iosystem, &
                    basepiotype=pio_double, &
                    dims=[nnodg], &
                    compdof=imap_nod_gh0, &
                    iodesc=my_iodesc &
                )
            else
                ! Writer processes end here
                call message_fini()
                stop
            end if
        end if
    end subroutine parallelio_init_parallel

    subroutine parallelio_open_files_parallel()
        !! Initialize all NetCDF output files for the simulation in parallel.
        !! Called immediately before the time stepping loop.

        use global, only : &
            noute, noutv, noutc, noutm, & ! Whether to write to NetCDF
            noutge, noutgv, noutgc, noutgw
        use sizes, only : myproc, nsubdom

        implicit none

        ! Just in case, filter only compute processes
        if (myproc < nsubdom) then

            ! Set counters to 0
            nscoue = 0
            nscouv = 0
            nscouc = 0
            nscoum = 0
            nscouge = 0
            nscougv = 0
            nscougc = 0
            nscougw = 0

            ! Open files

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
                call my_ele%open(my_iosystem, omode=pio_write)
                call my_maxele%open(my_iosystem, omode=pio_write)
            end if

            ! nodal velocity
            if (abs(noutgv) == 3) then
                call my_vel%open(my_iosystem, omode=pio_write)
            end if

            ! nodal concentration
            if (abs(noutgc) == 3) then
            end if

            ! nodal air pressure/wind velocity
            if (abs(noutgw) == 3) then
                call my_pr%open(my_iosystem, omode=pio_write)
                call my_wvel%open(my_iosystem, omode=pio_write)
            end if
        end if
    end subroutine parallelio_open_files_parallel

    subroutine parallelio_write_step_parallel(it, force_write)
        !! Write current state to all NetCDF output files for the simulation.

        use global, only : &
            eta2, etamax, np, time_a, uu2, vv2, &
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
                        scalar=eta2(1:np), &
                        piodesc=my_iodesc, &
                        sync=.true. &
                    )
                    call my_maxele%write_step( &
                        t=time_a, &
                        scalar=etamax(1:np), &
                        np=np, &
                        piodesc=my_iodesc, &
                        compute_max=.false., &
                        sync=.true. &
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
                        vector_u=uu2(1:np), &
                        vector_v=vv2(1:np), &
                        piodesc=my_iodesc, &
                        sync=.true. &
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
                        scalar=pr2(1:np), &
                        piodesc=my_iodesc, &
                        sync=.true. &
                    )
                    call my_wvel%write_step( &
                        t=time_a, &
                        vector_u=wvnxout(1:np), &
                        vector_v=wvnyout(1:np), &
                        piodesc=my_iodesc, &
                        sync=.true. &
                    )
                    nscougw = 0
                end if
            end if
        end if
    end subroutine parallelio_write_step_parallel

    subroutine parallelio_close_files_parallel()
        !! Called immediately after the time stepping loop.

        use global, only : &
            noute, noutv, noutc, noutm, &
            noutge, noutgv, noutgc, noutgw
        use sizes, only : imap_nod, imap_nod_gh0, resnode_pio

        implicit none

        integer :: piostat

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

        ! Deallocate associated arrays
        if (allocated(imap_nod)) deallocate(imap_nod)
        if (allocated(imap_nod_gh0)) deallocate(imap_nod_gh0)
        if (allocated(resnode_pio)) deallocate(resnode_pio)

        ! Finalize PIO
        call pio_finalize(my_iosystem, piostat)
        call piofile_check_error(piostat)
    end subroutine parallelio_close_files_parallel

end module parallelio_io
