module parallelio_io
    !! Module for interacting with NetCDF files in parallel.

    use mpi, only : mpi_comm_world
    use pio, only : iosystem_desc_t, io_desc_t, pio_init, pio_finalize, &
        pio_iotype_netcdf4p, pio_rearr_box, pio_rearr_subset, pio_write, &
        pio_double, pio_initdecomp
    use piodg, only : piodg_63, piodg_64, piodg_maxele
    use pioutil, only : piofile_check_error

    implicit none

    ! Committing the sin of global variables, for now

    ! PIODG variables
    type(piodg_63) :: my_63
    type(piodg_64) :: my_64
    type(piodg_maxele) :: my_maxele

    ! PIO variables
    type(iosystem_desc_t) :: my_iosystem
    type(iosystem_desc_t) :: my_iosystems(1) ! Only 1 kind of compute task
    type(io_desc_t) :: my_iodesc
    integer :: my_comm_comps(1) ! For nprocs > 2, and only 1 kind of program
    integer :: my_comm_io ! For nprocs > 2

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

        use sizes, only : mnproc, myproc, nsubdom, nnodg, imap_nod_gh0

        implicit none

        integer :: i ! For looping

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
            call pio_initdecomp( &
                iosystem=my_iosystem, &
                basepiotype=pio_double, &
                dims=[nnodg], &
                compdof=imap_nod_gh0, &
                iodesc=my_iodesc &
            )
        else if (mnproc > nsubdom) then
            if (myproc == 0) print *, "Surplus of ", mnproc - 2, &
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
                io_comm=my_comm_io &
            ) ! I/O processes do not return from this call
            my_iosystem = my_iosystems(1)
            if (myproc < nsubdom) then
                call pio_initdecomp( &
                    iosystem=my_iosystem, &
                    basepiotype=pio_double, &
                    dims=[nnodg], &
                    compdof=imap_nod_gh0, &
                    iodesc=my_iodesc &
                )
            end if
        end if
    end subroutine parallelio_init_parallel

    subroutine parallelio_open_files_parallel()
        !! Initialize all NetCDF output files for the simulation in parallel.
        !! Called immediately before the time stepping loop.

        use sizes, only : myproc, nsubdom

        implicit none

        ! Just in case, filter only compute processes
        if (myproc < nsubdom) then
            call my_63%open(my_iosystem, omode=pio_write)
            call my_64%open(my_iosystem, omode=pio_write)
            call my_maxele%open(my_iosystem, omode=pio_write)
        end if
    end subroutine parallelio_open_files_parallel

    subroutine parallelio_write_step_parallel()
        !! Write current state to all NetCDF output files for the simulation.

        use global, only : eta2, etamax, np, time_a, uu2, vv2

        implicit none

        call my_63%piodg_63_write_step( &
            t=time_a, &
            zeta=eta2, &
            piodesc=my_iodesc, &
            sync=.true. &
        )
        call my_64%piodg_64_write_step( &
            t=time_a, &
            u_vel=uu2, &
            v_vel=vv2, &
            piodesc=my_iodesc, &
            sync=.true. &
        )
        call my_maxele%piodg_maxele_write_step( &
            t=time_a, &
            zeta_max=etamax, &
            np=np, &
            piodesc=my_iodesc, &
            sync=.true. &
        )
    end subroutine parallelio_write_step_parallel

    subroutine parallelio_close_files_parallel()
        !! Called immediately after the time stepping loop.
        implicit none

        integer :: piostat

        ! Close files
        call my_63%close()
        call my_64%close()
        call my_maxele%close()

        ! Finalize PIO
        call pio_finalize(my_iosystem, piostat)
        call piofile_check_error(piostat)
    end subroutine parallelio_close_files_parallel

end module parallelio_io
