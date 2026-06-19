program nctest_parallel

    use ncdg, only: ncdg_63, ncdg_64, ncdg_maxele
    use mpi
    use pio

    implicit none

    ! MPI variables
    integer :: nprocs
    integer :: procno
    integer :: ierr

    type(ncdg_63) :: my_63
    type(ncdg_64) :: my_64
    type(ncdg_maxele) :: my_maxele

    ! Mesh (local)
    integer :: ics ! Cartesian
    integer :: np ! Number of nodes
    integer :: ne ! Number of elements
    integer :: nhy ! Number of nodes per element
    integer, allocatable :: nm(:, :) ! Node numbers for each element

    ! Coordinates and bathymetry
    real, allocatable :: x(:)
    real, allocatable :: y(:)
    real, allocatable :: dp(:)

    ! Elevation boundary
    integer :: nope ! Number of elevation segments
    integer :: neta ! Number of elevation nodes
    integer, allocatable :: ibtypee(:) ! Segment types
    integer, allocatable :: nvdll(:) ! Nodes per segment
    integer, allocatable :: nbdv(:, :) ! Node numbers for each segment

    ! Normal flow (discharge) boundary
    integer :: nbou ! Number of normal flow segments
    integer :: nvel ! Number of normal flow nodes
    integer, allocatable :: ibtype(:) ! Segment types
    integer, allocatable :: nvell(:) ! Nodes per segment
    integer, allocatable :: nbvv(:, :) ! Node numbers for each segment

    ! Parallell decomposition variables
    type(iosystem) :: my_iosystem
    type(iodesc) :: my_iodesc
    integer, allocatable :: imap_nod(:) ! Mapping to global node index
    logical, allocatable :: resnode(:) ! Whether local nodes are resident or ghost
    logical, allocatable :: resele(:) ! Whether local elements are resident or ghost

    ! Initialize MPI
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world, procno, ierr)
    call mpi_comm_size(mpi_comm_world, nprocs, ierr)

    if (procno == 0) print *, "Roll call!"
    call mpi_barrier(mpi_comm_world, ierr)
    print *, "Hello, I am processor ", procno, " of ", nprocs, "!"
    call mpi_barrier(mpi_comm_world, ierr)

    if (procno == 0) then
        print *, "Creating files and writing mesh like adcprep"
        ! Mesh (global)
        ics = 1 ! Cartesian
        np = 6 ! Number of nodes
        ne = 4 ! Number of elements
        nhy = 3 ! Number of nodes per element
        allocate(nm(ne, nhy)) ! Node numbers for each element
        nm(1, :) = [1, 3, 2]
        nm(2, :) = [3, 4, 2]
        nm(3, :) = [3, 5, 4]
        nm(4, :) = [5, 6, 4]

        ! Coordinates and bathymetry
        allocate(x(np))
        x(:) = [0.0, 0.0, 1.0, 1.0, 2.0, 2.0]
        allocate(y(np))
        y(:) = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
        allocate(dp(np)) ! Making this global node index
        dp(:) = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]

        ! Elevation boundary
        nope = 1 ! Number of elevation segments
        neta = 4 ! Number of elevation nodes
        allocate(ibtypee(nope)) ! Segment types
        ibtypee(:) = [0]
        allocate(nvdll(nope)) ! Nodes per segment
        nvdll(:) = [4]
        allocate(nbdv(nope, neta)) ! Node numbers for each segment
        nbdv(1, :) = [1, 2, 4, 6]

        ! Normal flow (discharge) boundary
        nbou = 1 ! Number of normal flow segments
        nvel = 4 ! Number of normal flow nodes
        allocate(ibtype(nbou)) ! Segment types
        ibtype(:) = [0]
        allocate(nvell(nbou)) ! Nodes per segment
        nvell(:) = [4]
        allocate(nbvv(nbou, nvel)) ! Node numbers for each segment
        nbvv(1, :) = [1, 3, 5, 6]

        ! Initialization (in serial)
        call my_63%init()
        call my_64%init()
        call my_maxele%init()

        ! Set metadata (in serial)
        call my_63%ncdg_63_set_metadata( &
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
        call my_64%ncdg_64_set_metadata( &
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
        call my_maxele%ncdg_maxele_set_metadata( &
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

        ! Write mesh (in serial)
        call my_63%write_mesh( &
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
            nbvv=nbvv&
        )
        call my_64%write_mesh( &
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
            nbvv=nbvv&
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
            nbvv=nbvv&
        )

        ! Close (in serial)
        call my_63%close()
        call my_64%close()
        call my_maxele%close()

        ! Deallocate
        deallocate(nm)
        deallocate(x)
        deallocate(y)
        deallocate(dp)
        deallocate(ibtypee)
        deallocate(nvdll)
        deallocate(nbdv)
        deallocate(ibtype)
        deallocate(nvell)
        deallocate(nbvv)
    end if

    ! Make sure that file creation is complete before letting any process move
    ! on to the next step
    call mpi_barrier(mpi_comm_world, ierr)

    ! Set values based on process
    select case (procno)
        case (0)
            ! Mesh (local)
            ics = 1 ! Cartesian
            np = 5 ! Number of nodes
            ne = 3 ! Number of elements
            nhy = 3 ! Number of nodes per element
            allocate(nm(ne, nhy)) ! Node numbers for each element
            nm(1, :) = [3, 2, 4]
            nm(2, :) = [5, 3, 4]
            nm(3, :) = [3, 1, 2]

            ! Coordinates and bathymetry
            allocate(x(np))
            x(:) = [2.0, 1.0, 1.0, 0.0, 0.0]
            allocate(y(np))
            y(:) = [0.0, 1.0, 0.0, 1.0, 0.0]
            allocate(dp(np)) ! Making this global node index
            dp(:) = [5.0, 4.0, 3.0, 2.0, 1.0]

            ! Elevation boundary
            nope = 1 ! Number of elevation segments
            neta = 3 ! Number of elevation nodes
            allocate(ibtypee(nope)) ! Segment types
            ibtypee(:) = [0]
            allocate(nvdll(nope)) ! Nodes per segment
            nvdll(:) = [3]
            allocate(nbdv(nope, neta)) ! Node numbers for each segment
            nbdv(1, :) = [5, 4, 2]

            ! Normal flow (discharge) boundary
            nbou = 1 ! Number of normal flow segments
            nvel = 3 ! Number of normal flow nodes
            allocate(ibtype(nbou)) ! Segment types
            ibtype(:) = [0]
            allocate(nvell(nbou)) ! Nodes per segment
            nvell(:) = [3]
            allocate(nbvv(nbou, nvel)) ! Node numbers for each segment
            nbvv(1, :) = [5, 3, 1]

            ! Parallell decomposition variables
            allocate(imap_nod(np))
            imap_nod(:) = [5, 4, 3, 2, 1]
            allocate(resnode(np)) ! Whether local nodes are resident or ghost
            resnode = [.false., .true., .true., .true., .true.]
            allocate(resele(ne)) ! Whether local elements are resident or ghost
            resele = [.true., .true., .false.]
        case (1)
            ! Mesh (local)
            ics = 1 ! Cartesian
            np = 5 ! Number of nodes
            ne = 3 ! Number of elements
            nhy = 3 ! Number of nodes per element
            allocate(nm(ne, nhy)) ! Node numbers for each element
            nm(1, :) = [3, 1, 4]
            nm(2, :) = [3, 4, 5]
            nm(3, :) = [1, 2, 4]

            ! Coordinates and bathymetry
            allocate(x(np))
            x(:) = [2.0, 2.0, 1.0, 1.0, 0.0]
            allocate(y(np))
            y(:) = [0.0, 1.0, 0.0, 1.0, 1.0]
            allocate(dp(np)) ! Making this global node index
            dp(:) = [5.0, 6.0, 3.0, 4.0, 2.0]

            ! Elevation boundary
            nope = 1 ! Number of elevation segments
            neta = 3 ! Number of elevation nodes
            allocate(ibtypee(nope)) ! Segment types
            ibtypee(:) = [0]
            allocate(nvdll(nope)) ! Nodes per segment
            nvdll(:) = [3]
            allocate(nbdv(nope, neta)) ! Node numbers for each segment
            nbdv(1, :) = [5, 4, 2]

            ! Normal flow (discharge) boundary
            nbou = 1 ! Number of normal flow segments
            nvel = 3 ! Number of normal flow nodes
            allocate(ibtype(nbou)) ! Segment types
            ibtype(:) = [0]
            allocate(nvell(nbou)) ! Nodes per segment
            nvell(:) = [3]
            allocate(nbvv(nbou, nvel)) ! Node numbers for each segment
            nbvv(1, :) = [3, 1, 2]

            ! Parallell decomposition variables
            allocate(imap_nod(np))
            imap_nod(:) = [5, 6, 3, 4, 2]
            allocate(resnode(np)) ! Whether local nodes are resident or ghost
            resnode(:) = [.true., .true., .false., .false., .false.]
            allocate(resele(ne)) ! Whether local elements are resident or ghost
            resele(:) = [.true., .false., .true.]
        case default
            ! Do nothing
    end select

    ! Initialize PIO and decompose
    select case (nprocs)
        case (:1)
            stop "Error: Fewer processes than domains."
        case (2)
            if (procno == 0) print *, "No surplus, so every process will write"
            call pio_init_intracomm( &
                comp_rank=procno, &
                comp_comm=mpi_comm_world, &ncells
                num_iotasks=2, &
                num_aggregator=2, &
                stride=1, &
                rearr=pio_rearr_box, &
                iosystem=my_iosystem &
            )
            call pio_initdecomp( &
                iosystem=my_iosystem, &
                basepiotype=pio_double, &
                dims=nnodg, &
                compdof=imap_nod, &
                iodesc=my_iodesc, &
            )
        case (3:)
            if (procno == 0) print *, "Surplus of ", nprocs - 2, " processes detected."
            stop "Not implemented yet"
        case default
            stop "Error: How did I get here?"
    end select

    select case (procno)
        case (0:1)
            ! Open files
            ! call my_63%pio_open()
            ! call my_64%pio_open()
            ! call my_maxele%pio_open()
            ! First time step
            ! Second time step
            ! Third time step
            ! Close
    end select

    ! Deallocate
    if (procno == 0 .or. procno == 1) then
        deallocate(nm)
        deallocate(x)
        deallocate(y)
        deallocate(dp)
        deallocate(ibtypee)
        deallocate(nvdll)
        deallocate(nbdv)
        deallocate(ibtype)
        deallocate(nvell)
        deallocate(nbvv)
        deallocate(resnode)
        deallocate(resele)
    end if

    ! Finalize PIO

    ! Finalize MPI
    call mpi_finalize(ierr)

end program nctest_parallel
