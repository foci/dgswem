program nctest_parallel

    use ncdg, only: ncdg_63, ncdg_64, ncdg_maxele
    use mpi
    use netcdf ! For now

    implicit none

    ! MPI variables
    integer :: nprocs
    integer :: procno
    integer :: ierr

    type(ncdg_63) :: my_63
    type(ncdg_64) :: my_64
    type(ncdg_maxele) :: my_maxele

    ! Bare-bones mesh consisting of a single triangle
    integer, parameter :: ics = 1 ! Cartesian
    integer, parameter :: np = 3 ! Number of nodes
    integer, parameter :: ne = 1 ! Number of elements
    integer, parameter :: nhy = 3 ! Number of nodes per element
    integer :: nm(ne, nhy) ! Node numbers for each element

    ! Coordinates and bathymetry
    real :: x(np) = [0.0, 1.0, 1.0]
    real :: y(np) = [0.0, 0.0, 1.0]
    real :: dp(np) = [0.0, 1.0, 2.0]

    ! Elevation boundary
    integer, parameter :: nope = 1 ! Number of elevation segments
    integer, parameter :: neta = 2 ! Number of elevation nodes
    integer :: ibtypee(nope) = [0] ! Segment types
    integer :: nvdll(nope) = [2] ! Nodes per segment
    integer :: nbdv(nope, neta) ! Node numbers for each segment

    ! Normal flow (discharge) boundary
    integer, parameter :: nbou = 1 ! Number of normal flow segments
    integer, parameter :: nvel = 3 ! Number of normal flow nodes
    integer :: ibtype(nbou) = [1] ! Segment types
    integer :: nvell(nbou) = [3] ! Nodes per segment
    integer :: nbvv(nbou, nvel) ! Node numbers for each segment

    ! Fill in arrays
    nm(1, :) = [1, 2, 3]
    nbdv(1, :) = [1, 2]
    nbvv(1, :) = [2, 3, 1]

    ! Initialize MPI
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world, procno, ierr)
    call mpi_comm_size(mpi_comm_world, nprocs, ierr)

    print *, "Hello, I am processor ", procno, " of ", nprocs, "!"

    ! Create files
    call my_63%pinit()
    call my_64%pinit()
    call my_maxele%pinit()

    ! Close
    call my_63%close()
    call my_64%close()
    call my_maxele%close()

    ! Open again
    ! call my_63%popen(mode=ior(nf90_write, nf90_mpiio))
    ! call my_64%popen(mode=ior(nf90_write, nf90_mpiio))
    ! call my_maxele%popen(mode=ior(nf90_write, nf90_mpiio))

    ! Write mesh
    ! First time step
    ! Second time step
    ! Third times step
    ! Close again
    ! call my_63%close()
    ! call my_64%close()
    ! call my_maxele%close()

    ! Finalize MPI
    call mpi_finalize(ierr)

end program nctest_parallel
