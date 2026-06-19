program piotest
    !! Test program for PIO development

    use mpi
    use pioutil, only : piofile

    implicit none

    ! MPI variables
    integer :: nprocs
    integer :: procno
    integer :: ierr

    ! PIO variables

    ! NCDG variables
    type(piofile) :: my_piofile

    ! Initialize MPI
    ! Initialize PIO

    call my_piofile%open()

    ! Finalize MPI
    ! Finalize PIO
end program piotest
