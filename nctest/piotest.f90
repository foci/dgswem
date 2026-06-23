program piotest
    !! Test program for PIO development

    use mpi
    use pioutil, only : piofile
    use pio

    implicit none

    ! MPI variables
    integer :: nprocs
    integer :: procno
    integer :: ierr

    ! PIO variables
    type(iosystem_desc_t) :: piosystem
    logical :: has_netcdf4p

    ! PIODG variables
    type(piofile) :: my_piofile
    type(piofile) :: my_piofile2

    ! Initialize MPI
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world, procno, ierr)
    call mpi_comm_size(mpi_comm_world, nprocs, ierr)

    if (procno == 0) print *, "Roll call!"
    call mpi_barrier(mpi_comm_world, ierr)
    print *, "Hello, I am processor ", procno, " of ", nprocs, "!"
    call mpi_barrier(mpi_comm_world, ierr)

    has_netcdf4p = PIO_iotype_available(PIO_iotype_netcdf4p)
    if (procno == 0) then
        print *, "Is parallel NetCDF-4 supported? ", has_netcdf4p
    end if
    call mpi_barrier(mpi_comm_world, ierr)
    ! DOES NOT WORK
    ! CONDA INSTALLATION HAS PREPROCESSOR FLAG BUG
    ! USE PNETCDF

    ! Initialize PIO
    ierr = pio_set_log_level(5_c_int)
    call pio_init( &
        comp_rank=procno, &
        comp_comm=mpi_comm_world, &
        num_iotasks=nprocs, &
        num_aggregator=nprocs, &
        stride=1, &
        rearr=pio_rearr_subset, &
        iosystem=piosystem &
    )
    call pio_seterrorhandling(piosystem, pio_return_error)

    call my_piofile%create(piosystem, path="file1.nc", &
        piotype=pio_iotype_netcdf4p, cmode=pio_clobber)
    call my_piofile%close()
    call my_piofile%open(piosystem)
    call my_piofile%close()
    call my_piofile%create(piosystem, path="file2.nc", &
        piotype=pio_iotype_netcdf, cmode=pio_clobber)
    call my_piofile%close()
    call my_piofile2%open(piosystem, path="file1.nc", &
        piotype=pio_iotype_netcdf, omode=pio_write)

    ! Finalize PIO
    call pio_finalize(piosystem, ierr)

    ! Finalize MPI
    call mpi_finalize(ierr)

end program piotest
