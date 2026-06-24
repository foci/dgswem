module parallelio_io
    !! Module for interacting with NetCDF files in parallel.

    use mpi, only : mpi_comm_world
    use pio, only : iosystem_desc_t, pio_init, pio_finalize, &
        pio_iotype_netcdf4p, pio_rearr_box, pio_rearr_subset, pio_write
    use piodg, only : piodg_63, piodg_64, piodg_maxele
    use pioutil, only : piofile_check_error

    implicit none

    ! Committing the sin of global variables, for now
    type(iosystem_desc_t) :: piosystem
    type(piodg_63) :: my_63
    type(piodg_64) :: my_64
    type(piodg_maxele) :: my_maxele

    contains

    subroutine parallelio_open_files_parallel()
        !! Initialize all NetCDF output files for the simulation in parallel.

        use sizes, only : mnproc, myproc

        implicit none

        ! Initialize PIO
        call pio_init( &
            comp_rank=myproc, &
            comp_comm=mpi_comm_world, &
            num_iotasks=mnproc, &
            num_aggregator=mnproc, &
            stride=1, &
            rearr=pio_rearr_subset, &
            iosystem=piosystem &
        )

        ! fort.63.nc
        call my_63%open( &
            piosystem, &
            piotype=pio_iotype_netcdf4p, &
            omode=pio_write &
        )

        ! fort.64.nc
        call my_64%open( &
            piosystem, &
            piotype=pio_iotype_netcdf4p, &
            omode=pio_write &
        )

        ! maxele.63.nc
        call my_maxele%open( &
            piosystem, &
            piotype=pio_iotype_netcdf4p, &
            omode=pio_write &
        )
    end subroutine parallelio_open_files_parallel

    subroutine parallelio_close_files_parallel()
        implicit none

        integer :: piostat

        ! Close files
        call my_63%close()
        call my_64%close()
        call my_maxele%close()

        ! Finalize PIO
        call pio_finalize(piosystem, piostat)
        call piofile_check_error(piostat)
    end subroutine parallelio_close_files_parallel

end module parallelio_io
