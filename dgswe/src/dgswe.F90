      PROGRAM swe_tri

!$    USE omp_lib      
      USE globals
      USE basis
      USE allocation
      USE messenger2
      USE read_dginp
      USE output
      USE quit
      USE area_qpts_mod
      USE edge_qpts_mod

      IMPLICIT NONE
      INTEGER :: it
      INTEGER :: cnt_sol,cnt_sta
      INTEGER :: i,j,et
      INTEGER :: tstep
      REAL(rp) :: t_start,t_end
      


      
      ! Initialize MPI and determine rank if appropriate
      CALL message_init()
      
      ! Read in the keyword (or fixed) format input file 
      CALL read_input(myrank,dirname)        

      ! Compute # of dofs and # of nodes straight/curved elements etc.
      CALL sizes()

      ! Read in grid file
      CALL read_grid()

      ! Read in forcing file
      CALL read_forcing()     
   
      ! Find edge connectivity
      CALL connect()

      ! Get quadrature points for area integration
      CALL area_qpts(myrank,p,ctp,nel_type,nqpta,mnqpta,wpta,qpta)

      ! Calculate basis function values at edge quadrature points
      CALL edge_qpts(myrank,p,ctp,nel_type,nqpte,mnqpte,wpte,qpte)   
      
      ! Calculate basis function and derivative values at area quadrature points
      CALL area_basis()  

      ! Calculate basis function and derivative values at edge quadrature points
      CALL edge_basis()            
      
      ! Compute element area, edge length, edge normals, and bathymetry derivatives
      CALL element_data()                 

      ! Compute initial condition
      CALL initial()

      ! Boundary forcing interpolation
      CALL interp_forcing()
      
      ! Partition domain for element/edge blocking
      CALL metis2(npart)
      
      ! Decompose domain and prep element/edge blocking
      CALL edge_partition2()
      
      ! Read the local-to-global element and message passing files
      CALL read_message_files()      
      
      ! Set up send/recieve buffers and edge data structures
      CALL message_setup()
      
      ! Initialize the MPI persistent message passing calls
      CALL communication_setup()    

!       OPEN(unit=195,file=trim(out_direc) // 'edge_check.d')        
! 
!       DO et = 1,2*nel_type
!         DO i = 1,mnnds
!           WRITE(195,"(100(ES24.17,1x))") (dpsids(i,j,et), j = 1,mnqpta+4*mnqpte)
!         ENDDO
!       WRITE(195,*) " "
!       ENDDO      
      
#ifdef openmp
      t_start = omp_get_wtime()
#else
      CALL CPU_TIME(t_start)
#endif
      

      
      tf = tf*86400d0
      tstep = int(tf/dt)           

      t = 0d0
      cnt_sol = 0
      cnt_sta = 0
 
      CALL output_solution(.true.)
      CALL output_stations(.true.)
      
      
      
      IF (myrank == 0) THEN
        PRINT "(A)", " "
        PRINT "(A)", "---------------------------------------------"
        PRINT "(A)", "               Time Stepping                 "
        PRINT "(A)", "---------------------------------------------"
        PRINT "(A)", " "

        PRINT "(A,e12.4)", "Time step: ",dt
        PRINT "(A,e12.4)", "Final time: ",tf

        PRINT "(A)", " "
      ENDIF           
      
#ifdef CMPI      
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)     
#endif      

      DO it = 1,tstep

        CALL rk()

         t = t + dt

         cnt_sol = cnt_sol + 1
         IF(sol_opt > 0 .and. cnt_sol == tskp_sol) THEN
           CALL output_solution(.false.)             
           cnt_sol = 0
         ENDIF
         
         cnt_sta = cnt_sta + 1
         IF(sta_opt > 0 .and. cnt_sta == tskp_sta) THEN
           CALL output_stations(.false.)             
           cnt_sta = 0
         ENDIF         

      ENDDO

      
      CALL end_time(t_start,nproc)
      
     

      
!       OPEN(UNIT=101, FILE='CPUtime.log', STATUS='OLD', POSITION='APPEND')
!       WRITE(101,*) "grid = ", grid_file
!       WRITE(101,*) "p = ", p
!       WRITE(101,*) "dt, tf = ", dt,tf
! #ifdef rk22
!       WRITE(101,*) "RK 22"
! #endif
! #ifdef rk33
!       WRITE(101,*) "RK 33"
! #endif
!       WRITE(101,*) "nblk, npart = ", nblk, npart
!       WRITE(101,*) "mnpartel,mnparted = ", mnpartel,mnparted   
!       WRITE(101,*) "CPU time = ", t_end-t_start
!       WRITE(101,*) " "
!       CLOSE(101)      
      
      
      
      CALL close_output()
      CALL finish(myrank)
      
      END PROGRAM swe_tri
