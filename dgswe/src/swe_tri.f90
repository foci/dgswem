      PROGRAM swe_tri

      USE globals
      USE basis

      IMPLICIT NONE
      INTEGER :: it,tskp,cnt
      REAL(pres) :: tf,tstep,lines,t_start,t_end

      OPEN(unit=63,file='../output/solution_H.d')
      OPEN(unit=641,file='../output/solution_Qx.d')
      OPEN(unit=642,file='../output/solution_Qy.d')
      OPEN(unit=17,file='../output/edge_connect.d')

      PRINT*, ' '

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       grid_file = "square_grid_split.grd"
!       forcing_file = "square_grid_split.bfr"
!       dt = 1d-4
!       tf = .5d0

      grid_file = "../grids/inlet1.grd"
      forcing_file = "../grids/inlet1.bfr"
      dt = 1d0!/2d0      ! 1d0 for p=1,2, .5d0 for p=3
      tf = 1d0*86400d0
      dramp = .5d0
      cf = .003d0

!       grid_file = "../grids/inlet2.grd"
!       forcing_file = "../grids/inlet2.bfr"
!       dt = .5d0
!       tf = 1d0*86400d0
!       dramp = .5d0
!       cf = .003d0

!       grid_file = "../grids/converge.grd"
!       forcing_file = "../grids/converge.bfr"
!       dt = .5d0
!       tf = 86400d0
!       dramp = .5d0
!       cf = .003d0

!       grid_file = "../grids/converge3.grd"
!       forcing_file = "../grids/converge3.bfr"
!       dt = .125d0
!       tf = .2d0*86400d0
!       dramp = .08d0
!       cf = .0025d0

!       grid_file = "../grids/converge4.grd"
!       forcing_file = "../grids/converge4.bfr"
!       dt = .0625d0
!       tf = .2d0*86400d0
!       dramp = .08d0
!       cf = .0025d0

      p = 1
    
      nsp = 5
      nsp2 = 10

      tstep = int(tf/dt)
      lines = 100d0
      tskp = int(tf/(lines*dt)) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ndof = (p+1)*(p+2)/2

      ! Read in grid file
      CALL read_grid()

      ! Read in forcing file
      CALL read_forcing()
   
      ! Find edge connectivity
      CALL connect()
      
      ! Compute element area, edge length, edge normals, and bathymetry derivatives
      CALL element_data()

      ! Get quadrature points for area integration
      CALL area_qpts()

      ! Calculate basis function and derivative values at area quadrature points
      CALL area_basis()

      ! Calculate basis function values at edge quadrature points
      CALL edge_qpts()

      ! Calculate basis function and derivative values at edge quadrature points
      CALL edge_basis()

      ! Allocate arrays needed in time-stepping and rhs evaluation
      CALL alloc_arrays()

      CALL ptr_arrays()
      
      ! Set up netcdf output files
!       CALL file_setup()

      ! Compute initial condition, boundary forcing interpolation
      CALL initial()

      CALL CPU_TIME(t_start)

      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", "               Time Stepping                 "
      PRINT "(A)", "---------------------------------------------"
      PRINT "(A)", " "

      PRINT "(A,e12.4)", "Time step: ",dt
      PRINT "(A,e12.4)", "Final time: ",tf

      PRINT "(A)", " "

      t = 0d0
      cnt = 0
      DO it = 1,tstep

        CALL rk()

         t = t + dt

         cnt = cnt + 1
         IF(cnt == tskp) THEN

           PRINT("(A,e15.8)"), 't = ', t

           WRITE(63,"(e24.17)") t
           DO dof = 1,ndof
             WRITE(63,"(16000(e24.17,1x))") H(:,dof)
           ENDDO

           WRITE(641,"(e24.17)") t
           DO dof = 1,ndof
             WRITE(641,"(16000(e24.17,1x))") Qx(:,dof)
           ENDDO

           WRITE(642,"(e24.17)") t
           DO dof = 1,ndof
             WRITE(642,"(16000(e24.17,1x))") Qy(:,dof)
           ENDDO
             
           cnt = 0

         ENDIF

      ENDDO

      CALL CPU_TIME(t_end)

      PRINT*, ' '
      PRINT("(A,F10.5,A)"), "CPU time = ",t_end-t_start," seconds"

      PRINT*, ' '
      END PROGRAM swe_tri
