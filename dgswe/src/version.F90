      MODULE version
      
      IMPLICIT NONE
      
      CHARACTER(:), ALLOCATABLE, PARAMETER :: gitBranch = "master" 
      CHARACTER(:), ALLOCATABLE, PARAMETER :: gitSHA = "2a0cbe96a8d0887ce123c10bf03b2b246c208e23 +" 
      CHARACTER(:), ALLOCATABLE, PARAMETER :: compiler_version = "mpifort for MPICH version 3.1.4 ifort version 14.0.0" 
      CHARACTER(:), ALLOCATABLE, PARAMETER :: compiler_flags = "-O3 -xHost -DCMPI -Iodir_dgswe_mpi/" 
      CHARACTER(:), ALLOCATABLE, PARAMETER :: modified_files = "../plot/plot_PE_grid.m ../plot/plot_grid.m ../plot/plot_stations.m ../plot/work/show_plots.py ../src/area_integration.F90 ../src/boundary_edge_elev.F90 ../src/boundary_edge_flow.F90 ../src/boundary_edge_land.F90 ../src/edge_integration.F90 ../src/linear_solve.F90 ../src/messenger2.F90 ../src/numerical_flux.F90" 
      CHARACTER(:), ALLOCATABLE, PARAMETER :: compile_date = "Mon Jan 16 23:18:11 EST 2017" 
      CHARACTER(:), ALLOCATABLE, PARAMETER :: host = "chl-tilos" 
            
      CONTAINS      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
            
      SUBROUTINE version_information(unit)

      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: unit
      INTEGER :: length,ind,sind,eind            

      
      length = LEN(TRIM(ADJUSTL(modified_files)))
      
            
      CALL directed_output("Version Information",unit)
      CALL directed_output("  Branch: "//TRIM(ADJUSTL(gitBranch)),unit)
      CALL directed_output("  SHA: "//TRIM(ADJUSTL(gitSHA)),unit)
      CALL directed_output("  Modified files: ",unit)

        
      sind = 1        
      ind = INDEX(modified_files," ")        
      eind = ind
      DO WHILE (ind > 0)
        CALL directed_output("    " //TRIM(ADJUSTL(modified_files(sind:eind))),unit)
        sind = eind+1
        ind = INDEX(modified_files(sind:length)," ")
        eind = eind + ind
      ENDDO
        
      CALL directed_output(" ",unit)
      CALL directed_output("Compiler Information",unit)
      CALL directed_output("  Version: "//TRIM(ADJUSTL(compiler_version)),unit)
      CALL directed_output("  Flags: "//TRIM(ADJUSTL(compiler_flags)),unit)
      CALL directed_output("  Date compiled: "//TRIM(ADJUSTL(compile_date)),unit)
      CALL directed_output("  Host: "//TRIM(ADJUSTL(host)),unit)
      CALL directed_output(" ",unit)
     
 
      END SUBROUTINE version_information
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
     

      
      END MODULE version