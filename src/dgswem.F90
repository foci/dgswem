!******************************************************************************
!     Dgswem RELEASE VERSION 1 11/2013   
!     
!******************************************************************************

PROGRAM DGSWEM
  USE GLOBAL
  USE DG
  USE NodalAttributes
  IMPLICIT NONE
  
  type (sizes_type), pointer :: s
  type (dg_type), pointer :: dg_here
  type (global_type), pointer :: global_here
  type (nodalattr_type), pointer :: nodalattr_here
  
  integer :: timestep,NT
  
  allocate(s)
  allocate(dg_here)
  allocate(global_here)
  allocate(nodalattr_here)
  
  !print*,"HPX defined"
  CALL DGSWEM_INIT(s,dg_here,global_here,nodalattr_here,NT)
  
  print*,"s%DIRNAME = ", s%DIRNAME

  DO timestep = global_here%ITHS+1,global_here%NT
     global_here%ITIME_A = timestep
     CALL DG_TIMESTEP(s,dg_here,global_here,nodalattr_here,timestep)
  END DO
  
  deallocate(s)
  deallocate(dg_here)
  deallocate(global_here)
  deallocate(nodalattr_here)
  
  
END PROGRAM DGSWEM
