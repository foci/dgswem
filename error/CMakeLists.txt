if (Fortran_COMPILER_NAME MATCHES "gfortran.*")  
  # FFLAGS1 DPRE IMODS O_DIR  
  set (CMAKE_Fortran_FLAGS "-DREAL8 -DLINUX -DGFORTRAN -ffree-line-length-none -llapack")
  set (CMAKE_Fortran_FLAGS_RELEASE   "-O3 -FI")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-w -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

add_executable(
  error
  error.F90
  evaluate.F90
  read_solution.F90
  version.F90
  allocation.f90
  area_qpts.f90
  basis.f90
  globals.f90
  read_grid.f90
  read_input.f90
)
