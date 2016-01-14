if (Fortran_COMPILER_NAME MATCHES "gfortran.*")  
  # FFLAGS1 DPRE IMODS O_DIR  
  set (CMAKE_Fortran_FLAGS_RELEASE   "-w -O3 -FI -DREAL8 -DLINUX -DGFORTRAN")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-w -g -DREAL8 -DLINUX -DGFORTRAN")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

add_executable(
  adcpost
  adcpost.F
  post.F
  compare.F
  diffmerge.F
  post_global.F
)