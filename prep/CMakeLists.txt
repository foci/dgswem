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

#METIS1_OBJ PREP_SRC -L ../metis -lmetis -lU77

#file(GLOB SOURCES RELATIVE "metis/Lib/*" )
#message("METIS_SOURCES:" ${SOURCES})
#add_library(metis SHARED ${SOURCES})

#add_executable(
#  adcprep
#  adcprep.F
#  decomp.F
#  read_global.F
#  prep.F
#  interp.F
#  machdep.F
#  metis.F
#  pre_global.F
#  presizes.F
#)