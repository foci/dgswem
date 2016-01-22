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

include_directories(../metis/Lib)
add_library(metis SHARED
../metis/Lib/metis.h
../metis/Lib/balance.c
../metis/Lib/bucketsort.c
../metis/Lib/ccgraph.c
../metis/Lib/coarsen.c
../metis/Lib/compress.c
../metis/Lib/debug.c
../metis/Lib/defs.h
../metis/Lib/estmem.c
../metis/Lib/fm.c
../metis/Lib/fortran.c
../metis/Lib/frename.c
../metis/Lib/graph.c
../metis/Lib/initpart.c
../metis/Lib/kmetis.c
../metis/Lib/kvmetis.c
../metis/Lib/kwayfm.c
../metis/Lib/kwayrefine.c
../metis/Lib/kwayvolfm.c
../metis/Lib/kwayvolrefine.c
../metis/Lib/macros.h
../metis/Lib/Makefile
../metis/Lib/match.c
../metis/Lib/mbalance2.c
../metis/Lib/mbalance.c
../metis/Lib/mcoarsen.c
../metis/Lib/memory.c
../metis/Lib/mesh.c
../metis/Lib/meshpart.c
../metis/Lib/mfm2.c
../metis/Lib/mfm.c
../metis/Lib/mincover.c
../metis/Lib/minitpart2.c
../metis/Lib/minitpart.c
../metis/Lib/mkmetis.c
../metis/Lib/mkwayfmh.c
../metis/Lib/mkwayrefine.c
../metis/Lib/mmatch.c
../metis/Lib/mmd.c
../metis/Lib/mpmetis.c
../metis/Lib/mrefine2.c
../metis/Lib/mrefine.c
../metis/Lib/mutil.c
../metis/Lib/myqsort.c
../metis/Lib/ometis.c
../metis/Lib/parmetis.c
../metis/Lib/pmetis.c
../metis/Lib/pqueue.c
../metis/Lib/proto.h
../metis/Lib/refine.c
../metis/Lib/rename.h
../metis/Lib/separator.c
../metis/Lib/sfm.c
../metis/Lib/srefine.c
../metis/Lib/stat.c
../metis/Lib/struct.h
../metis/Lib/subdomains.c
../metis/Lib/timing.c
../metis/Lib/util.c
)

add_executable(
  adcprep
  adcprep.F
  decomp.F
  read_global.F
  prep.F
  interp.F
  machdep.F
  metis.F
  pre_global.F
  presizes.F
)

target_link_libraries(adcprep metis)
