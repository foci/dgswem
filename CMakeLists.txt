cmake_minimum_required (VERSION 2.6)
project (dgswem_serial)
enable_language (Fortran)

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE DEBUG CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# default installation
get_filename_component (default_prefix "../install" ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in install."
      FORCE)

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-cpp -DREAL8 -DCSCA -DLINUX -DRKSSP -DSLOPE5 -ffixed-line-length-none -ffree-line-length-none -Waliasing -Wno-tabs -Wline-truncation -Wsurprising -fcheck=all -g -fbacktrace -fimplicit-none -I CMakeFiles/${PROJECT_NAME}.dir/src/ -I .")
#set (CMAKE_Fortran_FLAGS_DEBUG   "-O2 -cpp -DREAL8 -DCSCA -DLINUX -DRKSSP -DSLOPE5 -ffixed-line-length-none -ffree-line-length-none -Waliasing -Wno-tabs -Wline-truncation -Wsurprising -fcheck=all -g -fbacktrace -fimplicit-none")
  #elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  #  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  #  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
  #  message("Ifort not enabled now")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

#add_library ("modules" STATIC ../src/sizes.f95 ../src/wind.f95  ../src/dg.f95  ../src/global.f95  ../src/nodalattr.f95   ../src/read_fort_dg.f95   ../src/fort_dg_module.f95)

# Uncomment to remove -rdynamic
#SET(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS)

  
add_executable(
  dgswem_serial
  ../src/sizes.f95 
  ../src/global.f95  
  ../src/dg.f95  
  ../src/nodalattr.f95   
  ../src/fort_dg_module.f95
  ../src/wind.f95  
  ../src/dgswem.f95
  ../src/init_fileunits.f95
  ../src/read_input.f95
  ../src/cstart.f95
  ../src/hstart.f95
  ../src/radiation_edge_hydro.f95
  ../src/ibarrier_edge_hydro.f95
  ../src/ebarrier_edge_hydro.f95
  ../src/calc_normal.f95
  ../src/create_edge_data.f95
  ../src/DG_timestep.f95
  ../src/DG_hydro_timestep.f95
  ../src/flow_edge_hydro.f95
  ../src/internal_edge_hydro.f95
  ../src/land_edge_hydro.f95
  ../src/ocean_edge_hydro.f95
  ../src/orthobasis_area.f95
  ../src/orthobasis_edge.f95
  ../src/p_enrichment.f95
  ../src/prep_DG.f95
  ../src/quad_pts_area.f95
  ../src/quad_pts_edge.f95
  ../src/rhs_dg_hydro.f95
  ../src/numerical_flux.f95
  ../src/tidal_potential.f95
  ../src/write_results.f95
  ../src/slopelimiter.f95
  ../src/prep_slopelim.f95
  ../src/wetdry.f95
  ../src/LDG_hydro.f95
  ../src/met_forcing.f95
  ../src/sta_basis.f95
  ../src/make_dirname.f95
  ../src/nodalattr_subs.f95
  ../src/dgswem_init.f95
  ../src/cwrappers.f95
  ../src/read_fort_dg.f95
  )

#target_link_libraries(dgswem_serial modules)