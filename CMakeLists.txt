cmake_minimum_required (VERSION 2.6)
project (dgswem C CXX Fortran)
#enable_language (Fortran)

set(CMAKE_VERBOSE_MAKEFILE on)

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
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -DHPX")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-cpp -DHPX -DREAL8 -DCSCA -DLINUX -DRKSSP -DSLOPE5 -ffixed-line-length-none -ffree-line-length-none -Waliasing -Wno-tabs -Wline-truncation -Wsurprising -fcheck=all -g -fbacktrace -fimplicit-none -I CMakeFiles/${PROJECT_NAME}.dir/src/ -I .")
#set (CMAKE_Fortran_FLAGS_DEBUG   "-O2 -cpp -DREAL8 -DCSCA -DLINUX -DRKSSP -DSLOPE5 -ffixed-line-length-none -ffree-line-length-none -Waliasing -Wno-tabs -Wline-truncation -Wsurprising -fcheck=all -g -fbacktrace -fimplicit-none")
  #elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  #  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  #  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
  #  message("Ifort not enabled now")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

if(NOT NO_LGD)
  find_package(libgeodecomp REQUIRED)
  include_directories(${libgeodecomp_INCLUDE_DIRS})
  
  message("libgeodecomp_INCLUDE_DIRS:" ${libgeodecomp_INCLUDE_DIRS})

endif(NOT NO_LGD)

set (CMAKE_CXX_FLAGS "-std=c++11")

add_subdirectory(src)
