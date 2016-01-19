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

if(NOT NO_LGD)
  
  #LGD
  find_package(libgeodecomp REQUIRED)
  include_directories(${libgeodecomp_INCLUDE_DIRS})
  link_directories(${libgeodecomp_LIB_DIRS})

  #HPX
  find_package(HPX REQUIRED)
  include_directories(${HPX_INCLUDE_DIRS})
  link_directories(${HPX_LIBRARY_DIRS})
  message("HPX_INCLUDE_DIRS:" ${HPX_INCLUDE_DIRS})
  message("HPX_LIBRARY_DIR:" ${HPX_LIBRARY_DIR})
  

  #  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${libgeodecomp_FLAGS} ${HPX_FLAGS} -L/home/gentryx/local_install/lib -O3")

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${libgeodecomp_FLAGS} ${HPX_FLAGS}")
  
  #  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${libgeodecomp_FLAGS} -L/home/gentryx/local_install/lib -DVERBOSE")
  #message("libgeodecomp_INCLUDE_DIRS:" ${libgeodecomp_INCLUDE_DIRS})
  #message("libgeodecomp_FLAGS:" ${libgeodecomp_FLAGS})
endif(NOT NO_LGD)

#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

message("CMAKE_CXX_COMPILER:" ${CMAKE_CXX_COMPILER})
message("CMAKE_BUILD_TYPE:" ${CMAKE_BUILD_TYPE})

message("CMAKE_CXX_FLAGS:" ${CMAKE_CXX_FLAGS})
message("CMAKE_INCLUDE_DIRECTORIES:" ${CMAKE_INCLUDE_DIRECTORIES})


add_subdirectory(src)
add_subdirectory(prep)
add_subdirectory(error)
