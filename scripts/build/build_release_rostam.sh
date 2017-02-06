#!/bin/bash

# Currently, this build script assumes you want to build dgswem-hpx 
# on ROSTAM (rostam.cct.lsu.edu) in debug using gnu compilers for Fortran modules, 
# and using the pre-installed HPX and LGD libraries in Zach's home directory.
# but in the future this script will support building other configurations and on 
# other machines.

BUILD_TYPE="Release"
BUILD_DIR="build_release"

REPO_PATH=".."

# Machine specific

CXX_COMPILER="/opt/mn/gcc/6.2.0/bin/g++"
C_COMPILER="/opt/mn/gcc/6.2.0/bin/gcc"
Fortran_COMPILER="/opt/mn/gcc/6.2.0/bin/gfortran"
PREFIX_PATH="/home/zbyerly/local_install_rostam"
CMAKE="/opt/eb/software/cmake/3.6.2/bin/cmake"

echo "Build type: ${BUILD_TYPE}"
echo "Build directory: ${BUILD_DIR}"
echo "repo path = ${REPO_PATH}"

mkdir ${BUILD_DIR} || echo "${BUILD_DIR} already exists, continuing"
cd ${BUILD_DIR} || { 
    echo "could not cd into build directory ${BUILD_DIR}, quitting."
    exit
}

$CMAKE \
    -DCMAKE_PREFIX_PATH=${PREFIX_PATH} \
    -DCMAKE_CXX_COMPILER=${CXX_COMPILER} \
    -DCMAKE_C_COMPILER=${C_COMPILER} \
    -DCMAKE_Fortran_COMPILER=${Fortran_COMPILER} \
    -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
    ${REPO_PATH}
