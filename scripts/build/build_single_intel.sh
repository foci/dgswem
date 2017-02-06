#!/bin/bash

# This script is to build the dgswem_cpp_serial exe
# with ifort used to build the Fortran modules
# and gnu compilers for the c++

BUILD_TYPE="Release"
BUILD_DIR="build_release_single_intel"

REPO_PATH=".."

# Machine specific

#CXX_COMPILER="/opt/mn/gcc/6.2.0/bin/g++"
#C_COMPILER="/opt/mn/gcc/6.2.0/bin/gcc"
#Fortran_COMPILER="/opt/mn/gcc/6.2.0/bin/gfortran"
CXX_COMPILER=g++
C_COMPILER=gcc
Fortran_COMPILER=ifort
#PREFIX_PATH="/home/zbyerly/local_install_debug"
#CMAKE="/opt/eb/software/cmake/3.6.2/bin/cmake"
CMAKE=cmake
echo "Build type: ${BUILD_TYPE}"
echo "Build directory: ${BUILD_DIR}"
echo "repo path = ${REPO_PATH}"

mkdir ${BUILD_DIR} || echo "${BUILD_DIR} already exists, continuing"
cd ${BUILD_DIR} || { 
    echo "could not cd into build directory ${BUILD_DIR}, quitting."
    exit
}

$CMAKE \
    -DNO_LGD=true \
    -DCMAKE_CXX_COMPILER=${CXX_COMPILER} \
    -DCMAKE_C_COMPILER=${C_COMPILER} \
    -DCMAKE_Fortran_COMPILER=${Fortran_COMPILER} \
    -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
    -DPROFILE=true \
    ${REPO_PATH}

