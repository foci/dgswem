#!/bin/bash

BUILD_TYPE="RelWithDebInfo"
BUILD_DIR="build_${BUILD_TYPE}_gnu"

REPO_PATH=".."

# Machine specific

CXX_COMPILER="/opt/mn/gcc/6.3.0/bin/g++"
C_COMPILER="/opt/mn/gcc/6.3.0/bin/gcc"
Fortran_COMPILER="/opt/mn/gcc/6.3.0/bin/gfortran"
#Fortran_COMPILER="/opt/eb/software/ifort/2013.5.192/bin/intel64/ifort"
PREFIX_PATH="/home/zbyerly/local_install_${BUILD_TYPE}"
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
    -DPROFILE="true" \
    ${REPO_PATH}