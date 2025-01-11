#!/bin/sh
#

. ./env

# Setup
builddir="build"
if [ ! -d ./${builddir} ]; then mkdir ${builddir}; fi

cmake \
    $* \
    ${lapack_libs} \
    -D CMAKE_Fortran_COMPILER=gfortran-12 \
    -D CMAKE_BUILD_TYPE=Debug \
    -D BUILD_SHARED_LIBS=OFF \
    -D enable-openmp=ON \
    -D enable-sym-ev-routine=OFF \
    -D enable-single=OFF \
    -D enable-gpu-code=OFF \
    -D enable-cuda=ON \
    -D enable-gpu-double=OFF \
    -D enable-pod-caching=OFF \
    -D enable-optimal-omp=ON \
    -S . -B ${builddir}

