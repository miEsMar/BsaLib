#!/bin/sh
#

. ./env

# Setup
builddir="build"
if [ ! -d ./${builddir} ]; then mkdir ${builddir}; fi


cmake \
    $* \
    -D CMAKE_Fortran_COMPILER=gfortran \
    -D BUILD_SHARED_LIBS=OFF \
    -D enable-openmp=ON \
    -D enable-sym-ev-routine=OFF \
    -D enable-single=OFF \
    -D enable-gpu-code=OFF \
    -D enable-cuda=ON \
    -D enable-gpu-double=OFF \
    -S . -B ${builddir}

