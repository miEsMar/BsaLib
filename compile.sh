#!/bin/sh
#


# Check if we are in cluster
case "$HOSTNAME" in
    glogin*)
        module purge
        module load gcc/13.2.0
        module load cmake

        # For LAPACK
        module load EB/apps
        module load LAPACK/3.12.0-GCC-13.2.0
        ;;
    *)
        ;;
esac


# Setup
builddir="build"
if [ ! -d ./${builddir} ]; then mkdir ${builddir}; fi


# CLI loop
gen_args=""
bld_args=""
while [ "$1" != "" ]
do
    case "$1" in
        -c | --clean | --rebuild)
            bld_args="${bld_args}--clean-first "
            ;;
        *)
            gen_args="${gen_args}$1 "
    esac
    shift
done
echo ${gen_args}


# Generation
cmake \
    ${gen_args} \
    -D CMAKE_Fortran_COMPILER=gfortran \
    -D BUILD_SHARED_LIBS=OFF \
    -D enable-openmp=ON \
    -D enable-sym-ev-routine=OFF \
    -D enable-single=OFF \
    -D enable-gpu-code=OFF \
    -D enable-cuda=ON \
    -D enable-gpu-double=OFF \
    -S . -B ${builddir}
if [ "$?" != "0" ]; then exit 1; fi


# Build
cmake --build ${builddir} --config Debug ${bld_args}
if [ "$?" != "0" ]; then exit 1; fi

cmake --build ${builddir} --config Release ${bld_args}
if [ "$?" != "0" ]; then exit 1; fi


