#
# Copyright (c), Michele Esposito Marzino
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#      * Redistributions of source code must retain the above copyright notice, this
#        list of conditions and the following disclaimer.
#
#      * Redistributions in binary form must reproduce the above copyright notice,
#        this list of conditions and the following disclaimer in the documentation
#        and/or other materials provided with the distribution.
#
#      * Neither the name of the copyright holder nor the names of its
#        contributors may be used to endorse or promote products derived from
#        this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

cmake_minimum_required( VERSION 3.15 )


# # if CMAKE_BUILD_TYPE is specified use it; otherwise set the default build
# # type to "RelWithDebInfo" ("-O2 -g" with gcc) prior to calling project()
# if(DEFINED CMAKE_BUILD_TYPE)
#   set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose build type")
# else()
#    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose build type") # later, RelWithDebInfo
# endif()



###########################################################
# Main project name
###########################################################
file(STRINGS "${CMAKE_SOURCE_DIR}/VERSION" PROJECT_VERSION)
set( project BsaLib )
project( ${project} 
   LANGUAGES Fortran 
   VERSION ${PROJECT_VERSION} 
   DESCRIPTION "Library for Bispectral analysis of systems under non-Gaussian wind action"
)
include(CheckFortranSourceCompiles)
include(CheckFortranSourceRuns)


###########################################################
# Options
###########################################################
option(enable-sequential-code  "Disable any code parallelisation"                   OFF)
option(enable-single           "Enable single floating precision"                   OFF)
option(enable-heap-arrays      "Allocatable arrays are put on the heap"             ON)
option(enable-shared-runtime   "Link to shared Fortran runtime libraries"           OFF)
option(enable-openmp           "Enable OpenMP parallelisation"                      ON)
option(enable-mkl              "Enable Intel MKL libraries (Intel compilers only)"  ON)
option(enable-sym-ev-routine   "Enable usage of 'dsyev' instead of 'dgesvd'"        ON)
option(enable-safe-return      "Enable a safe return with error code instead of stopping execution" ON)

# experimental
option(enable-pod-caching     "Enable POD caching. EXPERIMENTAL!"                   OFF)
option(enable-post-mesh-omp   "Enable using OpenMP in post meshing. EXPERIMENTAL!"  OFF)

# GPU options
option(enable-gpu-code         "Enable GPU kernel code with OpenCL (defaulr) or CUDA backends"  OFF)
option(enable-const-ptr-def    "In GPU related code, declare internal pointers `const`"         ON)
option(enable-ocl-ptx-code-gen "Enable OpenCL PTX debug code generation"                        ON)
option(enable-cuda             "IF GPU code enabled, use CUDA instead of OpenCL (default)"      OFF)
option(enable-gpu-double       "Enable double precision in GPU kernel code (discouraged)"       OFF)


set( CMAKE_DEBUG_POSTFIX "d" )


###########################################################
# General section
###########################################################
if (enable-sequential-code)
   message( STATUS "Sequential code execution is enabled. Disabling parallelisation.")
   set( enable-openmp   OFF )
   set( enable-mkl      OFF )
   set( enable-gpu-code OFF )
endif()


# check if intel compiler or not
message( STATUS "NOTE: Compiler ID:  ${CMAKE_Fortran_COMPILER_ID}")
if (CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
   message( STATUS "NOTE: Using Intel compiler")
   set(intel-compiler ON)
else()
   set(intel-compiler OFF)
endif()


if ( NOT intel-compiler )
   set( enable-mkl OFF )
endif()




# Flags common to all configurations
if (WIN32)
   add_compile_options(
      $<$<COMPILE_LANGUAGE:Fortran>:-warn:all>
      $<$<COMPILE_LANGUAGE:Fortran>:-stand:f18>
   )
else()
   if (intel-compiler)
      add_compile_options(
         "$<$<COMPILE_LANGUAGE:Fortran>:-warn;all>"
         "$<$<COMPILE_LANGUAGE:Fortran>:-stand;f18>"
      )
   else()
      add_compile_options(
         $<$<COMPILE_LANGUAGE:Fortran>:-Wall>
         $<$<COMPILE_LANGUAGE:Fortran>:-Wimplicit-interface>
         $<$<COMPILE_LANGUAGE:Fortran>:-Wunused-parameter>
         $<$<COMPILE_LANGUAGE:Fortran>:-pedantic>
         $<$<COMPILE_LANGUAGE:Fortran>:-std=f2018>
      )
   endif()
endif()


if (enable-safe-return)
   add_compile_options(-DBSALIB_SAFE_RETURN)
endif ()


# Configuration specific options
if (enable-heap-arrays)
   if (intel-compiler)
      if (WIN32)
         add_compile_options( $<$<COMPILE_LANGUAGE:Fortran>:-heap-arrays0> )
      else()
         add_compile_options( "$<$<COMPILE_LANGUAGE:Fortran>:-heap-arrays;0>" )
      endif()
   else()
      if (NOT enable-openmp)
         add_compile_options( $<$<COMPILE_LANGUAGE:Fortran>:-fmax-stack-var-size=0> )
      endif()
   endif()
endif()

if (intel-compiler)
   add_compile_options( $<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANGUAGE:Fortran>>:-O3> )
   if (WIN32)
      add_compile_options( $<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANGUAGE:Fortran>>:-Ob1> )
      add_compile_options( $<$<COMPILE_LANGUAGE:Fortran>:-Qdiag-disable:8291> ) # field widths recommendation
      add_compile_options( $<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANGUAGE:Fortran>>:/traceback> )
   endif()
else()
   add_compile_options(
      $<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANGUAGE:Fortran>>:-Ofast>
   )
endif()


add_compile_definitions( $<$<CONFIG:DEBUG>:BSA_DEBUG> )


if (enable-single)
   message( STATUS "NOTE: enabling single precision floating point arithmetic")
   add_compile_definitions( BSA_SINGLE_FLOATING_PRECISION )
   set( enable-gpu-double OFF ) # make sure also GPU code is built in single
endif()


# Disable annoying warning of incompatible RT libraries
if (WIN32 AND enable-gpu-code)
   set( ignore_libs "MSVCRTD;MSVCRT;LIBCMT" )
   foreach( lib IN LISTS ignore_libs )
      add_link_options( "/NODEFAULTLIB:${lib}" )
   endforeach()
   unset( ignore_libs )
endif()



###########################################################
# Targets definition
###########################################################

add_custom_target( docs ford docs.md 
   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
   COMMENT "Building BsaLib documentation"
)

add_subdirectory( BsaLib )

# reset default flags, if ever changed (hopefully backed up!)
foreach(config ${CMAKE_CONFIGURATION_TYPES})
   string( TOUPPER ${config} config )
   if (DEFINED CMAKE_Fortran_FLAGS_${config}-default)
      set( CMAKE_Fortran_FLAGS_${config} ${CMAKE_Fortran_FLAGS_${config}-default} )
   endif()
endforeach()

add_subdirectory( bsa )


message( STATUS "NOTE: Build status:")
message( STATUS "     OMP: ${enable-openmp}")
message( STATUS "     MKL: ${enable-mkl}")
message( STATUS "     GPU: ${enable-gpu-code}")

