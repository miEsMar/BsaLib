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

set(SRC
   ${CMAKE_CURRENT_SOURCE_DIR}/src/BsaLib.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/CONSTANTS/CONSTANTS.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/IO/BsaIO.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/IO/BsaIOImpl.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/bsa/BsaLibImpl.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/bsa/classic/BsaClassicImpl.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/bsa/data/BsaLibData.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/bsa/data/BsaLibDataImpl.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/bsa/functions/functions.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/bsa/functions/functionsImpl.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/bsa/meshing/BsaMesherImpl.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/bsa/meshing/point/MPoint.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/bsa/meshing/policy/MPolicy.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/bsa/meshing/zones/M2DPolygZone.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/bsa/meshing/zones/MRectZone.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/bsa/meshing/zones/MRectZoneImpl.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/bsa/meshing/zones/MTriangZone.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/bsa/meshing/zones/MTriangZoneImpl.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/bsa/meshing/zones/MZone.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/settings/Settings.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/settings/SettingsImpl.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/structure/StructureData.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/structure/StructureImpl.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/timing/Timer.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/utils/utility.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/wind/WindPSDImpl.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/wind/WindSetImpl.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/wind/WindType.F90
)



# TODO: at some point, enable shared library generation
add_library( ${project} STATIC ${SRC} )
set_target_properties( ${project} PROPERTIES
   DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX}
)


if ( enable-pod-caching )
   target_compile_definitions( ${project} PRIVATE BSA_USE_POD_DATA_CACHING )
endif()


# BUG: using OMP in post-meshing phase will result in slower runtime 
#      due to critical sections synchronisation where accessing shared 
#      file I/O. For better performances, consider rethinking about the
#      logic related, otherwise better use 1 thread only.
if ( enable-post-mesh-omp )
   target_compile_definitions( ${project} PRIVATE BSA_USE_POST_MESH_OMP )
endif()


# Preprocessing flag
if (intel-compiler)
   add_compile_options( $<$<COMPILE_LANGUAGE:Fortran>:-fpp> )
else()
   add_compile_options( $<$<COMPILE_LANGUAGE:Fortran>:-cpp> )
endif()


# Link to MSVC runtime libraries.
# Ignored if compiler does not target MSVC!
# Taken from https://cmake.org/cmake/help/latest/prop_tgt/MSVC_RUNTIME_LIBRARY.html
# NOTE: if using shared, append DLL!
set_property(TARGET ${project} PROPERTY
   MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")


# TODO: see if we can manage this better
# NOTE: disabling annoying -O2 default release flag...
set( CMAKE_Fortran_FLAGS_RELEASE-default ${CMAKE_Fortran_FLAGS_RELEASE} PARENT_SCOPE )
set( CMAKE_Fortran_FLAGS_RELEASE "" )


if (intel-compiler AND enable-mkl)
   message( STATUS "NOTE: Using Intel MKL libraries." )
   if (WIN32)
      target_compile_options( ${project} INTERFACE "-Qmkl:parallel" )
   else()
      target_compile_options( ${project} INTERFACE "-qmkl=parallel" )
   endif()
   set( enable-openmp ON )
else()
   find_package( LAPACK REQUIRED )
   if (NOT LAPACK_FOUND )
      message( SEND_ERROR "Could not find LAPACK." )
   endif()
   target_link_libraries( ${project} INTERFACE ${LAPACK_LIBRARIES} )
endif()


if ( enable-openmp )
   # NOTE: find_package for OpenMP used to not work on Windows.
   #       See https://discourse.cmake.org/t/unique-way-to-enable-openmp-with-cmake/10350
   #       However, it now seems to be fixed using version 3.31.3.
   #       If using previous version where this still fails, then you need 
   #       to uncommend the if statement below.

   # # Uncomment this if find_package() fails on Windows.
   # if (WIN32)
   #    if (intel-compiler)
   #       target_compile_options( ${project} PRIVATE "-Qopenmp" )
   #    else()
   #       target_compile_options( ${project} PRIVATE "-fopenmp" )
   #    endif()
   # else()
      find_package( OpenMP REQUIRED COMPONENTS Fortran )
      target_link_libraries( ${project} PUBLIC OpenMP::OpenMP_Fortran )
   # endif()
else()
   if ( intel-compiler )
      target_compile_options( ${project} PRIVATE $<IF:${WIN32},-assume:nocc_omp,> )
      target_compile_options( ${project} PRIVATE $<IF:${WIN32},-Qparallel,-parallel> )
   endif()
endif()



if ( enable-shared-runtime )
   message( STATUS "NOTE: Using SHARED runtime libraries!" )
else()
   if (intel-compiler)
      if (WIN32)
         target_compile_options( ${project} PUBLIC "-libs:static" )
      else()
         target_compile_options( ${project} PUBLIC "-static" )
         target_compile_options( ${project} PUBLIC "-static-intel" )
         if (enable-openmp)
            target_compile_options( ${project} PUBLIC "-qopenmp-link=static" )
         endif()
      endif()
   else()
      message( WARNING "Check how to link to static runtime libraries for non MSVC!" )
   endif()
endif()


if (enable-sym-ev-routine)
   message( STATUS "NOTE: Using 'dsyev'  method" )
else()
   message( STATUS "NOTE: Using 'dgesvd' method" )
   target_compile_definitions( ${project} PRIVATE "-DBSA_USE_SVD_METHOD" )
endif()



if ( enable-gpu-code )

   if (enable-cuda)
      enable_language( CUDA )
      if (enable-ocl-ptx-code-gen)
         set( enable-ocl-ptx-code-gen OFF )
      endif()
   endif()
   add_subdirectory( BsaCL )

   target_compile_definitions( ${project} PRIVATE BSA_USE_GPU )
   if (enable-cuda)
      target_compile_definitions( ${project} PRIVATE BSA_USE_CUDA )
   endif()
   target_link_libraries( ${project} PRIVATE BsaCL_Fortran )
endif()


# Set module directory
set(LIB_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/modules)
if(NOT EXISTS "${LIB_MOD_DIR}")
  make_directory("${LIB_MOD_DIR}")
endif()
set_target_properties( ${project} PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_MOD_DIR} )

target_include_directories( ${project}
   PUBLIC
   $<BUILD_INTERFACE:${LIB_MOD_DIR}>
   $<INSTALL_INTERFACE:${CMAKE_INSTALL_MODULEDIR}>
)

