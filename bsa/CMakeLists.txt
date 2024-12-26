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
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.enable_language( Fortran )

add_executable( bsa bsa.F90 )

if (enable-openmp)
   message( STATUS "NOTE: setting bsa target name to 'bsa_omp'")
   set_target_properties( bsa PROPERTIES OUTPUT_NAME "bsa_omp" )
endif()

target_link_libraries( bsa BsaLib )

if (intel-compiler AND WIN32)
   target_compile_options( bsa PRIVATE -w )
endif()


if (GNU)
   target_compile_options( bsa <$<CONFIG:DEBUG>:-fbacktrace> )
endif()


if (enable-gpu-code)
   if (enable-cuda)
      set_target_properties(bsa PROPERTIES OUTPUT_NAME "bsa_cuda")
      set_property(TARGET bsa PROPERTY CUDA_SEPARABLE_COMPILATION ON)
      if(APPLE)
         # We need to add the path to the driver (libcuda.dylib) as an rpath,
         # so that the static cuda runtime can find it at runtime.
         set_property(TARGET bsa 
            PROPERTY
            BUILD_RPATH ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
      endif()
   else()
      set_target_properties(bsa PROPERTIES OUTPUT_NAME "bsa_cl")
   endif()
endif()


set( LIB_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/modules )
if(NOT EXISTS "${LIB_MOD_DIR}")
  make_directory("${LIB_MOD_DIR}")
endif()
set_target_properties( bsa
   PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_MOD_DIR}
)

target_include_directories( ${project} 
   PUBLIC
   $<BUILD_INTERFACE:${LIB_MOD_DIR}>
   $<INSTALL_INTERFACE:${CMAKE_INSTALL_MODULEDIR}>
)

unset( LIB_MOD_DIR )
