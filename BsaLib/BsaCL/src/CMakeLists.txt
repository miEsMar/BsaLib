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

# TODO: activate generation of shared library as well
add_library( BsaCL STATIC bsacl.c )

if (enable-gpu-double)
   target_compile_definitions( BsaCL PRIVATE BSACL_USE_DOUBLE_PRECISION )
endif()


if (enable-ocl-ptx-code-gen)
   message( STATUS "NOTE: Compiling PTX code ON")
   target_compile_definitions( BsaCL PRIVATE BSACL_GENERATE_PTX_BINARY )
endif()


if (enable-const-ptr-def)
   target_compile_definitions( BsaCL PRIVATE BSACL_USE_CONST_EXT_POINTERS )
endif()

target_compile_definitions( BsaCL PRIVATE BSACL_USE_RETRICTED_POINTERS )
target_compile_definitions( BsaCL PRIVATE BSACL_USE_FUSED_OP )


if (enable-cuda)
   set_target_properties( BsaCL PROPERTIES CUDA_ARCHITECTURES native )
   set_source_files_properties( bsacl.c PROPERTIES LANGUAGE CUDA )
   target_compile_options( BsaCL PRIVATE -diag-suppress 546 )
   target_compile_options( BsaCL PRIVATE
      $<$<CONFIG:Release>: --use_fast_math --extra-device-vectorization>
   )
   target_compile_options( BsaCL PRIVATE
      $<$<CONFIG:Debug>: -lineinfo>
   )
   target_compile_definitions( BsaCL PRIVATE BSACL_USE_CUDA__ )
   set_target_properties( BsaCL PROPERTIES POSITION_INDEPENDENT_CODE ON )
else()
   target_compile_definitions( BsaCL PRIVATE "BASE_DIRECTORY=\"${CMAKE_CURRENT_SOURCE_DIR}/\"" )
   target_link_libraries( BsaCL OpenCL::OpenCL )
endif()

target_include_directories( BsaCL BEFORE PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../include" )

if (MSVC AND NOT enable-shared-runtime)
   set_property(TARGET BsaCL PROPERTY
      MSVC_RUNTIME_LIBRARY "MultiThreadedDebug")
endif()
