/**
 * This file is part of BsaLib.
 * Copyright (C) 2024  Michele Esposito Marzino 
 *
 * BsaLib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BsaLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with BsaLib.  If not, see <https://www.gnu.org/licenses/>.
 * */
#ifndef BSACL_BASE_H_
#define BSACL_BASE_H_

/* This controls the type of incoming external data*/
#ifdef BSA_SINGLE_FLOATING_PRECISION
# define __real float
#else
# define __real double
#endif


/* This controls the type of internal data */
#ifdef BSACL_USE_DOUBLE_PRECISION
# define BSACL_REAL double
# define BSACL_REAL_MIN DBL_MIN
# ifdef BSACL_INCLUDE
#  ifdef cl_khr_fp64
#    pragma OPENCL EXTENSION cl_khr_fp64 : enable
#  elif defined(cl_amd_fp64)
#    pragma OPENCL EXTENSION cl_amd_fp64 : enable
#  else
#    error "Double precision floating point not supported by OpenCL implementation."
#  endif
#  define BSACL_REAL_IS_DOUBLE
# else
#  pragma message("   --- [NOTE]:  enabling GPU double (f64) floating point precision.")
# endif
#else
# define BSACL_REAL float
# define BSACL_REAL_MIN FLT_MIN
#endif


#ifdef BSACL_USE_CUDA__
# ifndef BSACL_KERNEL_ID
#  define BSACL_KERNEL_ID 2
# endif
#else
# if (BSACL_KERNEL_ID==2)
#  ifndef BSACL_PASS_PARAMS_BY_MACRO__
#   define BSACL_PASS_PARAMS_BY_MACRO__
#  endif
# endif
#endif


// NOTE: default define BSACL_WIpWG if not passed as argument when 
//       compiling on-the-fly this CL source.
#ifndef BSACL_WIpWG
#  define BSACL_WIpWG 256
#endif


#ifndef STRINGIFYMACRO_LITERAL
#  define STRINGIFYMACRO_LITERAL(X) #X
#  define xstr(s) STRINGIFYMACRO_LITERAL(s)
#  define STRINGIFYMACRO_VALUE(X) xstr(X)
#endif


#ifdef BSACL_USE_CONST_EXT_POINTERS
# define __EXT_PTR_CONST const
#else
# define __EXT_PTR_CONST
#endif


#ifdef BSACL_USE_RETRICTED_POINTERS
# ifdef __cplusplus
// from: https://stackoverflow.com/questions/5947564/whats-a-good-way-to-check-availability-of-restrict-keyword
#  if defined(__GNUC__) && ((__GNUC__ > 3) || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#    define __RESTRICT_PTR __restrict
#  elif defined(_MSC_VER) && _MSC_VER >= 1400
#    define __RESTRICT_PTR __restrict
#  else
#    define __RESTRICT_PTR
#  endif
# else // C compiler, is a keyword.
#  define __RESTRICT_PTR restrict
# endif
# define __RESTRICT_PTR_CL __restrict__
#else
# define __RESTRICT_PTR
# define __RESTRICT_PTR_CL
#endif


// Wind PSDs
#define BSACL_PSD_TYPE_UNKNOWN   0
#define BSACL_PSD_TYPE_VONKARMAN 1
#define BSACL_PSD_TYPE_KAIMAL    2
#define BSACL_PSD_TYPE_DAVENPORT 5

// BUG: control (user defined)
#define BSACL_CONV_PULSATION


#ifdef BSACL_USE_CUDA__

# define BSACL_PI (BSACL_REAL)M_PI
# define ierr_t cudaError_t
# define BSACL_SUCCESS cudaSuccess
# define BSACL_DEVICE_FREE_MEM(X) cudaFree((void *)X)
# define BSACL_INT int
# define BSACL_UINT unsigned int
# define BSACL_USHORT unsigned short
# define INT int
# define UINT unsigned int
# define USHORT unsigned short
# define BOOL BSACL_UINT
# define BSACL_MEM_CAST(x) (x)
# define BSACL_MEM void*
# define BSACL_MEM_UINT_T unsigned int*
# define BSACL_MEM_INT_T  int*
# define BSACL_MEM_REAL_T __real*
# define DEVICE __device__
# define PRIVATE
# define GLOBAL
# define LOCAL __shared__
# define CONSTANT
# define KERNEL __global__
# define EXTERN_C extern "C"
# define POW_F powf
# define POW_D pow
# define FABS_F fabsf
# define FABS_D fabs
# ifdef BSACL_USE_DOUBLE_PRECISION
#  define POW pow
#  define FABS fabs
# else
#  define POW powf
#  define FABS fabsf
# endif
# define LOCAL_ID_X_DIM0 threadIdx.x
# define LOCAL_ID_Y_DIM1 threadIdx.y
# define LOCAL_ID_Z_DIM2 threadIdx.z
# define N_BLOCKS_WGROUPS_X_DIM0 gridDim.x
# define N_BLOCKS_WGROUPS_Y_DIM1 gridDim.y
# define N_BLOCKS_WGROUPS_Z_DIM2 gridDim.z
# define BLOCK_ID_X_DIM0 blockIdx.x
# define BLOCK_ID_Y_DIM1 blockIdx.y
# define BLOCK_ID_Z_DIM2 blockIdx.z
# define GLOBAL_ID_X_DIM0 threadIdx.x + blockDim.x * blockIdx.x
# define GLOBAL_ID_Y_DIM1 threadIdx.y + blockDim.y * blockIdx.y
# define GLOBAL_ID_Z_DIM2 threadIdx.z + blockDim.z * blockIdx.z
# define LOCAL_WORKGROUP_BARRIER __syncthreads()

# define BSACL_MEM_READ_ONLY  0
# define BSACL_MEM_WRITE_ONLY 0
# define BSACL_MEM_READ_WRITE 0
# define BSACL_MEM_HOST_NO_ACCESS  0
# define BSACL_MEM_HOST_READ_ONLY  0
# define BSACL_MEM_HOST_WRITE_ONLY 0
# define BSACL_MEM_USE_HOST_PTR   0
# define BSACL_MEM_ALLOC_HOST_PTR 0
# define BSACL_MEM_COPY_HOST_PTR  0

#else  // OpenCL

# define BSACL_PI (BSACL_REAL)3.14159265358979323846
# define ierr_t cl_int
# define BSACL_SUCCESS CL_SUCCESS
# define BSACL_DEVICE_FREE_MEM(X) clReleaseMemObject(X)
# define BSACL_INT int
# define BSACL_UINT cl_uint
# define BSACL_USHORT ushort
# define INT int
# define UINT uint
# define USHORT ushort
# define BOOL bool
# define BSACL_MEM_CAST(x)
# define BSACL_MEM cl_mem
# define BSACL_MEM_UINT_T cl_mem
# define BSACL_MEM_INT_T  cl_mem
# define BSACL_MEM_REAL_T cl_mem
# define DEVICE
# define PRIVATE __private
# define GLOBAL __global
# define LOCAL __local
# define CONSTANT __constant
# define KERNEL __kernel
# define EXTERN_C
# define POW_F powr
# define POW_D powr
# define POW   powr
# define FABS_F fabs
# define FABS_D fabs
# define FABS   fabs
# define LOCAL_ID_X_DIM0 get_local_id(0)
# define LOCAL_ID_Y_DIM1 get_local_id(1)
# define LOCAL_ID_Z_DIM2 get_local_id(2)
# define N_BLOCKS_WGROUPS_X_DIM0 get_num_groups(0)
# define N_BLOCKS_WGROUPS_Y_DIM1 get_num_groups(1)
# define N_BLOCKS_WGROUPS_Z_DIM2 get_num_groups(2)
# define BLOCK_ID_X_DIM0  get_group_id(0)
# define BLOCK_ID_Y_DIM1  get_group_id(1)
# define BLOCK_ID_Z_DIM2  get_group_id(2)
# define GLOBAL_ID_X_DIM0 get_global_id(0)
# define GLOBAL_ID_Y_DIM1 get_global_id(1)
# define GLOBAL_ID_Z_DIM2 get_global_id(2)
# define LOCAL_WORKGROUP_BARRIER barrier(CLK_LOCAL_MEM_FENCE)

# define BSACL_MEM_READ_ONLY  CL_MEM_READ_ONLY
# define BSACL_MEM_WRITE_ONLY CL_MEM_WRITE_ONLY
# define BSACL_MEM_READ_WRITE CL_MEM_READ_WRITE
# define BSACL_MEM_HOST_NO_ACCESS  CL_MEM_HOST_NO_ACCESS
# define BSACL_MEM_HOST_READ_ONLY  CL_MEM_HOST_READ_ONLY
# define BSACL_MEM_HOST_WRITE_ONLY CL_MEM_HOST_WRITE_ONLY
# define BSACL_MEM_USE_HOST_PTR   CL_MEM_USE_HOST_PTR
# define BSACL_MEM_ALLOC_HOST_PTR CL_MEM_ALLOC_HOST_PTR
# define BSACL_MEM_COPY_HOST_PTR  CL_MEM_COPY_HOST_PTR
#endif


#define BSACL_ERROR_CODE(x) ((ierr_t)x)


#if (defined(BSACL_USE_CUDA__)) && (defined(BSACL_USE_FUSED_OP))
# pragma message("   --- [NOTE]  Using FUSED operations!")
# ifdef BSACL_USE_DOUBLE_PRECISION
#  define BSACL_FUSE_ADD(x, y) __dadd_rn((x), (y))
#  define BSACL_FUSE_MUL(x, y) __dmul_rn((x), (y))
#  define BSACL_FUSE_FMA(x, y, z) __fma_rn((x), (y), (z))
#  define BSACL_FUSE_RCP(x)    __drcp_rn((x))
#  define BSACL_FUSE_DIV(x, y) __ddiv_rn((x), (y))
#  define BSACL_FUSE_SQRT(x) __dsqrt_rn((x))
# else
#  define BSACL_FUSE_ADD(x, y) __fadd_rn((x), (y))
#  define BSACL_FUSE_MUL(x, y) __fmul_rn((x), (y))
#  define BSACL_FUSE_FMA(x, y, z) __fmaf_rn((x), (y), (z))
#  define BSACL_FUSE_RCP(x)    __frcp_rn((x))
#  define BSACL_FUSE_DIV(x, y) __fdiv_rn((x), (y))
#  define BSACL_FUSE_SQRT(x) __fsqrt_rn((x))
# endif
#else
# define BSACL_FUSE_ADD(x, y) ((x) + (y))
# define BSACL_FUSE_MUL(x, y) ((x) * (y))
# define BSACL_FUSE_FMA(x, y, z) ((x) * (y) + (z))
# define BSACL_FUSE_RCP(x)    (1.f / (x))
# define BSACL_FUSE_DIV(x, y) ((x) / (y))
# define BSACL_FUSE_SQRT(x) (sqrt(x))
#endif


#endif // BSACL_BASE_H_
