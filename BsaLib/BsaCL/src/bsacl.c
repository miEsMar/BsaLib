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

// #ifndef BSACL_OPT_N_WORK_GROUPS
// #  define BSACL_OPT_N_WORK_GROUPS 10
// #endif

// #ifndef MIN_N_OF_WORK_GROUPS
// #  define MIN_N_OF_WORK_GROUPS 2
// #endif

// #ifndef BSACL_MAX_GPU_COUNT
// #  define BSACL_MAX_GPU_COUNT 8
// #endif


#ifndef BASE_DIRECTORY
# define BASE_DIRECTORY "./"  // By default, current binary directory.
#endif

#ifdef BSACL_ENABLE_EVALFCT_PTR
typedef void (*evalFct_t)(int, int, const double*, int, double*);
#endif

#ifdef BSACL_USE_CUDA__
# define _USE_MATH_DEFINES
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "_msgtypes.h"
#include "_errtypes.h"
#include "_base.h"
#include "bsacl.h"

#ifdef BSACL_USE_CUDA__

# ifdef BSACL_PASS_PARAMS_BY_MACRO__
#  undef BSACL_PASS_PARAMS_BY_MACRO__  // in CUDA, we pass always by kernel param!
# endif
// # define __CUDACC_RTC__
# pragma message("   --- [NOTE]:  Using  CUDA  in place of  OpenCL  specification!")
# include <limits>
# include "bsacl.cl"
# ifdef NULL // CUDA does not like (void *) ptr assignment
#   undef  NULL
#   define NULL 0
# endif

#else // using OpenCL

# include "_devtypes.h"
# define CL_USE_DEPRECATED_OPENCL_1_2_APIS
# define __CL_ENABLE_EXCEPTIONS
# include <CL/cl.h>

#endif



#define ABORT_INTERNAL_RETURN_VOID(X, Y)  {abortInternal_(X, Y); return;}
#define ABORT_INTERNAL_GOTO_RET_(X, Y, Z) {abortInternal_(X, Y); goto Z;}
#define ABORT_INTERNAL_RETURN_IERR(X, Y)  {abortInternal_(X, Y); return (X);}



//-------------------------------------
// GLOBAL BSACL VARIABLES
//-------------------------------------
static unsigned int n_threads__ = 0U;

static unsigned int bsacl_init_called__        = 0U;
static unsigned int bsacl_init_devmem_called__ = 0U;


#ifndef BSACL_USE_CUDA__


static cl_uint         n_platforms__  = 0;
static cl_platform_id *platform_all__ = NULL; // array of all available platforms
static cl_platform_id  platform__     = NULL; // final selected platform

static cl_device_type device_type__ = CL_DEVICE_TYPE_DEFAULT;
static cl_uint        n_devices__   = 0;       // n. devices for finally selected platform
static cl_device_id   *devices__    = NULL;    // array of all available devices for selected platform

/* Max memory allocation (bytes) allowed for one single buffer/image */
static cl_ulong *dev_max_mem_alloc_size__  = NULL;
static cl_ulong *dev_glob_mem_size__       = NULL;
static cl_ulong *dev_glob_mem_cache_size__ = NULL;
static cl_ulong *dev_local_mem_size__      = NULL;
static cl_uint  *dev_n_compute_units__     = NULL;
static size_t   *dev_max_work_group_size__ = NULL;
static cl_uint  *dev_max_work_items_n_dims__ = NULL;
static size_t  **dev_max_work_item_sizes__   = NULL;

static cl_context context__        = NULL;
static cl_command_queue *cqueues__ = NULL; /* Array of devices' command queues */
static cl_program *programs__      = NULL;
static cl_kernel  *kernels__       = NULL;
static unsigned int *thread_set_shared_kernel_args_ = NULL;
static char       **clSourceStrings__ = NULL;
static char       **evalFctStrings__  = NULL;
static char       **thread_prog_opts__= NULL;
static char       prog_compile_opts__[BUFSIZ] = {' '};

#endif // BSACL_USE_CUDA__




static unsigned short has_halted__  = 0U;
static unsigned short has_cleaned__ = 0U;

static unsigned short kernel_id_ = 2U;
#ifndef BSACL_USE_CUDA__
static unsigned short pass_params_by_macro_ = 1U;
#endif


// BSACL device memory buffers
BSACL_MEM_UINT_T d_tc__          = NULL;
BSACL_MEM_UINT_T d_nodes_load__  = NULL;
BSACL_MEM_REAL_T d_phiTc__       = NULL;
BSACL_MEM_REAL_T d_nod_corr__    = NULL;
BSACL_MEM_REAL_T d_wind_nodal_vel__   = NULL;
BSACL_MEM_REAL_T d_wind_turb_scales__ = NULL;
BSACL_MEM_REAL_T d_wind_turb_std__    = NULL;
BSACL_MEM_INT_T  d_wind_nodal_windz__ = NULL;
BSACL_MEM_REAL_T d_Mg__ = NULL;
BSACL_MEM_REAL_T d_Cg__ = NULL;
BSACL_MEM_REAL_T d_Kg__ = NULL;

static unsigned  *nfi__, *nfj__;
static __real    **fi__ = NULL;
static __real    **fj__ = NULL;
static BSACL_MEM_REAL_T *d_fi__ = NULL;
static BSACL_MEM_REAL_T *d_fj__ = NULL;

// ------------------------------------------------------------------
// EXTERNAL DATA (to be acquired, i.e. no internal memory allocation)
// ------------------------------------------------------------------

typedef struct extdata_t {

   int      *nodes_load__;  // list of loaded nodes in the structural model
   unsigned NNODES_LOAD__;
   unsigned NN__;

   unsigned PSD_ID__;

   unsigned NLIBS__;
   unsigned NDOFS__;

   unsigned NMODES_EFF__;
   unsigned NDEGW__;

   unsigned NTC__;
   unsigned NNOD_CORR__; // n. of nodal correlation entries (per direction!)

   unsigned NMODES__;

   unsigned NWZ__;  // N. of wind zones

   int  *tc__;           // list of effective turbulent components
   int  *modes_eff__;    // list of effective modes used

   __real *modmat__;
   __real *natfreqs__;

   __real *Mg__;
   __real *Cg__;
   __real *Kg__;

   __real *wfc__;       // wind force coefficients
   __real *phi_T_c__;   // Phi x C matrix
   __real *nod_corr__;  // nodal correlation array (no duplicates)

   __real *wind_nodal_vel__;
   __real *wind_turb_scales__;
   __real *wind_turb_std__;
   int    *wind_nodal_windz__;
} extdata_t;


/**
 * @brief Releases references to external memory, to avoid 
 *        unwanted deallocations and so memory corruption.
 * */
static inline void freeExtData(extdata_t *extdata) {
   extdata->nodes_load__       = NULL;
   extdata->modmat__           = NULL;
   extdata->natfreqs__         = NULL;
   extdata->modes_eff__        = NULL;
   extdata->wfc__              = NULL;
   extdata->nod_corr__         = NULL;
   extdata->phi_T_c__          = NULL;
   extdata->wind_nodal_vel__   = NULL;
   extdata->wind_nodal_windz__ = NULL;
   extdata->wind_turb_scales__ = NULL;
   extdata->wind_turb_std__    = NULL;
   extdata->NWZ__         = 0;
   extdata->NNODES_LOAD__ = 0;
   extdata->NLIBS__       = 0;
   extdata->NDOFS__       = 0;
   extdata->NMODES__      = 0;
   extdata->NMODES_EFF__  = 0;
   extdata->NDEGW__       = 0;
   extdata->NTC__         = 0;
   extdata->NNOD_CORR__   = 0;
   extdata->PSD_ID__      = 0;
}

extdata_t extdata__;

// BUG: Deprecated
static __real *S_uvw__       = NULL;
static __real *S_uvw_fiPfj__ = NULL;






// -----------------------------------
//  LOCAL helper functions
// -----------------------------------

#ifdef __cplusplus
extern "C" {
#endif


/**
 * @brief (internal) releases BSACL local device memory.
 *
 */
static inline ierr_t releaseDeviceMem_() {
   unsigned ierr_ = 0u;

   ierr_ |= BSACL_DEVICE_FREE_MEM(d_nodes_load__);
   ierr_ |= BSACL_DEVICE_FREE_MEM(d_nod_corr__);
   ierr_ |= BSACL_DEVICE_FREE_MEM(d_phiTc__);
   ierr_ |= BSACL_DEVICE_FREE_MEM(d_tc__);

   ierr_ |= BSACL_DEVICE_FREE_MEM(d_wind_nodal_vel__);
   ierr_ |= BSACL_DEVICE_FREE_MEM(d_wind_nodal_windz__);
   ierr_ |= BSACL_DEVICE_FREE_MEM(d_wind_turb_scales__);
   ierr_ |= BSACL_DEVICE_FREE_MEM(d_wind_turb_std__);

   ierr_ |= BSACL_DEVICE_FREE_MEM(d_Mg__);
   ierr_ |= BSACL_DEVICE_FREE_MEM(d_Cg__);
   ierr_ |= BSACL_DEVICE_FREE_MEM(d_Kg__);

   return (ierr_t)ierr_;
}



/**
 * @brief (internal) releases BSACL local memory.
 *
 */
static inline void freeMem_(void) {

   if (has_cleaned__) return;

#ifndef BSACL_USE_CUDA__

   if (clSourceStrings__ != NULL) {
      free(*clSourceStrings__);
      clSourceStrings__ = NULL;
   }

   if (platform_all__ != NULL) free(platform_all__);
   platform_all__ = NULL;
   platform__     = NULL;

   if (devices__ != NULL) free(devices__);
   devices__ = NULL;

   if (dev_max_mem_alloc_size__ != NULL) free(dev_max_mem_alloc_size__);
   dev_max_mem_alloc_size__ = NULL;
   if (dev_glob_mem_size__ != NULL)      free(dev_glob_mem_size__);
   dev_glob_mem_size__ = NULL;
   if (dev_glob_mem_cache_size__ != NULL) free(dev_glob_mem_cache_size__);
   dev_glob_mem_cache_size__ = NULL;
   if (dev_local_mem_size__ != NULL)      free(dev_local_mem_size__);
   dev_local_mem_size__ = NULL;
   if (dev_n_compute_units__ != NULL)     free(dev_n_compute_units__);
   dev_n_compute_units__ = NULL;
   if (dev_max_work_group_size__ != NULL)   free(dev_max_work_group_size__);
   dev_max_work_group_size__ = NULL;
   if (dev_max_work_items_n_dims__ != NULL) free(dev_max_work_items_n_dims__);
   dev_max_work_items_n_dims__ = NULL;

   if (dev_max_work_item_sizes__ != NULL) {
      for (unsigned i = 0; i < n_devices__; ++i) {
         if (dev_max_work_item_sizes__[i] != NULL) free(dev_max_work_item_sizes__[i]);
      }
      free(dev_max_work_item_sizes__);
      n_devices__ = 0;
   }

   if (cqueues__ != NULL) free(cqueues__);
   cqueues__ = NULL;
#endif

   freeExtData(&extdata__);

#ifdef BSACL_ENABLE_EVALFCT_PTR
#if (BSACL_KERNEL_ID==0)
   evalFunc_ptr__ = NULL;
#endif
#endif


   for (unsigned i = 0; i < n_threads__; ++i) {

#ifndef BSACL_USE_CUDA__
      free(thread_prog_opts__[i]);
#endif

      fi__[i] = NULL;
      fj__[i] = NULL;
      (void) BSACL_DEVICE_FREE_MEM(d_fi__[i]);
      (void) BSACL_DEVICE_FREE_MEM(d_fj__[i]);
   }
#ifndef BSACL_USE_CUDA__
   thread_prog_opts__ = NULL;

   free(thread_set_shared_kernel_args_);
   thread_set_shared_kernel_args_ = NULL;
#endif

   free(fi__);
   free(fj__);
   free(nfi__);
   free(nfj__);
   free(d_fi__);
   free(d_fj__);
   fi__   = NULL;
   fj__   = NULL;
   nfi__  = NULL;
   nfj__  = NULL;
   d_fi__ = NULL;
   d_fj__ = NULL;


   has_cleaned__       = 1U;
   bsacl_init_called__ = 0U;
}



static inline void printMsg_(const char *const msgtype, const char *const msg) {
   printf("%s%s", msgtype, msg);
}
static inline void printMsgWithIerr_(const char *const msgtype, const char *const msg, const int ierr) {
   printMsg_(msgtype, msg);
   if (ierr != 0) printf(" Aborting  (%d)", ierr);
   printf("\n");
}
static void abortInternal_(const int ierr, const char *const emsg) {
   if (emsg != NULL && ierr != 0) printMsgWithIerr_(ERRR_MSG, emsg, ierr);
   freeMem_();
   if (ierr != 0) has_halted__ = 1U;
}





#ifndef BSACL_USE_CUDA__
/**
 * @brief Get the all OpenCL Platforms available
 * 
 * @param npltfms pointer to the total n. of platforms found
 * @return BSACL_INT : error code
 */
BSACL_INT getNumPlatforms_(cl_uint *const npltfms) {
   BSACL_INT ierr_;

#ifdef BSA_DEBUG
   printMsg_(INFO_MSG, "Init querying for n. of available platforms...");
#endif

   ierr_ = clGetPlatformIDs(0, NULL, npltfms);
   if (*npltfms <= 0) {
      printMsg_(ERRR_MSG, "N<=0. OpenCL platform could not be correctly found.");
      return BSACL_PLATFORM_FIND_ERROR_;
   }
#ifdef BSA_DEBUG
   printf("\n%sFound  %d  platforms.", INFO_MSG, *npltfms);
#endif
   return ierr_;
}




BSACL_INT getAllPlatforms_() {
   BSACL_INT ierr_;

   ierr_ = getNumPlatforms_(&n_platforms__);
   if (ierr_ != BSACL_SUCCESS) {
      printf("\n%sCannot inquire n. of available CL platforms. Aborting.", ERRR_MSG);
      return ierr_;
   }

   platform_all__ = (cl_platform_id *) malloc(sizeof(cl_platform_id) * n_platforms__);
   ierr_          = clGetPlatformIDs(n_platforms__, platform_all__, NULL);
   if (ierr_ != BSACL_SUCCESS) {
      printf("\n%sCannot acquire  %d  available platforms. Aborting.", ERRR_MSG, n_platforms__);
      return ierr_;
   }
#ifdef BSA_DEBUG
   printf("\n%s%d  platforms acquired -- ok.", INFO_MSG, n_platforms__);
#endif
   return ierr_;
}



/**
 * @brief Select a platform by given platform info
 * 
 * @param selected_platform pointer to the actual platform criteria
 * @param pinfo platform info to query
 * @param pinfo_val value of info to query (if needed)
 * @param pinfo_nbytes n bytes of pinfo value (if any)
 * @return BSACL_INT 
 */
BSACL_INT getCLPlatformFromInfo_(
      cl_platform_id *selected_platform, cl_platform_info pinfo,
      const char * const pinfo_val,      size_t pinfo_nbytes)
{
   BSACL_INT ierr_;
   size_t param_size_;
   char buf_[BUFSIZ];

   *selected_platform = NULL;  // default initialise final platform to null

   if (platform_all__ == NULL) {
      ierr_ = getAllPlatforms_();
      if (platform_all__ == NULL) ierr_ = (BSACL_INT)-1;
      if (ierr_ != BSACL_SUCCESS) return ierr_;
   }

   if (pinfo_val != NULL) { // Search for info requested (all char[] !!)

      for (unsigned int i = 0; i < n_platforms__; i++) {

         ierr_ = clGetPlatformInfo(
            platform_all__[i], pinfo, BUFSIZ, &buf_, &param_size_);
         if (ierr_ == BSACL_SUCCESS) {
            if (strstr(buf_, pinfo_val) != NULL) {
               printf("\n%sPlatform with info  %s  found.", NOTE_MSG, pinfo_val);
               *selected_platform = platform_all__[i];
               break;
            }
         }
      }
   }

   // If info NOT found, or val not provided
   // Acquire first by default
   if (*selected_platform == NULL) {

      *selected_platform = platform_all__[0];

      if (*selected_platform == NULL) {
         printMsg_(ERRR_MSG, "Cannot even acquire first default CL platform. Check.");
         return BSACL_DEF_PLATFORM_ACQUIRE_ERROR_;
      }
      ierr_ = clGetPlatformInfo(
         *selected_platform, CL_PLATFORM_NAME, 1024, &buf_, &param_size_);

      buf_[param_size_] = '\0'; // BUG: make sure string is NULL terminated
      printf("\n%sPlatform with info   %s   NOT FOUND.", WARN_MSG, pinfo_val);
      printf("\n%sSelecting first by default (name: %s)", CONT_MSG, buf_);
   }
   printf("\n");
   return BSACL_SUCCESS;
};



/**
 * @brief Get the Selected Devices critical info needed for optimal setting.
 * 
 * @param  devices  list of selected available devices
 * @return (cl_int) Error code 
 */
void getCLDevicesInfo_() {
   size_t sz_, info_len_;

   // Allocate memory
   sz_ = sizeof(cl_ulong) * n_devices__;
   dev_max_mem_alloc_size__  = (cl_ulong *) malloc(sz_);
   dev_glob_mem_size__       = (cl_ulong *) malloc(sz_);
   dev_glob_mem_cache_size__ = (cl_ulong *) malloc(sz_);
   dev_local_mem_size__      = (cl_ulong *) malloc(sz_);

   sz_ = sizeof(cl_uint)  * n_devices__;
   dev_n_compute_units__       = (cl_uint *) malloc(sz_);
   dev_max_work_items_n_dims__ = (cl_uint *) malloc(sz_);
   sz_ = sizeof(size_t)  * n_devices__;
   dev_max_work_group_size__   = (size_t *) malloc(sz_);

   sz_ = sizeof(size_t*) * n_devices__;
   dev_max_work_item_sizes__ = (size_t **) malloc(sz_);

   for (unsigned i = 0; i < n_devices__; ++i) {

#ifdef BSA_DEBUG
      printf("%sRelevant info for device   #%d\n", INFO_MSG, i+1);
#endif

      clGetDeviceInfo(*(devices__ + i), CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), (void*)&(dev_glob_mem_size__[i]), 
         &info_len_);
#ifdef BSA_DEBUG
      printf("%s - glob mem size       =   %llu\n", CONT_MSG, dev_glob_mem_size__[i]);
#endif

      clGetDeviceInfo(*(devices__ + i), CL_DEVICE_GLOBAL_MEM_CACHE_SIZE, sizeof(cl_ulong), 
         (void*)&dev_glob_mem_cache_size__[i], &info_len_);
#ifdef BSA_DEBUG
      printf("%s - glob mem cache size =   %llu\n", CONT_MSG, dev_glob_mem_cache_size__[i]);
#endif

      clGetDeviceInfo(*(devices__ + i), CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), 
         (void*)&dev_local_mem_size__[i], &info_len_);
#ifdef BSA_DEBUG
      printf("%s - local mem size      =   %llu\n", CONT_MSG, dev_local_mem_size__[i]);
#endif

      clGetDeviceInfo(*(devices__ + i), CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(cl_ulong), 
         (void*)&dev_max_mem_alloc_size__[i], &info_len_);
#ifdef BSA_DEBUG
      printf("%s - max mem alloc size  =   %llu\n", CONT_MSG, dev_max_mem_alloc_size__[i]);
#endif

      clGetDeviceInfo(*(devices__ + i), CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), 
         (void*)&dev_n_compute_units__[i], &info_len_);
#ifdef BSA_DEBUG
      printf("%s - max compute units   =   %u\n", CONT_MSG, dev_n_compute_units__[i]);
#endif

      clGetDeviceInfo(*(devices__ + i), CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), 
         (void*)&dev_max_work_group_size__[i], &info_len_);
#ifdef BSA_DEBUG
      printf("%s - max work-group size =   %llu\n", CONT_MSG, dev_max_work_group_size__[i]);
#endif

      clGetDeviceInfo(*(devices__ + i), CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(cl_uint), 
         (void*)&dev_max_work_items_n_dims__[i], &info_len_);
#ifdef BSA_DEBUG
      printf("%s - max work-item dims  =   %u\n", CONT_MSG, dev_max_work_items_n_dims__[i]);
#endif

      sz_ = sizeof(size_t) * dev_max_work_items_n_dims__[i];
      dev_max_work_item_sizes__[i] = (size_t *) malloc(sz_);
      memset(dev_max_work_item_sizes__[i], (size_t)0, sz_);
      clGetDeviceInfo(*(devices__ + i), CL_DEVICE_MAX_WORK_ITEM_SIZES, 
         sizeof(size_t) * dev_max_work_items_n_dims__[i], 
         (void*)dev_max_work_item_sizes__[i], &info_len_);
#ifdef BSA_DEBUG
      printf("%s - max work-item sizes =   ", CONT_MSG);
      for (unsigned int j = 0; j < dev_max_work_items_n_dims__[i] - 1; ++j)
         printf("%llu - ", dev_max_work_item_sizes__[i][j]);
      printf("%llu\n", dev_max_work_item_sizes__[i][dev_max_work_items_n_dims__[i] - 1]);
#endif
   } // n devices
}



/**
 * @brief Get the Devices By Type object
 * 
 * @param platform Chosen platform
 * @param device_type Type of desired devices to query for.
 * @param devices  Stores the actual devices (of type) found for platform
 * @return (cl_int) Number of found devices (of specified type)
 */
cl_int getCLDevicesByType_(
      cl_platform_id *platform, cl_device_type device_type, cl_device_id **devices) {

   cl_int  ierr_ = 0;
   cl_uint num_devices_of_type_ = 0;

#ifdef BSA_DEBUG
   printMsg_(INFO_MSG, "Init querying devices by type..\n");
#endif

   // Get all available devices of type, for given platform
   ierr_ = clGetDeviceIDs(*platform, device_type, 0, NULL, &num_devices_of_type_);
   if (ierr_ != BSACL_SUCCESS) {
      printf(
         "\n\n%sCannot query requested devices type for provided platform. Aborting.", ERRR_MSG);
      return -1;
   }
   *devices = (cl_device_id *)malloc(sizeof(cl_device_id) * num_devices_of_type_);
   if (*devices == NULL) {
      printf("\n\n%sFailed to allocate memory for type of devices. Aborting.", ERRR_MSG);
      return -2;
   }

   // get actual devices by given type
   ierr_ = clGetDeviceIDs(*platform, device_type, num_devices_of_type_, *devices, NULL);
   if (ierr_ != BSACL_SUCCESS || devices == NULL) {
      printf("\n\n%sFailed to acquire type of devices for platform. Aborting (%d).", ERRR_MSG, ierr_);
      return -3;
   }
#ifdef BSA_DEBUG
   printf("%sFound  %d  device(s) of given type.\n", INFO_MSG, num_devices_of_type_);
#endif
   return num_devices_of_type_;
}



size_t readFileContent_(char ** buf_, const char * fname_)
{
   FILE *fp;
   errno_t nbytes_;
   long fsize_ = 0L;
   char ffp_[256];

   fsize_ = (long)(strlen(BASE_DIRECTORY) + strlen(fname_));
   strcpy(ffp_, BASE_DIRECTORY);
   nbytes_ = strcat_s(ffp_, fsize_ + 1, fname_);

   nbytes_ = fopen_s(&fp, ffp_, (const char *)"r");
   if (nbytes_ != 0) return 0;

   rewind(fp);
   fseek(fp, 0L, SEEK_END);
   fsize_ = ftell(fp);
   rewind(fp);

   *buf_ = (char *)malloc(sizeof(char) * (fsize_+1));
   if (*buf_ == NULL) return 0;

   nbytes_ = (int)fread_s(*buf_, fsize_, 1, fsize_, fp);
   if (nbytes_ == 0LL) {
      free(*buf_);
      ABORT_INTERNAL_RETURN_IERR(-431, "Failed to read source file.");
   }
   (*buf_)[nbytes_] = '\0'; // NULL terminate string

   return (size_t)nbytes_;
}




size_t assembleFinalCLSource_(char **clBaseStrings, char **evalFctStrings)
{
   size_t fsize_ = 0, clsize_ = 0, evalsize_ = 0, fpos_ = 0;

   if (clBaseStrings != NULL)  clsize_   = strlen(*clBaseStrings);
   if (evalFctStrings != NULL) evalsize_ = strlen(*evalFctStrings);

   fsize_ = clsize_ + evalsize_;
   if (fsize_ == 0) return 0;

   fsize_ += 4;
   clSourceStrings__  = (char **)malloc(sizeof(char *));
   *clSourceStrings__ = (char  *)malloc(sizeof(char) * fsize_);
   strcpy(*clSourceStrings__, *clBaseStrings);
   fpos_ = clsize_;
   if (evalFctStrings && *evalFctStrings) {
      *(*clSourceStrings__ + fpos_) = '\n';
      fpos_++;
      *(*clSourceStrings__ + fpos_) = '\n';
      strcpy(*clSourceStrings__ + ++fpos_, *evalFctStrings);
      fpos_ += evalsize_+1;
   }
   *(*clSourceStrings__ + fpos_) = '\0';  // null terminate

   return fpos_;
}




void assembleProgramBuildOptsString_(const unsigned i_thread)
{
   char *buf  = prog_compile_opts__;
   char *tbuf = thread_prog_opts__[i_thread];

   strcpy(buf, "-D BSACL_INCLUDE ");
   buf += 17;

   strcpy(buf, "-D BSACL_BASE_DIR=\"");
   buf += 19;
   strcpy(buf, BASE_DIRECTORY);
   buf += strlen(BASE_DIRECTORY);
   strcpy(buf, "\" \0");
   buf += 2;

#ifdef BSACL_USE_DOUBLE_PRECISION
   strcpy(buf, "-D BSACL_USE_DOUBLE_PRECISION ");
   buf += 33;
#endif

   strcpy(buf, " -D BSACL_KERNEL_ID=");
   buf += 20;
   sprintf(buf, "%-1u", kernel_id_);
   ++buf;

   strcpy(buf, " -D BSACL_WIND_PSD_ID=");
   buf += 22;
   strcpy(buf, STRINGIFYMACRO_VALUE(BSACL_PSD_TYPE_DAVENPORT));
   ++buf;

   if (pass_params_by_macro_ == 1) {
      strcpy(buf, " -D BSACL_PASS_PARAMS_BY_MACRO__");
      buf += 32;

      strcpy(buf, " -D NTC__=");
      buf += 10;
      sprintf(buf, "%-1u", extdata__.NTC__);
      ++buf;

      strcpy(buf, " -D NNL__=");
      buf += 10;
      sprintf(buf, "%-10u", extdata__.NNODES_LOAD__);
      buf += 10;

      strcpy(buf, " -D NN__=");
      buf += 9;
      sprintf(buf, "%-10u", extdata__.NN__);
      buf += 10;

      strcpy(buf, " -D NM_EFF__=");
      buf += 13;
      sprintf(buf, "%-5u", extdata__.NMODES_EFF__);
      buf += 5;

      if (kernel_id_ == 2) {
            strcpy(tbuf, " -D NFI__=");
            tbuf += 10;
            sprintf(tbuf, "%-10u", nfi__[i_thread]);
            tbuf += 10;

            strcpy(tbuf, " -D NFJ__=");
            tbuf += 10;
            sprintf(tbuf, "%-10u", nfj__[i_thread]);
            tbuf += 10;
      }
   }

#ifdef BSA_DEBUG
   strcpy(buf, " -g ");
   buf += 4;
#else
   strcpy(buf, " -cl-fast-relaxed-math ");
   buf += 23;
#endif
   *buf = '\0';

   // merge common flags into thread specific ones.
   strcat(thread_prog_opts__[i_thread], prog_compile_opts__);

#ifdef BSA_DEBUG
   if (i_thread == 0) {
      printf("%sBuild options:\n", INFO_MSG);
      printf("%s\t%s\n", CONT_MSG, thread_prog_opts__[i_thread]);
   }
#endif
   return;
}





cl_program buildCLProgram_(
      const unsigned i_thread, const unsigned device_id, ierr_t *ierr)
{
   cl_program program_ = NULL;
   size_t size_;
   char *clSource_;

   /** Create and build program. */
   size_ = readFileContent_(&clSource_, "bsacl.cl");
   size_ = assembleFinalCLSource_(&clSource_, evalFctStrings__);
   if (clSourceStrings__ == NULL) {
      *ierr = -324;
      printMsg_(ERRR_MSG, "Failed to assemble CL source strings.");
      goto ret_;
   }
   free(clSource_), clSource_ = NULL;
   if (evalFctStrings__ != NULL) { free(*evalFctStrings__); evalFctStrings__ = NULL; }

   program_ = clCreateProgramWithSource(context__, 1, (const char **)clSourceStrings__, &size_, ierr);
   if (*ierr != BSACL_SUCCESS) {
      printMsg_(ERRR_MSG, "Failed to create CL program object.");
      goto ret_;
   }

   assembleProgramBuildOptsString_(i_thread);

   *ierr = clBuildProgram(program_, 1, &(devices__[device_id]), thread_prog_opts__[i_thread], NULL, NULL);
   if (*ierr != BSACL_SUCCESS) {
      char buffer[20480];
      clGetProgramBuildInfo(program_, devices__[device_id],
         CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &size_);
      printf("\n Build log:\n\n%s\n", buffer);
      printMsg_(ERRR_MSG, "Failed to build CL program.");
      program_ = NULL;
      goto ret_;
   }

#ifdef BSACL_GENERATE_PTX_BINARY
   FILE *ptxfile_ = NULL;
   char *buffer;

   size_  = 0llu;
   buffer = (char *)malloc(sizeof(char) * 1024 * 1000); // NOTE: need to allocate enough memory!
   *ierr  = clGetProgramInfo(program_, CL_PROGRAM_BINARIES, 20480, &buffer, &size_);
   if (*ierr != BSACL_SUCCESS) {
      printMsg_(ERRR_MSG, "Failed to query CL program binaries.");
      goto ret_;
   }
   ptxfile_ = fopen("program.ptx", "w");
   if (ptxfile_ == NULL) {
      *ierr = -544;
      printMsg_(ERRR_MSG, "Failed to open file for writing PTX.");
      goto ret_;
   }
   fprintf(ptxfile_, "%s", buffer);
   fclose(ptxfile_);
   free(buffer);
#endif

ret_:
   return program_;
}




ierr_t setKernelArgs_(const unsigned i_thread, BSACL_MEM_REAL_T *d_res)
{
   ierr_t ierr_ = 0;
   BSACL_UINT iarg_ = 0;

   if (thread_set_shared_kernel_args_[i_thread] == 0U) {
      if (pass_params_by_macro_ == 0U)
         ierr_ |= clSetKernelArg(kernels__[i_thread], iarg_++,  sizeof(unsigned int), &extdata__.NTC__);
      ierr_ |= clSetKernelArg(kernels__[i_thread],    iarg_++,  sizeof(cl_mem),       &d_tc__);

      if (pass_params_by_macro_ == 0U) {
         ierr_ |= clSetKernelArg(kernels__[i_thread], iarg_++, sizeof(unsigned int), &extdata__.NNODES_LOAD__);
         ierr_ |= clSetKernelArg(kernels__[i_thread], iarg_++, sizeof(unsigned int), &extdata__.NMODES_EFF__);
         ierr_ |= clSetKernelArg(kernels__[i_thread], iarg_++, sizeof(unsigned int), &extdata__.NDEGW__);
      }
      ierr_ |= clSetKernelArg(kernels__[i_thread], iarg_++, sizeof(cl_mem), &d_nodes_load__);
      ierr_ |= clSetKernelArg(kernels__[i_thread], iarg_++, sizeof(cl_mem), &d_phiTc__);

      if (pass_params_by_macro_ == 0U) {
         ierr_ |= clSetKernelArg(kernels__[i_thread], iarg_++, sizeof(unsigned int), &extdata__.NN__);
         ierr_ |= clSetKernelArg(kernels__[i_thread], iarg_++, sizeof(unsigned int), &extdata__.NNOD_CORR__);
      }
      ierr_ |= clSetKernelArg(kernels__[i_thread], iarg_++, sizeof(cl_mem), &d_nod_corr__);

      ierr_ |= clSetKernelArg(kernels__[i_thread], iarg_++, sizeof(cl_mem), &d_wind_nodal_vel__);
      ierr_ |= clSetKernelArg(kernels__[i_thread], iarg_++, sizeof(cl_mem), &d_wind_turb_scales__);
      ierr_ |= clSetKernelArg(kernels__[i_thread], iarg_++, sizeof(cl_mem), &d_wind_turb_std__);
      ierr_ |= clSetKernelArg(kernels__[i_thread], iarg_++, sizeof(cl_mem), &d_wind_nodal_windz__);

      ierr_ |= clSetKernelArg(kernels__[i_thread], iarg_++, sizeof(cl_mem), &d_Mg__);
      ierr_ |= clSetKernelArg(kernels__[i_thread], iarg_++, sizeof(cl_mem), &d_Cg__);
      ierr_ |= clSetKernelArg(kernels__[i_thread], iarg_++, sizeof(cl_mem), &d_Kg__);

      thread_set_shared_kernel_args_[i_thread] = 1U;
   } else {
      if (pass_params_by_macro_ == 0U) iarg_ = 17U;
      else iarg_ = 11U;
   }

   if (pass_params_by_macro_ == 0U)
      ierr_ |= clSetKernelArg(kernels__[i_thread], iarg_++,  sizeof(unsigned int), &(nfi__[i_thread]));
   ierr_ |= clSetKernelArg(kernels__[i_thread], iarg_++,     sizeof(cl_mem),       &(d_fi__[i_thread]));
   if (pass_params_by_macro_ == 0U)
      ierr_ |= clSetKernelArg(kernels__[i_thread], iarg_++,  sizeof(unsigned int), &(nfj__[i_thread]));
   ierr_ |= clSetKernelArg(kernels__[i_thread], iarg_++,     sizeof(cl_mem),       &(d_fj__[i_thread]));

   ierr_ |= clSetKernelArg(kernels__[i_thread], iarg_++, sizeof(cl_mem), d_res);

   return ierr_;
}
#endif // BSACL_USE_CUDA__





ierr_t initThreadMem_(const unsigned n_threads)
{
#ifndef BSACL_USE_CUDA__
   thread_set_shared_kernel_args_ = (unsigned int*)malloc(sizeof(unsigned int) * n_threads);
   if (thread_set_shared_kernel_args_ == NULL) return -1;

   thread_prog_opts__ = (char **)malloc(sizeof(char*) * n_threads);
   if (thread_prog_opts__ == NULL) return -1;

   for (unsigned i = 0; i < n_threads; ++i) {

      thread_set_shared_kernel_args_[i] = 0U;

      thread_prog_opts__[i] = (char *)malloc(sizeof(char) * BUFSIZ);
      if (thread_prog_opts__[i] == NULL) return -1;
      thread_prog_opts__[i][0] = ' ';
      memset(thread_prog_opts__[i] + 1, '\0', BUFSIZ-1);
   }
#endif

   nfi__  = (unsigned *)malloc(sizeof(unsigned) * n_threads);
   if (nfi__ == NULL) return BSACL_ERROR_CODE(-1);

   nfj__  = (unsigned *)malloc(sizeof(unsigned) * n_threads);
   if (nfj__ == NULL) return BSACL_ERROR_CODE(-1);

   fi__   = (__real **)malloc(sizeof(__real*) * n_threads);
   if (fi__ == NULL) return BSACL_ERROR_CODE(-1);

   fj__   = (__real **)malloc(sizeof(__real*) * n_threads);
   if (fj__ == NULL) return BSACL_ERROR_CODE(-1);

   d_fi__ = (BSACL_MEM_REAL_T *)malloc(sizeof(BSACL_MEM_REAL_T) * n_threads);
   if (d_fi__ == NULL) return BSACL_ERROR_CODE(-1);

   d_fj__ = (BSACL_MEM_REAL_T *)malloc(sizeof(BSACL_MEM_REAL_T) * n_threads);
   if (d_fj__ == NULL) return BSACL_ERROR_CODE(-1);

   for (unsigned i = 0; i < n_threads; ++i) {
      fi__[i]   = NULL;
      fj__[i]   = NULL;
      d_fi__[i] = NULL;
      d_fj__[i] = NULL;
   }

   return BSACL_SUCCESS;
}





ierr_t createDeviceBuffer_(BSACL_MEM *d_buf, const size_t sz, unsigned flags, void *h_ptr)
{
   ierr_t ierr_ = BSACL_SUCCESS;

   // free device buffer if exists.
   if (*d_buf != NULL) {
      ierr_ = BSACL_DEVICE_FREE_MEM(*d_buf);
      if (ierr_ != BSACL_SUCCESS) {
         printMsg_(ERRR_MSG, "Failed to release device buffer.");
         goto err_;
      }
   }

#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMalloc(d_buf, sz);
#else
   *(d_buf) = clCreateBuffer(context__, flags, sz, h_ptr, &ierr_);
#endif

   if (ierr_ != BSACL_SUCCESS) {
      printMsg_(ERRR_MSG, "Failed to (re)create device buffer");
      goto err_;
   }

#ifdef BSACL_USE_CUDA__
   if (h_ptr) {
      ierr_ = cudaMemcpy((void *) *d_buf, h_ptr, sz, cudaMemcpyHostToDevice);
      if (ierr_ != BSACL_SUCCESS) {
         printMsg_(ERRR_MSG, "Failed to copy memory from host to device.");
         goto err_;
      }
   }
#endif

   goto ret_; // correct flow
err_:
   if (*d_buf) ierr_ = BSACL_DEVICE_FREE_MEM(*d_buf);
   ierr_ = ierr_ == BSACL_SUCCESS ? BSACL_ERROR_CODE(-45) : ierr_;
ret_:
   return ierr_;
}



// ************************************************************************************************
// ************************************************************************************************
//
//                        SECTION:   EXPOSED library functions (bsacl.h)
//
// ************************************************************************************************
// ************************************************************************************************



/**
 * @brief Initialises BsaCL runtime
 * @note  External data has to be fully initialised for this function to work!
 *
 * @return (int) Error code
 */
int bsaclInit(unsigned n_threads)
{
   if (has_halted__ == 1U) { return -999; }

   BSACL_INT ierr_ = BSACL_SUCCESS;

#ifndef BSACL_USE_CUDA__

   // Acquire platform
   ierr_ = getCLPlatformFromInfo_(&platform__, CL_PLATFORM_NAME, "NVIDIA\0", 0);
   if (ierr_ != BSACL_SUCCESS) goto ret_;

   // Query devices by type
   n_devices__ = getCLDevicesByType_(&platform__, device_type__, &devices__);
   if (n_devices__ <= 0) {
      ierr_ = (BSACL_INT)n_devices__;
      goto ret_;
   }

   getCLDevicesInfo_();

   context__ = clCreateContext(NULL, n_devices__, devices__, NULL, NULL, &ierr_);
   if (ierr_ != BSACL_SUCCESS) ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to create context", ret_);

   // One command queue per device (for the moment)
   cqueues__ = (cl_command_queue *) malloc(n_devices__ * sizeof(cl_command_queue));
   if (cqueues__ == NULL) {
      ierr_ = (cl_int)BSACL_CQUEUES_CREATION_ERROR_;
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to create command queues.", ret_);
   }
   for (unsigned i = 0; i < n_devices__; ++i) {
      cqueues__[i] = clCreateCommandQueue(context__, devices__[i], NULL, &ierr_);
      if (ierr_ != BSACL_SUCCESS) {
         printf("\n\n%sFailed to create  %u-th  command queue. Aborting (%d).", ERRR_MSG, i, ierr_);
         abortInternal_(0, NULL);
         goto ret_;
      }
   }

   // one program and kernel per thread.
   programs__ = (cl_program *)malloc(n_threads * sizeof(cl_program));
   if (programs__ == NULL) {
      ierr_ = (cl_int)BSACL_PROGRAMS_CREATION_ERROR_;
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to create program objects.", ret_);
   }

   kernels__ = (cl_kernel *)malloc(n_threads * sizeof(cl_kernel));
   if (kernels__ == NULL) {
      ierr_ = (cl_int)BSACL_KERNELS_CREATION_ERROR_;
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to create kernel objects.\n", ret_);
   }
#endif

   ierr_ = initThreadMem_(n_threads);
   if (ierr_ != BSACL_SUCCESS) ABORT_INTERNAL_GOTO_RET_(-32, "Failed to initialise thread memory.", ret_);

   n_threads__         = n_threads;
   has_cleaned__       = 0U; // if reset, then we'd need to clean again.
   bsacl_init_called__ = 1U; // NOTE: put it here, we avoid set to 1 (true) if error occurs

ret_:
   return ierr_;
}




int bsaclInitDeviceMemory(void)
{
   ierr_t ierr_;
   size_t sz_;

   if (extdata__.tc__ == NULL) {
      ierr_ = BSACL_ERROR_CODE(-83);
      ABORT_INTERNAL_GOTO_RET_(ierr_, "List of turb. components was not acquired. Aborting.", ret_);
   }
   sz_ = extdata__.NTC__ * sizeof(unsigned int);
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMalloc((void**)&d_tc__, sz_);
#else
   d_tc__ = clCreateBuffer(
      context__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sz_, (void *)extdata__.tc__, &ierr_);
#endif
   if (ierr_ != BSACL_SUCCESS || d_tc__==NULL)
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to create  d_tc__  device buffer", ret_);
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMemcpy((void *)d_tc__, (void *)extdata__.tc__, sz_, cudaMemcpyHostToDevice);
   if (ierr_ != BSACL_SUCCESS) 
         ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to write device buffer  d_tc__.", ret_);
#endif


   if (extdata__.nodes_load__ == NULL) {
      ierr_ = BSACL_ERROR_CODE(-84);
      ABORT_INTERNAL_GOTO_RET_(ierr_, "List of loaded nodes was not acquired. Aborting.", ret_);
   }
   sz_ = extdata__.NNODES_LOAD__ * sizeof(unsigned int);
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMalloc((void**)&d_nodes_load__, sz_);
#else
   d_nodes_load__ = clCreateBuffer(
      context__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sz_, (void *)extdata__.nodes_load__, &ierr_);
#endif
   if (ierr_ != BSACL_SUCCESS || d_nodes_load__==NULL) 
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to create  d_nodes_load__  device buffer", ret_);
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMemcpy((void *)d_nodes_load__, (void *)extdata__.nodes_load__, sz_, cudaMemcpyHostToDevice);
   if (ierr_ != BSACL_SUCCESS) 
         ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to write device buffer  d_nodes_load__.", ret_);
#endif



   if (extdata__.phi_T_c__ == NULL) {
      ierr_ = BSACL_ERROR_CODE(-85);
      ABORT_INTERNAL_GOTO_RET_(ierr_, "phiTc was not acquired. Aborting.", ret_);
   }
   sz_ = extdata__.NMODES_EFF__ * extdata__.NNODES_LOAD__ * extdata__.NDEGW__ * sizeof(__real);
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMalloc((void**)&d_phiTc__, sz_);
#else
   d_phiTc__ = clCreateBuffer(
      context__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sz_, (void *)extdata__.phi_T_c__, &ierr_);
#endif
   if (ierr_ != BSACL_SUCCESS || d_phiTc__==NULL) 
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to create phi_T_c__ on device", ret_);
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMemcpy((void *)d_phiTc__, (void *)extdata__.phi_T_c__, sz_, cudaMemcpyHostToDevice);
   if (ierr_ != BSACL_SUCCESS) 
         ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to write device buffer  d_phiTc__.", ret_);
#endif

   if (extdata__.nod_corr__ == NULL) {
      ierr_ = BSACL_ERROR_CODE(-86);
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Nodal spatial correlation was not acquired. Aborting.", ret_);
   }
   sz_ = extdata__.NNOD_CORR__ * extdata__.NTC__ * sizeof(__real);
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMalloc((void**)&d_nod_corr__, sz_);
#else
   d_nod_corr__ = clCreateBuffer(
      context__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sz_, (void *)extdata__.nod_corr__, &ierr_);
#endif
   if (ierr_ != BSACL_SUCCESS || d_nod_corr__==NULL) 
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to create d_nod_corr__ on device", ret_);
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMemcpy((void *)d_nod_corr__, (void *)extdata__.nod_corr__, sz_, cudaMemcpyHostToDevice);
   if (ierr_ != BSACL_SUCCESS) 
         ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to write device buffer  d_nod_corr__.", ret_);
#endif


   sz_ = extdata__.NN__ * sizeof(__real);
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMalloc((void**)&d_wind_nodal_vel__, sz_);
#else
   d_wind_nodal_vel__ = clCreateBuffer(
      context__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sz_, (void *)extdata__.wind_nodal_vel__, &ierr_);
#endif
   if (ierr_ != BSACL_SUCCESS || d_wind_nodal_vel__==NULL)
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to create d_wind_nodal_vel__ on device", ret_);
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMemcpy((void *)d_wind_nodal_vel__, (void *)extdata__.wind_nodal_vel__, sz_, cudaMemcpyHostToDevice);
   if (ierr_ != BSACL_SUCCESS) 
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to write device buffer  d_wind_nodal_vel__.", ret_);
#endif


   sz_ = extdata__.NN__ * sizeof(int);
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMalloc((void**)&d_wind_nodal_windz__, sz_);
#else
   d_wind_nodal_windz__ = clCreateBuffer(
      context__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sz_, (void *)extdata__.wind_nodal_windz__, &ierr_);
#endif
   if (ierr_ != BSACL_SUCCESS || d_wind_nodal_windz__==NULL) 
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to create d_wind_nodal_windz__ on device", ret_);
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMemcpy((void *)d_wind_nodal_windz__, (void *)extdata__.wind_nodal_windz__, sz_, cudaMemcpyHostToDevice);
   if (ierr_ != BSACL_SUCCESS) 
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to write device buffer  d_wind_nodal_windz__.", ret_);
#endif


   sz_ = (size_t)(extdata__.NWZ__ * 3 * sizeof(__real));
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMalloc((void**)&d_wind_turb_std__, sz_);
#else
   d_wind_turb_std__ = clCreateBuffer(
      context__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sz_, (void *)extdata__.wind_turb_std__, &ierr_);
#endif
   if (ierr_ != BSACL_SUCCESS || d_wind_turb_std__==NULL) 
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to create d_wind_turb_std__ on device", ret_);
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMemcpy((void *)d_wind_turb_std__, (void *)extdata__.wind_turb_std__, sz_, cudaMemcpyHostToDevice);
   if (ierr_ != BSACL_SUCCESS) 
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to write device buffer  d_wind_turb_std__.", ret_);
#endif

   sz_ *= 3;
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMalloc((void**)&d_wind_turb_scales__, sz_);
#else
   d_wind_turb_scales__ = clCreateBuffer(
      context__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sz_, (void *)extdata__.wind_turb_scales__, &ierr_);
#endif
   if (ierr_ != BSACL_SUCCESS || d_wind_turb_scales__==NULL) 
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to create d_wind_turb_scales__ on device", ret_);
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMemcpy((void *)d_wind_turb_scales__, (void *)extdata__.wind_turb_scales__, sz_, cudaMemcpyHostToDevice);
   if (ierr_ != BSACL_SUCCESS) 
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to write device buffer  d_wind_turb_scales__.", ret_);
#endif




   sz_ = extdata__.NMODES_EFF__ * sizeof(__real);
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMalloc((void**)&d_Mg__, sz_);
#else
   d_Mg__ = clCreateBuffer(
      context__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sz_, (void *)extdata__.Mg__, &ierr_);
#endif
   if (ierr_ != BSACL_SUCCESS)
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to create   d_Mg__   device buffer", ret_);
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMemcpy((void *)d_Mg__, (void *)extdata__.Mg__, sz_, cudaMemcpyHostToDevice);
   if (ierr_ != BSACL_SUCCESS) 
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to write device buffer  d_Mg__.", ret_);
#endif

#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMalloc((void**)&d_Kg__, sz_);
#else
   d_Kg__ = clCreateBuffer(
      context__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sz_, (void *)extdata__.Kg__, &ierr_);
#endif
   if (ierr_ != BSACL_SUCCESS) 
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to create   d_Kg__   device buffer", ret_);
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMemcpy((void *)d_Kg__, (void *)extdata__.Kg__, sz_, cudaMemcpyHostToDevice);
   if (ierr_ != BSACL_SUCCESS) 
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to write device buffer  d_Kg__.", ret_);
#endif

   sz_ *= extdata__.NMODES_EFF__;
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMalloc((void**)&d_Cg__, sz_);
#else
   d_Cg__ = clCreateBuffer(
      context__, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sz_, (void *)extdata__.Cg__, &ierr_);
#endif
   if (ierr_ != BSACL_SUCCESS) 
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to create   d_Cg__   device buffer", ret_);
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMemcpy((void *)d_Cg__, (void *)extdata__.Cg__, sz_, cudaMemcpyHostToDevice);
   if (ierr_ != BSACL_SUCCESS)
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to write device buffer  d_Cg__.", ret_);
#endif

   bsacl_init_devmem_called__ = 1U;

ret_:
   return ierr_;
}





/**
 * @brief BSACL core function. Enqueue selected kernel to run on GPU(s).
 *
 */
int bsaclRun(unsigned i_thread, const size_t dim, __real * __RESTRICT_PTR __EXT_PTR_CONST h_res)
{
   ierr_t ierr_;
#ifndef BSACL_USE_CUDA__
   char* clSource_ = NULL;
#endif

   if (has_halted__ == 1U) {
      ierr_ = BSACL_ERROR_CODE(-1);
      goto ret_;
   }

   ierr_ = BSACL_SUCCESS;

   // Init BsaCL runtime, if needed. Or, validate it.
   if (bsacl_init_called__ == 0U) {
      n_threads__ = 1U;
      ierr_ = (ierr_t)bsaclInit(1); // use only 1 thread if called internally
      if (ierr_ != BSACL_SUCCESS) {
         printMsg_(ERRR_MSG, "Failed to initialise BsaCL."); 
         goto ret_;
      }
      if (bsacl_init_devmem_called__ == 0U) {
         ierr_ = (ierr_t)bsaclInitDeviceMemory();
         if (ierr_ != BSACL_SUCCESS) 
            ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to initialise BsaCL device memory.", ret_);
      }
   } else {
      /*
       * Init called externally before Run.
       * However, it might that requirements have updated, so we need to check 
       * this run's compatibility with runtime initialisation.
       * */
      if (i_thread >= n_threads__) {
         ierr_ = BSACL_ERROR_CODE(-10);
         ABORT_INTERNAL_GOTO_RET_(ierr_, "Runtime has to be re-initialised due to incompatible settings.", ret_);
      }
      if (bsacl_init_devmem_called__ == 0U) {
         ierr_ = BSACL_ERROR_CODE(-101);
         ABORT_INTERNAL_GOTO_RET_(ierr_, "Device memory has not been initialised. Aborting.", ret_);
      }
   }


#ifdef BSACL_USE_CUDA__
   kernel_id_ = BSACL_KERNEL_ID;
#else
   cl_event mainKerEv_;
#endif
   size_t size_ = 0ull;


   if (fi__ == NULL || fj__ == NULL || fi__[i_thread]==NULL || fj__[i_thread] == NULL)
      { ierr_ = BSACL_ERROR_CODE(-20); printMsg_(ERRR_MSG, "No computation frequencies acquired."); goto ret_; }

   __real dInfl_  = (*(fi__[i_thread] + 1) - *fi__[i_thread]);
   dInfl_ *= (*(fj__[i_thread] + 1) - *fj__[i_thread]);
#ifdef BSACL_CONV_PULSATION
   dInfl_ *= 4 * (__real)(BSACL_PI * BSACL_PI);
#endif


#ifdef BSA_DEBUG
   printf("%sUsing kernel ID  %u\n", INFO_MSG, kernel_id_);
#endif /* ifdef BSA_DEBUG */


   /**
    * Define Kernel dimensions. Inlne here, so that each CPU thread has its own copy.
    * */
   size_t ntwi_;
   BSACL_UINT nwg_;
#ifdef BSACL_USE_CUDA__
   dim3 cuda_nblocks__, cuda_nthreads_perblock__;
#else
   const BSACL_UINT n_work_dims__ = 2U;
#endif
   size_t global_dims_ie_NTWI__[3] = { 0 };
   size_t local_dims_ie_WIpWG__[3] = { 0 };
   size_t n_work_groups__[3] = { 0 };

   // NOTE: use power function ??
   ntwi_ = extdata__.NMODES_EFF__  * extdata__.NMODES_EFF__  * extdata__.NMODES_EFF__;

   global_dims_ie_NTWI__[1] = ntwi_; // NM^3
   if (kernel_id_ == 2) {

      // find how many WG of given size fit in nfreqs
      size_t nfreqs_ = nfi__[i_thread] * nfj__[i_thread];
      nwg_ = (unsigned)(nfreqs_ / BSACL_WIpWG);
      ++nwg_;

      global_dims_ie_NTWI__[0] = nwg_*BSACL_WIpWG;

      local_dims_ie_WIpWG__[0] = BSACL_WIpWG;
      local_dims_ie_WIpWG__[1] = 1;

   } else {

      ntwi_ = extdata__.NNODES_LOAD__ * extdata__.NNODES_LOAD__ * extdata__.NNODES_LOAD__;
      nwg_  = (BSACL_UINT)ceil((double)ntwi_ / BSACL_WIpWG);

      global_dims_ie_NTWI__[0] = BSACL_WIpWG * nwg_;

      local_dims_ie_WIpWG__[0] = BSACL_WIpWG;
      local_dims_ie_WIpWG__[1] = 1;
   }
   n_work_groups__[0] = nwg_;
   n_work_groups__[1] = global_dims_ie_NTWI__[1];

#ifdef BSACL_USE_CUDA__
   cuda_nblocks__.x           = (unsigned int)n_work_groups__[0];
   cuda_nblocks__.y           = (unsigned int)n_work_groups__[1];
   cuda_nblocks__.z           = 1;
   cuda_nthreads_perblock__.x = (unsigned int)local_dims_ie_WIpWG__[0];
   cuda_nthreads_perblock__.y = (unsigned int)local_dims_ie_WIpWG__[1];
   cuda_nthreads_perblock__.z = 1;
#endif


#ifndef BSACL_USE_CUDA__
   programs__[i_thread] = buildCLProgram_(i_thread, /* device id */ 0, &ierr_);
   if (ierr_ != BSACL_SUCCESS) 
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to create CL program object.", ret_);

   /** Create kernel */
   kernels__[i_thread] = clCreateKernel(programs__[i_thread], "bfm_kernel", &ierr_);
   if (ierr_ != BSACL_SUCCESS) 
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to create CL kernel object.", ret_);
#endif // not def BSACL_USE_CUDA__


   size_t ntot_ = 0;
   ntot_ = nfi__[i_thread] * sizeof(__real);
   ierr_ = createDeviceBuffer_(
      BSACL_MEM_CAST(void **) &(d_fi__[i_thread]), ntot_, BSACL_MEM_READ_ONLY | BSACL_MEM_COPY_HOST_PTR, fi__[i_thread] );
   if (ierr_ != BSACL_SUCCESS) 
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to create d_fi__ device buffer.\n", ret_);

   ntot_ = nfj__[i_thread] * sizeof(__real);
   ierr_ = createDeviceBuffer_(
      BSACL_MEM_CAST(void **) &(d_fj__[i_thread]), ntot_, BSACL_MEM_READ_ONLY | BSACL_MEM_COPY_HOST_PTR, fj__[i_thread] );
   if (ierr_ != BSACL_SUCCESS) 
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to create d_fj__ device buffer.\n", ret_);


   // create result buffer on device.
   if (h_res == NULL) {
      ierr_ = BSACL_ERROR_CODE(-123);
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Host result array pointer is an invalid pointer.", ret_);
   }
   BSACL_MEM_REAL_T d_res_ = NULL;
   ntot_ = dim * n_work_groups__[0] * sizeof(__real);
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMalloc((void**)&d_res_, ntot_);
#else
   d_res_ = clCreateBuffer(
      context__, BSACL_MEM_READ_WRITE, ntot_, NULL, &ierr_);
#endif
   if (ierr_ != BSACL_SUCCESS || d_res_ == NULL) {
      if (d_res_) {
         ierr_ = BSACL_DEVICE_FREE_MEM(d_res_);
      ierr_ = BSACL_ERROR_CODE(-56);
      }
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to create  result  device buffer", ret_);
   }


#ifndef BSACL_USE_CUDA__
   // make sure to set passing by macro to 1 BEFORE setting kernel arguments.
   if (kernel_id_ == 2 && pass_params_by_macro_ == 0) pass_params_by_macro_ = 1U;

   ierr_ = setKernelArgs_(i_thread, &d_res_);
   if (ierr_ != BSACL_SUCCESS)
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to set kernel arguments.", ret_);
#endif


// #ifdef BSA_DEBUG
   printf("\n%sEnqueueing kernel using:\n", INFO_MSG);
   printf("%s - local grid dims  : %llu - %llu - %llu\n",
      CONT_MSG, local_dims_ie_WIpWG__[0], local_dims_ie_WIpWG__[1], local_dims_ie_WIpWG__[2]);
   printf("%s - global grid dims : %llu - %llu - %llu\n",
      CONT_MSG, global_dims_ie_NTWI__[0], global_dims_ie_NTWI__[1], global_dims_ie_NTWI__[2]);
   printf("%s - n. of WG/TB      : %llu - %llu - %llu\n\n",
      CONT_MSG, n_work_groups__[0], n_work_groups__[1], 1llu);
// #endif


#ifdef BSACL_USE_CUDA__
         bfm_kernel<<<cuda_nblocks__, cuda_nthreads_perblock__>>>(\
            extdata__.PSD_ID__,
            extdata__.NTC__,
            d_tc__,
            extdata__.NNODES_LOAD__,
            extdata__.NMODES_EFF__,
            extdata__.NDEGW__,
            d_nodes_load__,
            d_phiTc__,
            extdata__.NN__,
            extdata__.NNOD_CORR__,
            d_nod_corr__,
            d_wind_nodal_vel__,
            d_wind_turb_scales__,
            d_wind_turb_std__,
            d_wind_nodal_windz__,
            d_Mg__,
            d_Cg__,
            d_Kg__,
            nfi__[i_thread],
            d_fi__[i_thread],
            nfj__[i_thread],
            d_fj__[i_thread],
            d_res_);
         ierr_ = cudaDeviceSynchronize();
         // ierr_ = cudaGetLastError();
#else
         ierr_ = clEnqueueNDRangeKernel(
            cqueues__[i_thread], kernels__[i_thread], n_work_dims__, NULL, 
               global_dims_ie_NTWI__, local_dims_ie_WIpWG__, 
                  0, NULL, &mainKerEv_);
#endif // BSACL_USE_CUDA__

         if (ierr_ != BSACL_SUCCESS)
            ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to enqueue kernel execution.", ret_);

#ifndef BSACL_USE_CUDA__
         clWaitForEvents(1, &mainKerEv_);
         // ierr_ = clFinish(cqueues__[i_thread]);
         // if (ierr_ != BSACL_SUCCESS) ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to finish command queue (2).", ret_);
#endif

   printf("\n");

   // clReleaseEvent(mainKerEv_);
   __real *rtmp_;
   size_ = sizeof(__real) * dim * n_work_groups__[0];
   rtmp_ = (__real *)malloc(size_);
#ifdef BSACL_USE_CUDA__
   ierr_ = cudaMemcpy(rtmp_, d_res_, size_, cudaMemcpyDeviceToHost);
#else
   ierr_ = clEnqueueReadBuffer(cqueues__[i_thread], d_res_, CL_TRUE, 0, size_, rtmp_, 0, NULL, NULL);
#endif

   if (ierr_ != BSACL_SUCCESS)
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to read results from device.", ret_);

   /** Cumulate all single WG contributions */
   if (kernel_id_ == 1) {
      for (BSACL_UINT iwgx_ = 0; iwgx_ < n_work_groups__[0]; iwgx_++) {
         for (BSACL_UINT i_ = 0; i_ < dim; i_++) {
            h_res[i_] += rtmp_[iwgx_*dim + i_];
         }
      }
      for (BSACL_UINT i_ = 0; i_ < dim; i_++) {
         h_res[i_] *= dInfl_;
      }
   } else {
      for (BSACL_UINT i_ = 0; i_ < dim; i_++) {
         h_res[i_] = 0.;
         for (BSACL_UINT iwgx_ = 0; iwgx_ < n_work_groups__[0]; iwgx_++) {
            h_res[i_] += rtmp_[iwgx_ + (i_*n_work_groups__[0])];
         }
         h_res[i_] *= dInfl_;
      }
   }
   free(rtmp_);
   rtmp_ = NULL;
   ierr_ = BSACL_DEVICE_FREE_MEM(d_res_);
   if (ierr_ != BSACL_SUCCESS) 
      ABORT_INTERNAL_GOTO_RET_(ierr_, "Failed to release result device memory.", ret_);


#ifndef BSACL_USE_CUDA__
   ierr_ = clFinish(cqueues__[i_thread]);
   if (ierr_ != BSACL_SUCCESS) 
      printMsg_(ERRR_MSG, "Failed to finish default command queue.");
#endif


ret_:
#ifndef BSACL_USE_CUDA__
   if (clSource_) free(clSource_);
#endif
   return ierr_;
}




BSACL_INT bsaclSetKernelID(unsigned kid)
{
   ierr_t ret = { BSACL_SUCCESS };
#ifndef BSACL_USE_CUDA__
   if (kid < 1 || kid > 2) {
      ret = BSACL_ERROR_CODE(1);
   } else {
      kernel_id_ = kid;
   }
#endif
   return ret;
}





void bsaclAcquirePSDId(const uint32_t psdid)
{
   if (has_halted__ == 1U) return;
   extdata__.PSD_ID__ = psdid;
}




/**
 * @brief Get pointer to memory containing Modal Matrix
 * 
 * @note For dimensions, no need to pass any reference. Copy.
 * 
 * @param modmat pointer to base memory address
 * @param natf   vector of modal natural frequencies
 * @param ndofs  n. of total degrees of freedom
 * @param nmodes n. of modes (used)
 */
void bsaclAcquireStructModMat(
      __real * __RESTRICT_PTR __EXT_PTR_CONST modmat, __real * __RESTRICT_PTR __EXT_PTR_CONST natf, const uint32_t ndofs, const uint32_t nmodes)
{
   if (has_halted__ == 1U) return;
   extdata__.modmat__   = modmat;
   extdata__.natfreqs__ = natf;
   extdata__.NDOFS__    = ndofs;
   extdata__.NMODES__   = nmodes;
// #ifdef BSA_DEBUG
//    for (unsigned i = 0; i < ndofs; ++i) {
//       printf("\n");
//       for (unsigned j = 0; j < nmodes; ++j) {
//          printf("%10.4g  ", extdata__.modmat__[j*ndofs + i]);
//       }
//    } 
// #endif
}



void bsaclAcquireModalMatrices(
      __real * __RESTRICT_PTR __EXT_PTR_CONST Mg, __real * __RESTRICT_PTR __EXT_PTR_CONST Cg, __real * __RESTRICT_PTR __EXT_PTR_CONST Kg)
{
   if (has_halted__ == 1U) return;
   extdata__.Mg__ = Mg;
   extdata__.Cg__ = Cg;
   extdata__.Kg__ = Kg;
}





/**
 * @brief Get ptr to memory containing loaded nodes list
 * 
 * @param nodes_load 
 * @param nnodes_l 
 */
void bsaclAcquireLoadedNodesList(int * __RESTRICT_PTR __EXT_PTR_CONST nodes_load, const uint32_t nnodes_l)
{
   if (has_halted__ == 1U) return;
   if (extdata__.NNODES_LOAD__ != 0 && nnodes_l != extdata__.NNODES_LOAD__)
      ABORT_INTERNAL_RETURN_VOID(BSACL_VALUE_MISMATCH_ERROR_, "N. of loaded nodes does not match.");
   extdata__.nodes_load__  = nodes_load;
   extdata__.NNODES_LOAD__ = nnodes_l;
// #ifdef BSA_DEBUG
//    printf("\n");
//    for (unsigned i = 0; i < nnodes_l; ++i) printf("  %4d", extdata__.nodes_load__[i]);
// #endif
}



void bsaclAcquireTotalNOfNodes(const uint32_t nnodes_tot)
{
   if (has_halted__ == 1U) return;
   if (extdata__.NN__ != 0 && nnodes_tot != extdata__.NN__)
      ABORT_INTERNAL_RETURN_VOID(BSACL_VALUE_MISMATCH_ERROR_, "N. of total nodes does not match.");
   extdata__.NN__  = nnodes_tot;
}


/**
 * @brief Get ptr to memory holding list of effective used modes
 * 
 * @param modes 
 * @param nmodes_eff 
 */
void bsaclAcquireUsedModesList(int * __RESTRICT_PTR __EXT_PTR_CONST modes, const uint32_t nmodes_eff)
{
   if (has_halted__ == 1U) return;
   if (extdata__.NMODES_EFF__ != 0 && nmodes_eff != extdata__.NMODES_EFF__)
      ABORT_INTERNAL_RETURN_VOID(BSACL_VALUE_MISMATCH_ERROR_, "N. of effective modes does not match.");
   extdata__.modes_eff__  = modes;
   extdata__.NMODES_EFF__ = nmodes_eff;
// #ifdef BSA_DEBUG
//    printf("\n");
//    for (unsigned i = 0; i < nmodes_eff; ++i) printf("  %4d", extdata__.modes_eff__[i]);
// #endif
}


/**
 * @brief Get ptr to memory holding Wind Force Coefficients data
 * 
 * @param wfc       wind force coefficient matrix
 * @param nnodes_l  n. of loaded nodes
 * @param nlibs     n. of degrees of freedom per node
 * @param ndegw     n. of coefficients per element (transformation degree)
 */
void bsaclAcquireWindCoeffs(
      __real * __RESTRICT_PTR __EXT_PTR_CONST wfc, const uint32_t nnodes_l, const uint32_t nlibs, const uint32_t ndegw)
{
   if (has_halted__ == 1U) return;
   if (extdata__.NNODES_LOAD__ != 0 && nnodes_l != extdata__.NNODES_LOAD__) 
      ABORT_INTERNAL_RETURN_VOID(BSACL_VALUE_MISMATCH_ERROR_, "N. of loaded nodes does not match.");
   if (extdata__.NLIBS__ != 0 && nlibs != extdata__.NLIBS__)
      ABORT_INTERNAL_RETURN_VOID(BSACL_VALUE_MISMATCH_ERROR_, "N. of LIBS does not match.");
   if (extdata__.NDEGW__ != 0 && ndegw != extdata__.NDEGW__)
      ABORT_INTERNAL_RETURN_VOID(BSACL_VALUE_MISMATCH_ERROR_, "DEGW does not match.");
   extdata__.NNODES_LOAD__ = nnodes_l;
   extdata__.NLIBS__       = nlibs;
   extdata__.NDEGW__       = ndegw;
   extdata__.wfc__         = wfc;
}


/**
 * @brief Get ptr to memory holding list of effective turbulent components considered
 *
 * @param tc   list of turbulent components
 * @param ntc  n. of turb. components in the list
 */
void bsaclAcquireTurbComponentsList(int * __RESTRICT_PTR __EXT_PTR_CONST tc, const uint32_t ntc)
{
   if (has_halted__ == 1U) return;
   extdata__.NTC__ = ntc;
   extdata__.tc__  = tc;
}


/**
 * @brief Get ptr to memory holding the product of Phi x C
 * 
 * @param phi_T_c    product of modal matrix times wind-force-coeffs matrix
 * @param nnodes_l   n. of loaded nodes
 * @param nmodes_eff n. of effective modes
 * @param ndegw      n. of coefficients per element (transformation degree)
 */
void bsaclAcquirePhiTimesCMat(
      __real * __RESTRICT_PTR __EXT_PTR_CONST phi_T_c, const uint32_t nmodes_eff, const uint32_t nnodes_l, const uint32_t ndegw)
{
   if (has_halted__ == 1U) return;
   if (extdata__.NNODES_LOAD__ != 0 && nnodes_l != extdata__.NNODES_LOAD__)
      ABORT_INTERNAL_RETURN_VOID(BSACL_VALUE_MISMATCH_ERROR_, "N. of loaded nodes does not match.");
   if (extdata__.NMODES_EFF__ != 0 && nmodes_eff != extdata__.NMODES_EFF__)
      ABORT_INTERNAL_RETURN_VOID(BSACL_VALUE_MISMATCH_ERROR_, "N. of effective modes does not match.");
   if (extdata__.NDEGW__ != 0 && ndegw != extdata__.NDEGW__)
      ABORT_INTERNAL_RETURN_VOID(BSACL_VALUE_MISMATCH_ERROR_, "DEGW does not match.");
   extdata__.NNODES_LOAD__ = nnodes_l;
   extdata__.NMODES_EFF__  = nmodes_eff;
   extdata__.NDEGW__       = ndegw;
   extdata__.phi_T_c__     = phi_T_c;
}




void bsaclAcquireNodalCorrelation(__real * __RESTRICT_PTR __EXT_PTR_CONST nod_corr, const uint32_t nnod_corr)
{
   if (has_halted__ == 1U) return;
   extdata__.nod_corr__  = nod_corr;
   extdata__.NNOD_CORR__ = nnod_corr;
}



void bsaclAcquireWindNodalVelocities(__real * __RESTRICT_PTR __EXT_PTR_CONST nod_vel)
{
   if (has_halted__ == 1U) return;
   extdata__.wind_nodal_vel__ = nod_vel;
}

void bsaclAcquireWindNodalWindZones(int * __RESTRICT_PTR __EXT_PTR_CONST nod_wz)
{
   if (has_halted__ == 1U) return;
   extdata__.wind_nodal_windz__ = nod_wz;
}

void bsaclAcquireWindTurbScales(__real * __RESTRICT_PTR __EXT_PTR_CONST wt_scl, const uint32_t nwz)
{
   if (has_halted__ == 1U) return;
   extdata__.NWZ__ = nwz;
   extdata__.wind_turb_scales__ = wt_scl;
}

void bsaclAcquireWindTurbStd(__real * __RESTRICT_PTR __EXT_PTR_CONST wt_std, const uint32_t nwz)
{
   if (has_halted__ == 1U) return;
   extdata__.NWZ__ = nwz;
   extdata__.wind_turb_std__ = wt_std;
}




#ifdef BSACL_ENABLE_EVALFCT_PTR
void bsaclAcquireEvalFunc(evalFct_t fct)
{
   if (has_halted__ == 1U) return;
#ifdef BSACL_ENABLE_EVALFCT_PTR
#if (BSACL_KERNEL_ID==0)
   if (fct == NULL) ABORT_INTERNAL_RETURN_VOID(-11, "Passed evaluation function pointer is NULL.");
   evalFunc_ptr__ = fct;
#endif
#endif
   return;
}


void bsaclAcquireEvalFuncByStrings(char ** __RESTRICT_PTR __EXT_PTR_CONST strings)
{
   if (has_halted__ == 1U) return;
#ifndef BSACL_USE_CUDA__
   size_t ilen_ = strlen(*strings);
   evalFctStrings__  = (char **)malloc(sizeof(char *));
   *evalFctStrings__ = (char  *)malloc(sizeof(char) * (ilen_ + 1));
   strcpy(*evalFctStrings__, *strings);
   *(*evalFctStrings__ + ilen_) = '\0';
#endif
   return;
}


void bsaclAcquireEvalFuncByFile(char * __RESTRICT_PTR __EXT_PTR_CONST filename)
{
   if (has_halted__ == 1U) return;
#ifndef BSACL_USE_CUDA__
   size_t ierr_   = 0;
   char *fstrings = NULL;
   ierr_ = readFileContent_(&fstrings, filename);
   bsaclAcquireEvalFuncByStrings(&fstrings);
   free(fstrings);
#endif
   return;
}
#endif




/**
 * @brief Acquires reference to computational frequencies arrays.
 *        I.e. desired spatial domain discretisation.
 *
 * */
void bsaclAcquireComputationFreqs(const unsigned i_thread,
      const uint32_t nfi, __real * __RESTRICT_PTR __EXT_PTR_CONST fi, const uint32_t nfj, __real * __RESTRICT_PTR __EXT_PTR_CONST fj) 
{
   if (has_halted__ == 1U) return;
   if (bsacl_init_called__ == 0U) {
      abortInternal_(-89, "This function must be called after BsaCL initialisation.");
      return;
   }
   if (n_threads__ <= i_thread) {
      abortInternal_(-103, "Current thread index exceeds max n. of threads.");
      return;
   }
   nfi__[i_thread] = nfi;
   fi__[i_thread]  = fi;
   nfj__[i_thread] = nfj;
   fj__[i_thread]  = fj;
   return;
}



/**
 * @brief Deprecated.
 *
 * */
void bsaclAcquireBaseWindTurbPSD(__real * __RESTRICT_PTR __EXT_PTR_CONST S_uvw)
{
   if (has_halted__ == 1U) return;
   S_uvw__ = S_uvw;
   return;
}




/**
 * @brief Specify desired device type to be used in computation.
 *        Valid types:
 *          - BSACL_DEVICE_TYPE_CPU_
 *          - BSACL_DEVICE_TYPE_GPU_
 *          - BSACL_DEVICE_TYPE_ACC_
 *          - BSACL_DEVICE_TYPE_DEF_
 *
 * */
void bsaclSetDeviceType(const uint32_t itype)
{
   if (has_halted__ == 1U) return;
#ifndef BSACL_USE_CUDA__
   switch (itype) {
      case BSACL_DEVICE_TYPE_CPU_:
         device_type__ = CL_DEVICE_TYPE_CPU;
         break;
      case BSACL_DEVICE_TYPE_GPU_:
         device_type__ = CL_DEVICE_TYPE_GPU;
         break;
      case BSACL_DEVICE_TYPE_ACC_:
         device_type__ = CL_DEVICE_TYPE_ACCELERATOR;
         break;
      case BSACL_DEVICE_TYPE_DEF_:
         break;
      default:
         abortInternal_(BSACL_INVALID_DEVICE_TYPE_, "Invalid device type requested.");
   }
#endif
   return;
}




/**
 * @brief Verifies if user can allocate `idim` n. of bytes in 
 *        selected device's memory. Returns 0 if not possible,
 *        1 otherwise.
 *
 * */
void bsaclVerifyMaxAllocCondition(size_t idim, unsigned int * __RESTRICT_PTR __EXT_PTR_CONST ican)
{
   if (has_halted__ == 1U) return;
#ifndef BSACL_USE_CUDA__
   if (dev_max_mem_alloc_size__[0] < idim * sizeof(BSACL_REAL)) {
      *ican = 0U;
      return;
   }
#endif
   *ican = 1U;
   return;
}



/**
 * @brief Aborts execution with error code.
 *
 * */
void bsaclAbort(const int ierr)
{
   abortInternal_(ierr, NULL);
}


/**
 * @brief Cleans BSACL runtime memory.
 *
 */
void bsaclFinalise(void)
{
   ierr_t iret_, ierr_;

   ierr_ = releaseDeviceMem_();
   if (ierr_ != BSACL_SUCCESS)
      printMsg_(ERRR_MSG, "Failed to release device memory.");
   iret_ = ierr_;

#ifndef BSACL_USE_CUDA__

   ierr_ = BSACL_SUCCESS;
   for (unsigned i = 0; i < n_devices__; ++i) {
      ierr_ |= clReleaseCommandQueue(cqueues__[i]);
   }
   if (ierr_ != BSACL_SUCCESS) 
      printMsg_(ERRR_MSG, "Failed to release command queue(s).");
   iret_ |= ierr_;

   ierr_ = BSACL_SUCCESS;
   for (unsigned i = 0; i < n_devices__; ++i) {
      ierr_ |= clReleaseKernel(kernels__[i]);
   }
   if (ierr_ != BSACL_SUCCESS) 
      printMsg_(ERRR_MSG, "Failed to release kernel.");
   iret_ |= ierr_;

   ierr_ = BSACL_SUCCESS;
   for (unsigned i = 0; i < n_devices__; ++i) {
      ierr_ = clReleaseProgram(programs__[i]);
   }
   if (ierr_ != BSACL_SUCCESS) 
      printMsg_(ERRR_MSG, "Failed to release program.");
   iret_ |= ierr_;
#endif

   abortInternal_(iret_, "Failed to release BsaCL memory.");
}



// extern "C" {
#ifdef __cplusplus
}
#endif

