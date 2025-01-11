
#ifdef _OPENMP
# define __export_POD_trunc_id__  omp_get_thread_num()+1
#else
# define __export_POD_trunc_id__  1
#endif


#ifdef BSA_USE_OPTIMISED_OMP
# define __compute_bfm_inout__ intent(inout),
# define BSA_USES_ZONES_ARRAY
# ifndef BSA_USE_POD_DATA_CACHING
#   define __postmesh_use_per_thread_bfmundump
# endif
#else
# define __compute_bfm_inout__
# define __dump_needs_sync
# ifndef BSA_USE_POD_DATA_CACHING
#  define __print_premesh_update
#  ifndef _OPENMP
#    define __post_mesh_uses_global_bfmundump
#  endif
# endif
#endif


#if defined(BSA_DEBUG)
#  define __compute_validate_zone_layout
#endif


#if (defined(BSA_USE_POD_DATA_CACHING)) || (!defined(_BSA_M3MF_ONLY_PREMESH))
# define __interp_updates_m3mf
#endif

#if (defined(BSA_USE_POD_DATA_CACHING)) || (defined(_OPENMP))
# define __interp_use_new_proc
#endif



