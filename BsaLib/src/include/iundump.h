
#ifdef _OPENMP
# define __iun_dump omp_get_thread_num()+1
#else
# define __iun_dump 1
#endif

