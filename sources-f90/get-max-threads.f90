program get_max_threads
use omp_lib
    maxnp = OMP_GET_MAX_THREADS()
    print *,maxnp
stop
end
