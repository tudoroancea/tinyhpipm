


#include <stdio.h>
#include <stdlib.h>

#include <blasfeo/blasfeo_block_size.h>
#include <blasfeo/blasfeo_stdlib.h>


void blasfeo_malloc(void** ptr, size_t size) {
    *ptr = malloc(size);
    return;
}


// allocate memory aligned to typical cache line size (64 bytes)
void blasfeo_malloc_align(void** ptr, size_t size) {
#if defined(OS_WINDOWS)
    *ptr = _aligned_malloc(size, CACHE_LINE_SIZE);
#elif defined(__DSPACE__)
    // XXX fix this hack !!! (Andrea?)
    *ptr = malloc(size);
#elif (defined __XILINX_NONE_ELF__ || defined __XILINX_ULTRASCALE_NONE_ELF_JAILHOUSE__)
    *ptr = memalign(CACHE_LINE_SIZE, size);
#else
    int err = posix_memalign(ptr, CACHE_LINE_SIZE, size);
    if (err != 0) {
        printf("Memory allocation error");
        exit(1);
    }
#endif
    return;
}


void blasfeo_free(void* ptr) {
    free(ptr);
    return;
}


void blasfeo_free_align(void* ptr) {
#if defined(OS_WINDOWS)
    _aligned_free(ptr);
#else
    free(ptr);
#endif
    return;
}
