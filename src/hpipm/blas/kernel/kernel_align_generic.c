#include "hpipm/blas/kernel.h"
#include <stdint.h>

void blasfeo_align_4096_byte(void* ptr, void** ptr_align) {
    *ptr_align = (void*) ((((uintptr_t) ptr) + 4095) / 4096 * 4096);
}

void blasfeo_align_64_byte(void* ptr, void** ptr_align) {
    *ptr_align = (void*) ((((uintptr_t) ptr) + 63) / 64 * 64);
}
