#if defined(TARGET_X64_AVX) || defined(TARGET_X64_AVX2)
#include <emmintrin.h>  // SSE2
#include <immintrin.h>  // AVX
#include <mmintrin.h>
#include <pmmintrin.h>  // SSE3
#include <smmintrin.h>  // SSE4
#include <xmmintrin.h>  // SSE
#endif

#include "hpipm/common.h"

void hpipm_zero_memset(hpipm_size_t memsize, void* mem) {
    hpipm_size_t ii;
    hpipm_size_t memsize_m8 = memsize / 8;  // sizeof(double) is 8
    hpipm_size_t memsize_r8 = memsize % 8;
    double* double_ptr = mem;
    ii = 0;
#if defined(TARGET_X64_AVX) || defined(TARGET_X64_AVX2)
    __m256d
            y_zeros;

    y_zeros = _mm256_setzero_pd();
    if (memsize_m8 > 7) {
        for (; ii < memsize_m8 - 7; ii += 8) {
            _mm256_storeu_pd(double_ptr + ii + 0, y_zeros);
            _mm256_storeu_pd(double_ptr + ii + 4, y_zeros);
        }
    }
#else
    if (memsize_m8 > 7) {
        for (; ii < memsize_m8 - 7; ii += 8) {
            double_ptr[ii + 0] = 0.0;
            double_ptr[ii + 1] = 0.0;
            double_ptr[ii + 2] = 0.0;
            double_ptr[ii + 3] = 0.0;
            double_ptr[ii + 4] = 0.0;
            double_ptr[ii + 5] = 0.0;
            double_ptr[ii + 6] = 0.0;
            double_ptr[ii + 7] = 0.0;
        }
    }
#endif
    for (; ii < memsize_m8; ii++) {
        double_ptr[ii] = 0.0;
    }
    char* char_ptr = (char*) (&double_ptr[ii]);
    for (ii = 0; ii < memsize_r8; ii++) {
        char_ptr[ii] = 0;
    }
    return;
}
