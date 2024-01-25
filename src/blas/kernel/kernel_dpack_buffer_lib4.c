#include "tinyhpipm/blas/kernel.h"

// full non-transposed
void kernel_dpack_buffer_fn(int m, int n, double* A, int lda, double* pA, int sda) {

    const int ps = 4;

    int ii;

    for (ii = 0; ii < n - 3; ii += 4) {
        kernel_dpack_tt_4_lib4(m, A + ii * lda, lda, pA + ii * ps, sda);
    }
    if (ii < n) {
        kernel_dpack_tt_4_vs_lib4(m, A + ii * lda, lda, pA + ii * ps, sda, n - ii);
    }
}


// full transposed
void kernel_dpack_buffer_ft(int m, int n, double* A, int lda, double* pA, int sda) {

    const int ps = 4;

    int ii;

    for (ii = 0; ii < n - 3; ii += 4) {
        kernel_dpack_tn_4_lib4(m, A + ii * lda, lda, pA + ii * sda);
    }
    if (ii < n) {
        kernel_dpack_tn_4_vs_lib4(m, A + ii * lda, lda, pA + ii * sda, n - ii);
    }
}


// lower non-transposed
void kernel_dpack_buffer_ln(int m, double* A, int lda, double* pA, int sda) {

    const int ps = 4;

    int ii;

    for (ii = 0; ii < m - 3; ii += 4) {
        kernel_dpack_tt_4_lib4(m - ii, A + ii + ii * lda, lda, pA + ii * ps + ii * sda, sda);
    }
    if (ii < m) {
        kernel_dpack_tt_4_vs_lib4(m - ii, A + ii + ii * lda, lda, pA + ii * ps + ii * sda, sda, m - ii);
    }
}


// lower transposed
void kernel_dpack_buffer_lt(int m, double* A, int lda, double* pA, int sda) {

    const int ps = 4;

    int ii;

    for (ii = 0; ii < m - 3; ii += 4) {
        kernel_dpack_tn_4_lib4(m - ii, A + ii + ii * lda, lda, pA + ii * ps + ii * sda);
    }
    if (ii < m) {
        kernel_dpack_tn_4_vs_lib4(m - ii, A + ii + ii * lda, lda, pA + ii * ps + ii * sda, m - ii);
    }
}


// upper transposed
void kernel_dpack_buffer_ut(int m, double* A, int lda, double* pA, int sda) {

    const int ps = 4;

    int ii;

    for (ii = 0; ii < m - 3; ii += 4) {
        kernel_dpack_tn_4_lib4(ii + 4, A + ii * lda, lda, pA + ii * sda);
    }
    if (ii < m) {
        kernel_dpack_tn_4_vs_lib4(m, A + ii * lda, lda, pA + ii * sda, m - ii);
    }
}
