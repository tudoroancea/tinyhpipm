#include "hpipm/blas/kernel.h"


void kernel_dgetr_tn_4_lib(int kmax, double* A, int lda, double* C, int ldc) {

    int ii;
    ii = 0;

    for (; ii < kmax - 3; ii += 4) {
        C[0 + ldc * 0] = A[0 + lda * 0];
        C[1 + ldc * 0] = A[0 + lda * 1];
        C[2 + ldc * 0] = A[0 + lda * 2];
        C[3 + ldc * 0] = A[0 + lda * 3];

        C[0 + ldc * 1] = A[1 + lda * 0];
        C[1 + ldc * 1] = A[1 + lda * 1];
        C[2 + ldc * 1] = A[1 + lda * 2];
        C[3 + ldc * 1] = A[1 + lda * 3];

        C[0 + ldc * 2] = A[2 + lda * 0];
        C[1 + ldc * 2] = A[2 + lda * 1];
        C[2 + ldc * 2] = A[2 + lda * 2];
        C[3 + ldc * 2] = A[2 + lda * 3];

        C[0 + ldc * 3] = A[3 + lda * 0];
        C[1 + ldc * 3] = A[3 + lda * 1];
        C[2 + ldc * 3] = A[3 + lda * 2];
        C[3 + ldc * 3] = A[3 + lda * 3];

        A += 4;
        C += 4 * ldc;
    }
    for (; ii < kmax; ii++) {
        C[0 + ldc * 0] = A[0 + lda * 0];
        C[1 + ldc * 0] = A[0 + lda * 1];
        C[2 + ldc * 0] = A[0 + lda * 2];
        C[3 + ldc * 0] = A[0 + lda * 3];

        A += 1;
        C += 1 * ldc;
    }
}


void kernel_dgetr_tn_4_vs_lib(int kmax, double* A, int lda, double* C, int ldc, int m1) {

    if (m1 <= 0) { return; }

    int ii;
    ii = 0;

    if (m1 >= 4) {
        kernel_dgetr_tn_4_lib(kmax, A, lda, C, ldc);

    } else if (m1 == 1) {
        goto l1;
    } else if (m1 == 2) {
        goto l2;
    } else  // if(m1==3)
    {
        goto l3;
    }


l1:
    ii = 0;
    for (; ii < kmax; ii++) {
        C[0 + ldc * 0] = A[0 + lda * 0];

        A += 1;
        C += 1 * ldc;
    }


l2:
    ii = 0;
    for (; ii < kmax; ii++) {
        C[0 + ldc * 0] = A[0 + lda * 0];
        C[1 + ldc * 0] = A[0 + lda * 1];

        A += 1;
        C += 1 * ldc;
    }


l3:
    ii = 0;
    for (; ii < kmax; ii++) {
        C[0 + ldc * 0] = A[0 + lda * 0];
        C[1 + ldc * 0] = A[0 + lda * 1];
        C[2 + ldc * 0] = A[0 + lda * 2];

        A += 1;
        C += 1 * ldc;
    }
}
