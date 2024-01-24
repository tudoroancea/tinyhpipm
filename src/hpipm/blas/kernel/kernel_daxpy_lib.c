#include "hpipm/blas/kernel.h"
#include <stdio.h>
#include <stdlib.h>


void kernel_daxpy_11_lib(int n, double* ptr_alpha, double* x, double* y) {
    int ii;
    double alpha = *ptr_alpha;
    ii = 0;
    for (; ii < n - 3; ii += 4) {
        y[0] += alpha * x[0];
        y[1] += alpha * x[1];
        y[2] += alpha * x[2];
        y[3] += alpha * x[3];
        x += 4;
        y += 4;
    }
    for (; ii < n; ii++) {
        y[0] += alpha * x[0];
        x += 1;
        y += 1;
    }
}
