#include "hpipm/blas/kernel.h"
#include <stdio.h>
#include <stdlib.h>


void kernel_ddot_11_lib(int n, double* x, double* y, double* res) {

    int ii;

    double tmp_res[4];
    tmp_res[0] = 0.0;
    tmp_res[1] = 0.0;
    tmp_res[2] = 0.0;
    tmp_res[3] = 0.0;

    ii = 0;
    for (; ii < n - 3; ii += 4) {
        tmp_res[0] += x[0] * y[0];
        tmp_res[1] += x[1] * y[1];
        tmp_res[2] += x[2] * y[2];
        tmp_res[3] += x[3] * y[3];
        x += 4;
        y += 4;
    }
    for (; ii < n; ii++) {
        tmp_res[0] += x[0] * y[0];
        x += 1;
        y += 1;
    }

    *res = tmp_res[0] + tmp_res[1] + tmp_res[2] + tmp_res[3];
}
