#include "hpipm/blas/kernel.h"

void kernel_dvecld_inc1(int k, double* x) {
    int ii;
    double tmp;
    for (ii = 0; ii < k; ii++) {
        tmp = x[ii];
    }
}


void kernel_dveccp_inc1(int k, double* x, double* y) {
    int ii;
    for (ii = 0; ii < k; ii++) {
        y[ii] = x[ii];
    }
}
