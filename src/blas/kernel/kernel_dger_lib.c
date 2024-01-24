#include "hpipm/blas/kernel.h"


void kernel_dger_4_libc(int kmax, double* alpha, double* x, double* y, double* C, int ldc, double* D, int ldd) {

    if (kmax <= 0) { return; }


    int ii;

    double
            x_0,
            y_0, y_1, y_2, y_3;

    y_0 = *alpha * y[0];
    y_1 = *alpha * y[1];
    y_2 = *alpha * y[2];
    y_3 = *alpha * y[3];

    ii = 0;
    for (; ii < kmax - 3; ii += 4) {

        x_0 = x[0];
        D[0 + ldd * 0] = C[0 + ldc * 0] + x_0 * y_0;
        D[0 + ldd * 1] = C[0 + ldc * 1] + x_0 * y_1;
        D[0 + ldd * 2] = C[0 + ldc * 2] + x_0 * y_2;
        D[0 + ldd * 3] = C[0 + ldc * 3] + x_0 * y_3;

        x_0 = x[1];
        D[1 + ldd * 0] = C[1 + ldc * 0] + x_0 * y_0;
        D[1 + ldd * 1] = C[1 + ldc * 1] + x_0 * y_1;
        D[1 + ldd * 2] = C[1 + ldc * 2] + x_0 * y_2;
        D[1 + ldd * 3] = C[1 + ldc * 3] + x_0 * y_3;

        x_0 = x[2];
        D[2 + ldd * 0] = C[2 + ldc * 0] + x_0 * y_0;
        D[2 + ldd * 1] = C[2 + ldc * 1] + x_0 * y_1;
        D[2 + ldd * 2] = C[2 + ldc * 2] + x_0 * y_2;
        D[2 + ldd * 3] = C[2 + ldc * 3] + x_0 * y_3;

        x_0 = x[3];
        D[3 + ldd * 0] = C[3 + ldc * 0] + x_0 * y_0;
        D[3 + ldd * 1] = C[3 + ldc * 1] + x_0 * y_1;
        D[3 + ldd * 2] = C[3 + ldc * 2] + x_0 * y_2;
        D[3 + ldd * 3] = C[3 + ldc * 3] + x_0 * y_3;

        C += 4;
        D += 4;
        x += 4;
    }
    for (; ii < kmax; ii++) {

        x_0 = x[0];
        D[0 + ldd * 0] = C[0 + ldc * 0] + x_0 * y_0;
        D[0 + ldd * 1] = C[0 + ldc * 1] + x_0 * y_1;
        D[0 + ldd * 2] = C[0 + ldc * 2] + x_0 * y_2;
        D[0 + ldd * 3] = C[0 + ldc * 3] + x_0 * y_3;

        C += 1;
        D += 1;
        x += 1;
    }
}


static void kernel_dger_3_libc(int kmax, double* alpha, double* x, double* y, double* C, int ldc, double* D, int ldd) {

    if (kmax <= 0) { return; }


    int ii;

    double
            x_0,
            y_0, y_1, y_2;

    y_0 = *alpha * y[0];
    y_1 = *alpha * y[1];
    y_2 = *alpha * y[2];

    ii = 0;
    for (; ii < kmax - 3; ii += 4) {

        x_0 = x[0];
        D[0 + ldd * 0] = C[0 + ldc * 0] + x_0 * y_0;
        D[0 + ldd * 1] = C[0 + ldc * 1] + x_0 * y_1;
        D[0 + ldd * 2] = C[0 + ldc * 2] + x_0 * y_2;

        x_0 = x[1];
        D[1 + ldd * 0] = C[1 + ldc * 0] + x_0 * y_0;
        D[1 + ldd * 1] = C[1 + ldc * 1] + x_0 * y_1;
        D[1 + ldd * 2] = C[1 + ldc * 2] + x_0 * y_2;

        x_0 = x[2];
        D[2 + ldd * 0] = C[2 + ldc * 0] + x_0 * y_0;
        D[2 + ldd * 1] = C[2 + ldc * 1] + x_0 * y_1;
        D[2 + ldd * 2] = C[2 + ldc * 2] + x_0 * y_2;

        x_0 = x[3];
        D[3 + ldd * 0] = C[3 + ldc * 0] + x_0 * y_0;
        D[3 + ldd * 1] = C[3 + ldc * 1] + x_0 * y_1;
        D[3 + ldd * 2] = C[3 + ldc * 2] + x_0 * y_2;

        C += 4;
        D += 4;
        x += 4;
    }
    for (; ii < kmax; ii++) {

        x_0 = x[0];
        D[0 + ldd * 0] = C[0 + ldc * 0] + x_0 * y_0;
        D[0 + ldd * 1] = C[0 + ldc * 1] + x_0 * y_1;
        D[0 + ldd * 2] = C[0 + ldc * 2] + x_0 * y_2;

        C += 1;
        D += 1;
        x += 1;
    }
}


static void kernel_dger_2_libc(int kmax, double* alpha, double* x, double* y, double* C, int ldc, double* D, int ldd) {

    if (kmax <= 0) { return; }


    int ii;

    double
            x_0,
            y_0, y_1;

    y_0 = *alpha * y[0];
    y_1 = *alpha * y[1];

    ii = 0;
    for (; ii < kmax - 3; ii += 4) {

        x_0 = x[0];
        D[0 + ldd * 0] = C[0 + ldc * 0] + x_0 * y_0;
        D[0 + ldd * 1] = C[0 + ldc * 1] + x_0 * y_1;

        x_0 = x[1];
        D[1 + ldd * 0] = C[1 + ldc * 0] + x_0 * y_0;
        D[1 + ldd * 1] = C[1 + ldc * 1] + x_0 * y_1;

        x_0 = x[2];
        D[2 + ldd * 0] = C[2 + ldc * 0] + x_0 * y_0;
        D[2 + ldd * 1] = C[2 + ldc * 1] + x_0 * y_1;

        x_0 = x[3];
        D[3 + ldd * 0] = C[3 + ldc * 0] + x_0 * y_0;
        D[3 + ldd * 1] = C[3 + ldc * 1] + x_0 * y_1;

        C += 4;
        D += 4;
        x += 4;
    }
    for (; ii < kmax; ii++) {

        x_0 = x[0];
        D[0 + ldd * 0] = C[0 + ldc * 0] + x_0 * y_0;
        D[0 + ldd * 1] = C[0 + ldc * 1] + x_0 * y_1;

        C += 1;
        D += 1;
        x += 1;
    }
}


static void kernel_dger_1_libc(int kmax, double* alpha, double* x, double* y, double* C, int ldc, double* D, int ldd) {

    if (kmax <= 0) { return; }


    int ii;

    double
            x_0,
            y_0;

    y_0 = *alpha * y[0];

    ii = 0;
    for (; ii < kmax - 3; ii += 4) {

        x_0 = x[0];
        D[0 + ldd * 0] = C[0 + ldc * 0] + x_0 * y_0;

        x_0 = x[1];
        D[1 + ldd * 0] = C[1 + ldc * 0] + x_0 * y_0;

        x_0 = x[2];
        D[2 + ldd * 0] = C[2 + ldc * 0] + x_0 * y_0;

        x_0 = x[3];
        D[3 + ldd * 0] = C[3 + ldc * 0] + x_0 * y_0;

        C += 4;
        D += 4;
        x += 4;
    }
    for (; ii < kmax; ii++) {

        x_0 = x[0];
        D[0 + ldd * 0] = C[0 + ldc * 0] + x_0 * y_0;

        C += 1;
        D += 1;
        x += 1;
    }
}


void kernel_dger_4_vs_libc(int kmax, double* alpha, double* x, double* y, double* C, int ldc, double* D, int ldd, int km) {

    if (kmax <= 0) { return; }


    if (km == 1) {
        kernel_dger_1_libc(kmax, alpha, x, y, C, ldc, D, ldd);
    } else if (km == 2) {
        kernel_dger_2_libc(kmax, alpha, x, y, C, ldc, D, ldd);
    } else if (km == 3) {
        kernel_dger_3_libc(kmax, alpha, x, y, C, ldc, D, ldd);
    } else {
        kernel_dger_4_libc(kmax, alpha, x, y, C, ldc, D, ldd);
    }
}
