#include "hpipm/blas/kernel.h"

// XXX copy and scale y_n into z_n outside the kernel !!!!!
void kernel_dgemv_n_4_libc(int kmax, double* alpha, double* A, int lda, double* x, double* z) {

    if (kmax <= 0) {return;}


        int k;

    double
            a_00,
            a_01, a_02, a_03,
            x_0, x_1, x_2, x_3, y_0;

    x_0 = alpha[0] * x[0];
    x_1 = alpha[0] * x[1];
    x_2 = alpha[0] * x[2];
    x_3 = alpha[0] * x[3];

    k = 0;
    for (; k < kmax - 3; k += 4) {

        // 0
        y_0 = z[0];

        a_00 = A[0 + lda * 0];
        a_01 = A[0 + lda * 1];
        a_02 = A[0 + lda * 2];
        a_03 = A[0 + lda * 3];

        y_0 += a_00 * x_0;
        y_0 += a_01 * x_1;
        y_0 += a_02 * x_2;
        y_0 += a_03 * x_3;

        z[0] = y_0;


        // 1
        y_0 = z[1];

        a_00 = A[1 + lda * 0];
        a_01 = A[1 + lda * 1];
        a_02 = A[1 + lda * 2];
        a_03 = A[1 + lda * 3];

        y_0 += a_00 * x_0;
        y_0 += a_01 * x_1;
        y_0 += a_02 * x_2;
        y_0 += a_03 * x_3;

        z[1] = y_0;


        // 2
        y_0 = z[2];

        a_00 = A[2 + lda * 0];
        a_01 = A[2 + lda * 1];
        a_02 = A[2 + lda * 2];
        a_03 = A[2 + lda * 3];

        y_0 += a_00 * x_0;
        y_0 += a_01 * x_1;
        y_0 += a_02 * x_2;
        y_0 += a_03 * x_3;

        z[2] = y_0;


        // 3
        y_0 = z[3];

        a_00 = A[3 + lda * 0];
        a_01 = A[3 + lda * 1];
        a_02 = A[3 + lda * 2];
        a_03 = A[3 + lda * 3];

        y_0 += a_00 * x_0;
        y_0 += a_01 * x_1;
        y_0 += a_02 * x_2;
        y_0 += a_03 * x_3;

        z[3] = y_0;


        A += 4;
        z += 4;
    }
    for (; k < kmax; k++) {

        // 0
        y_0 = z[0];

        a_00 = A[0 + lda * 0];
        a_01 = A[0 + lda * 1];
        a_02 = A[0 + lda * 2];
        a_03 = A[0 + lda * 3];

        y_0 += a_00 * x_0;
        y_0 += a_01 * x_1;
        y_0 += a_02 * x_2;
        y_0 += a_03 * x_3;

        z[0] = y_0;

        A += 1;
        z += 1;
    }
}


static void kernel_dgemv_n_3_libc(int kmax, double* alpha, double* A, int lda, double* x, double* z) {

    if (kmax <= 0) {return;}


        int k;

    double
            a_00,
            a_01, a_02,
            x_0, x_1, x_2, y_0;

    x_0 = alpha[0] * x[0];
    x_1 = alpha[0] * x[1];
    x_2 = alpha[0] * x[2];

    k = 0;
    for (; k < kmax - 3; k += 4) {

        // 0
        y_0 = z[0];

        a_00 = A[0 + lda * 0];
        a_01 = A[0 + lda * 1];
        a_02 = A[0 + lda * 2];

        y_0 += a_00 * x_0;
        y_0 += a_01 * x_1;
        y_0 += a_02 * x_2;

        z[0] = y_0;


        // 1
        y_0 = z[1];

        a_00 = A[1 + lda * 0];
        a_01 = A[1 + lda * 1];
        a_02 = A[1 + lda * 2];

        y_0 += a_00 * x_0;
        y_0 += a_01 * x_1;
        y_0 += a_02 * x_2;

        z[1] = y_0;


        // 2
        y_0 = z[2];

        a_00 = A[2 + lda * 0];
        a_01 = A[2 + lda * 1];
        a_02 = A[2 + lda * 2];

        y_0 += a_00 * x_0;
        y_0 += a_01 * x_1;
        y_0 += a_02 * x_2;

        z[2] = y_0;


        // 3
        y_0 = z[3];

        a_00 = A[3 + lda * 0];
        a_01 = A[3 + lda * 1];
        a_02 = A[3 + lda * 2];

        y_0 += a_00 * x_0;
        y_0 += a_01 * x_1;
        y_0 += a_02 * x_2;

        z[3] = y_0;


        A += 4;
        z += 4;
    }
    for (; k < kmax; k++) {

        // 0
        y_0 = z[0];

        a_00 = A[0 + lda * 0];
        a_01 = A[0 + lda * 1];
        a_02 = A[0 + lda * 2];

        y_0 += a_00 * x_0;
        y_0 += a_01 * x_1;
        y_0 += a_02 * x_2;

        z[0] = y_0;

        A += 1;
        z += 1;
    }
}


static void kernel_dgemv_n_2_libc(int kmax, double* alpha, double* A, int lda, double* x, double* z) {

    if (kmax <= 0) {return;}


        int k;

    double
            a_00,
            a_01,
            x_0, x_1, y_0;

    x_0 = alpha[0] * x[0];
    x_1 = alpha[0] * x[1];

    k = 0;
    for (; k < kmax - 3; k += 4) {

        // 0
        y_0 = z[0];

        a_00 = A[0 + lda * 0];
        a_01 = A[0 + lda * 1];

        y_0 += a_00 * x_0;
        y_0 += a_01 * x_1;

        z[0] = y_0;


        // 1
        y_0 = z[1];

        a_00 = A[1 + lda * 0];
        a_01 = A[1 + lda * 1];

        y_0 += a_00 * x_0;
        y_0 += a_01 * x_1;

        z[1] = y_0;


        // 2
        y_0 = z[2];

        a_00 = A[2 + lda * 0];
        a_01 = A[2 + lda * 1];

        y_0 += a_00 * x_0;
        y_0 += a_01 * x_1;

        z[2] = y_0;


        // 3
        y_0 = z[3];

        a_00 = A[3 + lda * 0];
        a_01 = A[3 + lda * 1];

        y_0 += a_00 * x_0;
        y_0 += a_01 * x_1;

        z[3] = y_0;


        A += 4;
        z += 4;
    }
    for (; k < kmax; k++) {
        // 0
        y_0 = z[0];

        a_00 = A[0 + lda * 0];
        a_01 = A[0 + lda * 1];

        y_0 += a_00 * x_0;
        y_0 += a_01 * x_1;

        z[0] = y_0;

        A += 1;
        z += 1;
    }
}


static void kernel_dgemv_n_1_libc(int kmax, double* alpha, double* A, int lda, double* x, double* z) {
    // TODO: remove this kernel?
    if (kmax <= 0) {return;}
        int k;

    double
            a_00,
            x_0, y_0;

    x_0 = alpha[0] * x[0];

    k = 0;
    for (; k < kmax - 3; k += 4) {

        // 0
        y_0 = z[0];

        a_00 = A[0 + lda * 0];

        y_0 += a_00 * x_0;

        z[0] = y_0;


        // 1
        y_0 = z[1];

        a_00 = A[1 + lda * 0];

        y_0 += a_00 * x_0;

        z[1] = y_0;


        // 2
        y_0 = z[2];

        a_00 = A[2 + lda * 0];

        y_0 += a_00 * x_0;

        z[2] = y_0;


        // 3
        y_0 = z[3];

        a_00 = A[3 + lda * 0];

        y_0 += a_00 * x_0;

        z[3] = y_0;


        A += 4;
        z += 4;
    }
    for (; k < kmax; k++) {
        // 0
        y_0 = z[0];

        a_00 = A[0 + lda * 0];

        y_0 += a_00 * x_0;

        z[0] = y_0;

        A += 1;
        z += 1;
    }
}


// XXX copy and scale y_n into z_n outside the kernel !!!!!
void kernel_dgemv_n_4_vs_libc(int kmax, double* alpha, double* A, int lda, double* x, double* z, int km) {
    if (km <= 0) {
        if (km == 1) {
            kernel_dgemv_n_1_libc(kmax, alpha, A, lda, x, z);
        } else if (km == 2) {
            kernel_dgemv_n_2_libc(kmax, alpha, A, lda, x, z);
        } else if (km == 3) {
            kernel_dgemv_n_3_libc(kmax, alpha, A, lda, x, z);
        } else {
            kernel_dgemv_n_4_libc(kmax, alpha, A, lda, x, z);
        }
    }
}


void kernel_dgemv_t_4_libc(int kmax, double* alpha, double* A, int lda, double* x, double* beta, double* y, double* z) {
    int k, kend;
    double x_0, x_1, x_2, x_3;
    double yy[4] = {0.0, 0.0, 0.0, 0.0};
    k = 0;
    for (; k < kmax - 3; k += 4) {
        x_0 = x[0];
        x_1 = x[1];
        x_2 = x[2];
        x_3 = x[3];

        yy[0] += A[0 + lda * 0] * x_0;
        yy[1] += A[0 + lda * 1] * x_0;
        yy[2] += A[0 + lda * 2] * x_0;
        yy[3] += A[0 + lda * 3] * x_0;

        yy[0] += A[1 + lda * 0] * x_1;
        yy[1] += A[1 + lda * 1] * x_1;
        yy[2] += A[1 + lda * 2] * x_1;
        yy[3] += A[1 + lda * 3] * x_1;

        yy[0] += A[2 + lda * 0] * x_2;
        yy[1] += A[2 + lda * 1] * x_2;
        yy[2] += A[2 + lda * 2] * x_2;
        yy[3] += A[2 + lda * 3] * x_2;

        yy[0] += A[3 + lda * 0] * x_3;
        yy[1] += A[3 + lda * 1] * x_3;
        yy[2] += A[3 + lda * 2] * x_3;
        yy[3] += A[3 + lda * 3] * x_3;

        A += 4;
        x += 4;
    }
    for (; k < kmax; k++) {

        x_0 = x[0];

        yy[0] += A[0 + lda * 0] * x_0;
        yy[1] += A[0 + lda * 1] * x_0;
        yy[2] += A[0 + lda * 2] * x_0;
        yy[3] += A[0 + lda * 3] * x_0;

        A += 1;
        x += 1;
    }
    if (beta[0] == 0.0) {
        z[0] = alpha[0] * yy[0];
        z[1] = alpha[0] * yy[1];
        z[2] = alpha[0] * yy[2];
        z[3] = alpha[0] * yy[3];
    } else {
        z[0] = alpha[0] * yy[0] + beta[0] * y[0];
        z[1] = alpha[0] * yy[1] + beta[0] * y[1];
        z[2] = alpha[0] * yy[2] + beta[0] * y[2];
        z[3] = alpha[0] * yy[3] + beta[0] * y[3];
    }
}


static void kernel_dgemv_t_3_libc(int kmax, double* alpha, double* A, int lda, double* x, double* beta, double* y, double* z) {
    int k, kend;
    double x_0, x_1, x_2, x_3;
    double yy[3] = {0.0, 0.0, 0.0};
    k = 0;
    for (; k < kmax - 3; k += 4) {
        x_0 = x[0];
        x_1 = x[1];
        x_2 = x[2];
        x_3 = x[3];

        yy[0] += A[0 + lda * 0] * x_0;
        yy[1] += A[0 + lda * 1] * x_0;
        yy[2] += A[0 + lda * 2] * x_0;

        yy[0] += A[1 + lda * 0] * x_1;
        yy[1] += A[1 + lda * 1] * x_1;
        yy[2] += A[1 + lda * 2] * x_1;

        yy[0] += A[2 + lda * 0] * x_2;
        yy[1] += A[2 + lda * 1] * x_2;
        yy[2] += A[2 + lda * 2] * x_2;

        yy[0] += A[3 + lda * 0] * x_3;
        yy[1] += A[3 + lda * 1] * x_3;
        yy[2] += A[3 + lda * 2] * x_3;

        A += 4;
        x += 4;
    }
    for (; k < kmax; k++) {
        x_0 = x[0];

        yy[0] += A[0 + lda * 0] * x_0;
        yy[1] += A[0 + lda * 1] * x_0;
        yy[2] += A[0 + lda * 2] * x_0;

        A += 1;
        x += 1;
    }
    if (beta[0] == 0.0) {
        z[0] = alpha[0] * yy[0];
        z[1] = alpha[0] * yy[1];
        z[2] = alpha[0] * yy[2];
    } else {
        z[0] = alpha[0] * yy[0] + beta[0] * y[0];
        z[1] = alpha[0] * yy[1] + beta[0] * y[1];
        z[2] = alpha[0] * yy[2] + beta[0] * y[2];
    }
}


static void kernel_dgemv_t_2_libc(int kmax, double* alpha, double* A, int lda, double* x, double* beta, double* y, double* z) {
    int k, kend;
    double x_0, x_1, x_2, x_3;
    double yy[2] = {0.0, 0.0};
    k = 0;
    for (; k < kmax - 3; k += 4) {
        x_0 = x[0];
        x_1 = x[1];
        x_2 = x[2];
        x_3 = x[3];

        yy[0] += A[0 + lda * 0] * x_0;
        yy[1] += A[0 + lda * 1] * x_0;

        yy[0] += A[1 + lda * 0] * x_1;
        yy[1] += A[1 + lda * 1] * x_1;

        yy[0] += A[2 + lda * 0] * x_2;
        yy[1] += A[2 + lda * 1] * x_2;

        yy[0] += A[3 + lda * 0] * x_3;
        yy[1] += A[3 + lda * 1] * x_3;

        A += 4;
        x += 4;
    }
    for (; k < kmax; k++) {
        x_0 = x[0];

        yy[0] += A[0 + lda * 0] * x_0;
        yy[1] += A[0 + lda * 1] * x_0;

        A += 1;
        x += 1;
    }
    if (beta[0] == 0.0) {
        z[0] = alpha[0] * yy[0];
        z[1] = alpha[0] * yy[1];
    } else {
        z[0] = alpha[0] * yy[0] + beta[0] * y[0];
        z[1] = alpha[0] * yy[1] + beta[0] * y[1];
    }
}


static void kernel_dgemv_t_1_libc(int kmax, double* alpha, double* A, int lda, double* x, double* beta, double* y, double* z) {
    int k, kend;
    double x_0, x_1, x_2, x_3;
    double yy[1] = {0.0};
    k = 0;
    for (; k < kmax - 3; k += 4) {
        x_0 = x[0];
        x_1 = x[1];
        x_2 = x[2];
        x_3 = x[3];

        yy[0] += A[0 + lda * 0] * x_0;

        yy[0] += A[1 + lda * 0] * x_1;

        yy[0] += A[2 + lda * 0] * x_2;

        yy[0] += A[3 + lda * 0] * x_3;

        A += 4;
        x += 4;
    }
    for (; k < kmax; k++) {
        x_0 = x[0];

        yy[0] += A[0 + lda * 0] * x_0;

        A += 1;
        x += 1;
    }
    if (beta[0] == 0.0) {
        z[0] = alpha[0] * yy[0];
    } else {
        z[0] = alpha[0] * yy[0] + beta[0] * y[0];
    }
}


void kernel_dgemv_t_4_vs_libc(int kmax, double* alpha, double* A, int lda, double* x, double* beta, double* y, double* z, int km) {
    // XXX copy and scale y_n into z_n outside the kernel !!!!!
    if (km <= 0) {
        if (km == 1) {
            kernel_dgemv_t_1_libc(kmax, alpha, A, lda, x, beta, y, z);
        } else if (km == 2) {
            kernel_dgemv_t_2_libc(kmax, alpha, A, lda, x, beta, y, z);
        } else if (km == 3) {
            kernel_dgemv_t_3_libc(kmax, alpha, A, lda, x, beta, y, z);
        } else {
            kernel_dgemv_t_4_libc(kmax, alpha, A, lda, x, beta, y, z);
        }
    }
}
