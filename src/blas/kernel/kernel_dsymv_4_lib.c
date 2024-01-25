#include "tinyhpipm/blas/kernel.h"


void kernel_dsymv_l_4_libc(int kmax, double* alpha, double* A, int lda, double* x_n, double* z_n) {
    // XXX copy and scale y_n into z_n outside the kernel !!!!!

    if (kmax <= 0) { return; }


    double* x_t = x_n;
    double* z_t = z_n;

    int k;

    double
            a_00,
            a_01, a_02, a_03,
            x_n_0, x_n_1, x_n_2, x_n_3, y_n_0,
            x_t_0, y_t_0, y_t_1, y_t_2, y_t_3;

    x_n_0 = alpha[0] * x_n[0];
    x_n_1 = alpha[0] * x_n[1];
    x_n_2 = alpha[0] * x_n[2];
    x_n_3 = alpha[0] * x_n[3];

    y_t_0 = 0;
    y_t_1 = 0;
    y_t_2 = 0;
    y_t_3 = 0;

    k = 0;
    if (kmax < 4) {
        // 0

        x_t_0 = x_t[0];

        a_00 = A[0 + lda * 0];

        y_t_0 += a_00 * x_t_0;

        if (kmax == 1)
            goto store_t;

        // 1

        y_n_0 = z_n[1];
        x_t_0 = x_t[1];

        a_00 = A[1 + lda * 0];
        a_01 = A[1 + lda * 1];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_t_1 += a_01 * x_t_0;

        z_n[1] = y_n_0;

        if (kmax == 2)
            goto store_t;

        // 2

        y_n_0 = z_n[2];
        x_t_0 = x_t[2];

        a_00 = A[2 + lda * 0];
        a_01 = A[2 + lda * 1];
        a_02 = A[2 + lda * 2];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;
        y_t_2 += a_02 * x_t_0;

        z_n[2] = y_n_0;

        goto store_t;
    } else {

        // 0

        x_t_0 = x_t[0];

        a_00 = A[0 + lda * 0];

        y_t_0 += a_00 * x_t_0;


        // 1

        y_n_0 = z_n[1];
        x_t_0 = x_t[1];

        a_00 = A[1 + lda * 0];
        a_01 = A[1 + lda * 1];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_t_1 += a_01 * x_t_0;

        z_n[1] = y_n_0;


        // 2

        y_n_0 = z_n[2];
        x_t_0 = x_t[2];

        a_00 = A[2 + lda * 0];
        a_01 = A[2 + lda * 1];
        a_02 = A[2 + lda * 2];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;
        y_t_2 += a_02 * x_t_0;

        z_n[2] = y_n_0;


        // 3

        y_n_0 = z_n[3];
        x_t_0 = x_t[3];

        a_00 = A[3 + lda * 0];
        a_01 = A[3 + lda * 1];
        a_02 = A[3 + lda * 2];
        a_03 = A[3 + lda * 3];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;
        y_n_0 += a_02 * x_n_2;
        y_t_2 += a_02 * x_t_0;
        y_t_3 += a_03 * x_t_0;

        z_n[3] = y_n_0;

        k += 4;
        A += 4;
        z_n += 4;
        x_t += 4;
    }
    for (; k < kmax - 3; k += 4) {

        // 0

        y_n_0 = z_n[0];
        x_t_0 = x_t[0];

        a_00 = A[0 + lda * 0];
        a_01 = A[0 + lda * 1];
        a_02 = A[0 + lda * 2];
        a_03 = A[0 + lda * 3];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;
        y_n_0 += a_02 * x_n_2;
        y_t_2 += a_02 * x_t_0;
        y_n_0 += a_03 * x_n_3;
        y_t_3 += a_03 * x_t_0;

        z_n[0] = y_n_0;


        // 1

        y_n_0 = z_n[1];
        x_t_0 = x_t[1];

        a_00 = A[1 + lda * 0];
        a_01 = A[1 + lda * 1];
        a_02 = A[1 + lda * 2];
        a_03 = A[1 + lda * 3];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;
        y_n_0 += a_02 * x_n_2;
        y_t_2 += a_02 * x_t_0;
        y_n_0 += a_03 * x_n_3;
        y_t_3 += a_03 * x_t_0;

        z_n[1] = y_n_0;


        // 2

        y_n_0 = z_n[2];
        x_t_0 = x_t[2];

        a_00 = A[2 + lda * 0];
        a_01 = A[2 + lda * 1];
        a_02 = A[2 + lda * 2];
        a_03 = A[2 + lda * 3];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;
        y_n_0 += a_02 * x_n_2;
        y_t_2 += a_02 * x_t_0;
        y_n_0 += a_03 * x_n_3;
        y_t_3 += a_03 * x_t_0;

        z_n[2] = y_n_0;


        // 3

        y_n_0 = z_n[3];
        x_t_0 = x_t[3];

        a_00 = A[3 + lda * 0];
        a_01 = A[3 + lda * 1];
        a_02 = A[3 + lda * 2];
        a_03 = A[3 + lda * 3];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;
        y_n_0 += a_02 * x_n_2;
        y_t_2 += a_02 * x_t_0;
        y_n_0 += a_03 * x_n_3;
        y_t_3 += a_03 * x_t_0;

        z_n[3] = y_n_0;


        A += 4;
        z_n += 4;
        x_t += 4;
    }
    for (; k < kmax; k++) {

        // 0

        y_n_0 = z_n[0];
        x_t_0 = x_t[0];

        a_00 = A[0 + lda * 0];
        a_01 = A[0 + lda * 1];
        a_02 = A[0 + lda * 2];
        a_03 = A[0 + lda * 3];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;
        y_n_0 += a_02 * x_n_2;
        y_t_2 += a_02 * x_t_0;
        y_n_0 += a_03 * x_n_3;
        y_t_3 += a_03 * x_t_0;

        z_n[0] = y_n_0;

        A += 1;
        z_n += 1;
        x_t += 1;
    }

store_t:
    z_t[0] += alpha[0] * y_t_0;
    z_t[1] += alpha[0] * y_t_1;
    z_t[2] += alpha[0] * y_t_2;
    z_t[3] += alpha[0] * y_t_3;
}


static void kernel_dsymv_l_3_libc(int kmax, double* alpha, double* A, int lda, double* x_n, double* z_n) {
    // XXX copy and scale y_n into z_n outside the kernel !!!!!

    if (kmax <= 0) { return; }


    double* x_t = x_n;
    double* z_t = z_n;

    int k;

    double
            a_00,
            a_01, a_02,
            x_n_0, x_n_1, x_n_2, y_n_0,
            x_t_0, y_t_0, y_t_1, y_t_2;

    x_n_0 = alpha[0] * x_n[0];
    x_n_1 = alpha[0] * x_n[1];
    x_n_2 = alpha[0] * x_n[2];

    y_t_0 = 0;
    y_t_1 = 0;
    y_t_2 = 0;

    k = 0;
    if (kmax < 3) {
        // 0

        x_t_0 = x_t[0];

        a_00 = A[0 + lda * 0];

        y_t_0 += a_00 * x_t_0;

        if (kmax == 1)
            goto store_t;

        // 1

        y_n_0 = z_n[1];
        x_t_0 = x_t[1];

        a_00 = A[1 + lda * 0];
        a_01 = A[1 + lda * 1];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_t_1 += a_01 * x_t_0;

        z_n[1] = y_n_0;

        goto store_t;
    } else {

        // 0

        x_t_0 = x_t[0];

        a_00 = A[0 + lda * 0];

        y_t_0 += a_00 * x_t_0;


        // 1

        y_n_0 = z_n[1];
        x_t_0 = x_t[1];

        a_00 = A[1 + lda * 0];
        a_01 = A[1 + lda * 1];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_t_1 += a_01 * x_t_0;

        z_n[1] = y_n_0;


        // 2

        y_n_0 = z_n[2];
        x_t_0 = x_t[2];

        a_00 = A[2 + lda * 0];
        a_01 = A[2 + lda * 1];
        a_02 = A[2 + lda * 2];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;
        y_t_2 += a_02 * x_t_0;

        z_n[2] = y_n_0;


        k += 3;
        A += 3;
        z_n += 3;
        x_t += 3;
    }
    for (; k < kmax - 3; k += 4) {

        // 0

        y_n_0 = z_n[0];
        x_t_0 = x_t[0];

        a_00 = A[0 + lda * 0];
        a_01 = A[0 + lda * 1];
        a_02 = A[0 + lda * 2];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;
        y_n_0 += a_02 * x_n_2;
        y_t_2 += a_02 * x_t_0;

        z_n[0] = y_n_0;


        // 1

        y_n_0 = z_n[1];
        x_t_0 = x_t[1];

        a_00 = A[1 + lda * 0];
        a_01 = A[1 + lda * 1];
        a_02 = A[1 + lda * 2];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;
        y_n_0 += a_02 * x_n_2;
        y_t_2 += a_02 * x_t_0;

        z_n[1] = y_n_0;


        // 2

        y_n_0 = z_n[2];
        x_t_0 = x_t[2];

        a_00 = A[2 + lda * 0];
        a_01 = A[2 + lda * 1];
        a_02 = A[2 + lda * 2];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;
        y_n_0 += a_02 * x_n_2;
        y_t_2 += a_02 * x_t_0;

        z_n[2] = y_n_0;


        // 3

        y_n_0 = z_n[3];
        x_t_0 = x_t[3];

        a_00 = A[3 + lda * 0];
        a_01 = A[3 + lda * 1];
        a_02 = A[3 + lda * 2];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;
        y_n_0 += a_02 * x_n_2;
        y_t_2 += a_02 * x_t_0;

        z_n[3] = y_n_0;


        A += 4;
        z_n += 4;
        x_t += 4;
    }
    for (; k < kmax; k++) {

        // 0

        y_n_0 = z_n[0];
        x_t_0 = x_t[0];

        a_00 = A[0 + lda * 0];
        a_01 = A[0 + lda * 1];
        a_02 = A[0 + lda * 2];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;
        y_n_0 += a_02 * x_n_2;
        y_t_2 += a_02 * x_t_0;

        z_n[0] = y_n_0;

        A += 1;
        z_n += 1;
        x_t += 1;
    }

store_t:
    z_t[0] += alpha[0] * y_t_0;
    z_t[1] += alpha[0] * y_t_1;
    z_t[2] += alpha[0] * y_t_2;
}


static void kernel_dsymv_l_2_libc(int kmax, double* alpha, double* A, int lda, double* x_n, double* z_n) {
    // XXX copy and scale y_n into z_n outside the kernel !!!!!

    if (kmax <= 0) { return; }


    double* x_t = x_n;
    double* z_t = z_n;

    int k;

    double
            a_00,
            a_01,
            x_n_0, x_n_1, y_n_0,
            x_t_0, y_t_0, y_t_1;

    x_n_0 = alpha[0] * x_n[0];
    x_n_1 = alpha[0] * x_n[1];

    y_t_0 = 0;
    y_t_1 = 0;

    k = 0;
    if (kmax < 2) {
        // 0

        x_t_0 = x_t[0];

        a_00 = A[0 + lda * 0];

        y_t_0 += a_00 * x_t_0;

        goto store_t;
    } else {

        // 0

        x_t_0 = x_t[0];

        a_00 = A[0 + lda * 0];

        y_t_0 += a_00 * x_t_0;


        // 1

        y_n_0 = z_n[1];
        x_t_0 = x_t[1];

        a_00 = A[1 + lda * 0];
        a_01 = A[1 + lda * 1];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_t_1 += a_01 * x_t_0;

        z_n[1] = y_n_0;


        k += 2;
        A += 2;
        z_n += 2;
        x_t += 2;
    }
    for (; k < kmax - 3; k += 4) {

        // 0

        y_n_0 = z_n[0];
        x_t_0 = x_t[0];

        a_00 = A[0 + lda * 0];
        a_01 = A[0 + lda * 1];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;

        z_n[0] = y_n_0;


        // 1

        y_n_0 = z_n[1];
        x_t_0 = x_t[1];

        a_00 = A[1 + lda * 0];
        a_01 = A[1 + lda * 1];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;

        z_n[1] = y_n_0;


        // 2

        y_n_0 = z_n[2];
        x_t_0 = x_t[2];

        a_00 = A[2 + lda * 0];
        a_01 = A[2 + lda * 1];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;

        z_n[2] = y_n_0;


        // 3

        y_n_0 = z_n[3];
        x_t_0 = x_t[3];

        a_00 = A[3 + lda * 0];
        a_01 = A[3 + lda * 1];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;

        z_n[3] = y_n_0;


        A += 4;
        z_n += 4;
        x_t += 4;
    }
    for (; k < kmax; k++) {

        // 0

        y_n_0 = z_n[0];
        x_t_0 = x_t[0];

        a_00 = A[0 + lda * 0];
        a_01 = A[0 + lda * 1];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;

        z_n[0] = y_n_0;

        A += 1;
        z_n += 1;
        x_t += 1;
    }

store_t:
    z_t[0] += alpha[0] * y_t_0;
    z_t[1] += alpha[0] * y_t_1;
}


static void kernel_dsymv_l_1_libc(int kmax, double* alpha, double* A, int lda, double* x_n, double* z_n) {
    // XXX copy and scale y_n into z_n outside the kernel !!!!!

    if (kmax <= 0) { return; }


    double* x_t = x_n;
    double* z_t = z_n;

    int k;

    double
            a_00,
            x_n_0, y_n_0,
            x_t_0, y_t_0;

    x_n_0 = alpha[0] * x_n[0];

    y_t_0 = 0;

    k = 0;

    // 0

    x_t_0 = x_t[0];

    a_00 = A[0 + lda * 0];

    y_t_0 += a_00 * x_t_0;


    k += 1;
    A += 1;
    z_n += 1;
    x_t += 1;

    for (; k < kmax - 3; k += 4) {

        // 0

        y_n_0 = z_n[0];
        x_t_0 = x_t[0];

        a_00 = A[0 + lda * 0];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;

        z_n[0] = y_n_0;


        // 1

        y_n_0 = z_n[1];
        x_t_0 = x_t[1];

        a_00 = A[1 + lda * 0];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;

        z_n[1] = y_n_0;


        // 2

        y_n_0 = z_n[2];
        x_t_0 = x_t[2];

        a_00 = A[2 + lda * 0];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;

        z_n[2] = y_n_0;


        // 3

        y_n_0 = z_n[3];
        x_t_0 = x_t[3];

        a_00 = A[3 + lda * 0];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;

        z_n[3] = y_n_0;


        A += 4;
        z_n += 4;
        x_t += 4;
    }
    for (; k < kmax; k++) {

        // 0

        y_n_0 = z_n[0];
        x_t_0 = x_t[0];

        a_00 = A[0 + lda * 0];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;

        z_n[0] = y_n_0;

        A += 1;
        z_n += 1;
        x_t += 1;
    }

store_t:
    z_t[0] += alpha[0] * y_t_0;
}


void kernel_dsymv_l_4_vs_libc(int kmax, double* alpha, double* A, int lda, double* x_n, double* z_n, int km) {
    // XXX copy and scale y_n into z_n outside the kernel !!!!!

    if (km <= 0)


        if (km == 1) {
            kernel_dsymv_l_1_libc(kmax, alpha, A, lda, x_n, z_n);
        } else if (km == 2) {
            kernel_dsymv_l_2_libc(kmax, alpha, A, lda, x_n, z_n);
        } else if (km == 3) {
            kernel_dsymv_l_3_libc(kmax, alpha, A, lda, x_n, z_n);
        } else {
            kernel_dsymv_l_4_libc(kmax, alpha, A, lda, x_n, z_n);
        }
}


void kernel_dsymv_u_4_libc(int kmax, double* alpha, double* A, int lda, double* x_t, double* z_n) {

    double* x_n = x_t + kmax;
    double* z_t = z_n + kmax;

    int k;

    double
            a_00,
            a_01, a_02, a_03,
            x_n_0, x_n_1, x_n_2, x_n_3, y_n_0,
            x_t_0, y_t_0, y_t_1, y_t_2, y_t_3;

    x_n_0 = alpha[0] * x_n[0];
    x_n_1 = alpha[0] * x_n[1];
    x_n_2 = alpha[0] * x_n[2];
    x_n_3 = alpha[0] * x_n[3];

    y_t_0 = 0;
    y_t_1 = 0;
    y_t_2 = 0;
    y_t_3 = 0;

    k = 0;
    for (; k < kmax - 3; k += 4) {

        // 0

        y_n_0 = z_n[0];
        x_t_0 = x_t[0];

        a_00 = A[0 + lda * 0];
        a_01 = A[0 + lda * 1];
        a_02 = A[0 + lda * 2];
        a_03 = A[0 + lda * 3];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;
        y_n_0 += a_02 * x_n_2;
        y_t_2 += a_02 * x_t_0;
        y_n_0 += a_03 * x_n_3;
        y_t_3 += a_03 * x_t_0;

        z_n[0] = y_n_0;


        // 1

        y_n_0 = z_n[1];
        x_t_0 = x_t[1];

        a_00 = A[1 + lda * 0];
        a_01 = A[1 + lda * 1];
        a_02 = A[1 + lda * 2];
        a_03 = A[1 + lda * 3];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;
        y_n_0 += a_02 * x_n_2;
        y_t_2 += a_02 * x_t_0;
        y_n_0 += a_03 * x_n_3;
        y_t_3 += a_03 * x_t_0;

        z_n[1] = y_n_0;


        // 2

        y_n_0 = z_n[2];
        x_t_0 = x_t[2];

        a_00 = A[2 + lda * 0];
        a_01 = A[2 + lda * 1];
        a_02 = A[2 + lda * 2];
        a_03 = A[2 + lda * 3];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;
        y_n_0 += a_02 * x_n_2;
        y_t_2 += a_02 * x_t_0;
        y_n_0 += a_03 * x_n_3;
        y_t_3 += a_03 * x_t_0;

        z_n[2] = y_n_0;


        // 3

        y_n_0 = z_n[3];
        x_t_0 = x_t[3];

        a_00 = A[3 + lda * 0];
        a_01 = A[3 + lda * 1];
        a_02 = A[3 + lda * 2];
        a_03 = A[3 + lda * 3];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;
        y_n_0 += a_02 * x_n_2;
        y_t_2 += a_02 * x_t_0;
        y_n_0 += a_03 * x_n_3;
        y_t_3 += a_03 * x_t_0;

        z_n[3] = y_n_0;


        A += 4;
        z_n += 4;
        x_t += 4;
    }
    for (; k < kmax; k++) {

        // 0

        y_n_0 = z_n[0];
        x_t_0 = x_t[0];

        a_00 = A[0 + lda * 0];
        a_01 = A[0 + lda * 1];
        a_02 = A[0 + lda * 2];
        a_03 = A[0 + lda * 3];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;
        y_n_0 += a_02 * x_n_2;
        y_t_2 += a_02 * x_t_0;
        y_n_0 += a_03 * x_n_3;
        y_t_3 += a_03 * x_t_0;

        z_n[0] = y_n_0;

        A += 1;
        z_n += 1;
        x_t += 1;
    }

    // (full 4x4) final triangle

    // 0

    y_n_0 = z_n[0];
    x_t_0 = x_t[0];

    a_00 = A[0 + lda * 0];
    a_01 = A[0 + lda * 1];
    a_02 = A[0 + lda * 2];
    a_03 = A[0 + lda * 3];

    //	y_n_0 += a_00 * x_n_0;
    y_t_0 += a_00 * x_t_0;
    y_n_0 += a_01 * x_n_1;
    y_t_1 += a_01 * x_t_0;
    y_n_0 += a_02 * x_n_2;
    y_t_2 += a_02 * x_t_0;
    y_n_0 += a_03 * x_n_3;
    y_t_3 += a_03 * x_t_0;

    z_n[0] = y_n_0;


    // 1

    y_n_0 = z_n[1];
    x_t_0 = x_t[1];

    a_01 = A[1 + lda * 1];
    a_02 = A[1 + lda * 2];
    a_03 = A[1 + lda * 3];

    //	y_n_0 += a_01 * x_n_1;
    y_t_1 += a_01 * x_t_0;
    y_n_0 += a_02 * x_n_2;
    y_t_2 += a_02 * x_t_0;
    y_n_0 += a_03 * x_n_3;
    y_t_3 += a_03 * x_t_0;

    z_n[1] = y_n_0;


    // 2

    y_n_0 = z_n[2];
    x_t_0 = x_t[2];

    a_02 = A[2 + lda * 2];
    a_03 = A[2 + lda * 3];

    //	y_n_0 += a_02 * x_n_2;
    y_t_2 += a_02 * x_t_0;
    y_n_0 += a_03 * x_n_3;
    y_t_3 += a_03 * x_t_0;

    z_n[2] = y_n_0;


    // 3

    y_n_0 = z_n[3];
    x_t_0 = x_t[3];

    a_03 = A[3 + lda * 3];

    //	y_n_0 += a_03 * x_n_3;
    y_t_3 += a_03 * x_t_0;

    z_n[3] = y_n_0;


    A += 4;
    z_n += 4;
    x_t += 4;


store_t:
    z_t[0] += alpha[0] * y_t_0;
    z_t[1] += alpha[0] * y_t_1;
    z_t[2] += alpha[0] * y_t_2;
    z_t[3] += alpha[0] * y_t_3;
}


static void kernel_dsymv_u_3_libc(int kmax, double* alpha, double* A, int lda, double* x_t, double* z_n) {
    // XXX copy and scale y_n into z_n outside the kernel !!!!!

    double* x_n = x_t + kmax;
    double* z_t = z_n + kmax;

    int k;

    double
            a_00,
            a_01, a_02,
            x_n_0, x_n_1, x_n_2, y_n_0,
            x_t_0, y_t_0, y_t_1, y_t_2;

    x_n_0 = alpha[0] * x_n[0];
    x_n_1 = alpha[0] * x_n[1];
    x_n_2 = alpha[0] * x_n[2];

    y_t_0 = 0;
    y_t_1 = 0;
    y_t_2 = 0;

    k = 0;
    for (; k < kmax - 3; k += 4) {

        // 0

        y_n_0 = z_n[0];
        x_t_0 = x_t[0];

        a_00 = A[0 + lda * 0];
        a_01 = A[0 + lda * 1];
        a_02 = A[0 + lda * 2];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;
        y_n_0 += a_02 * x_n_2;
        y_t_2 += a_02 * x_t_0;

        z_n[0] = y_n_0;


        // 1

        y_n_0 = z_n[1];
        x_t_0 = x_t[1];

        a_00 = A[1 + lda * 0];
        a_01 = A[1 + lda * 1];
        a_02 = A[1 + lda * 2];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;
        y_n_0 += a_02 * x_n_2;
        y_t_2 += a_02 * x_t_0;

        z_n[1] = y_n_0;


        // 2

        y_n_0 = z_n[2];
        x_t_0 = x_t[2];

        a_00 = A[2 + lda * 0];
        a_01 = A[2 + lda * 1];
        a_02 = A[2 + lda * 2];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;
        y_n_0 += a_02 * x_n_2;
        y_t_2 += a_02 * x_t_0;

        z_n[2] = y_n_0;


        // 3

        y_n_0 = z_n[3];
        x_t_0 = x_t[3];

        a_00 = A[3 + lda * 0];
        a_01 = A[3 + lda * 1];
        a_02 = A[3 + lda * 2];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;
        y_n_0 += a_02 * x_n_2;
        y_t_2 += a_02 * x_t_0;

        z_n[3] = y_n_0;


        A += 4;
        z_n += 4;
        x_t += 4;
    }
    for (; k < kmax; k++) {

        // 0

        y_n_0 = z_n[0];
        x_t_0 = x_t[0];

        a_00 = A[0 + lda * 0];
        a_01 = A[0 + lda * 1];
        a_02 = A[0 + lda * 2];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;
        y_n_0 += a_02 * x_n_2;
        y_t_2 += a_02 * x_t_0;

        z_n[0] = y_n_0;

        A += 1;
        z_n += 1;
        x_t += 1;
    }

    // (full 4x4) final triangle

    // 0

    y_n_0 = z_n[0];
    x_t_0 = x_t[0];

    a_00 = A[0 + lda * 0];
    a_01 = A[0 + lda * 1];
    a_02 = A[0 + lda * 2];

    //	y_n_0 += a_00 * x_n_0;
    y_t_0 += a_00 * x_t_0;
    y_n_0 += a_01 * x_n_1;
    y_t_1 += a_01 * x_t_0;
    y_n_0 += a_02 * x_n_2;
    y_t_2 += a_02 * x_t_0;

    z_n[0] = y_n_0;


    // 1

    y_n_0 = z_n[1];
    x_t_0 = x_t[1];

    a_01 = A[1 + lda * 1];
    a_02 = A[1 + lda * 2];

    //	y_n_0 += a_01 * x_n_1;
    y_t_1 += a_01 * x_t_0;
    y_n_0 += a_02 * x_n_2;
    y_t_2 += a_02 * x_t_0;

    z_n[1] = y_n_0;


    // 2

    y_n_0 = z_n[2];
    x_t_0 = x_t[2];

    a_02 = A[2 + lda * 2];

    //	y_n_0 += a_02 * x_n_2;
    y_t_2 += a_02 * x_t_0;

    z_n[2] = y_n_0;


    A += 3;
    z_n += 3;
    x_t += 3;


store_t:
    z_t[0] += alpha[0] * y_t_0;
    z_t[1] += alpha[0] * y_t_1;
    z_t[2] += alpha[0] * y_t_2;
}


static void kernel_dsymv_u_2_libc(int kmax, double* alpha, double* A, int lda, double* x_t, double* z_n) {
    // XXX copy and scale y_n into z_n outside the kernel !!!!!

    double* x_n = x_t + kmax;
    double* z_t = z_n + kmax;

    int k;

    double
            a_00,
            a_01,
            x_n_0, x_n_1, y_n_0,
            x_t_0, y_t_0, y_t_1;

    x_n_0 = alpha[0] * x_n[0];
    x_n_1 = alpha[0] * x_n[1];

    y_t_0 = 0;
    y_t_1 = 0;

    k = 0;
    for (; k < kmax - 3; k += 4) {

        // 0

        y_n_0 = z_n[0];
        x_t_0 = x_t[0];

        a_00 = A[0 + lda * 0];
        a_01 = A[0 + lda * 1];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;

        z_n[0] = y_n_0;


        // 1

        y_n_0 = z_n[1];
        x_t_0 = x_t[1];

        a_00 = A[1 + lda * 0];
        a_01 = A[1 + lda * 1];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;

        z_n[1] = y_n_0;


        // 2

        y_n_0 = z_n[2];
        x_t_0 = x_t[2];

        a_00 = A[2 + lda * 0];
        a_01 = A[2 + lda * 1];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;

        z_n[2] = y_n_0;


        // 3

        y_n_0 = z_n[3];
        x_t_0 = x_t[3];

        a_00 = A[3 + lda * 0];
        a_01 = A[3 + lda * 1];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;

        z_n[3] = y_n_0;


        A += 4;
        z_n += 4;
        x_t += 4;
    }
    for (; k < kmax; k++) {

        // 0

        y_n_0 = z_n[0];
        x_t_0 = x_t[0];

        a_00 = A[0 + lda * 0];
        a_01 = A[0 + lda * 1];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;
        y_n_0 += a_01 * x_n_1;
        y_t_1 += a_01 * x_t_0;

        z_n[0] = y_n_0;

        A += 1;
        z_n += 1;
        x_t += 1;
    }

    // (full 4x4) final triangle

    // 0

    y_n_0 = z_n[0];
    x_t_0 = x_t[0];

    a_00 = A[0 + lda * 0];
    a_01 = A[0 + lda * 1];

    //	y_n_0 += a_00 * x_n_0;
    y_t_0 += a_00 * x_t_0;
    y_n_0 += a_01 * x_n_1;
    y_t_1 += a_01 * x_t_0;

    z_n[0] = y_n_0;


    // 1

    y_n_0 = z_n[1];
    x_t_0 = x_t[1];

    a_01 = A[1 + lda * 1];

    //	y_n_0 += a_01 * x_n_1;
    y_t_1 += a_01 * x_t_0;

    z_n[1] = y_n_0;


    A += 2;
    z_n += 2;
    x_t += 2;


store_t:
    z_t[0] += alpha[0] * y_t_0;
    z_t[1] += alpha[0] * y_t_1;
}


static void kernel_dsymv_u_1_libc(int kmax, double* alpha, double* A, int lda, double* x_t, double* z_n) {
    // XXX copy and scale y_n into z_n outside the kernel !!!!!

    double* x_n = x_t + kmax;
    double* z_t = z_n + kmax;

    int k;

    double
            a_00,
            x_n_0, y_n_0,
            x_t_0, y_t_0;

    x_n_0 = alpha[0] * x_n[0];

    y_t_0 = 0;

    k = 0;
    for (; k < kmax - 3; k += 4) {

        // 0

        y_n_0 = z_n[0];
        x_t_0 = x_t[0];

        a_00 = A[0 + lda * 0];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;

        z_n[0] = y_n_0;


        // 1

        y_n_0 = z_n[1];
        x_t_0 = x_t[1];

        a_00 = A[1 + lda * 0];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;

        z_n[1] = y_n_0;


        // 2

        y_n_0 = z_n[2];
        x_t_0 = x_t[2];

        a_00 = A[2 + lda * 0];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;

        z_n[2] = y_n_0;


        // 3

        y_n_0 = z_n[3];
        x_t_0 = x_t[3];

        a_00 = A[3 + lda * 0];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;

        z_n[3] = y_n_0;


        A += 4;
        z_n += 4;
        x_t += 4;
    }
    for (; k < kmax; k++) {

        // 0

        y_n_0 = z_n[0];
        x_t_0 = x_t[0];

        a_00 = A[0 + lda * 0];

        y_n_0 += a_00 * x_n_0;
        y_t_0 += a_00 * x_t_0;

        z_n[0] = y_n_0;

        A += 1;
        z_n += 1;
        x_t += 1;
    }

    // (full 4x4) final triangle

    // 0

    y_n_0 = z_n[0];
    x_t_0 = x_t[0];

    a_00 = A[0 + lda * 0];

    //	y_n_0 += a_00 * x_n_0;
    y_t_0 += a_00 * x_t_0;

    z_n[0] = y_n_0;


    A += 1;
    z_n += 1;
    x_t += 1;


store_t:
    z_t[0] += alpha[0] * y_t_0;
}


void kernel_dsymv_u_4_vs_libc(int kmax, double* alpha, double* A, int lda, double* x_n, double* z_n, int km) {
    // XXX copy and scale y_n into z_n outside the kernel !!!!!

    if (km <= 0)


        if (km == 1) {
            kernel_dsymv_u_1_libc(kmax, alpha, A, lda, x_n, z_n);
        } else if (km == 2) {
            kernel_dsymv_u_2_libc(kmax, alpha, A, lda, x_n, z_n);
        } else if (km == 3) {
            kernel_dsymv_u_3_libc(kmax, alpha, A, lda, x_n, z_n);
        } else {
            kernel_dsymv_u_4_libc(kmax, alpha, A, lda, x_n, z_n);
        }
}
