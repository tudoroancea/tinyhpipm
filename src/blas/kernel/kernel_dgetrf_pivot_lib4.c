#include <math.h>
#include <stdio.h>

#include "tinyhpipm/blas/kernel.h"


// swap two rows
void kernel_drowsw_lib4(int kmax, double* pA, double* pC) {

    const int bs = 4;

    int ii;
    double tmp;

    for (ii = 0; ii < kmax - 3; ii += 4) {
        tmp = pA[0 + bs * 0];
        pA[0 + bs * 0] = pC[0 + bs * 0];
        pC[0 + bs * 0] = tmp;
        tmp = pA[0 + bs * 1];
        pA[0 + bs * 1] = pC[0 + bs * 1];
        pC[0 + bs * 1] = tmp;
        tmp = pA[0 + bs * 2];
        pA[0 + bs * 2] = pC[0 + bs * 2];
        pC[0 + bs * 2] = tmp;
        tmp = pA[0 + bs * 3];
        pA[0 + bs * 3] = pC[0 + bs * 3];
        pC[0 + bs * 3] = tmp;
        pA += 4 * bs;
        pC += 4 * bs;
    }
    for (; ii < kmax; ii++) {
        tmp = pA[0 + bs * 0];
        pA[0 + bs * 0] = pC[0 + bs * 0];
        pC[0 + bs * 0] = tmp;
        pA += 1 * bs;
        pC += 1 * bs;
    }
}


// C numbering, starting from 0
void didamax_lib4(int n, int offset, double* pA, int sda, int* p_idamax, double* p_amax) {

    int ii;
    double tmp, amax;

    int idamax = -1;

    p_idamax[0] = idamax;
    if (n < 1) { return; }


    const int bs = 4;

    int na = (bs - offset % bs) % bs;
    na = n < na ? n : na;

    amax = -1.0;
    ii = 0;
    if (na > 0) {
        for (; ii < na; ii++) {
            tmp = fabs(pA[0]);
            if (tmp > amax) {
                idamax = ii + 0;
                amax = tmp;
            }
            pA += 1;
        }
        pA += bs * (sda - 1);
    }
    for (; ii < n - 3; ii += 4) {
        tmp = fabs(pA[0]);
        if (tmp > amax) {
            idamax = ii + 0;
            amax = tmp;
        }
        tmp = fabs(pA[1]);
        if (tmp > amax) {
            idamax = ii + 1;
            amax = tmp;
        }
        tmp = fabs(pA[2]);
        if (tmp > amax) {
            idamax = ii + 2;
            amax = tmp;
        }
        tmp = fabs(pA[3]);
        if (tmp > amax) {
            idamax = ii + 3;
            amax = tmp;
        }
        pA += bs * sda;
    }
    for (; ii < n; ii++) {
        tmp = fabs(pA[0]);
        if (tmp > amax) {
            idamax = ii + 0;
            amax = tmp;
        }
        pA += 1;
    }

    p_amax[0] = amax;
    p_idamax[0] = idamax;
}


// C numering (starting from zero) in the ipiv
// it process m>=4 rows and 4 cols
void kernel_dgetrf_pivot_4_lib4(int m, double* pA, int sda, double* inv_diag_A, int* ipiv) {

    const int bs = 4;

    // assume m>=4
    int ma = m - 4;

    double
            tmp0,
            tmp1, tmp2, tmp3,
            u_00, u_01, u_02, u_03,
            u_11, u_12, u_13,
            u_22, u_23,
            u_33;

    double* pB;

    int
            k,
            idamax;

    // first column
    didamax_lib4(m - 0, 0, &pA[0 + bs * 0], sda, &idamax, &tmp0);
    ipiv[0] = idamax;
    if (tmp0 != 0.0) {
        if (ipiv[0] != 0)
            kernel_drowsw_lib4(4, pA + 0, pA + ipiv[0] / bs * bs * sda + ipiv[0] % bs);

        tmp0 = 1.0 / pA[0 + bs * 0];
        inv_diag_A[0] = tmp0;
        pA[1 + bs * 0] *= tmp0;
        pA[2 + bs * 0] *= tmp0;
        pA[3 + bs * 0] *= tmp0;
        pB = pA + bs * sda;
        for (k = 0; k < ma - 3; k += 4) {
            pB[0 + bs * 0] *= tmp0;
            pB[1 + bs * 0] *= tmp0;
            pB[2 + bs * 0] *= tmp0;
            pB[3 + bs * 0] *= tmp0;
            pB += bs * sda;
        }
        for (; k < ma; k++) {
            pB[0 + bs * 0] *= tmp0;
            pB += 1;
        }
    } else {
        inv_diag_A[0] = 0.0;
    }

    // second column
    u_01 = pA[0 + bs * 1];
    tmp1 = pA[1 + bs * 1];
    tmp2 = pA[2 + bs * 1];
    tmp3 = pA[3 + bs * 1];
    tmp1 -= pA[1 + bs * 0] * u_01;
    tmp2 -= pA[2 + bs * 0] * u_01;
    tmp3 -= pA[3 + bs * 0] * u_01;
    pA[1 + bs * 1] = tmp1;
    pA[2 + bs * 1] = tmp2;
    pA[3 + bs * 1] = tmp3;
    pB = pA + bs * sda;
    for (k = 0; k < ma - 3; k += 4) {
        tmp0 = pB[0 + bs * 1];
        tmp1 = pB[1 + bs * 1];
        tmp2 = pB[2 + bs * 1];
        tmp3 = pB[3 + bs * 1];
        tmp0 -= pB[0 + bs * 0] * u_01;
        tmp1 -= pB[1 + bs * 0] * u_01;
        tmp2 -= pB[2 + bs * 0] * u_01;
        tmp3 -= pB[3 + bs * 0] * u_01;
        pB[0 + bs * 1] = tmp0;
        pB[1 + bs * 1] = tmp1;
        pB[2 + bs * 1] = tmp2;
        pB[3 + bs * 1] = tmp3;
        pB += bs * sda;
    }
    for (; k < ma; k++) {
        tmp0 = pB[0 + bs * 1];
        tmp0 -= pB[0 + bs * 0] * u_01;
        pB[0 + bs * 1] = tmp0;
        pB += 1;
    }

    didamax_lib4(m - 1, 1, &pA[1 + bs * 1], sda, &idamax, &tmp1);
    ipiv[1] = idamax + 1;
    if (tmp1 != 0) {
        if (ipiv[1] != 1)
            kernel_drowsw_lib4(4, pA + 1, pA + ipiv[1] / bs * bs * sda + ipiv[1] % bs);

        tmp1 = 1.0 / pA[1 + bs * 1];
        inv_diag_A[1] = tmp1;
        pA[2 + bs * 1] *= tmp1;
        pA[3 + bs * 1] *= tmp1;
        pB = pA + bs * sda;
        for (k = 0; k < ma - 3; k += 4) {
            pB[0 + bs * 1] *= tmp1;
            pB[1 + bs * 1] *= tmp1;
            pB[2 + bs * 1] *= tmp1;
            pB[3 + bs * 1] *= tmp1;
            pB += bs * sda;
        }
        for (; k < ma; k++) {
            pB[0 + bs * 1] *= tmp1;
            pB += 1;
        }
    } else {
        inv_diag_A[1] = 0.0;
    }

    // third column
    u_02 = pA[0 + bs * 2];
    u_12 = pA[1 + bs * 2];
    u_12 -= pA[1 + bs * 0] * u_02;
    pA[1 + bs * 2] = u_12;
    tmp2 = pA[2 + bs * 2];
    tmp3 = pA[3 + bs * 2];
    tmp2 -= pA[2 + bs * 0] * u_02;
    tmp3 -= pA[3 + bs * 0] * u_02;
    tmp2 -= pA[2 + bs * 1] * u_12;
    tmp3 -= pA[3 + bs * 1] * u_12;
    pA[2 + bs * 2] = tmp2;
    pA[3 + bs * 2] = tmp3;
    pB = pA + bs * sda;
    for (k = 0; k < ma - 3; k += 4) {
        tmp0 = pB[0 + bs * 2];
        tmp1 = pB[1 + bs * 2];
        tmp2 = pB[2 + bs * 2];
        tmp3 = pB[3 + bs * 2];
        tmp0 -= pB[0 + bs * 0] * u_02;
        tmp1 -= pB[1 + bs * 0] * u_02;
        tmp2 -= pB[2 + bs * 0] * u_02;
        tmp3 -= pB[3 + bs * 0] * u_02;
        tmp0 -= pB[0 + bs * 1] * u_12;
        tmp1 -= pB[1 + bs * 1] * u_12;
        tmp2 -= pB[2 + bs * 1] * u_12;
        tmp3 -= pB[3 + bs * 1] * u_12;
        pB[0 + bs * 2] = tmp0;
        pB[1 + bs * 2] = tmp1;
        pB[2 + bs * 2] = tmp2;
        pB[3 + bs * 2] = tmp3;
        pB += bs * sda;
    }
    for (; k < ma; k++) {
        tmp0 = pB[0 + bs * 2];
        tmp0 -= pB[0 + bs * 0] * u_02;
        tmp0 -= pB[0 + bs * 1] * u_12;
        pB[0 + bs * 2] = tmp0;
        pB += 1;
    }

    didamax_lib4(m - 2, 2, &pA[2 + bs * 2], sda, &idamax, &tmp2);
    ipiv[2] = idamax + 2;
    if (tmp2 != 0) {
        if (ipiv[2] != 2)
            kernel_drowsw_lib4(4, pA + 2, pA + ipiv[2] / bs * bs * sda + ipiv[2] % bs);

        tmp2 = 1.0 / pA[2 + bs * 2];
        inv_diag_A[2] = tmp2;
        pA[3 + bs * 2] *= tmp2;
        pB = pA + bs * sda;
        for (k = 0; k < ma - 3; k += 4) {
            pB[0 + bs * 2] *= tmp2;
            pB[1 + bs * 2] *= tmp2;
            pB[2 + bs * 2] *= tmp2;
            pB[3 + bs * 2] *= tmp2;
            pB += bs * sda;
        }
        for (; k < ma; k++) {
            pB[0 + bs * 2] *= tmp2;
            pB += 1;
        }
    } else {
        inv_diag_A[2] = 0.0;
    }

    // fourth column
    u_03 = pA[0 + bs * 3];
    u_13 = pA[1 + bs * 3];
    u_13 -= pA[1 + bs * 0] * u_03;
    pA[1 + bs * 3] = u_13;
    u_23 = pA[2 + bs * 3];
    u_23 -= pA[2 + bs * 0] * u_03;
    u_23 -= pA[2 + bs * 1] * u_13;
    pA[2 + bs * 3] = u_23;
    tmp3 = pA[3 + bs * 3];
    tmp3 -= pA[3 + bs * 0] * u_03;
    tmp3 -= pA[3 + bs * 1] * u_13;
    tmp3 -= pA[3 + bs * 2] * u_23;
    pA[3 + bs * 3] = tmp3;
    pB = pA + bs * sda;
    for (k = 0; k < ma - 3; k += 4) {
        tmp0 = pB[0 + bs * 3];
        tmp1 = pB[1 + bs * 3];
        tmp2 = pB[2 + bs * 3];
        tmp3 = pB[3 + bs * 3];
        tmp0 -= pB[0 + bs * 0] * u_03;
        tmp1 -= pB[1 + bs * 0] * u_03;
        tmp2 -= pB[2 + bs * 0] * u_03;
        tmp3 -= pB[3 + bs * 0] * u_03;
        tmp0 -= pB[0 + bs * 1] * u_13;
        tmp1 -= pB[1 + bs * 1] * u_13;
        tmp2 -= pB[2 + bs * 1] * u_13;
        tmp3 -= pB[3 + bs * 1] * u_13;
        tmp0 -= pB[0 + bs * 2] * u_23;
        tmp1 -= pB[1 + bs * 2] * u_23;
        tmp2 -= pB[2 + bs * 2] * u_23;
        tmp3 -= pB[3 + bs * 2] * u_23;
        pB[0 + bs * 3] = tmp0;
        pB[1 + bs * 3] = tmp1;
        pB[2 + bs * 3] = tmp2;
        pB[3 + bs * 3] = tmp3;
        pB += bs * sda;
    }
    for (; k < ma; k++) {
        tmp0 = pB[0 + bs * 3];
        tmp0 -= pB[0 + bs * 0] * u_03;
        tmp0 -= pB[0 + bs * 1] * u_13;
        tmp0 -= pB[0 + bs * 2] * u_23;
        pB[0 + bs * 3] = tmp0;
        pB += 1;
    }

    didamax_lib4(m - 3, 3, &pA[3 + bs * 3], sda, &idamax, &tmp3);
    ipiv[3] = idamax + 3;
    if (tmp3 != 0) {
        if (ipiv[3] != 3)
            kernel_drowsw_lib4(4, pA + 3, pA + ipiv[3] / bs * bs * sda + ipiv[3] % bs);

        tmp3 = 1.0 / pA[3 + bs * 3];
        inv_diag_A[3] = tmp3;
        pB = pA + bs * sda;
        for (k = 0; k < ma - 3; k += 4) {
            pB[0 + bs * 3] *= tmp3;
            pB[1 + bs * 3] *= tmp3;
            pB[2 + bs * 3] *= tmp3;
            pB[3 + bs * 3] *= tmp3;
            pB += bs * sda;
        }
        for (; k < ma; k++) {
            pB[0 + bs * 3] *= tmp3;
            pB += 1;
        }
    } else {
        inv_diag_A[3] = 0.0;
    }
}


// it process m>0 rows and 0<n<=4 cols
void kernel_dgetrf_pivot_4_vs_lib4(int m, double* pA, int sda, double* inv_diag_A, int* ipiv, int n) {

    if (m <= 0 || n <= 0) { return; }


    const int bs = 4;

    // assume m>=4
    int ma = m - 4;

    double
            tmp0,
            tmp1, tmp2, tmp3,
            u_00, u_01, u_02, u_03,
            u_11, u_12, u_13,
            u_22, u_23,
            u_33;

    double* pB;

    int
            k,
            idamax;

    int n4 = n < 4 ? n : 4;

    // first column

    // find pivot & scale
    didamax_lib4(m - 0, 0, &pA[0 + bs * 0], sda, &idamax, &tmp0);
    ipiv[0] = idamax;
    if (tmp0 != 0.0) {
        if (ipiv[0] != 0)
            kernel_drowsw_lib4(n4, pA + 0, pA + ipiv[0] / bs * bs * sda + ipiv[0] % bs);

        tmp0 = 1.0 / pA[0 + bs * 0];
        inv_diag_A[0] = tmp0;
        if (m >= 4) {
            pA[1 + bs * 0] *= tmp0;
            pA[2 + bs * 0] *= tmp0;
            pA[3 + bs * 0] *= tmp0;
            pB = pA + bs * sda;
            for (k = 0; k < ma - 3; k += 4) {
                pB[0 + bs * 0] *= tmp0;
                pB[1 + bs * 0] *= tmp0;
                pB[2 + bs * 0] *= tmp0;
                pB[3 + bs * 0] *= tmp0;
                pB += bs * sda;
            }
            for (; k < ma; k++) {
                pB[0 + bs * 0] *= tmp0;
                pB += 1;
            }
        } else  // m = {1,2,3}
        {
            if (m > 1) {
                pA[1 + bs * 0] *= tmp0;
                if (m > 2)
                    pA[2 + bs * 0] *= tmp0;
            }
        }
    } else {
        inv_diag_A[0] = 0.0;
    }

    if (n == 1 || m == 1) { return; }  // XXX for the first row there is nothing to do, so we can return here


    // second column

    // correct
    if (m >= 4) {
        u_01 = pA[0 + bs * 1];
        tmp1 = pA[1 + bs * 1];
        tmp2 = pA[2 + bs * 1];
        tmp3 = pA[3 + bs * 1];
        tmp1 -= pA[1 + bs * 0] * u_01;
        tmp2 -= pA[2 + bs * 0] * u_01;
        tmp3 -= pA[3 + bs * 0] * u_01;
        pA[1 + bs * 1] = tmp1;
        pA[2 + bs * 1] = tmp2;
        pA[3 + bs * 1] = tmp3;
        pB = pA + bs * sda;
        for (k = 0; k < ma - 3; k += 4) {
            tmp0 = pB[0 + bs * 1];
            tmp1 = pB[1 + bs * 1];
            tmp2 = pB[2 + bs * 1];
            tmp3 = pB[3 + bs * 1];
            tmp0 -= pB[0 + bs * 0] * u_01;
            tmp1 -= pB[1 + bs * 0] * u_01;
            tmp2 -= pB[2 + bs * 0] * u_01;
            tmp3 -= pB[3 + bs * 0] * u_01;
            pB[0 + bs * 1] = tmp0;
            pB[1 + bs * 1] = tmp1;
            pB[2 + bs * 1] = tmp2;
            pB[3 + bs * 1] = tmp3;
            pB += bs * sda;
        }
        for (; k < ma; k++) {
            tmp0 = pB[0 + bs * 1];
            tmp0 -= pB[0 + bs * 0] * u_01;
            pB[0 + bs * 1] = tmp0;
            pB += 1;
        }
    } else {  // m = {2,3}
        u_01 = pA[0 + bs * 1];
        tmp1 = pA[1 + bs * 1];
        tmp1 -= pA[1 + bs * 0] * u_01;
        pA[1 + bs * 1] = tmp1;
        if (m > 2) {
            tmp2 = pA[2 + bs * 1];
            tmp2 -= pA[2 + bs * 0] * u_01;
            pA[2 + bs * 1] = tmp2;
        }
    }

    // find pivot & scale
    didamax_lib4(m - 1, 1, &pA[1 + bs * 1], sda, &idamax, &tmp1);
    ipiv[1] = idamax + 1;
    if (tmp1 != 0) {
        if (ipiv[1] != 1)
            kernel_drowsw_lib4(n4, pA + 1, pA + ipiv[1] / bs * bs * sda + ipiv[1] % bs);

        tmp1 = 1.0 / pA[1 + bs * 1];
        inv_diag_A[1] = tmp1;
        if (m >= 4) {
            pA[2 + bs * 1] *= tmp1;
            pA[3 + bs * 1] *= tmp1;
            pB = pA + bs * sda;
            for (k = 0; k < ma - 3; k += 4) {
                pB[0 + bs * 1] *= tmp1;
                pB[1 + bs * 1] *= tmp1;
                pB[2 + bs * 1] *= tmp1;
                pB[3 + bs * 1] *= tmp1;
                pB += bs * sda;
            }
            for (; k < ma; k++) {
                pB[0 + bs * 1] *= tmp1;
                pB += 1;
            }
        } else  // m = {2,3}
        {
            if (m > 2)
                pA[2 + bs * 1] *= tmp1;
        }
    } else {
        inv_diag_A[1] = 0.0;
    }

    if (n == 2)


        // third column

        // correct
        if (m >= 4) {
            u_02 = pA[0 + bs * 2];
            u_12 = pA[1 + bs * 2];
            u_12 -= pA[1 + bs * 0] * u_02;
            pA[1 + bs * 2] = u_12;
            tmp2 = pA[2 + bs * 2];
            tmp3 = pA[3 + bs * 2];
            tmp2 -= pA[2 + bs * 0] * u_02;
            tmp3 -= pA[3 + bs * 0] * u_02;
            tmp2 -= pA[2 + bs * 1] * u_12;
            tmp3 -= pA[3 + bs * 1] * u_12;
            pA[2 + bs * 2] = tmp2;
            pA[3 + bs * 2] = tmp3;
            pB = pA + bs * sda;
            for (k = 0; k < ma - 3; k += 4) {
                tmp0 = pB[0 + bs * 2];
                tmp1 = pB[1 + bs * 2];
                tmp2 = pB[2 + bs * 2];
                tmp3 = pB[3 + bs * 2];
                tmp0 -= pB[0 + bs * 0] * u_02;
                tmp1 -= pB[1 + bs * 0] * u_02;
                tmp2 -= pB[2 + bs * 0] * u_02;
                tmp3 -= pB[3 + bs * 0] * u_02;
                tmp0 -= pB[0 + bs * 1] * u_12;
                tmp1 -= pB[1 + bs * 1] * u_12;
                tmp2 -= pB[2 + bs * 1] * u_12;
                tmp3 -= pB[3 + bs * 1] * u_12;
                pB[0 + bs * 2] = tmp0;
                pB[1 + bs * 2] = tmp1;
                pB[2 + bs * 2] = tmp2;
                pB[3 + bs * 2] = tmp3;
                pB += bs * sda;
            }
            for (; k < ma; k++) {
                tmp0 = pB[0 + bs * 2];
                tmp0 -= pB[0 + bs * 0] * u_02;
                tmp0 -= pB[0 + bs * 1] * u_12;
                pB[0 + bs * 2] = tmp0;
                pB += 1;
            }
        } else {  // m = {2,3}
            u_02 = pA[0 + bs * 2];
            u_12 = pA[1 + bs * 2];
            u_12 -= pA[1 + bs * 0] * u_02;
            pA[1 + bs * 2] = u_12;
            if (m > 2) {
                tmp2 = pA[2 + bs * 2];
                tmp2 -= pA[2 + bs * 0] * u_02;
                tmp2 -= pA[2 + bs * 1] * u_12;
                pA[2 + bs * 2] = tmp2;
            }
        }

    // find pivot & scale
    if (m > 2) {
        didamax_lib4(m - 2, 2, &pA[2 + bs * 2], sda, &idamax, &tmp2);
        ipiv[2] = idamax + 2;
        if (tmp2 != 0) {
            if (ipiv[2] != 2)
                kernel_drowsw_lib4(n4, pA + 2, pA + ipiv[2] / bs * bs * sda + ipiv[2] % bs);

            tmp2 = 1.0 / pA[2 + bs * 2];
            inv_diag_A[2] = tmp2;
            if (m >= 4) {
                pA[3 + bs * 2] *= tmp2;
                pB = pA + bs * sda;
                for (k = 0; k < ma - 3; k += 4) {
                    pB[0 + bs * 2] *= tmp2;
                    pB[1 + bs * 2] *= tmp2;
                    pB[2 + bs * 2] *= tmp2;
                    pB[3 + bs * 2] *= tmp2;
                    pB += bs * sda;
                }
                for (; k < ma; k++) {
                    pB[0 + bs * 2] *= tmp2;
                    pB += 1;
                }
            }
        } else {
            inv_diag_A[2] = 0.0;
        }
    }

    if (n < 4)


        // fourth column

        // correct
        if (m >= 4) {
            u_03 = pA[0 + bs * 3];
            u_13 = pA[1 + bs * 3];
            u_13 -= pA[1 + bs * 0] * u_03;
            pA[1 + bs * 3] = u_13;
            u_23 = pA[2 + bs * 3];
            u_23 -= pA[2 + bs * 0] * u_03;
            u_23 -= pA[2 + bs * 1] * u_13;
            pA[2 + bs * 3] = u_23;
            tmp3 = pA[3 + bs * 3];
            tmp3 -= pA[3 + bs * 0] * u_03;
            tmp3 -= pA[3 + bs * 1] * u_13;
            tmp3 -= pA[3 + bs * 2] * u_23;
            pA[3 + bs * 3] = tmp3;
            pB = pA + bs * sda;
            for (k = 0; k < ma - 3; k += 4) {
                tmp0 = pB[0 + bs * 3];
                tmp1 = pB[1 + bs * 3];
                tmp2 = pB[2 + bs * 3];
                tmp3 = pB[3 + bs * 3];
                tmp0 -= pB[0 + bs * 0] * u_03;
                tmp1 -= pB[1 + bs * 0] * u_03;
                tmp2 -= pB[2 + bs * 0] * u_03;
                tmp3 -= pB[3 + bs * 0] * u_03;
                tmp0 -= pB[0 + bs * 1] * u_13;
                tmp1 -= pB[1 + bs * 1] * u_13;
                tmp2 -= pB[2 + bs * 1] * u_13;
                tmp3 -= pB[3 + bs * 1] * u_13;
                tmp0 -= pB[0 + bs * 2] * u_23;
                tmp1 -= pB[1 + bs * 2] * u_23;
                tmp2 -= pB[2 + bs * 2] * u_23;
                tmp3 -= pB[3 + bs * 2] * u_23;
                pB[0 + bs * 3] = tmp0;
                pB[1 + bs * 3] = tmp1;
                pB[2 + bs * 3] = tmp2;
                pB[3 + bs * 3] = tmp3;
                pB += bs * sda;
            }
            for (; k < ma; k++) {
                tmp0 = pB[0 + bs * 3];
                tmp0 -= pB[0 + bs * 0] * u_03;
                tmp0 -= pB[0 + bs * 1] * u_13;
                tmp0 -= pB[0 + bs * 2] * u_23;
                pB[0 + bs * 3] = tmp0;
                pB += 1;
            }
        } else {  // m = {2,3}
            u_03 = pA[0 + bs * 3];
            u_13 = pA[1 + bs * 3];
            u_13 -= pA[1 + bs * 0] * u_03;
            pA[1 + bs * 3] = u_13;
            if (m > 2) {
                u_23 = pA[2 + bs * 3];
                u_23 -= pA[2 + bs * 0] * u_03;
                u_23 -= pA[2 + bs * 1] * u_13;
                pA[2 + bs * 3] = u_23;
            }
        }

    if (m > 3) {
        // find pivot & scale
        didamax_lib4(m - 3, 3, &pA[3 + bs * 3], sda, &idamax, &tmp3);
        ipiv[3] = idamax + 3;
        if (tmp3 != 0) {
            if (ipiv[3] != 3)
                kernel_drowsw_lib4(n4, pA + 3, pA + ipiv[3] / bs * bs * sda + ipiv[3] % bs);

            tmp3 = 1.0 / pA[3 + bs * 3];
            inv_diag_A[3] = tmp3;
            pB = pA + bs * sda;
            for (k = 0; k < ma - 3; k += 4) {
                pB[0 + bs * 3] *= tmp3;
                pB[1 + bs * 3] *= tmp3;
                pB[2 + bs * 3] *= tmp3;
                pB[3 + bs * 3] *= tmp3;
                pB += bs * sda;
            }
            for (; k < ma; k++) {
                pB[0 + bs * 3] *= tmp3;
                pB += 1;
            }
        } else {
            inv_diag_A[3] = 0.0;
        }
    }
}
