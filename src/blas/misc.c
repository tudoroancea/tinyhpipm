#include "tinyhpipm/blas/misc.h"
#include "tinyhpipm/blas/kernel.h"
#include "tinyhpipm/blas/struct.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// add scaled vector to diagonal
void ddiaad_lib(int kmax, double alpha, double* x, int offset, double* pD, int sdd) {
    const int bs = D_PS;

    int kna = (bs - offset % bs) % bs;
    kna = kmax < kna ? kmax : kna;

    int jj, ll;

    if (kna > 0) {
        for (ll = 0; ll < kna; ll++) {
            pD[ll + bs * ll] += alpha * x[ll];
        }
        pD += kna + bs * (sdd - 1) + kna * bs;
        x += kna;
        kmax -= kna;
    }
    for (jj = 0; jj < kmax - 3; jj += 4) {
        pD[jj * sdd + (jj + 0) * bs + 0] += alpha * x[jj + 0];
        pD[jj * sdd + (jj + 1) * bs + 1] += alpha * x[jj + 1];
        pD[jj * sdd + (jj + 2) * bs + 2] += alpha * x[jj + 2];
        pD[jj * sdd + (jj + 3) * bs + 3] += alpha * x[jj + 3];
    }
    for (ll = 0; ll < kmax - jj; ll++) {
        pD[jj * sdd + (jj + ll) * bs + ll] += alpha * x[jj + ll];
    }
}


// insert vector to diagonal, sparse formulation
void ddiain_libsp(int kmax, int* idx, double alpha, double* x, double* pD, int sdd) {

    const int bs = D_PS;

    int ii, jj;

    for (jj = 0; jj < kmax; jj++) {
        ii = idx[jj];
        pD[ii / bs * bs * sdd + ii % bs + ii * bs] = alpha * x[jj];
    }
}


// add scaled vector to diagonal, sparse formulation
void ddiaad_libsp(int kmax, int* idx, double alpha, double* x, double* pD, int sdd) {

    const int bs = D_PS;

    int ii, jj;

    for (jj = 0; jj < kmax; jj++) {
        ii = idx[jj];
        pD[ii / bs * bs * sdd + ii % bs + ii * bs] += alpha * x[jj];
    }
}


// add scaled vector to another vector and insert to diagonal, sparse formulation
void ddiaadin_libsp(int kmax, int* idx, double alpha, double* x, double* y, double* pD, int sdd) {

    const int bs = D_PS;

    int ii, jj;

    for (jj = 0; jj < kmax; jj++) {
        ii = idx[jj];
        pD[ii / bs * bs * sdd + ii % bs + ii * bs] = y[jj] + alpha * x[jj];
    }
}


// add scaled vector to another vector and insert to row, sparse formulation
void drowadin_libsp(int kmax, int* idx, double alpha, double* x, double* y, double* pD) {

    const int bs = D_PS;

    int ii, jj;

    for (jj = 0; jj < kmax; jj++) {
        ii = idx[jj];
        pD[ii * bs] = y[jj] + alpha * x[jj];
    }
}


// insert vector to diagonal, sparse formulation
void dcolin_libsp(int kmax, int* idx, double* x, double* pD, int sdd) {

    const int bs = D_PS;

    int ii, jj;

    for (jj = 0; jj < kmax; jj++) {
        ii = idx[jj];
        pD[ii / bs * bs * sdd + ii % bs] = x[jj];
    }
}


// add scaled vector to diagonal, sparse formulation
void dcolad_libsp(int kmax, double alpha, int* idx, double* x, double* pD, int sdd) {

    const int bs = D_PS;

    int ii, jj;

    for (jj = 0; jj < kmax; jj++) {
        ii = idx[jj];
        pD[ii / bs * bs * sdd + ii % bs] += alpha * x[jj];
    }
}


// swaps two cols
void dcolsw_lib(int kmax, int offsetA, double* pA, int sda, int offsetC, double* pC, int sdc) {

    const int bs = D_PS;

    int ii;

    double tmp;

    if (offsetA == offsetC) {
        if (offsetA > 0) {
            ii = 0;
            for (; ii < bs - offsetA; ii++) {
                tmp = pA[0 + bs * 0];
                pA[0 + bs * 0] = pC[0 + bs * 0];
                pC[0 + bs * 0] = tmp;
                pA += 1;
                pC += 1;
            }
            pA += bs * (sda - 1);
            pC += bs * (sdc - 1);
            kmax -= bs - offsetA;
        }
        ii = 0;
        for (; ii < kmax - 3; ii += 4) {
            tmp = pA[0 + bs * 0];
            pA[0 + bs * 0] = pC[0 + bs * 0];
            pC[0 + bs * 0] = tmp;
            tmp = pA[1 + bs * 0];
            pA[1 + bs * 0] = pC[1 + bs * 0];
            pC[1 + bs * 0] = tmp;
            tmp = pA[2 + bs * 0];
            pA[2 + bs * 0] = pC[2 + bs * 0];
            pC[2 + bs * 0] = tmp;
            tmp = pA[3 + bs * 0];
            pA[3 + bs * 0] = pC[3 + bs * 0];
            pC[3 + bs * 0] = tmp;
            pA += bs * sda;
            pC += bs * sdc;
        }
        for (; ii < kmax; ii++) {
            tmp = pA[0 + bs * 0];
            pA[0 + bs * 0] = pC[0 + bs * 0];
            pC[0 + bs * 0] = tmp;
            pA += 1;
            pC += 1;
        }
    } else {
        fprintf(stderr, "\n\tdcolsw: feature not implemented yet: offsetA!=offsetC\n\n");
        exit(1);
    }
}


// insert vector to vector, sparse formulation
void vecin_libsp(int kmax, int* idx, double* x, double* y) {

    int jj;

    for (jj = 0; jj < kmax; jj++) {
        y[idx[jj]] = x[jj];
    }
}


// adds vector to vector, sparse formulation
void vecad_libsp(int kmax, int* idx, double alpha, double* x, double* y) {

    int jj;

    for (jj = 0; jj < kmax; jj++) {
        y[idx[jj]] += alpha * x[jj];
    }
}


// insert element into strmat
void dgein1(double a, struct mat* sA, int ai, int aj) {
    if (ai == aj) {
        // invalidate stored inverse diagonal
        sA->use_dA = 0;
    }

    const int bs = D_PS;
    int sda = sA->cn;
    double* pA = sA->pA + ai / bs * bs * sda + ai % bs + aj * bs;
    pA[0] = a;
}


// extract element from strmat
double dgeex1(struct mat* sA, int ai, int aj) {
    const int bs = D_PS;
    int sda = sA->cn;
    double* pA = sA->pA + ai / bs * bs * sda + ai % bs + aj * bs;
    return pA[0];
}


// insert element into strvec
void vecin1(double a, struct vec* sx, int xi) {
    double* x = sx->pa + xi;
    x[0] = a;
}


// extract element from strvec
double vecex1(struct vec* sx, int xi) {
    double* x = sx->pa + xi;
    return x[0];
}


// set all elements of a strmat to a value
void dgese(int m, int n, double alpha, struct mat* sA, int ai, int aj) {

    // invalidate stored inverse diagonal
    sA->use_dA = 0;

    const int bs = D_PS;
    int sda = sA->cn;
    double* pA = sA->pA + ai % bs + ai / bs * bs * sda + aj * bs;
    int m0 = m < (bs - ai % bs) % bs ? m : (bs - ai % bs) % bs;
    int ii, jj;
    if (m0 > 0) {
        for (ii = 0; ii < m0; ii++) {
            for (jj = 0; jj < n; jj++) {
                pA[jj * bs] = alpha;
            }
            pA += 1;
        }
        pA += bs * (sda - 1);
        m -= m0;
    }
    for (ii = 0; ii < m - 3; ii += 4) {
        for (jj = 0; jj < n; jj++) {
            pA[0 + jj * bs] = alpha;
            pA[1 + jj * bs] = alpha;
            pA[2 + jj * bs] = alpha;
            pA[3 + jj * bs] = alpha;
        }
        pA += bs * sda;
    }
    for (; ii < m; ii++) {
        for (jj = 0; jj < n; jj++) {
            pA[jj * bs] = alpha;
        }
        pA += 1;
    }
}


// set all elements of a strvec to a value
void dvecse(int m, double alpha, struct vec* sx, int xi) {
    double* x = sx->pa + xi;
    int ii;
    for (ii = 0; ii < m; ii++)
        x[ii] = alpha;
}


// insert a vector into diagonal
void ddiain(int kmax, double alpha, struct vec* sx, int xi, struct mat* sA, int ai, int aj) {

    // invalidate stored inverse diagonal
    sA->use_dA = 0;

    const int bs = D_PS;
    int sda = sA->cn;
    double* pA = sA->pA + ai / bs * bs * sda + ai % bs + aj * bs;
    double* x = sx->pa + xi;
    int offsetA = ai % bs;

    int kna = (bs - offsetA % bs) % bs;
    kna = kmax < kna ? kmax : kna;

    int jj, ll;

    if (kna > 0) {
        for (ll = 0; ll < kna; ll++) {
            pA[ll + bs * ll] = alpha * x[ll];
        }
        pA += kna + bs * (sda - 1) + kna * bs;
        x += kna;
        kmax -= kna;
    }
    for (jj = 0; jj < kmax - 3; jj += 4) {
        pA[jj * sda + (jj + 0) * bs + 0] = alpha * x[jj + 0];
        pA[jj * sda + (jj + 1) * bs + 1] = alpha * x[jj + 1];
        pA[jj * sda + (jj + 2) * bs + 2] = alpha * x[jj + 2];
        pA[jj * sda + (jj + 3) * bs + 3] = alpha * x[jj + 3];
    }
    for (ll = 0; ll < kmax - jj; ll++) {
        pA[jj * sda + (jj + ll) * bs + ll] = alpha * x[jj + ll];
    }
}


// add scalar to diagonal
void ddiare(int kmax, double alpha, struct mat* sA, int ai, int aj) {

    // invalidate stored inverse diagonal
    sA->use_dA = 0;

    const int bs = D_PS;
    int sda = sA->cn;
    double* pA = sA->pA + ai / bs * bs * sda + ai % bs + aj * bs;
    int offsetA = ai % bs;

    int kna = (bs - offsetA % bs) % bs;
    kna = kmax < kna ? kmax : kna;

    int jj, ll;

    if (kna > 0) {
        for (ll = 0; ll < kna; ll++) {
            pA[ll + bs * ll] += alpha;
        }
        pA += kna + bs * (sda - 1) + kna * bs;
        kmax -= kna;
    }
    for (jj = 0; jj < kmax - 3; jj += 4) {
        pA[jj * sda + (jj + 0) * bs + 0] += alpha;
        pA[jj * sda + (jj + 1) * bs + 1] += alpha;
        pA[jj * sda + (jj + 2) * bs + 2] += alpha;
        pA[jj * sda + (jj + 3) * bs + 3] += alpha;
    }
    for (ll = 0; ll < kmax - jj; ll++) {
        pA[jj * sda + (jj + ll) * bs + ll] += alpha;
    }
}


// swap two rows of two matrix structs
void drowsw(int kmax, struct mat* sA, int ai, int aj, struct mat* sC, int ci, int cj) {

    // invalidate stored inverse diagonal
    sA->use_dA = 0;
    sC->use_dA = 0;

    const int bs = D_PS;
    int sda = sA->cn;
    double* pA = sA->pA + ai / bs * bs * sda + ai % bs + aj * bs;
    int sdc = sC->cn;
    double* pC = sC->pA + ci / bs * bs * sdc + ci % bs + cj * bs;
    kernel_drowsw_lib4(kmax, pA, pC);
}


// permute the rows of a matrix struct
void drowpe(int kmax, int* ipiv, struct mat* sA) {

    // invalidate stored inverse diagonal
    sA->use_dA = 0;

    int ii;
    for (ii = 0; ii < kmax; ii++) {
        if (ipiv[ii] != ii)
            drowsw(sA->n, sA, ii, 0, sA, ipiv[ii], 0);
    }
}


// inverse permute the rows of a matrix struct
void drowpei(int kmax, int* ipiv, struct mat* sA) {

    // invalidate stored inverse diagonal
    sA->use_dA = 0;

    int ii;
    for (ii = kmax - 1; ii >= 0; ii--) {
        if (ipiv[ii] != ii)
            drowsw(sA->n, sA, ii, 0, sA, ipiv[ii], 0);
    }
}


// extract row to vector
void drowex_lib(int kmax, double alpha, double* pD, double* x) {

    const int bs = D_PS;

    int jj, ll;

    for (jj = 0; jj < kmax - 3; jj += 4) {
        x[jj + 0] = alpha * pD[(jj + 0) * bs];
        x[jj + 1] = alpha * pD[(jj + 1) * bs];
        x[jj + 2] = alpha * pD[(jj + 2) * bs];
        x[jj + 3] = alpha * pD[(jj + 3) * bs];
    }
    for (; jj < kmax; jj++) {
        x[jj] = alpha * pD[(jj) *bs];
    }
}
// extract a row int a vector
void drowex(int kmax, double alpha, struct mat* sA, int ai, int aj, struct vec* sx, int xi) {
    const int bs = D_PS;
    int sda = sA->cn;
    double* pA = sA->pA + ai / bs * bs * sda + ai % bs + aj * bs;
    double* x = sx->pa + xi;
    drowex_lib(kmax, alpha, pA, x);
}

// insert vector to row
void drowin_lib(int kmax, double alpha, double* x, double* pD) {

    const int bs = D_PS;

    int jj, ll;

    for (jj = 0; jj < kmax - 3; jj += 4) {
        pD[(jj + 0) * bs] = alpha * x[jj + 0];
        pD[(jj + 1) * bs] = alpha * x[jj + 1];
        pD[(jj + 2) * bs] = alpha * x[jj + 2];
        pD[(jj + 3) * bs] = alpha * x[jj + 3];
    }
    for (; jj < kmax; jj++) {
        pD[(jj) *bs] = alpha * x[jj];
    }
}


// insert a vector into a row
void drowin(int kmax, double alpha, struct vec* sx, int xi, struct mat* sA, int ai, int aj) {

    // invalidate stored inverse diagonal
    sA->use_dA = 0;

    const int bs = D_PS;
    int sda = sA->cn;
    double* pA = sA->pA + ai / bs * bs * sda + ai % bs + aj * bs;
    double* x = sx->pa + xi;
    drowin_lib(kmax, alpha, x, pA);
}

// add scaled vector to row
void drowad_lib(int kmax, double alpha, double* x, double* pD) {

    const int bs = D_PS;

    int jj, ll;

    for (jj = 0; jj < kmax - 3; jj += 4) {
        pD[(jj + 0) * bs] += alpha * x[jj + 0];
        pD[(jj + 1) * bs] += alpha * x[jj + 1];
        pD[(jj + 2) * bs] += alpha * x[jj + 2];
        pD[(jj + 3) * bs] += alpha * x[jj + 3];
    }
    for (; jj < kmax; jj++) {
        pD[(jj) *bs] += alpha * x[jj];
    }
}

// add a vector to a row
void drowad(int kmax, double alpha, struct vec* sx, int xi, struct mat* sA, int ai, int aj) {

    // invalidate stored inverse diagonal
    sA->use_dA = 0;

    const int bs = D_PS;
    int sda = sA->cn;
    double* pA = sA->pA + ai / bs * bs * sda + ai % bs + aj * bs;
    double* x = sx->pa + xi;
    drowad_lib(kmax, alpha, x, pA);
}

// extract vector from column
void dcolex_lib(int kmax, int offset, double* pD, int sdd, double* x) {

    const int bs = D_PS;

    int kna = (bs - offset % bs) % bs;
    kna = kmax < kna ? kmax : kna;

    int jj, ll;

    if (kna > 0) {
        for (ll = 0; ll < kna; ll++) {
            x[ll] = pD[ll];
        }

        pD += kna + bs * (sdd - 1);
        x += kna;
        kmax -= kna;
    }

    for (jj = 0; jj < kmax - 3; jj += 4) {
        x[jj + 0] = pD[jj * sdd + 0];
        x[jj + 1] = pD[jj * sdd + 1];
        x[jj + 2] = pD[jj * sdd + 2];
        x[jj + 3] = pD[jj * sdd + 3];
    }
    for (ll = 0; ll < kmax - jj; ll++) {
        x[jj + ll] = pD[jj * sdd + ll];
    }
}


// extract vector from column
void dcolex(int kmax, struct mat* sA, int ai, int aj, struct vec* sx, int xi) {

    // invalidate stored inverse diagonal
    sA->use_dA = 0;

    const int bs = D_PS;
    int sda = sA->cn;
    double* pA = sA->pA + ai / bs * bs * sda + ai % bs + aj * bs;
    double* x = sx->pa + xi;
    dcolex_lib(kmax, ai % bs, pA, sda, x);
}


// insert vector to column
void dcolin_lib(int kmax, double* x, int offset, double* pD, int sdd) {

    const int bs = D_PS;

    int kna = (bs - offset % bs) % bs;
    kna = kmax < kna ? kmax : kna;

    int jj, ll;

    if (kna > 0) {
        for (ll = 0; ll < kna; ll++) {
            pD[ll] = x[ll];
        }
        pD += kna + bs * (sdd - 1);
        x += kna;
        kmax -= kna;
    }
    for (jj = 0; jj < kmax - 3; jj += 4) {
        pD[jj * sdd + 0] = x[jj + 0];
        pD[jj * sdd + 1] = x[jj + 1];
        pD[jj * sdd + 2] = x[jj + 2];
        pD[jj * sdd + 3] = x[jj + 3];
    }
    for (ll = 0; ll < kmax - jj; ll++) {
        pD[jj * sdd + ll] = x[jj + ll];
    }
}


// insert as vector as a column
void dcolin(int kmax, struct vec* sx, int xi, struct mat* sA, int ai, int aj) {
    // invalidate stored inverse diagonal
    sA->use_dA = 0;
    const int bs = D_PS;
    int sda = sA->cn;
    double* pA = sA->pA + ai / bs * bs * sda + ai % bs + aj * bs;
    double* x = sx->pa + xi;
    dcolin_lib(kmax, x, ai % bs, pA, sda);
}

// add scaled vector to column
void dcolad_lib(int kmax, double alpha, double* x, int offset, double* pD, int sdd) {

    const int bs = D_PS;

    int kna = (bs - offset % bs) % bs;
    kna = kmax < kna ? kmax : kna;

    int jj, ll;

    if (kna > 0) {
        for (ll = 0; ll < kna; ll++) {
            pD[ll] += alpha * x[ll];
        }
        pD += kna + bs * (sdd - 1);
        x += kna;
        kmax -= kna;
    }
    for (jj = 0; jj < kmax - 3; jj += 4) {
        pD[jj * sdd + 0] += alpha * x[jj + 0];
        pD[jj * sdd + 1] += alpha * x[jj + 1];
        pD[jj * sdd + 2] += alpha * x[jj + 2];
        pD[jj * sdd + 3] += alpha * x[jj + 3];
    }
    for (ll = 0; ll < kmax - jj; ll++) {
        pD[jj * sdd + ll] += alpha * x[jj + ll];
    }
}


// add scaled vector to column
void dcolad(int kmax, double alpha, struct vec* sx, int xi, struct mat* sA, int ai, int aj) {
    // invalidate stored inverse diagonal
    sA->use_dA = 0;
    const int bs = D_PS;
    int sda = sA->cn;
    double* pA = sA->pA + ai / bs * bs * sda + ai % bs + aj * bs;
    double* x = sx->pa + xi;
    dcolad_lib(kmax, alpha, x, ai % bs, pA, sda);
}


// scale a column
void dcolsc(int kmax, double alpha, struct mat* sA, int ai, int aj) {
    // invalidate stored inverse diagonal
    sA->use_dA = 0;
    const int bs = D_PS;
    int sda = sA->cn;
    double* pA = sA->pA + ai / bs * bs * sda + ai % bs + aj * bs;
    int kna = (bs - ai % bs) % bs;
    kna = kmax < kna ? kmax : kna;
    int jj, ll;

    if (kna > 0) {
        for (ll = 0; ll < kna; ll++) {
            pA[ll] *= alpha;
        }
        pA += kna + bs * (sda - 1);
        kmax -= kna;
    }
    for (jj = 0; jj < kmax - 3; jj += 4) {
        pA[jj * sda + 0] *= alpha;
        pA[jj * sda + 1] *= alpha;
        pA[jj * sda + 2] *= alpha;
        pA[jj * sda + 3] *= alpha;
    }
    for (ll = 0; ll < kmax - jj; ll++) {
        pA[jj * sda + ll] *= alpha;
    }
}


// swap two cols of two matrix structs
void dcolsw(int kmax, struct mat* sA, int ai, int aj, struct mat* sC, int ci, int cj) {

    // invalidate stored inverse diagonal
    sA->use_dA = 0;
    sC->use_dA = 0;

    const int bs = D_PS;
    int sda = sA->cn;
    double* pA = sA->pA + ai / bs * bs * sda + ai % bs + aj * bs;
    int sdc = sC->cn;
    double* pC = sC->pA + ci / bs * bs * sdc + ci % bs + cj * bs;
    dcolsw_lib(kmax, ai % bs, pA, sda, ci % bs, pC, sdc);
}


// permute the cols of a matrix struct
void dcolpe(int kmax, int* ipiv, struct mat* sA) {

    // invalidate stored inverse diagonal
    sA->use_dA = 0;

    int ii;
    for (ii = 0; ii < kmax; ii++) {
        if (ipiv[ii] != ii)
            dcolsw(sA->m, sA, 0, ii, sA, 0, ipiv[ii]);
    }
}


// inverse permute the cols of a matrix struct
void dcolpei(int kmax, int* ipiv, struct mat* sA) {

    // invalidate stored inverse diagonal
    sA->use_dA = 0;

    int ii;
    for (ii = kmax - 1; ii >= 0; ii--) {
        if (ipiv[ii] != ii)
            dcolsw(sA->m, sA, 0, ii, sA, 0, ipiv[ii]);
    }
}


// copy a generic strmat into a generic strmat
void dgecp(int m, int n, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj) {
    // invalidate stored inverse diagonal
    sB->use_dA = 0;

    const int bs = D_PS;
    const double alpha = 1.0;

    // get submatrices
    int sda = sA->cn;
    double* A = sA->pA + ai / bs * bs * sda + ai % bs + aj * bs;
    int sdb = sB->cn;
    double* B = sB->pA + bi / bs * bs * sdb + bi % bs + bj * bs;

    if (m <= 0 || n <= 0) {
        return;
    }

    int mna, ii;

    // compute offset from closest panels
    int offA = ai % bs;
    int offB = bi % bs;

    // A at the beginning of the block
    A -= offA;
    // A at the beginning of the block
    B -= offB;

    // same alignment
    if (offA == offB) {
        ii = 0;
        // clean up at the beginning
        mna = (4 - offB) % bs;
        if (mna > 0) {
            if (m < mna)  // mna<=3  ==>  m = { 1, 2 }
            {
                if (m == 1) {
                    kernel_dgecpsc_1_0_lib4(0, n, alpha, A + offA, B + offB);
                    return;
                } else  // if(m==2 && mna==3)
                {
                    kernel_dgecpsc_2_0_lib4(0, n, alpha, A + offA, B + offB);
                    return;
                }
            }
            if (mna == 1) {
                kernel_dgecpsc_1_0_lib4(0, n, alpha, A + offA, B + offB);
                A += 4 * sda;
                B += 4 * sdb;
                ii += 1;
            } else if (mna == 2) {
                kernel_dgecpsc_2_0_lib4(0, n, alpha, A + offA, B + offB);
                A += 4 * sda;
                B += 4 * sdb;
                ii += 2;
            } else  // if(mna==3)
            {
                kernel_dgecpsc_3_0_lib4(0, n, alpha, A + offA, B + offB);
                A += 4 * sda;
                B += 4 * sdb;
                ii += 3;
            }
        }
        // main loop
        for (; ii < m - 3; ii += 4) {
            kernel_dgecp_4_0_lib4(0, n, A, B);
            A += 4 * sda;
            B += 4 * sdb;
        }
        // clean up at the end
        if (ii < m) {
            if (m - ii == 1)
                kernel_dgecpsc_1_0_lib4(0, n, alpha, A, B);
            else if (m - ii == 2)
                kernel_dgecpsc_2_0_lib4(0, n, alpha, A, B);
            else  // if(m-ii==3)
                kernel_dgecpsc_3_0_lib4(0, n, alpha, A, B);
        }
    }
    // skip one element of A
    else if (offA == (offB + 1) % bs) {
        ii = 0;
        // clean up at the beginning
        mna = (4 - offB) % bs;
        if (mna > 0) {
            if (m < mna)  // mna<=3  ==>  m = { 1, 2 }
            {
                if (m == 1) {
                    kernel_dgecpsc_1_0_lib4(0, n, alpha, A + offA, B + offB);
                    return;
                } else  // if(m==2 && mna==3)
                {
                    kernel_dgecpsc_2_0_lib4(0, n, alpha, A + offA, B + offB);
                    return;
                }
            }
            if (mna == 1) {
                kernel_dgecpsc_1_0_lib4(0, n, alpha, A + offA, B + offB);
                // A += 4*sda;
                B += 4 * sdb;
                ii += 1;
            } else if (mna == 2) {
                kernel_dgecpsc_2_3_lib4(0, n, alpha, A, sda, B + 2);
                A += 4 * sda;
                B += 4 * sdb;
                ii += 2;
            } else  // if(mna==3)
            {
                kernel_dgecpsc_3_2_lib4(0, n, alpha, A, sda, B + 1);
                A += 4 * sda;
                B += 4 * sdb;
                ii += 3;
            }
        }
        // main loop
        for (; ii < m - 3; ii += 4) {
            kernel_dgecpsc_4_1_lib4(0, n, alpha, A, sda, B);
            A += 4 * sda;
            B += 4 * sdb;
        }
        // clean up at the end
        if (ii < m) {
            if (m - ii == 1)
                kernel_dgecpsc_1_0_lib4(0, n, alpha, A + 1, B);
            else if (m - ii == 2)
                kernel_dgecpsc_2_0_lib4(0, n, alpha, A + 1, B);
            else  // if(m-ii==3)
                kernel_dgecpsc_3_0_lib4(0, n, alpha, A + 1, B);
        }
    }
    // skip 2 elements of A
    else if (offA == (offB + 2) % bs) {
        ii = 0;
        // clean up at the beginning
        mna = (4 - offB) % bs;
        if (mna > 0) {
            if (m < mna) {
                if (m == 1) {
                    kernel_dgecpsc_1_0_lib4(0, n, alpha, A + offA, B + offB);
                    return;
                } else  // if(m==2 && mna==3)
                {
                    kernel_dgecpsc_2_3_lib4(0, n, alpha, A, sda, B + 1);
                    return;
                }
            }
            if (mna == 1) {
                kernel_dgecpsc_1_0_lib4(0, n, alpha, A + 1, B + 3);
                // A += 4*sda;
                B += 4 * sdb;
                ii += 1;
            } else if (mna == 2) {
                kernel_dgecpsc_2_0_lib4(0, n, alpha, A, B + 2);
                // A += 4*sda;
                B += 4 * sdb;
                ii += 2;
            } else  // if(mna==3)
            {
                kernel_dgecpsc_3_3_lib4(0, n, alpha, A, sda, B + 1);
                A += 4 * sda;
                B += 4 * sdb;
                ii += 3;
            }
        }
        // main loop
        for (; ii < m - 3; ii += 4) {
            kernel_dgecpsc_4_2_lib4(0, n, alpha, A, sda, B);
            A += 4 * sda;
            B += 4 * sdb;
        }
        // clean up at the end
        if (ii < m) {
            if (m - ii == 1)
                kernel_dgecpsc_1_0_lib4(0, n, alpha, A + 2, B);
            else if (m - ii == 2)
                kernel_dgecpsc_2_0_lib4(0, n, alpha, A + 2, B);
            else  // if(m-ii==3)
                kernel_dgecpsc_3_2_lib4(0, n, alpha, A, sda, B);
        }
    }
    // skip 3 elements of A
    else  // if(offA==(offB+3)%bs)
    {
        ii = 0;
        // clean up at the beginning
        mna = (4 - offB) % bs;
        if (mna > 0) {
            if (m < mna) {
                if (m == 1) {
                    kernel_dgecpsc_1_0_lib4(0, n, alpha, A + offA, B + offB);
                    return;
                } else  // if(m==2 && mna==3)
                {
                    kernel_dgecpsc_2_0_lib4(0, n, alpha, A + offA, B + offB);
                    return;
                }
            }
            if (mna == 1) {
                kernel_dgecpsc_1_0_lib4(0, n, alpha, A + offA, B + offB);
                // A += 4*sda;
                B += 4 * sdb;
                ii += 1;
            } else if (mna == 2) {
                kernel_dgecpsc_2_0_lib4(0, n, alpha, A + offA, B + offB);
                // A += 4*sda;
                B += 4 * sdb;
                ii += 2;
            } else  // if(mna==3)
            {
                kernel_dgecpsc_3_0_lib4(0, n, alpha, A + offA, B + offB);
                // A += 4*sda;
                B += 4 * sdb;
                ii += 3;
            }
        }
        // main loop
        for (; ii < m - 3; ii += 4) {
            kernel_dgecpsc_4_3_lib4(0, n, alpha, A, sda, B);
            A += 4 * sda;
            B += 4 * sdb;
        }
        // clean up at the end
        if (ii < m) {
            if (m - ii == 1)
                kernel_dgecpsc_1_0_lib4(0, n, alpha, A + 3, B);
            else if (m - ii == 2)
                kernel_dgecpsc_2_3_lib4(0, n, alpha, A, sda, B);
            else  // if(m-ii==3)
                kernel_dgecpsc_3_3_lib4(0, n, alpha, A, sda, B);
        }
    }
}


// copy a lower triangular strmat into a lower triangular strmat
void dtrcp_l(int m, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj) {

    // invalidate stored inverse diagonal
    sB->use_dA = 0;

    const int bs = D_PS;
    const double alpha = 1;

    // get submatrices
    int sda = sA->cn;
    double* A = sA->pA + ai / bs * bs * sda + ai % bs + aj * bs;
    int sdb = sB->cn;
    double* B = sB->pA + bi / bs * bs * sdb + bi % bs + bj * bs;

    if (m <= 0) {
        return;
    }

    int n = m;

    int mna, ii;

    int offA = ai % bs;
    int offB = bi % bs;

    // A at the beginning of the block
    A -= offA;

    // A at the beginning of the block
    B -= offB;

    // same alignment
    if (offA == offB) {
        ii = 0;
        // clean up at the beginning
        mna = (4 - offB) % bs;
        if (mna > 0) {
            if (m < mna)  // mna<=3  ==>  m = { 1, 2 }
            {
                if (m == 1) {
                    kernel_dgecpsc_1_0_lib4(1, ii, alpha, A + offA, B + offB);
                    return;
                } else  // if(m==2 && mna==3)
                {
                    kernel_dgecpsc_2_0_lib4(1, ii, alpha, A + offA, B + offB);
                    return;
                }
            }
            if (mna == 1) {
                kernel_dgecpsc_1_0_lib4(1, ii, alpha, A + offA, B + offB);
                A += 4 * sda;
                B += 4 * sdb;
                ii += 1;
            } else if (mna == 2) {
                kernel_dgecpsc_2_0_lib4(1, ii, alpha, A + offA, B + offB);
                A += 4 * sda;
                B += 4 * sdb;
                ii += 2;
            } else  // if(mna==3)
            {
                kernel_dgecpsc_3_0_lib4(1, ii, alpha, A + offA, B + offB);
                A += 4 * sda;
                B += 4 * sdb;
                ii += 3;
            }
        }
        // main loop
        for (; ii < m - 3; ii += 4) {
            kernel_dgecp_4_0_lib4(1, ii, A, B);
            A += 4 * sda;
            B += 4 * sdb;
        }
        // clean up at the end
        if (ii < m) {
            if (m - ii == 1)
                kernel_dgecpsc_1_0_lib4(1, ii, alpha, A, B);
            else if (m - ii == 2)
                kernel_dgecpsc_2_0_lib4(1, ii, alpha, A, B);
            else  // if(m-ii==3)
                kernel_dgecpsc_3_0_lib4(1, ii, alpha, A, B);
        }
    }
    // skip one element of A
    else if (offA == (offB + 1) % bs) {
        ii = 0;
        // clean up at the beginning
        mna = (4 - offB) % bs;
        if (mna > 0) {
            if (m < mna)  // mna<=3  ==>  m = { 1, 2 }
            {
                if (m == 1) {
                    kernel_dgecpsc_1_0_lib4(1, ii, alpha, A + offA, B + offB);
                    return;
                } else  // if(m==2 && mna==3)
                {
                    kernel_dgecpsc_2_0_lib4(1, ii, alpha, A + offA, B + offB);
                    return;
                }
            }
            if (mna == 1) {
                kernel_dgecpsc_1_0_lib4(1, ii, alpha, A + offA, B + offB);
                // A += 4*sda;
                B += 4 * sdb;
                ii += 1;
            } else if (mna == 2) {
                kernel_dgecpsc_2_3_lib4(1, ii, alpha, A, sda, B + 2);
                A += 4 * sda;
                B += 4 * sdb;
                ii += 2;
            } else  // if(mna==3)
            {
                kernel_dgecpsc_3_2_lib4(1, ii, alpha, A, sda, B + 1);
                A += 4 * sda;
                B += 4 * sdb;
                ii += 3;
            }
        }
        // main loop
        for (; ii < m - 3; ii += 4) {
            kernel_dgecpsc_4_1_lib4(1, ii, alpha, A, sda, B);
            A += 4 * sda;
            B += 4 * sdb;
        }
        // clean up at the end
        if (ii < m) {
            if (m - ii == 1)
                kernel_dgecpsc_1_0_lib4(1, ii, alpha, A + 1, B);
            else if (m - ii == 2)
                kernel_dgecpsc_2_0_lib4(1, ii, alpha, A + 1, B);
            else  // if(m-ii==3)
                kernel_dgecpsc_3_0_lib4(1, ii, alpha, A + 1, B);
        }
    }
    // skip 2 elements of A
    else if (offA == (offB + 2) % bs) {
        ii = 0;
        // clean up at the beginning
        mna = (4 - offB) % bs;
        if (mna > 0) {
            if (m < mna) {
                if (m == 1) {
                    kernel_dgecpsc_1_0_lib4(1, ii, alpha, A + offA, B + offB);
                    return;
                } else  // if(m==2 && mna==3)
                {
                    kernel_dgecpsc_2_3_lib4(1, ii, alpha, A, sda, B + 1);
                    return;
                }
            }
            if (mna == 1) {
                kernel_dgecpsc_1_0_lib4(1, ii, alpha, A + 1, B + 3);
                // A += 4*sda;
                B += 4 * sdb;
                ii += 1;
            } else if (mna == 2) {
                kernel_dgecpsc_2_0_lib4(1, ii, alpha, A, B + 2);
                // A += 4*sda;
                B += 4 * sdb;
                ii += 2;
            } else  // if(mna==3)
            {
                kernel_dgecpsc_3_3_lib4(1, ii, alpha, A, sda, B + 1);
                A += 4 * sda;
                B += 4 * sdb;
                ii += 3;
            }
        }
        // main loop
        for (; ii < m - 3; ii += 4) {
            kernel_dgecpsc_4_2_lib4(1, ii, alpha, A, sda, B);
            A += 4 * sda;
            B += 4 * sdb;
        }
        // clean up at the end
        if (ii < m) {
            if (m - ii == 1)
                kernel_dgecpsc_1_0_lib4(1, ii, alpha, A + 2, B);
            else if (m - ii == 2)
                kernel_dgecpsc_2_0_lib4(1, ii, alpha, A + 2, B);
            else  // if(m-ii==3)
                kernel_dgecpsc_3_2_lib4(1, ii, alpha, A, sda, B);
        }
    }
    // skip 3 elements of A
    else  // if(offA==(offB+3)%bs)
    {
        ii = 0;
        // clean up at the beginning
        mna = (4 - offB) % bs;
        if (mna > 0) {
            if (m < mna) {
                if (m == 1) {
                    kernel_dgecpsc_1_0_lib4(1, ii, alpha, A + offA, B + offB);
                    return;
                } else  // if(m==2 && mna==3)
                {
                    kernel_dgecpsc_2_0_lib4(1, ii, alpha, A + offA, B + offB);
                    return;
                }
            }
            if (mna == 1) {
                kernel_dgecpsc_1_0_lib4(1, ii, alpha, A + offA, B + offB);
                // A += 4*sda;
                B += 4 * sdb;
                ii += 1;
            } else if (mna == 2) {
                kernel_dgecpsc_2_0_lib4(1, ii, alpha, A + offA, B + offB);
                // A += 4*sda;
                B += 4 * sdb;
                ii += 2;
            } else  // if(mna==3)
            {
                kernel_dgecpsc_3_0_lib4(1, ii, alpha, A + offA, B + offB);
                // A += 4*sda;
                B += 4 * sdb;
                ii += 3;
            }
        }
        // main loop
        for (; ii < m - 3; ii += 4) {
            kernel_dgecpsc_4_3_lib4(1, ii, alpha, A, sda, B);
            A += 4 * sda;
            B += 4 * sdb;
        }
        // clean up at the end
        if (ii < m) {
            if (m - ii == 1)
                kernel_dgecpsc_1_0_lib4(1, ii, alpha, A + 3, B);
            else if (m - ii == 2)
                kernel_dgecpsc_2_3_lib4(1, ii, alpha, A, sda, B);
            else  // if(m-ii==3)
                kernel_dgecpsc_3_3_lib4(1, ii, alpha, A, sda, B);
        }
    }

    /* } */
}


// copy and scale a generic strmat into a generic strmat
void dgecpsc(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj) {

    // invalidate stored inverse diagonal
    sB->use_dA = 0;

    const int bs = D_PS;

    // extract dimension
    int sda = sA->cn;
    int sdb = sB->cn;

    // extract submatrix
    double* A = sA->pA + ai / bs * bs * sda + ai % bs + aj * bs;
    double* B = sB->pA + bi / bs * bs * sdb + bi % bs + bj * bs;

    if (m <= 0 || n <= 0) {
        return;
    }

    int mna, ii;

    int offA = ai % bs;
    int offB = bi % bs;

    // A at the beginning of the block
    A -= offA;

    // A at the beginning of the block
    B -= offB;

    // same alignment
    if (offA == offB) {
        ii = 0;
        // clean up at the beginning
        mna = (4 - offB) % bs;
        if (mna > 0) {
            if (m < mna)  // mna<=3  ==>  m = { 1, 2 }
            {
                if (m == 1) {
                    kernel_dgecpsc_1_0_lib4(0, n, alpha, A + offA, B + offB);
                    return;
                } else  // if(m==2 && mna==3)
                {
                    kernel_dgecpsc_2_0_lib4(0, n, alpha, A + offA, B + offB);
                    return;
                }
            }
            if (mna == 1) {
                kernel_dgecpsc_1_0_lib4(0, n, alpha, A + offA, B + offB);
                A += 4 * sda;
                B += 4 * sdb;
                ii += 1;
            } else if (mna == 2) {
                kernel_dgecpsc_2_0_lib4(0, n, alpha, A + offA, B + offB);
                A += 4 * sda;
                B += 4 * sdb;
                ii += 2;
            } else  // if(mna==3)
            {
                kernel_dgecpsc_3_0_lib4(0, n, alpha, A + offA, B + offB);
                A += 4 * sda;
                B += 4 * sdb;
                ii += 3;
            }
        }
        // main loop
        for (; ii < m - 3; ii += 4) {
            kernel_dgecpsc_4_0_lib4(0, n, alpha, A, B);
            A += 4 * sda;
            B += 4 * sdb;
        }
        // clean up at the end
        if (ii < m) {
            if (m - ii == 1)
                kernel_dgecpsc_1_0_lib4(0, n, alpha, A, B);
            else if (m - ii == 2)
                kernel_dgecpsc_2_0_lib4(0, n, alpha, A, B);
            else  // if(m-ii==3)
                kernel_dgecpsc_3_0_lib4(0, n, alpha, A, B);
        }
    }
    // skip one element of A
    else if (offA == (offB + 1) % bs) {
        ii = 0;
        // clean up at the beginning
        mna = (4 - offB) % bs;
        if (mna > 0) {
            if (m < mna)  // mna<=3  ==>  m = { 1, 2 }
            {
                if (m == 1) {
                    kernel_dgecpsc_1_0_lib4(0, n, alpha, A + offA, B + offB);
                    return;
                } else  // if(m==2 && mna==3)
                {
                    kernel_dgecpsc_2_0_lib4(0, n, alpha, A + offA, B + offB);
                    return;
                }
            }
            if (mna == 1) {
                kernel_dgecpsc_1_0_lib4(0, n, alpha, A + offA, B + offB);
                // A += 4*sda;
                B += 4 * sdb;
                ii += 1;
            } else if (mna == 2) {
                kernel_dgecpsc_2_3_lib4(0, n, alpha, A, sda, B + 2);
                A += 4 * sda;
                B += 4 * sdb;
                ii += 2;
            } else  // if(mna==3)
            {
                kernel_dgecpsc_3_2_lib4(0, n, alpha, A, sda, B + 1);
                A += 4 * sda;
                B += 4 * sdb;
                ii += 3;
            }
        }
        // main loop
        for (; ii < m - 3; ii += 4) {
            kernel_dgecpsc_4_1_lib4(0, n, alpha, A, sda, B);
            A += 4 * sda;
            B += 4 * sdb;
        }
        // clean up at the end
        if (ii < m) {
            if (m - ii == 1)
                kernel_dgecpsc_1_0_lib4(0, n, alpha, A + 1, B);
            else if (m - ii == 2)
                kernel_dgecpsc_2_0_lib4(0, n, alpha, A + 1, B);
            else  // if(m-ii==3)
                kernel_dgecpsc_3_0_lib4(0, n, alpha, A + 1, B);
        }
    }
    // skip 2 elements of A
    else if (offA == (offB + 2) % bs) {
        ii = 0;
        // clean up at the beginning
        mna = (4 - offB) % bs;
        if (mna > 0) {
            if (m < mna) {
                if (m == 1) {
                    kernel_dgecpsc_1_0_lib4(0, n, alpha, A + offA, B + offB);
                    return;
                } else  // if(m==2 && mna==3)
                {
                    kernel_dgecpsc_2_3_lib4(0, n, alpha, A, sda, B + 1);
                    return;
                }
            }
            if (mna == 1) {
                kernel_dgecpsc_1_0_lib4(0, n, alpha, A + 1, B + 3);
                // A += 4*sda;
                B += 4 * sdb;
                ii += 1;
            } else if (mna == 2) {
                kernel_dgecpsc_2_0_lib4(0, n, alpha, A, B + 2);
                // A += 4*sda;
                B += 4 * sdb;
                ii += 2;
            } else  // if(mna==3)
            {
                kernel_dgecpsc_3_3_lib4(0, n, alpha, A, sda, B + 1);
                A += 4 * sda;
                B += 4 * sdb;
                ii += 3;
            }
        }
        // main loop
        for (; ii < m - 3; ii += 4) {
            kernel_dgecpsc_4_2_lib4(0, n, alpha, A, sda, B);
            A += 4 * sda;
            B += 4 * sdb;
        }
        // clean up at the end
        if (ii < m) {
            if (m - ii == 1)
                kernel_dgecpsc_1_0_lib4(0, n, alpha, A + 2, B);
            else if (m - ii == 2)
                kernel_dgecpsc_2_0_lib4(0, n, alpha, A + 2, B);
            else  // if(m-ii==3)
                kernel_dgecpsc_3_2_lib4(0, n, alpha, A, sda, B);
        }
    }
    // skip 3 elements of A
    else  // if(offA==(offB+3)%bs)
    {
        ii = 0;
        // clean up at the beginning
        mna = (4 - offB) % bs;
        if (mna > 0) {
            if (m < mna) {
                if (m == 1) {
                    kernel_dgecpsc_1_0_lib4(0, n, alpha, A + offA, B + offB);
                    return;
                } else  // if(m==2 && mna==3)
                {
                    kernel_dgecpsc_2_0_lib4(0, n, alpha, A + offA, B + offB);
                    return;
                }
            }
            if (mna == 1) {
                kernel_dgecpsc_1_0_lib4(0, n, alpha, A + offA, B + offB);
                // A += 4*sda;
                B += 4 * sdb;
                ii += 1;
            } else if (mna == 2) {
                kernel_dgecpsc_2_0_lib4(0, n, alpha, A + offA, B + offB);
                // A += 4*sda;
                B += 4 * sdb;
                ii += 2;
            } else  // if(mna==3)
            {
                kernel_dgecpsc_3_0_lib4(0, n, alpha, A + offA, B + offB);
                // A += 4*sda;
                B += 4 * sdb;
                ii += 3;
            }
        }
        // main loop
        for (; ii < m - 3; ii += 4) {
            kernel_dgecpsc_4_3_lib4(0, n, alpha, A, sda, B);
            A += 4 * sda;
            B += 4 * sdb;
        }
        // clean up at the end
        if (ii < m) {
            if (m - ii == 1)
                kernel_dgecpsc_1_0_lib4(0, n, alpha, A + 3, B);
            else if (m - ii == 2)
                kernel_dgecpsc_2_3_lib4(0, n, alpha, A, sda, B);
            else  // if(m-ii==3)
                kernel_dgecpsc_3_3_lib4(0, n, alpha, A, sda, B);
        }
    }
}


// copy  and scale a lower triangular strmat into a lower triangular strmat
void dtrcpsc_l(int m, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj) {

    // invalidate stored inverse diagonal
    sB->use_dA = 0;

    const int bs = D_PS;
    int sda = sA->cn;
    double* A = sA->pA + ai / bs * bs * sda + ai % bs + aj * bs;
    int sdb = sB->cn;
    double* B = sB->pA + bi / bs * bs * sdb + bi % bs + bj * bs;

    /* dtrcp_l_lib(m, 1.0, ai%bs, pA, sda, ci%bs, pC, sdc); */
    /* // copies a lower triangular packed matrix into a lower triangular packed matrix */
    /* void dtrcp_l_lib(int m, double alpha, int offsetA, double *A, int sda, int offsetB, double *B, int sdb) */
    /* { */

    if (m <= 0) {
        return;
    }

    int n = m;

    int mna, ii;

    int offA = ai % bs;
    int offB = bi % bs;

    // A at the beginning of the block
    A -= offA;

    // A at the beginning of the block
    B -= offB;

    // same alignment
    if (offA == offB) {
        ii = 0;
        // clean up at the beginning
        mna = (4 - offB) % bs;
        if (mna > 0) {
            if (m < mna)  // mna<=3  ==>  m = { 1, 2 }
            {
                if (m == 1) {
                    kernel_dgecpsc_1_0_lib4(1, ii, alpha, A + offA, B + offB);
                    return;
                } else  // if(m==2 && mna==3)
                {
                    kernel_dgecpsc_2_0_lib4(1, ii, alpha, A + offA, B + offB);
                    return;
                }
            }
            if (mna == 1) {
                kernel_dgecpsc_1_0_lib4(1, ii, alpha, A + offA, B + offB);
                A += 4 * sda;
                B += 4 * sdb;
                ii += 1;
            } else if (mna == 2) {
                kernel_dgecpsc_2_0_lib4(1, ii, alpha, A + offA, B + offB);
                A += 4 * sda;
                B += 4 * sdb;
                ii += 2;
            } else  // if(mna==3)
            {
                kernel_dgecpsc_3_0_lib4(1, ii, alpha, A + offA, B + offB);
                A += 4 * sda;
                B += 4 * sdb;
                ii += 3;
            }
        }
        // main loop
        for (; ii < m - 3; ii += 4) {
            kernel_dgecpsc_4_0_lib4(1, ii, alpha, A, B);
            A += 4 * sda;
            B += 4 * sdb;
        }
        // clean up at the end
        if (ii < m) {
            if (m - ii == 1)
                kernel_dgecpsc_1_0_lib4(1, ii, alpha, A, B);
            else if (m - ii == 2)
                kernel_dgecpsc_2_0_lib4(1, ii, alpha, A, B);
            else  // if(m-ii==3)
                kernel_dgecpsc_3_0_lib4(1, ii, alpha, A, B);
        }
    }
    // skip one element of A
    else if (offA == (offB + 1) % bs) {
        ii = 0;
        // clean up at the beginning
        mna = (4 - offB) % bs;
        if (mna > 0) {
            if (m < mna)  // mna<=3  ==>  m = { 1, 2 }
            {
                if (m == 1) {
                    kernel_dgecpsc_1_0_lib4(1, ii, alpha, A + offA, B + offB);
                    return;
                } else  // if(m==2 && mna==3)
                {
                    kernel_dgecpsc_2_0_lib4(1, ii, alpha, A + offA, B + offB);
                    return;
                }
            }
            if (mna == 1) {
                kernel_dgecpsc_1_0_lib4(1, ii, alpha, A + offA, B + offB);
                // A += 4*sda;
                B += 4 * sdb;
                ii += 1;
            } else if (mna == 2) {
                kernel_dgecpsc_2_3_lib4(1, ii, alpha, A, sda, B + 2);
                A += 4 * sda;
                B += 4 * sdb;
                ii += 2;
            } else  // if(mna==3)
            {
                kernel_dgecpsc_3_2_lib4(1, ii, alpha, A, sda, B + 1);
                A += 4 * sda;
                B += 4 * sdb;
                ii += 3;
            }
        }
        // main loop
        for (; ii < m - 3; ii += 4) {
            kernel_dgecpsc_4_1_lib4(1, ii, alpha, A, sda, B);
            A += 4 * sda;
            B += 4 * sdb;
        }
        // clean up at the end
        if (ii < m) {
            if (m - ii == 1)
                kernel_dgecpsc_1_0_lib4(1, ii, alpha, A + 1, B);
            else if (m - ii == 2)
                kernel_dgecpsc_2_0_lib4(1, ii, alpha, A + 1, B);
            else  // if(m-ii==3)
                kernel_dgecpsc_3_0_lib4(1, ii, alpha, A + 1, B);
        }
    }
    // skip 2 elements of A
    else if (offA == (offB + 2) % bs) {
        ii = 0;
        // clean up at the beginning
        mna = (4 - offB) % bs;
        if (mna > 0) {
            if (m < mna) {
                if (m == 1) {
                    kernel_dgecpsc_1_0_lib4(1, ii, alpha, A + offA, B + offB);
                    return;
                } else  // if(m==2 && mna==3)
                {
                    kernel_dgecpsc_2_3_lib4(1, ii, alpha, A, sda, B + 1);
                    return;
                }
            }
            if (mna == 1) {
                kernel_dgecpsc_1_0_lib4(1, ii, alpha, A + 1, B + 3);
                // A += 4*sda;
                B += 4 * sdb;
                ii += 1;
            } else if (mna == 2) {
                kernel_dgecpsc_2_0_lib4(1, ii, alpha, A, B + 2);
                // A += 4*sda;
                B += 4 * sdb;
                ii += 2;
            } else  // if(mna==3)
            {
                kernel_dgecpsc_3_3_lib4(1, ii, alpha, A, sda, B + 1);
                A += 4 * sda;
                B += 4 * sdb;
                ii += 3;
            }
        }
        // main loop
        for (; ii < m - 3; ii += 4) {
            kernel_dgecpsc_4_2_lib4(1, ii, alpha, A, sda, B);
            A += 4 * sda;
            B += 4 * sdb;
        }
        // clean up at the end
        if (ii < m) {
            if (m - ii == 1)
                kernel_dgecpsc_1_0_lib4(1, ii, alpha, A + 2, B);
            else if (m - ii == 2)
                kernel_dgecpsc_2_0_lib4(1, ii, alpha, A + 2, B);
            else  // if(m-ii==3)
                kernel_dgecpsc_3_2_lib4(1, ii, alpha, A, sda, B);
        }
    }
    // skip 3 elements of A
    else  // if(offA==(offB+3)%bs)
    {
        ii = 0;
        // clean up at the beginning
        mna = (4 - offB) % bs;
        if (mna > 0) {
            if (m < mna) {
                if (m == 1) {
                    kernel_dgecpsc_1_0_lib4(1, ii, alpha, A + offA, B + offB);
                    return;
                } else  // if(m==2 && mna==3)
                {
                    kernel_dgecpsc_2_0_lib4(1, ii, alpha, A + offA, B + offB);
                    return;
                }
            }
            if (mna == 1) {
                kernel_dgecpsc_1_0_lib4(1, ii, alpha, A + offA, B + offB);
                // A += 4*sda;
                B += 4 * sdb;
                ii += 1;
            } else if (mna == 2) {
                kernel_dgecpsc_2_0_lib4(1, ii, alpha, A + offA, B + offB);
                // A += 4*sda;
                B += 4 * sdb;
                ii += 2;
            } else  // if(mna==3)
            {
                kernel_dgecpsc_3_0_lib4(1, ii, alpha, A + offA, B + offB);
                // A += 4*sda;
                B += 4 * sdb;
                ii += 3;
            }
        }
        // main loop
        for (; ii < m - 3; ii += 4) {
            kernel_dgecpsc_4_3_lib4(1, ii, alpha, A, sda, B);
            A += 4 * sda;
            B += 4 * sdb;
        }
        // clean up at the end
        if (ii < m) {
            if (m - ii == 1)
                kernel_dgecpsc_1_0_lib4(1, ii, alpha, A + 3, B);
            else if (m - ii == 2)
                kernel_dgecpsc_2_3_lib4(1, ii, alpha, A, sda, B);
            else  // if(m-ii==3)
                kernel_dgecpsc_3_3_lib4(1, ii, alpha, A, sda, B);
        }
    }

    /* } */
}


// scale a generic strmat
void dgesc(int m, int n, double alpha, struct mat* sA, int ai, int aj) {
    // invalidate stored inverse diagonal
    sA->use_dA = 0;

    dgecpsc(m, n, alpha, sA, ai, aj, sA, ai, aj);
}


// scale a triangular strmat
void dtrsc_l(int m, double alpha, struct mat* sA, int ai, int aj) {
    // invalidate stored inverse diagonal
    sA->use_dA = 0;

    dtrcpsc_l(m, alpha, sA, ai, aj, sA, ai, aj);
}


// copy a strvec into a strvec
void dveccp(int m, struct vec* sa, int ai, struct vec* sc, int ci) {
    double* pa = sa->pa + ai;
    double* pc = sc->pa + ci;
    kernel_dveccp_inc1(m, pa, pc);
}


// scale a strvec
void dvecsc(int m, double alpha, struct vec* sa, int ai) {
    double* pa = sa->pa + ai;
    int ii;
    ii = 0;
    for (; ii < m - 3; ii += 4) {
        pa[ii + 0] *= alpha;
        pa[ii + 1] *= alpha;
        pa[ii + 2] *= alpha;
        pa[ii + 3] *= alpha;
    }
    for (; ii < m; ii++) {
        pa[ii + 0] *= alpha;
    }
}


// copy and scale a strvec into a strvec
void dveccpsc(int m, double alpha, struct vec* sa, int ai, struct vec* sc, int ci) {
    double* pa = sa->pa + ai;
    double* pc = sc->pa + ci;
    int ii;
    ii = 0;
    for (; ii < m - 3; ii += 4) {
        pc[ii + 0] = alpha * pa[ii + 0];
        pc[ii + 1] = alpha * pa[ii + 1];
        pc[ii + 2] = alpha * pa[ii + 2];
        pc[ii + 3] = alpha * pa[ii + 3];
    }
    for (; ii < m; ii++) {
        pc[ii + 0] = alpha * pa[ii + 0];
    }
}

// scales and adds a packed matrix into a packed matrix: B = B + alpha*A
void dgead_lib(int m, int n, double alpha, int offsetA, double* A, int sda, int offsetB, double* B, int sdb) {
    if (m <= 0 || n <= 0) {
        return;
    }

    const int bs = D_PS;

    int mna, ii;

    int offA = offsetA % bs;
    int offB = offsetB % bs;

    // A at the beginning of the block
    A -= offA;

    // A at the beginning of the block
    B -= offB;

    // same alignment
    if (offA == offB) {
        ii = 0;
        // clean up at the beginning
        mna = (4 - offB) % bs;
        if (mna > 0) {
            if (m < mna)  // mna<=3  ==>  m = { 1, 2 }
            {
                if (m == 1) {
                    kernel_dgead_1_0_lib4(n, alpha, A + offA, B + offB);
                    return;
                } else  // if(m==2 && mna==3)
                {
                    kernel_dgead_2_0_lib4(n, alpha, A + offA, B + offB);
                    return;
                }
            }
            if (mna == 1) {
                kernel_dgead_1_0_lib4(n, alpha, A + offA, B + offB);
                A += 4 * sda;
                B += 4 * sdb;
                ii += 1;
            } else if (mna == 2) {
                kernel_dgead_2_0_lib4(n, alpha, A + offA, B + offB);
                A += 4 * sda;
                B += 4 * sdb;
                ii += 2;
            } else  // if(mna==3)
            {
                kernel_dgead_3_0_lib4(n, alpha, A + offA, B + offB);
                A += 4 * sda;
                B += 4 * sdb;
                ii += 3;
            }
        }
        // main loop
        for (; ii < m - 3; ii += 4) {
            kernel_dgead_4_0_lib4(n, alpha, A, B);
            A += 4 * sda;
            B += 4 * sdb;
        }
        // clean up at the end
        if (ii < m) {
            if (m - ii == 1)
                kernel_dgead_1_0_lib4(n, alpha, A, B);
            else if (m - ii == 2)
                kernel_dgead_2_0_lib4(n, alpha, A, B);
            else  // if(m-ii==3)
                kernel_dgead_3_0_lib4(n, alpha, A, B);
        }
    }
    // skip one element of A
    else if (offA == (offB + 1) % bs) {
        ii = 0;
        // clean up at the beginning
        mna = (4 - offB) % bs;
        if (mna > 0) {
            if (m < mna)  // mna<=3  ==>  m = { 1, 2 }
            {
                if (m == 1) {
                    kernel_dgead_1_0_lib4(n, alpha, A + offA, B + offB);
                    return;
                } else  // if(m==2 && mna==3)
                {
                    kernel_dgead_2_0_lib4(n, alpha, A + offA, B + offB);
                    return;
                }
            }
            if (mna == 1) {
                kernel_dgead_1_0_lib4(n, alpha, A + offA, B + offB);
                // A += 4*sda;
                B += 4 * sdb;
                ii += 1;
            } else if (mna == 2) {
                kernel_dgead_2_3_lib4(n, alpha, A, sda, B + 2);
                A += 4 * sda;
                B += 4 * sdb;
                ii += 2;
            } else  // if(mna==3)
            {
                kernel_dgead_3_2_lib4(n, alpha, A, sda, B + 1);
                A += 4 * sda;
                B += 4 * sdb;
                ii += 3;
            }
        }
        // main loop
        for (; ii < m - 3; ii += 4) {
            kernel_dgead_4_1_lib4(n, alpha, A, sda, B);
            A += 4 * sda;
            B += 4 * sdb;
        }
        // clean up at the end
        if (ii < m) {
            if (m - ii == 1)
                kernel_dgead_1_0_lib4(n, alpha, A + 1, B);
            else if (m - ii == 2)
                kernel_dgead_2_0_lib4(n, alpha, A + 1, B);
            else  // if(m-ii==3)
                kernel_dgead_3_0_lib4(n, alpha, A + 1, B);
        }
    }
    // skip 2 elements of A
    else if (offA == (offB + 2) % bs) {
        ii = 0;
        // clean up at the beginning
        mna = (4 - offB) % bs;
        if (mna > 0) {
            if (m < mna) {
                if (m == 1) {
                    kernel_dgead_1_0_lib4(n, alpha, A + offA, B + offB);
                    return;
                } else  // if(m==2 && mna==3)
                {
                    kernel_dgead_2_3_lib4(n, alpha, A, sda, B + 1);
                    return;
                }
            }
            if (mna == 1) {
                kernel_dgead_1_0_lib4(n, alpha, A + 1, B + 3);
                // A += 4*sda;
                B += 4 * sdb;
                ii += 1;
            } else if (mna == 2) {
                kernel_dgead_2_0_lib4(n, alpha, A, B + 2);
                // A += 4*sda;
                B += 4 * sdb;
                ii += 2;
            } else  // if(mna==3)
            {
                kernel_dgead_3_3_lib4(n, alpha, A, sda, B + 1);
                A += 4 * sda;
                B += 4 * sdb;
                ii += 3;
            }
        }
        // main loop
        for (; ii < m - 3; ii += 4) {
            kernel_dgead_4_2_lib4(n, alpha, A, sda, B);
            A += 4 * sda;
            B += 4 * sdb;
        }
        // clean up at the end
        if (ii < m) {
            if (m - ii == 1)
                kernel_dgead_1_0_lib4(n, alpha, A + 2, B);
            else if (m - ii == 2)
                kernel_dgead_2_0_lib4(n, alpha, A + 2, B);
            else  // if(m-ii==3)
                kernel_dgead_3_2_lib4(n, alpha, A, sda, B);
        }
    }
    // skip 3 elements of A
    else  // if(offA==(offB+3)%bs)
    {
        ii = 0;
        // clean up at the beginning
        mna = (4 - offB) % bs;
        if (mna > 0) {
            if (m < mna) {
                if (m == 1) {
                    kernel_dgead_1_0_lib4(n, alpha, A + offA, B + offB);
                    return;
                } else  // if(m==2 && mna==3)
                {
                    kernel_dgead_2_0_lib4(n, alpha, A + offA, B + offB);
                    return;
                }
            }
            if (mna == 1) {
                kernel_dgead_1_0_lib4(n, alpha, A + offA, B + offB);
                // A += 4*sda;
                B += 4 * sdb;
                ii += 1;
            } else if (mna == 2) {
                kernel_dgead_2_0_lib4(n, alpha, A + offA, B + offB);
                // A += 4*sda;
                B += 4 * sdb;
                ii += 2;
            } else  // if(mna==3)
            {
                kernel_dgead_3_0_lib4(n, alpha, A + offA, B + offB);
                // A += 4*sda;
                B += 4 * sdb;
                ii += 3;
            }
        }
        // main loop
        for (; ii < m - 3; ii += 4) {
            kernel_dgead_4_3_lib4(n, alpha, A, sda, B);
            A += 4 * sda;
            B += 4 * sdb;
        }
        // clean up at the end
        if (ii < m) {
            if (m - ii == 1)
                kernel_dgead_1_0_lib4(n, alpha, A + 3, B);
            else if (m - ii == 2)
                kernel_dgead_2_3_lib4(n, alpha, A, sda, B);
            else  // if(m-ii==3)
                kernel_dgead_3_3_lib4(n, alpha, A, sda, B);
        }
    }
}
// scale and add a generic strmat into a generic strmat
void dgead(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sC, int ci, int cj) {
    // invalidate stored inverse diagonal
    sC->use_dA = 0;

    const int bs = D_PS;
    int sda = sA->cn;
    double* pA = sA->pA + ai / bs * bs * sda + ai % bs + aj * bs;
    int sdc = sC->cn;
    double* pC = sC->pA + ci / bs * bs * sdc + ci % bs + cj * bs;
    dgead_lib(m, n, alpha, ai % bs, pA, sda, ci % bs, pC, sdc);
}


// scales and adds a strvec into a strvec
void dvecad(int m, double alpha, struct vec* sa, int ai, struct vec* sc, int ci) {
    double* pa = sa->pa + ai;
    double* pc = sc->pa + ci;
    int ii;
    ii = 0;
    for (; ii < m - 3; ii += 4) {
        pc[ii + 0] += alpha * pa[ii + 0];
        pc[ii + 1] += alpha * pa[ii + 1];
        pc[ii + 2] += alpha * pa[ii + 2];
        pc[ii + 3] += alpha * pa[ii + 3];
    }
    for (; ii < m; ii++) {
        pc[ii + 0] += alpha * pa[ii + 0];
    }
}

// transpose general matrix; m and n are referred to the original matrix
void dgetr_lib(int m, int n, double alpha, int offsetA, double* pA, int sda, int offsetC, double* pC, int sdc) {

    /*

    m = 5
    n = 3
    offsetA = 1
    offsetC = 2

    A =
     x x x
     -
     x x x
     x x x
     x x x
     x x x

    C =
     x x x x x
     x x x x x
     -
     x x x x x

    */


    if (m <= 0 || n <= 0) {
        return;
    }

    const int bs = D_PS;

    int mna = (bs - offsetA % bs) % bs;
    mna = m < mna ? m : mna;
    int nna = (bs - offsetC % bs) % bs;
    nna = n < nna ? n : nna;

    int ii;

    ii = 0;

    if (mna > 0) {
        if (mna == 1)
            kernel_dgetr_1_lib4(0, n, nna, alpha, pA, pC, sdc);
        else if (mna == 2)
            kernel_dgetr_2_lib4(0, n, nna, alpha, pA, pC, sdc);
        else  // if(mna==3)
            kernel_dgetr_3_lib4(0, n, nna, alpha, pA, pC, sdc);
        ii += mna;
        pA += mna + bs * (sda - 1);
        pC += mna * bs;
    }
    for (; ii < m - 3; ii += 4)
    //	for( ; ii<m; ii+=4)
    {
        kernel_dgetr_4_lib4(0, n, nna, alpha, pA, pC, sdc);
        pA += bs * sda;
        pC += bs * bs;
    }

    // clean-up at the end using smaller kernels
    if (ii == m) {
        return;
    }

    if (m - ii == 1)
        kernel_dgetr_1_lib4(0, n, nna, alpha, pA, pC, sdc);
    else if (m - ii == 2)
        kernel_dgetr_2_lib4(0, n, nna, alpha, pA, pC, sdc);
    else if (m - ii == 3)
        kernel_dgetr_3_lib4(0, n, nna, alpha, pA, pC, sdc);
}
// copy and transpose a generic strmat into a generic strmat
void dgetr(int m, int n, struct mat* sA, int ai, int aj, struct mat* sC, int ci, int cj) {
    // invalidate stored inverse diagonal
    sC->use_dA = 0;

    const int bs = D_PS;
    int sda = sA->cn;
    double* pA = sA->pA + ai / bs * bs * sda + ai % bs + aj * bs;
    int sdc = sC->cn;
    double* pC = sC->pA + ci / bs * bs * sdc + ci % bs + cj * bs;
    dgetr_lib(m, n, 1.0, ai % bs, pA, sda, ci % bs, pC, sdc);  // TODO remove alpha !!!
}

// transpose lower triangular matrix
void dtrtr_l_lib(int m, double alpha, int offsetA, double* pA, int sda, int offsetC, double* pC, int sdc) {

    /*

    A =
     x
     x x
     x x x
     x x x x

     x x x x x
     x x x x x x
     x x x x x x x
     x x x x x x x x

    C =
     x x x x x x x x

       x x x x x x x
         x x x x x x
           x x x x x
             x x x x

               x x x
                 x x
                   x

    */

    int n = m;

    if (m <= 0 || n <= 0) {
        return;
    }

    const int bs = D_PS;

    int mna = (bs - offsetA % bs) % bs;
    mna = m < mna ? m : mna;
    int nna = (bs - offsetC % bs) % bs;
    nna = n < nna ? n : nna;

    int ii;

    ii = 0;

    if (mna > 0) {
        if (mna == 1) {
            pC[0] = alpha * pA[0];
        } else if (mna == 2) {
            if (nna == 1) {
                pC[0 + bs * 0] = alpha * pA[0 + bs * 0];
                pC[0 + bs * 1] = alpha * pA[1 + bs * 0];
                pC[1 + bs * (0 + sdc)] = alpha * pA[1 + bs * 1];
            } else {
                pC[0 + bs * 0] = alpha * pA[0 + bs * 0];
                pC[0 + bs * 1] = alpha * pA[1 + bs * 0];
                pC[1 + bs * 1] = alpha * pA[1 + bs * 1];
            }
        } else  // if(mna==3)
        {
            if (nna == 1) {
                pC[0 + bs * 0] = alpha * pA[0 + bs * 0];
                pC[0 + bs * 1] = alpha * pA[1 + bs * 0];
                pC[0 + bs * 2] = alpha * pA[2 + bs * 0];
                pC[1 + bs * (0 + sdc)] = alpha * pA[1 + bs * 1];
                pC[1 + bs * (1 + sdc)] = alpha * pA[2 + bs * 1];
                pC[2 + bs * (1 + sdc)] = alpha * pA[2 + bs * 2];
            } else if (nna == 2) {
                pC[0 + bs * 0] = alpha * pA[0 + bs * 0];
                pC[0 + bs * 1] = alpha * pA[1 + bs * 0];
                pC[0 + bs * 2] = alpha * pA[2 + bs * 0];
                pC[1 + bs * 1] = alpha * pA[1 + bs * 1];
                pC[1 + bs * 2] = alpha * pA[2 + bs * 1];
                pC[2 + bs * (1 + sdc)] = alpha * pA[2 + bs * 2];
            } else {
                pC[0 + bs * 0] = alpha * pA[0 + bs * 0];
                pC[0 + bs * 1] = alpha * pA[1 + bs * 0];
                pC[0 + bs * 2] = alpha * pA[2 + bs * 0];
                pC[1 + bs * 1] = alpha * pA[1 + bs * 1];
                pC[1 + bs * 2] = alpha * pA[2 + bs * 1];
                pC[2 + bs * 2] = alpha * pA[2 + bs * 2];
            }
        }
        ii += mna;
        pA += mna + bs * (sda - 1);
        pC += mna * bs;
    }
    for (; ii < m - 3; ii += 4) {
        kernel_dgetr_4_lib4(1, ii, nna, alpha, pA, pC, sdc);
        pA += bs * sda;
        pC += bs * bs;
    }

    // clean-up at the end using smaller kernels
    if (ii == m) {
        return;
    }

    if (m - ii == 1)
        kernel_dgetr_1_lib4(1, ii, nna, alpha, pA, pC, sdc);
    else if (m - ii == 2)
        kernel_dgetr_2_lib4(1, ii, nna, alpha, pA, pC, sdc);
    else if (m - ii == 3)
        kernel_dgetr_3_lib4(1, ii, nna, alpha, pA, pC, sdc);
}
// copy and transpose a lower triangular strmat into an upper triangular strmat
void dtrtr_l(int m, struct mat* sA, int ai, int aj, struct mat* sC, int ci, int cj) {
    // invalidate stored inverse diagonal
    sC->use_dA = 0;

    const int bs = D_PS;
    int sda = sA->cn;
    double* pA = sA->pA + ai / bs * bs * sda + ai % bs + aj * bs;
    int sdc = sC->cn;
    double* pC = sC->pA + ci / bs * bs * sdc + ci % bs + cj * bs;
    dtrtr_l_lib(m, 1.0, ai % bs, pA, sda, ci % bs, pC, sdc);  // TODO remove alpha !!!
}

// transpose an aligned upper triangular matrix into an aligned lower triangular matrix
void dtrtr_u_lib(int m, double alpha, int offsetA, double* pA, int sda, int offsetC, double* pC, int sdc) {

    /*

    A =
     x x x x x x x x
       x x x x x x x

         x x x x x x
           x x x x x
             x x x x
               x x x
                 x x
                   x

    C =
     x

     x x
     x x x
     x x x x
     x x x x x
     x x x x x x
     x x x x x x x
     x x x x x x x x

    */

    int n = m;

    if (m <= 0 || n <= 0) {
        return;
    }

    const int bs = D_PS;

    int mna = (bs - offsetA % bs) % bs;
    mna = m < mna ? m : mna;
    int nna = (bs - offsetC % bs) % bs;
    nna = n < nna ? n : nna;
    int tna = nna;

    int ii;

    ii = 0;

    if (mna > 0) {
        if (mna == 1) {
            kernel_dgetr_1_lib4(0, n, nna, alpha, pA, pC, sdc);
            if (nna != 1) {
                //				pC[0+bs*0] = alpha * pA[0+bs*0];
                pA += 1 * bs;
                pC += 1;
                tna = (bs - (offsetC + 1) % bs) % bs;
            } else  // if(nna==1)
            {
                //				pC[0+bs*0] = alpha * pA[0+bs*0];
                pA += 1 * bs;
                pC += 1 + (sdc - 1) * bs;
                tna = 0;  //(bs-(offsetC+1)%bs)%bs;
            }
            //			kernel_dgetr_1_lib4(0, n-1, tna, alpha, pA, pC, sdc);
        } else if (mna == 2) {
            if (nna == 0 || nna == 3) {
                pC[0 + bs * 0] = alpha * pA[0 + bs * 0];
                pC[1 + bs * 0] = alpha * pA[0 + bs * 1];
                pC[1 + bs * 1] = alpha * pA[1 + bs * 1];
                pA += 2 * bs;
                pC += 2;
                tna = (bs - (offsetC + 2) % bs) % bs;
                kernel_dgetr_2_lib4(0, n - 2, tna, alpha, pA, pC, sdc);
            } else if (nna == 1) {
                pC[0 + bs * 0] = alpha * pA[0 + bs * 0];
                pA += 1 * bs;
                pC += 1 + (sdc - 1) * bs;
                //				pC[0+bs*0] = alpha * pA[0+bs*0];
                //				pC[0+bs*1] = alpha * pA[1+bs*0];
                kernel_dgetr_2_lib4(0, n - 1, 0, alpha, pA, pC, sdc);
                pA += 1 * bs;
                pC += 1;
                tna = 3;  //(bs-(offsetC+2)%bs)%bs;
                //				kernel_dgetr_2_lib4(0, n-2, tna, alpha, pA, pC, sdc);
            } else if (nna == 2) {
                pC[0 + bs * 0] = alpha * pA[0 + bs * 0];
                pC[1 + bs * 0] = alpha * pA[0 + bs * 1];
                pC[1 + bs * 1] = alpha * pA[1 + bs * 1];
                pA += 2 * bs;
                pC += 2 + (sdc - 1) * bs;
                tna = 0;  //(bs-(offsetC+2)%bs)%bs;
                kernel_dgetr_2_lib4(0, n - 2, tna, alpha, pA, pC, sdc);
            }
        } else  // if(mna==3)
        {
            if (nna == 0) {
                pC[0 + bs * 0] = alpha * pA[0 + bs * 0];
                pC[1 + bs * 0] = alpha * pA[0 + bs * 1];
                pC[1 + bs * 1] = alpha * pA[1 + bs * 1];
                pC[2 + bs * 0] = alpha * pA[0 + bs * 2];
                pC[2 + bs * 1] = alpha * pA[1 + bs * 2];
                pC[2 + bs * 2] = alpha * pA[2 + bs * 2];
                pA += 3 * bs;
                pC += 3;
                tna = 1;
                kernel_dgetr_3_lib4(0, n - 3, tna, alpha, pA, pC, sdc);
            } else if (nna == 1) {
                pC[0 + bs * 0] = alpha * pA[0 + bs * 0];
                pA += bs;
                pC += 1 + (sdc - 1) * bs;
                pC[0 + bs * 0] = alpha * pA[0 + bs * 0];
                pC[0 + bs * 1] = alpha * pA[1 + bs * 0];
                pC[1 + bs * 0] = alpha * pA[0 + bs * 1];
                pC[1 + bs * 1] = alpha * pA[1 + bs * 1];
                pC[1 + bs * 2] = alpha * pA[2 + bs * 1];
                pA += 2 * bs;
                pC += 2;
                tna = 2;
                kernel_dgetr_3_lib4(0, n - 3, tna, alpha, pA, pC, sdc);
            } else if (nna == 2) {
                pC[0 + bs * 0] = alpha * pA[0 + bs * 0];
                pC[1 + bs * 0] = alpha * pA[0 + bs * 1];
                pC[1 + bs * 1] = alpha * pA[1 + bs * 1];
                pA += 2 * bs;
                pC += 2 + (sdc - 1) * bs;
                //				pC[0+bs*0] = alpha * pA[0+bs*0];
                //				pC[0+bs*1] = alpha * pA[1+bs*0];
                //				pC[0+bs*2] = alpha * pA[2+bs*0];
                kernel_dgetr_3_lib4(0, n - 2, 0, alpha, pA, pC, sdc);
                pA += 1 * bs;
                pC += 1;
                tna = 3;
                //				kernel_dgetr_3_lib4(0, n-3, tna, alpha, pA, pC, sdc);
            } else  // if(nna==3)
            {
                pC[0 + bs * 0] = alpha * pA[0 + bs * 0];
                pC[1 + bs * 0] = alpha * pA[0 + bs * 1];
                pC[1 + bs * 1] = alpha * pA[1 + bs * 1];
                pC[2 + bs * 0] = alpha * pA[0 + bs * 2];
                pC[2 + bs * 1] = alpha * pA[1 + bs * 2];
                pC[2 + bs * 2] = alpha * pA[2 + bs * 2];
                pA += 3 * bs;
                pC += 3 + (sdc - 1) * bs;
                tna = 0;
                kernel_dgetr_3_lib4(0, n - 3, tna, alpha, pA, pC, sdc);
            }
        }
        ii += mna;
        pA += mna + bs * (sda - 1);
        pC += mna * bs;
    }
    for (; ii < m - 3; ii += 4) {
        if (tna == 0) {
            pC[0 + bs * 0] = alpha * pA[0 + bs * 0];
            pC[1 + bs * 0] = alpha * pA[0 + bs * 1];
            pC[1 + bs * 1] = alpha * pA[1 + bs * 1];
            pC[2 + bs * 0] = alpha * pA[0 + bs * 2];
            pC[2 + bs * 1] = alpha * pA[1 + bs * 2];
            pC[2 + bs * 2] = alpha * pA[2 + bs * 2];
            pC[3 + bs * 0] = alpha * pA[0 + bs * 3];
            pC[3 + bs * 1] = alpha * pA[1 + bs * 3];
            pC[3 + bs * 2] = alpha * pA[2 + bs * 3];
            pC[3 + bs * 3] = alpha * pA[3 + bs * 3];
            pA += 4 * bs;
            pC += sdc * bs;
            kernel_dgetr_4_lib4(0, n - ii - 4, 0, alpha, pA, pC, sdc);
        } else if (tna == 1) {
            pC[0 + bs * 0] = alpha * pA[0 + bs * 0];
            pA += bs;
            pC += 1 + (sdc - 1) * bs;
            pC[0 + bs * 0] = alpha * pA[0 + bs * 0];
            pC[0 + bs * 1] = alpha * pA[1 + bs * 0];
            pC[1 + bs * 0] = alpha * pA[0 + bs * 1];
            pC[1 + bs * 1] = alpha * pA[1 + bs * 1];
            pC[1 + bs * 2] = alpha * pA[2 + bs * 1];
            pC[2 + bs * 0] = alpha * pA[0 + bs * 2];
            pC[2 + bs * 1] = alpha * pA[1 + bs * 2];
            pC[2 + bs * 2] = alpha * pA[2 + bs * 2];
            pC[2 + bs * 3] = alpha * pA[3 + bs * 2];
            pA += 3 * bs;
            pC += 3;
            kernel_dgetr_4_lib4(0, n - ii - 4, 1, alpha, pA, pC, sdc);
        } else if (tna == 2) {
            pC[0 + bs * 0] = alpha * pA[0 + bs * 0];
            pC[1 + bs * 0] = alpha * pA[0 + bs * 1];
            pC[1 + bs * 1] = alpha * pA[1 + bs * 1];
            pA += 2 * bs;
            pC += 2 + (sdc - 1) * bs;
            pC[0 + bs * 0] = alpha * pA[0 + bs * 0];
            pC[0 + bs * 1] = alpha * pA[1 + bs * 0];
            pC[0 + bs * 2] = alpha * pA[2 + bs * 0];
            pC[1 + bs * 0] = alpha * pA[0 + bs * 1];
            pC[1 + bs * 1] = alpha * pA[1 + bs * 1];
            pC[1 + bs * 2] = alpha * pA[2 + bs * 1];
            pC[1 + bs * 3] = alpha * pA[3 + bs * 1];
            pA += 2 * bs;
            pC += 2;
            kernel_dgetr_4_lib4(0, n - ii - 4, 2, alpha, pA, pC, sdc);
        } else  // if(tna==3)
        {
            pC[0 + bs * 0] = alpha * pA[0 + bs * 0];
            pC[1 + bs * 0] = alpha * pA[0 + bs * 1];
            pC[1 + bs * 1] = alpha * pA[1 + bs * 1];
            pC[2 + bs * 0] = alpha * pA[0 + bs * 2];
            pC[2 + bs * 1] = alpha * pA[1 + bs * 2];
            pC[2 + bs * 2] = alpha * pA[2 + bs * 2];
            pA += 3 * bs;
            pC += 3 + (sdc - 1) * bs;
            kernel_dgetr_4_lib4(0, n - ii - 3, 0, alpha, pA, pC, sdc);
            //			pC[0+bs*0] = alpha * pA[0+bs*0];
            //			pC[0+bs*1] = alpha * pA[1+bs*0];
            //			pC[0+bs*2] = alpha * pA[2+bs*0];
            //			pC[0+bs*3] = alpha * pA[3+bs*0];
            pA += bs;
            pC += 1;
            //			kernel_dgetr_4_lib4(0, n-ii-4, tna, alpha, pA, pC, sdc);
        }
        pA += bs * sda;
        pC += bs * bs;
    }

    // clean-up at the end
    if (ii == m) {
        return;
    }

    if (m - ii == 1) {
        pC[0 + bs * 0] = alpha * pA[0 + bs * 0];
    } else if (m - ii == 2) {
        if (tna != 1) {
            pC[0 + bs * 0] = alpha * pA[0 + bs * 0];
            pC[1 + bs * 0] = alpha * pA[0 + bs * 1];
            pC[1 + bs * 1] = alpha * pA[1 + bs * 1];
        } else  // if(tna==1)
        {
            pC[0 + bs * 0] = alpha * pA[0 + bs * 0];
            pA += bs;
            pC += 1 + (sdc - 1) * bs;
            pC[0 + bs * 0] = alpha * pA[0 + bs * 0];
            pC[0 + bs * 1] = alpha * pA[1 + bs * 0];
        }
    } else if (m - ii == 3) {
        if (tna == 0 || tna == 3) {
            pC[0 + bs * 0] = alpha * pA[0 + bs * 0];
            pC[1 + bs * 0] = alpha * pA[0 + bs * 1];
            pC[1 + bs * 1] = alpha * pA[1 + bs * 1];
            pC[2 + bs * 0] = alpha * pA[0 + bs * 2];
            pC[2 + bs * 1] = alpha * pA[1 + bs * 2];
            pC[2 + bs * 2] = alpha * pA[2 + bs * 2];
        } else if (tna == 1) {
            pC[0 + bs * 0] = alpha * pA[0 + bs * 0];
            pA += bs;
            pC += 1 + (sdc - 1) * bs;
            pC[0 + bs * 0] = alpha * pA[0 + bs * 0];
            pC[0 + bs * 1] = alpha * pA[1 + bs * 0];
            pC[1 + bs * 0] = alpha * pA[0 + bs * 0];
            pC[1 + bs * 1] = alpha * pA[1 + bs * 1];
            pC[1 + bs * 2] = alpha * pA[2 + bs * 1];
        } else  // if(tna==2)
        {
            pC[0 + bs * 0] = alpha * pA[0 + bs * 0];
            pC[1 + bs * 0] = alpha * pA[0 + bs * 1];
            pC[1 + bs * 1] = alpha * pA[1 + bs * 1];
            pA += 2 * bs;
            pC += 2 + (sdc - 1) * bs;
            pC[0 + bs * 0] = alpha * pA[0 + bs * 0];
            pC[0 + bs * 1] = alpha * pA[1 + bs * 0];
            pC[0 + bs * 2] = alpha * pA[2 + bs * 0];
        }
    }
}
// copy and transpose an upper triangular strmat into a lower triangular strmat
void dtrtr_u(int m, struct mat* sA, int ai, int aj, struct mat* sC, int ci, int cj) {
    // invalidate stored inverse diagonal
    sC->use_dA = 0;

    const int bs = D_PS;
    int sda = sA->cn;
    double* pA = sA->pA + ai / bs * bs * sda + ai % bs + aj * bs;
    int sdc = sC->cn;
    double* pC = sC->pA + ci / bs * bs * sdc + ci % bs + cj * bs;
    dtrtr_u_lib(m, 1.0, ai % bs, pA, sda, ci % bs, pC, sdc);  // TODO remove alpha !!!
}


// insert a strvec to diagonal of strmat, sparse formulation
void ddiain_sp(int kmax, double alpha, struct vec* sx, int xi, int* idx, struct mat* sD, int di, int dj) {
    // invalidate stored inverse diagonal
    sD->use_dA = 0;

    const int bs = D_PS;
    double* x = sx->pa + xi;
    int sdd = sD->cn;
    double* pD = sD->pA;
    int ii, jj;
    for (jj = 0; jj < kmax; jj++) {
        ii = idx[jj];
        pD[(ii + di) / bs * bs * sdd + (ii + di) % bs + (ii + dj) * bs] = alpha * x[jj];
    }
}

// extract diagonal to vector
void ddiaex_lib(int kmax, double alpha, int offset, double* pD, int sdd, double* x) {

    const int bs = D_PS;

    int kna = (bs - offset % bs) % bs;
    kna = kmax < kna ? kmax : kna;

    int jj, ll;

    if (kna > 0) {
        for (ll = 0; ll < kna; ll++) {
            x[ll] = alpha * pD[ll + bs * ll];
        }
        pD += kna + bs * (sdd - 1) + kna * bs;
        x += kna;
        kmax -= kna;
    }
    for (jj = 0; jj < kmax - 3; jj += 4) {
        x[jj + 0] = alpha * pD[jj * sdd + (jj + 0) * bs + 0];
        x[jj + 1] = alpha * pD[jj * sdd + (jj + 1) * bs + 1];
        x[jj + 2] = alpha * pD[jj * sdd + (jj + 2) * bs + 2];
        x[jj + 3] = alpha * pD[jj * sdd + (jj + 3) * bs + 3];
    }
    for (ll = 0; ll < kmax - jj; ll++) {
        x[jj + ll] = alpha * pD[jj * sdd + (jj + ll) * bs + ll];
    }
}
// extract a vector from diagonal
void ddiaex(int kmax, double alpha, struct mat* sA, int ai, int aj, struct vec* sx, int xi) {
    const int bs = D_PS;
    int sda = sA->cn;
    double* pA = sA->pA + ai / bs * bs * sda + ai % bs + aj * bs;
    double* x = sx->pa + xi;
    int offsetA = ai % bs;

    int kna = (bs - offsetA % bs) % bs;
    kna = kmax < kna ? kmax : kna;

    int jj, ll;

    if (kna > 0) {
        for (ll = 0; ll < kna; ll++) {
            x[ll] = alpha * pA[ll + bs * ll];
        }
        pA += kna + bs * (sda - 1) + kna * bs;
        x += kna;
        kmax -= kna;
    }
    for (jj = 0; jj < kmax - 3; jj += 4) {
        x[jj + 0] = alpha * pA[jj * sda + (jj + 0) * bs + 0];
        x[jj + 1] = alpha * pA[jj * sda + (jj + 1) * bs + 1];
        x[jj + 2] = alpha * pA[jj * sda + (jj + 2) * bs + 2];
        x[jj + 3] = alpha * pA[jj * sda + (jj + 3) * bs + 3];
    }
    for (ll = 0; ll < kmax - jj; ll++) {
        x[jj + ll] = alpha * pA[jj * sda + (jj + ll) * bs + ll];
    }
}

// extract diagonal to vector, sparse formulation
void ddiaex_libsp(int kmax, int* idx, double alpha, double* pD, int sdd, double* x) {

    const int bs = D_PS;

    int ii, jj;

    for (jj = 0; jj < kmax; jj++) {
        ii = idx[jj];
        x[jj] = alpha * pD[ii / bs * bs * sdd + ii % bs + ii * bs];
    }
}
// extract the diagonal of a strmat to a strvec, sparse formulation
void ddiaex_sp(int kmax, double alpha, int* idx, struct mat* sD, int di, int dj, struct vec* sx, int xi) {
    const int bs = D_PS;
    double* x = sx->pa + xi;
    int sdd = sD->cn;
    double* pD = sD->pA;
    int ii, jj;
    for (jj = 0; jj < kmax; jj++) {
        ii = idx[jj];
        x[jj] = alpha * pD[(ii + di) / bs * bs * sdd + (ii + di) % bs + (ii + dj) * bs];
    }
}


// add a vector to diagonal
void ddiaad(int kmax, double alpha, struct vec* sx, int xi, struct mat* sA, int ai, int aj) {

    // invalidate stored inverse diagonal
    sA->use_dA = 0;

    const int bs = D_PS;
    int sda = sA->cn;
    double* pA = sA->pA + ai / bs * bs * sda + ai % bs + aj * bs;
    double* x = sx->pa + xi;
    int offsetA = ai % bs;

    int kna = (bs - offsetA % bs) % bs;
    kna = kmax < kna ? kmax : kna;

    int jj, ll;

    if (kna > 0) {
        for (ll = 0; ll < kna; ll++) {
            pA[ll + bs * ll] += alpha * x[ll];
        }
        pA += kna + bs * (sda - 1) + kna * bs;
        x += kna;
        kmax -= kna;
    }
    for (jj = 0; jj < kmax - 3; jj += 4) {
        pA[jj * sda + (jj + 0) * bs + 0] += alpha * x[jj + 0];
        pA[jj * sda + (jj + 1) * bs + 1] += alpha * x[jj + 1];
        pA[jj * sda + (jj + 2) * bs + 2] += alpha * x[jj + 2];
        pA[jj * sda + (jj + 3) * bs + 3] += alpha * x[jj + 3];
    }
    for (ll = 0; ll < kmax - jj; ll++) {
        pA[jj * sda + (jj + ll) * bs + ll] += alpha * x[jj + ll];
    }
}


// add scaled strvec to diagonal of strmat, sparse formulation
void ddiaad_sp(int kmax, double alpha, struct vec* sx, int xi, int* idx, struct mat* sD, int di, int dj) {

    // invalidate stored inverse diagonal
    sD->use_dA = 0;

    const int bs = D_PS;
    double* x = sx->pa + xi;
    int sdd = sD->cn;
    double* pD = sD->pA;
    int ii, jj;
    for (jj = 0; jj < kmax; jj++) {
        ii = idx[jj];
        pD[(ii + di) / bs * bs * sdd + (ii + di) % bs + (ii + dj) * bs] += alpha * x[jj];
    }
}


// add scaled strvec to another strvec and insert to diagonal of strmat, sparse formulation
void ddiaadin_sp(int kmax, double alpha, struct vec* sx, int xi, struct vec* sy, int yi, int* idx, struct mat* sD, int di, int dj) {

    // invalidate stored inverse diagonal
    sD->use_dA = 0;

    const int bs = D_PS;
    double* x = sx->pa + xi;
    double* y = sy->pa + yi;
    int sdd = sD->cn;
    double* pD = sD->pA;
    int ii, jj;
    for (jj = 0; jj < kmax; jj++) {
        ii = idx[jj];
        pD[(ii + di) / bs * bs * sdd + (ii + di) % bs + (ii + dj) * bs] = y[jj] + alpha * x[jj];
    }
}

// add scaled vector to row, sparse formulation
void drowad_libsp(int kmax, int* idx, double alpha, double* x, double* pD) {
    const int bs = D_PS;
    int ii, jj;

    for (jj = 0; jj < kmax; jj++) {
        ii = idx[jj];
        pD[ii * bs] += alpha * x[jj];
    }
}
// add scaled strvec to row of strmat, sparse formulation
void drowad_sp(int kmax, double alpha, struct vec* sx, int xi, int* idx, struct mat* sD, int di, int dj) {

    // invalidate stored inverse diagonal
    sD->use_dA = 0;

    const int bs = D_PS;
    double* x = sx->pa + xi;
    int sdd = sD->cn;
    double* pD = sD->pA + di / bs * bs * sdd + di % bs + dj * bs;
    drowad_libsp(kmax, idx, alpha, x, pD);
}


// add scaled strvec to strvec, sparse formulation
void dvecad_sp(int m, double alpha, struct vec* sx, int xi, int* idx, struct vec* sz, int zi) {
    double* x = sx->pa + xi;
    double* z = sz->pa + zi;
    int ii;
    for (ii = 0; ii < m; ii++)
        z[idx[ii]] += alpha * x[ii];
}


// insert scaled strvec to strvec, sparse formulation
void dvecin_sp(int m, double alpha, struct vec* sx, int xi, int* idx, struct vec* sz, int zi) {
    double* x = sx->pa + xi;
    double* z = sz->pa + zi;
    int ii;
    for (ii = 0; ii < m; ii++)
        z[idx[ii]] = alpha * x[ii];
}


// extract scaled strvec to strvec, sparse formulation
void dvecex_sp(int m, double alpha, int* idx, struct vec* sx, int xi, struct vec* sz, int zi) {
    double* x = sx->pa + xi;
    double* z = sz->pa + zi;
    int ii;
    for (ii = 0; ii < m; ii++)
        z[ii] = alpha * x[idx[ii]];
}


// z += alpha * x[idx]
void dvecexad_sp(int m, double alpha, int* idx, struct vec* sx, int xi, struct vec* sz, int zi) {
    double* x = sx->pa + xi;
    double* z = sz->pa + zi;
    int ii;
    for (ii = 0; ii < m; ii++)
        z[ii] += alpha * x[idx[ii]];
}


// clip strvec between two strvec
void dveccl(int m, struct vec* sxm, int xim, struct vec* sx, int xi, struct vec* sxp, int xip, struct vec* sz, int zi) {

    double* xm = sxm->pa + xim;
    double* x = sx->pa + xi;
    double* xp = sxp->pa + xip;
    double* z = sz->pa + zi;

    int ii;

    for (ii = 0; ii < m; ii++) {
        if (x[ii] >= xp[ii]) {
            z[ii] = xp[ii];
        } else if (x[ii] <= xm[ii]) {
            z[ii] = xm[ii];
        } else {
            z[ii] = x[ii];
        }
    }
}


// clip strvec between two strvec, with mask
void dveccl_mask(int m, struct vec* sxm, int xim, struct vec* sx, int xi, struct vec* sxp, int xip, struct vec* sz, int zi, struct vec* sm, int mi) {

    double* xm = sxm->pa + xim;
    double* x = sx->pa + xi;
    double* xp = sxp->pa + xip;
    double* z = sz->pa + zi;
    double* mask = sm->pa + mi;

    int ii;

    for (ii = 0; ii < m; ii++) {
        if (x[ii] >= xp[ii]) {
            z[ii] = xp[ii];
            mask[ii] = 1.0;
        } else if (x[ii] <= xm[ii]) {
            z[ii] = xm[ii];
            mask[ii] = -1.0;
        } else {
            z[ii] = x[ii];
            mask[ii] = 0.0;
        }
    }
}


// zero out strvec to strvec with mask
void dvecze(int m, struct vec* sm, int mi, struct vec* sv, int vi, struct vec* se, int ei) {
    double* mask = sm->pa + mi;
    double* v = sv->pa + vi;
    double* e = se->pa + ei;

    int ii;

    for (ii = 0; ii < m; ii++) {
        if (mask[ii] == 0) {
            e[ii] = v[ii];
        } else {
            e[ii] = 0;
        }
    }
}


// compute inf norm of vector
void dvecnrm_inf(int m, struct vec* sx, int xi, double* ptr_norm) {
    int ii;
    double* x = sx->pa + xi;
    double norm = 0.0;
    int is_nan = 0;
    double tmp;
    for (ii = 0; ii < m; ii++) {
        norm = fmax(norm, fabs(x[ii]));
        is_nan |= x[ii] != x[ii];
    }
    *ptr_norm = is_nan == 0 ? norm : NAN;
}


// compute 2 norm of vector
void dvecnrm_2(int m, struct vec* sx, int xi, double* ptr_norm) {
    int ii;
    double* x = sx->pa + xi;
    double norm = 0.0;
    for (ii = 0; ii < m; ii++) {
        norm += x[ii] * x[ii];
    }
    norm = sqrt(norm);
    *ptr_norm = norm;
}


// permute elements of a vector struct
void dvecpe(int kmax, int* ipiv, struct vec* sx, int xi) {
    int ii;
    double tmp;
    double* x = sx->pa + xi;
    for (ii = 0; ii < kmax; ii++) {
        if (ipiv[ii] != ii) {
            tmp = x[ipiv[ii]];
            x[ipiv[ii]] = x[ii];
            x[ii] = tmp;
        }
    }
}


// inverse permute elements of a vector struct
void dvecpei(int kmax, int* ipiv, struct vec* sx, int xi) {
    int ii;
    double tmp;
    double* x = sx->pa + xi;
    for (ii = kmax - 1; ii >= 0; ii--) {
        if (ipiv[ii] != ii) {
            tmp = x[ipiv[ii]];
            x[ipiv[ii]] = x[ii];
            x[ii] = tmp;
        }
    }
}
