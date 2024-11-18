#include "tinyhpipm/blas/kernel.h"
#include "tinyhpipm/blas/misc.h"
#include "tinyhpipm/blas/struct.h"
#include <stdio.h>
#include <stdlib.h>

void dpotrf_l(int m, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj) {
    if (m <= 0) {
        return;
    }

    if (ci != 0 | di != 0) {
        fprintf(stderr, "\n\tdpotrf_l: feature not implemented yet: ci=%d, di=%d\n", ci, di);
        exit(1);
    }

    const int ps = D_PS;

    double alpha = 1.0;

    int sdc = sC->cn;
    int sdd = sD->cn;
    double* pC = sC->pA + cj * ps;
    double* pD = sD->pA + dj * ps;
    double* dD = sD->dA;

    if (di == 0 & dj == 0)  // XXX what to do if di and dj are not zero
        sD->use_dA = m;
    else
        sD->use_dA = 0;

    int i, j, l;

    i = 0;
    for (; i < m - 3; i += 4) {
        j = 0;
        for (; j < i; j += 4) {
            kernel_dtrsm_nt_rl_inv_4x4_lib4(j, &pD[i * sdd], &pD[j * sdd], &alpha, &pC[j * ps + i * sdc], &pD[j * ps + i * sdd], &pD[j * ps + j * sdd], &dD[j]);
        }
        kernel_dpotrf_nt_l_4x4_lib4(j, &pD[i * sdd], &pD[j * sdd], &pC[j * ps + j * sdc], &pD[j * ps + j * sdd], &dD[j]);
    }
    if (m > i) {
        // clean up loops definitions
        j = 0;
        if (m - i == 4) {
            for (; j < i; j += 4) {
                kernel_dtrsm_nt_rl_inv_4x4_lib4(j, &pD[i * sdd], &pD[j * sdd], &alpha, &pC[j * ps + i * sdc], &pD[j * ps + i * sdd], &pD[j * ps + j * sdd], &dD[j]);
            }
            kernel_dpotrf_nt_l_4x4_lib4(j, &pD[i * sdd], &pD[j * sdd], &pC[j * ps + j * sdc], &pD[j * ps + j * sdd], &dD[j]);
        } else {
            for (; j < i; j += 4) {
                kernel_dtrsm_nt_rl_inv_4x4_vs_lib4(j, &pD[i * sdd], &pD[j * sdd], &alpha, &pC[j * ps + i * sdc], &pD[j * ps + i * sdd], &pD[j * ps + j * sdd], &dD[j], m - i, m - j);
            }
            kernel_dpotrf_nt_l_4x4_vs_lib4(j, &pD[i * sdd], &pD[j * sdd], &pC[j * ps + j * sdc], &pD[j * ps + j * sdd], &dD[j], m - i, m - j);
        }
    }
}


// dpotrf
void dpotrf_l_mn(int m, int n, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj) {

    if (m <= 0 || n <= 0)
        return;

    if (ci != 0 | di != 0) {
        fprintf(stderr, "\n\tdpotrf_l_mn: feature not implemented yet: ci=%d, di=%d\n", ci, di);
        exit(1);
    }

    const int ps = D_PS;

    double alpha = 1.0;

    int sdc = sC->cn;
    int sdd = sD->cn;
    double* pC = sC->pA + cj * ps;
    double* pD = sD->pA + dj * ps;
    double* dD = sD->dA;

    if (di == 0 & dj == 0)  // XXX what to do if di and dj are not zero
        sD->use_dA = 1;
    else
        sD->use_dA = 0;

    int i, j, l;

    i = 0;
    for (; i < m - 3; i += 4) {
        j = 0;
        for (; j < i & j < n - 3; j += 4) {
            kernel_dtrsm_nt_rl_inv_4x4_lib4(j, &pD[i * sdd], &pD[j * sdd], &alpha, &pC[j * ps + i * sdc], &pD[j * ps + i * sdd], &pD[j * ps + j * sdd], &dD[j]);
        }
        if (j < n) {
            if (j < i)  // dtrsm
            {
                kernel_dtrsm_nt_rl_inv_4x4_vs_lib4(j, &pD[i * sdd], &pD[j * sdd], &alpha, &pC[j * ps + i * sdc], &pD[j * ps + i * sdd], &pD[j * ps + j * sdd], &dD[j], m - i, n - j);
            } else  // dpotrf
            {
                if (j < n - 3) {
                    kernel_dpotrf_nt_l_4x4_lib4(j, &pD[i * sdd], &pD[j * sdd], &pC[j * ps + j * sdc], &pD[j * ps + j * sdd], &dD[j]);
                } else {
                    kernel_dpotrf_nt_l_4x4_vs_lib4(j, &pD[i * sdd], &pD[j * sdd], &pC[j * ps + j * sdc], &pD[j * ps + j * sdd], &dD[j], m - i, n - j);
                }
            }
        }
    }
    if (m > i) {
        // clean up loops definitions
        j = 0;
        for (; j < i & j < n; j += 4) {
            kernel_dtrsm_nt_rl_inv_4x4_vs_lib4(j, &pD[i * sdd], &pD[j * sdd], &alpha, &pC[j * ps + i * sdc], &pD[j * ps + i * sdd], &pD[j * ps + j * sdd], &dD[j], m - i, n - j);
        }
        if (j < n) {
            kernel_dpotrf_nt_l_4x4_vs_lib4(j, &pD[i * sdd], &pD[j * sdd], &pC[j * ps + j * sdc], &pD[j * ps + j * sdd], &dD[j], m - i, n - j);
        }
    }
}


void dsyrk_dpotrf_ln(int m, int k, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj) {
    //	dsyrk_dpotrf_ln_mn(m, m, k, sA, ai, aj, sB, bi, bj, sC, ci, cj, sD, di, dj);
    //	return;

    if (m <= 0)
        return;

    if (ai != 0 | bi != 0 | ci != 0 | di != 0) {
        fprintf(stderr, "\ndsyrk_dpotrf_ln: feature not implemented yet: ai=%d, bi=%d, ci=%d, di=%d\n", ai, bi, ci, di);
        exit(1);
    }

    const int ps = 4;

    int sda = sA->cn;
    int sdb = sB->cn;
    int sdc = sC->cn;
    int sdd = sD->cn;
    double* pA = sA->pA + aj * ps;
    double* pB = sB->pA + bj * ps;
    double* pC = sC->pA + cj * ps;
    double* pD = sD->pA + dj * ps;
    double* dD = sD->dA;  // XXX what to do if di and dj are not zero

    if (di == 0 & dj == 0)
        sD->use_dA = 1;
    else
        sD->use_dA = 0;

    int i, j, l;

    i = 0;

    for (; i < m - 3; i += 4) {
        j = 0;
        for (; j < i; j += 4) {
            kernel_dgemm_dtrsm_nt_rl_inv_4x4_lib4(k, &pA[i * sda], &pB[j * sdb], j, &pD[i * sdd], &pD[j * sdd], &pC[j * ps + i * sdc], &pD[j * ps + i * sdd], &pD[j * ps + j * sdd], &dD[j]);
        }
        kernel_dsyrk_dpotrf_nt_l_4x4_lib4(k, &pA[i * sda], &pB[j * sdb], j, &pD[i * sdd], &pD[j * sdd], &pC[j * ps + j * sdc], &pD[j * ps + j * sdd], &dD[j]);
    }
    if (m > i) {
        j = 0;
        if (m - i == 4) {
            for (; j < i; j += 4) {
                kernel_dgemm_dtrsm_nt_rl_inv_4x4_lib4(k, &pA[i * sda], &pB[j * sdb], j, &pD[i * sdd], &pD[j * sdd], &pC[j * ps + i * sdc], &pD[j * ps + i * sdd], &pD[j * ps + j * sdd], &dD[j]);
            }
            kernel_dsyrk_dpotrf_nt_l_4x4_lib4(k, &pA[i * sda], &pB[j * sdb], j, &pD[i * sdd], &pD[j * sdd], &pC[j * ps + j * sdc], &pD[j * ps + j * sdd], &dD[j]);
        } else {
            for (; j < i; j += 4) {
                kernel_dgemm_dtrsm_nt_rl_inv_4x4_vs_lib4(k, &pA[i * sda], &pB[j * sdb], j, &pD[i * sdd], &pD[j * sdd], &pC[j * ps + i * sdc], &pD[j * ps + i * sdd], &pD[j * ps + j * sdd], &dD[j], m - i, m - j);
            }
            kernel_dsyrk_dpotrf_nt_l_4x4_vs_lib4(k, &pA[i * sda], &pB[j * sdb], j, &pD[i * sdd], &pD[j * sdd], &pC[j * ps + j * sdc], &pD[j * ps + j * sdd], &dD[j], m - i, m - j);
        }
    }
}

// dsyrk dpotrf
void dsyrk_dpotrf_ln_mn(int m, int n, int k, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj) {

    if (m <= 0 || n <= 0) {
        return;
    }

    if (ai != 0 | bi != 0 | ci != 0 | di != 0) {
        fprintf(stderr, "\ndsyrk_dpotrf_ln_mn: feature not implemented yet: ai=%d, bi=%d, ci=%d, di=%d\n", ai, bi, ci, di);
        exit(1);
    }

    const int ps = 4;

    int sda = sA->cn;
    int sdb = sB->cn;
    int sdc = sC->cn;
    int sdd = sD->cn;
    double* pA = sA->pA + aj * ps;
    double* pB = sB->pA + bj * ps;
    double* pC = sC->pA + cj * ps;
    double* pD = sD->pA + dj * ps;
    double* dD = sD->dA;  // XXX what to do if di and dj are not zero

    if (di == 0 & dj == 0)
        sD->use_dA = 1;
    else
        sD->use_dA = 0;

    int i, j, l;

    i = 0;

    for (; i < m - 3; i += 4) {
        j = 0;
        for (; j < i & j < n - 3; j += 4) {
            kernel_dgemm_dtrsm_nt_rl_inv_4x4_lib4(k, &pA[i * sda], &pB[j * sdb], j, &pD[i * sdd], &pD[j * sdd], &pC[j * ps + i * sdc], &pD[j * ps + i * sdd], &pD[j * ps + j * sdd], &dD[j]);
        }
        if (j < n) {
            if (j < i)  // dgemm
            {
                kernel_dgemm_dtrsm_nt_rl_inv_4x4_vs_lib4(k, &pA[i * sda], &pB[j * sdb], j, &pD[i * sdd], &pD[j * sdd], &pC[j * ps + i * sdc], &pD[j * ps + i * sdd], &pD[j * ps + j * sdd], &dD[j], m - i, n - j);
            } else  // dsyrk
            {
                if (j < n - 3) {
                    kernel_dsyrk_dpotrf_nt_l_4x4_lib4(k, &pA[i * sda], &pB[j * sdb], j, &pD[i * sdd], &pD[j * sdd], &pC[j * ps + j * sdc], &pD[j * ps + j * sdd], &dD[j]);
                } else {
                    kernel_dsyrk_dpotrf_nt_l_4x4_vs_lib4(k, &pA[i * sda], &pB[j * sdb], j, &pD[i * sdd], &pD[j * sdd], &pC[j * ps + j * sdc], &pD[j * ps + j * sdd], &dD[j], m - i, n - j);
                }
            }
        }
    }
    if (m > i) {
        j = 0;
        for (; j < i & j < n; j += 4) {
            kernel_dgemm_dtrsm_nt_rl_inv_4x4_vs_lib4(k, &pA[i * sda], &pB[j * sdb], j, &pD[i * sdd], &pD[j * sdd], &pC[j * ps + i * sdc], &pD[j * ps + i * sdd], &pD[j * ps + j * sdd], &dD[j], m - i, n - j);
        }
        if (j < n) {
            kernel_dsyrk_dpotrf_nt_l_4x4_vs_lib4(k, &pA[i * sda], &pB[j * sdb], j, &pD[i * sdd], &pD[j * sdd], &pC[j * ps + j * sdc], &pD[j * ps + j * sdd], &dD[j], m - i, n - j);
        }
    }
}

int dorglq_worksize(int m, int n, int k) {
    return 0;
}


void dorglq(int m, int n, int k, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj, void* work) {
    if (m <= 0 | n <= 0)
        return;

    // TODO check that k <= m <= n

    if (di != 0) {
        fprintf(stderr, "\ndorglq: feature not implemented yet: di=%d\n", di);
        exit(1);
    }

    // invalidate stored inverse diagonal of result matrix
    sD->use_dA = 0;

    const int ps = 4;

    // extract dimensions
    int sdc = sC->cn;
    int sdd = sD->cn;

    // go to submatrix
    double* pC = &(MATEL(sC, ci, cj));
    double* dC = sC->dA + ci;
    double* pD = &(MATEL(sD, di, dj));

    double pT[144] = {0};  // XXX smaller ?
    double pK[96] = {0};  // XXX smaller ?

    int ii, jj, ll, idx;

    // set result matrix to the identity
    dgese(m, n, 0.0, sD, di, dj);
    ddiare(m, 1.0, sD, di, dj);

    int kr4 = k % 4;
    int km4 = k - kr4;

    // clear out the end
    if (kr4 > 0) {
        if (kr4 == 1) {
            kernel_dlarft_1_lib4(n - km4, pC + km4 * sdc + km4 * ps, dC + km4, pT);
            for (jj = 0; jj < m - km4 - 3; jj += 4) {
                kernel_dlarfb1_rt_4_lib4(n - km4, pC + km4 * sdc + km4 * ps, pT, pD + (km4 + jj) * sdd + km4 * ps);
            }
            for (ll = 0; ll < m - km4 - jj; ll++) {
                kernel_dlarfb1_rt_1_lib4(n - km4, pC + km4 * sdc + km4 * ps, pT, pD + (km4 + jj) * sdd + km4 * ps + ll);
            }
        } else if (kr4 == 2) {
            kernel_dlarft_2_lib4(n - km4, pC + km4 * sdc + km4 * ps, dC + km4, pT);
            for (jj = 0; jj < m - km4 - 3; jj += 4) {
                kernel_dlarfb2_rt_4_lib4(n - km4, pC + km4 * sdc + km4 * ps, pT, pD + (km4 + jj) * sdd + km4 * ps);
            }
            for (ll = 0; ll < m - km4 - jj; ll++) {
                kernel_dlarfb2_rt_1_lib4(n - km4, pC + km4 * sdc + km4 * ps, pT, pD + (km4 + jj) * sdd + km4 * ps + ll);
            }
        } else  // kr4==3
        {
            kernel_dlarft_3_lib4(n - km4, pC + km4 * sdc + km4 * ps, dC + km4, pT);
            for (jj = 0; jj < m - km4 - 3; jj += 4) {
                kernel_dlarfb3_rt_4_lib4(n - km4, pC + km4 * sdc + km4 * ps, pT, pD + (km4 + jj) * sdd + km4 * ps);
            }
            for (ll = 0; ll < m - km4 - jj; ll++) {
                kernel_dlarfb3_rt_1_lib4(n - km4, pC + km4 * sdc + km4 * ps, pT, pD + (km4 + jj) * sdd + km4 * ps + ll);
            }
        }
    }
    // main loop
    for (ii = 0; ii < km4; ii += 4) {
        idx = km4 - ii - 4;
        kernel_dlarft_4_lib4(n - idx, pC + idx * sdc + idx * ps, dC + idx, pT);
        for (jj = 0; jj < m - idx - 3; jj += 4) {
            kernel_dlarfb4_rt_4_lib4(n - idx, pC + idx * sdc + idx * ps, pT, pD + (idx + jj) * sdd + idx * ps);
        }
        for (ll = 0; ll < m - idx - jj; ll++) {
            kernel_dlarfb4_rt_1_lib4(n - idx, pC + idx * sdc + idx * ps, pT, pD + (idx + jj) * sdd + idx * ps + ll);
        }
    }
}


int dgelqf_worksize(int m, int n) {
    return 0;
}


void dgelqf(int m, int n, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj, void* work) {
    if (m <= 0 | n <= 0)
        return;

    // invalidate stored inverse diagonal of result matrix
    sD->use_dA = 0;

    const int ps = 4;

    // extract dimensions
    int sdc = sC->cn;
    int sdd = sD->cn;

    // go to submatrix
    double* pC = &(MATEL(sC, ci, cj));
    double* pD = &(MATEL(sD, di, dj));

    double* dD = sD->dA + di;
    double pT[144] = {0};  // XXX smaller ?
    double pK[96] = {0};  // XXX smaller ?
    /* if(pC!=pD) */
    /* dgecp_lib(m, n, 1.0, ci&(ps-1), pC, sdc, di&(ps-1), pD, sdd); */
    /* // where ci&(ps-1) == ci%ps */

    if (pC != pD)
        // copy strmat submatrix
        dgecp(m, n, sC, ci, cj, sD, di, dj);

    int ii, jj, ll;
    int imax0 = (ps - (di & (ps - 1))) & (ps - 1);
    int imax = m < n ? m : n;
    imax0 = imax < imax0 ? imax : imax0;
    if (imax0 > 0) {
        kernel_dgelqf_vs_lib4(m, n, imax0, di & (ps - 1), pD, sdd, dD);
        pD += imax0 - ps + ps * sdd + imax0 * ps;
        dD += imax0;
        m -= imax0;
        n -= imax0;
        imax -= imax0;
    }
    ii = 0;
    for (ii = 0; ii < imax - 4; ii += 4) {
        //		kernel_dgelqf_vs_lib4(4, n-ii, 4, 0, pD+ii*sdd+ii*ps, sdd, dD+ii);
        //		kernel_dgelqf_4_lib4(n-ii, pD+ii*sdd+ii*ps, dD+ii);
        //		kernel_dlarft_4_lib4(n-ii, pD+ii*sdd+ii*ps, dD+ii, pT);
        kernel_dgelqf_dlarft4_4_lib4(n - ii, pD + ii * sdd + ii * ps, dD + ii, pT);
        jj = ii + 4;
        for (; jj < m - 3; jj += 4) {
            kernel_dlarfb4_rn_4_lib4(n - ii, pD + ii * sdd + ii * ps, pT, pD + jj * sdd + ii * ps);
        }
        for (ll = 0; ll < m - jj; ll++) {
            kernel_dlarfb4_rn_1_lib4(n - ii, pD + ii * sdd + ii * ps, pT, pD + ll + jj * sdd + ii * ps);
        }
    }
    if (ii < imax) {
        if (ii == imax - 4) {
            kernel_dgelqf_4_lib4(n - ii, pD + ii * sdd + ii * ps, dD + ii);
        } else {
            kernel_dgelqf_vs_lib4(m - ii, n - ii, imax - ii, ii & (ps - 1), pD + ii * sdd + ii * ps, sdd, dD + ii);
        }
    }
}

void dgetrf_np(int m, int n, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj) {

    if (ci != 0 | di != 0) {
        fprintf(stderr, "\n\tdgetf_np: feature not implemented yet: ci=%d, di=%d\n", ci, di);
        exit(1);
    }

    const int ps = 4;

    int sdc = sC->cn;
    int sdd = sD->cn;
    double* pC = sC->pA + cj * ps;
    double* pD = sD->pA + dj * ps;
    double* dD = sD->dA;  // XXX what to do if di and dj are not zero

    if (di == 0 && dj == 0)
        sD->use_dA = 1;
    else
        sD->use_dA = 0;

    if (m <= 0 | n <= 0)
        return;

    double d1 = 1.0;

    int ii, jj, ie;

    // main loop
    ii = 0;
    for (; ii < m - 3; ii += 4) {
        jj = 0;
        // solve lower
        ie = n < ii ? n : ii;  // ie is multiple of 4
        for (; jj < ie - 3; jj += 4) {
            kernel_dtrsm_nn_ru_inv_4x4_lib4(jj, &pD[ii * sdd], &pD[jj * ps], sdd, &d1, &pC[jj * ps + ii * sdc], &pD[jj * ps + ii * sdd], &pD[jj * ps + jj * sdd], &dD[jj]);
        }
        if (jj < ie) {
            kernel_dtrsm_nn_ru_inv_4x4_vs_lib4(jj, &pD[ii * sdd], &pD[jj * ps], sdd, &d1, &pC[jj * ps + ii * sdc], &pD[jj * ps + ii * sdd], &pD[jj * ps + jj * sdd], &dD[jj], m - ii, ie - jj);
            jj += 4;
        }
        // factorize
        if (jj < n - 3) {
            kernel_dgetrf_nn_4x4_lib4(jj, &pD[ii * sdd], &pD[jj * ps], sdd, &pC[jj * ps + ii * sdc], &pD[jj * ps + ii * sdd], &dD[jj]);
            jj += 4;
        } else if (jj < n) {
            kernel_dgetrf_nn_4x4_vs_lib4(jj, &pD[ii * sdd], &pD[jj * ps], sdd, &pC[jj * ps + ii * sdc], &pD[jj * ps + ii * sdd], &dD[jj], m - ii, n - jj);
            jj += 4;
        }
        // solve upper
        for (; jj < n - 3; jj += 4) {
            kernel_dtrsm_nn_ll_one_4x4_lib4(ii, &pD[ii * sdd], &pD[jj * ps], sdd, &d1, &pC[jj * ps + ii * sdc], &pD[jj * ps + ii * sdd], &pD[ii * ps + ii * sdd]);
        }
        if (jj < n) {
            kernel_dtrsm_nn_ll_one_4x4_vs_lib4(ii, &pD[ii * sdd], &pD[jj * ps], sdd, &d1, &pC[jj * ps + ii * sdc], &pD[jj * ps + ii * sdd], &pD[ii * ps + ii * sdd], m - ii, n - jj);
        }
    }
    jj = 0;
    // solve lower
    ie = n < ii ? n : ii;  // ie is multiple of 4
    for (; jj < ie; jj += 4) {
        kernel_dtrsm_nn_ru_inv_4x4_vs_lib4(jj, &pD[ii * sdd], &pD[jj * ps], sdd, &d1, &pC[jj * ps + ii * sdc], &pD[jj * ps + ii * sdd], &pD[jj * ps + jj * sdd], &dD[jj], m - ii, ie - jj);
    }
    // factorize
    if (jj < n) {
        kernel_dgetrf_nn_4x4_vs_lib4(jj, &pD[ii * sdd], &pD[jj * ps], sdd, &pC[jj * ps + ii * sdc], &pD[jj * ps + ii * sdd], &dD[jj], m - ii, n - jj);
        jj += 4;
    }
    // solve upper
    for (; jj < n; jj += 4) {
        kernel_dtrsm_nn_ll_one_4x4_vs_lib4(ii, &pD[ii * sdd], &pD[jj * ps], sdd, &d1, &pC[jj * ps + ii * sdc], &pD[jj * ps + ii * sdd], &pD[ii * ps + ii * sdd], m - ii, n - jj);
    }
}

void dgetrf_rp(int m, int n, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj, int* ipiv) {
    if (ci != 0 | di != 0) {
        printf("\ndgetrf_rp: feature not implemented yet: ci=%d, di=%d\n", ci, di);
        exit(1);
    }

    const int ps = D_PS;

    int sdc = sC->cn;
    int sdd = sD->cn;
    double* pC = sC->pA + cj * ps;
    double* pD = sD->pA + dj * ps;
    double* dD = sD->dA;  // XXX what to do if di and dj are not zero

    if (di == 0 && dj == 0)
        sD->use_dA = 1;
    else
        sD->use_dA = 0;

    if (m <= 0 | n <= 0)
        return;

    int ii, jj, i0, i1, j0, ll, p;

    double d1 = 1.0;
    double dm1 = -1.0;

    // needs to perform row-excanges on the yet-to-be-factorized matrix too
    if (pC != pD)
        dgecp(m, n, sC, ci, cj, sD, di, dj);

    // minimum matrix size
    p = n < m ? n : m;  // XXX

    // main loop
    // 4 columns at a time
    jj = 0;
    for (; jj < p - 3; jj += 4)  // XXX
    {
        // pivot & factorize & solve lower
        ii = jj;
        i0 = ii;
        for (; ii < m - 3; ii += 4) {
            kernel_dgemm_nn_4x4_lib4(jj, &dm1, &pD[ii * sdd], 0, &pD[jj * ps], sdd, &d1, &pD[jj * ps + ii * sdd], &pD[jj * ps + ii * sdd]);
        }
        if (m - ii > 0) {
            kernel_dgemm_nn_4x4_vs_lib4(jj, &dm1, &pD[ii * sdd], 0, &pD[jj * ps], sdd, &d1, &pD[jj * ps + ii * sdd], &pD[jj * ps + ii * sdd], m - ii, 4);
        }
        kernel_dgetrf_pivot_4_lib4(m - i0, &pD[jj * ps + i0 * sdd], sdd, &dD[jj], &ipiv[i0]);
        ipiv[i0 + 0] += i0;
        if (ipiv[i0 + 0] != i0 + 0) {
            kernel_drowsw_lib4(jj, pD + (i0 + 0) / ps * ps * sdd + (i0 + 0) % ps, pD + (ipiv[i0 + 0]) / ps * ps * sdd + (ipiv[i0 + 0]) % ps);
            kernel_drowsw_lib4(n - jj - 4, pD + (i0 + 0) / ps * ps * sdd + (i0 + 0) % ps + (jj + 4) * ps, pD + (ipiv[i0 + 0]) / ps * ps * sdd + (ipiv[i0 + 0]) % ps + (jj + 4) * ps);
        }
        ipiv[i0 + 1] += i0;
        if (ipiv[i0 + 1] != i0 + 1) {
            kernel_drowsw_lib4(jj, pD + (i0 + 1) / ps * ps * sdd + (i0 + 1) % ps, pD + (ipiv[i0 + 1]) / ps * ps * sdd + (ipiv[i0 + 1]) % ps);
            kernel_drowsw_lib4(n - jj - 4, pD + (i0 + 1) / ps * ps * sdd + (i0 + 1) % ps + (jj + 4) * ps, pD + (ipiv[i0 + 1]) / ps * ps * sdd + (ipiv[i0 + 1]) % ps + (jj + 4) * ps);
        }
        ipiv[i0 + 2] += i0;
        if (ipiv[i0 + 2] != i0 + 2) {
            kernel_drowsw_lib4(jj, pD + (i0 + 2) / ps * ps * sdd + (i0 + 2) % ps, pD + (ipiv[i0 + 2]) / ps * ps * sdd + (ipiv[i0 + 2]) % ps);
            kernel_drowsw_lib4(n - jj - 4, pD + (i0 + 2) / ps * ps * sdd + (i0 + 2) % ps + (jj + 4) * ps, pD + (ipiv[i0 + 2]) / ps * ps * sdd + (ipiv[i0 + 2]) % ps + (jj + 4) * ps);
        }
        ipiv[i0 + 3] += i0;
        if (ipiv[i0 + 3] != i0 + 3) {
            kernel_drowsw_lib4(jj, pD + (i0 + 3) / ps * ps * sdd + (i0 + 3) % ps, pD + (ipiv[i0 + 3]) / ps * ps * sdd + (ipiv[i0 + 3]) % ps);
            kernel_drowsw_lib4(n - jj - 4, pD + (i0 + 3) / ps * ps * sdd + (i0 + 3) % ps + (jj + 4) * ps, pD + (ipiv[i0 + 3]) / ps * ps * sdd + (ipiv[i0 + 3]) % ps + (jj + 4) * ps);
        }

        // solve upper
        ll = jj + 4;
        for (; ll < n - 3; ll += 4) {
            kernel_dtrsm_nn_ll_one_4x4_lib4(i0, &pD[i0 * sdd], &pD[ll * ps], sdd, &d1, &pD[ll * ps + i0 * sdd], &pD[ll * ps + i0 * sdd], &pD[i0 * ps + i0 * sdd]);
        }
        if (n - ll > 0) {
            kernel_dtrsm_nn_ll_one_4x4_vs_lib4(i0, &pD[i0 * sdd], &pD[ll * ps], sdd, &d1, &pD[ll * ps + i0 * sdd], &pD[ll * ps + i0 * sdd], &pD[i0 * ps + i0 * sdd], 4, n - ll);
        }
    }
    if (m >= n) {
        if (n - jj > 0) {
            // 1-4 columns at a time
            // pivot & factorize & solve lower
            ii = jj;
            i0 = ii;
            for (; ii < m; ii += 4) {
                kernel_dgemm_nn_4x4_vs_lib4(jj, &dm1, &pD[ii * sdd], 0, &pD[jj * ps], sdd, &d1, &pD[jj * ps + ii * sdd], &pD[jj * ps + ii * sdd], m - ii, n - jj);
            }
            kernel_dgetrf_pivot_4_vs_lib4(m - i0, &pD[jj * ps + i0 * sdd], sdd, &dD[jj], &ipiv[i0], n - jj);
            ipiv[i0 + 0] += i0;
            if (ipiv[i0 + 0] != i0 + 0) {
                kernel_drowsw_lib4(jj, pD + (i0 + 0) / ps * ps * sdd + (i0 + 0) % ps, pD + (ipiv[i0 + 0]) / ps * ps * sdd + (ipiv[i0 + 0]) % ps);
                kernel_drowsw_lib4(n - jj - 4, pD + (i0 + 0) / ps * ps * sdd + (i0 + 0) % ps + (jj + 4) * ps, pD + (ipiv[i0 + 0]) / ps * ps * sdd + (ipiv[i0 + 0]) % ps + (jj + 4) * ps);
            }
            if (n - jj > 1) {
                ipiv[i0 + 1] += i0;
                if (ipiv[i0 + 1] != i0 + 1) {
                    kernel_drowsw_lib4(jj, pD + (i0 + 1) / ps * ps * sdd + (i0 + 1) % ps, pD + (ipiv[i0 + 1]) / ps * ps * sdd + (ipiv[i0 + 1]) % ps);
                    kernel_drowsw_lib4(n - jj - 4, pD + (i0 + 1) / ps * ps * sdd + (i0 + 1) % ps + (jj + 4) * ps, pD + (ipiv[i0 + 1]) / ps * ps * sdd + (ipiv[i0 + 1]) % ps + (jj + 4) * ps);
                }
                if (n - jj > 2) {
                    ipiv[i0 + 2] += i0;
                    if (ipiv[i0 + 2] != i0 + 2) {
                        kernel_drowsw_lib4(jj, pD + (i0 + 2) / ps * ps * sdd + (i0 + 2) % ps, pD + (ipiv[i0 + 2]) / ps * ps * sdd + (ipiv[i0 + 2]) % ps);
                        kernel_drowsw_lib4(n - jj - 4, pD + (i0 + 2) / ps * ps * sdd + (i0 + 2) % ps + (jj + 4) * ps, pD + (ipiv[i0 + 2]) / ps * ps * sdd + (ipiv[i0 + 2]) % ps + (jj + 4) * ps);
                    }
                    if (n - jj > 3) {
                        ipiv[i0 + 3] += i0;
                        if (ipiv[i0 + 3] != i0 + 3) {
                            kernel_drowsw_lib4(jj, pD + (i0 + 3) / ps * ps * sdd + (i0 + 3) % ps, pD + (ipiv[i0 + 3]) / ps * ps * sdd + (ipiv[i0 + 3]) % ps);
                            kernel_drowsw_lib4(n - jj - 4, pD + (i0 + 3) / ps * ps * sdd + (i0 + 3) % ps + (jj + 4) * ps, pD + (ipiv[i0 + 3]) / ps * ps * sdd + (ipiv[i0 + 3]) % ps + (jj + 4) * ps);
                        }
                    }
                }
            }

            // solve upper
            if (0)  // there is no upper
            {
                ll = jj + 4;
                for (; ll < n; ll += 4) {
                    kernel_dtrsm_nn_ll_one_4x4_vs_lib4(i0, &pD[i0 * sdd], &pD[ll * ps], sdd, &d1, &pD[ll * ps + i0 * sdd], &pD[ll * ps + i0 * sdd], &pD[i0 * ps + i0 * sdd], m - i0, n - ll);
                }
            }
        }
    } else {
        if (m - jj > 0) {
            // 1-4 rows at a time
            // pivot & factorize & solve lower
            ii = jj;
            i0 = ii;
            kernel_dgemm_nn_4x4_vs_lib4(jj, &dm1, &pD[ii * sdd], 0, &pD[jj * ps], sdd, &d1, &pD[jj * ps + ii * sdd], &pD[jj * ps + ii * sdd], m - ii, n - jj);
            kernel_dgetrf_pivot_4_vs_lib4(m - i0, &pD[jj * ps + i0 * sdd], sdd, &dD[jj], &ipiv[i0], n - jj);
            ipiv[i0 + 0] += i0;
            if (ipiv[i0 + 0] != i0 + 0) {
                kernel_drowsw_lib4(jj, pD + (i0 + 0) / ps * ps * sdd + (i0 + 0) % ps, pD + (ipiv[i0 + 0]) / ps * ps * sdd + (ipiv[i0 + 0]) % ps);
                kernel_drowsw_lib4(n - jj - 4, pD + (i0 + 0) / ps * ps * sdd + (i0 + 0) % ps + (jj + 4) * ps, pD + (ipiv[i0 + 0]) / ps * ps * sdd + (ipiv[i0 + 0]) % ps + (jj + 4) * ps);
            }
            if (m - i0 > 1) {
                ipiv[i0 + 1] += i0;
                if (ipiv[i0 + 1] != i0 + 1) {
                    kernel_drowsw_lib4(jj, pD + (i0 + 1) / ps * ps * sdd + (i0 + 1) % ps, pD + (ipiv[i0 + 1]) / ps * ps * sdd + (ipiv[i0 + 1]) % ps);
                    kernel_drowsw_lib4(n - jj - 4, pD + (i0 + 1) / ps * ps * sdd + (i0 + 1) % ps + (jj + 4) * ps, pD + (ipiv[i0 + 1]) / ps * ps * sdd + (ipiv[i0 + 1]) % ps + (jj + 4) * ps);
                }
                if (m - i0 > 2) {
                    ipiv[i0 + 2] += i0;
                    if (ipiv[i0 + 2] != i0 + 2) {
                        kernel_drowsw_lib4(jj, pD + (i0 + 2) / ps * ps * sdd + (i0 + 2) % ps, pD + (ipiv[i0 + 2]) / ps * ps * sdd + (ipiv[i0 + 2]) % ps);
                        kernel_drowsw_lib4(n - jj - 4, pD + (i0 + 2) / ps * ps * sdd + (i0 + 2) % ps + (jj + 4) * ps, pD + (ipiv[i0 + 2]) / ps * ps * sdd + (ipiv[i0 + 2]) % ps + (jj + 4) * ps);
                    }
                    if (m - i0 > 3) {
                        ipiv[i0 + 3] += i0;
                        if (ipiv[i0 + 3] != i0 + 3) {
                            kernel_drowsw_lib4(jj, pD + (i0 + 3) / ps * ps * sdd + (i0 + 3) % ps, pD + (ipiv[i0 + 3]) / ps * ps * sdd + (ipiv[i0 + 3]) % ps);
                            kernel_drowsw_lib4(n - jj - 4, pD + (i0 + 3) / ps * ps * sdd + (i0 + 3) % ps + (jj + 4) * ps, pD + (ipiv[i0 + 3]) / ps * ps * sdd + (ipiv[i0 + 3]) % ps + (jj + 4) * ps);
                        }
                    }
                }
            }

            // solve upper
            ll = jj + 4;
            for (; ll < n; ll += 4) {
                kernel_dtrsm_nn_ll_one_4x4_vs_lib4(i0, &pD[i0 * sdd], &pD[ll * ps], sdd, &d1, &pD[ll * ps + i0 * sdd], &pD[ll * ps + i0 * sdd], &pD[i0 * ps + i0 * sdd], m - i0, n - ll);
            }
            return;
        }
    }
}