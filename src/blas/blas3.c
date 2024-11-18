
#include "tinyhpipm/blas/kernel.h"
#include "tinyhpipm/blas/misc.h"
#include "tinyhpipm/blas/struct.h"
#include <stdio.h>
#include <stdlib.h>

void dgemm_nn(int m, int n, int k, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, double beta, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj) {
    if (m <= 0 || n <= 0) {
        return;
    }

    // invalidate stored inverse diagonal of result matrix
    sD->use_dA = 0;

    const int ps = D_PS;

    int sda = sA->cn;
    int sdb = sB->cn;
    int sdc = sC->cn;
    int sdd = sD->cn;

    int air = ai & (ps - 1);
    int bir = bi & (ps - 1);

    // pA, pB point to panels edges
    double* pA = sA->pA + aj * ps + (ai - air) * sda;
    double* pB = sB->pA + bj * ps + (bi - bir) * sdb;
    double* pC = sC->pA + cj * ps;
    double* pD = sD->pA + dj * ps;

    int offsetB = bir;

    int ci0 = ci - air;
    int di0 = di - air;
    int offsetC;
    int offsetD;
    if (ci0 >= 0) {
        pC += ci0 / ps * ps * sdc;
        offsetC = ci0 % ps;
    } else {
        pC += -ps * sdc;
        offsetC = ps + ci0;
    }

    if (di0 >= 0) {
        pD += di0 / ps * ps * sdd;
        offsetD = di0 % ps;
    } else {
        pD += -ps * sdd;
        offsetD = ps + di0;
    }

    int i, j;


    // algorithm scheme
    if (air != 0) {
        // clean up at the beginning
        j = 0;
        for (; j < n; j += 4) {
            kernel_dgemm_nn_4x4_gen_lib4(k, &alpha, &pA[0], offsetB, &pB[j * ps], sdb, &beta, offsetC, &pC[j * ps], sdc, offsetD, &pD[j * ps], sdd, air, air + m, 0, n - j);
        }
        m -= 1 * ps - air;
        pA += 1 * ps * sda;
        pC += 1 * ps * sdc;
        pD += 1 * ps * sdd;
    }
    if (offsetC == 0 & offsetD == 0) {
        // main loop aligned
        i = 0;
        for (; i < m - 3; i += 4) {
            j = 0;
            for (; j < n - 3; j += 4) {
                kernel_dgemm_nn_4x4_lib4(k, &alpha, &pA[i * sda], offsetB, &pB[j * ps], sdb, &beta, &pC[j * ps + i * sdc], &pD[j * ps + i * sdd]);
            }
            if (j < n) {
                kernel_dgemm_nn_4x4_vs_lib4(k, &alpha, &pA[i * sda], offsetB, &pB[j * ps], sdb, &beta, &pC[j * ps + i * sdc], &pD[j * ps + i * sdd], m - i, n - j);
            }
        }
        if (m > i) {
            // clean up loops definitions
            j = 0;
            for (; j < n; j += 4) {
                kernel_dgemm_nn_4x4_vs_lib4(k, &alpha, &pA[i * sda], offsetB, &pB[j * ps], sdb, &beta, &pC[j * ps + i * sdc], &pD[j * ps + i * sdd], m - i, n - j);
            }
        }
    } else {
        // main loop C, D not aligned
        i = 0;
        for (; i < m; i += 4) {
            j = 0;
            for (; j < n; j += 4) {
                kernel_dgemm_nn_4x4_gen_lib4(k, &alpha, &pA[i * sda], offsetB, &pB[j * ps], sdb, &beta, offsetC, &pC[j * ps + i * sdc], sdc, offsetD, &pD[j * ps + i * sdd], sdd, 0, m - i, 0, n - j);
            }
        }
    }
}


void dgemm_nt(int m, int n, int k, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, double beta, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj) {
    if (m <= 0 | n <= 0) {
        return;
    }
    // invalidate stored inverse diagonal of result matrix
    sD->use_dA = 0;
    const int ps = D_PS;

    int sda = sA->cn;
    int sdb = sB->cn;
    int sdc = sC->cn;
    int sdd = sD->cn;
    int air = ai & (ps - 1);
    int bir = bi & (ps - 1);
    double* pA = sA->pA + aj * ps + (ai - air) * sda;
    double* pB = sB->pA + bj * ps + (bi - bir) * sdb;
    double* pC = sC->pA + cj * ps;
    double* pD = sD->pA + dj * ps;

    int ci0 = ci - air;
    int di0 = di - air;
    int offsetC;
    int offsetD;
    if (ci0 >= 0) {
        pC += ci0 / ps * ps * sdc;
        offsetC = ci0 % ps;
    } else {
        pC += -ps * sdc;
        offsetC = ps + ci0;
    }
    if (di0 >= 0) {
        pD += di0 / ps * ps * sdd;
        offsetD = di0 % ps;
    } else {
        pD += -ps * sdd;
        offsetD = ps + di0;
    }

    int i, j, idxB;
    // algorithm scheme
    if (air != 0) {
        // clean up at the beginning
        j = 0;
        idxB = 0;
        // clean up at the beginning
        if (bir != 0) {
            kernel_dgemm_nt_4x4_gen_lib4(k, &alpha, &pA[0], &pB[idxB * sdb], &beta, offsetC, &pC[j * ps] - bir * ps, sdc, offsetD, &pD[j * ps] - bir * ps, sdd, air, air + m, bir, bir + n - j);
            j += ps - bir;
            idxB += 4;
        }
        // main loop
        for (; j < n; j += 4, idxB += 4) {
            kernel_dgemm_nt_4x4_gen_lib4(k, &alpha, &pA[0], &pB[idxB * sdb], &beta, offsetC, &pC[j * ps], sdc, offsetD, &pD[j * ps], sdd, air, air + m, 0, n - j);
        }
        m -= ps - air;
        pA += ps * sda;
        pC += ps * sdc;
        pD += ps * sdd;
        // TODO instaed use buffer to align A !!!
    }
    if (offsetC == 0 & offsetD == 0) {
        // main loop aligned
        i = 0;
        for (; i < m - 3; i += 4) {
            j = 0;
            idxB = 0;
            // clean up at the beginning
            if (bir != 0) {
                kernel_dgemm_nt_4x4_gen_lib4(k, &alpha, &pA[i * sda], &pB[idxB * sdb], &beta, offsetC, &pC[j * ps + i * sdc] - bir * ps, sdc, offsetD, &pD[j * ps + i * sdd] - bir * ps, sdd, 0, m - i, bir, bir + n - j);
                j += ps - bir;
                idxB += 4;
            }
            // main loop
            for (; j < n - 3; j += 4, idxB += 4) {
                kernel_dgemm_nt_4x4_lib4(k, &alpha, &pA[i * sda], &pB[idxB * sdb], &beta, &pC[j * ps + i * sdc], &pD[j * ps + i * sdd]);
            }
            if (j < n) {
                kernel_dgemm_nt_4x4_vs_lib4(k, &alpha, &pA[i * sda], &pB[idxB * sdb], &beta, &pC[j * ps + i * sdc], &pD[j * ps + i * sdd], m - i, n - j);
            }
        }
        if (m > i) {
            // clean up loops definitions
            j = 0;
            idxB = 0;
            // clean up at the beginning
            if (bir != 0) {
                kernel_dgemm_nt_4x4_gen_lib4(k, &alpha, &pA[i * sda], &pB[idxB * sdb], &beta, offsetC, &pC[j * ps + i * sdc] - bir * ps, sdc, offsetD, &pD[j * ps + i * sdd] - bir * ps, sdd, 0, m - i, bir, bir + n - j);
                j += ps - bir;
                idxB += 4;
            }
            // main loop
            for (; j < n; j += 4, idxB += 4) {
                kernel_dgemm_nt_4x4_vs_lib4(k, &alpha, &pA[i * sda], &pB[idxB * sdb], &beta, &pC[j * ps + i * sdc], &pD[j * ps + i * sdd], m - i, n - j);
            }
        }
    } else {
        // main loop C, D not aligned
        i = 0;

        for (; i < m; i += 4) {
            j = 0;
            idxB = 0;
            // clean up at the beginning
            if (bir != 0) {
                kernel_dgemm_nt_4x4_gen_lib4(k, &alpha, &pA[i * sda], &pB[idxB * sdb], &beta, offsetC, &pC[j * ps + i * sdc] - bir * ps, sdc, offsetD, &pD[j * ps + i * sdd] - bir * ps, sdd, 0, m - i, bir, bir + n - j);
                j += ps - bir;
                idxB += 4;
            }
            // main loop
            for (; j < n; j += 4, idxB += 4) {
                kernel_dgemm_nt_4x4_gen_lib4(k, &alpha, &pA[i * sda], &pB[idxB * sdb], &beta, offsetC, &pC[j * ps + i * sdc], sdc, offsetD, &pD[j * ps + i * sdd], sdd, 0, m - i, 0, n - j);
            }
        }
    }
}

void dgemm_dn(int m, int n, double alpha, struct vec* sA, int ai, struct mat* sB, int bi, int bj, double beta, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj) {
    if (m <= 0 | n <= 0) {
        return;
    }
    if (bi != 0 | ci != 0 | di != 0) {
        fprintf(stderr, "\n\tdgemm_dn: feature not implemented yet: bi=%d, ci=%d, di=%d\n", bi, ci, di);
        exit(1);
    }

    // invalidate stored inverse diagonal of result matrix
    sD->use_dA = 0;

    const int bs = 4;
    int sdb = sB->cn;
    int sdc = sC->cn;
    int sdd = sD->cn;
    double* dA = sA->pa + ai;
    double* pB = sB->pA + bj * bs;
    double* pC = sC->pA + cj * bs;
    double* pD = sD->pA + dj * bs;

    int ii;

    ii = 0;
    if (beta == 0.0) {
        for (; ii < m - 3; ii += 4) {
            kernel_dgemm_diag_left_4_a0_lib4(n, &alpha, &dA[ii], &pB[ii * sdb], &pD[ii * sdd]);
        }
    } else {
        for (; ii < m - 3; ii += 4) {
            kernel_dgemm_diag_left_4_lib4(n, &alpha, &dA[ii], &pB[ii * sdb], &beta, &pC[ii * sdc], &pD[ii * sdd]);
        }
    }
    if (m - ii > 0) {
        if (m - ii == 1)
            kernel_dgemm_diag_left_1_lib4(n, &alpha, &dA[ii], &pB[ii * sdb], &beta, &pC[ii * sdc], &pD[ii * sdd]);
        else if (m - ii == 2)
            kernel_dgemm_diag_left_2_lib4(n, &alpha, &dA[ii], &pB[ii * sdb], &beta, &pC[ii * sdc], &pD[ii * sdd]);
        else  // if(m-ii==3)
            kernel_dgemm_diag_left_3_lib4(n, &alpha, &dA[ii], &pB[ii * sdb], &beta, &pC[ii * sdc], &pD[ii * sdd]);
    }
}


// dgemm with B diagonal matrix (stored as strvec)
void dgemm_nd(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct vec* sB, int bi, double beta, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj) {
    if (m <= 0 | n <= 0) {
        return;
    }
    if (ai != 0 | ci != 0 | di != 0) {
        fprintf(stderr, "\n\tdgemm_nd: feature not implemented yet: ai=%d, ci=%d, di=%d\n", ai, ci, di);
        exit(1);
    }

    // invalidate stored inverse diagonal of result matrix
    sD->use_dA = 0;

    const int bs = 4;
    int sda = sA->cn;
    int sdc = sC->cn;
    int sdd = sD->cn;
    double* pA = sA->pA + aj * bs;
    double* dB = sB->pa + bi;
    double* pC = sC->pA + cj * bs;
    double* pD = sD->pA + dj * bs;

    int ii;

    ii = 0;
    if (beta == 0.0) {
        for (; ii < n - 3; ii += 4) {
            kernel_dgemm_diag_right_4_a0_lib4(m, &alpha, &pA[ii * bs], sda, &dB[ii], &pD[ii * bs], sdd);
        }
    } else {
        for (; ii < n - 3; ii += 4) {
            kernel_dgemm_diag_right_4_lib4(m, &alpha, &pA[ii * bs], sda, &dB[ii], &beta, &pC[ii * bs], sdc, &pD[ii * bs], sdd);
        }
    }
    if (n - ii > 0) {
        if (n - ii == 1)
            kernel_dgemm_diag_right_1_lib4(m, &alpha, &pA[ii * bs], sda, &dB[ii], &beta, &pC[ii * bs], sdc, &pD[ii * bs], sdd);
        else if (n - ii == 2)
            kernel_dgemm_diag_right_2_lib4(m, &alpha, &pA[ii * bs], sda, &dB[ii], &beta, &pC[ii * bs], sdc, &pD[ii * bs], sdd);
        else  // if(n-ii==3)
            kernel_dgemm_diag_right_3_lib4(m, &alpha, &pA[ii * bs], sda, &dB[ii], &beta, &pC[ii * bs], sdc, &pD[ii * bs], sdd);
    }
}

// dtrsm_right_lower_transposed_notunit
void dtrsm_rltn(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj) {
    if (m <= 0 || n <= 0) {
        return;
    }

    const int ps = D_PS;

    // invalidate stored inverse diagonal of result matrix
    sD->use_dA = 0;

    // TODO alpha !!!!!

    int sda = sA->cn;
    int sdb = sB->cn;
    int sdd = sD->cn;
    int bir = bi & (ps - 1);
    int dir = di & (ps - 1);
    double* pA = sA->pA + aj * ps;
    double* pB = sB->pA + bj * ps + (bi - bir) * sdb;
    double* pD = sD->pA + dj * ps + (di - dir) * sdd;
    double* dA = sA->dA;

    if (ai != 0 | bir != 0 | dir != 0 | alpha != 1.0) {
        fprintf(stderr, "\n\tdtrsm_rltn: feature not implemented yet: ai=%d, bi=%d, di=%d, alpha=%f\n", ai, bi, di, alpha);
        exit(1);
    }

    int i, j;
    // TODO to avoid touching A, better temporarely use sD.dA ?????
    if (ai == 0 & aj == 0) {
        if (sA->use_dA < n) {
            ddiaex_lib(n, 1.0, ai, pA, sda, dA);
            for (i = 0; i < n; i++)
                dA[i] = 1.0 / dA[i];
            sA->use_dA = n;
        }
    } else {
        ddiaex_lib(n, 1.0, ai, pA, sda, dA);
        for (i = 0; i < n; i++)
            dA[i] = 1.0 / dA[i];
        sA->use_dA = 0;
    }

    i = 0;
    for (; i < m - 3; i += 4) {
        j = 0;
        for (; j < n - 3; j += 4) {
            kernel_dtrsm_nt_rl_inv_4x4_lib4(j, &pD[i * sdd], &pA[j * sda], &alpha, &pB[j * ps + i * sdb], &pD[j * ps + i * sdd], &pA[j * ps + j * sda], &dA[j]);
        }
        if (j < n) {
            kernel_dtrsm_nt_rl_inv_4x4_vs_lib4(j, &pD[i * sdd], &pA[j * sda], &alpha, &pB[j * ps + i * sdb], &pD[j * ps + i * sdd], &pA[j * ps + j * sda], &dA[j], m - i, n - j);
        }
    }
    if (m > i) {
        j = 0;
        for (; j < n; j += 4) {
            kernel_dtrsm_nt_rl_inv_4x4_vs_lib4(j, &pD[i * sdd], &pA[j * sda], &alpha, &pB[j * ps + i * sdb], &pD[j * ps + i * sdd], &pA[j * ps + j * sda], &dA[j], m - i, n - j);
        }
    }
}

// dtrmm_right_lower_nottransposed_notunit (B, i.e. the first matrix, is triangular !!!)
void dtrmm_rlnn(int m, int n, double alpha, struct mat* sB, int bi, int bj, struct mat* sA, int ai, int aj, struct mat* sD, int di, int dj) {

    const int ps = 4;

    int sda = sA->cn;
    int sdb = sB->cn;
    int sdd = sD->cn;
    int air = ai & (ps - 1);
    int bir = bi & (ps - 1);
    double* pA = sA->pA + aj * ps + (ai - air) * sda;
    double* pB = sB->pA + bj * ps + (bi - bir) * sdb;
    double* pD = sD->pA + dj * ps;

    int offsetB = bir;

    int di0 = di - air;
    int offsetD;

    // invalidate stored inverse diagonal of result matrix
    sD->use_dA = 0;

    if (di0 >= 0) {
        pD += di0 / ps * ps * sdd;
        offsetD = di0 % ps;
    } else {
        pD += -ps * sdd;
        offsetD = ps + di0;
    }

    int ii, jj;

    // algorithm scheme
    if (air != 0) {
        jj = 0;
        for (; jj < n; jj += 4) {
            kernel_dtrmm_nn_rl_4x4_gen_lib4(n - jj, &alpha, &pA[jj * ps], offsetB, &pB[jj * sdb + jj * ps], sdb, offsetD, &pD[jj * ps], sdd, air, air + m, 0, n - jj);
        }
        m -= ps - air;
        pA += ps * sda;
        pD += ps * sdd;
        goto select_loop;
    }
select_loop:
    if (offsetD == 0) {
        ii = 0;
        for (; ii < m - 3; ii += 4) {
            jj = 0;
            for (; jj < n - 5; jj += 4) {
                kernel_dtrmm_nn_rl_4x4_lib4(n - jj, &alpha, &pA[ii * sda + jj * ps], offsetB, &pB[jj * sdb + jj * ps], sdb, &pD[ii * sdd + jj * ps]);
            }
            for (; jj < n; jj += 4) {
                kernel_dtrmm_nn_rl_4x4_vs_lib4(n - jj, &alpha, &pA[ii * sda + jj * ps], offsetB, &pB[jj * sdb + jj * ps], sdb, &pD[ii * sdd + jj * ps], 4, n - jj);
            }
        }
    } else {
        // main loop C, D not aligned
        ii = 0;
        for (; ii < m; ii += 4) {
            jj = 0;
            for (; jj < n; jj += 4) {
                kernel_dtrmm_nn_rl_4x4_gen_lib4(n - jj, &alpha, &pA[ii * sda + jj * ps], offsetB, &pB[jj * sdb + jj * ps], sdb, offsetD, &pD[ii * sdd + jj * ps], sdd, 0, m - ii, 0, n - jj);
            }
        }
    }


    // clear
    if (ii < m) {
        jj = 0;
        for (; jj < n; jj += 4) {
            kernel_dtrmm_nn_rl_4x4_vs_lib4(n - jj, &alpha, &pA[ii * sda + jj * ps], offsetB, &pB[jj * sdb + jj * ps], sdb, &pD[ii * sdd + jj * ps], m - ii, n - jj);
        }
    }
}

void dsyrk_ln(int m, int k, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, double beta, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj) {
    // fast return
    if (m <= 0) {
        return;
    }
    // invalidate stored inverse diagonal of result matrix
    sD->use_dA = 0;
    const int ps = 4;

    int sda = sA->cn;
    int sdb = sB->cn;
    int sdc = sC->cn;
    int sdd = sD->cn;
    int air = ai & (ps - 1);
    int bir = bi & (ps - 1);
    double* pA = sA->pA + aj * ps + (ai - air) * sda;
    double* pB = sB->pA + bj * ps + (bi - bir) * sdb;
    double* pC = sC->pA + cj * ps;
    double* pD = sD->pA + dj * ps;

    int ci0 = ci;  //-air;
    int di0 = di;  //-air;
    int offsetC;
    int offsetD;
    if (ci0 >= 0) {
        pC += ci0 / ps * ps * sdc;
        offsetC = ci0 % ps;
    } else {
        pC += -ps * sdc;
        offsetC = ps + ci0;
    }
    if (di0 >= 0) {
        pD += di0 / ps * ps * sdd;
        offsetD = di0 % ps;
    } else {
        pD += -ps * sdd;
        offsetD = ps + di0;
    }

    struct mat sAt;
    int sAt_size;
    void* mem;
    char* mem_align;

    double *pU, *pA2;

    double pU0[1 * 4 * K_MAX_STACK];
    int sdu0 = (k + 3) / 4 * 4;
    sdu0 = sdu0 < K_MAX_STACK ? sdu0 : K_MAX_STACK;

    // allocate memory
    if (k > K_MAX_STACK) {
        sAt_size = memsize_mat(12, k);
        mem = malloc(sAt_size + 64);
        align_64_byte(mem, (void**) &mem_align);
        create_mat(12, k, &sAt, (void*) mem_align);
        pU = sAt.pA;
    } else {
        pU = pU0;
    }


    int i, j;

    int idxB;


    // algorithm scheme
    if (offsetC == 0 & offsetD == 0) {
        if (bir == 0) {
            // main loop aligned
            i = 0;
            for (; i < m - 3; i += 4) {
                if (air == 0) {
                    pA2 = pA + i * sda;
                } else {
                    kernel_dpacp_nn_4_lib4(k, air, pA + i * sda, sda, pU);
                    pA2 = pU;
                }
                j = 0;
                // main loop
                for (; j < i; j += 4) {
                    kernel_dgemm_nt_4x4_lib4(k, &alpha, pA2, &pB[j * sdb], &beta, &pC[j * ps + i * sdc], &pD[j * ps + i * sdd]);
                }
                kernel_dsyrk_nt_l_4x4_lib4(k, &alpha, pA2, &pB[j * sdb], &beta, &pC[j * ps + i * sdc], &pD[j * ps + i * sdd]);
            }
            if (m > i) {
                if (air == 0) {
                    pA2 = pA + i * sda;
                } else {
                    kernel_dpacp_nn_4_vs_lib4(k, air, pA + i * sda, sda, pU, m - i);
                    pA2 = pU;
                }
                j = 0;
                // main loop
                for (; j < i; j += 4) {
                    kernel_dgemm_nt_4x4_vs_lib4(k, &alpha, pA2, &pB[j * sdb], &beta, &pC[j * ps + i * sdc], &pD[j * ps + i * sdd], m - i, m - j);
                }
                kernel_dsyrk_nt_l_4x4_vs_lib4(k, &alpha, pA2, &pB[j * sdb], &beta, &pC[j * ps + i * sdc], &pD[j * ps + i * sdd], m - i, m - j);
            }
        } else {
            // main loop aligned
            i = 0;
            for (; i < m - 3; i += 4) {
                if (air == 0) {
                    pA2 = pA + i * sda;
                } else {
                    kernel_dpacp_nn_4_lib4(k, air, pA + i * sda, sda, pU);
                    pA2 = pU;
                }
                j = 0;
                idxB = 0;
                if (j < i) {
                    kernel_dgemm_nt_4x4_gen_lib4(k, &alpha, pA2, &pB[idxB * sdb], &beta, 0, &pC[j * ps + i * sdc] - bir * ps, sdc, 0, &pD[j * ps + i * sdd] - bir * ps, sdd, 0, m - i, bir, 4);
                    j += ps - bir;
                    idxB += 4;
                    // main loop
                    for (; j < i + (ps - bir) - ps; j += 4, idxB += 4) {
                        kernel_dgemm_nt_4x4_lib4(k, &alpha, pA2, &pB[idxB * sdb], &beta, &pC[j * ps + i * sdc], &pD[j * ps + i * sdd]);
                    }
                    kernel_dgemm_nt_4x4_vs_lib4(k, &alpha, pA2, &pB[idxB * sdb], &beta, &pC[j * ps + i * sdc], &pD[j * ps + i * sdd], m - i, bir);
                    j += bir;
                }
                kernel_dsyrk_nt_l_4x4_gen_lib4(k, &alpha, pA2, &pB[j * sdb], &beta, 0, &pC[j * ps + i * sdc] - bir * ps, sdc, 0, &pD[j * ps + i * sdd] - bir * ps, sdd, 0, m - i, bir, 4);
                kernel_dsyrk_nt_l_4x4_gen_lib4(k, &alpha, pA2, &pB[(j + 4) * sdb], &beta, 0, &pC[j * ps + i * sdc] + (ps - bir) * ps, sdc, 0, &pD[j * ps + i * sdd] + (ps - bir) * ps, sdd, ps - bir, m - i, 0, 4);
            }
            if (m > i) {
                j = 0;
                if (air == 0) {
                    pA2 = pA + i * sda;
                } else {
                    kernel_dpacp_nn_4_vs_lib4(k, air, pA + i * sda, sda, pU, m - i);
                    pA2 = pU;
                }
                if (bir != 0) {
                    idxB = 0;
                    if (j < i) {
                        kernel_dgemm_nt_4x4_gen_lib4(k, &alpha, pA2, &pB[idxB * sdb], &beta, offsetC, &pC[j * ps + i * sdc] - bir * ps, sdc, offsetD, &pD[j * ps + i * sdd] - bir * ps, sdd, 0, m - i, bir, m - j);
                        j += ps - bir;
                        idxB += 4;
                        // main loop
                        for (; j < i + (ps - bir) - ps; j += 4, idxB += 4) {
                            kernel_dgemm_nt_4x4_gen_lib4(k, &alpha, pA2, &pB[idxB * sdb], &beta, offsetC, &pC[j * ps + i * sdc], sdc, offsetD, &pD[j * ps + i * sdd], sdd, 0, m - i, 0, m - j);
                        }
                        kernel_dgemm_nt_4x4_gen_lib4(k, &alpha, pA2, &pB[idxB * sdb], &beta, offsetC, &pC[j * ps + i * sdc], sdc, offsetD, &pD[j * ps + i * sdd], sdd, 0, m - i, 0, bir);  // XXX n1
                        j += bir;
                    }
                    kernel_dsyrk_nt_l_4x4_gen_lib4(k, &alpha, pA2, &pB[j * sdb], &beta, offsetC, &pC[j * ps + i * sdc] - bir * ps, sdc, offsetD, &pD[j * ps + i * sdd] - bir * ps, sdd, 0, m - i, bir, bir + m - j);
                    kernel_dsyrk_nt_l_4x4_gen_lib4(k, &alpha, pA2, &pB[(j + 4) * sdb], &beta, offsetC, &pC[j * ps + i * sdc] + (ps - bir) * ps, sdc, offsetD, &pD[j * ps + i * sdd] + (ps - bir) * ps, sdd, ps - bir, m - i, 0, m - j);
                } else {
                    // main loop
                    for (; j < i; j += 4) {
                        kernel_dgemm_nt_4x4_gen_lib4(k, &alpha, pA2, &pB[j * sdb], &beta, offsetC, &pC[j * ps + i * sdc], sdc, offsetD, &pD[j * ps + i * sdd], sdd, 0, m - i, 0, m - j);
                    }
                    kernel_dsyrk_nt_l_4x4_gen_lib4(k, &alpha, pA2, &pB[j * sdb], &beta, offsetC, &pC[j * ps + i * sdc], sdc, offsetD, &pD[j * ps + i * sdd], sdd, 0, m - i, 0, m - j);
                }
            }
        }
    } else {
        if (bir == 0) {
            // main loop C, D not aligned
            i = 0;
            for (; i < m; i += 4) {
                if (air == 0) {
                    pA2 = pA + i * sda;
                } else {
                    kernel_dpacp_nn_4_vs_lib4(k, air, pA + i * sda, sda, pU, m - i);
                    pA2 = pU;
                }
                j = 0;
                // main loop
                for (; j < i; j += 4) {
                    kernel_dgemm_nt_4x4_gen_lib4(k, &alpha, pA2, &pB[j * sdb], &beta, offsetC, &pC[j * ps + i * sdc], sdc, offsetD, &pD[j * ps + i * sdd], sdd, 0, m - i, 0, m - j);
                }
                kernel_dsyrk_nt_l_4x4_gen_lib4(k, &alpha, pA2, &pB[j * sdb], &beta, offsetC, &pC[j * ps + i * sdc], sdc, offsetD, &pD[j * ps + i * sdd], sdd, 0, m - i, 0, m - j);
            }
        } else {
            // main loop aligned
            i = 0;
            for (; i < m; i += 4) {
                if (air == 0) {
                    pA2 = pA + i * sda;
                } else {
                    kernel_dpacp_nn_4_vs_lib4(k, air, pA + i * sda, sda, pU, m - i);
                    pA2 = pU;
                }
                j = 0;
                idxB = 0;
                if (j < i) {
                    kernel_dgemm_nt_4x4_gen_lib4(k, &alpha, pA2, &pB[idxB * sdb], &beta, offsetC, &pC[j * ps + i * sdc] - bir * ps, sdc, offsetD, &pD[j * ps + i * sdd] - bir * ps, sdd, 0, m - i, bir, m - j);
                    j += ps - bir;
                    idxB += 4;
                    // main loop
                    for (; j < i + (ps - bir) - ps; j += 4, idxB += 4) {
                        kernel_dgemm_nt_4x4_gen_lib4(k, &alpha, pA2, &pB[idxB * sdb], &beta, offsetC, &pC[j * ps + i * sdc], sdc, offsetD, &pD[j * ps + i * sdd], sdd, 0, m - i, 0, m - j);
                    }
                    kernel_dgemm_nt_4x4_gen_lib4(k, &alpha, pA2, &pB[idxB * sdb], &beta, offsetC, &pC[j * ps + i * sdc], sdc, offsetD, &pD[j * ps + i * sdd], sdd, 0, m - i, 0, bir);  // XXX n1
                    j += bir;
                }
                kernel_dsyrk_nt_l_4x4_gen_lib4(k, &alpha, pA2, &pB[j * sdb], &beta, offsetC, &pC[j * ps + i * sdc] - bir * ps, sdc, offsetD, &pD[j * ps + i * sdd] - bir * ps, sdd, 0, m - i, bir, bir + m - j);
                kernel_dsyrk_nt_l_4x4_gen_lib4(k, &alpha, pA2, &pB[(j + 4) * sdb], &beta, offsetC, &pC[j * ps + i * sdc] + (ps - bir) * ps, sdc, offsetD, &pD[j * ps + i * sdd] + (ps - bir) * ps, sdd, ps - bir, m - i, 0, m - j);
            }
        }
    }
    if (k > K_MAX_STACK) {
        free(mem);
    }
}


void dsyrk_ln_mn(int m, int n, int k, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, double beta, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj) {
    if (m <= 0 | n <= 0)
        return;

    if (ai != 0 | bi != 0) {
        fprintf(stderr, "\n\tdsyrk_ln_mn: feature not implemented yet: ai=%d, bi=%d\n", ai, bi);
        exit(1);
    }

    // invalidate stored inverse diagonal of result matrix
    sD->use_dA = 0;

    const int ps = 4;

    int i, j;

    int sda = sA->cn;
    int sdb = sB->cn;
    int sdc = sC->cn;
    int sdd = sD->cn;
    double* pA = sA->pA + aj * ps;
    double* pB = sB->pA + bj * ps;
    double* pC = sC->pA + cj * ps + (ci - (ci & (ps - 1))) * sdc;
    double* pD = sD->pA + dj * ps + (di - (di & (ps - 1))) * sdd;

    // TODO ai and bi
    int offsetC;
    int offsetD;
    offsetC = ci & (ps - 1);
    offsetD = di & (ps - 1);


    // algorithm scheme
    if (offsetC == 0 & offsetD == 0) {
        // main loop aligned
        i = 0;

        for (; i < m - 3; i += 4) {
            j = 0;
            for (; j < i & j < n - 3; j += 4) {
                kernel_dgemm_nt_4x4_lib4(k, &alpha, &pA[i * sda], &pB[j * sdb], &beta, &pC[j * ps + i * sdc], &pD[j * ps + i * sdd]);
            }
            if (j < n) {
                if (j < i)  // dgemm
                {
                    kernel_dgemm_nt_4x4_vs_lib4(k, &alpha, &pA[i * sda], &pB[j * sdb], &beta, &pC[j * ps + i * sdc], &pD[j * ps + i * sdd], m - i, n - j);
                } else  // dsyrk
                {
                    if (j < n - 3) {
                        kernel_dsyrk_nt_l_4x4_lib4(k, &alpha, &pA[i * sda], &pB[j * sdb], &beta, &pC[j * ps + i * sdc], &pD[j * ps + i * sdd]);
                    } else {
                        kernel_dsyrk_nt_l_4x4_vs_lib4(k, &alpha, &pA[i * sda], &pB[j * sdb], &beta, &pC[j * ps + i * sdc], &pD[j * ps + i * sdd], m - i, n - j);
                    }
                }
            }
        }
        if (m > i) {
            // clean up loops definitions
            j = 0;
            for (; j < i & j < n; j += 4) {
                kernel_dgemm_nt_4x4_vs_lib4(k, &alpha, &pA[i * sda], &pB[j * sdb], &beta, &pC[j * ps + i * sdc], &pD[j * ps + i * sdd], m - i, n - j);
            }
            if (j < n) {
                kernel_dsyrk_nt_l_4x4_vs_lib4(k, &alpha, &pA[i * sda], &pB[j * sdb], &beta, &pC[j * ps + i * sdc], &pD[j * ps + i * sdd], m - i, n - j);
            }
        }
    } else {
        // main loop C, D not aligned
        i = 0;
        for (; i < m; i += 4) {
            j = 0;
            for (; j < i & j < n; j += 4) {
                kernel_dgemm_nt_4x4_gen_lib4(k, &alpha, &pA[i * sda], &pB[j * sdb], &beta, offsetC, &pC[j * ps + i * sdc], sdc, offsetD, &pD[j * ps + i * sdd], sdd, 0, m - i, 0, n - j);
            }
            if (j < n) {
                kernel_dsyrk_nt_l_4x4_gen_lib4(k, &alpha, &pA[i * sda], &pB[j * sdb], &beta, offsetC, &pC[j * ps + i * sdc], sdc, offsetD, &pD[j * ps + i * sdd], sdd, 0, m - i, 0, n - j);
            }
        }
    }
}
