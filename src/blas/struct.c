#include "tinyhpipm/blas/struct.h"
#include <stdio.h>
#include <stdlib.h>

size_t memsize_mat(int m, int n) {
    const int al = D_PS * D_PLD;
    int pm = (m + D_PS - 1) & ~(D_PS - 1);  // require  pm % D_PS == 0
    int cn = (n + D_PLD - 1) & ~(D_PLD - 1);  // require  cn % D_PLD == 0
    int tmp = m < n ? (m + al - 1) & ~(al - 1) : (n + al - 1) & ~(al - 1);  // al(min(m,n)) // XXX max ???
    size_t memsize = (pm * cn + tmp) * sizeof(double);
    // TODO: check if we don't need to include cache line alignment
    return memsize;
}
size_t memsize_diag_mat(int m, int n) {
    const int al = D_PS * D_PLD;
    int tmp = m < n ? (m + al - 1) / al * al : (n + al - 1) / al * al;  // al(min(m,n)) // XXX max ???
    size_t memsize = tmp * sizeof(double);
    return memsize;
}

size_t memsize_vec(int m) {
    int pm = (m + D_PS - 1) & ~(D_PS - 1);  // require  pm % D_PS == 0
    size_t memsize = pm * sizeof(double);
    return memsize;
}

void create_mat(int m, int n, struct mat* sA, void* memory) {
    sA->mem = memory;  // the memory has been created using mem_size
    int al = D_PS * D_PLD;
    sA->m = m;
    sA->n = n;
    int pm = (m + D_PS - 1) & ~(D_PS - 1);  // require  pm % D_PS == 0
    int cn = (n + D_PLD - 1) & ~(D_PLD - 1);  // require  cn % D_PLD == 0
    sA->pm = pm;
    sA->cn = cn;
    double* ptr = (double*) memory;
    sA->pA = ptr;
    ptr += pm * cn;
    int tmp = m < n ? (m + al - 1) / al * al : (n + al - 1) / al * al;  // al(min(m,n)) // XXX max ???
    sA->dA = ptr;
    ptr += tmp;
    sA->memsize = (pm * cn + tmp) * sizeof(double);
    sA->use_dA = 0;  // invalidate stored inverse diagonal
}
void create_vec(int m, struct vec* sa, void* memory) {
    sa->mem = memory;
    sa->m = m;
    int pm = (m + D_PS - 1) & ~(D_PS - 1);  // require  pm % D_PS == 0
    sa->pm = pm;
    double* ptr = (double*) memory;
    sa->pa = ptr;
    //	ptr += pm;
    sa->memsize = pm * sizeof(double);
}

/*
void allocate_mat(int m, int n, struct mat* sA) {
    size_t size = memsize_mat(m, n);
    void* mem = aligned_alloc(CACHE_LINE_SIZE, size);
    create_mat(m, n, sA, mem);
}

void allocate_vec(int m, struct vec* sa) {
    size_t size = memsize_vec(m);
    void* mem = aligned_alloc(CACHE_LINE_SIZE, size);
    create_vec(m, sa, mem);
}
*/


void free_mat(struct mat* sA) {
    free(sA->mem);
}
void free_vec(struct vec* sa) {
    free(sa->mem);
}


// convert a matrix into a matrix structure
void pack_mat(int m, int n, double* A, int lda, struct mat* sA, int ai, int aj) {
    if (m <= 0 || n <= 0) {
        return;
    }

    // invalidate stored inverse diagonal
    sA->use_dA = 0;

    const int bs = 4;
    int sda = sA->cn;
    double* pA = sA->pA + aj * bs + ai / bs * bs * sda + ai % bs;
    int i, ii, j, jj, m0, m1, m2;
    double *B, *pB;
    sA->use_dA = 0;

    // row vector in sA
    if (m == 1) {
        for (jj = 0; jj < n; jj++) {
            pA[jj * bs] = A[jj * lda];
        }
        return;
    }

    m0 = (bs - ai % bs) % bs;
    if (m0 > m)
        m0 = m;
    m1 = m - m0;
    jj = 0;
    for (; jj < n - 3; jj += 4) {
        B = A + jj * lda;
        pB = pA + jj * bs;
        ii = 0;
        if (m0 > 0) {
            for (; ii < m0; ii++) {
                pB[ii + bs * 0] = B[ii + lda * 0];
                pB[ii + bs * 1] = B[ii + lda * 1];
                pB[ii + bs * 2] = B[ii + lda * 2];
                pB[ii + bs * 3] = B[ii + lda * 3];
            }
            B += m0;
            pB += m0 + bs * (sda - 1);
        }

        for (; ii < m - 3; ii += 4) {
            // col 0
            pB[0 + bs * 0] = B[0 + lda * 0];
            pB[1 + bs * 0] = B[1 + lda * 0];
            pB[2 + bs * 0] = B[2 + lda * 0];
            pB[3 + bs * 0] = B[3 + lda * 0];
            // col 1
            pB[0 + bs * 1] = B[0 + lda * 1];
            pB[1 + bs * 1] = B[1 + lda * 1];
            pB[2 + bs * 1] = B[2 + lda * 1];
            pB[3 + bs * 1] = B[3 + lda * 1];
            // col 2
            pB[0 + bs * 2] = B[0 + lda * 2];
            pB[1 + bs * 2] = B[1 + lda * 2];
            pB[2 + bs * 2] = B[2 + lda * 2];
            pB[3 + bs * 2] = B[3 + lda * 2];
            // col 3
            pB[0 + bs * 3] = B[0 + lda * 3];
            pB[1 + bs * 3] = B[1 + lda * 3];
            pB[2 + bs * 3] = B[2 + lda * 3];
            pB[3 + bs * 3] = B[3 + lda * 3];
            // update
            B += 4;
            pB += bs * sda;
        }
        for (; ii < m; ii++) {
            // col 0
            pB[0 + bs * 0] = B[0 + lda * 0];
            // col 1
            pB[0 + bs * 1] = B[0 + lda * 1];
            // col 2
            pB[0 + bs * 2] = B[0 + lda * 2];
            // col 3
            pB[0 + bs * 3] = B[0 + lda * 3];
            // update
            B += 1;
            pB += 1;
        }
    }
    for (; jj < n; jj++) {

        B = A + jj * lda;
        pB = pA + jj * bs;

        ii = 0;
        if (m0 > 0) {
            for (; ii < m0; ii++) {
                pB[ii + bs * 0] = B[ii + lda * 0];
            }
            B += m0;
            pB += m0 + bs * (sda - 1);
        }
        for (; ii < m - 3; ii += 4) {
            // col 0
            pB[0 + bs * 0] = B[0 + lda * 0];
            pB[1 + bs * 0] = B[1 + lda * 0];
            pB[2 + bs * 0] = B[2 + lda * 0];
            pB[3 + bs * 0] = B[3 + lda * 0];
            // update
            B += 4;
            pB += bs * sda;
        }
        for (; ii < m; ii++) {
            // col 0
            pB[0 + bs * 0] = B[0 + lda * 0];
            // update
            B += 1;
            pB += 1;
        }
    }
}


// convert a lower triangular matrix into a matrix structure
// TODO	vectorize triangle in the 4x4 diag blocks
void pack_l_mat(int m, int n, double* A, int lda, struct mat* sA, int ai, int aj) {
    if (m <= 0 || n <= 0) {
        return;
    }

    // invalidate stored inverse diagonal
    sA->use_dA = 0;

    const int bs = 4;
    int sda = sA->cn;
    double* pA = sA->pA + aj * bs + ai / bs * bs * sda + ai % bs;
    int i, ii, j, jj, m0, m1, m2;
    double *B, *pB;
    sA->use_dA = 0;

    // row vector in sA
    if (m == 1) {
        for (jj = 0; jj < n; jj++) {
            pA[jj * bs] = A[jj * lda];
        }
        return;
    }


    m0 = (bs - ai % bs) % bs;
    if (m0 > m)
        m0 = m;
    m1 = m - m0;
    if (m0 > 0) {
        fprintf(stderr, "\n\tpack_l_mat: feature not implemented yet: ai!=0\n");
        exit(1);
    }
    jj = 0;
    for (; jj < n - 3; jj += 4) {
        B = A + jj + jj * lda;
        pB = pA + jj * bs + jj * sda;
        ii = jj;
        // col 0
        pB[0 + bs * 0] = B[0 + lda * 0];
        pB[1 + bs * 0] = B[1 + lda * 0];
        pB[2 + bs * 0] = B[2 + lda * 0];
        pB[3 + bs * 0] = B[3 + lda * 0];
        // col 1
        //		pB[0+bs*1] = B[0+lda*1];
        pB[1 + bs * 1] = B[1 + lda * 1];
        pB[2 + bs * 1] = B[2 + lda * 1];
        pB[3 + bs * 1] = B[3 + lda * 1];
        // col 2
        //		pB[0+bs*2] = B[0+lda*2];
        //		pB[1+bs*2] = B[1+lda*2];
        pB[2 + bs * 2] = B[2 + lda * 2];
        pB[3 + bs * 2] = B[3 + lda * 2];
        // col 3
        //		pB[0+bs*3] = B[0+lda*3];
        //		pB[1+bs*3] = B[1+lda*3];
        //		pB[2+bs*3] = B[2+lda*3];
        pB[3 + bs * 3] = B[3 + lda * 3];
        B += 4;
        pB += bs * sda;
        ii += 4;

        for (; ii < m - 3; ii += 4) {
            // col 0
            pB[0 + bs * 0] = B[0 + lda * 0];
            pB[1 + bs * 0] = B[1 + lda * 0];
            pB[2 + bs * 0] = B[2 + lda * 0];
            pB[3 + bs * 0] = B[3 + lda * 0];
            // col 1
            pB[0 + bs * 1] = B[0 + lda * 1];
            pB[1 + bs * 1] = B[1 + lda * 1];
            pB[2 + bs * 1] = B[2 + lda * 1];
            pB[3 + bs * 1] = B[3 + lda * 1];
            // col 2
            pB[0 + bs * 2] = B[0 + lda * 2];
            pB[1 + bs * 2] = B[1 + lda * 2];
            pB[2 + bs * 2] = B[2 + lda * 2];
            pB[3 + bs * 2] = B[3 + lda * 2];
            // col 3
            pB[0 + bs * 3] = B[0 + lda * 3];
            pB[1 + bs * 3] = B[1 + lda * 3];
            pB[2 + bs * 3] = B[2 + lda * 3];
            pB[3 + bs * 3] = B[3 + lda * 3];
            // update
            B += 4;
            pB += bs * sda;
        }
        for (; ii < m; ii++) {
            // col 0
            pB[0 + bs * 0] = B[0 + lda * 0];
            // col 1
            pB[0 + bs * 1] = B[0 + lda * 1];
            // col 2
            pB[0 + bs * 2] = B[0 + lda * 2];
            // col 3
            pB[0 + bs * 3] = B[0 + lda * 3];
            // update
            B += 1;
            pB += 1;
        }
    }
    if (jj < n) {
        B = A + jj + jj * lda;
        pB = pA + jj * bs + jj * sda;
        if (n - jj == 1) {
            // col 0
            pB[0 + bs * 0] = B[0 + lda * 0];
        } else if (n - jj == 2) {
            // col 0
            pB[0 + bs * 0] = B[0 + lda * 0];
            pB[1 + bs * 0] = B[1 + lda * 0];
            // col 1
            //			pB[0+bs*1] = B[0+lda*1];
            pB[1 + bs * 1] = B[1 + lda * 1];
        } else  // if(n-jj==3)
        {
            // col 0
            pB[0 + bs * 0] = B[0 + lda * 0];
            pB[1 + bs * 0] = B[1 + lda * 0];
            pB[2 + bs * 0] = B[2 + lda * 0];
            // col 1
            //			pB[0+bs*1] = B[0+lda*1];
            pB[1 + bs * 1] = B[1 + lda * 1];
            pB[2 + bs * 1] = B[2 + lda * 1];
            // col 2
            //			pB[0+bs*2] = B[0+lda*2];
            //			pB[1+bs*2] = B[1+lda*2];
            pB[2 + bs * 2] = B[2 + lda * 2];
        }
    }
}


// convert a upper triangular matrix into a matrix structure
// TODO	vectorize triangle in the 4x4 diag blocks
void pack_u_mat(int m, int n, double* A, int lda, struct mat* sA, int ai, int aj) {
    if (m <= 0 || n <= 0) {
        return;
    }

    // invalidate stored inverse diagonal
    sA->use_dA = 0;

    const int bs = 4;
    int sda = sA->cn;
    double* pA = sA->pA + aj * bs + ai / bs * bs * sda + ai % bs;
    int i, ii, j, jj, m0, m1, m2;
    double *B, *pB;
    sA->use_dA = 0;

    // row vector in sA
    if (m == 1) {
        for (jj = 0; jj < n; jj++) {
            pA[jj * bs] = A[jj * lda];
        }
        return;
    }


    m0 = (bs - ai % bs) % bs;
    if (m0 > m)
        m0 = m;
    m1 = m - m0;
    jj = 0;
    for (; jj < n - 3; jj += 4) {
        B = A + jj * lda;
        pB = pA + jj * bs;
        ii = 0;
        if (m0 > 0) {
            for (; ii < m0; ii++) {
                pB[ii + bs * 0] = B[ii + lda * 0];
                pB[ii + bs * 1] = B[ii + lda * 1];
                pB[ii + bs * 2] = B[ii + lda * 2];
                pB[ii + bs * 3] = B[ii + lda * 3];
            }
            B += m0;
            pB += m0 + bs * (sda - 1);
        }
        for (; ii < jj - 3; ii += 4) {
            // col 0
            pB[0 + bs * 0] = B[0 + lda * 0];
            pB[1 + bs * 0] = B[1 + lda * 0];
            pB[2 + bs * 0] = B[2 + lda * 0];
            pB[3 + bs * 0] = B[3 + lda * 0];
            // col 1
            pB[0 + bs * 1] = B[0 + lda * 1];
            pB[1 + bs * 1] = B[1 + lda * 1];
            pB[2 + bs * 1] = B[2 + lda * 1];
            pB[3 + bs * 1] = B[3 + lda * 1];
            // col 2
            pB[0 + bs * 2] = B[0 + lda * 2];
            pB[1 + bs * 2] = B[1 + lda * 2];
            pB[2 + bs * 2] = B[2 + lda * 2];
            pB[3 + bs * 2] = B[3 + lda * 2];
            // col 3
            pB[0 + bs * 3] = B[0 + lda * 3];
            pB[1 + bs * 3] = B[1 + lda * 3];
            pB[2 + bs * 3] = B[2 + lda * 3];
            pB[3 + bs * 3] = B[3 + lda * 3];
            // update
            B += 4;
            pB += bs * sda;
        }
        // col 0
        pB[0 + bs * 0] = B[0 + lda * 0];
        //		pB[1+bs*0] = B[1+lda*0];
        //		pB[2+bs*0] = B[2+lda*0];
        //		pB[3+bs*0] = B[3+lda*0];
        // col 1
        pB[0 + bs * 1] = B[0 + lda * 1];
        pB[1 + bs * 1] = B[1 + lda * 1];
        //		pB[2+bs*1] = B[2+lda*1];
        //		pB[3+bs*1] = B[3+lda*1];
        // col 2
        pB[0 + bs * 2] = B[0 + lda * 2];
        pB[1 + bs * 2] = B[1 + lda * 2];
        pB[2 + bs * 2] = B[2 + lda * 2];
        //		pB[3+bs*2] = B[3+lda*2];
        // col 3
        pB[0 + bs * 3] = B[0 + lda * 3];
        pB[1 + bs * 3] = B[1 + lda * 3];
        pB[2 + bs * 3] = B[2 + lda * 3];
        pB[3 + bs * 3] = B[3 + lda * 3];
    }
    for (; jj < n; jj++) {

        B = A + jj * lda;
        pB = pA + jj * bs;

        ii = 0;
        if (m0 > 0) {
            for (; ii < m0; ii++) {
                pB[ii + bs * 0] = B[ii + lda * 0];
            }
            B += m0;
            pB += m0 + bs * (sda - 1);
        }
        for (; ii < jj - 3; ii += 4) {
            // col 0
            pB[0 + bs * 0] = B[0 + lda * 0];
            pB[1 + bs * 0] = B[1 + lda * 0];
            pB[2 + bs * 0] = B[2 + lda * 0];
            pB[3 + bs * 0] = B[3 + lda * 0];
            // update
            B += 4;
            pB += bs * sda;
        }
        for (; ii <= jj; ii++) {
            // col 0
            pB[0 + bs * 0] = B[0 + lda * 0];
            // update
            B += 1;
            pB += 1;
        }
    }
}


// convert and transpose a matrix into a matrix structure
void pack_tran_mat(int m, int n, double* A, int lda, struct mat* sA, int ai, int aj) {

    // invalidate stored inverse diagonal
    sA->use_dA = 0;

    const int bs = 4;
    int sda = sA->cn;
    double* pA = sA->pA + aj * bs + ai / bs * bs * sda + ai % bs;
    int i, ii, j, m0, m1, m2;
    double *B, *pB;
    sA->use_dA = 0;

    // row vector in sA
    if (n == 1) {
        for (ii = 0; ii < m; ii++) {
            pA[ii * bs] = A[ii];
        }
        return;
    }


    m0 = (bs - ai % bs) % bs;
    if (m0 > n)
        m0 = n;
    m1 = n - m0;
    ii = 0;
    if (m0 > 0) {
        for (j = 0; j < m; j++) {
            for (i = 0; i < m0; i++) {
                pA[i + j * bs + ii * sda] = A[j + (i + ii) * lda];
            }
        }
        A += m0 * lda;
        pA += m0 + bs * (sda - 1);
    }
    ii = 0;
    for (; ii < m1 - 3; ii += bs) {
        j = 0;
        B = A + ii * lda;
        pB = pA + ii * sda;
        for (; j < m - 3; j += 4) {
            // unroll 0
            pB[0 + 0 * bs] = B[0 + 0 * lda];
            pB[1 + 0 * bs] = B[0 + 1 * lda];
            pB[2 + 0 * bs] = B[0 + 2 * lda];
            pB[3 + 0 * bs] = B[0 + 3 * lda];
            // unroll 1
            pB[0 + 1 * bs] = B[1 + 0 * lda];
            pB[1 + 1 * bs] = B[1 + 1 * lda];
            pB[2 + 1 * bs] = B[1 + 2 * lda];
            pB[3 + 1 * bs] = B[1 + 3 * lda];
            // unroll 2
            pB[0 + 2 * bs] = B[2 + 0 * lda];
            pB[1 + 2 * bs] = B[2 + 1 * lda];
            pB[2 + 2 * bs] = B[2 + 2 * lda];
            pB[3 + 2 * bs] = B[2 + 3 * lda];
            // unroll 3
            pB[0 + 3 * bs] = B[3 + 0 * lda];
            pB[1 + 3 * bs] = B[3 + 1 * lda];
            pB[2 + 3 * bs] = B[3 + 2 * lda];
            pB[3 + 3 * bs] = B[3 + 3 * lda];
            B += 4;
            pB += 4 * bs;
        }
        for (; j < m; j++) {
            // unroll 0
            pB[0 + 0 * bs] = B[0 + 0 * lda];
            pB[1 + 0 * bs] = B[0 + 1 * lda];
            pB[2 + 0 * bs] = B[0 + 2 * lda];
            pB[3 + 0 * bs] = B[0 + 3 * lda];
            B += 1;
            pB += 1 * bs;
        }
    }
    if (ii < m1) {
        m2 = m1 - ii;
        if (bs < m2) m2 = bs;
        for (j = 0; j < m; j++) {
            for (i = 0; i < m2; i++) {
                pA[i + j * bs + ii * sda] = A[j + (i + ii) * lda];
            }
        }
    }
}


// convert a vector into a vector structure
void pack_vec(int m, double* x, int xi, struct vec* sa, int ai) {
    double* pa = sa->pa + ai;
    int ii;
    if (xi == 1) {
        for (ii = 0; ii < m; ii++)
            pa[ii] = x[ii];
    } else {
        for (ii = 0; ii < m; ii++)
            pa[ii] = x[ii * xi];
    }
}


// convert a matrix structure into a matrix
void unpack_mat(int m, int n, struct mat* sA, int ai, int aj, double* A, int lda) {
    const int bs = 4;
    int sda = sA->cn;
    double* pA = sA->pA + aj * bs + ai / bs * bs * sda + ai % bs;
    int i, ii, jj;
    int m0 = (bs - ai % bs) % bs;
    if (m0 > m)
        m0 = m;
    double* ptr_pA;


    jj = 0;
    for (; jj < n - 3; jj += 4) {
        ptr_pA = pA + jj * bs;
        ii = 0;
        if (m0 > 0) {
            for (; ii < m0; ii++) {
                // unroll 0
                A[ii + lda * (jj + 0)] = ptr_pA[0 + bs * 0];
                // unroll 1
                A[ii + lda * (jj + 1)] = ptr_pA[0 + bs * 1];
                // unroll 2
                A[ii + lda * (jj + 2)] = ptr_pA[0 + bs * 2];
                // unroll 3
                A[ii + lda * (jj + 3)] = ptr_pA[0 + bs * 3];
                ptr_pA++;
            }
            ptr_pA += (sda - 1) * bs;
        }
        for (; ii < m - bs + 1; ii += bs) {
            // unroll 0
            A[0 + ii + lda * (jj + 0)] = ptr_pA[0 + bs * 0];
            A[1 + ii + lda * (jj + 0)] = ptr_pA[1 + bs * 0];
            A[2 + ii + lda * (jj + 0)] = ptr_pA[2 + bs * 0];
            A[3 + ii + lda * (jj + 0)] = ptr_pA[3 + bs * 0];
            // unroll 0
            A[0 + ii + lda * (jj + 1)] = ptr_pA[0 + bs * 1];
            A[1 + ii + lda * (jj + 1)] = ptr_pA[1 + bs * 1];
            A[2 + ii + lda * (jj + 1)] = ptr_pA[2 + bs * 1];
            A[3 + ii + lda * (jj + 1)] = ptr_pA[3 + bs * 1];
            // unroll 0
            A[0 + ii + lda * (jj + 2)] = ptr_pA[0 + bs * 2];
            A[1 + ii + lda * (jj + 2)] = ptr_pA[1 + bs * 2];
            A[2 + ii + lda * (jj + 2)] = ptr_pA[2 + bs * 2];
            A[3 + ii + lda * (jj + 2)] = ptr_pA[3 + bs * 2];
            // unroll 0
            A[0 + ii + lda * (jj + 3)] = ptr_pA[0 + bs * 3];
            A[1 + ii + lda * (jj + 3)] = ptr_pA[1 + bs * 3];
            A[2 + ii + lda * (jj + 3)] = ptr_pA[2 + bs * 3];
            A[3 + ii + lda * (jj + 3)] = ptr_pA[3 + bs * 3];
            ptr_pA += sda * bs;
        }
        for (; ii < m; ii++) {
            // unroll 0
            A[ii + lda * (jj + 0)] = ptr_pA[0 + bs * 0];
            // unroll 1
            A[ii + lda * (jj + 1)] = ptr_pA[0 + bs * 1];
            // unroll 2
            A[ii + lda * (jj + 2)] = ptr_pA[0 + bs * 2];
            // unroll 3
            A[ii + lda * (jj + 3)] = ptr_pA[0 + bs * 3];
            ptr_pA++;
        }
    }
    for (; jj < n; jj++) {
        ptr_pA = pA + jj * bs;
        ii = 0;
        if (m0 > 0) {
            for (; ii < m0; ii++) {
                A[ii + lda * jj] = ptr_pA[0];
                ptr_pA++;
            }
            ptr_pA += (sda - 1) * bs;
        }
        for (; ii < m - bs + 1; ii += bs) {
            A[0 + ii + lda * jj] = ptr_pA[0];
            A[1 + ii + lda * jj] = ptr_pA[1];
            A[2 + ii + lda * jj] = ptr_pA[2];
            A[3 + ii + lda * jj] = ptr_pA[3];
            ptr_pA += sda * bs;
        }
        for (; ii < m; ii++) {
            A[ii + lda * jj] = ptr_pA[0];
            ptr_pA++;
        }
    }
}


// convert and transpose a matrix structure into a matrix
void unpack_tran_mat(int m, int n, struct mat* sA, int ai, int aj, double* A, int lda) {
    const int bs = 4;
    int sda = sA->cn;
    double* pA = sA->pA + aj * bs + ai / bs * bs * sda + ai % bs;
    int i, ii, jj;
    int m0 = (bs - ai % bs) % bs;
    if (m0 > m)
        m0 = m;
    double* ptr_pA;


    jj = 0;
    for (; jj < n - 3; jj += 4) {
        ptr_pA = pA + jj * bs;
        ii = 0;
        if (m0 > 0) {
            for (; ii < m0; ii++) {
                // unroll 0
                A[jj + 0 + lda * ii] = ptr_pA[0 + bs * 0];
                // unroll 1
                A[jj + 1 + lda * ii] = ptr_pA[0 + bs * 1];
                // unroll 2
                A[jj + 2 + lda * ii] = ptr_pA[0 + bs * 2];
                // unroll 3
                A[jj + 3 + lda * ii] = ptr_pA[0 + bs * 3];
                ptr_pA++;
            }
            ptr_pA += (sda - 1) * bs;
        }
        for (; ii < m - bs + 1; ii += bs) {
            // unroll 0
            A[jj + 0 + lda * (ii + 0)] = ptr_pA[0 + bs * 0];
            A[jj + 0 + lda * (ii + 1)] = ptr_pA[1 + bs * 0];
            A[jj + 0 + lda * (ii + 2)] = ptr_pA[2 + bs * 0];
            A[jj + 0 + lda * (ii + 3)] = ptr_pA[3 + bs * 0];
            // unroll 1
            A[jj + 1 + lda * (ii + 0)] = ptr_pA[0 + bs * 1];
            A[jj + 1 + lda * (ii + 1)] = ptr_pA[1 + bs * 1];
            A[jj + 1 + lda * (ii + 2)] = ptr_pA[2 + bs * 1];
            A[jj + 1 + lda * (ii + 3)] = ptr_pA[3 + bs * 1];
            // unroll 2
            A[jj + 2 + lda * (ii + 0)] = ptr_pA[0 + bs * 2];
            A[jj + 2 + lda * (ii + 1)] = ptr_pA[1 + bs * 2];
            A[jj + 2 + lda * (ii + 2)] = ptr_pA[2 + bs * 2];
            A[jj + 2 + lda * (ii + 3)] = ptr_pA[3 + bs * 2];
            // unroll 3
            A[jj + 3 + lda * (ii + 0)] = ptr_pA[0 + bs * 3];
            A[jj + 3 + lda * (ii + 1)] = ptr_pA[1 + bs * 3];
            A[jj + 3 + lda * (ii + 2)] = ptr_pA[2 + bs * 3];
            A[jj + 3 + lda * (ii + 3)] = ptr_pA[3 + bs * 3];
            ptr_pA += sda * bs;
        }
        for (; ii < m; ii++) {
            // unroll 0
            A[jj + 0 + lda * ii] = ptr_pA[0 + bs * 0];
            // unroll 1
            A[jj + 1 + lda * ii] = ptr_pA[0 + bs * 1];
            // unroll 2
            A[jj + 2 + lda * ii] = ptr_pA[0 + bs * 2];
            // unroll 3
            A[jj + 3 + lda * ii] = ptr_pA[0 + bs * 3];
            ptr_pA++;
        }
    }
    for (; jj < n; jj++) {
        ptr_pA = pA + jj * bs;
        ii = 0;
        if (m0 > 0) {
            for (; ii < m0; ii++) {
                A[jj + lda * ii] = ptr_pA[0];
                ptr_pA++;
            }
            ptr_pA += (sda - 1) * bs;
        }
        for (; ii < m - bs + 1; ii += bs) {
            i = 0;
            for (; i < bs; i++) {
                A[jj + lda * (i + ii)] = ptr_pA[0];
                ptr_pA++;
            }
            ptr_pA += (sda - 1) * bs;
        }
        for (; ii < m; ii++) {
            A[jj + lda * ii] = ptr_pA[0];
            ptr_pA++;
        }
    }
}


// convert a vector structure into a vector
void unpack_vec(int m, struct vec* sa, int ai, double* x, int xi) {
    double* pa = sa->pa + ai;
    int ii;
    if (xi == 1) {
        for (ii = 0; ii < m; ii++)
            x[ii] = pa[ii];
    } else {
        for (ii = 0; ii < m; ii++)
            x[ii * xi] = pa[ii];
    }
}
