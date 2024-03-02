#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void naive_dgemv_n(int m, int n, double* A, int lda, double* x, double* z) {
    for (int ii = 0; ii < m; ii++) {
        z[ii] = 0.0;
        for (int jj = 0; jj < n; jj++) {
            z[ii] += A[ii + lda * jj] * x[jj];
        }
    }
}


void naive_dgemm_nn(int m, int n, int k, double* A, int lda, double* B, int ldb, double* C, int ldc) {
    for (int jj = 0; jj < n; jj++) {
        for (int ii = 0; ii < m; ii++) {
            C[ii + ldc * jj] = 0;
            for (int kk = 0; kk < k; kk++) {
                C[ii + ldc * jj] += A[ii + lda * kk] * B[kk + ldb * jj];
            }
        }
    }
}


void naive_daxpy(int n, double da, double* dx, double* dy) {
    for (int i = 0; i < n; i++) {
        dy[i] += da * dx[i];
    }
}


void naive_dscal(int n, double da, double* dx) {
    for (int i = 0; i < n; i++) {
        dx[i] *= da;
    }
}

void naive_dzeros(int m, int n, double* A) {
    for (int i; i < m * n; i++) {
        A[i] = 0.0;
    }
}


void naive_dmcopy(int row, int col, double* A, int lda, double* B, int ldb) {
    int i, j;
    for (j = 0; j < col; j++) {
        for (i = 0; i < row; i++) {
            B[i + j * ldb] = A[i + j * lda];
        }
    }
}


int idamax_3l(int n, double* x) {
    if (n <= 0) return 0;
    if (n == 1) return 0;
    double dabs;
    double dmax = (x[0] > 0 ? x[0] : -x[0]);
    int idmax = 0;
    int jj;
    for (jj = 1; jj < n; jj++) {
        dabs = (x[jj] > 0 ? x[jj] : -x[jj]);
        if (dabs > dmax) {
            dmax = dabs;
            idmax = jj;
        }
    }
    return idmax;
}


void dswap_3l(int n, double* x, int incx, double* y, int incy) {
    if (n <= 0) return;
    double temp;
    int jj;
    for (jj = 0; jj < n; jj++) {
        temp = x[0];
        x[0] = y[0];
        y[0] = temp;
        x += incx;
        y += incy;
    }
}


void dger_3l(int m, int n, double alpha, double* x, int incx, double* y, int incy, double* A, int lda) {
    if (m == 0 || n == 0 || alpha == 0.0) return;
    int i, j;
    double *px, *py, temp;
    py = y;
    for (j = 0; j < n; j++) {
        temp = alpha * py[0];
        px = x;
        for (i = 0; i < m; i++) {
            A[i + lda * j] += px[0] * temp;
            px += incx;
        }
        py += incy;
    }
}


void dgetf2_3l(int m, int n, double* A, int lda, int* ipiv, int* info) {
    if (m <= 0 || n <= 0) return;
    int i, j, jp;
    double Ajj;
    int size_min = (m < n ? m : n);
    for (j = 0; j < size_min; j++) {  // find the pivot and test for singularity
        jp = j + idamax_3l(m - j, &A[j + lda * j]);
        ipiv[j] = jp;
        if (A[jp + lda * j] != 0) {
            // apply the interchange to columns 0:n-1
            if (jp != j) {
                dswap_3l(n, &A[j], lda, &A[jp], lda);
            }
            // compute elements j+1:m-1 of j-th column
            if (j < m - 1) {
                Ajj = A[j + lda * j];
                if ((Ajj > 0 ? Ajj : -Ajj) >= 2.22e-16) {
                    naive_dscal(m - j - 1, 1.0 / Ajj, &A[j + 1 + lda * j]);
                } else {
                    for (i = j + 1; i < m; i++) {
                        A[i + lda * j] /= Ajj;
                    }
                }
            }
        } else if (*info == 0) {
            *info = j + 1;
        }

        if (j < size_min) {
            // update trailing submatrix
            dger_3l(m - j - 1, n - j - 1, -1.0, &A[j + 1 + lda * j], 1, &A[j + lda * (j + 1)], lda, &A[j + 1 + lda * (j + 1)], lda);
        }
    }
}


void dlaswp_3l(int n, double* A, int lda, int k1, int k2, int* ipiv) {
    int i, j, k, ix, ix0, i1, i2, n32, ip;
    double temp;
    ix0 = k1;
    i1 = k1;
    i2 = k2;
    n32 = (n / 32) * 32;
    if (n32 != 0) {
        for (j = 0; j < n32; j += 32) {
            ix = ix0;
            for (i = i1; i < i2; i++) {
                ip = ipiv[ix];
                if (ip != i) {
                    for (k = j; k < j + 32; k++) {
                        temp = A[i + lda * k];
                        A[i + lda * k] = A[ip + lda * k];
                        A[ip + lda * k] = temp;
                    }
                }
                ix++;
            }
        }
    }
    if (n32 != n) {
        ix = ix0;
        for (i = i1; i < i2; i++) {
            ip = ipiv[ix];
            if (ip != i) {
                for (k = n32; k < n; k++) {
                    temp = A[i + lda * k];
                    A[i + lda * k] = A[ip + lda * k];
                    A[ip + lda * k] = temp;
                }
            }
            ix++;
        }
    }
}


// left lower no-transp unit
void dtrsm_l_l_n_u_3l(int m, int n, double* A, int lda, double* B, int ldb) {
    if (m == 0 || n == 0) return;
    int i, j, k;
    for (j = 0; j < n; j++) {
        for (k = 0; k < m; k++) {
            for (i = k + 1; i < m; i++) {
                B[i + ldb * j] -= B[k + ldb * j] * A[i + lda * k];
            }
        }
    }
}


// left upper no-transp non-unit
void dtrsm_l_u_n_n_3l(int m, int n, double* A, int lda, double* B, int ldb) {
    if (m == 0 || n == 0) return;
    int i, j, k;
    for (j = 0; j < n; j++) {
        for (k = m - 1; k >= 0; k--) {
            B[k + ldb * j] /= A[k + lda * k];
            for (i = 0; i < k; i++) {
                B[i + ldb * j] -= B[k + ldb * j] * A[i + lda * k];
            }
        }
    }
}


void dgetrs_3l(int n, int nrhs, double* A, int lda, int* ipiv, double* B, int ldb, int* info) {
    if (n == 0 || nrhs == 0) return;
    // solve A * X = B
    // apply row interchanges to the rhs
    dlaswp_3l(nrhs, B, ldb, 0, n, ipiv);
    // solve L*X = B, overwriting B with X
    dtrsm_l_l_n_u_3l(n, nrhs, A, lda, B, ldb);
    // solve U*X = B, overwriting B with X
    dtrsm_l_u_n_n_3l(n, nrhs, A, lda, B, ldb);
}


void naive_dgesv(int n, int nrhs, double* A, int lda, int* ipiv, double* B, int ldb, int* info) {
    *info = 0;
    // compute the LU factorization of A
    dgetf2_3l(n, n, A, lda, ipiv, info);
    if (*info == 0) {
        // solve the system A*X = B, overwriting B with X
        dgetrs_3l(n, nrhs, A, lda, ipiv, B, ldb, info);
    }
}


/* one norm of a matrix */
double onenorm(int row, int col, double* ptrA) {
    double max, temp;
    int i, j;
    temp = 0;
    for (j = 0; j < col; j++) {
        temp = fabs(*(ptrA + j * row));
        for (i = 1; i < row; i++) {
            temp += fabs(*(ptrA + j * row + i));
        }
        if (j == 0) max = temp;
        else if (max > temp)
            temp = max;
    }
    return temp;
}


// computes the Pade approximation of degree m of the matrix A, used in expm
void padeapprox(int m, int row, double* A) {
    int row2 = row * row;
    /*	int i1 = 1;*/
    /*	double d0 = 0;*/
    /*	double d1 = 1;*/
    /*	double dm1 = -1;*/
    int ii;
    double* U = malloc(row * row * sizeof(double));
    double* V = malloc(row * row * sizeof(double));
    if (m == 3) {
        double c[] = {120, 60, 12, 1};
        double* A0 = malloc(row * row * sizeof(double));
        for (ii = 0; ii < row2; ii++) A0[ii] = 0.0;
        for (ii = 0; ii < row; ii++) A0[ii * (row + 1)] = 1.0;
        double* A2 = malloc(row * row * sizeof(double));
        //		char ta = 'n'; double alpha = 1; double beta = 0;
        //		dgemm_(&ta, &ta, &row, &row, &row, &alpha, A, &row, A, &row, &beta, A2, &row);
        naive_dgemm_nn(row, row, row, A, row, A, row, A2, row);
        double* temp = malloc(row * row * sizeof(double));
        //		dscal_(&row2, &d0, temp, &i1);
        naive_dscal(row2, 0, temp);
        //		daxpy_(&row2, &c[3], A2, &i1, temp, &i1);
        naive_daxpy(row2, c[3], A2, temp);
        //		daxpy_(&row2, &c[1], A0, &i1, temp, &i1);
        naive_daxpy(row2, c[1], A0, temp);
        //		dgemm_(&ta, &ta, &row, &row, &row, &alpha, A, &row, temp, &row, &beta, U, &row);
        naive_dgemm_nn(row, row, row, A, row, temp, row, U, row);
        //		dscal_(&row2, &d0, V, &i1);
        naive_dscal(row2, 0, V);
        //		daxpy_(&row2, &c[2], A2, &i1, V, &i1);
        naive_daxpy(row2, c[2], A2, V);
        //		daxpy_(&row2, &c[0], A0, &i1, V, &i1);
        naive_daxpy(row2, c[0], A0, V);
        free(A0);
        free(A2);
        free(temp);
    } else if (m == 5) {
        double c[] = {30240, 15120, 3360, 420, 30, 1};
        double* A0 = malloc(row * row * sizeof(double));
        for (ii = 0; ii < row2; ii++) A0[ii] = 0.0;
        for (ii = 0; ii < row; ii++) A0[ii * (row + 1)] = 1.0;
        double* A2 = malloc(row * row * sizeof(double));
        double* A4 = malloc(row * row * sizeof(double));
        //		char ta = 'n'; double alpha = 1; double beta = 0;
        //		dgemm_(&ta, &ta, &row, &row, &row, &alpha, A, &row, A, &row, &beta, A2, &row);
        naive_dgemm_nn(row, row, row, A, row, A, row, A2, row);
        //		dgemm_(&ta, &ta, &row, &row, &row, &alpha, A2, &row, A2, &row, &beta, A4, &row);
        naive_dgemm_nn(row, row, row, A2, row, A2, row, A4, row);
        naive_dmcopy(row, row, A4, row, V, row);
        double* temp = malloc(row * row * sizeof(double));
        naive_dmcopy(row, row, A4, row, temp, row);
        //		daxpy_(&row2, &c[3], A2, &i1, temp, &i1);
        naive_daxpy(row2, c[3], A2, temp);
        //		daxpy_(&row2, &c[1], A0, &i1, temp, &i1);
        naive_daxpy(row2, c[1], A0, temp);
        //		dgemm_(&ta, &ta, &row, &row, &row, &alpha, A, &row, temp, &row, &beta, U, &row);
        naive_dgemm_nn(row, row, row, A, row, temp, row, U, row);
        //		dscal_(&row2, &c[4], V, &i1);
        naive_dscal(row2, c[4], V);
        //		daxpy_(&row2, &c[2], A2, &i1, V, &i1);
        naive_daxpy(row2, c[2], A2, V);
        //		daxpy_(&row2, &c[0], A0, &i1, V, &i1);
        naive_daxpy(row2, c[0], A0, V);
        free(A0);
        free(A2);
        free(A4);
        free(temp);
    } else if (m == 7) {
        double c[] = {17297280, 8648640, 1995840, 277200, 25200, 1512, 56, 1};
        double* A0 = malloc(row * row * sizeof(double));
        for (ii = 0; ii < row2; ii++) A0[ii] = 0.0;
        for (ii = 0; ii < row; ii++) A0[ii * (row + 1)] = 1.0;
        double* A2 = malloc(row * row * sizeof(double));
        double* A4 = malloc(row * row * sizeof(double));
        double* A6 = malloc(row * row * sizeof(double));
        //		char ta = 'n'; double alpha = 1; double beta = 1;
        //		dgemm_(&ta, &ta, &row, &row, &row, &alpha, A, &row, A, &row, &beta, A2, &row);
        naive_dgemm_nn(row, row, row, A, row, A, row, A2, row);
        //		dgemm_(&ta, &ta, &row, &row, &row, &alpha, A2, &row, A2, &row, &beta, A4, &row);
        naive_dgemm_nn(row, row, row, A2, row, A2, row, A4, row);
        //		dgemm_(&ta, &ta, &row, &row, &row, &alpha, A4, &row, A2, &row, &beta, A6, &row);
        naive_dgemm_nn(row, row, row, A4, row, A2, row, A6, row);
        double* temp = malloc(row * row * sizeof(double));
        //		dscal_(&row2, &d0, temp, &i1);
        naive_dscal(row2, 0, temp);
        //		daxpy_(&row2, &c[3], A2, &i1, temp, &i1);
        naive_daxpy(row2, c[3], A2, temp);
        //		daxpy_(&row2, &c[1], A0, &i1, temp, &i1);
        naive_daxpy(row2, c[1], A0, temp);
        //		daxpy_(&row2, &c[5], A4, &i1, temp, &i1);
        naive_daxpy(row2, c[5], A4, temp);
        //		daxpy_(&row2, &c[7], A6, &i1, temp, &i1);
        naive_daxpy(row2, c[7], A6, temp);
        //		dgemm_(&ta, &ta, &row, &row, &row, &alpha, A, &row, temp, &row, &beta, U, &row);
        naive_dgemm_nn(row, row, row, A, row, temp, row, U, row);
        //		dscal_(&row2, &d0, V, &i1);
        naive_dscal(row2, 0, V);
        //		daxpy_(&row2, &c[2], A2, &i1, V, &i1);
        naive_daxpy(row2, c[2], A2, V);
        //		daxpy_(&row2, &c[0], A0, &i1, V, &i1);
        naive_daxpy(row2, c[0], A0, V);
        //		daxpy_(&row2, &c[4], A4, &i1, V, &i1);
        naive_daxpy(row2, c[4], A4, V);
        //		daxpy_(&row2, &c[6], A6, &i1, V, &i1);
        naive_daxpy(row2, c[6], A6, V);
        free(A0);
        free(A2);
        free(A4);
        free(A6);
        free(temp);
    } else if (m == 9) {
        double c[] = {17643225600, 8821612800, 2075673600, 302702400, 30270240, 2162160, 110880, 3960, 90, 1};
        double* A0 = malloc(row * row * sizeof(double));
        for (ii = 0; ii < row2; ii++) A0[ii] = 0.0;
        for (ii = 0; ii < row; ii++) A0[ii * (row + 1)] = 1.0;
        double* A2 = malloc(row * row * sizeof(double));
        double* A4 = malloc(row * row * sizeof(double));
        double* A6 = malloc(row * row * sizeof(double));
        double* A8 = malloc(row * row * sizeof(double));
        //		char ta = 'n'; double alpha = 1; double beta = 0;
        //		dgemm_(&ta, &ta, &row, &row, &row, &alpha, A, &row, A, &row, &beta, A2, &row);
        naive_dgemm_nn(row, row, row, A, row, A, row, A2, row);
        //		dgemm_(&ta, &ta, &row, &row, &row, &alpha, A2, &row, A2, &row, &beta, A4, &row);
        naive_dgemm_nn(row, row, row, A2, row, A2, row, A4, row);
        //		dgemm_(&ta, &ta, &row, &row, &row, &alpha, A4, &row, A2, &row, &beta, A6, &row);
        naive_dgemm_nn(row, row, row, A4, row, A2, row, A6, row);
        //		dgemm_(&ta, &ta, &row, &row, &row, &alpha, A6, &row, A2, &row, &beta, A8, &row);
        naive_dgemm_nn(row, row, row, A6, row, A2, row, A8, row);
        naive_dmcopy(row, row, A8, row, V, row);
        double* temp = malloc(row * row * sizeof(double));
        naive_dmcopy(row, row, A8, row, temp, row);
        //		daxpy_(&row2, &c[3], A2, &i1, temp, &i1);
        naive_daxpy(row2, c[3], A2, temp);
        //		daxpy_(&row2, &c[1], A0, &i1, temp, &i1);
        naive_daxpy(row2, c[1], A0, temp);
        //		daxpy_(&row2, &c[5], A4, &i1, temp, &i1);
        naive_daxpy(row2, c[5], A4, temp);
        //		daxpy_(&row2, &c[7], A6, &i1, temp, &i1);
        naive_daxpy(row2, c[7], A6, temp);
        //		dgemm_(&ta, &ta, &row, &row, &row, &alpha, A, &row, temp, &row, &beta, U, &row);
        naive_dgemm_nn(row, row, row, A, row, temp, row, U, row);
        //		dscal_(&row2, &c[8], V, &i1);
        naive_dscal(row2, c[8], V);
        //		daxpy_(&row2, &c[2], A2, &i1, V, &i1);
        naive_daxpy(row2, c[2], A2, V);
        //		daxpy_(&row2, &c[0], A0, &i1, V, &i1);
        naive_daxpy(row2, c[0], A0, V);
        //		daxpy_(&row2, &c[4], A4, &i1, V, &i1);
        naive_daxpy(row2, c[4], A4, V);
        //		daxpy_(&row2, &c[6], A6, &i1, V, &i1);
        naive_daxpy(row2, c[6], A6, V);
        free(A0);
        free(A2);
        free(A4);
        free(A6);
        free(A8);
        free(temp);
    } else if (m == 13)  // tested
    {
        double c[] = {64764752532480000, 32382376266240000, 7771770303897600, 1187353796428800, 129060195264000, 10559470521600, 670442572800, 33522128640, 1323241920, 40840800, 960960, 16380, 182, 1};
        double* A0 = malloc(row * row * sizeof(double));
        for (ii = 0; ii < row2; ii++) A0[ii] = 0.0;
        for (ii = 0; ii < row; ii++) A0[ii * (row + 1)] = 1.0;
        double* A2 = malloc(row * row * sizeof(double));
        double* A4 = malloc(row * row * sizeof(double));
        double* A6 = malloc(row * row * sizeof(double));
        //		char ta = 'n'; double alpha = 1; double beta = 0;
        //		dgemm_(&ta, &ta, &row, &row, &row, &alpha, A, &row, A, &row, &beta, A2, &row);
        naive_dgemm_nn(row, row, row, A, row, A, row, A2, row);
        //		dgemm_(&ta, &ta, &row, &row, &row, &alpha, A2, &row, A2, &row, &beta, A4, &row);
        naive_dgemm_nn(row, row, row, A2, row, A2, row, A4, row);
        //		dgemm_(&ta, &ta, &row, &row, &row, &alpha, A4, &row, A2, &row, &beta, A6, &row);
        naive_dgemm_nn(row, row, row, A4, row, A2, row, A6, row);
        naive_dmcopy(row, row, A2, row, U, row);
        double* temp = malloc(row * row * sizeof(double));
        //		dscal_(&row2, &c[9], U, &i1);
        naive_dscal(row2, c[9], U);
        //		daxpy_(&row2, &c[11], A4, &i1, U, &i1);
        naive_daxpy(row2, c[11], A4, U);
        //		daxpy_(&row2, &c[13], A6, &i1, U, &i1);
        naive_daxpy(row2, c[13], A6, U);
        //		dgemm_(&ta, &ta, &row, &row, &row, &alpha, A6, &row, U, &row, &beta, temp, &row);
        naive_dgemm_nn(row, row, row, A6, row, U, row, temp, row);
        //		daxpy_(&row2, &c[7], A6, &i1, temp, &i1);
        naive_daxpy(row2, c[7], A6, temp);
        //		daxpy_(&row2, &c[5], A4, &i1, temp, &i1);
        naive_daxpy(row2, c[5], A4, temp);
        //		daxpy_(&row2, &c[3], A2, &i1, temp, &i1);
        naive_daxpy(row2, c[3], A2, temp);
        //		daxpy_(&row2, &c[1], A0, &i1, temp, &i1);
        naive_daxpy(row2, c[1], A0, temp);
        //		dgemm_(&ta, &ta, &row, &row, &row, &alpha, A, &row, temp, &row, &beta, U, &row);
        naive_dgemm_nn(row, row, row, A, row, temp, row, U, row);
        naive_dmcopy(row, row, A2, row, temp, row);
        //		dscal_(&row2, &c[8], V, &i1);
        naive_dscal(row2, c[8], V);
        //		daxpy_(&row2, &c[12], A6, &i1, temp, &i1);
        naive_daxpy(row2, c[12], A6, temp);
        //		daxpy_(&row2, &c[10], A4, &i1, temp, &i1);
        naive_daxpy(row2, c[10], A4, temp);
        //		dgemm_(&ta, &ta, &row, &row, &row, &alpha, A6, &row, temp, &row, &beta, V, &row);
        naive_dgemm_nn(row, row, row, A6, row, temp, row, V, row);
        //		daxpy_(&row2, &c[6], A6, &i1, V, &i1);
        naive_daxpy(row2, c[6], A6, V);
        //		daxpy_(&row2, &c[4], A4, &i1, V, &i1);
        naive_daxpy(row2, c[4], A4, V);
        //		daxpy_(&row2, &c[2], A2, &i1, V, &i1);
        naive_daxpy(row2, c[2], A2, V);
        //		daxpy_(&row2, &c[0], A0, &i1, V, &i1);
        naive_daxpy(row2, c[0], A0, V);
        free(A0);
        free(A2);
        free(A4);
        free(A6);
        free(temp);
    } else {
        printf("%s\n", "Wrong Pade approximatin degree");
        exit(1);
    }
    double* D = malloc(row * row * sizeof(double));
    //	dcopy_(&row2, V, &i1, A, &i1);
    naive_dmcopy(row, row, V, row, A, row);
    //	daxpy_(&row2, &d1, U, &i1, A, &i1);
    naive_daxpy(row2, 1.0, U, A);
    //	dcopy_(&row2, V, &i1, D, &i1);
    naive_dmcopy(row, row, V, row, D, row);
    //	daxpy_(&row2, &dm1, U, &i1, D, &i1);
    naive_daxpy(row2, -1.0, U, D);
    int* ipiv = (int*) malloc(row * sizeof(int));
    int info = 0;
    //	dgesv_(&row, &row, D, &row, ipiv, A, &row, &info);
    naive_dgesv(row, row, D, row, ipiv, A, row, &info);
    free(ipiv);
    free(D);
    free(U);
    free(V);
}


void expm(int row, double* A) {
    int i;
    int m_vals[] = {3, 5, 7, 9, 13};
    double theta[] = {0.01495585217958292, 0.2539398330063230, 0.9504178996162932, 2.097847961257068, 5.371920351148152};
    int lentheta = 5;
    double normA = onenorm(row, row, A);
    if (normA <= theta[4]) {
        for (i = 0; i < lentheta; i++) {
            if (normA <= theta[i]) {
                padeapprox(m_vals[i], row, A);
                break;
            }
        }
    } else {
        int s;
        double t = frexp(normA / (theta[4]), &s);
        s = s - (t == 0.5);
        t = pow(2, -s);
        int row2 = row * row;
        /*		int i1 = 1;*/
        //		dscal_(&row2, &t, A, &i1);
        naive_dscal(row2, t, A);
        padeapprox(m_vals[4], row, A);
        double* temp = malloc(row * row * sizeof(double));
        //		char ta = 'n'; double alpha = 1; double beta = 0;
        for (i = 0; i < s; i++) {
            //			dgemm_(&ta, &ta, &row, &row, &row, &alpha, A, &row, A, &row, &beta, temp, &row);
            naive_dgemm_nn(row, row, row, A, row, A, row, temp, row);
            naive_dmcopy(row, row, temp, row, A, row);
        }
        free(temp);
    }
}
