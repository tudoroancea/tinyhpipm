#ifndef TINYHPIPM_TEST_UTILS_NAIVE_BLAS_H
#define TINYHPIPM_TEST_UTILS_NAIVE_BLAS_H

void naive_dgemv_n(int m, int n, double* A, int lda, double* x, double* z);
void naive_dgemm_nn(int m, int n, int k, double* A, int lda, double* B, int ldb, double* C, int ldc);
void naive_daxpy(int n, double da, double* dx, double* dy);
void naive_dscal(int n, double da, double* dx);

void naive_dzeros(int m, int n, double* A);

void naive_dmcopy(int row, int col, double* A, int lda, double* B, int ldb);

void naive_dgesv(int n, int nrhs, double* A, int lda, int* ipiv, double* B, int ldb, int* info);

void expm(int row, double* A);

#endif  // TINYHPIPM_TEST_UTILS_NAIVE_BLAS_H
