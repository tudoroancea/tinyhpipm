#ifndef TINYHPIPM_TEST_UTILS_NAIVE_BLAS_H
#define TINYHPIPM_TEST_UTILS_NAIVE_BLAS_H


/*
 * @brief Linear combination z = a * x + y
 *
 * @param[in] n number of elements of x and y
 * @param[in] a scalar a
 * @param[in] x vector x
 * @param[in] y vector y
 * @param[out] dz vector z
 */
void naive_daxpy(int n, double a, const double* x, double* y);

/*
 * @brief Scaling of a vector x by a scalar a
 *
 * @param[in] n number of elements of x
 * @param[in] da scalar a
 * @param[in] dx vector x
 * @param[out] dy vector y
 */
void naive_dscal(int n, double da, double* dx);

/*
 * @brief Sets all elements of a matrix to zero
 *
 * @param[in] m number of rows of A
 * @param[in] n number of columns of A
 * @param[out] A matrix A, row-major order
 */
void naive_dzeros(int m, int n, double* A);

/*
 * @brief Copies a matrix A to B
 *
 * @param[in] row number of rows of A and B
 * @param[in] col number of columns of A and B
 * @param[in] A matrix A, row-major order
 * @param[in] lda leading dimension of A
 * @param[out] B matrix B, row-major order
 * @param[in] ldb leading dimension of B
 */
void naive_dmcopy(int row, int col, const double* A, int lda, double* B, int ldb);

/*
 * @brief General matrix-vector multiplication z = A * x, where A is not transposed.
 *
 * @param[in] m number of rows of A
 * @param[in] n number of columns of A
 * @param[in] A matrix A, not transposed and in row-major order
 * @param[in] lda leading dimension of A
 * @param[in] x vector x
 * @param[out] z vector z
 */
void naive_dgemv_n(int m, int n, const double* A, int lda, const double* x, double* z);

/*
 * @brief General matrix-matrix multiplication C = A * B, where both A and B are not transposed.
 *
 * @param[in] m number of rows of A
 * @param[in] n number of columns of A
 * @param[in] k number of columns of B
 * @param[in] A matrix A, not transposed and in row-major order
 * @param[in] lda leading dimension of A
 * @param[in] B matrix B, not transposed and in row-major order
 * @param[in] ldb leading dimension of B
 * @param[out] C matrix C, not transposed and in row-major order
 * @param[in] ldc leading dimension of C
 */
void naive_dgemm_nn(int m, int n, int k, const double* A, int lda, const double* B, int ldb, double* C, int ldc);

/*
 * @brief Solves a linear system of equations A * X = B
 *
 * @param[in] n number of rows and columns of A and B
 * @param[in] nrhs number of columns of B
 * @param[in] A matrix A
 * @param[in] lda leading dimension of A
 * @param[out] ipiv pivot indices
 * @param[in] B matrix B
 * @param[in] ldb leading dimension of B
 * @param[out] info status of the computation
 */
void naive_dgesv(int n, int nrhs, double* A, int lda, int* ipiv, double* B, int ldb, int* info);

/*
 * @brief Computes the matrix exponential of a matrix A
 *
 * @param[in] row number of rows of A
 * @param[in, out] A matrix A
 */
void expm(int row, double* A);

#endif  // TINYHPIPM_TEST_UTILS_NAIVE_BLAS_H
