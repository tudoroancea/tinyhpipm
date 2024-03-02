#ifndef TINYHPIPM_BLAS_BLAS2_H
#define TINYHPIPM_BLAS_BLAS2_H

#include "tinyhpipm/blas/struct.h"

// z <= beta * y + alpha * A * x
void dgemv_n(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct vec* sx, int xi, double beta, struct vec* sy, int yi, struct vec* sz, int zi);  // to keep
// z <= beta * y + alpha * A^T * x
void dgemv_t(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct vec* sx, int xi, double beta, struct vec* sy, int yi, struct vec* sz, int zi);  // to keep
// z <= beta * y + alpha * A * x, A diagonal
void dgemv_d(int m, double alpha, struct vec* sA, int ai, struct vec* sx, int xi, double beta, struct vec* sy, int yi, struct vec* sz, int zi);  // to keep
// z <= inv( A ) * x, A (m)x(n)
void dtrsv_lnn_mn(int m, int n, struct mat* sA, int ai, int aj, struct vec* sx, int xi, struct vec* sz, int zi);  // to keep
// z <= inv( A^T ) * x, A (m)x(n)
void dtrsv_ltn_mn(int m, int n, struct mat* sA, int ai, int aj, struct vec* sx, int xi, struct vec* sz, int zi);  // to keep
// z <= inv( A ) * x, A (m)x(m) lower, not_transposed
void dtrsv_lnn(int m, struct mat* sA, int ai, int aj, struct vec* sx, int xi, struct vec* sz, int zi);  // to keep
// z <= inv( A ) * x, A (m)x(m) lower, not_transposed, assuming unit diagonal
// void dtrsv_lnu(int m, struct mat* sA, int ai, int aj, struct vec* sx, int xi, struct vec* sz, int zi);
// z <= inv( A^T ) * x, A (m)x(m) lower, transposed
void dtrsv_ltn(int m, struct mat* sA, int ai, int aj, struct vec* sx, int xi, struct vec* sz, int zi);  // to keep
// z <= inv( A^T ) * x, A (m)x(m) lower, transposed, assuming unit diagonal
// void dtrsv_ltu(int m, struct mat* sA, int ai, int aj, struct vec* sx, int xi, struct vec* sz, int zi);
// z <= inv( A^T ) * x, A (m)x(m) upper, not_transposed
// void dtrsv_unn(int m, struct mat* sA, int ai, int aj, struct vec* sx, int xi, struct vec* sz, int zi);
// z <= inv( A^T ) * x, A (m)x(m) upper, transposed
// void dtrsv_utn(int m, struct mat* sA, int ai, int aj, struct vec* sx, int xi, struct vec* sz, int zi);
// z <= A * x ; A lower triangular
void dtrmv_lnn(int m, struct mat* sA, int ai, int aj, struct vec* sx, int xi, struct vec* sz, int zi);  // to keep
// z <= A * x ; A lower triangular, assuming unit diagonal
// void dtrmv_lnu(int m, struct mat* sA, int ai, int aj, struct vec* sx, int xi, struct vec* sz, int zi);
// z <= A^T * x ; A lower triangular
void dtrmv_ltn(int m, struct mat* sA, int ai, int aj, struct vec* sx, int xi, struct vec* sz, int zi);  // to keep
// z <= A^T * x ; A lower triangular, assuming unit diagonal
// void dtrmv_ltu(int m, struct mat* sA, int ai, int aj, struct vec* sx, int xi, struct vec* sz, int zi);
// z <= beta * y + alpha * A * x ; A upper triangular
// void dtrmv_unn(int m, struct mat* sA, int ai, int aj, struct vec* sx, int xi, struct vec* sz, int zi);
// z <= A^T * x ; A upper triangular
// void dtrmv_utn(int m, struct mat* sA, int ai, int aj, struct vec* sx, int xi, struct vec* sz, int zi);
// z_n <= beta_n * y_n + alpha_n * A  * x_n
// z_t <= beta_t * y_t + alpha_t * A^T * x_t
void dgemv_nt(int m, int n, double alpha_n, double alpha_t, struct mat* sA, int ai, int aj, struct vec* sx_n, int xi_n, struct vec* sx_t, int xi_t, double beta_n, double beta_t, struct vec* sy_n, int yi_n, struct vec* sy_t, int yi_t, struct vec* sz_n, int zi_n, struct vec* sz_t, int zi_t);  // to keep
// z <= beta * y + alpha * A * x, where A is symmetric and only the lower triangular patr of A is accessed
void dsymv_l(int m, double alpha, struct mat* sA, int ai, int aj, struct vec* sx, int xi, double beta, struct vec* sy, int yi, struct vec* sz, int zi);  // to keep
// void dsymv_l_mn(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct vec* sx, int xi, double beta, struct vec* sy, int yi, struct vec* sz, int zi);


#endif  // TINYHPIPM_BLAS_BLAS2_H
