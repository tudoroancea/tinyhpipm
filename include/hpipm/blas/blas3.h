#ifndef TINYHPIPM_BLAS_BLAS3_H
#define TINYHPIPM_BLAS_BLAS3_H

#include "hpipm/blas/struct.h"

/************************************************
 * general matrix-matrix multiplication
 ************************************************/
// D <= beta * C + alpha * A * B
void dgemm_nn(int m, int n, int k, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, double beta, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj);  // to keep
// D <= beta * C + alpha * A * B^T
void dgemm_nt(int m, int n, int k, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, double beta, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj);  // to keep
// D <= beta * C + alpha * A^T * B
// void dgemm_tn(int m, int n, int k, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, double beta, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj);
// D <= beta * C + alpha * A^T * B^T
// void dgemm_tt(int m, int n, int k, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, double beta, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj);


/*****************************************************
 * diagonal matrix-matrix multiplication
 * ***************************************************/
// D <= alpha * A * B + beta * C, with A diagonal (stored as strvec)
// void dgemm_diag_left_lib(int m, int n, double alpha, double* dA, double* pB, int sdb, double beta, double* pC, int sdc, double* pD, int sdd);
void dgemm_dn(int m, int n, double alpha, struct vec* sA, int ai, struct mat* sB, int bi, int bj, double beta, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj);  // to keep
// D <= alpha * A * B + beta * C, with B diagonal (stored as strvec)
void dgemm_nd(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct vec* sB, int bi, double beta, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj);  // to keep


/*****************************************************
 * triangular matrix-matrix multiplication
 * ***************************************************/
// D <= alpha * A * B ; A lower triangular
// void dtrmm_llnn(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);
// D <= alpha * A * B ; A lower triangular assuming unit diagonal
// void dtrmm_llnu(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);
// D <= alpha * A^T * B ; A lower triangular
// void dtrmm_lltn(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);
// D <= alpha * A^T * B ; A lower triangular assuming unit diagonal
// void dtrmm_lltu(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);
// D <= alpha * A * B ; A upper triangular
// void dtrmm_lunn(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);
// D <= alpha * A * B ; A upper triangular assuming unit diagonal
// void dtrmm_lunu(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);
// D <= alpha * A^T * B ; A upper triangular
// void dtrmm_lutn(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);
// D <= alpha * A^T * B ; A upper triangular assuming unit diagonal
// void dtrmm_lutu(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);
// D <= alpha * B * A ; A lower triangular
void dtrmm_rlnn(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);  // to keep
// D <= alpha * B * A ; A lower triangular assuming unit diagonal
// void dtrmm_rlnu(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);
// D <= alpha * B * A^T ; A lower triangular
// void dtrmm_rltn(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);
// D <= alpha * B * A^T ; A lower triangular assuming unit diagonal
// void dtrmm_rltu(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);
// D <= alpha * B * A ; A upper triangular
// void dtrmm_runn(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);
// D <= alpha * B * A ; A upper triangular assuming unit diagonal
// void dtrmm_runu(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);
// D <= alpha * B * A^T ; A upper triangular
// void dtrmm_rutn(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);
// D <= alpha * B * A^T ; A upper triangular assuming unit diagonal
// void dtrmm_rutu(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);


/*****************************************************
 * triangular matrix solve
 * ***************************************************/
// D <= alpha * A^{-1} * B , with A lower triangular employing explicit inverse of diagonal
// void dtrsm_llnn(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);
// D <= alpha * A^{-1} * B , with A lower triangular assuming unit diagonal
// void dtrsm_llnu(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);
// D <= alpha * A^{-T} * B , with A lower triangular employing explicit inverse of diagonal
// void dtrsm_lltn(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);
// D <= alpha * A^{-T} * B , with A lower triangular assuming unit diagonal
// void dtrsm_lltu(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);
// D <= alpha * A^{-1} * B , with A upper triangular employing explicit inverse of diagonal
// void dtrsm_lunn(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);
// D <= alpha * A^{-1} * B , with A upper triangular assuming unit diagonal
// void dtrsm_lunu(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);
// D <= alpha * A^{-T} * B , with A upper triangular employing explicit inverse of diagonal
// void dtrsm_lutn(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);
// D <= alpha * A^{-T} * B , with A upper triangular assuming unit diagonal
// void dtrsm_lutu(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);
// D <= alpha * B * A^{-1} , with A lower triangular employing explicit inverse of diagonal
// void dtrsm_rlnn(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);  // to keep
// D <= alpha * B * A^{-1} , with A lower triangular assuming unit diagonal
// void dtrsm_rlnu(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);
// D <= alpha * B * A^{-T} , with A lower triangular employing explicit inverse of diagonal
void dtrsm_rltn(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);  // to keep
// D <= alpha * B * A^{-T} , with A lower triangular assuming unit diagonal
// void dtrsm_rltu(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);
// D <= alpha * B * A^{-1} , with A upper triangular employing explicit inverse of diagonal
// void dtrsm_runn(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);
// D <= alpha * B * A^{-1} , with A upper triangular assuming unit diagonal
// void dtrsm_runu(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);
// D <= alpha * B * A^{-T} , with A upper triangular employing explicit inverse of diagonal
// void dtrsm_rutn(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);
// D <= alpha * B * A^{-T} , with A upper triangular assuming unit diagonal
// void dtrsm_rutu(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sD, int di, int dj);


/*****************************************************
 * rank-1 update
 * ***************************************************/
// D <= beta * C + alpha * A * B^T ; C, D lower triangular
void dsyrk_ln(int m, int k, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, double beta, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj);  // to keep
void dsyrk_ln_mn(int m, int n, int k, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, double beta, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj);  // to keep
// D <= beta * C + alpha * A^T * B ; C, D lower triangular
// void dsyrk_lt(int m, int k, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, double beta, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj);
// D <= beta * C + alpha * A * B^T ; C, D upper triangular
// void dsyrk_un(int m, int k, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, double beta, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj);
// D <= beta * C + alpha * A^T * B ; C, D upper triangular
// void dsyrk_ut(int m, int k, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, double beta, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj);

/*****************************************************
 * rank-2 update
 * ***************************************************/
// D <= beta * C + alpha * A * B^T + alpha * B * A^T; C, D lower triangular
// void dsyr2k_ln(int m, int k, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, double beta, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj);
// D <= beta * C + alpha * A^T * B + alpha * B^T * A; C, D lower triangular
// void dsyr2k_lt(int m, int k, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, double beta, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj);
// D <= beta * C + alpha * A * B^T + alpha * B * A^T; C, D upper triangular
// void dsyr2k_un(int m, int k, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, double beta, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj);
// D <= beta * C + alpha * A^T * B + alpha * B^T * A; C, D upper triangular
// void dsyr2k_ut(int m, int k, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, double beta, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj);


#endif
