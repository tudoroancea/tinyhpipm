#ifndef TINYHPIPM_BLAS_MISC_H
#define TINYHPIPM_BLAS_MISC_H
#include "tinyhpipm/blas/struct.h"

/***************************************************************************************
 *  get and set
 ***************************************************************************************/
// sx[xi] <= a
// void dvecin1(double a, struct vec* sx, int xi);
// <= sx[xi]
// double dvecex1(struct vec* sx, int xi);
// z[idx] <= alpha * x
// void dvecin_sp(int m, double alpha, struct vec* sx, int xi, int* idx, struct vec* sz, int zi);
/*
 *@brief extract part of a vector x into a vector y based on a sparsity pattern and scale it
 *       y[idx] = alpha * x[idx]
 *
 * @param[in] m number of elements in the vector
 * @param[in] alpha scaling factor
 * @param[in] idx indices of the elements to extract, all relative to the starting index of the vector
 * @param[in] sx vector struct
 * @param[in] xi starting index of the vector
 * @param[out] sz vector struct
 * @param[in] zi starting index of the vector
 */
void dvecex_sp(int m, double alpha, int* idx, struct vec* sx, int xi, struct vec* sz, int zi);
// sA[ai, aj] <= a
// void dgein1(double a, struct mat* sA, int ai, int aj);
// <= sA[ai, aj]
// double dgeex1(struct mat* sA, int ai, int aj);
/*
 * @brief insert a (scaled) vector into a row (i.e. copies it into the row, no row is added)
 *        A[ai, aj:aj+kmax] = alpha * x[xi:xi+kmax]
 *
 * @param[in] kmax number of elements in the vector
 * @param[in] alpha scaling factor
 * @param[in] sx vector struct
 * @param[in] xi starting index of the vector
 * @param[in,out] sA matrix struct
 * @param[in] ai row index
 * @param[in] aj starting column index
 */
void drowin(int kmax, double alpha, struct vec* sx, int xi, struct mat* sA, int ai, int aj);
/*
 * @brief extract a (scaled) vector from a row
 *        x[xi:xi+kmax] = alpha * A[ai, aj:aj+kmax]
 *
 * @param[in] kmax number of elements in the vector
 * @param[in] alpha scaling factor
 * @param[in,out] sA matrix struct
 * @param[in] ai row index
 * @param[in] aj starting column index
 * @param[in] sx vector struct
 * @param[in] xi starting index of the vector
 */
void drowex(int kmax, double alpha, struct mat* sA, int ai, int aj, struct vec* sx, int xi);
void dcolex(int kmax, struct mat* sA, int ai, int aj, struct vec* sx, int xi);
void dcolin(int kmax, struct vec* sx, int xi, struct mat* sA, int ai, int aj);
// diag(A) <= alpha*x
// void ddiain(int kmax, double alpha, struct vec* sx, int xi, struct mat* sA, int ai, int aj);
// diag(A)[idx] <= alpha*x
// void ddiain_sp(int kmax, double alpha, struct vec* sx, int xi, int* idx, struct mat* sD, int di, int dj);
// x <= alpha * diag(A)
void ddiaex_lib(int kmax, double alpha, int offset, double* pD, int sdd, double* x);
void ddiaex(int kmax, double alpha, struct mat* sA, int ai, int aj, struct vec* sx, int xi);
// x <= diag(A)[idx]
// void ddiaex_sp(int kmax, double alpha, int* idx, struct mat* sD, int di, int dj, struct vec* sx, int xi);

/*
 * @brief sets all elements of a vector x to a constant alpha
 *        x[xi:xi+m] <= alpha
 *
 * @param[in] m number of elements
 * @param[in] alpha constant
 * @param[in,out] sx vector struct
 * @param[in] xi starting index
 */
void dvecse(int m, double alpha, struct vec* sx, int xi);
// zero out strvec to strvec with mask
// void dvecze(int m, struct vec* sm, int mi, struct vec* sv, int vi, struct vec* se, int ei);
/*
 * @brief sets all elements of a matrix A to a constant alpha
 *        A <= alpha
 *
 * @param[in] m number of rows
 * @param[in] n number of columns
 * @param[in] alpha constant
 * @param[in,out] sA matrix struct
 * @param[in] ai starting row index
 * @param[in] aj starting column index
 */
void dgese(int m, int n, double alpha, struct mat* sA, int ai, int aj);

/***************************************************************************************
 *  copy
 ***************************************************************************************/
/*
 * @brief Copy part of a vector into another vector.
 *        y[yi:yi+m] = x[xi:xi+m]
 *
 * @param[in] m number of elements to copy
 * @param[in] sx source vector struct
 * @param[in] xi starting index of source vector
 * @param[out] sy destination vector struct
 * @param[in] yi starting index of destination vector
 */
void dveccp(int m, struct vec* sx, int xi, struct vec* sy, int yi);
/*
 * @brief copy part of a general matrix A into another matrix B
 *        B[bi:bi+m, bj:bj+n] = A[ai:ai+m, aj:aj+n]
 *
 * @param[in] m number of rows to copy
 * @param[in] n number of columns to copy
 * @param[in,out] sA source matrix struct
 * @param[in] ai starting row index
 * @param[in] aj starting column index
 * @param[in,out] sB destination matrix struct
 * @param[in] bi starting row index
 * @param[in] bj starting column index
 */
void dgecp(int m, int n, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj);
/*
 * @brief copy part of the lower triangular part of a matrix A into another matrix B
 *        B[bi:bi+m, bj:bj+n] = A[ai:ai+m, aj:aj+n]
 *
 * @param[in] m number of rows to copy
 * @param[in] n number of columns to copy
 * @param[in,out] sA source matrix struct
 * @param[in] ai starting row index
 * @param[in] aj starting column index
 * @param[in,out] sB destination matrix struct
 * @param[in] bi starting row index
 * @param[in] bj starting column index
 */
void dtrcp_l(int m, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj);

/***************************************************************************************
 *  transpositions
 ***************************************************************************************/
// B <= A'
void dgetr(int m, int n, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj);
// B <= A', A lower triangular
void dtrtr_l(int m, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj);
// B <= A', A upper triangular
// void dtrtr_u(int m, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj);

/********************************************************************
 * Swap
 ********************************************************************/
/*
 * @brief swap two rows of two matrix structs
 *        B[bi, bj:bj+kmax] = A[ai, aj:aj+kmax]
 *        A[ai, aj:aj+kmax] = B[bi, bj:bj+kmax]
 *
 * @param[in] kmax number of elements
 * @param[in,out] sA matrix struct
 * @param[in] ai row index
 * @param[in] aj starting column index
 * @param[in,out] sB matrix struct
 * @param[in] bi row index
 * @param[in] bj starting column index
 */
void drowsw(int kmax, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj);
// void dcolsw(int kmax, struct mat* sA, int ai, int aj, struct mat* sC, int ci, int cj);


/***************************************************************************************
 * extended blas level 1 routines
 ***************************************************************************************/
// B <= B + alpha*A
void dgead(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj);
// y <= y + alpha*x
// void dvecad(int m, double alpha, struct vec* sx, int xi, struct vec* sy, int yi);
// diag(A) += alpha*x
// void ddiaad(int kmax, double alpha, struct vec* sx, int xi, struct mat* sA, int ai, int aj);
// diag(A)[idx] += alpha*x
void ddiaad_sp(int kmax, double alpha, struct vec* sx, int xi, const int* idx, struct mat* sD, int di, int dj);
// diag(A)[idx] = y + alpha*x
// void ddiaadin_sp(int kmax, double alpha, struct vec* sx, int xi, struct vec* sy, int yi, int* idx, struct mat* sD, int di, int dj);
// diag(A) += alpha
void ddiare(int kmax, double alpha, struct mat* sA, int ai, int aj);
/*
 * @brief add a (scaled) vector to a row
 *        A[ai, aj:aj+kmax] += alpha * x[xi:xi+kmax]
 *
 * @param[in] kmax number of elements in the vector
 * @param[in] alpha scaling factor
 * @param[in] sx vector struct
 * @param[in] xi starting index of the vector
 * @param[in,out] sA matrix struct
 * @param[in] ai row index
 * @param[in] aj starting column index
 */
void drowad(int kmax, double alpha, struct vec* sx, int xi, struct mat* sA, int ai, int aj);
/*
 * @brief add a (scaled) vector to a row, sparse formulation
 *        A[ai, aj:aj+kmax] += alpha * x[xi:xi+kmax]
 *
 * @param[in] kmax number of elements in the vector
 * @param[in] alpha scaling factor
 * @param[in] sx vector struct
 * @param[in] xi starting index of the vector
 * @param[in] idx indices of the elements to add
 * @param[in,out] sA diagonal struct
 * @param[in] ai starting index of the diagonal
 * @param[in] aj starting index of the diagonal
 */
void drowad_sp(int kmax, double alpha, struct vec* sx, int xi, int* idx, struct mat* sA, int ai, int aj);
void dcolad(int kmax, double alpha, struct vec* sx, int xi, struct mat* sA, int ai, int aj);

// z[idx] += alpha * x
void dvecad_sp(int m, double alpha, struct vec* sx, int xi, int* idx, struct vec* sz, int zi);
// z += alpha * x[idx]
// void dvecexad_sp(int m, double alpha, int* idx, struct vec* sx, int xi, struct vec* sz, int zi);

// x <= alpha*x
void dvecsc(int m, double alpha, struct vec* sx, int xi);
// y <= alpha*x
void dveccpsc(int m, double alpha, struct vec* sx, int xi, struct vec* sy, int yi);
// A <= alpha*A
// void dgesc(int m, int n, double alpha, struct mat* sA, int ai, int aj);
void dcolsc(int kmax, double alpha, struct mat* sA, int ai, int aj);
// B <= alpha*A
// void dgecpsc(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj);
// A <= alpha*A, (A lower triangular)
// void dtrsc_l(int m, double alpha, struct mat* sA, int ai, int aj);
// B <= alpha*A, (A lower triangular)
// void dtrcpsc_l(int m, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj);


/***************************************************************************************
 *  clipping
 ***************************************************************************************/
// void dveccl(int m, struct vec* sxm, int xim, struct vec* sx, int xi, struct vec* sxp, int xip, struct vec* sz, int zi);
// void dveccl_mask(int m, struct vec* sxm, int xim, struct vec* sx, int xi, struct vec* sxp, int xip, struct vec* sz, int zi, struct vec* sm, int mi);

/***************************************************************************************
 *  norms
 ***************************************************************************************/
/*
 * @brief compute inf norm of vector
 *        norm = max_{i=0,...,m-1} |x[xi+i]|
 *
 * @param[in] m number of elements
 * @param[in,out] sx vector struct
 * @param[in] xi starting index
 * @param[out] ptr_norm pointer to the computed norm
 */
void dvecnrm_inf(int m, struct vec* sx, int xi, double* ptr_norm);

// void dvecnrm_2(int m, struct vec* sx, int xi, double* ptr_norm);

/***************************************************************************************
 * permutations
 ***************************************************************************************/
// void dvecpe(int kmax, int* ipiv, struct vec* sx, int xi);
// void dvecpei(int kmax, int* ipiv, struct vec* sx, int xi);
// void drowpe(int kmax, int* ipiv, struct mat* sA);
// void drowpei(int kmax, int* ipiv, struct mat* sA);
// void dcolpe(int kmax, int* ipiv, struct mat* sA);
// void dcolpei(int kmax, int* ipiv, struct mat* sA);

#endif  // TINYHPIPM_BLAS_MISC_H
