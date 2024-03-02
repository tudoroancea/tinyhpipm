#ifndef TINYHPIPM_BLAS_MISC_H
#define TINYHPIPM_BLAS_MISC_H
#include "tinyhpipm/blas/struct.h"

// --- insert/extract
//
// sA[ai, aj] <= a
// void dgein1(double a, struct mat* sA, int ai, int aj);
// <= sA[ai, aj]
// double dgeex1(struct mat* sA, int ai, int aj);

// --- set
// A <= alpha
void dgese(int m, int n, double alpha, struct mat* sA, int ai, int aj);  // to keep

// --- copy / scale
// B <= A
void dgecp(int m, int n, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj);  // to keep
// A <= alpha*A
// void dgesc(int m, int n, double alpha, struct mat* sA, int ai, int aj);
// B <= alpha*A
// void dgecpsc(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj);
// B <= A, A lower triangular
void dtrcp_l(int m, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj);  // to keep
// void dtrcpsc_l(int m, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj);
// void dtrsc_l(int m, double alpha, struct mat* sA, int ai, int aj);

// --- sum
// B <= B + alpha*A
void dgead_lib(int m, int n, double alpha, int offsetA, double* A, int sda, int offsetB, double* B, int sdb);
void dgead(int m, int n, double alpha, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj);  // to keep
// y <= y + alpha*x
// void dvecad(int m, double alpha, struct vec* sx, int xi, struct vec* sy, int yi);

// --- traspositions
// B <= A'
void dgetr(int m, int n, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj);  // to keep

// B <= A', A lower triangular
void dtrtr_l(int m, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj);  // to keep

// B <= A', A upper triangular
// void dtrtr_u(int m, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj);

/********************************************************************
 * operations on diagonal
 ********************************************************************/

// diag(A) += alpha
void ddiare(int kmax, double alpha, struct mat* sA, int ai, int aj);  // to keep

// diag(A) <= alpha*x
// void ddiain(int kmax, double alpha, struct vec* sx, int xi, struct mat* sA, int ai, int aj);

// diag(A)[idx] <= alpha*x
// void ddiain_sp(int kmax, double alpha, struct vec* sx, int xi, int* idx, struct mat* sD, int di, int dj);

// x <= diag(A)
void ddiaex_lib(int kmax, double alpha, int offset, double* pD, int sdd, double* x);
void ddiaex(int kmax, double alpha, struct mat* sA, int ai, int aj, struct vec* sx, int xi);  // to keep

// x <= diag(A)[idx]
// void ddiaex_sp(int kmax, double alpha, int* idx, struct mat* sD, int di, int dj, struct vec* sx, int xi);

// diag(A) += alpha*x
// void ddiaad(int kmax, double alpha, struct vec* sx, int xi, struct mat* sA, int ai, int aj);

// diag(A)[idx] += alpha*x
void ddiaex_libsp(int kmax, int* idx, double alpha, double* pD, int sdd, double* x);
void ddiaad_sp(int kmax, double alpha, struct vec* sx, int xi, int* idx, struct mat* sD, int di, int dj);  // to keep

// diag(A)[idx] = y + alpha*x
// void ddiaadin_sp(int kmax, double alpha, struct vec* sx, int xi, struct vec* sy, int yi, int* idx, struct mat* sD, int di, int dj);

/********************************************************************
 * operations on rows
 ********************************************************************/

void drowin(int kmax, double alpha, struct vec* sx, int xi, struct mat* sA, int ai, int aj);  // to keep
void drowex(int kmax, double alpha, struct mat* sA, int ai, int aj, struct vec* sx, int xi);  // to keep
void drowad(int kmax, double alpha, struct vec* sx, int xi, struct mat* sA, int ai, int aj);  // to keep
void drowad_sp(int kmax, double alpha, struct vec* sx, int xi, int* idx, struct mat* sD, int di, int dj);  // to keep
void drowsw(int kmax, struct mat* sA, int ai, int aj, struct mat* sC, int ci, int cj);  // to keep
// void drowpe(int kmax, int* ipiv, struct mat* sA);
// void drowpei(int kmax, int* ipiv, struct mat* sA);

/********************************************************************
 * operations on columns
 ********************************************************************/

void dcolex(int kmax, struct mat* sA, int ai, int aj, struct vec* sx, int xi);  // to keep
void dcolin(int kmax, struct vec* sx, int xi, struct mat* sA, int ai, int aj);  // to keep
void dcolad(int kmax, double alpha, struct vec* sx, int xi, struct mat* sA, int ai, int aj);  // to keep
void dcolsc(int kmax, double alpha, struct mat* sA, int ai, int aj);  // to keep
// void dcolsw(int kmax, struct mat* sA, int ai, int aj, struct mat* sC, int ci, int cj);
// void dcolpe(int kmax, int* ipiv, struct mat* sA);
// void dcolpei(int kmax, int* ipiv, struct mat* sA);

/********************************************************************
 * operations on vectors
 ********************************************************************/

// a <= alpha
void dvecse(int m, double alpha, struct vec* sx, int xi);  // to keep
// sx[xi] <= a
// void dvecin1(double a, struct vec* sx, int xi);
// <= sx[xi]
// double dvecex1(struct vec* sx, int xi);
// y <= x
void dveccp(int m, struct vec* sx, int xi, struct vec* sy, int yi);  // to keep
// x <= alpha*x
void dvecsc(int m, double alpha, struct vec* sx, int xi);  // to keep
// y <= alpha*x
void dveccpsc(int m, double alpha, struct vec* sx, int xi, struct vec* sy, int yi);  // to keep
// z[idx] += alpha * x
void dvecad_sp(int m, double alpha, struct vec* sx, int xi, int* idx, struct vec* sz, int zi);  // to keep
// z[idx] <= alpha * x
// void dvecin_sp(int m, double alpha, struct vec* sx, int xi, int* idx, struct vec* sz, int zi);
// z <= alpha * x[idx]
void dvecex_sp(int m, double alpha, int* idx, struct vec* sx, int xi, struct vec* sz, int zi);  // to keep
// z += alpha * x[idx]
// void dvecexad_sp(int m, double alpha, int* idx, struct vec* sx, int xi, struct vec* sz, int zi);

// void dveccl(int m, struct vec* sxm, int xim, struct vec* sx, int xi, struct vec* sxp, int xip, struct vec* sz, int zi);
// void dveccl_mask(int m, struct vec* sxm, int xim, struct vec* sx, int xi, struct vec* sxp, int xip, struct vec* sz, int zi, struct vec* sm, int mi);

// zero out strvec to strvec with mask
// void dvecze(int m, struct vec* sm, int mi, struct vec* sv, int vi, struct vec* se, int ei);

// compute inf norm of vector
void dvecnrm_inf(int m, struct vec* sx, int xi, double* ptr_norm);  // to keep

// void dvecnrm_2(int m, struct vec* sx, int xi, double* ptr_norm);

// void dvecpe(int kmax, int* ipiv, struct vec* sx, int xi);

// void dvecpei(int kmax, int* ipiv, struct vec* sx, int xi);

#endif  // TINYHPIPM_BLAS_MISC_H
