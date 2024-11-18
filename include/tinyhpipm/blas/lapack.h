#ifndef TINYHPIPM_BLAS_LAPACK_H
#define TINYHPIPM_BLAS_LAPACK_H


#include "tinyhpipm/blas/struct.h"

/*
 * LAPACK
 */

// D <= chol( C ) ; C, D lower triangular
void dpotrf_l(int m, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj);  // to keep
void dpotrf_l_mn(int m, int n, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj);  // to keep
// D <= chol( C ) ; C, D upper triangular
// void dpotrf_u(int m, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj);
// D <= chol( C + A * B' ) ; C, D lower triangular
// D <= chol( C + A * B^T ) ; C, D lower triangular
void dsyrk_dpotrf_ln(int m, int k, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj);  // to keep
void dsyrk_dpotrf_ln_mn(int m, int n, int k, struct mat* sA, int ai, int aj, struct mat* sB, int bi, int bj, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj);  // to keep
// D <= lu( C ) ; no pivoting
// void dgetrf_np(int m, int n, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj);
// D <= lu( C ) ; row pivoting
void dgetrf_rp(int m, int n, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj, int* ipiv);  // to keep
// D <= qr( C )
// int dgeqrf_worksize(int m, int n);  // in bytes
// void dgeqrf(int m, int n, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj, void* work);
// D <= Q factor, where C is the output of the LQ factorization
int dorglq_worksize(int m, int n, int k);  // in bytes // to keep
void dorglq(int m, int n, int k, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj, void* work);  // to keep
// D <= lq( C )
void dgelqf(int m, int n, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj, void* work);  // to keep
int dgelqf_worksize(int m, int n);  // in bytes // to keep
// D <= lq( C ), positive diagonal elements
// void dgelqf_pd(int m, int n, struct mat* sC, int ci, int cj, struct mat* sD, int di, int dj, void* work);
// [L, A] <= lq( [L, A] ), positive diagonal elements, array of matrices, with
// L lower triangular, of size (m)x(m)
// A full, of size (m)x(n1)
// void dgelqf_pd_la(int m, int n1, struct mat* sL, int li, int lj, struct mat* sA, int ai, int aj, void* work);
// [L, L, A] <= lq( [L, L, A] ), positive diagonal elements, array of matrices, with:
// L lower triangular, of size (m)x(m)
// A full, of size (m)x(n1)
// void dgelqf_pd_lla(int m, int n1, struct mat* sL0, int l0i, int l0j, struct mat* sL1, int l1i, int l1j, struct mat* sA, int ai, int aj, void* work);

#endif  // TINYHPIPM_BLAS_LAPACK_H
