#ifndef TINYHPIPM_BLAS_BLAS1_H
#define TINYHPIPM_BLAS_BLAS1_H

#include "hpipm/blas/struct.h"

// z = y + alpha*x
// z[zi:zi+n] = alpha*x[xi:xi+n] + y[yi:yi+n]
// NB: Different arguments semantic compare to equivalent standard BLAS routine
void daxpy(int kmax, double alpha, struct vec* sx, int xi, struct vec* sy, int yi, struct vec* sz, int zi);
// z = beta*y + alpha*x
void daxpby(int kmax, double alpha, struct vec* sx, int xi, double beta, struct vec* sy, int yi, struct vec* sz, int zi);
// z = x .* y
void dvecmul(int m, struct vec* sx, int xi, struct vec* sy, int yi, struct vec* sz, int zi);
// z += x .* y
void dvecmulacc(int m, struct vec* sx, int xi, struct vec* sy, int yi, struct vec* sz, int zi);
// z = x .* y, return sum(z) = x^T * y
double dvecmuldot(int m, struct vec* sx, int xi, struct vec* sy, int yi, struct vec* sz, int zi);
// return x^T * y
double ddot(int m, struct vec* sx, int xi, struct vec* sy, int yi);
// construct givens plane rotation
void drotg(double a, double b, double* c, double* s);
// apply plane rotation [a b] [c -s; s; c] to the aj0 and aj1 columns of A at row index ai
void dcolrot(int m, struct mat* sA, int ai, int aj0, int aj1, double c, double s);
// apply plane rotation [c s; -s c] [a; b] to the ai0 and ai1 rows of A at column index aj
void drowrot(int m, struct mat* sA, int ai0, int ai1, int aj, double c, double s);

#endif  // TINYHPIPM_BLAS_BLAS1_H
