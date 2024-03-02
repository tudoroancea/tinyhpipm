#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "tinyhpipm/blas/blas1.h"
#include "tinyhpipm/blas/kernel.h"
#include "tinyhpipm/blas/struct.h"


void daxpy(int m, double alpha, struct vec* sx, int xi, struct vec* sy, int yi, struct vec* sz, int zi) {

    if (m <= 0) { return; }


    double* x = sx->pa + xi;
    double* y = sy->pa + yi;
    double* z = sz->pa + zi;

    int ii;

    ii = 0;
    for (; ii < m - 3; ii += 4) {
        z[ii + 0] = y[ii + 0] + alpha * x[ii + 0];
        z[ii + 1] = y[ii + 1] + alpha * x[ii + 1];
        z[ii + 2] = y[ii + 2] + alpha * x[ii + 2];
        z[ii + 3] = y[ii + 3] + alpha * x[ii + 3];
    }
    for (; ii < m; ii++) {
        z[ii + 0] = y[ii + 0] + alpha * x[ii + 0];
    }
}


void daxpby(int m, double alpha, struct vec* sx, int xi, double beta, struct vec* sy, int yi, struct vec* sz, int zi) {

    if (m <= 0) { return; }


    double* x = sx->pa + xi;
    double* y = sy->pa + yi;
    double* z = sz->pa + zi;

    int ii;

    ii = 0;
    for (; ii < m - 3; ii += 4) {
        z[ii + 0] = beta * y[ii + 0] + alpha * x[ii + 0];
        z[ii + 1] = beta * y[ii + 1] + alpha * x[ii + 1];
        z[ii + 2] = beta * y[ii + 2] + alpha * x[ii + 2];
        z[ii + 3] = beta * y[ii + 3] + alpha * x[ii + 3];
    }
    for (; ii < m; ii++) {
        z[ii + 0] = beta * y[ii + 0] + alpha * x[ii + 0];
    }
}


// multiply two vectors
void dvecmul(int m, struct vec* sx, int xi, struct vec* sy, int yi, struct vec* sz, int zi) {

    if (m <= 0) { return; }


    double* x = sx->pa + xi;
    double* y = sy->pa + yi;
    double* z = sz->pa + zi;
    int ii;

    ii = 0;

    for (; ii < m; ii++) {
        z[ii + 0] = x[ii + 0] * y[ii + 0];
    }
}


// multiply two vectors and add result to another vector
void dvecmulacc(int m, struct vec* sx, int xi, struct vec* sy, int yi, struct vec* sz, int zi) {

    if (m <= 0) { return; }


    double* x = sx->pa + xi;
    double* y = sy->pa + yi;
    double* z = sz->pa + zi;
    int ii;

    ii = 0;
    for (; ii < m; ii++) {
        z[ii + 0] += x[ii + 0] * y[ii + 0];
    }
}


// multiply two vectors and compute ddot product
double dvecmuldot(int m, struct vec* sx, int xi, struct vec* sy, int yi, struct vec* sz, int zi) {
    if (m <= 0) { return 0.0; }

    double* x = sx->pa + xi;
    double* y = sy->pa + yi;
    double* z = sz->pa + zi;
    int ii;
    double ddot = 0.0;
    ii = 0;
    for (; ii < m; ii++) {
        z[ii + 0] = x[ii + 0] * y[ii + 0];
        ddot += z[ii + 0];
    }
    return ddot;
}


// compute ddot product of two vectors
double ddot(int m, struct vec* sx, int xi, struct vec* sy, int yi) {

    if (m <= 0) { return 0.0; }

    double* x = sx->pa + xi;
    double* y = sy->pa + yi;
    int ii;
    double ddot = 0.0;

    ii = 0;
    for (; ii < m - 3; ii += 4) {
        ddot += x[ii + 0] * y[ii + 0];
        ddot += x[ii + 1] * y[ii + 1];
        ddot += x[ii + 2] * y[ii + 2];
        ddot += x[ii + 3] * y[ii + 3];
    }
    for (; ii < m; ii++) {
        ddot += x[ii + 0] * y[ii + 0];
    }
    return ddot;
}


void drotg(double a, double b, double* c, double* s) {
    double aa = fabs(a);
    double bb = fabs(b);
    double roe = aa >= bb ? a : b;
    double scale = aa + bb;
    double r;
    if (scale == 0) {
        *c = 1.0;
        *s = 0.0;
    } else {
        aa = a / scale;
        bb = b / scale;
        r = scale * sqrt(aa * aa + bb * bb);
        r = r * (roe >= 0 ? 1 : -1);
        *c = a / r;
        *s = b / r;
    }
}


void dcolrot(int m, struct mat* sA, int ai, int aj0, int aj1, double c, double s) {
    const int ps = 4;
    int sda = sA->cn;
    double* px = sA->pA + ai / ps * ps * sda + ai % ps + aj0 * ps;
    double* py = sA->pA + ai / ps * ps * sda + ai % ps + aj1 * ps;
    int mna = (ps - ai % ps) % ps;
    int ii;
    double d_tmp;
    ii = 0;
    if (mna > 0) {
        for (; ii < mna; ii++) {
            d_tmp = c * px[0] + s * py[0];
            py[0] = c * py[0] - s * px[0];
            px[0] = d_tmp;
            px++;
            py++;
        }
        px += ps * (sda - 1);
        py += ps * (sda - 1);
    }
    for (; ii < m - 3; ii += 4) {
        //
        d_tmp = c * px[0] + s * py[0];
        py[0] = c * py[0] - s * px[0];
        px[0] = d_tmp;
        //
        d_tmp = c * px[1] + s * py[1];
        py[1] = c * py[1] - s * px[1];
        px[1] = d_tmp;
        //
        d_tmp = c * px[2] + s * py[2];
        py[2] = c * py[2] - s * px[2];
        px[2] = d_tmp;
        //
        d_tmp = c * px[3] + s * py[3];
        py[3] = c * py[3] - s * px[3];
        px[3] = d_tmp;
        //
        px += ps * sda;
        py += ps * sda;
    }
    for (; ii < m; ii++) {
        //
        d_tmp = c * px[0] + s * py[0];
        py[0] = c * py[0] - s * px[0];
        px[0] = d_tmp;
        //
        px++;
        py++;
    }
}


void drowrot(int m, struct mat* sA, int ai0, int ai1, int aj, double c, double s) {
    const int ps = 4;
    int sda = sA->cn;
    double* px = sA->pA + ai0 / ps * ps * sda + ai0 % ps + aj * ps;
    double* py = sA->pA + ai1 / ps * ps * sda + ai1 % ps + aj * ps;
    int ii;
    double d_tmp;
    ii = 0;
    for (; ii < m - 3; ii += 4) {
        //
        d_tmp = c * px[0 * ps] + s * py[0 * ps];
        py[0 * ps] = c * py[0 * ps] - s * px[0 * ps];
        px[0 * ps] = d_tmp;
        //
        d_tmp = c * px[1 * ps] + s * py[1 * ps];
        py[1 * ps] = c * py[1 * ps] - s * px[1 * ps];
        px[1 * ps] = d_tmp;
        //
        d_tmp = c * px[2 * ps] + s * py[2 * ps];
        py[2 * ps] = c * py[2 * ps] - s * px[2 * ps];
        px[2 * ps] = d_tmp;
        //
        d_tmp = c * px[3 * ps] + s * py[3 * ps];
        py[3 * ps] = c * py[3 * ps] - s * px[3 * ps];
        px[3 * ps] = d_tmp;
        //
        px += 4 * ps;
        py += 4 * ps;
    }
    for (; ii < m; ii++) {
        //
        d_tmp = c * px[0 * ps] + s * py[0 * ps];
        py[0 * ps] = c * py[0 * ps] - s * px[0 * ps];
        px[0 * ps] = d_tmp;
        //
        px += 1 * ps;
        py += 1 * ps;
    }
}
