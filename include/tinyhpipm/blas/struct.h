#ifndef TINYHPIPM_BLAS_STRUCT_H
#define TINYHPIPM_BLAS_STRUCT_H

#ifdef __cplusplus
extern "C" {
#endif

#include "tinyhpipm/common.h"

#define D_PS 4  // panel size
#define D_PLD 4  // GCD of panel length
#define CACHE_LINE_SIZE 64  // bytes, holds 8 double precision numbers (on ARM Cortex A76 we have 64B, on Apple M2 we have 128B)
#define K_MAX_STACK 300  // Maximum inner product length K for buffer allocation on stack (decrease this value if stack size is exceeded)


/************************************************
 * matrix and vector struct
 ************************************************/

/*
 * matrix struct in Panel-Major format. The data is stored by sets of rows (panels)
 * and each panel is stored in column-major format. This allows the elements to be stored in the same
 * order they are accessed by gemm kerels.
 *
 * For alignment and memory contiguity reasons, we actually store more rows and columns
 * than needed to have a multiple of D_PS rows and D_PLD columns.
 */
struct mat {
    double* mem;  // pointer to passed chunk of memory
    double* pA;  // pointer to a pm*pn array of doubles, the first is aligned to cache line size
    double* dA;  // pointer to a min(m,n) (or max???) array of doubles
    int m;  // rows
    int n;  // cols
    int pm;  // packed number or rows, multiple of D_PS
    int cn;  // packed number or cols, multiplie of D_PLD
    int use_dA;  // flag to tell if dA can be used
    int memsize;  // size of needed memory
};
struct vec {
    double* mem;  // pointer to passed chunk of memory
    double* pa;  // pointer to a pm array of doubles, the first is aligned to cache line size
    int m;  // size
    int pm;  // packed size
    int memsize;  // size of needed memory
};

/************************************************
 * matrix and vector indexing
 ************************************************/

#define MATEL(sA, ai, aj) ((sA)->pA[((ai) - ((ai) & (D_PS - 1))) * (sA)->cn + (aj) * D_PS + ((ai) & (D_PS - 1))])
#define VECEL(sa, ai) ((sa)->pa[ai])


/************************************************
 * matrix and vector creation and destruction
 ************************************************/

// computes the memory needed for the matrix content, with the necessary padding,
// but with all the other variables in the struct (pointers, ints for dimensions, etc.),
// i.e. only what has the be allocated dynamically
hpipm_size_t memsize_mat(int m, int n);
// computes the memory needed for the vector content, with the necessary padding,
// but with all the other variables in the struct (pointers, ints for dimensions, etc.),
// i.e. only what has the be allocated dynamically
hpipm_size_t memsize_vec(int m);

// takes in the
void create_mat(int m, int n, struct mat* sA, void* memory);
// create a vec struct of size m by using memory passed by a pointer (pointer is not updated)
void create_vec(int m, struct vec* sA, void* memory);

/*
// combines (in a single call) the create_mat() and allocate_mat() functions (pointer is updated)
void allocate_mat(int m, int n, struct mat* sA);
// combines (in a single call) the create_vec() and allocate_vec() functions (pointer is updated)
void allocate_vec(int m, struct vec* sa);
*/

// free memory of a mat structure
void free_mat(struct mat* sA);
// free memory of a vec structure
void free_vec(struct vec* sa);


/************************************************
 * conversion between mat/vec and double*
 ************************************************/

// pack the column-major matrix A (with leading dimension lda) into the matrix struct B (at row and col offsets bi and bj)
void pack_mat(int m, int n, double* A, int lda, struct mat* sB, int bi, int bj);
// pack the lower-triangular column-major matrix A (with leading dimension lda) into the matrix struct B (at row and col offsets bi and bj)
void pack_l_mat(int m, int n, double* A, int lda, struct mat* sB, int bi, int bj);
// pack the upper-triangular column-major matrix A (with leading dimension lda) into the matrix struct B (at row and col offsets bi and bj)
void pack_u_mat(int m, int n, double* A, int lda, struct mat* sB, int bi, int bj);
// transpose and pack the column-major matrix A (with leading dimension lda) into the matrix struct B (at row and col offsets bi and bj)
void pack_tran_mat(int m, int n, double* A, int lda, struct mat* sB, int bi, int bj);
// pack the vector x (using increment incx) into the vector structure y (at offset yi)
void pack_vec(int m, double* x, int incx, struct vec* sy, int yi);
// unpack the matrix structure A (at row and col offsets ai and aj) into the column-major matrix B (with leading dimension ldb)
void unpack_mat(int m, int n, struct mat* sA, int ai, int aj, double* B, int ldb);
// transpose and unpack the matrix structure A (at row and col offsets ai and aj) into the column-major matrix B (with leading dimension ldb)
void unpack_tran_mat(int m, int n, struct mat* sA, int ai, int aj, double* B, int ldb);
// unpack the vector structure x (at offset xi) into the vector y (using increment incy)
void unpack_vec(int m, struct vec* sx, int xi, double* y, int incy);

#ifdef __cplusplus
}
#endif

#endif  // TINYHPIPM_BLAS_STRUCT_H
