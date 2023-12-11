#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <blasfeo/blasfeo_block_size.h>
#include <blasfeo/blasfeo_common.h>
#include <blasfeo/blasfeo_d_aux.h>


/* Explicitly panel-major */

// return the memory size (in bytes) needed for a strmat
size_t blasfeo_pm_memsize_dmat(int ps, int m, int n) {
    int nc = 4;  // hard coded !!! D_PLD
    int al = ps * nc;
    int pm = (m + ps - 1) / ps * ps;
    int cn = (n + nc - 1) / nc * nc;
    int tmp = m < n ? (m + al - 1) / al * al : (n + al - 1) / al * al;  // al(min(m,n)) // XXX max ???
    size_t memsize = (pm * cn + tmp) * sizeof(double);
    return memsize;
}


// create a matrix structure for a matrix of size m*n by using memory passed by a pointer
void blasfeo_pm_create_dmat(int ps, int m, int n, struct blasfeo_pm_dmat* sA, void* memory) {
    sA->mem = memory;
    sA->ps = ps;
    int nc = 4;  // hard coded !!! D_PLD
    int al = ps * nc;
    sA->m = m;
    sA->n = n;
    int pm = (m + ps - 1) / ps * ps;
    int cn = (n + nc - 1) / nc * nc;
    sA->pm = pm;
    sA->cn = cn;
    double* ptr = (double*) memory;
    sA->pA = ptr;
    ptr += pm * cn;
    int tmp = m < n ? (m + al - 1) / al * al : (n + al - 1) / al * al;  // al(min(m,n)) // XXX max ???
    sA->dA = ptr;
    ptr += tmp;
    sA->memsize = (pm * cn + tmp) * sizeof(double);
    sA->use_dA = 0;  // invalidate stored inverse diagonal
    return;
}


// print a matrix structure
#if defined(EXT_DEP)
void blasfeo_pm_print_dmat(int m, int n, struct blasfeo_pm_dmat* sA, int ai, int aj) {
    int ii, jj;
    for (ii = 0; ii < m; ii++) {
        for (jj = 0; jj < n; jj++) {
            printf("%9.5f ", BLASFEO_PM_DMATEL(sA, ai + ii, aj + jj));
        }
        printf("\n");
    }
    printf("\n");
    return;
}
#endif


/* Explicitly column-major */

// return the memory size (in bytes) needed for a strmat
size_t blasfeo_cm_memsize_dmat(int m, int n) {
    int tmp = m < n ? m : n;  // al(min(m,n)) // XXX max ???
    size_t memsize = (m * n + tmp) * sizeof(double);
    return memsize;
}


// create a matrix structure for a matrix of size m*n by using memory passed by a pointer
void blasfeo_cm_create_dmat(int m, int n, struct blasfeo_pm_dmat* sA, void* memory) {
    sA->mem = memory;
    sA->m = m;
    sA->n = n;
    sA->use_dA = 0;  // invalidate stored inverse diagonal
    double* ptr = (double*) memory;
    sA->pA = ptr;
    ptr += m * n;
    int tmp = m < n ? m : n;  // al(min(m,n)) // XXX max ???
    sA->dA = ptr;
    ptr += tmp;
    sA->use_dA = 0;
    sA->memsize = (m * n + tmp) * sizeof(double);
    return;
}
