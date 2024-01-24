#include <stdio.h>
#include <stdlib.h>

#include "hpipm/blas.h"
#include "hpipm/common.h"
#include "hpipm/dense/d_dense_qcqp.h"
#include "hpipm/dense/d_dense_qcqp_dim.h"
#include "hpipm/dense/d_dense_qcqp_sol.h"


hpipm_size_t d_dense_qcqp_sol_strsize() {
    return sizeof(struct d_dense_qcqp_sol);
}


hpipm_size_t d_dense_qcqp_sol_memsize(struct d_dense_qcqp_dim* dim) {

    int nv = dim->nv;
    int ne = dim->ne;
    int nb = dim->nb;
    int ng = dim->ng;
    int nq = dim->nq;
    int ns = dim->ns;

    hpipm_size_t size = 0;

    size += 4 * sizeof(struct vec);  // v pi lam t

    size += 1 * memsize_vec(nv + 2 * ns);  // ux
    size += 1 * memsize_vec(ne);  // pi
    size += 2 * memsize_vec(2 * nb + 2 * ng + 2 * nq + 2 * ns);  // lam t

    size = (size + 63) / 64 * 64;  // make multiple of typical cache line size
    size += 64;  // align to typical cache line size

    return size;
}


void d_dense_qcqp_sol_create(struct d_dense_qcqp_dim* dim, struct d_dense_qcqp_sol* qp_sol, void* mem) {

    // loop index
    int ii;

    // zero memory (to avoid corrupted memory like e.g. NaN)
    hpipm_size_t memsize = d_dense_qcqp_sol_memsize(dim);
    hpipm_zero_memset(memsize, mem);

    // extract dim
    int nv = dim->nv;
    int ne = dim->ne;
    int nb = dim->nb;
    int ng = dim->ng;
    int nq = dim->nq;
    int ns = dim->ns;


    // vector struct stuff
    struct vec* sv_ptr = (struct vec*) mem;

    qp_sol->v = sv_ptr;
    sv_ptr += 1;
    qp_sol->pi = sv_ptr;
    sv_ptr += 1;
    qp_sol->lam = sv_ptr;
    sv_ptr += 1;
    qp_sol->t = sv_ptr;
    sv_ptr += 1;


    // align to typical cache line size
    hpipm_size_t l_ptr = (hpipm_size_t) sv_ptr;
    l_ptr = (l_ptr + 63) / 64 * 64;


    // double stuff
    char* c_ptr;
    c_ptr = (char*) l_ptr;

    char* tmp_ptr;

    // v
    create_vec(nv + 2 * ns, qp_sol->v, c_ptr);
    c_ptr += qp_sol->v->memsize;
    // pi
    create_vec(ne, qp_sol->pi, c_ptr);
    c_ptr += qp_sol->pi->memsize;
    // lam
    create_vec(2 * nb + 2 * ng + 2 * nq + 2 * ns, qp_sol->lam, c_ptr);
    c_ptr += qp_sol->lam->memsize;
    // t
    create_vec(2 * nb + 2 * ng + 2 * nq + 2 * ns, qp_sol->t, c_ptr);
    c_ptr += qp_sol->t->memsize;


    qp_sol->dim = dim;

    qp_sol->memsize = d_dense_qcqp_sol_memsize(dim);


#if defined(RUNTIME_CHECKS)
    if (c_ptr > ((char*) mem) + qp_sol->memsize) {
        printf("\nCreate_dense_qp_sol: outsize memory bounds!\n\n");
        exit(1);
    }
#endif
}


void d_dense_qcqp_sol_get(char* field, struct d_dense_qcqp_sol* qp_sol, void* value) {
    if (hpipm_strcmp(field, "v")) {
        d_dense_qcqp_sol_get_v(qp_sol, value);
    } else if (hpipm_strcmp(field, "pi")) {
        d_dense_qcqp_sol_get_pi(qp_sol, value);
    } else if (hpipm_strcmp(field, "lam_lb")) {
        d_dense_qcqp_sol_get_lam_lb(qp_sol, value);
    } else if (hpipm_strcmp(field, "lam_ub")) {
        d_dense_qcqp_sol_get_lam_ub(qp_sol, value);
    } else if (hpipm_strcmp(field, "lam_lg")) {
        d_dense_qcqp_sol_get_lam_lg(qp_sol, value);
    } else if (hpipm_strcmp(field, "lam_ug")) {
        d_dense_qcqp_sol_get_lam_ug(qp_sol, value);
    } else if (hpipm_strcmp(field, "lam_uq")) {
        d_dense_qcqp_sol_get_lam_ug(qp_sol, value);
    } else {
        printf("error: d_dense_qcqp_sol_get: wrong field name '%s'. Exiting.\n", field);
        exit(1);
    }
}


void d_dense_qcqp_sol_get_v(struct d_dense_qcqp_sol* qp_sol, double* v) {
    int nv = qp_sol->dim->nv;
    unpack_vec(nv, qp_sol->v, 0, v, 1);
}


void d_dense_qcqp_sol_get_pi(struct d_dense_qcqp_sol* qp_sol, double* pi) {
    int ne = qp_sol->dim->ne;
    unpack_vec(ne, qp_sol->pi, 0, pi, 1);
}


void d_dense_qcqp_sol_get_lam_lb(struct d_dense_qcqp_sol* qp_sol, double* lam_lb) {
    int nb = qp_sol->dim->nb;
    unpack_vec(nb, qp_sol->lam, 0, lam_lb, 1);
}


void d_dense_qcqp_sol_get_lam_ub(struct d_dense_qcqp_sol* qp_sol, double* lam_ub) {
    int nb = qp_sol->dim->nb;
    int ng = qp_sol->dim->ng;
    int nq = qp_sol->dim->nq;
    unpack_vec(nb, qp_sol->lam, nb + ng + nq, lam_ub, 1);
}


void d_dense_qcqp_sol_get_lam_lg(struct d_dense_qcqp_sol* qp_sol, double* lam_lg) {
    int nb = qp_sol->dim->nb;
    int ng = qp_sol->dim->ng;
    unpack_vec(ng, qp_sol->lam, nb, lam_lg, 1);
}


void d_dense_qcqp_sol_get_lam_ug(struct d_dense_qcqp_sol* qp_sol, double* lam_ug) {
    int nb = qp_sol->dim->nb;
    int ng = qp_sol->dim->ng;
    int nq = qp_sol->dim->nq;
    unpack_vec(ng, qp_sol->lam, 2 * nb + ng + nq, lam_ug, 1);
}


void d_dense_qcqp_sol_get_lam_uq(struct d_dense_qcqp_sol* qp_sol, double* lam_uq) {
    int nb = qp_sol->dim->nb;
    int ng = qp_sol->dim->ng;
    int nq = qp_sol->dim->nq;
    unpack_vec(nq, qp_sol->lam, 2 * nb + 2 * ng + nq, lam_uq, 1);
}


void d_dense_qcqp_sol_set_v(double* v, struct d_dense_qcqp_sol* qp_sol) {
    int nv = qp_sol->dim->nv;
    pack_vec(nv, v, 1, qp_sol->v, 0);
}
