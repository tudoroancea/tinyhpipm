#include <stdio.h>
#include <stdlib.h>

#include "tinyhpipm/blas.h"
#include "tinyhpipm/common.h"
#include "tinyhpipm/dense/d_dense_qcqp.h"
#include "tinyhpipm/dense/d_dense_qcqp_dim.h"
#include "tinyhpipm/dense/d_dense_qcqp_res.h"
#include "tinyhpipm/dense/d_dense_qcqp_sol.h"


hpipm_size_t d_dense_qcqp_res_memsize(struct d_dense_qcqp_dim* dim) {

    // loop index
    int ii;

    // extract ocp qp size
    int nv = dim->nv;
    int ne = dim->ne;
    int nb = dim->nb;
    int ng = dim->ng;
    int nq = dim->nq;
    int ns = dim->ns;

    hpipm_size_t size = 0;

    size += 4 * sizeof(struct vec);  // res_g res_b res_d res_m

    size += 1 * memsize_vec(nv + 2 * ns);  // res_g
    size += 1 * memsize_vec(ne);  // res_b
    size += 2 * memsize_vec(2 * nb + 2 * ng + 2 * nq + 2 * ns);  // res_d res_m

    size = (size + 63) / 64 * 64;  // make multiple of typical cache line size
    size += 1 * 64;  // align once to typical cache line size

    return size;
}


void d_dense_qcqp_res_create(struct d_dense_qcqp_dim* dim, struct d_dense_qcqp_res* res, void* mem) {

    // loop index
    int ii;

    // zero memory (to avoid corrupted memory like e.g. NaN)
    hpipm_size_t memsize = d_dense_qcqp_res_memsize(dim);
    hpipm_zero_memset(memsize, mem);

    // extract ocp qp size
    int nv = dim->nv;
    int ne = dim->ne;
    int nb = dim->nb;
    int ng = dim->ng;
    int nq = dim->nq;
    int ns = dim->ns;


    // vector struct
    struct vec* sv_ptr = (struct vec*) mem;

    res->res_g = sv_ptr;
    sv_ptr += 1;
    res->res_b = sv_ptr;
    sv_ptr += 1;
    res->res_d = sv_ptr;
    sv_ptr += 1;
    res->res_m = sv_ptr;
    sv_ptr += 1;


    // align to typical cache line size
    hpipm_size_t s_ptr = (hpipm_size_t) sv_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;


    // void stuf
    char* c_ptr = (char*) s_ptr;

    create_vec(nv + 2 * ns, res->res_g, c_ptr);
    c_ptr += (res->res_g)->memsize;

    create_vec(ne, res->res_b, c_ptr);
    c_ptr += (res->res_b)->memsize;

    create_vec(2 * nb + 2 * ng + 2 * nq + 2 * ns, res->res_d, c_ptr);
    c_ptr += (res->res_d)->memsize;

    create_vec(2 * nb + 2 * ng + 2 * nq + 2 * ns, res->res_m, c_ptr);
    c_ptr += (res->res_m)->memsize;

    res->dim = dim;

    res->memsize = d_dense_qcqp_res_memsize(dim);


#if defined(RUNTIME_CHECKS)
    if (c_ptr > ((char*) mem) + res->memsize) {
        printf("\ncreate_dense_qcpp_res: outsize memory bounds!\n\n");
        exit(1);
    }
#endif
}


hpipm_size_t d_dense_qcqp_res_ws_memsize(struct d_dense_qcqp_dim* dim) {

    // loop index
    int ii;

    // extract ocp qp size
    int nv = dim->nv;
    int ne = dim->ne;
    int nb = dim->nb;
    int ng = dim->ng;
    int nq = dim->nq;
    int ns = dim->ns;

    hpipm_size_t size = 0;

    size += 7 * sizeof(struct vec);  // 2*tmp_nv 2*tmp_nbgq tmp_ns q_fun q_adj

    size += 3 * memsize_vec(nv);  // 2*tmp_nv q_adj
    size += 2 * memsize_vec(nb + ng + nq);  // 2*tmp_nbgq
    size += 1 * memsize_vec(ns);  // tmp_ns
    size += 1 * memsize_vec(nq);  // q_fun

    size = (size + 63) / 64 * 64;  // make multiple of typical cache line size
    size += 1 * 64;  // align once to typical cache line size

    return size;
}


void d_dense_qcqp_res_ws_create(struct d_dense_qcqp_dim* dim, struct d_dense_qcqp_res_ws* ws, void* mem) {

    // loop index
    int ii;

    // zero memory (to avoid corrupted memory like e.g. NaN)
    hpipm_size_t memsize = d_dense_qcqp_res_ws_memsize(dim);
    hpipm_zero_memset(memsize, mem);

    // extract ocp qp size
    int nv = dim->nv;
    int ne = dim->ne;
    int nb = dim->nb;
    int ng = dim->ng;
    int nq = dim->nq;
    int ns = dim->ns;


    // vector struct
    struct vec* sv_ptr = (struct vec*) mem;

    ws->tmp_nv = sv_ptr;
    sv_ptr += 2;
    ws->tmp_nbgq = sv_ptr;
    sv_ptr += 2;
    ws->tmp_ns = sv_ptr;
    sv_ptr += 1;
    ws->q_fun = sv_ptr;
    sv_ptr += 1;
    ws->q_adj = sv_ptr;
    sv_ptr += 1;


    // align to typical cache line size
    hpipm_size_t s_ptr = (hpipm_size_t) sv_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;


    // void stuf
    char* c_ptr = (char*) s_ptr;


    create_vec(nv, ws->tmp_nv + 0, c_ptr);
    c_ptr += (ws->tmp_nv + 0)->memsize;

    create_vec(nv, ws->tmp_nv + 1, c_ptr);
    c_ptr += (ws->tmp_nv + 1)->memsize;

    create_vec(nb + ng + nq, ws->tmp_nbgq + 0, c_ptr);
    c_ptr += (ws->tmp_nbgq + 0)->memsize;

    create_vec(nb + ng + nq, ws->tmp_nbgq + 1, c_ptr);
    c_ptr += (ws->tmp_nbgq + 1)->memsize;

    create_vec(ns, ws->tmp_ns + 0, c_ptr);
    c_ptr += (ws->tmp_ns + 0)->memsize;

    create_vec(nq, ws->q_fun, c_ptr);
    c_ptr += (ws->q_fun)->memsize;

    create_vec(nv, ws->q_adj, c_ptr);
    c_ptr += (ws->q_adj)->memsize;

    ws->use_q_fun = 0;
    ws->use_q_adj = 0;

    ws->memsize = d_dense_qcqp_res_ws_memsize(dim);


#if defined(RUNTIME_CHECKS)
    if (c_ptr > ((char*) mem) + ws->memsize) {
        printf("\ncreate_dense_qp_res_workspace: outsize memory bounds!\n\n");
        exit(1);
    }
#endif
}


void d_dense_qcqp_res_compute(struct d_dense_qcqp* qp, struct d_dense_qcqp_sol* qp_sol, struct d_dense_qcqp_res* res, struct d_dense_qcqp_res_ws* ws) {

    int ii, idx;

    int nv = qp->dim->nv;
    int ne = qp->dim->ne;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;
    int nq = qp->dim->nq;
    int ns = qp->dim->ns;

    int nvt = nv + 2 * ns;
    int net = ne;
    int nct = 2 * nb + 2 * ng + 2 * nq + 2 * ns;

    double nct_inv = 1.0 / nct;

    struct mat* Hg = qp->Hv;
    struct mat* A = qp->A;
    struct mat* Ct = qp->Ct;
    struct mat* Hq = qp->Hq;
    struct vec* gz = qp->gz;
    struct vec* b = qp->b;
    struct vec* d = qp->d;
    //	struct vec *gq = qp->gq;
    struct vec* m = qp->m;
    int* idxb = qp->idxb;
    struct vec* Z = qp->Z;
    int* idxs_rev = qp->idxs_rev;

    struct vec* v = qp_sol->v;
    struct vec* pi = qp_sol->pi;
    struct vec* lam = qp_sol->lam;
    struct vec* t = qp_sol->t;

    struct vec* res_g = res->res_g;
    struct vec* res_b = res->res_b;
    struct vec* res_d = res->res_d;
    struct vec* res_m = res->res_m;

    struct vec* tmp_nv = ws->tmp_nv;
    struct vec* tmp_nbgq = ws->tmp_nbgq;
    struct vec* tmp_ns = ws->tmp_ns;

    double mu, tmp;
    res->obj = 0.0;

    // res g
    //	dsymv_l(nv, 1.0, Hg, 0, 0, v, 0, 1.0, gz, 0, res_g, 0);
    dsymv_l(nv, 1.0, Hg, 0, 0, v, 0, 2.0, gz, 0, res_g, 0);
    res->obj += 0.5 * ddot(nv, res_g, 0, v, 0);
    daxpy(nv, -1.0, gz, 0, res_g, 0, res_g, 0);

    if (nb + ng + nq > 0) {
        daxpy(nb + ng + nq, -1.0, lam, 0, lam, nb + ng + nq, tmp_nbgq + 0, 0);
        //		daxpy(nb+ng,  1.0, d, 0, t, 0, res_d, 0);
        //		daxpy(nb+ng,  1.0, d, nb+ng, t, nb+ng, res_d, nb+ng);
        daxpy(2 * nb + 2 * ng + 2 * nq, 1.0, d, 0, t, 0, res_d, 0);
        // box
        if (nb > 0) {
            dvecad_sp(nb, 1.0, tmp_nbgq + 0, 0, idxb, res_g, 0);
            dvecex_sp(nb, 1.0, idxb, v, 0, tmp_nbgq + 1, 0);
        }
        // general
        if (ng > 0) {
            dgemv_nt(nv, ng, 1.0, 1.0, Ct, 0, 0, tmp_nbgq + 0, nb, v, 0, 1.0, 0.0, res_g, 0, tmp_nbgq + 1, nb, res_g, 0, tmp_nbgq + 1, nb);
        }
        // quadratic
        if (nq > 0) {
            //			daxpy(nq,  1.0, d, 2*nb+2*ng+2*ns, t, 2*nb+2*ng+2*ns, res_d, 2*nb+2*ng+2*ns);
            if (ws->use_q_fun & ws->use_q_adj) {
                dveccp(nq, ws->q_fun, 0, tmp_nbgq + 1, nb + ng);
                daxpy(nv, 1.0, ws->q_adj, 0, res_g, 0, res_g, 0);
            } else {
                for (ii = 0; ii < nq; ii++) {
                    dsymv_l(nv, 1.0, Hq + ii, 0, 0, v, 0, 0.0, tmp_nv + 0, 0, tmp_nv + 0, 0);
                    tmp = VECEL(tmp_nbgq + 0, nb + ng + ii);
                    daxpy(nv, tmp, tmp_nv + 0, 0, res_g, 0, res_g, 0);
                    dcolex(nv, Ct, 0, ng + ii, tmp_nv + 1, 0);
                    daxpy(nv, tmp, tmp_nv + 1, 0, res_g, 0, res_g, 0);
                    daxpy(nv, 0.5, tmp_nv + 0, 0, tmp_nv + 1, 0, tmp_nv + 0, 0);
                    tmp = ddot(nv, tmp_nv + 0, 0, v, 0);
                    VECEL(tmp_nbgq + 1, nb + ng + ii) = tmp;
                }
            }
        }
        daxpy(nb + ng + nq, -1.0, tmp_nbgq + 1, 0, res_d, 0, res_d, 0);
        daxpy(nb + ng + nq, 1.0, tmp_nbgq + 1, 0, res_d, nb + ng + nq, res_d, nb + ng + nq);
    }
    if (ns > 0) {
        // res_g
        //		dgemv_d(2*ns, 1.0, Z, 0, v, nv, 1.0, gz, nv, res_g, nv);
        dgemv_d(2 * ns, 1.0, Z, 0, v, nv, 2.0, gz, nv, res_g, nv);
        res->obj += 0.5 * ddot(2 * ns, res_g, nv, v, nv);
        daxpy(2 * ns, -1.0, gz, nv, res_g, nv, res_g, nv);

        daxpy(2 * ns, -1.0, lam, 2 * nb + 2 * ng + 2 * nq, res_g, nv, res_g, nv);
        for (ii = 0; ii < nb + ng + nq; ii++) {
            idx = idxs_rev[ii];
            if (idx != -1) {
                VECEL(res_g, nv + idx) -= VECEL(lam, ii);
                VECEL(res_g, nv + ns + idx) -= VECEL(lam, nb + ng + nq + ii);
                // res_d
                VECEL(res_d, ii) -= VECEL(v, nv + idx);
                VECEL(res_d, nb + ng + nq + ii) -= VECEL(v, nv + ns + idx);
            }
        }
        // res_d
        daxpy(2 * ns, -1.0, v, nv, t, 2 * nb + 2 * ng + 2 * nq, res_d, 2 * nb + 2 * ng + 2 * nq);
        daxpy(2 * ns, 1.0, d, 2 * nb + 2 * ng + 2 * nq, res_d, 2 * nb + 2 * ng + 2 * nq, res_d, 2 * nb + 2 * ng + 2 * nq);
    }

    // res b, res g
    if (ne > 0)
        dgemv_nt(ne, nv, -1.0, -1.0, A, 0, 0, v, 0, pi, 0, 1.0, 1.0, b, 0, res_g, 0, res_b, 0, res_g, 0);

    // res_m res_mu
    mu = dvecmuldot(nct, lam, 0, t, 0, res_m, 0);
    daxpy(nct, -1.0, m, 0, res_m, 0, res_m, 0);
    res->res_mu = mu * nct_inv;
}


void d_dense_qcqp_res_compute_inf_norm(struct d_dense_qcqp_res* res) {

    int nv = res->dim->nv;
    int ne = res->dim->ne;
    int nb = res->dim->nb;
    int ng = res->dim->ng;
    int nq = res->dim->nq;
    int ns = res->dim->ns;

    int nvt = nv + 2 * ns;
    int net = ne;
    int nct = 2 * nb + 2 * ng + 2 * nq + 2 * ns;

    struct vec* res_g = res->res_g;
    struct vec* res_b = res->res_b;
    struct vec* res_d = res->res_d;
    struct vec* res_m = res->res_m;

    // compute infinity norm
    dvecnrm_inf(nvt, res_g, 0, res->res_max + 0);
    dvecnrm_inf(net, res_b, 0, res->res_max + 1);
    dvecnrm_inf(nct, res_d, 0, res->res_max + 2);
    dvecnrm_inf(nct, res_m, 0, res->res_max + 3);
}
