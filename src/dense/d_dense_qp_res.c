#include <stdio.h>
#include <stdlib.h>

#include "tinyhpipm/blas.h"
#include "tinyhpipm/dense/d_dense_qp.h"
#include "tinyhpipm/dense/d_dense_qp_dim.h"
#include "tinyhpipm/dense/d_dense_qp_res.h"
#include "tinyhpipm/dense/d_dense_qp_sol.h"


hpipm_size_t d_dense_qp_res_memsize(struct d_dense_qp_dim* dim) {

    // loop index
    int ii;

    // extract ocp qp size
    int nv = dim->nv;
    int ne = dim->ne;
    int nb = dim->nb;
    int ng = dim->ng;
    int ns = dim->ns;

    hpipm_size_t size = 0;

    size += 4 * sizeof(struct vec);  // res_g res_b res_d res_m

    size += 1 * memsize_vec(nv + 2 * ns);  // res_g
    size += 1 * memsize_vec(ne);  // res_b
    size += 2 * memsize_vec(2 * nb + 2 * ng + 2 * ns);  // res_d res_m

    size = (size + 63) / 64 * 64;  // make multiple of typical cache line size
    size += 1 * 64;  // align once to typical cache line size

    return size;
}


void d_dense_qp_res_create(struct d_dense_qp_dim* dim, struct d_dense_qp_res* res, void* mem) {

    // loop index
    int ii;

    // zero memory (to avoid corrupted memory like e.g. NaN)
    hpipm_size_t memsize = d_dense_qp_res_memsize(dim);
    hpipm_zero_memset(memsize, mem);

    // extract ocp qp size
    int nv = dim->nv;
    int ne = dim->ne;
    int nb = dim->nb;
    int ng = dim->ng;
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

    create_vec(2 * nb + 2 * ng + 2 * ns, res->res_d, c_ptr);
    c_ptr += (res->res_d)->memsize;

    create_vec(2 * nb + 2 * ng + 2 * ns, res->res_m, c_ptr);
    c_ptr += (res->res_m)->memsize;

    res->dim = dim;

    res->memsize = d_dense_qp_res_memsize(dim);


#if defined(RUNTIME_CHECKS)
    if (c_ptr > ((char*) mem) + res->memsize) {
        printf("\ncreate_dense_qp_res: outsize memory bounds!\n\n");
        exit(1);
    }
#endif
}


hpipm_size_t d_dense_qp_res_ws_memsize(struct d_dense_qp_dim* dim) {

    // loop index
    int ii;

    // extract ocp qp size
    int nv = dim->nv;
    int ne = dim->ne;
    int nb = dim->nb;
    int ng = dim->ng;
    int ns = dim->ns;

    hpipm_size_t size = 0;

    size += 3 * sizeof(struct vec);  // 2*tmp_nbg tmp_ns

    size += 2 * memsize_vec(nb + ng);  // tmp_nbg
    size += 1 * memsize_vec(ns);  // tmp_ns

    size = (size + 63) / 64 * 64;  // make multiple of typical cache line size
    size += 1 * 64;  // align once to typical cache line size

    return size;
}


void d_dense_qp_res_ws_create(struct d_dense_qp_dim* dim, struct d_dense_qp_res_ws* ws, void* mem) {

    // loop index
    int ii;

    // zero memory (to avoid corrupted memory like e.g. NaN)
    hpipm_size_t memsize = d_dense_qp_res_ws_memsize(dim);
    hpipm_size_t memsize_m8 = memsize / 8;  // sizeof(double) is 8
    //	hpipm_size_t memsize_r8 = memsize - 8*memsize_m8;
    double* double_ptr = mem;
    // XXX exploit that it is multiple of 64 bytes !!!!!
    if (memsize_m8 > 7)
        for (ii = 0; ii < memsize_m8 - 7; ii += 8) {
            double_ptr[ii + 0] = 0.0;
            double_ptr[ii + 1] = 0.0;
            double_ptr[ii + 2] = 0.0;
            double_ptr[ii + 3] = 0.0;
            double_ptr[ii + 4] = 0.0;
            double_ptr[ii + 5] = 0.0;
            double_ptr[ii + 6] = 0.0;
            double_ptr[ii + 7] = 0.0;
        }
    //	for(; ii<memsize_m8; ii++)
    //		{
    //		double_ptr[ii] = 0.0;
    //		}
    //	char *char_ptr = (char *) (&double_ptr[ii]);
    //	for(ii=0; ii<memsize_r8; ii++)
    //		{
    //		char_ptr[ii] = 0;
    //		}

    // extract ocp qp size
    int nv = dim->nv;
    int ne = dim->ne;
    int nb = dim->nb;
    int ng = dim->ng;
    int ns = dim->ns;


    // vector struct
    struct vec* sv_ptr = (struct vec*) mem;

    ws->tmp_nbg = sv_ptr;
    sv_ptr += 2;
    ws->tmp_ns = sv_ptr;
    sv_ptr += 1;


    // align to typical cache line size
    hpipm_size_t s_ptr = (hpipm_size_t) sv_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;


    // void stuf
    char* c_ptr = (char*) s_ptr;


    create_vec(nb + ng, ws->tmp_nbg + 0, c_ptr);
    c_ptr += (ws->tmp_nbg + 0)->memsize;

    create_vec(nb + ng, ws->tmp_nbg + 1, c_ptr);
    c_ptr += (ws->tmp_nbg + 1)->memsize;

    create_vec(ns, ws->tmp_ns + 0, c_ptr);
    c_ptr += (ws->tmp_ns + 0)->memsize;

    ws->memsize = d_dense_qp_res_ws_memsize(dim);


#if defined(RUNTIME_CHECKS)
    if (c_ptr > ((char*) mem) + ws->memsize) {
        printf("\ncreate_dense_qp_res_workspace: outsize memory bounds!\n\n");
        exit(1);
    }
#endif
}


void d_dense_qp_res_compute(struct d_dense_qp* qp, struct d_dense_qp_sol* qp_sol, struct d_dense_qp_res* res, struct d_dense_qp_res_ws* ws) {

    int ii, idx;

    int nv = qp->dim->nv;
    int ne = qp->dim->ne;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;
    int ns = qp->dim->ns;

    int nvt = nv + 2 * ns;
    int net = ne;
    int nct = 2 * nb + 2 * ng + 2 * ns;

    // TODO use nc_mask_inv from cws if available !!!!!
    double nct_inv = 1.0 / nct;

    struct mat* Hg = qp->Hv;
    struct mat* A = qp->A;
    struct mat* Ct = qp->Ct;
    struct vec* gz = qp->gz;
    struct vec* b = qp->b;
    struct vec* d = qp->d;
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

    struct vec* tmp_nbg = ws->tmp_nbg;
    struct vec* tmp_ns = ws->tmp_ns;

    double mu, tmp;
    res->obj = 0.0;

    // res g
    //	dsymv_l(nv, 1.0, Hg, 0, 0, v, 0, 1.0, gz, 0, res_g, 0);
    dsymv_l(nv, 1.0, Hg, 0, 0, v, 0, 2.0, gz, 0, res_g, 0);
    res->obj += 0.5 * ddot(nv, res_g, 0, v, 0);
    daxpy(nv, -1.0, gz, 0, res_g, 0, res_g, 0);

    if (nb + ng > 0) {
        daxpy(nb + ng, -1.0, lam, 0, lam, nb + ng, tmp_nbg + 0, 0);
        //		daxpy(nb+ng,  1.0, d, 0, t, 0, res_d, 0);
        //		daxpy(nb+ng,  1.0, d, nb+ng, t, nb+ng, res_d, nb+ng);
        daxpy(2 * nb + 2 * ng, 1.0, d, 0, t, 0, res_d, 0);
        // box
        if (nb > 0) {
            dvecad_sp(nb, 1.0, tmp_nbg + 0, 0, idxb, res_g, 0);
            dvecex_sp(nb, 1.0, idxb, v, 0, tmp_nbg + 1, 0);
        }
        // general
        if (ng > 0) {
            dgemv_nt(nv, ng, 1.0, 1.0, Ct, 0, 0, tmp_nbg + 0, nb, v, 0, 1.0, 0.0, res_g, 0, tmp_nbg + 1, nb, res_g, 0, tmp_nbg + 1, nb);
        }
        daxpy(nb + ng, -1.0, tmp_nbg + 1, 0, res_d, 0, res_d, 0);
        daxpy(nb + ng, 1.0, tmp_nbg + 1, 0, res_d, nb + ng, res_d, nb + ng);
    }
    if (ns > 0) {
        // res_g
        //		dgemv_d(2*ns, 1.0, Z, 0, v, nv, 1.0, gz, nv, res_g, nv);
        dgemv_d(2 * ns, 1.0, Z, 0, v, nv, 2.0, gz, nv, res_g, nv);
        res->obj += 0.5 * ddot(2 * ns, res_g, nv, v, nv);
        daxpy(2 * ns, -1.0, gz, nv, res_g, nv, res_g, nv);

        daxpy(2 * ns, -1.0, lam, 2 * nb + 2 * ng, res_g, nv, res_g, nv);
        for (ii = 0; ii < nb + ng; ii++) {
            idx = idxs_rev[ii];
            if (idx != -1) {
                VECEL(res_g, nv + idx) -= VECEL(lam, ii);
                VECEL(res_g, nv + ns + idx) -= VECEL(lam, nb + ng + ii);
                // res_d
                VECEL(res_d, ii) -= VECEL(v, nv + idx);
                VECEL(res_d, nb + ng + ii) -= VECEL(v, nv + ns + idx);
            }
        }
        // res_d
        daxpy(2 * ns, -1.0, v, nv, t, 2 * nb + 2 * ng, res_d, 2 * nb + 2 * ng);
        daxpy(2 * ns, 1.0, d, 2 * nb + 2 * ng, res_d, 2 * nb + 2 * ng, res_d, 2 * nb + 2 * ng);
    }

    // res b, res g
    if (ne > 0)
        dgemv_nt(ne, nv, -1.0, -1.0, A, 0, 0, v, 0, pi, 0, 1.0, 1.0, b, 0, res_g, 0, res_b, 0, res_g, 0);

    // res_m res_mu
    mu = dvecmuldot(nct, lam, 0, t, 0, res_m, 0);
    daxpy(nct, -1.0, m, 0, res_m, 0, res_m, 0);
    // TODO use nc_mask_inv from cws if available !!!!!
    res->res_mu = mu * nct_inv;
}


void d_dense_qp_res_compute_lin(struct d_dense_qp* qp, struct d_dense_qp_sol* qp_sol, struct d_dense_qp_sol* qp_step, struct d_dense_qp_res* res, struct d_dense_qp_res_ws* ws) {

    int ii, idx;

    int nv = qp->dim->nv;
    int ne = qp->dim->ne;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;
    int ns = qp->dim->ns;

    int nvt = nv + 2 * ns;
    int net = ne;
    int nct = 2 * nb + 2 * ng + 2 * ns;

    double nct_inv = 1.0 / nct;

    struct mat* Hg = qp->Hv;
    struct mat* A = qp->A;
    struct mat* Ct = qp->Ct;
    struct vec* gz = qp->gz;
    struct vec* b = qp->b;
    struct vec* d = qp->d;
    struct vec* m = qp->m;
    int* idxb = qp->idxb;
    struct vec* Z = qp->Z;
    int* idxs_rev = qp->idxs_rev;

    struct vec* v = qp_step->v;
    struct vec* pi = qp_step->pi;
    struct vec* lam = qp_step->lam;
    struct vec* t = qp_step->t;

    struct vec* Lam = qp_sol->lam;
    struct vec* T = qp_sol->t;

    struct vec* res_g = res->res_g;
    struct vec* res_b = res->res_b;
    struct vec* res_d = res->res_d;
    struct vec* res_m = res->res_m;

    struct vec* tmp_nbg = ws->tmp_nbg;
    struct vec* tmp_ns = ws->tmp_ns;

    double mu, tmp;

    // res g
    dsymv_l(nv, 1.0, Hg, 0, 0, v, 0, 1.0, gz, 0, res_g, 0);

    if (nb + ng > 0) {
        daxpy(nb + ng, -1.0, lam, 0, lam, nb + ng, tmp_nbg + 0, 0);
        //		daxpy(nb+ng,  1.0, d, 0, t, 0, res_d, 0);
        //		daxpy(nb+ng,  1.0, d, nb+ng, t, nb+ng, res_d, nb+ng);
        daxpy(2 * nb + 2 * ng, 1.0, d, 0, t, 0, res_d, 0);
        // box
        if (nb > 0) {
            dvecad_sp(nb, 1.0, tmp_nbg + 0, 0, idxb, res_g, 0);
            dvecex_sp(nb, 1.0, idxb, v, 0, tmp_nbg + 1, 0);
        }
        // general
        if (ng > 0) {
            dgemv_nt(nv, ng, 1.0, 1.0, Ct, 0, 0, tmp_nbg + 0, nb, v, 0, 1.0, 0.0, res_g, 0, tmp_nbg + 1, nb, res_g, 0, tmp_nbg + 1, nb);
        }
        daxpy(nb + ng, -1.0, tmp_nbg + 1, 0, res_d, 0, res_d, 0);
        daxpy(nb + ng, 1.0, tmp_nbg + 1, 0, res_d, nb + ng, res_d, nb + ng);
    }
    if (ns > 0) {
        // res_g
        dgemv_d(2 * ns, 1.0, Z, 0, v, nv, 1.0, gz, nv, res_g, nv);
        daxpy(2 * ns, -1.0, lam, 2 * nb + 2 * ng, res_g, nv, res_g, nv);
        for (ii = 0; ii < nb + ng; ii++) {
            idx = idxs_rev[ii];
            if (idx != -1) {
                VECEL(res_g, nv + idx) -= VECEL(lam, ii);
                VECEL(res_g, nv + ns + idx) -= VECEL(lam, nb + ng + ii);
                // res_d
                VECEL(res_d, ii) -= VECEL(v, nv + idx);
                VECEL(res_d, nb + ng + ii) -= VECEL(v, nv + ns + idx);
            }
        }
        // res_d
        daxpy(2 * ns, -1.0, v, nv, t, 2 * nb + 2 * ng, res_d, 2 * nb + 2 * ng);
        daxpy(2 * ns, 1.0, d, 2 * nb + 2 * ng, res_d, 2 * nb + 2 * ng, res_d, 2 * nb + 2 * ng);
    }

    // res b, res g
    if (ne > 0)
        dgemv_nt(ne, nv, -1.0, -1.0, A, 0, 0, v, 0, pi, 0, 1.0, 1.0, b, 0, res_g, 0, res_b, 0, res_g, 0);

    // res_m res_mu
    //	dveccpsc(nct, -1.0, m, 0, res_m, 0);
    dveccp(nct, m, 0, res_m, 0);  // TODO scale by -1 ?????
    dvecmulacc(nct, Lam, 0, t, 0, res_m, 0);
    dvecmulacc(nct, lam, 0, T, 0, res_m, 0);
    //	for(ii=0; ii<nct; ii++) (res_m->pa)[ii] += 1e-3;
    //	mu = dvecmuldot(nct, lam, 0, t, 0, res_m, 0);
    //	daxpy(nct, -1.0, m, 0, res_m, 0, res_m, 0);
    //	res->res_mu = mu*nct_inv;
    // TODO use nc_mask_inv from cws if available !!!!!
}


void d_dense_qp_res_compute_inf_norm(struct d_dense_qp_res* res) {

    int nv = res->dim->nv;
    int ne = res->dim->ne;
    int nb = res->dim->nb;
    int ng = res->dim->ng;
    int ns = res->dim->ns;

    int nvt = nv + 2 * ns;
    int net = ne;
    int nct = 2 * nb + 2 * ng + 2 * ns;

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


void d_dense_qp_res_get_all(struct d_dense_qp_res* res, double* res_g, double* res_ls, double* res_us, double* res_b, double* res_d_lb, double* res_d_ub, double* res_d_lg, double* res_d_ug, double* res_d_ls, double* res_d_us, double* res_m_lb, double* res_m_ub, double* res_m_lg, double* res_m_ug, double* res_m_ls, double* res_m_us) {

    int nv = res->dim->nv;
    int ne = res->dim->ne;
    int nb = res->dim->nb;
    int ng = res->dim->ng;
    int ns = res->dim->ns;

    unpack_vec(nv, res->res_g, 0, res_g, 1);
    unpack_vec(ns, res->res_g, nv, res_ls, 1);
    unpack_vec(ns, res->res_g, nv + ns, res_us, 1);

    unpack_vec(ne, res->res_b, 0, res_b, 1);
    unpack_vec(nb, res->res_d, 0, res_d_lb, 1);
    unpack_vec(ng, res->res_d, nb, res_d_lg, 1);
    unpack_vec(nb, res->res_d, nb + ng, res_d_ub, 1);
    unpack_vec(ng, res->res_d, 2 * nb + ng, res_d_ug, 1);
    unpack_vec(ns, res->res_d, 2 * nb + 2 * ng, res_d_ls, 1);
    unpack_vec(ns, res->res_d, 2 * nb + 2 * ng + ns, res_d_us, 1);
    unpack_vec(nb, res->res_m, 0, res_m_lb, 1);
    unpack_vec(ng, res->res_m, nb, res_m_lg, 1);
    unpack_vec(nb, res->res_m, nb + ng, res_m_ub, 1);
    unpack_vec(ng, res->res_m, 2 * nb + ng, res_m_ug, 1);
    unpack_vec(ns, res->res_m, 2 * nb + 2 * ng, res_m_ls, 1);
    unpack_vec(ns, res->res_m, 2 * nb + 2 * ng + ns, res_m_us, 1);
}
