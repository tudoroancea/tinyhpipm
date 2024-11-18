#include <stdio.h>
#include <stdlib.h>

#include "tinyhpipm/blas.h"
#include "tinyhpipm/common.h"
#include "tinyhpipm/ocp/d_ocp_qcqp_dim.h"
#include "tinyhpipm/ocp/d_ocp_qcqp_res.h"


hpipm_size_t d_ocp_qcqp_res_memsize(struct d_ocp_qcqp_dim* dim) {

    // loop index
    int ii;

    // extract ocp qp size
    int N = dim->N;
    int* nx = dim->nx;
    int* nu = dim->nu;
    int* nb = dim->nb;
    int* ng = dim->ng;
    int* nq = dim->nq;
    int* ns = dim->ns;

    // compute core qp size
    int nvt = 0;
    int net = 0;
    int nct = 0;
    for (ii = 0; ii < N; ii++) {
        nvt += nx[ii] + nu[ii] + 2 * ns[ii];
        net += nx[ii + 1];
        nct += 2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii];
    }
    nvt += nx[ii] + nu[ii] + 2 * ns[ii];
    nct += 2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii];

    hpipm_size_t size = 0;

    size += 3 * (N + 1) * sizeof(struct vec);  // res_g res_d res_m
    size += 3 * N * sizeof(struct vec);  // res_b

    size += 1 * memsize_vec(nvt);  // res_g
    size += 1 * memsize_vec(net);  // res_b
    size += 2 * memsize_vec(nct);  // res_d res_m

    size = (size + 63) / 64 * 64;  // make multiple of typical cache line size
    size += 1 * 64;  // align once to typical cache line size

    return size;
}


void d_ocp_qcqp_res_create(struct d_ocp_qcqp_dim* dim, struct d_ocp_qcqp_res* res, void* mem) {

    // loop index
    int ii;

    // zero memory (to avoid corrupted memory like e.g. NaN)
    hpipm_size_t memsize = d_ocp_qcqp_res_memsize(dim);
    hpipm_zero_memset(memsize, mem);

    // extract ocp qp size
    int N = dim->N;
    int* nx = dim->nx;
    int* nu = dim->nu;
    int* nb = dim->nb;
    int* ng = dim->ng;
    int* nq = dim->nq;
    int* ns = dim->ns;

    // compute core qp size
    int nvt = 0;
    int net = 0;
    int nct = 0;
    for (ii = 0; ii < N; ii++) {
        nvt += nx[ii] + nu[ii] + 2 * ns[ii];
        net += nx[ii + 1];
        nct += 2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii];
    }
    nvt += nx[ii] + nu[ii] + 2 * ns[ii];
    nct += 2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii];


    // vector struct
    struct vec* sv_ptr = (struct vec*) mem;

    res->res_g = sv_ptr;
    sv_ptr += N + 1;
    res->res_b = sv_ptr;
    sv_ptr += N;
    res->res_d = sv_ptr;
    sv_ptr += N + 1;
    res->res_m = sv_ptr;
    sv_ptr += N + 1;


    // align to typical cache line size
    hpipm_size_t s_ptr = (hpipm_size_t) sv_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;


    // void stuf
    char* c_ptr = (char*) s_ptr;

    create_vec(nvt, res->res_g, c_ptr);
    c_ptr += memsize_vec(nvt);

    create_vec(net, res->res_b, c_ptr);
    c_ptr += memsize_vec(net);

    create_vec(nct, res->res_d, c_ptr);
    c_ptr += memsize_vec(nct);

    create_vec(nct, res->res_m, c_ptr);
    c_ptr += memsize_vec(nct);

    // alias
    //
    c_ptr = (char*) res->res_g->pa;
    for (ii = 0; ii <= N; ii++) {
        create_vec(nu[ii] + nx[ii] + 2 * ns[ii], res->res_g + ii, c_ptr);
        c_ptr += nu[ii] * sizeof(double);
        c_ptr += nx[ii] * sizeof(double);
        c_ptr += ns[ii] * sizeof(double);
        c_ptr += ns[ii] * sizeof(double);
    }
    //
    c_ptr = (char*) res->res_b->pa;
    for (ii = 0; ii < N; ii++) {
        create_vec(nx[ii + 1], res->res_b + ii, c_ptr);
        c_ptr += (nx[ii + 1]) * sizeof(double);
    }
    //
    c_ptr = (char*) res->res_d->pa;
    for (ii = 0; ii <= N; ii++) {
        create_vec(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], res->res_d + ii, c_ptr);
        c_ptr += nb[ii] * sizeof(double);
        c_ptr += ng[ii] * sizeof(double);
        c_ptr += nq[ii] * sizeof(double);
        c_ptr += nb[ii] * sizeof(double);
        c_ptr += ng[ii] * sizeof(double);
        c_ptr += nq[ii] * sizeof(double);
        c_ptr += ns[ii] * sizeof(double);
        c_ptr += ns[ii] * sizeof(double);
    }
    //
    c_ptr = (char*) res->res_m->pa;
    for (ii = 0; ii <= N; ii++) {
        create_vec(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], res->res_m + ii, c_ptr);
        c_ptr += nb[ii] * sizeof(double);
        c_ptr += ng[ii] * sizeof(double);
        c_ptr += nq[ii] * sizeof(double);
        c_ptr += nb[ii] * sizeof(double);
        c_ptr += ng[ii] * sizeof(double);
        c_ptr += nq[ii] * sizeof(double);
        c_ptr += ns[ii] * sizeof(double);
        c_ptr += ns[ii] * sizeof(double);
    }


    res->dim = dim;

    res->memsize = d_ocp_qcqp_res_memsize(dim);


#if defined(RUNTIME_CHECKS)
    if (c_ptr > ((char*) mem) + res->memsize) {
        printf("\ncreate_ocp_qcqp_res: outside memory bounds!\n\n");
        exit(1);
    }
#endif
}


hpipm_size_t d_ocp_qcqp_res_ws_memsize(struct d_ocp_qcqp_dim* dim) {

    // loop index
    int ii;

    // extract ocp qp size
    int N = dim->N;
    int* nx = dim->nx;
    int* nu = dim->nu;
    int* nb = dim->nb;
    int* ng = dim->ng;
    int* nq = dim->nq;
    int* ns = dim->ns;

    // compute core qp size and max size
    int nuM = 0;
    int nxM = 0;
    int nbM = 0;
    int ngM = 0;
    int nqM = 0;
    int nsM = 0;
    for (ii = 0; ii <= N; ii++) {
        nuM = nu[ii] > nuM ? nu[ii] : nuM;
        nxM = nx[ii] > nxM ? nx[ii] : nxM;
        nbM = nb[ii] > nbM ? nb[ii] : nbM;
        ngM = ng[ii] > ngM ? ng[ii] : ngM;
        nqM = nq[ii] > nqM ? nq[ii] : nqM;
        nsM = ns[ii] > nsM ? ns[ii] : nsM;
    }

    hpipm_size_t size = 0;

    size += (5 + 2 * (N + 1)) * sizeof(struct vec);  // 2*tmp_nuxM 2*tmp_nbgqM tmp_nsM q_fun q_adj

    size += 2 * memsize_vec(nuM + nxM);  // 2*tmp_nuxM
    size += 2 * memsize_vec(nbM + ngM + nqM);  // tmp_nbgqM
    size += 1 * memsize_vec(nsM);  // tmp_nsM
    for (ii = 0; ii <= N; ii++) {
        size += 1 * memsize_vec(nq[ii]);  // q_fun
        size += 1 * memsize_vec(nu[ii] + nx[ii]);  // q_adj
    }

    size = (size + 63) / 64 * 64;  // make multiple of typical cache line size
    size += 1 * 64;  // align once to typical cache line size

    return size;
}


void d_ocp_qcqp_res_ws_create(struct d_ocp_qcqp_dim* dim, struct d_ocp_qcqp_res_ws* ws, void* mem) {

    // loop index
    int ii;

    // zero memory (to avoid corrupted memory like e.g. NaN)
    hpipm_size_t memsize = d_ocp_qcqp_res_ws_memsize(dim);
    hpipm_zero_memset(memsize, mem);

    // extract ocp qp size
    int N = dim->N;
    int* nx = dim->nx;
    int* nu = dim->nu;
    int* nb = dim->nb;
    int* ng = dim->ng;
    int* nq = dim->nq;
    int* ns = dim->ns;


    // compute core qp size and max size
    int nuM = 0;
    int nxM = 0;
    int nbM = 0;
    int ngM = 0;
    int nqM = 0;
    int nsM = 0;
    for (ii = 0; ii <= N; ii++) {
        nuM = nu[ii] > nuM ? nu[ii] : nuM;
        nxM = nx[ii] > nxM ? nx[ii] : nxM;
        nbM = nb[ii] > nbM ? nb[ii] : nbM;
        ngM = ng[ii] > ngM ? ng[ii] : ngM;
        nqM = nq[ii] > nqM ? nq[ii] : nqM;
        nsM = ns[ii] > nsM ? ns[ii] : nsM;
    }


    // vector struct
    struct vec* sv_ptr = (struct vec*) mem;

    ws->tmp_nuxM = sv_ptr;
    sv_ptr += 2;
    ws->tmp_nbgqM = sv_ptr;
    sv_ptr += 2;
    ws->tmp_nsM = sv_ptr;
    sv_ptr += 1;
    ws->q_fun = sv_ptr;
    sv_ptr += N + 1;
    ws->q_adj = sv_ptr;
    sv_ptr += N + 1;


    // align to typical cache line size
    hpipm_size_t s_ptr = (hpipm_size_t) sv_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;


    // void stuf
    char* c_ptr = (char*) s_ptr;


    create_vec(nuM + nxM, ws->tmp_nuxM + 0, c_ptr);
    c_ptr += (ws->tmp_nuxM + 0)->memsize;
    create_vec(nuM + nxM, ws->tmp_nuxM + 1, c_ptr);
    c_ptr += (ws->tmp_nuxM + 1)->memsize;

    create_vec(nbM + ngM + nqM, ws->tmp_nbgqM + 0, c_ptr);
    c_ptr += (ws->tmp_nbgqM + 0)->memsize;
    create_vec(nbM + ngM + nqM, ws->tmp_nbgqM + 1, c_ptr);
    c_ptr += (ws->tmp_nbgqM + 1)->memsize;

    create_vec(nsM, ws->tmp_nsM, c_ptr);
    c_ptr += (ws->tmp_nsM)->memsize;

    for (ii = 0; ii <= N; ii++) {
        create_vec(nq[ii], ws->q_fun + ii, c_ptr);
        c_ptr += (ws->q_fun + ii)->memsize;
    }

    for (ii = 0; ii <= N; ii++) {
        create_vec(nu[ii] + nx[ii], ws->q_adj + ii, c_ptr);
        c_ptr += (ws->q_adj + ii)->memsize;
    }

    ws->use_q_fun = 0;
    ws->use_q_adj = 0;

    ws->memsize = memsize;  // d_ocp_qcqp_res_ws_memsize(dim);


#if defined(RUNTIME_CHECKS)
    if (c_ptr > ((char*) mem) + ws->memsize) {
        printf("\ncreate_ocp_qcqp_res_workspace: outside memory bounds!\n\n");
        exit(1);
    }
#endif
}


void d_ocp_qcqp_res_compute(struct d_ocp_qcqp* qp, struct d_ocp_qcqp_sol* qp_sol, struct d_ocp_qcqp_res* res, struct d_ocp_qcqp_res_ws* ws) {

    // loop index
    int ii, jj, idx;

    //
    int N = qp->dim->N;
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;
    int* ns = qp->dim->ns;

    int nct = 0;
    for (ii = 0; ii <= N; ii++)
        nct += 2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii];

    double nct_inv = 1.0 / nct;

    struct mat* BAbt = qp->BAbt;
    struct mat* RSQrq = qp->RSQrq;
    struct mat* DCt = qp->DCt;
    struct vec* b = qp->b;
    struct vec* rqz = qp->rqz;
    struct vec* d = qp->d;
    struct vec* m = qp->m;
    int** idxb = qp->idxb;
    struct vec* Z = qp->Z;
    int** idxs_rev = qp->idxs_rev;
    struct mat** Hq = qp->Hq;

    struct vec* ux = qp_sol->ux;
    struct vec* pi = qp_sol->pi;
    struct vec* lam = qp_sol->lam;
    struct vec* t = qp_sol->t;

    struct vec* res_g = res->res_g;
    struct vec* res_b = res->res_b;
    struct vec* res_d = res->res_d;
    struct vec* res_m = res->res_m;

    struct vec* tmp_nuxM = ws->tmp_nuxM;
    struct vec* tmp_nbgqM = ws->tmp_nbgqM;
    struct vec* tmp_nsM = ws->tmp_nsM;

    struct vec* q_fun = ws->q_fun;
    struct vec* q_adj = ws->q_adj;

    int nx0, nx1, nu0, nu1, nb0, ng0, nq0, ns0;

    //
    double tmp;
    double mu = 0.0;
    res->obj = 0.0;

    // loop over stages
    for (ii = 0; ii <= N; ii++) {

        nx0 = nx[ii];
        nu0 = nu[ii];
        nb0 = nb[ii];
        ng0 = ng[ii];
        nq0 = nq[ii];
        ns0 = ns[ii];

        //		dsymv_l(nu0+nx0, 1.0, RSQrq+ii, 0, 0, ux+ii, 0, 1.0, rqz+ii, 0, res_g+ii, 0);
        dsymv_l(nu0 + nx0, 1.0, RSQrq + ii, 0, 0, ux + ii, 0, 2.0, rqz + ii, 0, res_g + ii, 0);
        res->obj += 0.5 * ddot(nu0 + nx0, res_g + ii, 0, ux + ii, 0);
        daxpy(nu0 + nx0, -1.0, rqz + ii, 0, res_g + ii, 0, res_g + ii, 0);

        if (ii > 0)
            daxpy(nx0, -1.0, pi + (ii - 1), 0, res_g + ii, nu0, res_g + ii, nu0);

        if (nb0 + ng0 + nq0 > 0) {
            daxpy(nb0 + ng0 + nq0, -1.0, lam + ii, 0, lam + ii, nb0 + ng0 + nq0, tmp_nbgqM + 0, 0);
            daxpy(2 * nb0 + 2 * ng0 + 2 * nq0, 1.0, d + ii, 0, t + ii, 0, res_d + ii, 0);
            // box
            if (nb0 > 0) {
                dvecad_sp(nb0, 1.0, tmp_nbgqM + 0, 0, idxb[ii], res_g + ii, 0);
                dvecex_sp(nb0, 1.0, idxb[ii], ux + ii, 0, tmp_nbgqM + 1, 0);
            }
            // general
            if (ng0 > 0) {
                dgemv_nt(nu0 + nx0, ng0, 1.0, 1.0, DCt + ii, 0, 0, tmp_nbgqM + 0, nb[ii], ux + ii, 0, 1.0, 0.0, res_g + ii, 0, tmp_nbgqM + 1, nb0, res_g + ii, 0, tmp_nbgqM + 1, nb0);
            }
            // quadratic
            if (nq0 > 0) {
                //				daxpy(nq0,  1.0, d, 2*nb0+2*ng0+2*ns, t, 2*nb+2*ng+2*ns, res_d, 2*nb+2*ng+2*ns);
                if (ws->use_q_fun & ws->use_q_adj) {
                    dveccp(nq0, ws->q_fun + ii, 0, tmp_nbgqM + 1, nb0 + ng0);
                    daxpy(nu0 + nx0, 1.0, ws->q_adj + ii, 0, res_g + ii, 0, res_g + ii, 0);
                } else {
                    for (jj = 0; jj < nq0; jj++) {
                        dsymv_l(nu0 + nx0, 1.0, &Hq[ii][jj], 0, 0, ux + ii, 0, 0.0, tmp_nuxM, 0, tmp_nuxM, 0);
                        tmp = VECEL(tmp_nbgqM + 0, nb0 + ng0 + jj);
                        daxpy(nu0 + nx0, tmp, tmp_nuxM, 0, res_g + ii, 0, res_g + ii, 0);
                        dcolex(nu0 + nx0, DCt + ii, 0, ng0 + jj, tmp_nuxM + 1, 0);
                        daxpy(nu0 + nx0, tmp, tmp_nuxM + 1, 0, res_g + ii, 0, res_g + ii, 0);
                        daxpy(nu0 + nx0, 0.5, tmp_nuxM, 0, tmp_nuxM + 1, 0, tmp_nuxM, 0);
                        tmp = ddot(nu0 + nx0, tmp_nuxM, 0, ux + ii, 0);
                        VECEL(tmp_nbgqM + 1, nb0 + ng0 + jj) = tmp;
                    }
                }
            }
            daxpy(nb0 + ng0 + nq0, -1.0, tmp_nbgqM + 1, 0, res_d + ii, 0, res_d + ii, 0);
            daxpy(nb0 + ng0 + nq0, 1.0, tmp_nbgqM + 1, 0, res_d + ii, nb0 + ng0 + nq0, res_d + ii, nb0 + ng0 + nq0);
        }
        if (ns0 > 0) {
            // res_g
            //			dgemv_d(2*ns0, 1.0, Z+ii, 0, ux+ii, nu0+nx0, 1.0, rqz+ii, nu0+nx0, res_g+ii, nu0+nx0);
            dgemv_d(2 * ns0, 1.0, Z + ii, 0, ux + ii, nu0 + nx0, 2.0, rqz + ii, nu0 + nx0, res_g + ii, nu0 + nx0);
            res->obj += 0.5 * ddot(2 * ns0, res_g + ii, nu0 + nx0, ux + ii, nu0 + nx0);
            daxpy(2 * ns0, -1.0, rqz + ii, nu0 + nx0, res_g + ii, nu0 + nx0, res_g + ii, nu0 + nx0);

            daxpy(2 * ns0, -1.0, lam + ii, 2 * nb0 + 2 * ng0 + 2 * nq0, res_g + ii, nu0 + nx0, res_g + ii, nu0 + nx0);
            for (jj = 0; jj < nb0 + ng0 + nq0; jj++) {
                idx = idxs_rev[ii][jj];
                if (idx != -1) {
                    VECEL(res_g + ii, nu0 + nx0 + idx) -= VECEL(lam + ii, jj);
                    VECEL(res_g + ii, nu0 + nx0 + ns0 + idx) -= VECEL(lam + ii, nb0 + ng0 + nq0 + jj);
                    // res_d
                    VECEL(res_d + ii, jj) -= VECEL(ux + ii, nu0 + nx0 + idx);
                    VECEL(res_d + ii, nb0 + ng0 + nq0 + jj) -= VECEL(ux + ii, nu0 + nx0 + ns0 + idx);
                }
            }
            // res_d
            daxpy(2 * ns0, -1.0, ux + ii, nu0 + nx0, t + ii, 2 * nb0 + 2 * ng0 + 2 * nq0, res_d + ii, 2 * nb0 + 2 * ng0 + 2 * nq0);
            daxpy(2 * ns0, 1.0, d + ii, 2 * nb0 + 2 * ng0 + 2 * nq0, res_d + ii, 2 * nb0 + 2 * ng0 + 2 * nq0, res_d + ii, 2 * nb0 + 2 * ng0 + 2 * nq0);
        }

        if (ii < N) {

            nu1 = nu[ii + 1];
            nx1 = nx[ii + 1];

            daxpy(nx1, -1.0, ux + (ii + 1), nu1, b + ii, 0, res_b + ii, 0);

            dgemv_nt(nu0 + nx0, nx1, 1.0, 1.0, BAbt + ii, 0, 0, pi + ii, 0, ux + ii, 0, 1.0, 1.0, res_g + ii, 0, res_b + ii, 0, res_g + ii, 0, res_b + ii, 0);
        }

        mu += dvecmuldot(2 * nb0 + 2 * ng0 + 2 * nq0 + 2 * ns0, lam + ii, 0, t + ii, 0, res_m + ii, 0);
        daxpy(2 * nb0 + 2 * ng0 + 2 * nq0 + 2 * ns0, -1.0, m + ii, 0, res_m + ii, 0, res_m + ii, 0);
    }

    res->res_mu = mu * nct_inv;
}


void d_ocp_qcqp_res_compute_inf_norm(struct d_ocp_qcqp_res* res) {

    struct d_ocp_qcqp_dim* dim = res->dim;
    int N = dim->N;
    int* nx = dim->nx;
    int* nu = dim->nu;
    int* nb = dim->nb;
    int* ng = dim->ng;
    int* nq = dim->nq;
    int* ns = dim->ns;

    int ii;

    int nv = 0;
    int ne = 0;
    int nc = 0;

    for (ii = 0; ii < N; ii++) {
        nv += nu[ii] + nx[ii] + 2 * ns[ii];
        ne += nx[ii + 1];
        nc += 2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii];
    }
    ii = N;
    nv += nu[ii] + nx[ii] + 2 * ns[ii];
    nc += 2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii];

    // compute infinity norm
    dvecnrm_inf(nv, res->res_g, 0, res->res_max + 0);
    dvecnrm_inf(ne, res->res_b, 0, res->res_max + 1);
    dvecnrm_inf(nc, res->res_d, 0, res->res_max + 2);
    dvecnrm_inf(nc, res->res_m, 0, res->res_max + 3);
}


void d_ocp_qcqp_res_get_max_res_stat(struct d_ocp_qcqp_res* res, double* value) {

    *value = res->res_max[0];
}


void d_ocp_qcqp_res_get_max_res_eq(struct d_ocp_qcqp_res* res, double* value) {

    *value = res->res_max[1];
}


void d_ocp_qcqp_res_get_max_res_ineq(struct d_ocp_qcqp_res* res, double* value) {

    *value = res->res_max[2];
}


void d_ocp_qcqp_res_get_max_res_comp(struct d_ocp_qcqp_res* res, double* value) {

    *value = res->res_max[3];
}
