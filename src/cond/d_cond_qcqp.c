#include <stdio.h>
#include <stdlib.h>

#include "tinyhpipm/blas.h"
#include "tinyhpipm/common.h"

#include "tinyhpipm/cond/d_cond.h"
#include "tinyhpipm/cond/d_cond_aux.h"
#include "tinyhpipm/cond/d_cond_qcqp.h"
#include "tinyhpipm/dense/d_dense_qcqp.h"
#include "tinyhpipm/dense/d_dense_qcqp_dim.h"
#include "tinyhpipm/dense/d_dense_qcqp_sol.h"
#include "tinyhpipm/ocp/d_ocp_qcqp.h"
#include "tinyhpipm/ocp/d_ocp_qcqp_dim.h"
#include "tinyhpipm/ocp/d_ocp_qcqp_sol.h"

void d_cond_qcqp_compute_dim(struct d_ocp_qcqp_dim* ocp_dim, struct d_dense_qcqp_dim* dense_dim) {
    int N = ocp_dim->N;
    int* nx = ocp_dim->nx;
    int* nu = ocp_dim->nu;
    int* nbx = ocp_dim->nbx;
    int* nbu = ocp_dim->nbu;
    int* ng = ocp_dim->ng;
    int* nq = ocp_dim->nq;
    int* ns = ocp_dim->ns;
    int* nsbx = ocp_dim->nsbx;
    int* nsbu = ocp_dim->nsbu;
    int* nsg = ocp_dim->nsg;
    int* nsq = ocp_dim->nsq;

    int ii;

    int nvc = 0;
    int nec = 0;
    int nbc = 0;
    int ngc = 0;
    int nqc = 0;
    int nsc = 0;
    int nsbc = 0;
    int nsgc = 0;
    int nsqc = 0;

    // first stage
    nvc += nx[0] + nu[0];
    nbc += nbx[0] + nbu[0];
    ngc += ng[0];
    nqc += nq[0];
    nsc += ns[0];
    nsbc += nsbx[0] + nsbu[0];
    nsgc += nsg[0];
    nsqc += nsq[0];
    // remaining stages
    for (ii = 1; ii <= N; ii++) {
        nvc += nu[ii];
        nbc += nbu[ii];
        ngc += nbx[ii] + ng[ii];
        nqc += nq[ii];
        nsc += ns[ii];
        nsbc += nsbu[ii];
        nsgc += nsbx[ii] + nsg[ii];
        nsqc += nsq[ii];
    }

    // XXX must use setters to correctly set qp ones too !
    d_dense_qcqp_dim_set_nv(nvc, dense_dim);
    d_dense_qcqp_dim_set_ne(nec, dense_dim);
    d_dense_qcqp_dim_set_nb(nbc, dense_dim);
    d_dense_qcqp_dim_set_ng(ngc, dense_dim);
    d_dense_qcqp_dim_set_nq(nqc, dense_dim);
    d_dense_qcqp_dim_set_ns(nsc, dense_dim);
    d_dense_qcqp_dim_set_nsb(nsbc, dense_dim);
    d_dense_qcqp_dim_set_nsg(nsgc, dense_dim);
    d_dense_qcqp_dim_set_nsq(nsqc, dense_dim);
}


hpipm_size_t d_cond_qcqp_arg_memsize() {

    hpipm_size_t size = 0;

    size += 1 * sizeof(struct d_cond_qp_arg);
    size += 1 * d_cond_qp_arg_memsize();

    size = (size + 63) / 64 * 64;  // make multiple of typical cache line size
    size += 1 * 64;  // align once to typical cache line size

    return size;
}


void d_cond_qcqp_arg_create(struct d_cond_qcqp_arg* cond_arg, void* mem) {

    // loop index
    int ii;

    // zero memory (to avoid corrupted memory like e.g. NaN)
    hpipm_size_t memsize = d_cond_qcqp_arg_memsize();
    hpipm_zero_memset(memsize, mem);

    struct d_cond_qp_arg* arg_ptr = mem;

    cond_arg->qp_arg = arg_ptr;
    arg_ptr += 1;

    // align to typical cache line size
    hpipm_size_t s_ptr = (hpipm_size_t) arg_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;

    // void
    char* c_ptr = (char*) s_ptr;

    d_cond_qp_arg_create(cond_arg->qp_arg, c_ptr);
    c_ptr += cond_arg->qp_arg->memsize;


    cond_arg->memsize = d_cond_qcqp_arg_memsize();

#if defined(RUNTIME_CHECKS)
    if (c_ptr > ((char*) mem) + cond_arg->memsize) {
        printf("\nerror: d_cond_qcqp_arg_create: outside memory bounds!\n\n");
        exit(1);
    }
#endif
}


void d_cond_qcqp_arg_set_default(struct d_cond_qcqp_arg* cond_arg) {

    cond_arg->cond_last_stage = 1;  // condense last stage
    cond_arg->comp_prim_sol = 1;  // compute primal solution (v)
    cond_arg->comp_dual_sol_eq = 1;  // compute dual solution equality constr (pi)
    cond_arg->comp_dual_sol_ineq = 1;  // compute dual solution inequality constr (lam t)
    cond_arg->square_root_alg = 1;  // square root algorithm (faster but requires RSQ>0)

    // set arg of qp struct
    cond_arg->qp_arg->cond_last_stage = cond_arg->cond_last_stage;
    cond_arg->qp_arg->comp_prim_sol = cond_arg->comp_prim_sol;
    cond_arg->qp_arg->comp_dual_sol_eq = cond_arg->comp_dual_sol_eq;
    cond_arg->qp_arg->comp_dual_sol_ineq = cond_arg->comp_dual_sol_ineq;
    cond_arg->qp_arg->square_root_alg = cond_arg->square_root_alg;
}


void d_cond_qcqp_arg_set_ric_alg(int ric_alg, struct d_cond_qcqp_arg* cond_arg) {

    cond_arg->square_root_alg = ric_alg;

    // set arg of qp struct
    cond_arg->qp_arg->square_root_alg = cond_arg->square_root_alg;
}


void d_cond_qcqp_arg_set_cond_last_stage(int cond_last_stage, struct d_cond_qcqp_arg* cond_arg) {

    cond_arg->cond_last_stage = cond_last_stage;

    // set arg of qp struct
    cond_arg->qp_arg->cond_last_stage = cond_arg->cond_last_stage;
}


hpipm_size_t d_cond_qcqp_ws_memsize(struct d_ocp_qcqp_dim* ocp_dim, struct d_cond_qcqp_arg* cond_arg) {

    hpipm_size_t size = 0;

    size += 1 * sizeof(struct d_cond_qp_ws);
    size += 1 * d_cond_qp_ws_memsize(ocp_dim->qp_dim, cond_arg->qp_arg);

    int ii;

    int N = ocp_dim->N;
    int* nx = ocp_dim->nx;
    int* nu = ocp_dim->nu;
    int* nb = ocp_dim->nb;
    int* ng = ocp_dim->ng;
    int* nq = ocp_dim->nq;

    // compute core qp size and max size
    //	int nvt = 0;
    //	int net = 0;
    //	int nbt = 0;
    //	int ngt = 0;
    int nxM = 0;
    int nuM = 0;
    //	int nbM = 0;
    //	int ngM = 0;

    for (ii = 0; ii < N; ii++) {
        //		nvt += nx[ii]+nu[ii];
        //		net += nx[ii+1];
        //		nbt += nb[ii];
        //		ngt += ng[ii];
        nxM = nx[ii] > nxM ? nx[ii] : nxM;
        nuM = nu[ii] > nuM ? nu[ii] : nuM;
        //		nbM = nb[ii]>nbM ? nb[ii] : nbM;
        //		ngM = ng[ii]>ngM ? ng[ii] : ngM;
    }
    ii = N;
    //	nvt += nx[ii]+nu[ii];
    //	nbt += nb[ii];
    //	ngt += ng[ii];
    nxM = nx[ii] > nxM ? nx[ii] : nxM;
    nuM = nu[ii] > nuM ? nu[ii] : nuM;
    //	nbM = nb[ii]>nbM ? nb[ii] : nbM;
    //	ngM = ng[ii]>ngM ? ng[ii] : ngM;

    int nvc = 0;
    nvc += nx[0] + nu[0];
    for (ii = 1; ii <= N; ii++)
        nvc += nu[ii];

    size += 3 * (N + 1) * sizeof(struct mat);  // hess_array GammaQ tmp_DCt
    size += 1 * sizeof(struct mat);  // zero_hess
    size += 1 * (N + 1) * sizeof(struct vec);  // grad_array
    size += 3 * sizeof(struct vec);  // zero_grad tmp_nvc tmp_nuxM

    size += 1 * memsize_mat(nuM + nxM + 1, nuM + nxM);  // zero_hess
    size += 1 * memsize_vec(nuM + nxM);  // zero_grad
    size += 1 * memsize_vec(nvc);  // tmp_nvc
    size += 1 * memsize_vec(nuM + nxM);  // tmp_nuxM
    int nu_tmp = 0;
    for (ii = 0; ii < N; ii++) {
        nu_tmp += nu[ii];
        size += memsize_mat(nu_tmp + nx[0] + 1, nx[ii + 1]);  // GammaQ
        //		size += memsize_mat(nu_tmp+nx[0], nx[ii+1]); // GammaQ TODO without the +1 ???
    }
    for (ii = 0; ii <= N; ii++) {
        size += memsize_mat(nu[ii] + nx[ii], ng[ii] + nq[ii]);
    }

    size = (size + 63) / 64 * 64;  // make multiple of typical cache line size
    size += 1 * 64;  // align once to typical cache line size

    return size;
}


void d_cond_qcqp_ws_create(struct d_ocp_qcqp_dim* ocp_dim, struct d_cond_qcqp_arg* cond_arg, struct d_cond_qcqp_ws* cond_ws, void* mem) {

    // loop index
    int ii;

    // zero memory (to avoid corrupted memory like e.g. NaN)
    hpipm_size_t memsize = d_cond_qcqp_ws_memsize(ocp_dim, cond_arg);
    hpipm_zero_memset(memsize, mem);

    int N = ocp_dim->N;
    int* nx = ocp_dim->nx;
    int* nu = ocp_dim->nu;
    int* nb = ocp_dim->nb;
    int* ng = ocp_dim->ng;
    int* nq = ocp_dim->nq;

    // compute core qp dim and max dim
    //	int nvt = 0;
    //	int net = 0;
    //	int nbt = 0;
    //	int ngt = 0;
    int nxM = 0;
    int nuM = 0;
    //	int nbM = 0;
    //	int ngM = 0;

    for (ii = 0; ii < N; ii++) {
        //		nvt += nx[ii]+nu[ii];
        //		net += nx[ii+1];
        //		nbt += nb[ii];
        //		ngt += ng[ii];
        nxM = nx[ii] > nxM ? nx[ii] : nxM;
        nuM = nu[ii] > nuM ? nu[ii] : nuM;
        //		nbM = nb[ii]>nbM ? nb[ii] : nbM;
        //		ngM = ng[ii]>ngM ? ng[ii] : ngM;
    }
    ii = N;
    //	nvt += nx[ii]+nu[ii];
    //	nbt += nb[ii];
    //	ngt += ng[ii];
    nxM = nx[ii] > nxM ? nx[ii] : nxM;
    nuM = nu[ii] > nuM ? nu[ii] : nuM;
    //	nbM = nb[ii]>nbM ? nb[ii] : nbM;
    //	ngM = ng[ii]>ngM ? ng[ii] : ngM;

    int nvc = 0;
    nvc += nx[0] + nu[0];
    for (ii = 1; ii <= N; ii++)
        nvc += nu[ii];


    // cond qp ws struct
    struct d_cond_qp_ws* ws_ptr = mem;

    cond_ws->qp_ws = ws_ptr;
    ws_ptr += 1;

    // matrix struct
    struct mat* sm_ptr = (struct mat*) ws_ptr;

    cond_ws->hess_array = sm_ptr;
    sm_ptr += N + 1;
    cond_ws->zero_hess = sm_ptr;
    sm_ptr += 1;
    cond_ws->GammaQ = sm_ptr;
    sm_ptr += N + 1;
    cond_ws->tmp_DCt = sm_ptr;
    sm_ptr += N + 1;

    // vector struct
    struct vec* sv_ptr = (struct vec*) sm_ptr;

    cond_ws->grad_array = sv_ptr;
    sv_ptr += N + 1;
    cond_ws->zero_grad = sv_ptr;
    sv_ptr += 1;
    cond_ws->tmp_nvc = sv_ptr;
    sv_ptr += 1;
    cond_ws->tmp_nuxM = sv_ptr;
    sv_ptr += 1;

    // align to typical cache line size
    hpipm_size_t s_ptr = (hpipm_size_t) sv_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;

    // void stuf
    char* c_ptr = (char*) s_ptr;
    char* c_tmp;

    //
    int nu_tmp = 0;
    for (ii = 0; ii < N; ii++) {
        nu_tmp += nu[ii];
        create_mat(nu_tmp + nx[0] + 1, nx[ii + 1], cond_ws->GammaQ + ii, c_ptr);
        c_ptr += (cond_ws->GammaQ + ii)->memsize;
    }
    //
    for (ii = 0; ii <= N; ii++) {
        create_mat(nu[ii] + nx[ii], ng[ii] + nq[ii], cond_ws->tmp_DCt + ii, c_ptr);
        c_ptr += (cond_ws->tmp_DCt + ii)->memsize;
    }
    //
    d_cond_qp_ws_create(ocp_dim->qp_dim, cond_arg->qp_arg, cond_ws->qp_ws, c_ptr);
    c_ptr += cond_ws->qp_ws->memsize;
    //
    create_mat(nuM + nxM + 1, nuM + nxM, cond_ws->zero_hess, c_ptr);
    c_ptr += cond_ws->zero_hess->memsize;
    //
    create_vec(nuM + nxM, cond_ws->zero_grad, c_ptr);
    c_ptr += cond_ws->zero_grad->memsize;
    //
    create_vec(nvc, cond_ws->tmp_nvc, c_ptr);
    c_ptr += cond_ws->tmp_nvc->memsize;
    //
    create_vec(nuM + nxM, cond_ws->tmp_nuxM, c_ptr);
    c_ptr += cond_ws->tmp_nuxM->memsize;

    cond_ws->memsize = memsize;  // d_cond_qcqp_ws_memsize(ocp_dim, cond_arg);

#if defined(RUNTIME_CHECKS)
    if (c_ptr > ((char*) mem) + cond_ws->memsize) {
        printf("\nerror: d_cond_qcqp_ws_create: outside memory bounds!\n\n");
        exit(1);
    }
#endif
}


void d_cond_qcqp_qc(struct d_ocp_qcqp* ocp_qp, struct mat* Hq2, int* Hq_nzero2, struct mat* Ct2, struct vec* d2, struct d_cond_qcqp_arg* cond_arg, struct d_cond_qcqp_ws* cond_ws) {

    // cond quadr constr
    int ii, jj, kk;

    int N = ocp_qp->dim->N;
    if (cond_arg->cond_last_stage == 0)
        N -= 1;

    // early return
    if (N < 0) { return; }


    int* nu = ocp_qp->dim->nu;
    int* nx = ocp_qp->dim->nx;
    int* nbu = ocp_qp->dim->nbu;
    int* nbx = ocp_qp->dim->nbx;
    int* ng = ocp_qp->dim->ng;
    int* nq = ocp_qp->dim->nq;

    // early return
    if (N == 0) {
        for (jj = 0; jj < nq[0]; jj++) {
            dgecp(nu[0] + nx[0], nu[0] + nx[0], ocp_qp->Hq[0] + jj, 0, 0, Hq2 + jj, 0, 0);
        }
    }

    int nvc = nu[0] + nx[0];
    int nbc = nbu[0] + nbx[0];
    int ngc = ng[0];
    int nqc = nq[0];
    for (ii = 1; ii <= N; ii++) {
        nvc += nu[ii];
        nbc += nbu[ii];
        ngc += nbx[ii] + ng[ii];
        nqc += nq[ii];
    }

    int nu_tmp, nq_tmp, nu_tot_tmp, nq_tot_tmp;

    double rho;

    dgese(nvc, nqc, 0.0, Ct2, 0, ngc);

    nu_tmp = 0;
    nq_tmp = 0;
    nu_tot_tmp = nvc - nx[0];
    nq_tot_tmp = nqc;
    for (kk = 0; kk <= N; kk++) {

        nu_tot_tmp -= nu[kk];
        nq_tot_tmp -= nq[kk];

        nq_tmp = nq_tot_tmp;

        for (jj = 0; jj < nq[kk]; jj++) {

            dgese(nvc + 1, nvc, 0.0, Hq2 + nq_tmp, 0, 0);
            Hq_nzero2[nq_tmp] = 0;
            rho = 0.0;

            if (kk == 0) {
                dtrcp_l(nu[0] + nx[0], ocp_qp->Hq[kk] + nq[jj], 0, 0, Hq2 + nq_tmp, nu_tot_tmp, nu_tot_tmp);
                if (nx[0] > 0) {
                    if (ocp_qp->Hq_nzero[kk][jj] % 2 == 1)  // Q nzero
                        Hq_nzero2[nq_tmp] |= 1;
                    if ((ocp_qp->Hq_nzero[kk][jj] >> 1) % 2 == 1)  // S nzero
                        Hq_nzero2[nq_tmp] |= 2;
                }
                if ((ocp_qp->Hq_nzero[kk][jj] >> 2) % 2 == 1)  // R nzero
                    Hq_nzero2[nq_tmp] |= 4;
            } else {
                // XXX make Q full or use SYMM
                // hessian
                if (ocp_qp->Hq_nzero[kk][jj] % 2 == 1)
                // if((ocp_qp->Hq_nzero[kk][jj]&0x1)!=0) // faster ???
                {
                    //					dtrtr_l(nx[kk], ocp_qp->Hq[kk]+jj, nu[kk], nu[kk], ocp_qp->Hq[kk]+jj, nu[kk], nu[kk]); // buggy ???
                    dtrtr_l(nu[kk] + nx[kk], ocp_qp->Hq[kk] + jj, 0, 0, ocp_qp->Hq[kk] + jj, 0, 0);
                    dgemm_nn(nu_tmp + nx[0] + 1, nx[kk], nx[kk], 1.0, cond_ws->qp_ws->Gamma + kk - 1, 0, 0, ocp_qp->Hq[kk] + jj, nu[kk], nu[kk], 0.0, cond_ws->GammaQ + kk - 1, 0, 0, cond_ws->GammaQ + kk - 1, 0, 0);  // TODO check for diagonal Q to reduce this cost
                    drowex(nx[kk], 1.0, cond_ws->GammaQ + kk - 1, nu_tmp + nx[0], 0, cond_ws->tmp_nuxM, 0);
                    //					dsymv_l(nx[kk], 1.0, ocp_qp->Hq[kk]+jj, nu[kk], nu[kk], cond_ws->qp_ws->Gammab+kk-1, 0, 0.0, cond_ws->tmp_nuxM, 0, cond_ws->tmp_nuxM, 0);
                    rho = 0.5 * ddot(nx[kk], cond_ws->tmp_nuxM, 0, cond_ws->qp_ws->Gammab + kk - 1, 0);
                    dsyrk_ln_mn(nu_tmp + nx[0] + 1, nu_tmp + nx[0], nx[kk], 1.0, cond_ws->qp_ws->Gamma + kk - 1, 0, 0, cond_ws->GammaQ + kk - 1, 0, 0, 0.0, Hq2 + nq_tmp, nu_tot_tmp + nu[kk], nu_tot_tmp + nu[kk], Hq2 + nq_tmp, nu_tot_tmp + nu[kk], nu_tot_tmp + nu[kk]);
                    if (nx[0] > 0)
                        Hq_nzero2[nq_tmp] |= 7;
                    else
                        Hq_nzero2[nq_tmp] |= 4;
                }
                if ((ocp_qp->Hq_nzero[kk][jj] >> 2) % 2 == 1)
                // if((ocp_qp->Hq_nzero[kk][jj]&0x4)!=0) // faster ???
                {
                    dgead(nu[kk], nu[kk], 1.0, ocp_qp->Hq[kk] + jj, 0, 0, Hq2 + nq_tmp, nu_tot_tmp, nu_tot_tmp);
                    Hq_nzero2[nq_tmp] |= 4;
                }
                if ((ocp_qp->Hq_nzero[kk][jj] >> 1) % 2 == 1)
                // if((ocp_qp->Hq_nzero[kk][jj]&0x2)!=0) // faster ???
                {
                    dgemm_nn(nu_tmp + nx[0] + 1, nu[kk], nx[kk], 1.0, cond_ws->qp_ws->Gamma + kk - 1, 0, 0, ocp_qp->Hq[kk] + jj, nu[kk], 0, 1.0, Hq2 + nq_tmp, nu_tot_tmp + nu[kk], nu_tot_tmp, Hq2 + nq_tmp, nu_tot_tmp + nu[kk], nu_tot_tmp);
                    if (nx[0] > 0)
                        Hq_nzero2[nq_tmp] |= 6;
                    else
                        Hq_nzero2[nq_tmp] |= 4;
                }
                // gradient
                drowex(nu_tmp + nx[0], 1.0, Hq2 + nq_tmp, nvc, nu_tot_tmp + nu[kk], cond_ws->tmp_nvc, 0);
                dcolad(nu_tmp + nx[0], 1.0, cond_ws->tmp_nvc, 0, Ct2, nu_tot_tmp + nu[kk], ngc + nq_tmp);
            }

            //			printf("\nrho %d %f\n", kk, rho);
            VECEL(d2, nbc + ngc + nq_tmp) -= rho;
            VECEL(d2, 2 * nbc + 2 * ngc + nqc + nq_tmp) += rho;

            nq_tmp++;
        }

        nu_tmp += nu[kk];
    }
}


void d_cond_qcqp_qc_lhs(struct d_ocp_qcqp* ocp_qp, struct mat* Hq2, int* Hq_nzero2, struct mat* Ct2, struct d_cond_qcqp_arg* cond_arg, struct d_cond_qcqp_ws* cond_ws) {

    // cond quadr constr
    int ii, jj, kk;

    int N = ocp_qp->dim->N;
    if (cond_arg->cond_last_stage == 0)
        N -= 1;

    // early return
    if (N < 0) { return; }


    int* nu = ocp_qp->dim->nu;
    int* nx = ocp_qp->dim->nx;
    int* nbu = ocp_qp->dim->nbu;
    int* nbx = ocp_qp->dim->nbx;
    int* ng = ocp_qp->dim->ng;
    int* nq = ocp_qp->dim->nq;

    // early return
    if (N == 0) {
        for (jj = 0; jj < nq[0]; jj++) {
            dgecp(nu[0] + nx[0], nu[0] + nx[0], ocp_qp->Hq[0] + jj, 0, 0, Hq2 + jj, 0, 0);
        }
    }

    int nvc = nu[0] + nx[0];
    int nbc = nbu[0] + nbx[0];
    int ngc = ng[0];
    int nqc = nq[0];
    for (ii = 1; ii <= N; ii++) {
        nvc += nu[ii];
        nbc += nbu[ii];
        ngc += nbx[ii] + ng[ii];
        nqc += nq[ii];
    }

    int nu_tmp, nq_tmp, nu_tot_tmp, nq_tot_tmp;

    double rho;

    dgese(nvc, nqc, 0.0, Ct2, 0, ngc);

    nu_tmp = 0;
    nq_tmp = 0;
    nu_tot_tmp = nvc - nx[0];
    nq_tot_tmp = nqc;
    for (kk = 0; kk <= N; kk++) {

        nu_tot_tmp -= nu[kk];
        nq_tot_tmp -= nq[kk];

        nq_tmp = nq_tot_tmp;

        for (jj = 0; jj < nq[kk]; jj++) {

            dgese(nvc + 1, nvc, 0.0, Hq2 + nq_tmp, 0, 0);
            Hq_nzero2[nq_tmp] = 0;
            rho = 0.0;

            if (kk == 0) {
                dtrcp_l(nu[0] + nx[0], ocp_qp->Hq[kk] + nq[jj], 0, 0, Hq2 + nq_tmp, nu_tot_tmp, nu_tot_tmp);
                if (nx[0] > 0) {
                    if (ocp_qp->Hq_nzero[kk][jj] % 2 == 1)  // Q nzero
                        Hq_nzero2[nq_tmp] |= 1;
                    if ((ocp_qp->Hq_nzero[kk][jj] >> 1) % 2 == 1)  // S nzero
                        Hq_nzero2[nq_tmp] |= 2;
                }
                if ((ocp_qp->Hq_nzero[kk][jj] >> 2) % 2 == 1)  // R nzero
                    Hq_nzero2[nq_tmp] |= 4;
            } else {
                // XXX make Q full or use SYMM
                // hessian
                if (ocp_qp->Hq_nzero[kk][jj] % 2 == 1) {
                    //					dtrtr_l(nx[kk], ocp_qp->Hq[kk]+jj, nu[kk], nu[kk], ocp_qp->Hq[kk]+jj, nu[kk], nu[kk]); // buggy ???
                    dtrtr_l(nu[kk] + nx[kk], ocp_qp->Hq[kk] + jj, 0, 0, ocp_qp->Hq[kk] + jj, 0, 0);
                    dgemm_nn(nu_tmp + nx[0] + 1, nx[kk], nx[kk], 1.0, cond_ws->qp_ws->Gamma + kk - 1, 0, 0, ocp_qp->Hq[kk] + jj, nu[kk], nu[kk], 0.0, cond_ws->GammaQ + kk - 1, 0, 0, cond_ws->GammaQ + kk - 1, 0, 0);
                    //					drowex(nx[kk], 1.0, cond_ws->GammaQ+kk-1, nu_tmp+nx[0], 0, cond_ws->tmp_nuxM, 0);
                    //					dsymv_l(nx[kk], 1.0, ocp_qp->Hq[kk]+jj, nu[kk], nu[kk], cond_ws->qp_ws->Gammab+kk-1, 0, 0.0, cond_ws->tmp_nuxM, 0, cond_ws->tmp_nuxM, 0);
                    //					rho = 0.5*ddot(nx[kk], cond_ws->tmp_nuxM, 0, cond_ws->qp_ws->Gammab+kk-1, 0);
                    dsyrk_ln_mn(nu_tmp + nx[0] + 1, nu_tmp + nx[0], nx[kk], 1.0, cond_ws->qp_ws->Gamma + kk - 1, 0, 0, cond_ws->GammaQ + kk - 1, 0, 0, 0.0, Hq2 + nq_tmp, nu_tot_tmp + nu[kk], nu_tot_tmp + nu[kk], Hq2 + nq_tmp, nu_tot_tmp + nu[kk], nu_tot_tmp + nu[kk]);
                    if (nx[0] > 0)
                        Hq_nzero2[nq_tmp] |= 7;
                    else
                        Hq_nzero2[nq_tmp] |= 4;
                }
                if ((ocp_qp->Hq_nzero[kk][jj] >> 2) % 2 == 1) {
                    dgead(nu[kk], nu[kk], 1.0, ocp_qp->Hq[kk] + jj, 0, 0, Hq2 + nq_tmp, nu_tot_tmp, nu_tot_tmp);
                    Hq_nzero2[nq_tmp] |= 4;
                }
                if ((ocp_qp->Hq_nzero[kk][jj] >> 1) % 2 == 1) {
                    dgemm_nn(nu_tmp + nx[0] + 1, nu[kk], nx[kk], 1.0, cond_ws->qp_ws->Gamma + kk - 1, 0, 0, ocp_qp->Hq[kk] + jj, nu[kk], 0, 1.0, Hq2 + nq_tmp, nu_tot_tmp + nu[kk], nu_tot_tmp, Hq2 + nq_tmp, nu_tot_tmp + nu[kk], nu_tot_tmp);
                    if (nx[0] > 0)
                        Hq_nzero2[nq_tmp] |= 6;
                    else
                        Hq_nzero2[nq_tmp] |= 4;
                }
                // gradient
                drowex(nu_tmp + nx[0], 1.0, Hq2 + nq_tmp, nvc, nu_tot_tmp + nu[kk], cond_ws->tmp_nvc, 0);
                dcolad(nu_tmp + nx[0], 1.0, cond_ws->tmp_nvc, 0, Ct2, nu_tot_tmp + nu[kk], ngc + nq_tmp);
            }

            nq_tmp++;
        }

        nu_tmp += nu[kk];
    }
}


void d_cond_qcqp_qc_rhs(struct d_ocp_qcqp* ocp_qp, struct vec* d2, struct d_cond_qcqp_arg* cond_arg, struct d_cond_qcqp_ws* cond_ws) {

    // cond quadr constr
    int ii, jj, kk;

    int N = ocp_qp->dim->N;
    if (cond_arg->cond_last_stage == 0)
        N -= 1;

    // early return
    if (N < 0) { return; }


    int* nu = ocp_qp->dim->nu;
    int* nx = ocp_qp->dim->nx;
    int* nbu = ocp_qp->dim->nbu;
    int* nbx = ocp_qp->dim->nbx;
    int* ng = ocp_qp->dim->ng;
    int* nq = ocp_qp->dim->nq;

    int nvc = nu[0] + nx[0];
    int nbc = nbu[0] + nbx[0];
    int ngc = ng[0];
    int nqc = nq[0];
    for (ii = 1; ii <= N; ii++) {
        nvc += nu[ii];
        nbc += nbu[ii];
        ngc += nbx[ii] + ng[ii];
        nqc += nq[ii];
    }

    int nu_tmp, nq_tmp, nu_tot_tmp, nq_tot_tmp;

    double rho;

    nu_tmp = 0;
    nq_tmp = 0;
    nu_tot_tmp = nvc - nx[0];
    nq_tot_tmp = nqc;
    for (kk = 0; kk <= N; kk++) {

        nu_tot_tmp -= nu[kk];
        nq_tot_tmp -= nq[kk];

        nq_tmp = nq_tot_tmp;

        for (jj = 0; jj < nq[kk]; jj++) {

            rho = 0.0;

            if (kk == 0) {
                // nothing to do
            } else {
                if (ocp_qp->Hq_nzero[kk][jj] % 2 == 1) {
                    //					dsymv_l(nx[kk], 1.0, ocp_qp->Hq[kk]+jj, nu[kk], nu[kk], cond_ws->qp_ws->Gammab+kk-1, 0, 0.0, cond_ws->tmp_nuxM, 0, cond_ws->tmp_nuxM, 0); // XXX buggy !!!
                    dtrtr_l(nu[kk] + nx[kk], ocp_qp->Hq[kk] + jj, 0, 0, ocp_qp->Hq[kk] + jj, 0, 0);
                    dgemv_n(nx[kk], nx[kk], 1.0, ocp_qp->Hq[kk] + jj, nu[kk], nu[kk], cond_ws->qp_ws->Gammab + kk - 1, 0, 0.0, cond_ws->tmp_nuxM, 0, cond_ws->tmp_nuxM, 0);
                    rho = 0.5 * ddot(nx[kk], cond_ws->tmp_nuxM, 0, cond_ws->qp_ws->Gammab + kk - 1, 0);
                    //					dgemv_n(nu_tmp+nx[0], nx[kk], 1.0, cond_ws->qp_ws->Gamma+kk-1, 0, 0, cond_ws->tmp_nuxM, 0, 0.0, cond_ws->tmp_nvc, 0, cond_ws->tmp_nvc, 0);
                }

                // XXX this does not affect rho (i.e. the RHS)
                //				if((ocp_qp->Hq_nzero[kk][jj]>>1)%2==1)
                //					{
                //					// TODO S
                //					}

                //				dcolad(nu_tmp+nx[0], 1.0, cond_ws->tmp_nvc, 0, Ct2, nu_tot_tmp+nu[kk], ngc+nq_tmp);
            }

            VECEL(d2, nbc + ngc + nq_tmp) -= rho;
            VECEL(d2, 2 * nbc + 2 * ngc + nqc + nq_tmp) += rho;

            nq_tmp++;
        }

        nu_tmp += nu[kk];
    }
}


void d_cond_qcqp_cond(struct d_ocp_qcqp* ocp_qp, struct d_dense_qcqp* dense_qp, struct d_cond_qcqp_arg* cond_arg, struct d_cond_qcqp_ws* cond_ws) {

    // create tmp QP
    struct d_ocp_qp tmp_ocp_qp;

    // alias
    tmp_ocp_qp.dim = ocp_qp->dim->qp_dim;
    tmp_ocp_qp.idxb = ocp_qp->idxb;
    tmp_ocp_qp.BAbt = ocp_qp->BAbt;
    tmp_ocp_qp.b = ocp_qp->b;
    tmp_ocp_qp.RSQrq = ocp_qp->RSQrq;
    tmp_ocp_qp.rqz = ocp_qp->rqz;
    tmp_ocp_qp.DCt = ocp_qp->DCt;
    tmp_ocp_qp.d = ocp_qp->d;
    tmp_ocp_qp.d_mask = ocp_qp->d_mask;
    tmp_ocp_qp.Z = ocp_qp->Z;
    tmp_ocp_qp.idxs_rev = ocp_qp->idxs_rev;

    d_cond_BAbt(&tmp_ocp_qp, NULL, NULL, cond_arg->qp_arg, cond_ws->qp_ws);

    d_cond_RSQrq(&tmp_ocp_qp, dense_qp->Hv, dense_qp->gz, cond_arg->qp_arg, cond_ws->qp_ws);

    d_cond_DCtd(&tmp_ocp_qp, dense_qp->idxb, dense_qp->Ct, dense_qp->d, dense_qp->d_mask, dense_qp->idxs_rev, dense_qp->Z, dense_qp->gz, cond_arg->qp_arg, cond_ws->qp_ws);

    d_cond_qcqp_qc(ocp_qp, dense_qp->Hq, dense_qp->Hq_nzero, dense_qp->Ct, dense_qp->d, cond_arg, cond_ws);
}


void d_cond_qcqp_cond_lhs(struct d_ocp_qcqp* ocp_qp, struct d_dense_qcqp* dense_qp, struct d_cond_qcqp_arg* cond_arg, struct d_cond_qcqp_ws* cond_ws) {

    // create tmp QP
    struct d_ocp_qp tmp_ocp_qp;

    // alias
    tmp_ocp_qp.dim = ocp_qp->dim->qp_dim;
    tmp_ocp_qp.idxb = ocp_qp->idxb;
    tmp_ocp_qp.BAbt = ocp_qp->BAbt;
    tmp_ocp_qp.b = ocp_qp->b;
    tmp_ocp_qp.RSQrq = ocp_qp->RSQrq;
    tmp_ocp_qp.rqz = ocp_qp->rqz;
    tmp_ocp_qp.DCt = ocp_qp->DCt;
    tmp_ocp_qp.d = ocp_qp->d;
    tmp_ocp_qp.d_mask = ocp_qp->d_mask;
    tmp_ocp_qp.Z = ocp_qp->Z;
    tmp_ocp_qp.idxs_rev = ocp_qp->idxs_rev;

    d_cond_BAt(&tmp_ocp_qp, NULL, cond_arg->qp_arg, cond_ws->qp_ws);

    d_cond_RSQ(&tmp_ocp_qp, dense_qp->Hv, cond_arg->qp_arg, cond_ws->qp_ws);

    d_cond_DCt(&tmp_ocp_qp, dense_qp->idxb, dense_qp->Ct, dense_qp->idxs_rev, dense_qp->Z, cond_arg->qp_arg, cond_ws->qp_ws);

    d_cond_qcqp_qc_lhs(ocp_qp, dense_qp->Hq, dense_qp->Hq_nzero, dense_qp->Ct, cond_arg, cond_ws);
}


void d_cond_qcqp_cond_rhs(struct d_ocp_qcqp* ocp_qp, struct d_dense_qcqp* dense_qp, struct d_cond_qcqp_arg* cond_arg, struct d_cond_qcqp_ws* cond_ws) {

    // create tmp QP
    struct d_ocp_qp tmp_ocp_qp;

    // alias
    tmp_ocp_qp.dim = ocp_qp->dim->qp_dim;
    tmp_ocp_qp.idxb = ocp_qp->idxb;
    tmp_ocp_qp.BAbt = ocp_qp->BAbt;
    tmp_ocp_qp.b = ocp_qp->b;
    tmp_ocp_qp.RSQrq = ocp_qp->RSQrq;
    tmp_ocp_qp.rqz = ocp_qp->rqz;
    tmp_ocp_qp.DCt = ocp_qp->DCt;
    tmp_ocp_qp.d = ocp_qp->d;
    tmp_ocp_qp.d_mask = ocp_qp->d_mask;
    tmp_ocp_qp.Z = ocp_qp->Z;
    tmp_ocp_qp.idxs_rev = ocp_qp->idxs_rev;

    d_cond_b(&tmp_ocp_qp, NULL, cond_arg->qp_arg, cond_ws->qp_ws);

    d_cond_rq(&tmp_ocp_qp, dense_qp->gz, cond_arg->qp_arg, cond_ws->qp_ws);

    d_cond_d(&tmp_ocp_qp, dense_qp->d, dense_qp->d_mask, dense_qp->gz, cond_arg->qp_arg, cond_ws->qp_ws);

    d_cond_qcqp_qc_rhs(ocp_qp, dense_qp->d, cond_arg, cond_ws);
}


void d_cond_qcqp_expand_sol(struct d_ocp_qcqp* ocp_qp, struct d_dense_qcqp_sol* dense_qp_sol, struct d_ocp_qcqp_sol* ocp_qp_sol, struct d_cond_qcqp_arg* cond_arg, struct d_cond_qcqp_ws* cond_ws) {

    int ii, jj;

    // extract dim
    int N = ocp_qp->dim->N;
    int* nu = ocp_qp->dim->nu;
    int* nx = ocp_qp->dim->nx;
    int* ng = ocp_qp->dim->ng;
    int* nq = ocp_qp->dim->nq;

    // create tmp QP
    struct d_ocp_qp tmp_ocp_qp;

    // alias
    tmp_ocp_qp.dim = ocp_qp->dim->qp_dim;
    tmp_ocp_qp.idxb = ocp_qp->idxb;
    tmp_ocp_qp.BAbt = ocp_qp->BAbt;
    tmp_ocp_qp.b = ocp_qp->b;
    tmp_ocp_qp.RSQrq = ocp_qp->RSQrq;
    tmp_ocp_qp.rqz = ocp_qp->rqz;
    tmp_ocp_qp.DCt = ocp_qp->DCt;
    tmp_ocp_qp.d = ocp_qp->d;
    tmp_ocp_qp.d_mask = ocp_qp->d_mask;
    tmp_ocp_qp.Z = ocp_qp->Z;
    tmp_ocp_qp.idxs_rev = ocp_qp->idxs_rev;
    // TODO d_mask ????????


    // create tmp QP
    struct d_ocp_qp_sol tmp_ocp_qp_sol;

    // alias
    tmp_ocp_qp_sol.dim = ocp_qp_sol->dim->qp_dim;
    tmp_ocp_qp_sol.ux = ocp_qp_sol->ux;
    tmp_ocp_qp_sol.pi = ocp_qp_sol->pi;
    tmp_ocp_qp_sol.lam = ocp_qp_sol->lam;
    tmp_ocp_qp_sol.t = ocp_qp_sol->t;


    // create tmp QP
    struct d_dense_qp_sol tmp_dense_qp_sol;

    // alias
    tmp_dense_qp_sol.dim = dense_qp_sol->dim->qp_dim;
    tmp_dense_qp_sol.v = dense_qp_sol->v;
    tmp_dense_qp_sol.pi = dense_qp_sol->pi;
    tmp_dense_qp_sol.lam = dense_qp_sol->lam;
    tmp_dense_qp_sol.t = dense_qp_sol->t;

    int bkp_comp_prim_sol = cond_arg->qp_arg->comp_prim_sol;
    int bkp_comp_dual_sol_eq = cond_arg->qp_arg->comp_dual_sol_eq;
    int bkp_comp_dual_sol_ineq = cond_arg->qp_arg->comp_dual_sol_ineq;

    cond_arg->qp_arg->comp_prim_sol = 1 & bkp_comp_prim_sol;
    cond_arg->qp_arg->comp_dual_sol_eq = 0 & bkp_comp_dual_sol_eq;
    cond_arg->qp_arg->comp_dual_sol_ineq = 1 & bkp_comp_dual_sol_ineq;

    //	cond_arg->comp_dual_sol = 0; // compute dual solution
    //	cond_arg->qp_arg->comp_dual_sol = 0; // compute dual solution

    d_expand_sol(&tmp_ocp_qp, &tmp_dense_qp_sol, &tmp_ocp_qp_sol, cond_arg->qp_arg, cond_ws->qp_ws);

    // linearize quadr constr
    for (ii = 0; ii <= N; ii++) {
        dgecp(nu[ii] + nx[ii], ng[ii] + nq[ii], ocp_qp->DCt + ii, 0, 0, cond_ws->tmp_DCt + ii, 0, 0);
        for (jj = 0; jj < nq[ii]; jj++) {
            dsymv_l(nu[ii] + nx[ii], 1.0, ocp_qp->Hq[ii] + jj, 0, 0, ocp_qp_sol->ux + ii, 0, 0.0, cond_ws->tmp_nuxM, 0, cond_ws->tmp_nuxM, 0);
            dcolad(nu[ii] + nx[ii], 1.0, cond_ws->tmp_nuxM, 0, cond_ws->tmp_DCt + ii, 0, ng[ii] + jj);
        }
    }

    tmp_ocp_qp.DCt = cond_ws->tmp_DCt;

    cond_arg->qp_arg->comp_prim_sol = 0 & bkp_comp_prim_sol;
    ;
    cond_arg->qp_arg->comp_dual_sol_eq = 1 & bkp_comp_dual_sol_eq;
    cond_arg->qp_arg->comp_dual_sol_ineq = 0 & bkp_comp_dual_sol_ineq;

    d_expand_sol(&tmp_ocp_qp, &tmp_dense_qp_sol, &tmp_ocp_qp_sol, cond_arg->qp_arg, cond_ws->qp_ws);

    cond_arg->qp_arg->comp_prim_sol = bkp_comp_prim_sol;
    cond_arg->qp_arg->comp_dual_sol_eq = bkp_comp_dual_sol_eq;
    cond_arg->qp_arg->comp_dual_sol_ineq = bkp_comp_dual_sol_ineq;
}
