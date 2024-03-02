#include <math.h>
#include <stdlib.h>

#include "tinyhpipm/blas.h"
#include "tinyhpipm/common.h"
#include "tinyhpipm/dense/d_dense_qp.h"
#include "tinyhpipm/dense/d_dense_qp_dim.h"
#include "tinyhpipm/dense/d_dense_qp_ipm.h"
#include "tinyhpipm/dense/d_dense_qp_res.h"
#include "tinyhpipm/dense/d_dense_qp_sol.h"
#include "tinyhpipm/ipm_core/d_core_qp_ipm.h"
#include "tinyhpipm/ipm_core/d_core_qp_ipm_aux.h"

// range-space (Schur complement) method
void d_fact_solve_kkt_unconstr_dense_qp(struct d_dense_qp* qp, struct d_dense_qp_sol* qp_sol, struct d_dense_qp_ipm_arg* arg, struct d_dense_qp_ipm_ws* ws) {

    int ii, jj;

    int nv = qp->dim->nv;
    int ne = qp->dim->ne;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;

    struct mat* Hg = qp->Hv;
    struct mat* A = qp->A;
    struct vec* gz = qp->gz;
    struct vec* b = qp->b;

    struct vec* v = qp_sol->v;
    struct vec* pi = qp_sol->pi;

    struct mat* Lv = ws->Lv;
    struct mat* Le = ws->Le;
    struct mat* Ctx = ws->Ctx;
    struct mat* AL = ws->AL;
    struct vec* lv = ws->lv;

    // null space
    struct mat* A_LQ = ws->A_LQ;
    struct mat* A_Q = ws->A_Q;
    struct mat* Zt = ws->Zt;
    struct mat* ZtH = ws->ZtH;
    struct mat* ZtHZ = ws->ZtHZ;
    struct vec* xy = ws->xy;
    struct vec* Yxy = ws->Yxy;
    struct vec* xz = ws->xz;
    struct vec* tmp_nv = ws->tmp_nv;
    void* lq_work_null = ws->lq_work_null;
    void* orglq_work_null = ws->orglq_work_null;

    if (ne > 0) {
        if (arg->kkt_fact_alg == 0)  // null space method
        {
            // TODO check to cache LQ across IPM iterations !!!!!
            dgelqf(ne, nv, A, 0, 0, A_LQ, 0, 0, lq_work_null);
            //			printf("\nA_LQ\n");
            //			blasfeo_print_dmat(ne, nv, A_LQ, 0, 0);

            // TODO cache dA containing tau into another vector !!!!!

            // TODO change dorglq API to pass tau explicitly as a vector !!!!!
            // TODO allocate its dedicated workspace !!!!!
            dorglq(nv, nv, ne, A_LQ, 0, 0, A_Q, 0, 0, orglq_work_null);
            //			printf("\nA_Q\n");
            //			blasfeo_print_dmat(nv, nv, A_Q, 0, 0);

            dgecp(nv - ne, nv, A_Q, ne, 0, Zt, 0, 0);
            //			printf("\nZt\n");
            //			blasfeo_print_dmat(nv-ne, nv, Zt, 0, 0);

#if 0
//printf("\nA\n");
//blasfeo_print_dmat(ne, nv, A, 0, 0);
//printf("\nA_LQ\n");
//blasfeo_print_dmat(ne, nv, A_LQ, 0, 0);
//printf("\nA_Q\n");
//blasfeo_print_dmat(nv, nv, A_Q, 0, 0);

for(ii=0; ii<ne; ii++)
	for(jj=ii+1; jj<nv; jj++)
		MATEL(A_LQ, ii, jj) = 0;
blasfeo_dgemm_nn(ne, nv, nv, 1.0, A_LQ, 0, 0, A_Q, 0, 0, -1.0, A, 0, 0, AL, 0, 0);
double max_err = 0.0;
double tmp;
for(ii=0; ii<ne; ii++)
	for(jj=ii+1; jj<nv; jj++)
		{
		tmp = MATEL(AL, ii, jj);
		max_err = fabs(tmp)>max_err ? fabs(tmp) : max_err;
		}
printf("\nA_LQ * A_Q - A max err %e\n", max_err);
//blasfeo_print_exp_dmat(ne, nv, AL, 0, 0);
//exit(1);
#endif

            dgemm_nt(nv - ne, nv, nv, 1.0, Zt, 0, 0, Hg, 0, 0, 0.0, ZtH, 0, 0, ZtH, 0, 0);
            //			printf("\nZtH\n");
            //			blasfeo_print_dmat(nv-ne, nv, ZtH, 0, 0);

            dsyrk_ln(nv - ne, nv, 1.0, ZtH, 0, 0, Zt, 0, 0, 0.0, ZtHZ, 0, 0, ZtHZ, 0, 0);
            //			printf("\nZtHZ\n");
            //			blasfeo_print_dmat(nv-ne, nv-ne, ZtHZ, 0, 0);

            dpotrf_l(nv - ne, ZtHZ, 0, 0, ZtHZ, 0, 0);
            //			printf("\nZtHZ\n");
            //			blasfeo_print_dmat(nv-ne, nv-ne, ZtHZ, 0, 0);

            dtrsv_lnn(ne, A_LQ, 0, 0, b, 0, xy, 0);
            //			printf("\nxy\n");
            //			blasfeo_print_dvec(ne, xy, 0);

            dgemv_t(ne, nv, 1.0, A_Q, 0, 0, xy, 0, 0.0, Yxy, 0, Yxy, 0);
            //			printf("\nYxy\n");
            //			blasfeo_print_dvec(nv, Yxy, 0);

            dgemv_n(nv - ne, nv, -1.0, ZtH, 0, 0, Yxy, 0, 0.0, xz, 0, xz, 0);
            //			printf("\nxz\n");
            //			blasfeo_print_dvec(nv-ne, xz, 0);

            dgemv_n(nv - ne, nv, -1.0, Zt, 0, 0, gz, 0, 1.0, xz, 0, xz, 0);
            //			printf("\nxz\n");
            //			blasfeo_print_dvec(nv-ne, xz, 0);

            dtrsv_lnn(nv - ne, ZtHZ, 0, 0, xz, 0, xz, 0);
            //			printf("\nxz\n");
            //			blasfeo_print_dvec(nv-ne, xz, 0);

            dtrsv_ltn(nv - ne, ZtHZ, 0, 0, xz, 0, xz, 0);
            //			printf("\nxz\n");
            //			blasfeo_print_dvec(nv-ne, xz, 0);

            dgemv_t(nv - ne, nv, 1.0, Zt, 0, 0, xz, 0, 1.0, Yxy, 0, v, 0);
            //			printf("\nv\n");
            //			blasfeo_print_dvec(nv, v, 0);

            dsymv_l(nv, 1.0, Hg, 0, 0, v, 0, 1.0, gz, 0, tmp_nv, 0);
            //			printf("\ntmp_nv\n");
            //			blasfeo_print_dvec(nv, tmp_nv, 0);

            dgemv_n(ne, nv, 1.0, A_Q, 0, 0, tmp_nv, 0, 0.0, pi, 0, pi, 0);
            //			printf("\npi\n");
            //			blasfeo_print_dvec(ne, pi, 0);

            dtrsv_ltn(ne, A_LQ, 0, 0, pi, 0, pi, 0);
            //			printf("\npi\n");
            //			blasfeo_print_dvec(ne, pi, 0);

            // TODO
            //			printf("\ndone!\n");
            //			exit(1);
        } else  // range space method
        {
            dpotrf_l(nv, Hg, 0, 0, Lv, 0, 0);

            //			dgecp(ne, nv, A, 0, 0, AL, 0, 0);
            dtrsm_rltn(ne, nv, 1.0, Lv, 0, 0, A, 0, 0, AL, 0, 0);

            dgese(ne, ne, 0.0, Le, 0, 0);
            dsyrk_dpotrf_ln(ne, nv, AL, 0, 0, AL, 0, 0, Le, 0, 0, Le, 0, 0);

            dtrsv_lnn(nv, Lv, 0, 0, gz, 0, lv, 0);

            dgemv_n(ne, nv, 1.0, AL, 0, 0, lv, 0, 1.0, b, 0, pi, 0);

            dtrsv_lnn(ne, Le, 0, 0, pi, 0, pi, 0);
            dtrsv_ltn(ne, Le, 0, 0, pi, 0, pi, 0);

            dgemv_t(ne, nv, 1.0, A, 0, 0, pi, 0, -1.0, gz, 0, v, 0);

            dtrsv_lnn(nv, Lv, 0, 0, v, 0, v, 0);
            dtrsv_ltn(nv, Lv, 0, 0, v, 0, v, 0);
        }
    } else {
#if 0
		dpotrf_l(nv, Hg, 0, 0, Lv, 0, 0);

		dveccp(nv, gz, 0, v, 0);
		dvecsc(nv, -1.0, v, 0);

		dtrsv_lnn(nv, Lv, 0, 0, v, 0, v, 0);
		dtrsv_ltn(nv, Lv, 0, 0, v, 0, v, 0);
#else
        drowin(nv, 1.0, gz, 0, Hg, nv, 0);
        dpotrf_l_mn(nv + 1, nv, Hg, 0, 0, Lv, 0, 0);

        drowex(nv, -1.0, Lv, nv, 0, v, 0);
        dtrsv_ltn(nv, Lv, 0, 0, v, 0, v, 0);
#endif
    }
}


static void COND_SLACKS_FACT_solve(struct d_dense_qp* qp, struct d_dense_qp_sol* qp_sol, struct d_dense_qp_ipm_arg* arg, struct d_dense_qp_ipm_ws* ws) {

    int ii, idx;

    int nv = qp->dim->nv;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;
    int ns = qp->dim->ns;

    struct vec* Z = qp->Z;
    int* idxs_rev = qp->idxs_rev;

    //	struct vec *dv = ws->sol_step->v;
    struct vec* dv = qp_sol->v;

    //	struct vec *res_g = ws->res->res_g;
    struct vec* res_g = qp->gz;

    struct vec* Gamma = ws->Gamma;
    struct vec* gamma = ws->gamma;
    struct vec* Zs_inv = ws->Zs_inv;
    struct vec* tmp_nbg = ws->tmp_nbg;

    double* ptr_Gamma = Gamma->pa;
    double* ptr_gamma = gamma->pa;
    double* ptr_Z = Z->pa;
    double* ptr_Zs_inv = Zs_inv->pa;
    double* ptr_dv = dv->pa;
    double* ptr_res_g = res_g->pa;
    double* ptr_tmp0 = (tmp_nbg + 0)->pa;
    double* ptr_tmp1 = (tmp_nbg + 1)->pa;
    double* ptr_tmp2 = (tmp_nbg + 2)->pa;
    double* ptr_tmp3 = (tmp_nbg + 3)->pa;

    double tmp0, tmp1;

    dveccp(nb + ng, Gamma, 0, tmp_nbg + 0, 0);
    dveccp(nb + ng, Gamma, nb + ng, tmp_nbg + 1, 0);
    dveccp(nb + ng, gamma, 0, tmp_nbg + 2, 0);
    dveccp(nb + ng, gamma, nb + ng, tmp_nbg + 3, 0);

    // idxs_rev
    for (ii = 0; ii < nb + ng; ii++) {
        idx = idxs_rev[ii];
        if (idx != -1) {
            // ii  <= constr index
            // idx <= slack index
            ptr_Zs_inv[0 + idx] = ptr_Z[0 + idx] + arg->reg_prim + ptr_Gamma[0 + ii] + ptr_Gamma[2 * nb + 2 * ng + idx];
            ptr_Zs_inv[ns + idx] = ptr_Z[ns + idx] + arg->reg_prim + ptr_Gamma[nb + ng + ii] + ptr_Gamma[2 * nb + 2 * ng + ns + idx];
            ptr_dv[nv + idx] = ptr_res_g[nv + idx] + ptr_gamma[0 + ii] + ptr_gamma[2 * nb + 2 * ng + idx];
            ptr_dv[nv + ns + idx] = ptr_res_g[nv + ns + idx] + ptr_gamma[nb + ng + ii] + ptr_gamma[2 * nb + 2 * ng + ns + idx];
            ptr_Zs_inv[0 + idx] = 1.0 / ptr_Zs_inv[0 + idx];
            ptr_Zs_inv[ns + idx] = 1.0 / ptr_Zs_inv[ns + idx];
            tmp0 = ptr_dv[nv + idx] * ptr_Zs_inv[0 + idx];
            tmp1 = ptr_dv[nv + ns + idx] * ptr_Zs_inv[ns + idx];
            ptr_tmp0[ii] = ptr_tmp0[ii] - ptr_tmp0[ii] * ptr_Zs_inv[0 + idx] * ptr_tmp0[ii];
            ptr_tmp1[ii] = ptr_tmp1[ii] - ptr_tmp1[ii] * ptr_Zs_inv[ns + idx] * ptr_tmp1[ii];
            ptr_tmp2[ii] = ptr_tmp2[ii] - ptr_Gamma[0 + ii] * tmp0;
            ptr_tmp3[ii] = ptr_tmp3[ii] - ptr_Gamma[nb + ng + ii] * tmp1;
        }
    }

    daxpy(nb + ng, 1.0, tmp_nbg + 1, 0, tmp_nbg + 0, 0, tmp_nbg + 0, 0);
    daxpy(nb + ng, -1.0, tmp_nbg + 3, 0, tmp_nbg + 2, 0, tmp_nbg + 1, 0);
}


static void COND_SLACKS_solve(struct d_dense_qp* qp, struct d_dense_qp_sol* qp_sol, struct d_dense_qp_ipm_ws* ws) {

    int ii, idx;

    int nv = qp->dim->nv;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;
    int ns = qp->dim->ns;

    int* idxs_rev = qp->idxs_rev;

    //	struct vec *dv = ws->sol_step->v;
    struct vec* dv = qp_sol->v;

    //	struct vec *res_g = ws->res->res_g;
    struct vec* res_g = qp->gz;

    struct vec* Gamma = ws->Gamma;
    struct vec* gamma = ws->gamma;
    struct vec* Zs_inv = ws->Zs_inv;
    struct vec* tmp_nbg = ws->tmp_nbg;

    double* ptr_Gamma = Gamma->pa;
    double* ptr_gamma = gamma->pa;
    double* ptr_Zs_inv = Zs_inv->pa;
    double* ptr_dv = dv->pa;
    double* ptr_res_g = res_g->pa;
    double* ptr_tmp2 = (tmp_nbg + 2)->pa;
    double* ptr_tmp3 = (tmp_nbg + 3)->pa;

    double tmp0, tmp1;

    dveccp(nb + ng, gamma, 0, tmp_nbg + 2, 0);
    dveccp(nb + ng, gamma, nb + ng, tmp_nbg + 3, 0);

    // idxs_rev
    for (ii = 0; ii < nb + ng; ii++) {
        idx = idxs_rev[ii];
        if (idx != -1) {
            // ii  <= constr index
            // idx <= slack index
            ptr_dv[nv + idx] = ptr_res_g[nv + idx] + ptr_gamma[0 + ii] + ptr_gamma[2 * nb + 2 * ng + idx];
            ptr_dv[nv + ns + idx] = ptr_res_g[nv + ns + idx] + ptr_gamma[nb + ng + ii] + ptr_gamma[2 * nb + 2 * ng + ns + idx];
            tmp0 = ptr_dv[nv + idx] * ptr_Zs_inv[0 + idx];
            tmp1 = ptr_dv[nv + ns + idx] * ptr_Zs_inv[ns + idx];
            ptr_tmp2[ii] = ptr_tmp2[ii] - ptr_Gamma[0 + ii] * tmp0;
            ptr_tmp3[ii] = ptr_tmp3[ii] - ptr_Gamma[nb + ng + ii] * tmp1;
        }
    }

    daxpy(nb + ng, -1.0, tmp_nbg + 3, 0, tmp_nbg + 2, 0, tmp_nbg + 1, 0);
}


static void EXPAND_SLACKS(struct d_dense_qp* qp, struct d_dense_qp_sol* qp_sol, struct d_dense_qp_ipm_ws* ws) {

    int ii, idx;

    int nv = qp->dim->nv;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;
    int ns = qp->dim->ns;

    int* idxs_rev = qp->idxs_rev;

    struct vec* dv = qp_sol->v;
    struct vec* dt = qp_sol->t;

    struct vec* Gamma = ws->Gamma;
    struct vec* Zs_inv = ws->Zs_inv;

    double* ptr_Gamma = Gamma->pa;
    double* ptr_dv = dv->pa;
    double* ptr_dt = dt->pa;
    double* ptr_Zs_inv = Zs_inv->pa;

    // idxs_rev
    for (ii = 0; ii < nb + ng; ii++) {
        idx = idxs_rev[ii];
        if (idx != -1) {
            // ii  <= constr index
            // idx <= slack index
            ptr_dv[nv + idx] = -ptr_Zs_inv[0 + idx] * (ptr_dv[nv + idx] + ptr_dt[ii] * ptr_Gamma[ii]);
            ptr_dv[nv + ns + idx] = -ptr_Zs_inv[ns + idx] * (ptr_dv[nv + ns + idx] + ptr_dt[nb + ng + ii] * ptr_Gamma[nb + ng + ii]);
            ptr_dt[2 * nb + 2 * ng + idx] = ptr_dv[nv + idx];
            ptr_dt[2 * nb + 2 * ng + ns + idx] = ptr_dv[nv + ns + idx];
            ptr_dt[0 + ii] = ptr_dt[0 + ii] + ptr_dv[nv + idx];
            ptr_dt[nb + ng + ii] = ptr_dt[nb + ng + ii] + ptr_dv[nv + ns + idx];
        }
    }
}


// range-space (Schur complement) method
void d_fact_solve_kkt_step_dense_qp(struct d_dense_qp* qp, struct d_dense_qp_sol* qp_sol, struct d_dense_qp_ipm_arg* arg, struct d_dense_qp_ipm_ws* ws) {

    int ii, jj;

    int nv = qp->dim->nv;
    int ne = qp->dim->ne;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;
    int ns = qp->dim->ns;

    struct mat* Hg = qp->Hv;
    struct mat* A = qp->A;
    struct mat* Ct = qp->Ct;
    int* idxb = qp->idxb;

    struct vec* res_g = qp->gz;
    struct vec* res_b = qp->b;

    struct vec* dv = qp_sol->v;
    struct vec* dpi = qp_sol->pi;
    struct vec* dt = qp_sol->t;

    struct mat* Lv = ws->Lv;
    struct mat* Le = ws->Le;
    struct mat* Ctx = ws->Ctx;
    struct mat* AL = ws->AL;
    struct vec* lv = ws->lv;
    struct vec* sv = ws->sv;
    struct vec* se = ws->se;
    struct vec* Gamma = ws->Gamma;
    struct vec* gamma = ws->gamma;
    struct vec* tmp_nbg = ws->tmp_nbg;

    // null space
    struct mat* A_LQ = ws->A_LQ;
    struct mat* A_Q = ws->A_Q;
    struct mat* Zt = ws->Zt;
    struct mat* ZtH = ws->ZtH;
    struct mat* ZtHZ = ws->ZtHZ;
    struct vec* xy = ws->xy;
    struct vec* Yxy = ws->Yxy;
    struct vec* xz = ws->xz;
    struct vec* tmp_nv = ws->tmp_nv;
    void* lq_work_null = ws->lq_work_null;
    void* orglq_work_null = ws->orglq_work_null;

    double tmp;

    struct d_core_qp_ipm_workspace* cws = ws->core_workspace;

    if (nb + ng > 0) {
        d_compute_Gamma_gamma_qp(qp->d->pa, qp->m->pa, cws);
    }

    if (ne > 0) {

        if (arg->kkt_fact_alg == 0)  // null space method
        {

            if (ws->use_A_fact == 0) {
                dgelqf(ne, nv, A, 0, 0, A_LQ, 0, 0, lq_work_null);

                // TODO cache dA containing tau into another vector !!!!!

                // TODO change dorglq API to pass tau explicitly as a vector !!!!!
                // TODO allocate its dedicated workspace !!!!!
                dorglq(nv, nv, ne, A_LQ, 0, 0, A_Q, 0, 0, orglq_work_null);

                dgecp(nv - ne, nv, A_Q, ne, 0, Zt, 0, 0);

#if 0
printf("\nA\n");
blasfeo_print_dmat(ne, nv, A, 0, 0);
printf("\nA_LQ\n");
blasfeo_print_dmat(ne, nv, A_LQ, 0, 0);
printf("\nA_Q\n");
blasfeo_print_dmat(nv, nv, A_Q, 0, 0);

for(ii=0; ii<ne; ii++)
	for(jj=ii+1; jj<nv; jj++)
		MATEL(A_LQ, ii, jj) = 0;
blasfeo_dgemm_nn(ne, nv, nv, 1.0, A_LQ, 0, 0, A_Q, 0, 0, -1.0, A, 0, 0, AL, 0, 0);
double max_err = 0.0;
double tmp;
for(ii=0; ii<ne; ii++)
	for(jj=ii+1; jj<nv; jj++)
		{
		tmp = MATEL(AL, ii, jj);
		max_err = fabs(tmp)>max_err ? fabs(tmp) : max_err;
		}
printf("\nA_LQ * A_Q - A max err %e\n", max_err);
//blasfeo_print_exp_dmat(ne, nv, AL, 0, 0);
//exit(1);
#endif
                ws->use_A_fact = 1;
            }


            dtrcp_l(nv, Hg, 0, 0, Lv, 0, 0);
            //			dgecp(nv, nv, Hg, 0, 0, Lv, 0, 0);
            //			ddiare(nv, arg->reg_prim, Lv, 0, 0); // XXX leave out ???

            dveccp(nv, res_g, 0, lv, 0);

            if (ns > 0) {
                COND_SLACKS_FACT_solve(qp, qp_sol, arg, ws);
            } else if (nb + ng > 0) {
                daxpy(nb + ng, 1.0, Gamma, nb + ng, Gamma, 0, tmp_nbg + 0, 0);
                daxpy(nb + ng, -1.0, gamma, nb + ng, gamma, 0, tmp_nbg + 1, 0);
            }
            if (nb > 0) {
                ddiaad_sp(nb, 1.0, tmp_nbg + 0, 0, idxb, Lv, 0, 0);
                dvecad_sp(nb, 1.0, tmp_nbg + 1, 0, idxb, lv, 0);
            }
            if (ng > 0) {
                dgemv_n(nv, ng, 1.0, Ct, 0, 0, tmp_nbg + 1, nb, 1.0, lv, 0, lv, 0);
                dgemm_nd(nv, ng, 1.0, Ct, 0, 0, tmp_nbg + 0, nb, 0.0, Ctx, 0, 0, Ctx, 0, 0);
                dsyrk_ln(nv, ng, 1.0, Ctx, 0, 0, Ct, 0, 0, 1.0, Lv, 0, 0, Lv, 0, 0);
            }

            dtrtr_l(nv, Lv, 0, 0, Lv, 0, 0);

            dgemm_nt(nv - ne, nv, nv, 1.0, Zt, 0, 0, Lv, 0, 0, 0.0, ZtH, 0, 0, ZtH, 0, 0);
            dsyrk_ln(nv - ne, nv, 1.0, ZtH, 0, 0, Zt, 0, 0, 0.0, ZtHZ, 0, 0, ZtHZ, 0, 0);
            ddiare(nv - ne, arg->reg_prim, ZtHZ, 0, 0);  // XXX leave in ???
            dpotrf_l(nv - ne, ZtHZ, 0, 0, ZtHZ, 0, 0);

            dtrsv_lnn(ne, A_LQ, 0, 0, res_b, 0, xy, 0);
            // printf("\nxy\n");
            // blasfeo_print_tran_dvec(ne, xy, 0);
            dgemv_t(ne, nv, 1.0, A_Q, 0, 0, xy, 0, 0.0, Yxy, 0, Yxy, 0);

            dgemv_n(nv - ne, nv, -1.0, ZtH, 0, 0, Yxy, 0, 0.0, xz, 0, xz, 0);
            dgemv_n(nv - ne, nv, -1.0, Zt, 0, 0, lv, 0, 1.0, xz, 0, xz, 0);
            dtrsv_lnn(nv - ne, ZtHZ, 0, 0, xz, 0, xz, 0);
            dtrsv_ltn(nv - ne, ZtHZ, 0, 0, xz, 0, xz, 0);

            dgemv_t(nv - ne, nv, 1.0, Zt, 0, 0, xz, 0, 1.0, Yxy, 0, dv, 0);

            dsymv_l(nv, 1.0, Lv, 0, 0, dv, 0, 1.0, lv, 0, tmp_nv, 0);
            dgemv_n(ne, nv, 1.0, A_Q, 0, 0, tmp_nv, 0, 0.0, dpi, 0, dpi, 0);
            dtrsv_ltn(ne, A_LQ, 0, 0, dpi, 0, dpi, 0);

        } else  // schur-complement method
        {

            if (arg->scale) {

                //			dtrcp_l(nv, Hg, 0, 0, Lv, 0, 0);
                dgecp(nv, nv, Hg, 0, 0, Lv, 0, 0);

                dveccp(nv, res_g, 0, lv, 0);

                if (ns > 0) {
                    COND_SLACKS_FACT_solve(qp, qp_sol, arg, ws);
                } else if (nb + ng > 0) {
                    daxpy(nb + ng, 1.0, Gamma, nb + ng, Gamma, 0, tmp_nbg + 0, 0);
                    daxpy(nb + ng, -1.0, gamma, nb + ng, gamma, 0, tmp_nbg + 1, 0);
                }
                if (nb > 0) {
                    ddiaad_sp(nb, 1.0, tmp_nbg + 0, 0, idxb, Lv, 0, 0);
                    dvecad_sp(nb, 1.0, tmp_nbg + 1, 0, idxb, lv, 0);
                }
                if (ng > 0) {
                    dgemv_n(nv, ng, 1.0, Ct, 0, 0, tmp_nbg + 1, nb, 1.0, lv, 0, lv, 0);
                    dgemm_nd(nv, ng, 1.0, Ct, 0, 0, tmp_nbg + 0, nb, 0.0, Ctx, 0, 0, Ctx, 0, 0);
                    dsyrk_ln(nv, ng, 1.0, Ctx, 0, 0, Ct, 0, 0, 1.0, Lv, 0, 0, Lv, 0, 0);
                }

                ddiaex(nv, 1.0, Lv, 0, 0, sv, 0);
                for (ii = 0; ii < nv; ii++) {
                    tmp = sqrt(sv->pa[ii]);
                    //				tmp = sqrt(tmp);
                    //				tmp = sqrt(sv->pa[ii]+tmp);
                    //				tmp = 1.0;
                    sv->pa[ii] = tmp == 0 ? 1.0 : 1.0 / tmp;
                }

                dgemm_dn(nv, nv, 1.0, sv, 0, Lv, 0, 0, 0.0, Lv, 0, 0, Lv, 0, 0);
                dgemm_nd(nv, nv, 1.0, Lv, 0, 0, sv, 0, 0.0, Lv, 0, 0, Lv, 0, 0);
                ddiare(nv, arg->reg_prim, Lv, 0, 0);
                dpotrf_l(nv, Lv, 0, 0, Lv, 0, 0);

                dgemv_d(nv, 1.0, sv, 0, lv, 0, 0.0, lv, 0, lv, 0);
                dveccp(nv, lv, 0, dv, 0);

                dgecp(ne, nv, A, 0, 0, AL, 0, 0);
                dgemm_nd(ne, nv, 1.0, AL, 0, 0, sv, 0, 0.0, AL, 0, 0, AL, 0, 0);
                dtrsm_rltn(ne, nv, 1.0, Lv, 0, 0, AL, 0, 0, AL, 0, 0);

                dtrsv_lnn(nv, Lv, 0, 0, lv, 0, lv, 0);

                dgese(ne, ne, 0.0, Le, 0, 0);
                dsyrk_ln(ne, nv, 1.0, AL, 0, 0, AL, 0, 0, 1.0, Le, 0, 0, Le, 0, 0);

                ddiaex(ne, 1.0, Le, 0, 0, se, 0);
                for (ii = 0; ii < ne; ii++) {
                    tmp = sqrt(se->pa[ii]);
                    //				tmp = sqrt(tmp);
                    //				tmp = sqrt(se->pa[ii]+tmp);
                    //				tmp = 1.0;
                    se->pa[ii] = tmp == 0 ? 1.0 : 1.0 / tmp;
                }

                dgemm_dn(ne, ne, 1.0, se, 0, Le, 0, 0, 0.0, Le, 0, 0, Le, 0, 0);
                dgemm_nd(ne, ne, 1.0, Le, 0, 0, se, 0, 0.0, Le, 0, 0, Le, 0, 0);
                ddiare(ne, arg->reg_dual, Le, 0, 0);
                dpotrf_l(ne, Le, 0, 0, Le, 0, 0);

                dgemv_n(ne, nv, 1.0, AL, 0, 0, lv, 0, 1.0, res_b, 0, dpi, 0);

                dgemv_d(ne, 1.0, se, 0, dpi, 0, 0.0, dpi, 0, dpi, 0);
                dtrsv_lnn(ne, Le, 0, 0, dpi, 0, dpi, 0);
                dtrsv_ltn(ne, Le, 0, 0, dpi, 0, dpi, 0);
                dgemv_d(ne, 1.0, se, 0, dpi, 0, 0.0, dpi, 0, dpi, 0);

                dgemv_t(ne, nv, 1.0, A, 0, 0, dpi, 0, 0.0, lv, 0, lv, 0);
                dgemv_d(nv, 1.0, sv, 0, lv, 0, -1.0, dv, 0, dv, 0);

                dtrsv_lnn(nv, Lv, 0, 0, dv, 0, dv, 0);
                dtrsv_ltn(nv, Lv, 0, 0, dv, 0, dv, 0);
                dgemv_d(nv, 1.0, sv, 0, dv, 0, 0.0, dv, 0, dv, 0);

            } else  // no scale
            {

                //			dtrcp_l(nv, Hg, 0, 0, Lv, 0, 0);
                dgecp(nv, nv, Hg, 0, 0, Lv, 0, 0);
                ddiare(nv, arg->reg_prim, Lv, 0, 0);

                dveccp(nv, res_g, 0, lv, 0);

                if (ns > 0) {
                    COND_SLACKS_FACT_solve(qp, qp_sol, arg, ws);
                } else if (nb + ng > 0) {
                    daxpy(nb + ng, 1.0, Gamma, nb + ng, Gamma, 0, tmp_nbg + 0, 0);
                    daxpy(nb + ng, -1.0, gamma, nb + ng, gamma, 0, tmp_nbg + 1, 0);
                }
                if (nb > 0) {
                    ddiaad_sp(nb, 1.0, tmp_nbg + 0, 0, idxb, Lv, 0, 0);
                    dvecad_sp(nb, 1.0, tmp_nbg + 1, 0, idxb, lv, 0);
                }
                if (ng > 0) {
                    dgemv_n(nv, ng, 1.0, Ct, 0, 0, tmp_nbg + 1, nb, 1.0, lv, 0, lv, 0);
                    dgemm_nd(nv, ng, 1.0, Ct, 0, 0, tmp_nbg + 0, nb, 0.0, Ctx, 0, 0, Ctx, 0, 0);
                    dsyrk_dpotrf_ln(nv, ng, Ctx, 0, 0, Ct, 0, 0, Lv, 0, 0, Lv, 0, 0);
                } else {
                    dpotrf_l(nv, Lv, 0, 0, Lv, 0, 0);
                }
                // int pd = 1;
                // for(ii=0; ii<nv; ii++)
                //	if(Lv->dA[ii]==0.0)
                //		pd = 0;
                // printf(" chol pd %d\n", pd);

                dveccp(nv, lv, 0, dv, 0);

                dtrsm_rltn(ne, nv, 1.0, Lv, 0, 0, A, 0, 0, AL, 0, 0);

                dtrsv_lnn(nv, Lv, 0, 0, lv, 0, lv, 0);

                dgemv_n(ne, nv, 1.0, AL, 0, 0, lv, 0, 1.0, res_b, 0, dpi, 0);

                dgese(ne, ne, 0.0, Le, 0, 0);
                ddiare(ne, arg->reg_dual, Le, 0, 0);
                dsyrk_dpotrf_ln(ne, nv, AL, 0, 0, AL, 0, 0, Le, 0, 0, Le, 0, 0);

                dtrsv_lnn(ne, Le, 0, 0, dpi, 0, dpi, 0);
                dtrsv_ltn(ne, Le, 0, 0, dpi, 0, dpi, 0);

                dgemv_t(ne, nv, 1.0, A, 0, 0, dpi, 0, -1.0, dv, 0, dv, 0);

                dtrsv_lnn(nv, Lv, 0, 0, dv, 0, dv, 0);
                dtrsv_ltn(nv, Lv, 0, 0, dv, 0, dv, 0);


            }  // scale
        }

    } else  // ne==0
    {

        if (arg->scale) {

            dtrcp_l(nv, Hg, 0, 0, Lv, 0, 0);
            dveccp(nv, res_g, 0, lv, 0);

            if (ns > 0) {
                COND_SLACKS_FACT_solve(qp, qp_sol, arg, ws);
            } else if (nb + ng > 0) {
                daxpy(nb + ng, 1.0, Gamma, nb + ng, Gamma, 0, tmp_nbg + 0, 0);
                daxpy(nb + ng, -1.0, gamma, nb + ng, gamma, 0, tmp_nbg + 1, 0);
            }
            if (nb > 0) {
                ddiaad_sp(nb, 1.0, tmp_nbg + 0, 0, idxb, Lv, 0, 0);
                dvecad_sp(nb, 1.0, tmp_nbg + 1, 0, idxb, lv, 0);
            }
            if (ng > 0) {
                dgemm_nd(nv, ng, 1.0, Ct, 0, 0, tmp_nbg + 0, nb, 0.0, Ctx, 0, 0, Ctx, 0, 0);
                dgemv_n(nv, ng, 1.0, Ct, 0, 0, tmp_nbg + 1, nb, 1.0, lv, 0, lv, 0);
                dsyrk_ln(nv, ng, 1.0, Ctx, 0, 0, Ct, 0, 0, 1.0, Lv, 0, 0, Lv, 0, 0);
            }

            ddiaex(nv, 1.0, Lv, 0, 0, sv, 0);
            for (ii = 0; ii < nv; ii++) {
                tmp = sqrt(sv->pa[ii]);
                //				tmp = sqrt(tmp);
                //				tmp = sqrt(sv->pa[ii]+tmp);
                //				tmp = 1.0;
                sv->pa[ii] = tmp == 0 ? 1.0 : 1.0 / tmp;
            }

            dgemm_dn(nv, nv, 1.0, sv, 0, Lv, 0, 0, 0.0, Lv, 0, 0, Lv, 0, 0);
            dgemm_nd(nv, nv, 1.0, Lv, 0, 0, sv, 0, 0.0, Lv, 0, 0, Lv, 0, 0);
            ddiare(nv, arg->reg_prim, Lv, 0, 0);
            dpotrf_l_mn(nv, nv, Lv, 0, 0, Lv, 0, 0);

            dveccp(nv, lv, 0, dv, 0);
            dvecsc(nv, -1.0, dv, 0);

            dgemv_d(nv, 1.0, sv, 0, dv, 0, 0.0, dv, 0, dv, 0);
            dtrsv_lnn(nv, Lv, 0, 0, dv, 0, dv, 0);
            dtrsv_ltn(nv, Lv, 0, 0, dv, 0, dv, 0);
            dgemv_d(nv, 1.0, sv, 0, dv, 0, 0.0, dv, 0, dv, 0);
        } else  // no scale
        {
            //		dtrcp_l(nv, Hg, 0, 0, Lv, 0, 0);
            dgecp(nv, nv, Hg, 0, 0, Lv, 0, 0);
            drowin(nv, 1.0, res_g, 0, Lv, nv, 0);
            ddiare(nv, arg->reg_prim, Lv, 0, 0);

            if (ns > 0) {
                COND_SLACKS_FACT_solve(qp, qp_sol, arg, ws);
            } else if (nb + ng > 0) {
                daxpy(nb + ng, 1.0, Gamma, nb + ng, Gamma, 0, tmp_nbg + 0, 0);
                daxpy(nb + ng, -1.0, gamma, nb + ng, gamma, 0, tmp_nbg + 1, 0);
            }
            if (nb > 0) {
                ddiaad_sp(nb, 1.0, tmp_nbg + 0, 0, idxb, Lv, 0, 0);
                drowad_sp(nb, 1.0, tmp_nbg + 1, 0, idxb, Lv, nv, 0);
            }
            if (ng > 0) {
                dgemm_nd(nv, ng, 1.0, Ct, 0, 0, tmp_nbg + 0, nb, 0.0, Ctx, 0, 0, Ctx, 0, 0);
                drowin(ng, 1.0, tmp_nbg + 1, nb, Ctx, nv, 0);
                dsyrk_dpotrf_ln_mn(nv + 1, nv, ng, Ctx, 0, 0, Ct, 0, 0, Lv, 0, 0, Lv, 0, 0);
            } else {
                dpotrf_l_mn(nv + 1, nv, Lv, 0, 0, Lv, 0, 0);
            }

            drowex(nv, -1.0, Lv, nv, 0, dv, 0);
            dtrsv_ltn(nv, Lv, 0, 0, dv, 0, dv, 0);

        }  // scale

    }  // ne>0

    // blasfeo_print_tran_dvec(nv, dv, 0);
    // blasfeo_print_tran_dvec(ne, dpi, 0);
    // exit(1);
    if (nb + ng > 0) {
        if (nb > 0)
            dvecex_sp(nb, 1.0, idxb, dv, 0, dt, 0);

        if (ng > 0)
            dgemv_t(nv, ng, 1.0, Ct, 0, 0, dv, 0, 0.0, dt, nb, dt, nb);

        dveccp(nb + ng, dt, 0, dt, nb + ng);
        dvecsc(nb + ng, -1.0, dt, nb + ng);

        if (ns > 0)
            EXPAND_SLACKS(qp, qp_sol, ws);

        d_compute_lam_t_qp(qp->d->pa, qp->m->pa, qp_sol->lam->pa, qp_sol->t->pa, cws);
    }
}


void d_fact_lq_solve_kkt_step_dense_qp(struct d_dense_qp* qp, struct d_dense_qp_sol* qp_sol, struct d_dense_qp_ipm_arg* arg, struct d_dense_qp_ipm_ws* ws) {

    int ii;

    int nv = qp->dim->nv;
    int ne = qp->dim->ne;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;
    int ns = qp->dim->ns;

    struct mat* Hg = qp->Hv;
    struct mat* A = qp->A;
    struct mat* Ct = qp->Ct;
    int* idxb = qp->idxb;
    struct vec* res_g = qp->gz;
    struct vec* res_b = qp->b;

    struct vec* dv = qp_sol->v;
    struct vec* dpi = qp_sol->pi;
    struct vec* dt = qp_sol->t;

    struct mat* Lv = ws->Lv;
    struct mat* Le = ws->Le;
    struct mat* Ctx = ws->Ctx;
    struct mat* AL = ws->AL;
    struct vec* lv = ws->lv;
    struct vec* sv = ws->sv;
    struct vec* se = ws->se;
    struct vec* Gamma = ws->Gamma;
    struct vec* gamma = ws->gamma;
    struct vec* tmp_nbg = ws->tmp_nbg;
    void* lq_work0 = ws->lq_work0;
    void* lq_work1 = ws->lq_work1;
    struct mat* lq0 = ws->lq0;
    struct mat* lq1 = ws->lq1;

    // null space
    struct mat* A_LQ = ws->A_LQ;
    struct mat* A_Q = ws->A_Q;
    struct mat* Zt = ws->Zt;
    struct mat* ZtH = ws->ZtH;
    struct mat* ZtHZ = ws->ZtHZ;
    struct vec* xy = ws->xy;
    struct vec* Yxy = ws->Yxy;
    struct vec* xz = ws->xz;
    struct vec* tmp_nv = ws->tmp_nv;
    void* lq_work_null = ws->lq_work_null;
    void* orglq_work_null = ws->orglq_work_null;

    double tmp;

    struct d_core_qp_ipm_workspace* cws = ws->core_workspace;

    ws->scale = 0;

    if (nb + ng > 0) {
        d_compute_Gamma_gamma_qp(qp->d->pa, qp->m->pa, cws);
    }

    if (ne > 0) {

        if (arg->kkt_fact_alg == 0)  // null space method
        {

            if (ws->use_A_fact == 0) {
                dgelqf(ne, nv, A, 0, 0, A_LQ, 0, 0, lq_work_null);

                // TODO cache dA containing tau into another vector !!!!!

                // TODO change dorglq API to pass tau explicitly as a vector !!!!!
                // TODO allocate its dedicated workspace !!!!!
                dorglq(nv, nv, ne, A_LQ, 0, 0, A_Q, 0, 0, orglq_work_null);

                dgecp(nv - ne, nv, A_Q, ne, 0, Zt, 0, 0);

                ws->use_A_fact = 1;
            }


            dtrcp_l(nv, Hg, 0, 0, Lv, 0, 0);
            //			dgecp(nv, nv, Hg, 0, 0, Lv, 0, 0);

            dveccp(nv, res_g, 0, lv, 0);

            if (ns > 0) {
                COND_SLACKS_FACT_solve(qp, qp_sol, arg, ws);
            } else if (nb + ng > 0) {
                daxpy(nb + ng, 1.0, Gamma, nb + ng, Gamma, 0, tmp_nbg + 0, 0);
                daxpy(nb + ng, -1.0, gamma, nb + ng, gamma, 0, tmp_nbg + 1, 0);
            }
            if (nb > 0) {
                ddiaad_sp(nb, 1.0, tmp_nbg + 0, 0, idxb, Lv, 0, 0);
                dvecad_sp(nb, 1.0, tmp_nbg + 1, 0, idxb, lv, 0);
            }
            if (ng > 0) {
                dgemv_n(nv, ng, 1.0, Ct, 0, 0, tmp_nbg + 1, nb, 1.0, lv, 0, lv, 0);
                dgemm_nd(nv, ng, 1.0, Ct, 0, 0, tmp_nbg + 0, nb, 0.0, Ctx, 0, 0, Ctx, 0, 0);
                dsyrk_ln(nv, ng, 1.0, Ctx, 0, 0, Ct, 0, 0, 1.0, Lv, 0, 0, Lv, 0, 0);
            }

            dtrtr_l(nv, Lv, 0, 0, Lv, 0, 0);

            dgemm_nt(nv - ne, nv, nv, 1.0, Zt, 0, 0, Lv, 0, 0, 0.0, ZtH, 0, 0, ZtH, 0, 0);
            dsyrk_ln(nv - ne, nv, 1.0, ZtH, 0, 0, Zt, 0, 0, 0.0, ZtHZ, 0, 0, ZtHZ, 0, 0);
            dpotrf_l(nv - ne, ZtHZ, 0, 0, ZtHZ, 0, 0);

            dtrsv_lnn(ne, A_LQ, 0, 0, res_b, 0, xy, 0);
            dgemv_t(ne, nv, 1.0, A_Q, 0, 0, xy, 0, 0.0, Yxy, 0, Yxy, 0);

            dgemv_n(nv - ne, nv, -1.0, ZtH, 0, 0, Yxy, 0, 0.0, xz, 0, xz, 0);
            dgemv_n(nv - ne, nv, -1.0, Zt, 0, 0, lv, 0, 1.0, xz, 0, xz, 0);
            dtrsv_lnn(nv - ne, ZtHZ, 0, 0, xz, 0, xz, 0);
            dtrsv_ltn(nv - ne, ZtHZ, 0, 0, xz, 0, xz, 0);

            dgemv_t(nv - ne, nv, 1.0, Zt, 0, 0, xz, 0, 1.0, Yxy, 0, dv, 0);

            dsymv_l(nv, 1.0, Lv, 0, 0, dv, 0, 1.0, lv, 0, tmp_nv, 0);
            dgemv_n(ne, nv, 1.0, A_Q, 0, 0, tmp_nv, 0, 0.0, dpi, 0, dpi, 0);
            dtrsv_ltn(ne, A_LQ, 0, 0, dpi, 0, dpi, 0);

        } else  // schur-complement method
        {

            // XXX needed ???
            dgese(nv, nv + nv + ng, 0.0, lq1, 0, 0);  // TODO not the first part for HP and RF

            if (ws->use_hess_fact == 0) {
                dtrcp_l(nv, Hg, 0, 0, Lv + 1, 0, 0);
                ddiare(nv, arg->reg_prim, Lv + 1, 0, 0);
                dpotrf_l(nv, Lv + 1, 0, 0, Lv + 1, 0, 0);
                ws->use_hess_fact = 1;
            }
            // int pd = 1;
            // for(ii=0; ii<nv; ii++)
            //	if((Lv+1)->dA[ii]==0.0)
            //		pd = 0;
            // printf(" lq pd %d\n", pd);

            dveccp(nv, res_g, 0, lv, 0);

            if (ns > 0) {
                COND_SLACKS_FACT_solve(qp, qp_sol, arg, ws);
            } else if (nb + ng > 0) {
                daxpy(nb + ng, 1.0, Gamma, nb + ng, Gamma, 0, tmp_nbg + 0, 0);
                daxpy(nb + ng, -1.0, gamma, nb + ng, gamma, 0, tmp_nbg + 1, 0);
            }
            if (nb > 0) {
                for (ii = 0; ii < nb; ii++) {
                    tmp = VECEL(tmp_nbg + 0, ii);
                    tmp = tmp >= 0.0 ? tmp : 0.0;
                    tmp = sqrt(tmp);
                    MATEL(lq1, idxb[ii], nv + idxb[ii]) = tmp > 0.0 ? tmp : 0.0;
                }
                dvecad_sp(nb, 1.0, tmp_nbg + 1, 0, idxb, lv, 0);
            }
            if (ng > 0) {
                for (ii = 0; ii < ng; ii++) {
                    tmp = VECEL(tmp_nbg + 0, nb + ii);
                    tmp = tmp >= 0.0 ? tmp : 0.0;
                    tmp = sqrt(tmp);
                    VECEL(tmp_nbg + 0, nb + ii) = tmp;
                }
                dgemm_nd(nv, ng, 1.0, Ct, 0, 0, tmp_nbg + 0, nb, 0.0, lq1, 0, nv + nv, lq1, 0, nv + nv);
                dgemv_n(nv, ng, 1.0, Ct, 0, 0, tmp_nbg + 1, nb, 1.0, lv, 0, lv, 0);
            }

            //		ddiare(nv, arg->reg_prim, lq1, 0, nv);

            // blasfeo_print_dmat(nv, nv, lq1, 0, 0);
            // blasfeo_print_dmat(nv, nv, lq1, 0, nv);
            // blasfeo_print_dmat(nv, ng, lq1, 0, nv+ng);

#if defined(LA_HIGH_PERFORMANCE) | defined(LA_REFERENCE)
            //		dtrcp_l(nv, Lv+1, 0, 0, lq1, 0, 0);
            //		dgelqf_pd(nv, nv+nv+ng, lq1, 0, 0, lq1, 0, 0, lq_work1);
            //		dgelqf_pd_la(nv, nv+ng, lq1, 0, 0, lq1, 0, nv, lq_work1);
            //		dgelqf_pd_lla(nv, ng, lq1, 0, 0, lq1, 0, nv, lq1, 0, 2*nv, lq_work1);
            //		dtrcp_l(nv, lq1, 0, 0, Lv, 0, 0);
            dtrcp_l(nv, Lv + 1, 0, 0, Lv, 0, 0);
            dgelqf_pd_lla(nv, ng, Lv, 0, 0, lq1, 0, nv, lq1, 0, 2 * nv, lq_work1);  // TODO reduce lq1 size !!!
#else  // LA_BLAS_WRAPPER
            dtrcp_l(nv, Lv + 1, 0, 0, lq1, 0, 0);
            dgelqf(nv, nv + nv + ng, lq1, 0, 0, lq1, 0, 0, lq_work1);
            dtrcp_l(nv, lq1, 0, 0, Lv, 0, 0);
            for (ii = 0; ii < nv; ii++)
                if (MATEL(Lv, ii, ii) < 0)
                    dcolsc(nv - ii, -1.0, Lv, ii, ii);
#endif

            // blasfeo_print_dmat(nv, nv, Lv, 0, 0);

            dveccp(nv, lv, 0, dv, 0);

            dtrsm_rltn(ne, nv, 1.0, Lv, 0, 0, A, 0, 0, AL, 0, 0);

            dtrsv_lnn(nv, Lv, 0, 0, lv, 0, lv, 0);

            dgemv_n(ne, nv, 1.0, AL, 0, 0, lv, 0, 1.0, res_b, 0, dpi, 0);

            dgecp(ne, nv, AL, 0, 0, lq0, 0, ne);

#if defined(LA_HIGH_PERFORMANCE)
            //		dgese(ne, ne, 0.0, lq0, 0, 0);
            //		ddiare(ne, arg->reg_dual, lq0, 0, 0);
            //		dgelqf_pd(ne, ne+nv, lq0, 0, 0, lq0, 0, 0, lq_work0);
            //		dgelqf_pd_la(ne, nv, lq0, 0, 0, lq0, 0, ne, lq_work0);
            //		dtrcp_l(ne, lq0, 0, 0, Le, 0, 0);
            dgese(ne, ne, 0.0, Le, 0, 0);
            ddiare(ne, arg->reg_dual, Le, 0, 0);
            dgelqf_pd_la(ne, nv, Le, 0, 0, lq0, 0, ne, lq_work0);  // TODO reduce lq0 size !!!
#elif defined(LA_REFERENCE)
            //		dgese(ne, ne, 0.0, lq0, 0, 0);
            //		ddiare(ne, arg->reg_dual, lq0, 0, 0);
            //		dgelqf_pd(ne, ne+nv, lq0, 0, 0, lq0, 0, 0, lq_work0);
            //		dgelqf_pd_la(ne, nv, lq0, 0, 0, lq0, 0, ne, lq_work0);
            //		dtrcp_l(ne, lq0, 0, 0, Le, 0, 0);
            dgese(ne, ne, 0.0, Le, 0, 0);
            ddiare(ne, arg->reg_dual, Le, 0, 0);
            dgelqf_pd_la(ne, nv, Le, 0, 0, lq0, 0, ne, lq_work0);  // TODO reduce lq0 size !!!
#else  // LA_BLAS_WRAPPER
            dgese(ne, ne, 0.0, lq0, 0, 0);
            ddiare(ne, arg->reg_dual, lq0, 0, 0);
            dgelqf(ne, ne + nv, lq0, 0, 0, lq0, 0, 0, lq_work0);
            dtrcp_l(ne, lq0, 0, 0, Le, 0, 0);
            for (ii = 0; ii < ne; ii++)
                if (MATEL(Le, ii, ii) < 0)
                    dcolsc(ne - ii, -1.0, Le, ii, ii);
#endif

            //		blasfeo_print_dmat(ne, ne, Le, 0, 0);

            dtrsv_lnn(ne, Le, 0, 0, dpi, 0, dpi, 0);
            dtrsv_ltn(ne, Le, 0, 0, dpi, 0, dpi, 0);

            dgemv_t(ne, nv, 1.0, A, 0, 0, dpi, 0, -1.0, dv, 0, dv, 0);

            dtrsv_lnn(nv, Lv, 0, 0, dv, 0, dv, 0);
            dtrsv_ltn(nv, Lv, 0, 0, dv, 0, dv, 0);

        }  // schur-complement method

    } else  // ne==0
    {

        // XXX needed ???
        dgese(nv, nv + nv + ng, 0.0, lq1, 0, 0);  // TODO not the first part for HP and RF

        if (ws->use_hess_fact == 0) {
            dtrcp_l(nv, Hg, 0, 0, Lv + 1, 0, 0);
            ddiare(nv, arg->reg_prim, Lv + 1, 0, 0);
            dpotrf_l(nv, Lv + 1, 0, 0, Lv + 1, 0, 0);
            ws->use_hess_fact = 1;
        }

        dveccp(nv, res_g, 0, lv, 0);

        if (ns > 0) {
            COND_SLACKS_FACT_solve(qp, qp_sol, arg, ws);
        } else if (nb + ng > 0) {
            daxpy(nb + ng, 1.0, Gamma, nb + ng, Gamma, 0, tmp_nbg + 0, 0);
            daxpy(nb + ng, -1.0, gamma, nb + ng, gamma, 0, tmp_nbg + 1, 0);
        }
        if (nb > 0) {
            for (ii = 0; ii < nb; ii++) {
                tmp = VECEL(tmp_nbg + 0, ii);
                tmp = tmp >= 0.0 ? tmp : 0.0;
                tmp = sqrt(tmp);
                MATEL(lq1, idxb[ii], nv + idxb[ii]) = tmp > 0.0 ? tmp : 0.0;
            }
            dvecad_sp(nb, 1.0, tmp_nbg + 1, 0, idxb, lv, 0);
        }
        if (ng > 0) {
            for (ii = 0; ii < ng; ii++) {
                tmp = VECEL(tmp_nbg + 0, nb + ii);
                tmp = tmp >= 0.0 ? tmp : 0.0;
                tmp = sqrt(tmp);
                VECEL(tmp_nbg + 0, nb + ii) = tmp;
            }
            dgemm_nd(nv, ng, 1.0, Ct, 0, 0, tmp_nbg + 0, nb, 0.0, lq1, 0, nv + nv, lq1, 0, nv + nv);
            dgemv_n(nv, ng, 1.0, Ct, 0, 0, tmp_nbg + 1, nb, 1.0, lv, 0, lv, 0);
        }

        //		ddiare(nv, arg->reg_prim, lq1, 0, nv);

#if defined(LA_HIGH_PERFORMANCE)
        //		dtrcp_l(nv, Lv+1, 0, 0, lq1, 0, 0);
        //		dgelqf_pd(nv, nv+nv+ng, lq1, 0, 0, lq1, 0, 0, lq_work1);
        //		dgelqf_pd_la(nv, nv+ng, lq1, 0, 0, lq1, 0, nv, lq_work1);
        //		dgelqf_pd_lla(nv, ng, lq1, 0, 0, lq1, 0, nv, lq1, 0, 2*nv, lq_work1);
        //		dtrcp_l(nv, lq1, 0, 0, Lv, 0, 0);
        dtrcp_l(nv, Lv + 1, 0, 0, Lv, 0, 0);
        dgelqf_pd_lla(nv, ng, Lv, 0, 0, lq1, 0, nv, lq1, 0, 2 * nv, lq_work1);  // TODO reduce lq1 size !!!
#elif defined(LA_REFERENCE)
        //		dtrcp_l(nv, Lv+1, 0, 0, lq1, 0, 0);
        //		dgelqf_pd(nv, nv+nv+ng, lq1, 0, 0, lq1, 0, 0, lq_work1);
        //		dgelqf_pd_la(nv, nv+ng, lq1, 0, 0, lq1, 0, nv, lq_work1);
        //		dgelqf_pd_lla(nv, ng, lq1, 0, 0, lq1, 0, nv, lq1, 0, 2*nv, lq_work1);
        //		dtrcp_l(nv, lq1, 0, 0, Lv, 0, 0);
        dtrcp_l(nv, Lv + 1, 0, 0, Lv, 0, 0);
        dgelqf_pd_lla(nv, ng, Lv, 0, 0, lq1, 0, nv, lq1, 0, 2 * nv, lq_work1);  // TODO reduce lq1 size !!!
#else  // LA_BLAS_WRAPPER
        dtrcp_l(nv, Lv + 1, 0, 0, lq1, 0, 0);
        dgelqf(nv, nv + nv + ng, lq1, 0, 0, lq1, 0, 0, lq_work1);
        dtrcp_l(nv, lq1, 0, 0, Lv, 0, 0);
        for (ii = 0; ii < nv; ii++)
            if (MATEL(Lv, ii, ii) < 0)
                dcolsc(nv - ii, -1.0, Lv, ii, ii);
#endif

#if 0
if(nv<30)
{
//	blasfeo_print_dmat(nv, nv+nb+ng, lq1, 0, 0);
blasfeo_print_dmat(nv, nv, Lv, 0, 0);
exit(1);
}
#endif

        dveccp(nv, lv, 0, dv, 0);
        dvecsc(nv, -1.0, dv, 0);

        dtrsv_lnn(nv, Lv, 0, 0, dv, 0, dv, 0);
        dtrsv_ltn(nv, Lv, 0, 0, dv, 0, dv, 0);

    }  // ne>0

    if (nb + ng > 0) {
        if (nb > 0)
            dvecex_sp(nb, 1.0, idxb, dv, 0, dt, 0);

        dvecse(ng, 0.0, dt, nb);
        if (ng > 0)
            dgemv_t(nv, ng, 1.0, Ct, 0, 0, dv, 0, 0.0, dt, nb, dt, nb);

        dveccp(nb + ng, dt, 0, dt, nb + ng);
        dvecsc(nb + ng, -1.0, dt, nb + ng);

        if (ns > 0)
            EXPAND_SLACKS(qp, qp_sol, ws);

        d_compute_lam_t_qp(qp->d->pa, qp->m->pa, qp_sol->lam->pa, qp_sol->t->pa, cws);
    }
}


#if 0
void d_fact_solve_lu_kkt_step_dense_qp(struct d_dense_qp *qp, struct d_dense_qp_sol *qp_sol, struct d_dense_qp_ipm_arg *arg, struct d_dense_qp_ipm_ws *ws)
	{

	int ii;

	int nv = qp->dim->nv;
	int ne = qp->dim->ne;
	int nb = qp->dim->nb;
	int ng = qp->dim->ng;
	int ns = qp->dim->ns;

	struct mat *Hg = qp->Hv;
	struct mat *A = qp->A;
	struct mat *Ct = qp->Ct;
	int *idxb = qp->idxb;

	struct vec *res_g = qp->gz;
	struct vec *res_b = qp->b;

	struct vec *dv = qp_sol->v;
	struct vec *dpi = qp_sol->pi;
	struct vec *dt = qp_sol->t;

	struct mat *Lv = ws->Lv;
	struct mat *Le = ws->Le;
	struct mat *Ctx = ws->Ctx;
	struct mat *AL = ws->AL;
	struct vec *lv = ws->lv;
	struct vec *sv = ws->sv;
	struct vec *se = ws->se;
	struct vec *Gamma = ws->Gamma;
	struct vec *gamma = ws->gamma;
	struct vec *tmp_nbg = ws->tmp_nbg;
	int *ipiv_v = ws->ipiv_v;
	int *ipiv_e = ws->ipiv_e;

	double tmp;

	struct d_core_qp_ipm_workspace *cws = ws->core_workspace;

	if(nb+ng>0)
		{
		d_compute_Gamma_gamma_qp(qp->d->pa, qp->m->pa, cws);
		}

	if(ne>0)
		{

		if(arg->scale)
			{

//			dtrcp_l(nv, Hg, 0, 0, Lv, 0, 0);
			dgecp(nv, nv, Hg, 0, 0, Lv, 0, 0);

			dveccp(nv, res_g, 0, lv, 0);

			if(ns>0)
				{
				COND_SLACKS_FACT_solve(qp, qp_sol, arg, ws);
				}
			else if(nb+ng>0)
				{
				daxpy(nb+ng,  1.0, Gamma, nb+ng, Gamma, 0, tmp_nbg+0, 0);
				daxpy(nb+ng, -1.0, gamma, nb+ng, gamma, 0, tmp_nbg+1, 0);
				}
			if(nb>0)
				{
				ddiaad_sp(nb, 1.0, tmp_nbg+0, 0, idxb, Lv, 0, 0);
				dvecad_sp(nb, 1.0, tmp_nbg+1, 0, idxb, lv, 0);
				}
			if(ng>0)
				{
				dgemv_n(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+1, nb, 1.0, lv, 0, lv, 0);
				dgemm_nd(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+0, nb, 0.0, Ctx, 0, 0, Ctx, 0, 0);
				dsyrk_ln(nv, ng, 1.0, Ctx, 0, 0, Ct, 0, 0, 1.0, Lv, 0, 0, Lv, 0, 0);
				}

			ddiaex(nv, 1.0, Lv, 0, 0, sv, 0);
			for(ii=0; ii<nv; ii++)
				{
				tmp = sqrt(sv->pa[ii]);
//				tmp = sqrt(tmp);
//				tmp = sqrt(sv->pa[ii]+tmp);
//				tmp = 1.0;
				sv->pa[ii] = tmp==0 ? 1.0 : 1.0/tmp;
				}

			dgemm_dn(nv, nv, 1.0, sv, 0, Lv, 0, 0, 0.0, Lv, 0, 0, Lv, 0, 0);
			dgemm_nd(nv, nv, 1.0, Lv, 0, 0, sv, 0, 0.0, Lv, 0, 0, Lv, 0, 0);
			ddiare(nv, arg->reg_prim, Lv, 0, 0);
			dpotrf_l(nv, Lv, 0, 0, Lv, 0, 0);

			dgemv_d(nv, 1.0, sv, 0, lv, 0, 0.0, lv, 0, lv, 0);
			dveccp(nv, lv, 0, dv, 0);

			dgecp(ne, nv, A, 0, 0, AL, 0, 0);
			dgemm_nd(ne, nv, 1.0, AL, 0, 0, sv, 0, 0.0, AL, 0, 0, AL, 0, 0);
			dtrsm_rltn(ne, nv, 1.0, Lv, 0, 0, AL, 0, 0, AL, 0, 0);

			dtrsv_lnn(nv, Lv, 0, 0, lv, 0, lv, 0);

			dgese(ne, ne, 0.0, Le, 0, 0);
			dsyrk_ln(ne, nv, 1.0, AL, 0, 0, AL, 0, 0, 1.0, Le, 0, 0, Le, 0, 0);

			ddiaex(ne, 1.0, Le, 0, 0, se, 0);
			for(ii=0; ii<ne; ii++)
				{
				tmp = sqrt(se->pa[ii]);
//				tmp = sqrt(tmp);
//				tmp = sqrt(se->pa[ii]+tmp);
//				tmp = 1.0;
				se->pa[ii] = tmp==0 ? 1.0 : 1.0/tmp;
				}

			dgemm_dn(ne, ne, 1.0, se, 0, Le, 0, 0, 0.0, Le, 0, 0, Le, 0, 0);
			dgemm_nd(ne, ne, 1.0, Le, 0, 0, se, 0, 0.0, Le, 0, 0, Le, 0, 0);
			ddiare(ne, arg->reg_prim, Le, 0, 0);
			dpotrf_l(ne, Le, 0, 0, Le, 0, 0);

			dgemv_n(ne, nv, 1.0, AL, 0, 0, lv, 0, 1.0, res_b, 0, dpi, 0);

			dgemv_d(ne, 1.0, se, 0, dpi, 0, 0.0, dpi, 0, dpi, 0);
			dtrsv_lnn(ne, Le, 0, 0, dpi, 0, dpi, 0);
			dtrsv_ltn(ne, Le, 0, 0, dpi, 0, dpi, 0);
			dgemv_d(ne, 1.0, se, 0, dpi, 0, 0.0, dpi, 0, dpi, 0);

			dgemv_t(ne, nv, 1.0, A, 0, 0, dpi, 0, 0.0, lv, 0, lv, 0);
			dgemv_d(nv, 1.0, sv, 0, lv, 0, -1.0, dv, 0, dv, 0);

			dtrsv_lnn(nv, Lv, 0, 0, dv, 0, dv, 0);
			dtrsv_ltn(nv, Lv, 0, 0, dv, 0, dv, 0);
			dgemv_d(nv, 1.0, sv, 0, dv, 0, 0.0, dv, 0, dv, 0);

			}
		else // no scale
			{

//			dtrcp_l(nv, Hg, 0, 0, Lv, 0, 0);
			dgecp(nv, nv, Hg, 0, 0, Lv, 0, 0);

			dveccp(nv, res_g, 0, lv, 0);

			if(ns>0)
				{
				COND_SLACKS_FACT_solve(qp, qp_sol, arg, ws);
				}
			else if(nb+ng>0)
				{
				daxpy(nb+ng,  1.0, Gamma, nb+ng, Gamma, 0, tmp_nbg+0, 0);
				daxpy(nb+ng, -1.0, gamma, nb+ng, gamma, 0, tmp_nbg+1, 0);
				}
			if(nb>0)
				{
				ddiaad_sp(nb, 1.0, tmp_nbg+0, 0, idxb, Lv, 0, 0);
				dvecad_sp(nb, 1.0, tmp_nbg+1, 0, idxb, lv, 0);
				}
			if(ng>0)
				{
				dgemv_n(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+1, nb, 1.0, lv, 0, lv, 0);
				dgemm_nd(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+0, nb, 0.0, Ctx, 0, 0, Ctx, 0, 0);
				dsyrk_ln(nv, ng, 1.0, Ctx, 0, 0, Ct, 0, 0, 1.0, Lv, 0, 0, Lv, 0, 0);
				}

			dtrtr_l(nv, Lv, 0, 0, Lv, 0, 0);
			dgetrF(nv, nv, Lv, 0, 0, Lv, 0, 0, ipiv_v);

			dveccp(nv, lv, 0, dv, 0);
			dvecpe(nv, ipiv_v, lv, 0);
			dtrsv_lnu(nv, Lv, 0, 0, lv, 0, lv, 0);

			dgecp(ne, nv, A, 0, 0, AL+1, 0, 0);
			dcolpe(nv, ipiv_v, AL+1);
			dtrsm_rltu(ne, nv, 1.0, Lv, 0, 0, AL+1, 0, 0, AL+1, 0, 0);
//			dtrsm_runn(ne, nv, 1.0, Lv, 0, 0, A, 0, 0, AL+0, 0, 0);
			dtrtr_u(nv, Lv, 0, 0, Lv+1, 0, 0);
			dtrsm_rltn(ne, nv, 1.0, Lv+1, 0, 0, A, 0, 0, AL+0, 0, 0);

			dgemv_n(ne, nv, 1.0, AL+0, 0, 0, lv, 0, 1.0, res_b, 0, dpi, 0);

			dsyrk_ln(ne, nv, 1.0, AL+0, 0, 0, AL+1, 0, 0, 0.0, Le, 0, 0, Le, 0, 0);

#if 0
			dpotrf_l(ne, Le, 0, 0, Le, 0, 0);

			dtrsv_lnn(ne, Le, 0, 0, dpi, 0, dpi, 0);
			dtrsv_ltn(ne, Le, 0, 0, dpi, 0, dpi, 0);
#else
			dtrtr_l(ne, Le, 0, 0, Le, 0, 0);
			dgetrF(ne, ne, Le, 0, 0, Le, 0, 0, ipiv_e);

			dvecpe(ne, ipiv_e, dpi, 0);
			dtrsv_lnu(ne, Le, 0, 0, dpi, 0, dpi, 0);
			dtrsv_unn(ne, Le, 0, 0, dpi, 0, dpi, 0);
#endif

			dgemv_t(ne, nv, 1.0, A, 0, 0, dpi, 0, -1.0, dv, 0, dv, 0);

			dvecpe(nv, ipiv_v, dv, 0);
			dtrsv_lnu(nv, Lv, 0, 0, dv, 0, dv, 0);
			dtrsv_unn(nv, Lv, 0, 0, dv, 0, dv, 0);

			} // scale

		}
	else // ne==0
		{

		if(arg->scale)
			{

			dtrcp_l(nv, Hg, 0, 0, Lv, 0, 0);
			dveccp(nv, res_g, 0, lv, 0);

			if(ns>0)
				{
				COND_SLACKS_FACT_solve(qp, qp_sol, arg, ws);
				}
			else if(nb+ng>0)
				{
				daxpy(nb+ng,  1.0, Gamma, nb+ng, Gamma, 0, tmp_nbg+0, 0);
				daxpy(nb+ng, -1.0, gamma, nb+ng, gamma, 0, tmp_nbg+1, 0);
				}
			if(nb>0)
				{
				ddiaad_sp(nb, 1.0, tmp_nbg+0, 0, idxb, Lv, 0, 0);
				dvecad_sp(nb, 1.0, tmp_nbg+1, 0, idxb, lv, 0);
				}
			if(ng>0)
				{
				dgemm_nd(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+0, nb, 0.0, Ctx, 0, 0, Ctx, 0, 0);
				dgemv_n(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+1, nb, 1.0, lv, 0, lv, 0);
				dsyrk_ln(nv, ng, 1.0, Ctx, 0, 0, Ct, 0, 0, 1.0, Lv, 0, 0, Lv, 0, 0);
				}

			ddiaex(nv, 1.0, Lv, 0, 0, sv, 0);
			for(ii=0; ii<nv; ii++)
				{
				tmp = sqrt(sv->pa[ii]);
//				tmp = sqrt(tmp);
//				tmp = sqrt(sv->pa[ii]+tmp);
//				tmp = 1.0;
				sv->pa[ii] = tmp==0 ? 1.0 : 1.0/tmp;
				}

			dgemm_dn(nv, nv, 1.0, sv, 0, Lv, 0, 0, 0.0, Lv, 0, 0, Lv, 0, 0);
			dgemm_nd(nv, nv, 1.0, Lv, 0, 0, sv, 0, 0.0, Lv, 0, 0, Lv, 0, 0);
			ddiare(nv, arg->reg_prim, Lv, 0, 0);
			dpotrf_l_mn(nv, nv, Lv, 0, 0, Lv, 0, 0);

			dveccp(nv, lv, 0, dv, 0);
			dvecsc(nv, -1.0, dv, 0);

			dgemv_d(nv, 1.0, sv, 0, dv, 0, 0.0, dv, 0, dv, 0);
			dtrsv_lnn(nv, Lv, 0, 0, dv, 0, dv, 0);
			dtrsv_ltn(nv, Lv, 0, 0, dv, 0, dv, 0);
			dgemv_d(nv, 1.0, sv, 0, dv, 0, 0.0, dv, 0, dv, 0);
			}
		else // no scale
			{
	//		dtrcp_l(nv, Hg, 0, 0, Lv, 0, 0);
			dgecp(nv, nv, Hg, 0, 0, Lv, 0, 0);
			drowin(nv, 1.0, res_g, 0, Lv, nv, 0);

			if(ns>0)
				{
				COND_SLACKS_FACT_solve(qp, qp_sol, arg, ws);
				}
			else if(nb+ng>0)
				{
				daxpy(nb+ng,  1.0, Gamma, nb+ng, Gamma, 0, tmp_nbg+0, 0);
				daxpy(nb+ng, -1.0, gamma, nb+ng, gamma, 0, tmp_nbg+1, 0);
				}
			if(nb>0)
				{
				ddiaad_sp(nb, 1.0, tmp_nbg+0, 0, idxb, Lv, 0, 0);
				drowad_sp(nb, 1.0, tmp_nbg+1, 0, idxb, Lv, nv, 0);
				}
			if(ng>0)
				{
				dgemm_nd(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+0, nb, 0.0, Ctx, 0, 0, Ctx, 0, 0);
				drowin(ng, 1.0, tmp_nbg+1, nb, Ctx, nv, 0);
				dsyrk_ln_mn(nv+1, nv, ng, 1.0, Ctx, 0, 0, Ct, 0, 0, 1.0, Lv, 0, 0, Lv, 0, 0);
				}

			drowex(nv, -1.0, Lv, nv, 0, dv, 0);

			dtrtr_l(nv, Lv, 0, 0, Lv, 0, 0);
			dgetrF(nv, nv, Lv, 0, 0, Lv, 0, 0, ipiv_v);

			dvecpe(nv, ipiv_v, dv, 0);
			dtrsv_lnu(nv, Lv, 0, 0, dv, 0, dv, 0);
			dtrsv_unn(nv, Lv, 0, 0, dv, 0, dv, 0);

			} // scale

		} // ne>0

	if(nb+ng>0)
		{
		if(nb>0)
			dvecex_sp(nb, 1.0, idxb, dv, 0, dt, 0);

		if(ng>0)
			dgemv_t(nv, ng, 1.0, Ct, 0, 0, dv, 0, 0.0, dt, nb, dt, nb);

		dveccp(nb+ng, dt, 0, dt, nb+ng);
		dvecsc(nb+ng, -1.0, dt, nb+ng);

		if(ns>0)
			EXPAND_SLACKS(qp, qp_sol, ws);

		d_compute_lam_t_qp(qp->d->pa, qp->m->pa, qp_sol->lam->pa, qp_sol->t->pa, cws);
		}



	}
#endif


// range-space (Schur complement) method
void d_solve_kkt_step_dense_qp(struct d_dense_qp* qp, struct d_dense_qp_sol* qp_sol, struct d_dense_qp_ipm_arg* arg, struct d_dense_qp_ipm_ws* ws) {

    int nv = qp->dim->nv;
    int ne = qp->dim->ne;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;
    int ns = qp->dim->ns;

    struct mat* A = qp->A;
    struct mat* Ct = qp->Ct;
    int* idxb = qp->idxb;
    //	struct vec *res_g = ws->res->res_g;
    //	struct vec *res_b = ws->res->res_b;
    struct vec* res_g = qp->gz;
    struct vec* res_b = qp->b;

    struct vec* dv = qp_sol->v;
    struct vec* dpi = qp_sol->pi;
    struct vec* dt = qp_sol->t;

    struct mat* Lv = ws->Lv;
    struct mat* Le = ws->Le;
    struct mat* Ctx = ws->Ctx;
    struct mat* AL = ws->AL;
    struct vec* lv = ws->lv;
    struct vec* sv = ws->sv;
    struct vec* se = ws->se;
    struct vec* gamma = ws->gamma;
    struct vec* tmp_nbg = ws->tmp_nbg;

    // null space
    struct mat* A_LQ = ws->A_LQ;
    struct mat* A_Q = ws->A_Q;
    struct mat* Zt = ws->Zt;
    struct mat* ZtH = ws->ZtH;
    struct mat* ZtHZ = ws->ZtHZ;
    struct vec* xy = ws->xy;
    struct vec* Yxy = ws->Yxy;
    struct vec* xz = ws->xz;
    struct vec* tmp_nv = ws->tmp_nv;
    void* lq_work = ws->lq_work_null;

    struct d_core_qp_ipm_workspace* cws = ws->core_workspace;

    if (nb > 0 | ng > 0) {
        d_compute_gamma_qp(qp->d->pa, qp->m->pa, cws);
    }

    dveccp(nv, res_g, 0, lv, 0);

    if (ns > 0) {
        COND_SLACKS_solve(qp, qp_sol, ws);
    } else if (nb + ng > 0) {
        daxpy(nb + ng, -1.0, gamma, nb + ng, gamma, 0, tmp_nbg + 1, 0);
    }
    if (nb > 0) {
        dvecad_sp(nb, 1.0, tmp_nbg + 1, 0, idxb, lv, 0);
    }
    if (ng > 0) {
        dgemv_n(nv, ng, 1.0, Ct, 0, 0, tmp_nbg + 1, nb, 1.0, lv, 0, lv, 0);
    }

    if (ne > 0) {

        if (arg->kkt_fact_alg == 0)  // null space method
        {

            dtrsv_lnn(ne, A_LQ, 0, 0, res_b, 0, xy, 0);
            dgemv_t(ne, nv, 1.0, A_Q, 0, 0, xy, 0, 0.0, Yxy, 0, Yxy, 0);

            dgemv_n(nv - ne, nv, -1.0, ZtH, 0, 0, Yxy, 0, 0.0, xz, 0, xz, 0);
            dgemv_n(nv - ne, nv, -1.0, Zt, 0, 0, lv, 0, 1.0, xz, 0, xz, 0);
            dtrsv_lnn(nv - ne, ZtHZ, 0, 0, xz, 0, xz, 0);
            dtrsv_ltn(nv - ne, ZtHZ, 0, 0, xz, 0, xz, 0);

            dgemv_t(nv - ne, nv, 1.0, Zt, 0, 0, xz, 0, 1.0, Yxy, 0, dv, 0);

            dsymv_l(nv, 1.0, Lv, 0, 0, dv, 0, 1.0, lv, 0, tmp_nv, 0);
            dgemv_n(ne, nv, 1.0, A_Q, 0, 0, tmp_nv, 0, 0.0, dpi, 0, dpi, 0);
            dtrsv_ltn(ne, A_LQ, 0, 0, dpi, 0, dpi, 0);

        } else  // schur-complement method
        {

            if (ws->scale) {

                dgemv_d(nv, 1.0, sv, 0, lv, 0, 0.0, lv, 0, lv, 0);
                dveccp(nv, lv, 0, dv, 0);

                dtrsv_lnn(nv, Lv, 0, 0, lv, 0, lv, 0);

                dgemv_n(ne, nv, 1.0, AL, 0, 0, lv, 0, 1.0, res_b, 0, dpi, 0);

                dgemv_d(ne, 1.0, se, 0, dpi, 0, 0.0, dpi, 0, dpi, 0);
                dtrsv_lnn(ne, Le, 0, 0, dpi, 0, dpi, 0);
                dtrsv_ltn(ne, Le, 0, 0, dpi, 0, dpi, 0);
                dgemv_d(ne, 1.0, se, 0, dpi, 0, 0.0, dpi, 0, dpi, 0);

                dgemv_t(ne, nv, 1.0, A, 0, 0, dpi, 0, 0.0, lv, 0, lv, 0);
                dgemv_d(nv, 1.0, sv, 0, lv, 0, -1.0, dv, 0, dv, 0);

                dtrsv_lnn(nv, Lv, 0, 0, dv, 0, dv, 0);
                dtrsv_ltn(nv, Lv, 0, 0, dv, 0, dv, 0);
                dgemv_d(nv, 1.0, sv, 0, dv, 0, 0.0, dv, 0, dv, 0);

            } else  // no scale
            {

                dveccp(nv, lv, 0, dv, 0);

                dtrsv_lnn(nv, Lv, 0, 0, lv, 0, lv, 0);

                dgemv_n(ne, nv, 1.0, AL, 0, 0, lv, 0, 1.0, res_b, 0, dpi, 0);

                dtrsv_lnn(ne, Le, 0, 0, dpi, 0, dpi, 0);
                dtrsv_ltn(ne, Le, 0, 0, dpi, 0, dpi, 0);

                dgemv_t(ne, nv, 1.0, A, 0, 0, dpi, 0, -1.0, dv, 0, dv, 0);

                dtrsv_lnn(nv, Lv, 0, 0, dv, 0, dv, 0);
                dtrsv_ltn(nv, Lv, 0, 0, dv, 0, dv, 0);

            }  // scale

        }  // schur-complemen method

    } else  // ne==0
    {

        if (ws->scale) {

            dveccp(nv, lv, 0, dv, 0);
            dvecsc(nv, -1.0, dv, 0);

            dgemv_d(nv, 1.0, sv, 0, dv, 0, 0.0, dv, 0, dv, 0);
            dtrsv_lnn(nv, Lv, 0, 0, dv, 0, dv, 0);
            dtrsv_ltn(nv, Lv, 0, 0, dv, 0, dv, 0);
            dgemv_d(nv, 1.0, sv, 0, dv, 0, 0.0, dv, 0, dv, 0);

        } else  // no scale
        {

            dveccp(nv, lv, 0, dv, 0);
            dvecsc(nv, -1.0, dv, 0);

            dtrsv_lnn(nv, Lv, 0, 0, dv, 0, dv, 0);
            dtrsv_ltn(nv, Lv, 0, 0, dv, 0, dv, 0);

        }  // scale

    }  // ne>0

    if (nb + ng > 0) {
        if (nb > 0)
            dvecex_sp(nb, 1.0, idxb, dv, 0, dt, 0);

        if (ng > 0)
            dgemv_t(nv, ng, 1.0, Ct, 0, 0, dv, 0, 0.0, dt, nb, dt, nb);

        dveccp(nb + ng, dt, 0, dt, nb + ng);
        dvecsc(nb + ng, -1.0, dt, nb + ng);

        if (ns > 0)
            EXPAND_SLACKS(qp, qp_sol, ws);

        d_compute_lam_t_qp(qp->d->pa, qp->m->pa, qp_sol->lam->pa, qp_sol->t->pa, cws);
    }
}


void d_dense_qp_remove_lin_dep_eq(struct d_dense_qp* qp, struct d_dense_qp_ipm_arg* arg, struct d_dense_qp_ipm_ws* ws) {

    int ii, jj, ll;
    int stop_jj, jj0;

    int nv = qp->dim->nv;
    int ne = qp->dim->ne;

    struct mat* A = qp->A;
    struct vec* b = qp->b;

    struct mat* A_li = ws->A_li;
    struct vec* b_li = ws->b_li;
    struct mat* Ab_LU = ws->Ab_LU;
    void* lq_work_null = ws->lq_work_null;
    int* ipiv_v = ws->ipiv_v;
    int* ipiv_e = ws->ipiv_e;
    int* ipiv_e1 = ws->ipiv_e1;

    int ne_li = 0;

    // TODO tuning, single precision
    double thr = 1e-14;

    double tmp_diag, pivot, tmp, tmp_b, tmp_max;
    int idx0_max, idx1_max, tmp_int;

    ws->status = SUCCESS;

    if (ne > 0) {
        // augment A with b
        dgecp(ne, nv, A, 0, 0, Ab_LU, 0, 0);
        dcolin(ne, b, 0, Ab_LU, 0, nv);

        // row-pivot LU factorization
        dgetrf_rp(ne, nv + 1, Ab_LU, 0, 0, Ab_LU, 0, 0, ipiv_e);

        // get pivot in absolute form
        for (ll = 0; ll < ne; ll++) {
            ipiv_e1[ll] = ll;
        }
        for (ll = 0; ll < ne; ll++) {
            tmp_int = ipiv_e1[ll];
            ipiv_e1[ll] = ipiv_e1[ipiv_e[ll]];
            ipiv_e1[ipiv_e[ll]] = tmp_int;
        }

        jj0 = 0;
        for (ii = 0; ii < ne; ii++) {
            pivot = MATEL(Ab_LU, ii, ii);
            if (fabs(pivot) <= thr) {
                jj = ii + 1 > jj0 ? ii + 1 : jj0;
                stop_jj = 0;
                while (stop_jj == 0 & jj < nv) {
                    tmp_max = thr;
                    idx0_max = -1;
                    idx1_max = -1;
                    for (ll = ii; ll < ne & ll <= jj; ll++) {
                        tmp = fabs(MATEL(Ab_LU, ll, jj));
                        if (tmp > tmp_max) {
                            tmp_max = tmp;
                            idx0_max = ll;
                            idx1_max = jj;
                            jj0 = jj;
                            stop_jj = 1;
                        }
                    }
                    jj++;
                }
                if (stop_jj == 1) {
                    // swap rows
                    if (tmp_max > thr & idx0_max != ii) {
                        drowsw(nv + 1 - idx1_max, Ab_LU, ii, idx1_max, Ab_LU, idx0_max, idx1_max);
                        tmp_int = ipiv_e1[ii];
                        ipiv_e1[ii] = ipiv_e1[idx0_max];
                        ipiv_e1[idx0_max] = tmp_int;
                    }
                    // copy li eq
                    dgecp(1, nv, A, ipiv_e1[ii], 0, A_li, ne_li, 0);
                    dveccp(1, b, ipiv_e1[ii], b_li, ne_li);
                    ne_li++;
                    // pivot
                    pivot = MATEL(Ab_LU, ii, idx1_max);
                    // clear below TODO implement using level 2 BLAS !!!
                    for (ll = ii + 1; ll < ne & ll <= idx1_max; ll++) {
                        tmp = fabs(MATEL(Ab_LU, ll, idx1_max));
                        if (tmp != 0.0) {
                            tmp = -tmp / pivot;
                            dgead(1, nv + 1 - idx1_max, tmp, Ab_LU, ii, idx1_max, Ab_LU, ll, idx1_max);
                        }
                    }
                } else {
                    // all remaining matrix is zero: check b and return
                    for (ll = ii; ll < ne; ll++) {
                        tmp = fabs(MATEL(Ab_LU, ll, nv));
                        if (tmp > thr) {
                            ws->status = INCONS_EQ;
                        }
                    }
                    goto swap_A_b;
                }
            } else {
                // copy li eq
                dgecp(1, nv, A, ipiv_e1[ii], 0, A_li, ne_li, 0);
                dveccp(1, b, ipiv_e1[ii], b_li, ne_li);
                ne_li++;
            }
        }

    swap_A_b:
        if (ne_li < ne) {
            //			printf("\nne %d, ne_li %d\n", ne, ne_li);
            ws->ne_bkp = qp->dim->ne;
            qp->dim->ne = ne_li;
            ws->A_bkp = qp->A;
            qp->A = A_li;
            ws->b_bkp = qp->b;
            qp->b = b_li;
        }
    }

    // printf("\nne %d ne_li %d\n", ne, ne_li);
    // printf("\nA_li\n");
    // blasfeo_print_dmat(ne_li, nv, A_li, 0, 0);
    // printf("\nb_li\n");
    // blasfeo_print_tran_dvec(ne_li, b_li, 0);
}


void d_dense_qp_restore_lin_dep_eq(struct d_dense_qp* qp, struct d_dense_qp_ipm_arg* arg, struct d_dense_qp_ipm_ws* ws) {

    int ii, jj;

    int nv = qp->dim->nv;
    int ne = qp->dim->ne;

    struct mat* A = qp->A;
    struct vec* b = qp->b;

    struct mat* A_li = ws->A_li;
    struct vec* b_li = ws->b_li;
    void* lq_work_null = ws->lq_work_null;
    int* ipiv_v = ws->ipiv_v;

    if (ne > 0) {
        if (ne < ws->ne_bkp) {
            qp->dim->ne = ws->ne_bkp;
            qp->A = ws->A_bkp;
            qp->b = ws->b_bkp;
        }
    }

    // printf("\nne %d\n", ne);
    // printf("\nA\n");
    // blasfeo_print_dmat(ne, nv, A, 0, 0);
    // printf("\nb\n");
    // blasfeo_print_tran_dvec(ne, b, 0);
}


void d_dense_qp_compute_obj(struct d_dense_qp* qp, struct d_dense_qp_sol* qp_sol, struct d_dense_qp_ipm_arg* arg, struct d_dense_qp_ipm_ws* ws) {

    int nv = qp->dim->nv;
    int ns = qp->dim->ns;

    struct mat* Hg = qp->Hv;
    struct vec* Z = qp->Z;
    struct vec* gz = qp->gz;

    struct vec* v = qp_sol->v;

    // TODO soft constraints !!!!!!!!!!!!!!!!!!!!!!!

    struct vec* tmp_nv = ws->tmp_nv;
    struct vec* tmp_2ns = ws->tmp_2ns;

    dsymv_l(nv, 0.5, Hg, 0, 0, v, 0, 1.0, gz, 0, tmp_nv, 0);
    qp_sol->obj = ddot(nv, tmp_nv, 0, v, 0);

    dgemv_d(2 * ns, 0.5, Z, 0, v, nv, 1.0, gz, nv, tmp_2ns, 0);
    qp_sol->obj += ddot(2 * ns, tmp_2ns, 0, v, nv);

    qp_sol->valid_obj = 1;
}
