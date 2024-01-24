#include <math.h>
#include <stdlib.h>

#include "hpipm/blas.h"
#include "hpipm/common.h"
#include "hpipm/ipm_core/d_core_qp_ipm.h"
#include "hpipm/ipm_core/d_core_qp_ipm_aux.h"
#include "hpipm/ocp/d_ocp_qp.h"
#include "hpipm/ocp/d_ocp_qp_dim.h"
#include "hpipm/ocp/d_ocp_qp_ipm.h"
#include "hpipm/ocp/d_ocp_qp_sol.h"


// backward Riccati recursion
void d_ocp_qp_fact_solve_kkt_unconstr(struct d_ocp_qp* qp, struct d_ocp_qp_sol* qp_sol, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws) {

    int ii;

    int N = qp->dim->N;
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;

    struct mat* BAbt = qp->BAbt;
    struct mat* RSQrq = qp->RSQrq;
    struct vec* b = qp->b;
    struct vec* rqz = qp->rqz;

    struct vec* ux = qp_sol->ux;
    struct vec* pi = qp_sol->pi;

    struct mat* L = ws->L;
    struct vec* l = ws->l;
    struct mat* AL = ws->AL;

    struct vec* tmp_nuxM = ws->tmp_nuxM;

    if (ws->square_root_alg) {
        ws->valid_ric_p = 0;

        // factorization and backward substitution

        // last stage
        drowin(nu[N] + nx[N], 1.0, rqz + N, 0, RSQrq + N, nu[N] + nx[N], 0);
        ddiare(nu[N] + nx[N], arg->reg_prim, RSQrq + N, 0, 0);
        dpotrf_l_mn(nu[N] + nx[N] + 1, nu[N] + nx[N], RSQrq + N, 0, 0, L + N, 0, 0);
        ddiare(nu[N] + nx[N], -arg->reg_prim, RSQrq + N, 0, 0);

        // middle stages
        for (ii = 0; ii < N; ii++) {
            drowin(nx[N - ii], 1.0, b + N - ii - 1, 0, BAbt + N - ii - 1, nu[N - ii - 1] + nx[N - ii - 1], 0);
            dtrmm_rlnn(nu[N - ii - 1] + nx[N - ii - 1] + 1, nx[N - ii], 1.0, L + (N - ii), nu[N - ii], nu[N - ii], BAbt + N - ii - 1, 0, 0, AL, 0, 0);
            dgead(1, nx[N - ii], 1.0, L + N - ii, nu[N - ii] + nx[N - ii], nu[N - ii], AL, nu[N - ii - 1] + nx[N - ii - 1], 0);

            drowin(nu[N - ii - 1] + nx[N - ii - 1], 1.0, rqz + N - ii - 1, 0, RSQrq + N - ii - 1, nu[N - ii - 1] + nx[N - ii - 1], 0);
            ddiare(nu[N - ii - 1] + nx[N - ii - 1], arg->reg_prim, RSQrq + N - ii - 1, 0, 0);
            dsyrk_dpotrf_ln_mn(nu[N - ii - 1] + nx[N - ii - 1] + 1, nu[N - ii - 1] + nx[N - ii - 1], nx[N - ii], AL, 0, 0, AL, 0, 0, RSQrq + (N - ii - 1), 0, 0, L + (N - ii - 1), 0, 0);
            ddiare(nu[N - ii - 1] + nx[N - ii - 1], -arg->reg_prim, RSQrq + N - ii - 1, 0, 0);
        }

        // forward substitution

        // first stage
        ii = 0;
        drowex(nu[ii] + nx[ii], -1.0, L + ii, nu[ii] + nx[ii], 0, ux + ii, 0);
        dveccpsc(nu[ii] + nx[ii], -1.0, ux + ii, 0, l + ii, 0);
        drowex(nx[ii + 1], 1.0, L + ii + 1, nu[ii + 1] + nx[ii + 1], nu[ii + 1], tmp_nuxM, 0);
        dveccp(nx[ii + 1], tmp_nuxM, 0, l + ii + 1, nu[ii + 1]);  // TODO remove tmp_nuxM and use l instead !!!!!
        dtrsv_ltn(nu[ii] + nx[ii], L + ii, 0, 0, ux + ii, 0, ux + ii, 0);
        dgemv_t(nu[ii] + nx[ii], nx[ii + 1], 1.0, BAbt + ii, 0, 0, ux + ii, 0, 1.0, b + ii, 0, ux + ii + 1, nu[ii + 1]);
        if (arg->comp_dual_sol_eq) {
            dtrmv_ltn(nx[ii + 1], L + (ii + 1), nu[ii + 1], nu[ii + 1], ux + (ii + 1), nu[ii + 1], pi + ii, 0);
            daxpy(nx[ii + 1], 1.0, tmp_nuxM, 0, pi + ii, 0, pi + ii, 0);
            dtrmv_lnn(nx[ii + 1], L + (ii + 1), nu[ii + 1], nu[ii + 1], pi + ii, 0, pi + ii, 0);
        }

        // middle stages
        for (ii = 1; ii < N; ii++) {
            drowex(nu[ii], -1.0, L + ii, nu[ii] + nx[ii], 0, ux + ii, 0);
            dveccpsc(nu[ii], -1.0, ux + ii, 0, l + ii, 0);
            drowex(nx[ii + 1], 1.0, L + (ii + 1), nu[ii + 1] + nx[ii + 1], nu[ii + 1], tmp_nuxM, 0);
            dveccp(nx[ii + 1], tmp_nuxM, 0, l + ii + 1, nu[ii + 1]);  // TODO remove tmp_nuxM and use l instead !!!!!
            dtrsv_ltn_mn(nu[ii] + nx[ii], nu[ii], L + ii, 0, 0, ux + ii, 0, ux + ii, 0);
            dgemv_t(nu[ii] + nx[ii], nx[ii + 1], 1.0, BAbt + ii, 0, 0, ux + ii, 0, 1.0, b + ii, 0, ux + (ii + 1), nu[ii + 1]);
            if (arg->comp_dual_sol_eq) {
                dtrmv_ltn(nx[ii + 1], L + (ii + 1), nu[ii + 1], nu[ii + 1], ux + (ii + 1), nu[ii + 1], pi + ii, 0);
                daxpy(nx[ii + 1], 1.0, tmp_nuxM, 0, pi + ii, 0, pi + ii, 0);
                dtrmv_lnn(nx[ii + 1], L + (ii + 1), nu[ii + 1], nu[ii + 1], pi + ii, 0, pi + ii, 0);
            }
        }

        // last stage
        ii = N;
        drowex(nu[ii], -1.0, L + ii, nu[ii] + nx[ii], 0, ux + ii, 0);
        dveccpsc(nu[ii], -1.0, ux + ii, 0, l + ii, 0);
        dtrsv_ltn_mn(nu[ii] + nx[ii], nu[ii], L + ii, 0, 0, ux + ii, 0, ux + ii, 0);

    } else {
        ws->valid_ric_p = 1;

        struct mat* P = ws->P;
        struct mat* Ls = ws->Ls;

        // factorization and backward substitution

        // last stage
        drowin(nu[N] + nx[N], 1.0, rqz + N, 0, RSQrq + N, nu[N] + nx[N], 0);
        ddiare(nu[N], arg->reg_prim, RSQrq + N, 0, 0);
        dpotrf_l_mn(nu[N] + nx[N] + 1, nu[N], RSQrq + N, 0, 0, L + N, 0, 0);
        ddiare(nu[N], -arg->reg_prim, RSQrq + N, 0, 0);
        dgecp(nx[N] + 1, nu[N], L + N, nu[N], 0, Ls, 0, 0);
        dsyrk_ln_mn(nx[N] + 1, nx[N], nu[N], -1.0, Ls, 0, 0, Ls, 0, 0, 1.0, RSQrq + N, nu[N], nu[N], P + N, 0, 0);
        dtrtr_l(nx[N], P + N, 0, 0, P + N, 0, 0);

        // middle statges
        for (ii = 0; ii < N - 1; ii++) {
            drowin(nx[N - ii], 1.0, b + N - ii - 1, 0, BAbt + N - ii - 1, nu[N - ii - 1] + nx[N - ii - 1], 0);
            dgemm_nt(nu[N - ii - 1] + nx[N - ii - 1] + 1, nx[N - ii], nx[N - ii], 1.0, BAbt + N - ii - 1, 0, 0, P + N - ii, 0, 0, 0.0, AL, 0, 0, AL, 0, 0);  // TODO symm
            dgead(1, nx[N - ii], 1.0, P + N - ii, nx[N - ii], 0, AL, nu[N - ii - 1] + nx[N - ii - 1], 0);
            drowin(nu[N - ii - 1] + nx[N - ii - 1], 1.0, rqz + N - ii - 1, 0, RSQrq + N - ii - 1, nu[N - ii - 1] + nx[N - ii - 1], 0);
            ddiare(nu[N - ii - 1], arg->reg_prim, RSQrq + N - ii - 1, 0, 0);
            dsyrk_ln_mn(nu[N - ii - 1] + nx[N - ii - 1] + 1, nu[N - ii - 1] + nx[N - ii - 1], nx[N - ii], 1.0, AL, 0, 0, BAbt + N - ii - 1, 0, 0, 1.0, RSQrq + N - ii - 1, 0, 0, L + N - ii - 1, 0, 0);
            ddiare(nu[N - ii - 1], -arg->reg_prim, RSQrq + N - ii - 1, 0, 0);
            dpotrf_l_mn(nu[N - ii - 1] + nx[N - ii - 1] + 1, nu[N - ii - 1], L + N - ii - 1, 0, 0, L + N - ii - 1, 0, 0);
            dgecp(nx[N - ii - 1] + 1, nu[N - ii - 1], L + N - ii - 1, nu[N - ii - 1], 0, Ls, 0, 0);
            dsyrk_ln_mn(nx[N - ii - 1] + 1, nx[N - ii - 1], nu[N - ii - 1], -1.0, Ls, 0, 0, Ls, 0, 0, 1.0, L + N - ii - 1, nu[N - ii - 1], nu[N - ii - 1], P + N - ii - 1, 0, 0);
            dtrtr_l(nx[N - ii - 1], P + N - ii - 1, 0, 0, P + N - ii - 1, 0, 0);
        }

        // first stage: factorize P in L too
        if (N > 0) {
            drowin(nx[N - ii], 1.0, b + N - ii - 1, 0, BAbt + N - ii - 1, nu[N - ii - 1] + nx[N - ii - 1], 0);
            dgemm_nt(nu[N - ii - 1] + nx[N - ii - 1] + 1, nx[N - ii], nx[N - ii], 1.0, BAbt + N - ii - 1, 0, 0, P + N - ii, 0, 0, 0.0, AL, 0, 0, AL, 0, 0);  // TODO symm
            dgead(1, nx[N - ii], 1.0, P + N - ii, nx[N - ii], 0, AL, nu[N - ii - 1] + nx[N - ii - 1], 0);
            drowin(nu[N - ii - 1] + nx[N - ii - 1], 1.0, rqz + N - ii - 1, 0, RSQrq + N - ii - 1, nu[N - ii - 1] + nx[N - ii - 1], 0);
            ddiare(nu[N - ii - 1], arg->reg_prim, RSQrq + N - ii - 1, 0, 0);
            dsyrk_dpotrf_ln_mn(nu[N - ii - 1] + nx[N - ii - 1] + 1, nu[N - ii - 1] + nx[N - ii - 1], nx[N - ii], AL, 0, 0, BAbt + N - ii - 1, 0, 0, RSQrq + N - ii - 1, 0, 0, L + N - ii - 1, 0, 0);
            ddiare(nu[N - ii - 1], -arg->reg_prim, RSQrq + N - ii - 1, 0, 0);
        }

        // forward substitution

        // first stage
        ii = 0;
        drowex(nu[ii] + nx[ii], -1.0, L + ii, nu[ii] + nx[ii], 0, ux + ii, 0);
        dveccpsc(nu[ii] + nx[ii], -1.0, ux + ii, 0, l + ii, 0);
        drowex(nx[ii + 1], 1.0, P + ii + 1, nx[ii + 1], 0, tmp_nuxM, 0);
        dveccp(nx[ii + 1], tmp_nuxM, 0, l + ii + 1, nu[ii + 1]);  // TODO remove tmp_nuxM and use l instead !!!!!
        dtrsv_ltn(nu[ii] + nx[ii], L + ii, 0, 0, ux + ii, 0, ux + ii, 0);
        dgemv_t(nu[ii] + nx[ii], nx[ii + 1], 1.0, BAbt + ii, 0, 0, ux + ii, 0, 1.0, b + ii, 0, ux + ii + 1, nu[ii + 1]);
        if (arg->comp_dual_sol_eq) {
            dgemv_n(nx[ii + 1], nx[ii + 1], 1.0, P + ii + 1, 0, 0, ux + ii + 1, nu[ii + 1], 1.0, tmp_nuxM, 0, pi + ii, 0);
        }

        // middle stages
        for (ii = 1; ii < N; ii++) {
            drowex(nu[ii], -1.0, L + ii, nu[ii] + nx[ii], 0, ux + ii, 0);
            dveccpsc(nu[ii], -1.0, ux + ii, 0, l + ii, 0);
            drowex(nx[ii + 1], 1.0, P + ii + 1, nx[ii + 1], 0, tmp_nuxM, 0);
            dveccp(nx[ii + 1], tmp_nuxM, 0, l + ii + 1, nu[ii + 1]);  // TODO remove tmp_nuxM and use l instead !!!!!
            dtrsv_ltn_mn(nu[ii] + nx[ii], nu[ii], L + ii, 0, 0, ux + ii, 0, ux + ii, 0);
            dgemv_t(nu[ii] + nx[ii], nx[ii + 1], 1.0, BAbt + ii, 0, 0, ux + ii, 0, 1.0, b + ii, 0, ux + (ii + 1), nu[ii + 1]);
            if (arg->comp_dual_sol_eq) {
                dgemv_n(nx[ii + 1], nx[ii + 1], 1.0, P + ii + 1, 0, 0, ux + ii + 1, nu[ii + 1], 1.0, tmp_nuxM, 0, pi + ii, 0);
            }
        }

        // last stage
        ii = N;
        drowex(nu[ii], -1.0, L + ii, nu[ii] + nx[ii], 0, ux + ii, 0);
        dveccpsc(nu[ii], -1.0, ux + ii, 0, l + ii, 0);
        dtrsv_ltn_mn(nu[ii] + nx[ii], nu[ii], L + ii, 0, 0, ux + ii, 0, ux + ii, 0);
    }
}


static void COND_SLACKS_FACT_solve(int ss, struct d_ocp_qp* qp, struct d_ocp_qp_sol* qp_sol, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws) {

    int ii, idx;

    int nx0 = qp->dim->nx[ss];
    int nu0 = qp->dim->nu[ss];
    int nb0 = qp->dim->nb[ss];
    int ng0 = qp->dim->ng[ss];
    int ns0 = qp->dim->ns[ss];

    struct vec* Z = qp->Z + ss;
    int* idxs_rev0 = qp->idxs_rev[ss];

    //	struct vec *res_g = ws->res->res_g+ss; // TODO !!!
    struct vec* res_g = qp->rqz + ss;

    //	struct vec *dux = ws->sol_step->ux+ss; // TODO !!!
    struct vec* dux = qp_sol->ux + ss;

    struct vec* Gamma = ws->Gamma + ss;
    struct vec* gamma = ws->gamma + ss;
    struct vec* Zs_inv = ws->Zs_inv + ss;
    struct vec* tmp_nbgM = ws->tmp_nbgM;

    double* ptr_Gamma = Gamma->pa;
    double* ptr_gamma = gamma->pa;
    double* ptr_Z = Z->pa;
    double* ptr_Zs_inv = Zs_inv->pa;
    double* ptr_dux = dux->pa;
    double* ptr_res_g = res_g->pa;
    double* ptr_tmp0 = (tmp_nbgM + 0)->pa;
    double* ptr_tmp1 = (tmp_nbgM + 1)->pa;
    double* ptr_tmp2 = (tmp_nbgM + 2)->pa;
    double* ptr_tmp3 = (tmp_nbgM + 3)->pa;

    // idxs_rev
    for (ii = 0; ii < nb0 + ng0; ii++) {
        idx = idxs_rev0[ii];
        if (idx != -1) {
            // ii   constr index
            // idx <= slack index
            ptr_Zs_inv[0 + idx] = ptr_Z[0 + idx] + arg->reg_prim + ptr_Gamma[0 + ii] + ptr_Gamma[2 * nb0 + 2 * ng0 + idx];
            ptr_Zs_inv[ns0 + idx] = ptr_Z[ns0 + idx] + arg->reg_prim + ptr_Gamma[nb0 + ng0 + ii] + ptr_Gamma[2 * nb0 + 2 * ng0 + ns0 + idx];
            ptr_dux[nu0 + nx0 + idx] = ptr_res_g[nu0 + nx0 + idx] + ptr_gamma[0 + ii] + ptr_gamma[2 * nb0 + 2 * ng0 + idx];
            ptr_dux[nu0 + nx0 + ns0 + idx] = ptr_res_g[nu0 + nx0 + ns0 + idx] + ptr_gamma[nb0 + ng0 + ii] + ptr_gamma[2 * nb0 + 2 * ng0 + ns0 + idx];
            ptr_Zs_inv[0 + idx] = 1.0 / ptr_Zs_inv[0 + idx];
            ptr_Zs_inv[ns0 + idx] = 1.0 / ptr_Zs_inv[ns0 + idx];
            ptr_tmp0[ii] = ptr_Gamma[0 + ii] - ptr_Gamma[0 + ii] * ptr_Zs_inv[0 + idx] * ptr_Gamma[0 + ii];
            ptr_tmp1[ii] = ptr_Gamma[nb0 + ng0 + ii] - ptr_Gamma[nb0 + ng0 + ii] * ptr_Zs_inv[ns0 + idx] * ptr_Gamma[nb0 + ng0 + ii];
            ptr_tmp2[ii] = ptr_gamma[0 + ii] - ptr_Gamma[0 + ii] * ptr_Zs_inv[0 + idx] * ptr_dux[nu0 + nx0 + idx];
            ptr_tmp3[ii] = ptr_gamma[nb0 + ng0 + ii] - ptr_Gamma[nb0 + ng0 + ii] * ptr_Zs_inv[ns0 + idx] * ptr_dux[nu0 + nx0 + ns0 + idx];  // + ???
        } else {
            ptr_tmp0[ii] = ptr_Gamma[0 + ii];
            ptr_tmp1[ii] = ptr_Gamma[nb0 + ng0 + ii];
            ptr_tmp2[ii] = ptr_gamma[0 + ii];
            ptr_tmp3[ii] = ptr_gamma[nb0 + ng0 + ii];
        }
    }

    daxpy(nb0 + ng0, 1.0, tmp_nbgM + 1, 0, tmp_nbgM + 0, 0, tmp_nbgM + 0, 0);
    daxpy(nb0 + ng0, -1.0, tmp_nbgM + 3, 0, tmp_nbgM + 2, 0, tmp_nbgM + 1, 0);
}


static void COND_SLACKS_solve(int ss, struct d_ocp_qp* qp, struct d_ocp_qp_sol* qp_sol, struct d_ocp_qp_ipm_ws* ws) {

    int ii, idx;

    int nx0 = qp->dim->nx[ss];
    int nu0 = qp->dim->nu[ss];
    int nb0 = qp->dim->nb[ss];
    int ng0 = qp->dim->ng[ss];
    int ns0 = qp->dim->ns[ss];

    int* idxs_rev0 = qp->idxs_rev[ss];

    //	struct vec *res_g = ws->res->res_g+ss; // TODO !!!
    struct vec* res_g = qp->rqz + ss;

    //	struct vec *dux = ws->sol_step->ux+ss; // TODO !!!
    struct vec* dux = qp_sol->ux + ss;

    struct vec* Gamma = ws->Gamma + ss;
    struct vec* gamma = ws->gamma + ss;
    struct vec* Zs_inv = ws->Zs_inv + ss;
    struct vec* tmp_nbgM = ws->tmp_nbgM;

    double* ptr_Gamma = Gamma->pa;
    double* ptr_gamma = gamma->pa;
    double* ptr_Zs_inv = Zs_inv->pa;
    double* ptr_dux = dux->pa;
    double* ptr_res_g = res_g->pa;
    double* ptr_tmp2 = (tmp_nbgM + 2)->pa;
    double* ptr_tmp3 = (tmp_nbgM + 3)->pa;

    // idxs_rev
    for (ii = 0; ii < nb0 + ng0; ii++) {
        idx = idxs_rev0[ii];
        if (idx != -1) {
            // ii  <= constr index
            // idx <= slack index
            ptr_dux[nu0 + nx0 + idx] = ptr_res_g[nu0 + nx0 + idx] + ptr_gamma[0 + ii] + ptr_gamma[2 * nb0 + 2 * ng0 + idx];
            ptr_dux[nu0 + nx0 + ns0 + idx] = ptr_res_g[nu0 + nx0 + ns0 + idx] + ptr_gamma[nb0 + ng0 + ii] + ptr_gamma[2 * nb0 + 2 * ng0 + ns0 + idx];
            ptr_tmp2[ii] = ptr_gamma[0 + ii] - ptr_Gamma[0 + ii] * ptr_Zs_inv[0 + idx] * ptr_dux[nu0 + nx0 + idx];
            ptr_tmp3[ii] = ptr_gamma[nb0 + ng0 + ii] - ptr_Gamma[nb0 + ng0 + ii] * ptr_Zs_inv[ns0 + idx] * ptr_dux[nu0 + nx0 + ns0 + idx];  // + ???
        } else {
            ptr_tmp2[ii] = ptr_gamma[0 + ii];
            ptr_tmp3[ii] = ptr_gamma[nb0 + ng0 + ii];
        }
    }

    daxpy(nb0 + ng0, -1.0, tmp_nbgM + 3, 0, tmp_nbgM + 2, 0, tmp_nbgM + 1, 0);
}


static void EXPAND_SLACKS(int ss, struct d_ocp_qp* qp, struct d_ocp_qp_sol* qp_sol, struct d_ocp_qp_ipm_ws* ws) {

    int ii, idx;

    int nx0 = qp->dim->nx[ss];
    int nu0 = qp->dim->nu[ss];
    int nb0 = qp->dim->nb[ss];
    int ng0 = qp->dim->ng[ss];
    int ns0 = qp->dim->ns[ss];

    int* idxs_rev0 = qp->idxs_rev[ss];

    struct vec* dux = qp_sol->ux + ss;
    struct vec* dt = qp_sol->t + ss;

    struct vec* Gamma = ws->Gamma + ss;
    struct vec* Zs_inv = ws->Zs_inv + ss;

    double* ptr_Gamma = Gamma->pa;
    double* ptr_dux = dux->pa;
    double* ptr_dt = dt->pa;
    double* ptr_Zs_inv = Zs_inv->pa;

    // idxs_rev
    for (ii = 0; ii < nb0 + ng0; ii++) {
        idx = idxs_rev0[ii];
        if (idx != -1) {
            // ii  <= constr index
            // idx <= slack index
            ptr_dux[nu0 + nx0 + idx] = -ptr_Zs_inv[0 + idx] * (ptr_dux[nu0 + nx0 + idx] + ptr_Gamma[0 + ii] * ptr_dt[0 + ii]);
            ptr_dux[nu0 + nx0 + ns0 + idx] = -ptr_Zs_inv[ns0 + idx] * (ptr_dux[nu0 + nx0 + ns0 + idx] + ptr_Gamma[nb0 + ng0 + ii] * ptr_dt[nb0 + ng0 + ii]);
            ptr_dt[2 * nb0 + 2 * ng0 + idx] = ptr_dux[nu0 + nx0 + idx];
            ptr_dt[2 * nb0 + 2 * ng0 + ns0 + idx] = ptr_dux[nu0 + nx0 + ns0 + idx];
            ptr_dt[0 + ii] = ptr_dt[0 + ii] + ptr_dux[nu0 + nx0 + idx];
            ptr_dt[nb0 + ng0 + ii] = ptr_dt[nb0 + ng0 + ii] + ptr_dux[nu0 + nx0 + ns0 + idx];
        }
    }
}


// backward Riccati recursion
void d_ocp_qp_fact_solve_kkt_step(struct d_ocp_qp* qp, struct d_ocp_qp_sol* qp_sol, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws) {

    int N = qp->dim->N;
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;
    int* ns = qp->dim->ns;

    struct mat* BAbt = qp->BAbt;
    struct mat* RSQrq = qp->RSQrq;
    struct mat* DCt = qp->DCt;
    struct vec* Z = qp->Z;
    struct vec* res_g = qp->rqz;
    struct vec* res_b = qp->b;
    struct vec* res_d = qp->d;
    struct vec* res_m = qp->m;
    int** idxb = qp->idxb;

    struct vec* dux = qp_sol->ux;
    struct vec* dpi = qp_sol->pi;
    struct vec* dlam = qp_sol->lam;
    struct vec* dt = qp_sol->t;

    struct mat* L = ws->L;
    struct vec* l = ws->l;
    struct mat* AL = ws->AL;
    struct vec* Gamma = ws->Gamma;
    struct vec* gamma = ws->gamma;
    struct vec* Pb = ws->Pb;
    struct vec* Zs_inv = ws->Zs_inv;
    struct vec* tmp_nuxM = ws->tmp_nuxM;
    struct vec* tmp_nbgM = ws->tmp_nbgM;

    double *ptr0, *ptr1, *ptr2, *ptr3;

    //
    int ii, nn, ss, idx;

    struct d_core_qp_ipm_workspace* cws = ws->core_workspace;

    d_compute_Gamma_gamma_qp(res_d[0].pa, res_m[0].pa, cws);

    if (ws->square_root_alg) {
        ws->valid_ric_p = 0;

        // factorization and backward substitution

        // last stage
        ss = N;
        dtrcp_l(nu[ss] + nx[ss], RSQrq + ss, 0, 0, L + ss, 0, 0);  // TODO blasfeo_dtrcp_l with m and n, for m>=n
        // dgecp(nu[ss] + nx[ss], nu[ss] + nx[ss], RSQrq + ss, 0, 0, L + ss, 0, 0);  // TODO blasfeo_dtrcp_l with m and n, for m>=n
        ddiare(nu[ss] + nx[ss], arg->reg_prim, L + ss, 0, 0);
        drowin(nu[ss] + nx[ss], 1.0, res_g + ss, 0, L + ss, nu[ss] + nx[ss], 0);

        if (ns[ss] > 0) {
            COND_SLACKS_FACT_solve(ss, qp, qp_sol, arg, ws);
        } else if (nb[ss] + ng[ss] > 0) {
            daxpy(nb[ss] + ng[ss], 1.0, Gamma + ss, nb[ss] + ng[ss], Gamma + ss, 0, tmp_nbgM + 0, 0);
            daxpy(nb[ss] + ng[ss], -1.0, gamma + ss, nb[ss] + ng[ss], gamma + ss, 0, tmp_nbgM + 1, 0);
        }
        if (nb[ss] > 0) {
            ddiaad_sp(nb[ss], 1.0, tmp_nbgM + 0, 0, idxb[ss], L + ss, 0, 0);
            drowad_sp(nb[ss], 1.0, tmp_nbgM + 1, 0, idxb[ss], L + ss, nu[ss] + nx[ss], 0);
        }
        if (ng[ss] > 0) {
            dgemm_nd(nu[ss] + nx[ss], ng[ss], 1.0, DCt + ss, 0, 0, tmp_nbgM + 0, nb[ss], 0.0, AL + 0, 0, 0, AL + 0, 0, 0);
            drowin(ng[ss], 1.0, tmp_nbgM + 1, nb[ss], AL + 0, nu[ss] + nx[ss], 0);
            dsyrk_dpotrf_ln_mn(nu[ss] + nx[ss] + 1, nu[ss] + nx[ss], ng[ss], AL + 0, 0, 0, DCt + ss, 0, 0, L + ss, 0, 0, L + ss, 0, 0);
        } else {
            dpotrf_l_mn(nu[ss] + nx[ss] + 1, nu[ss] + nx[ss], L + ss, 0, 0, L + ss, 0, 0);
        }

        // middle stages
        for (nn = 0; nn < N; nn++) {
            ss = N - nn - 1;
            drowin(nx[ss + 1], 1.0, res_b + ss, 0, BAbt + ss, nu[ss] + nx[ss], 0);
            dtrmm_rlnn(nu[ss] + nx[ss] + 1, nx[ss + 1], 1.0, L + ss + 1, nu[ss + 1], nu[ss + 1], BAbt + ss, 0, 0, AL, 0, 0);
            drowex(nx[ss + 1], 1.0, AL, nu[ss] + nx[ss], 0, tmp_nuxM, 0);
            dtrmv_lnn(nx[ss + 1], L + ss + 1, nu[ss + 1], nu[ss + 1], tmp_nuxM, 0, Pb + ss, 0);
            dgead(1, nx[ss + 1], 1.0, L + ss + 1, nu[ss + 1] + nx[ss + 1], nu[ss + 1], AL, nu[ss] + nx[ss], 0);

            dtrcp_l(nu[ss] + nx[ss], RSQrq + ss, 0, 0, L + ss, 0, 0);
            // dgecp(nu[ss] + nx[ss], nu[ss] + nx[ss], RSQrq + ss, 0, 0, L + ss, 0, 0);
            ddiare(nu[ss] + nx[ss], arg->reg_prim, L + ss, 0, 0);
            drowin(nu[ss] + nx[ss], 1.0, res_g + ss, 0, L + ss, nu[ss] + nx[ss], 0);

            if (ns[ss] > 0) {
                COND_SLACKS_FACT_solve(ss, qp, qp_sol, arg, ws);
            } else if (nb[ss] + ng[ss] > 0) {
                daxpy(nb[ss] + ng[ss], 1.0, Gamma + ss, nb[ss] + ng[ss], Gamma + ss, 0, tmp_nbgM + 0, 0);
                daxpy(nb[ss] + ng[ss], -1.0, gamma + ss, nb[ss] + ng[ss], gamma + ss, 0, tmp_nbgM + 1, 0);
            }
            if (nb[ss] > 0) {
                ddiaad_sp(nb[ss], 1.0, tmp_nbgM + 0, 0, idxb[ss], L + ss, 0, 0);
                drowad_sp(nb[ss], 1.0, tmp_nbgM + 1, 0, idxb[ss], L + ss, nu[ss] + nx[ss], 0);
            }
            if (ng[ss] > 0) {
                dgemm_nd(nu[ss] + nx[ss], ng[ss], 1.0, DCt + ss, 0, 0, tmp_nbgM + 0, nb[ss], 0.0, AL + 0, 0, nx[ss + 1], AL + 0, 0, nx[ss + 1]);
                drowin(ng[ss], 1.0, tmp_nbgM + 1, nb[ss], AL + 0, nu[ss] + nx[ss], nx[ss + 1]);
                dgecp(nu[ss] + nx[ss], nx[ss + 1], AL + 0, 0, 0, AL + 1, 0, 0);
                dgecp(nu[ss] + nx[ss], ng[ss], DCt + ss, 0, 0, AL + 1, 0, nx[ss + 1]);
                dsyrk_dpotrf_ln_mn(nu[ss] + nx[ss] + 1, nu[ss] + nx[ss], nx[ss + 1] + ng[ss], AL + 0, 0, 0, AL + 1, 0, 0, L + ss, 0, 0, L + ss, 0, 0);
            } else {
                dsyrk_dpotrf_ln_mn(nu[ss] + nx[ss] + 1, nu[ss] + nx[ss], nx[ss + 1], AL, 0, 0, AL, 0, 0, L + ss, 0, 0, L + ss, 0, 0);
            }
        }

        // forward substitution

        // first stage
        ss = 0;
        drowex(nu[ss] + nx[ss], -1.0, L + ss, nu[ss] + nx[ss], 0, dux + ss, 0);
        dveccpsc(nu[ss] + nx[ss], -1.0, dux + ss, 0, l + ss, 0);
        drowex(nx[ss + 1], 1.0, L + ss + 1, nu[ss + 1] + nx[ss + 1], nu[ss + 1], tmp_nuxM, 0);
        dveccp(nx[ss + 1], tmp_nuxM, 0, l + ss + 1, nu[ss + 1]);  // TODO remove tmp_nuxM and use l instead !!!!!
        dtrsv_ltn(nu[ss] + nx[ss], L + ss, 0, 0, dux + ss, 0, dux + ss, 0);
        dgemv_t(nu[ss] + nx[ss], nx[ss + 1], 1.0, BAbt + ss, 0, 0, dux + ss, 0, 1.0, res_b + ss, 0, dux + ss + 1, nu[ss + 1]);
        if (arg->comp_dual_sol_eq) {
            dtrmv_ltn(nx[ss + 1], L + ss + 1, nu[ss + 1], nu[ss + 1], dux + ss + 1, nu[ss + 1], dpi + ss, 0);
            daxpy(nx[ss + 1], 1.0, tmp_nuxM, 0, dpi + ss, 0, dpi + ss, 0);
            dtrmv_lnn(nx[ss + 1], L + ss + 1, nu[ss + 1], nu[ss + 1], dpi + ss, 0, dpi + ss, 0);
        }

        // middle stages
        for (ss = 1; ss < N; ss++) {
            drowex(nu[ss], -1.0, L + ss, nu[ss] + nx[ss], 0, dux + ss, 0);
            dveccpsc(nu[ss], -1.0, dux + ss, 0, l + ss, 0);
            drowex(nx[ss + 1], 1.0, L + ss + 1, nu[ss + 1] + nx[ss + 1], nu[ss + 1], tmp_nuxM, 0);
            dveccp(nx[ss + 1], tmp_nuxM, 0, l + ss + 1, nu[ss + 1]);  // TODO remove tmp_nuxM and use l instead !!!!!
            dtrsv_ltn_mn(nu[ss] + nx[ss], nu[ss], L + ss, 0, 0, dux + ss, 0, dux + ss, 0);
            dgemv_t(nu[ss] + nx[ss], nx[ss + 1], 1.0, BAbt + ss, 0, 0, dux + ss, 0, 1.0, res_b + ss, 0, dux + (ss + 1), nu[ss + 1]);
            if (arg->comp_dual_sol_eq) {
                dtrmv_ltn(nx[ss + 1], L + ss + 1, nu[ss + 1], nu[ss + 1], dux + ss + 1, nu[ss + 1], dpi + ss, 0);
                daxpy(nx[ss + 1], 1.0, tmp_nuxM, 0, dpi + ss, 0, dpi + ss, 0);
                dtrmv_lnn(nx[ss + 1], L + ss + 1, nu[ss + 1], nu[ss + 1], dpi + ss, 0, dpi + ss, 0);
            }
        }

        ss = N;
        drowex(nu[ss], -1.0, L + ss, nu[ss] + nx[ss], 0, dux + ss, 0);
        dveccpsc(nu[ss], -1.0, dux + ss, 0, l + ss, 0);
        dtrsv_ltn_mn(nu[ss] + nx[ss], nu[ss], L + ss, 0, 0, dux + ss, 0, dux + ss, 0);

    } else  // classical algorithm
    {
        ws->valid_ric_p = 1;

        struct mat* P = ws->P;
        struct mat* Ls = ws->Ls;

        // factorization and backward substitution

        // last stage
        ss = N;
        dtrcp_l(nu[ss] + nx[ss], RSQrq + ss, 0, 0, L + ss, 0, 0);  // TODO blasfeo_dtrcp_l with m and n, for m>=n
        // dgecp(nu[ss] + nx[ss], nu[ss] + nx[ss], RSQrq + ss, 0, 0, L + ss, 0, 0);  // TODO blasfeo_dtrcp_l with m and n, for m>=n
        ddiare(nu[ss] + nx[ss], arg->reg_prim, L + ss, 0, 0);
        drowin(nu[ss] + nx[ss], 1.0, res_g + ss, 0, L + ss, nu[ss] + nx[ss], 0);

        if (ns[ss] > 0) {
            COND_SLACKS_FACT_solve(ss, qp, qp_sol, arg, ws);
        } else if (nb[ss] + ng[ss] > 0) {
            daxpy(nb[ss] + ng[ss], 1.0, Gamma + ss, nb[ss] + ng[ss], Gamma + ss, 0, tmp_nbgM + 0, 0);
            daxpy(nb[ss] + ng[ss], -1.0, gamma + ss, nb[ss] + ng[ss], gamma + ss, 0, tmp_nbgM + 1, 0);
        }
        if (nb[ss] > 0) {
            ddiaad_sp(nb[ss], 1.0, tmp_nbgM + 0, 0, idxb[ss], L + ss, 0, 0);
            drowad_sp(nb[ss], 1.0, tmp_nbgM + 1, 0, idxb[ss], L + ss, nu[ss] + nx[ss], 0);
        }
        if (ng[ss] > 0) {
            dgemm_nd(nu[ss] + nx[ss], ng[ss], 1.0, DCt + ss, 0, 0, tmp_nbgM + 0, nb[ss], 0.0, AL + 0, 0, 0, AL + 0, 0, 0);
            drowin(ng[ss], 1.0, tmp_nbgM + 1, nb[ss], AL + 0, nu[ss] + nx[ss], 0);
            dsyrk_ln_mn(nu[ss] + nx[ss] + 1, nu[ss] + nx[ss], ng[ss], 1.0, AL + 0, 0, 0, DCt + ss, 0, 0, 1.0, L + ss, 0, 0, L + ss, 0, 0);
        }
        dpotrf_l_mn(nu[ss] + nx[ss] + 1, nu[ss], L + ss, 0, 0, L + ss, 0, 0);
        dgecp(nx[ss] + 1, nu[ss], L + ss, nu[ss], 0, Ls, 0, 0);
        dsyrk_ln_mn(nx[ss] + 1, nx[ss], nu[ss], -1.0, Ls, 0, 0, Ls, 0, 0, 1.0, L + ss, nu[ss], nu[ss], P + ss, 0, 0);
        dtrtr_l(nx[ss], P + ss, 0, 0, P + ss, 0, 0);

        // middle stages
        for (nn = 0; nn < N - 1; nn++) {
            ss = N - nn - 1;
            drowin(nx[ss + 1], 1.0, res_b + ss, 0, BAbt + ss, nu[ss] + nx[ss], 0);
            dgemm_nt(nu[ss] + nx[ss] + 1, nx[ss + 1], nx[ss + 1], 1.0, BAbt + ss, 0, 0, P + ss + 1, 0, 0, 0.0, AL, 0, 0, AL, 0, 0);  // TODO symm
            drowex(nx[ss + 1], 1.0, AL, nu[ss] + nx[ss], 0, Pb + ss, 0);
            dgead(1, nx[ss + 1], 1.0, P + ss + 1, nx[ss + 1], 0, AL, nu[ss] + nx[ss], 0);

            dtrcp_l(nu[ss] + nx[ss], RSQrq + ss, 0, 0, L + ss, 0, 0);
            // dgecp(nu[ss] + nx[ss], nu[ss] + nx[ss], RSQrq + ss, 0, 0, L + ss, 0, 0);
            ddiare(nu[ss] + nx[ss], arg->reg_prim, L + ss, 0, 0);
            drowin(nu[ss] + nx[ss], 1.0, res_g + ss, 0, L + ss, nu[ss] + nx[ss], 0);

            if (ns[ss] > 0) {
                COND_SLACKS_FACT_solve(ss, qp, qp_sol, arg, ws);
            } else if (nb[ss] + ng[ss] > 0) {
                daxpy(nb[ss] + ng[ss], 1.0, Gamma + ss, nb[ss] + ng[ss], Gamma + ss, 0, tmp_nbgM + 0, 0);
                daxpy(nb[ss] + ng[ss], -1.0, gamma + ss, nb[ss] + ng[ss], gamma + ss, 0, tmp_nbgM + 1, 0);
            }
            if (nb[ss] > 0) {
                ddiaad_sp(nb[ss], 1.0, tmp_nbgM + 0, 0, idxb[ss], L + ss, 0, 0);
                drowad_sp(nb[ss], 1.0, tmp_nbgM + 1, 0, idxb[ss], L + ss, nu[ss] + nx[ss], 0);
            }
            if (ng[ss] > 0) {
                dgemm_nd(nu[ss] + nx[ss], ng[ss], 1.0, DCt + ss, 0, 0, tmp_nbgM + 0, nb[ss], 0.0, AL + 0, 0, nx[ss + 1], AL + 0, 0, nx[ss + 1]);
                drowin(ng[ss], 1.0, tmp_nbgM + 1, nb[ss], AL + 0, nu[ss] + nx[ss], nx[ss + 1]);
                dsyrk_ln_mn(nu[ss] + nx[ss] + 1, nu[ss] + nx[ss], ng[ss], 1.0, AL + 0, 0, nx[ss + 1], DCt + ss, 0, 0, 1.0, L + ss, 0, 0, L + ss, 0, 0);
            }
            dsyrk_ln_mn(nu[ss] + nx[ss] + 1, nu[ss] + nx[ss], nx[ss + 1], 1.0, AL, 0, 0, BAbt + ss, 0, 0, 1.0, L + ss, 0, 0, L + ss, 0, 0);
            dpotrf_l_mn(nu[ss] + nx[ss] + 1, nu[ss], L + ss, 0, 0, L + ss, 0, 0);
            dgecp(nx[ss] + 1, nu[ss], L + ss, nu[ss], 0, Ls, 0, 0);
            dsyrk_ln_mn(nx[ss] + 1, nx[ss], nu[ss], -1.0, Ls, 0, 0, Ls, 0, 0, 1.0, L + ss, nu[ss], nu[ss], P + ss, 0, 0);
            dtrtr_l(nx[ss], P + ss, 0, 0, P + ss, 0, 0);
        }

        // first stage: factorize P in L too
        if (N > 0) {
            ss = N - nn - 1;
            drowin(nx[ss + 1], 1.0, res_b + ss, 0, BAbt + ss, nu[ss] + nx[ss], 0);
            dgemm_nt(nu[ss] + nx[ss] + 1, nx[ss + 1], nx[ss + 1], 1.0, BAbt + ss, 0, 0, P + ss + 1, 0, 0, 0.0, AL, 0, 0, AL, 0, 0);  // TODO symm
            drowex(nx[ss + 1], 1.0, AL, nu[ss] + nx[ss], 0, Pb + ss, 0);
            dgead(1, nx[ss + 1], 1.0, P + ss + 1, nx[ss + 1], 0, AL, nu[ss] + nx[ss], 0);

            dtrcp_l(nu[ss] + nx[ss], RSQrq + ss, 0, 0, L + ss, 0, 0);
            // dgecp(nu[ss] + nx[ss], nu[ss] + nx[ss], RSQrq + ss, 0, 0, L + ss, 0, 0);
            ddiare(nu[ss] + nx[ss], arg->reg_prim, L + ss, 0, 0);
            drowin(nu[ss] + nx[ss], 1.0, res_g + ss, 0, L + ss, nu[ss] + nx[ss], 0);

            if (ns[ss] > 0) {
                COND_SLACKS_FACT_solve(ss, qp, qp_sol, arg, ws);
            } else if (nb[ss] + ng[ss] > 0) {
                daxpy(nb[ss] + ng[ss], 1.0, Gamma + ss, nb[ss] + ng[ss], Gamma + ss, 0, tmp_nbgM + 0, 0);
                daxpy(nb[ss] + ng[ss], -1.0, gamma + ss, nb[ss] + ng[ss], gamma + ss, 0, tmp_nbgM + 1, 0);
            }
            if (nb[ss] > 0) {
                ddiaad_sp(nb[ss], 1.0, tmp_nbgM + 0, 0, idxb[ss], L + ss, 0, 0);
                drowad_sp(nb[ss], 1.0, tmp_nbgM + 1, 0, idxb[ss], L + ss, nu[ss] + nx[ss], 0);
            }
            if (ng[ss] > 0) {
                dgemm_nd(nu[ss] + nx[ss], ng[ss], 1.0, DCt + ss, 0, 0, tmp_nbgM + 0, nb[ss], 0.0, AL + 0, 0, nx[ss + 1], AL + 0, 0, nx[ss + 1]);
                drowin(ng[ss], 1.0, tmp_nbgM + 1, nb[ss], AL + 0, nu[ss] + nx[ss], nx[ss + 1]);
                dsyrk_ln_mn(nu[ss] + nx[ss] + 1, nu[ss] + nx[ss], ng[ss], 1.0, AL + 0, 0, nx[ss + 1], DCt + ss, 0, 0, 1.0, L + ss, 0, 0, L + ss, 0, 0);
            }
            dsyrk_dpotrf_ln_mn(nu[ss] + nx[ss] + 1, nu[ss] + nx[ss], nx[ss + 1], AL, 0, 0, BAbt + ss, 0, 0, L + ss, 0, 0, L + ss, 0, 0);
        }

        // forward substitution

        // first stage
        ss = 0;
        drowex(nu[ss] + nx[ss], -1.0, L + ss, nu[ss] + nx[ss], 0, dux + ss, 0);
        // XXX P+0 is empty !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //		dveccpsc(nu[ss], -1.0, dux+ss, 0, l+ss, 0);
        //		drowex(nx[ss], 1.0, P+ss, nx[ss], 0, l+ss, nu[ss]);
        dveccpsc(nu[ss] + nx[ss], -1.0, dux + ss, 0, l + ss, 0);
        drowex(nx[ss + 1], 1.0, P + ss + 1, nx[ss + 1], 0, tmp_nuxM, 0);
        dveccp(nx[ss + 1], tmp_nuxM, 0, l + ss + 1, nu[ss + 1]);  // TODO remove tmp_nuxM and use l instead !!!!!
        dtrsv_ltn(nu[ss] + nx[ss], L + ss, 0, 0, dux + ss, 0, dux + ss, 0);
        dgemv_t(nu[ss] + nx[ss], nx[ss + 1], 1.0, BAbt + ss, 0, 0, dux + ss, 0, 1.0, res_b + ss, 0, dux + ss + 1, nu[ss + 1]);
        if (arg->comp_dual_sol_eq) {
            dgemv_n(nx[ss + 1], nx[ss + 1], 1.0, P + ss + 1, 0, 0, dux + ss + 1, nu[ss + 1], 1.0, tmp_nuxM, 0, dpi + ss, 0);
        }

        // middle stages
        for (ss = 1; ss < N; ss++) {
            drowex(nu[ss], -1.0, L + ss, nu[ss] + nx[ss], 0, dux + ss, 0);
            dveccpsc(nu[ss], -1.0, dux + ss, 0, l + ss, 0);
            drowex(nx[ss + 1], 1.0, P + ss + 1, nx[ss + 1], 0, tmp_nuxM, 0);
            dveccp(nx[ss + 1], tmp_nuxM, 0, l + ss + 1, nu[ss + 1]);  // TODO remove tmp_nuxM and use l instead !!!!!
            dtrsv_ltn_mn(nu[ss] + nx[ss], nu[ss], L + ss, 0, 0, dux + ss, 0, dux + ss, 0);
            dgemv_t(nu[ss] + nx[ss], nx[ss + 1], 1.0, BAbt + ss, 0, 0, dux + ss, 0, 1.0, res_b + ss, 0, dux + (ss + 1), nu[ss + 1]);
            if (arg->comp_dual_sol_eq) {
                dgemv_n(nx[ss + 1], nx[ss + 1], 1.0, P + ss + 1, 0, 0, dux + ss + 1, nu[ss + 1], 1.0, tmp_nuxM, 0, dpi + ss, 0);
            }
        }

        ss = N;
        drowex(nu[ss], -1.0, L + ss, nu[ss] + nx[ss], 0, dux + ss, 0);
        dveccpsc(nu[ss], -1.0, dux + ss, 0, l + ss, 0);
        dtrsv_ltn_mn(nu[ss] + nx[ss], nu[ss], L + ss, 0, 0, dux + ss, 0, dux + ss, 0);
    }

    for (ss = 0; ss <= N; ss++)
        dvecex_sp(nb[ss], 1.0, idxb[ss], dux + ss, 0, dt + ss, 0);
    for (ss = 0; ss <= N; ss++)
        dgemv_t(nu[ss] + nx[ss], ng[ss], 1.0, DCt + ss, 0, 0, dux + ss, 0, 0.0, dt + ss, nb[ss], dt + ss, nb[ss]);

    for (ss = 0; ss <= N; ss++) {
        dveccp(nb[ss] + ng[ss], dt + ss, 0, dt + ss, nb[ss] + ng[ss]);
        dvecsc(nb[ss] + ng[ss], -1.0, dt + ss, nb[ss] + ng[ss]);
    }

    for (ss = 0; ss <= N; ss++) {
        if (ns[ss] > 0)
            EXPAND_SLACKS(ss, qp, qp_sol, ws);
    }

    d_compute_lam_t_qp(res_d[0].pa, res_m[0].pa, dlam[0].pa, dt[0].pa, cws);
}


void d_ocp_qp_fact_lq_solve_kkt_step(struct d_ocp_qp* qp, struct d_ocp_qp_sol* qp_sol, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws) {

    // TODO find something better ???
    if (!ws->square_root_alg) {
        d_ocp_qp_fact_solve_kkt_step(qp, qp_sol, arg, ws);
    }

    int N = qp->dim->N;
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;
    int* ns = qp->dim->ns;

    struct mat* BAbt = qp->BAbt;
    struct mat* RSQrq = qp->RSQrq;
    struct mat* DCt = qp->DCt;
    struct vec* Z = qp->Z;
    struct vec* res_g = qp->rqz;
    struct vec* res_b = qp->b;
    struct vec* res_d = qp->d;
    struct vec* res_m = qp->m;
    int** idxb = qp->idxb;

    struct vec* dux = qp_sol->ux;
    struct vec* dpi = qp_sol->pi;
    struct vec* dlam = qp_sol->lam;
    struct vec* dt = qp_sol->t;

    struct mat* L = ws->L;
    struct vec* l = ws->l;
    struct mat* Lh = ws->Lh;
    struct mat* AL = ws->AL;
    struct vec* Gamma = ws->Gamma;
    struct vec* gamma = ws->gamma;
    struct vec* Pb = ws->Pb;
    struct vec* Zs_inv = ws->Zs_inv;
    struct vec* tmp_nuxM = ws->tmp_nuxM;
    struct vec* tmp_nbgM = ws->tmp_nbgM;
    struct mat* lq0 = ws->lq0;
    void* lq_work0 = ws->lq_work0;

    double *ptr0, *ptr1, *ptr2, *ptr3;

    double tmp;

    //
    int ii, nn, ss, idx;

    struct d_core_qp_ipm_workspace* cws = ws->core_workspace;

    ws->valid_ric_p = 1;

    d_compute_Gamma_gamma_qp(res_d[0].pa, res_m[0].pa, cws);

    // factorization and backward substitution

    // last stage
    ss = N;

    dveccp(nu[ss] + nx[ss], res_g + ss, 0, dux + ss, 0);

#if defined(LA_HIGH_PERFORMANCE) | defined(LA_REFERENCE)
    dgese(nu[ss] + nx[ss], nu[ss] + nx[ss] + ng[ss], 0.0, lq0, 0, nu[ss] + nx[ss]);
#else
    dgese(nu[ss] + nx[ss], 2 * nu[ss] + 2 * nx[ss] + ng[ss], 0.0, lq0, 0, 0);
#endif

    if (ns[ss] > 0) {
        COND_SLACKS_FACT_solve(ss, qp, qp_sol, arg, ws);
    } else if (nb[ss] + ng[ss] > 0) {
        daxpy(nb[ss] + ng[ss], 1.0, Gamma + ss, nb[ss] + ng[ss], Gamma + ss, 0, tmp_nbgM + 0, 0);
        daxpy(nb[ss] + ng[ss], -1.0, gamma + ss, nb[ss] + ng[ss], gamma + ss, 0, tmp_nbgM + 1, 0);
    }
    if (nb[ss] > 0) {
        for (ii = 0; ii < nb[ss]; ii++) {
            tmp = VECEL(tmp_nbgM + 0, ii);
            tmp = tmp >= 0.0 ? tmp : 0.0;
            tmp = sqrt(tmp);
            MATEL(lq0, idxb[ss][ii], nu[ss] + nx[ss] + idxb[ss][ii]) = tmp > 0.0 ? tmp : 0.0;
        }
        dvecad_sp(nb[ss], 1.0, tmp_nbgM + 1, 0, idxb[ss], dux + ss, 0);
    }
    if (ng[ss] > 0) {
        for (ii = 0; ii < ng[ss]; ii++) {
            tmp = VECEL(tmp_nbgM + 0, nb[ss] + ii);
            tmp = tmp >= 0.0 ? tmp : 0.0;
            tmp = sqrt(tmp);
            VECEL(tmp_nbgM + 0, nb[ss] + ii) = tmp;
        }
        dgemm_nd(nu[ss] + nx[ss], ng[ss], 1.0, DCt + ss, 0, 0, tmp_nbgM + 0, nb[ss], 0.0, lq0, 0, 2 * nu[ss] + 2 * nx[ss], lq0, 0, 2 * nu[ss] + 2 * nx[ss]);
        dgemv_n(nu[ss] + nx[ss], ng[ss], 1.0, DCt + ss, 0, 0, tmp_nbgM + 1, nb[ss], 1.0, dux + ss, 0, dux + ss, 0);
    }

    if (ws->use_hess_fact[ss] == 0) {
        dtrcp_l(nu[ss] + nx[ss], RSQrq + ss, 0, 0, Lh + ss, 0, 0);
        ddiare(nu[ss] + nx[ss], arg->reg_prim, Lh + ss, 0, 0);
        dpotrf_l(nu[ss] + nx[ss], Lh + ss, 0, 0, Lh + ss, 0, 0);
        ws->use_hess_fact[ss] = 1;
    }

#if defined(LA_HIGH_PERFORMANCE) | defined(LA_REFERENCE)
    dtrcp_l(nu[ss] + nx[ss], Lh + ss, 0, 0, L + ss, 0, 0);
    //	dgelqf_pd(nu[ss]+nx[ss], 2*nu[ss]+2*nx[ss]+ng[ss], lq0, 0, 0, lq0, 0, 0, lq_work0);
    dgelqf_pd_lla(nu[ss] + nx[ss], ng[ss], L + ss, 0, 0, lq0, 0, nu[ss] + nx[ss], lq0, 0, 2 * nu[ss] + 2 * nx[ss], lq_work0);  // TODO reduce lq1 size !!!
#else
    dtrcp_l(nu[ss] + nx[ss], Lh + ss, 0, 0, lq0, 0, 0);
    dgelqf(nu[ss] + nx[ss], 2 * nu[ss] + 2 * nx[ss] + ng[ss], lq0, 0, 0, lq0, 0, 0, lq_work0);
    dtrcp_l(nu[ss] + nx[ss], lq0, 0, 0, L + ss, 0, 0);
    for (ii = 0; ii < nu[ss] + nx[ss]; ii++)
        if (MATEL(L + ss, ii, ii) < 0)
            dcolsc(nu[ss] + nx[ss] - ii, -1.0, L + ss, ii, ii);
#endif

    dtrsv_lnn_mn(nu[ss] + nx[ss], nu[ss], L + ss, 0, 0, dux + ss, 0, dux + ss, 0);


    // middle stages
    for (nn = 0; nn < N - 1; nn++) {
        ss = N - nn - 1;

#if defined(LA_HIGH_PERFORMANCE) | defined(LA_REFERENCE)
        dgese(nu[ss] + nx[ss], nu[ss] + nx[ss] + ng[ss], 0.0, lq0, 0, nu[ss] + nx[ss]);
#else
        dgese(nu[ss] + nx[ss], 2 * nu[ss] + 2 * nx[ss] + nx[ss + 1] + ng[ss], 0.0, lq0, 0, 0);
#endif

        dtrmm_rlnn(nu[ss] + nx[ss], nx[ss + 1], 1.0, L + ss + 1, nu[ss + 1], nu[ss + 1], BAbt + ss, 0, 0, lq0, 0, 2 * nu[ss] + 2 * nx[ss] + ng[ss]);
        dtrmv_ltn(nx[ss + 1], L + ss + 1, nu[ss + 1], nu[ss + 1], res_b + ss, 0, Pb + ss, 0);
        dtrmv_lnn(nx[ss + 1], L + ss + 1, nu[ss + 1], nu[ss + 1], Pb + ss, 0, Pb + ss, 0);

        dveccp(nu[ss] + nx[ss], res_g + ss, 0, dux + ss, 0);
        daxpy(nx[ss + 1], 1.0, dux + ss + 1, nu[ss + 1], Pb + ss, 0, tmp_nuxM, 0);
        dgemv_n(nu[ss] + nx[ss], nx[ss + 1], 1.0, BAbt + ss, 0, 0, tmp_nuxM, 0, 1.0, dux + ss, 0, dux + ss, 0);

        if (ns[ss] > 0) {
            COND_SLACKS_FACT_solve(ss, qp, qp_sol, arg, ws);
        } else if (nb[ss] + ng[ss] > 0) {
            daxpy(nb[ss] + ng[ss], 1.0, Gamma + ss, nb[ss] + ng[ss], Gamma + ss, 0, tmp_nbgM + 0, 0);
            daxpy(nb[ss] + ng[ss], -1.0, gamma + ss, nb[ss] + ng[ss], gamma + ss, 0, tmp_nbgM + 1, 0);
        }
        if (nb[ss] > 0) {
            for (ii = 0; ii < nb[ss]; ii++) {
                tmp = VECEL(tmp_nbgM + 0, ii);
                tmp = tmp >= 0.0 ? tmp : 0.0;
                tmp = sqrt(tmp);
                MATEL(lq0, idxb[ss][ii], nu[ss] + nx[ss] + idxb[ss][ii]) = tmp > 0.0 ? tmp : 0.0;
            }
            dvecad_sp(nb[ss], 1.0, tmp_nbgM + 1, 0, idxb[ss], dux + ss, 0);
        }
        if (ng[ss] > 0) {
            for (ii = 0; ii < ng[ss]; ii++) {
                tmp = VECEL(tmp_nbgM + 0, nb[ss] + ii);
                tmp = tmp >= 0.0 ? tmp : 0.0;
                tmp = sqrt(tmp);
                VECEL(tmp_nbgM + 0, nb[ss] + ii) = tmp;
            }
            dgemm_nd(nu[ss] + nx[ss], ng[ss], 1.0, DCt + ss, 0, 0, tmp_nbgM + 0, nb[ss], 0.0, lq0, 0, 2 * nu[ss] + 2 * nx[ss], lq0, 0, 2 * nu[ss] + 2 * nx[ss]);
            dgemv_n(nu[ss] + nx[ss], ng[ss], 1.0, DCt + ss, 0, 0, tmp_nbgM + 1, nb[ss], 1.0, dux + ss, 0, dux + ss, 0);
        }

        if (ws->use_hess_fact[ss] == 0) {
            dtrcp_l(nu[ss] + nx[ss], RSQrq + ss, 0, 0, Lh + ss, 0, 0);
            ddiare(nu[ss] + nx[ss], arg->reg_prim, Lh + ss, 0, 0);
            dpotrf_l(nu[ss] + nx[ss], Lh + ss, 0, 0, Lh + ss, 0, 0);
            ws->use_hess_fact[ss] = 1;
        }

#if defined(LA_HIGH_PERFORMANCE) | defined(LA_REFERENCE)
        dtrcp_l(nu[ss] + nx[ss], Lh + ss, 0, 0, L + ss, 0, 0);
        //		dgelqf_pd(nu[ss]+nx[ss], 2*nu[ss]+2*nx[ss]+nx[ss+1]+ng[ss], lq0, 0, 0, lq0, 0, 0, lq_work0);
        //		dgelqf_pd_la(nu[ss]+nx[ss], nu[ss]+nx[ss]+nx[ss+1]+ng[ss], L+ss, 0, 0, lq0, 0, nu[ss]+nx[ss], lq_work0); // TODO reduce lq1 size !!!
        dgelqf_pd_lla(nu[ss] + nx[ss], nx[ss + 1] + ng[ss], L + ss, 0, 0, lq0, 0, nu[ss] + nx[ss], lq0, 0, 2 * nu[ss] + 2 * nx[ss], lq_work0);  // TODO reduce lq1 size !!!
#else
        dtrcp_l(nu[ss] + nx[ss], Lh + ss, 0, 0, lq0, 0, 0);
        dgelqf(nu[ss] + nx[ss], 2 * nu[ss] + 2 * nx[ss] + nx[ss + 1] + ng[ss], lq0, 0, 0, lq0, 0, 0, lq_work0);
        dtrcp_l(nu[ss] + nx[ss], lq0, 0, 0, L + ss, 0, 0);
        for (ii = 0; ii < nu[ss] + nx[ss]; ii++)
            if (MATEL(L + ss, ii, ii) < 0)
                dcolsc(nu[ss] + nx[ss] - ii, -1.0, L + ss, ii, ii);
#endif

        dtrsv_lnn_mn(nu[ss] + nx[ss], nu[ss], L + ss, 0, 0, dux + ss, 0, dux + ss, 0);
    }

    // first stage
    nn = N - 1;
    ss = N - nn - 1;

#if defined(LA_HIGH_PERFORMANCE) | defined(LA_REFERENCE)
    dgese(nu[ss] + nx[ss], nu[ss] + nx[ss] + ng[ss], 0.0, lq0, 0, nu[ss] + nx[ss]);
#else
    dgese(nu[ss] + nx[ss], 2 * nu[ss] + 2 * nx[ss] + nx[ss + 1] + ng[ss], 0.0, lq0, 0, 0);
#endif

    dtrmm_rlnn(nu[ss] + nx[ss], nx[ss + 1], 1.0, L + ss + 1, nu[ss + 1], nu[ss + 1], BAbt + ss, 0, 0, lq0, 0, 2 * nu[ss] + 2 * nx[ss] + ng[ss]);
    dtrmv_ltn(nx[ss + 1], L + ss + 1, nu[ss + 1], nu[ss + 1], res_b + ss, 0, Pb + ss, 0);
    dtrmv_lnn(nx[ss + 1], L + ss + 1, nu[ss + 1], nu[ss + 1], Pb + ss, 0, Pb + ss, 0);

    dveccp(nu[ss] + nx[ss], res_g + ss, 0, dux + ss, 0);
    daxpy(nx[ss + 1], 1.0, dux + ss + 1, nu[ss + 1], Pb + ss, 0, tmp_nuxM, 0);
    dgemv_n(nu[ss] + nx[ss], nx[ss + 1], 1.0, BAbt + ss, 0, 0, tmp_nuxM, 0, 1.0, dux + ss, 0, dux + ss, 0);

    if (ns[ss] > 0) {
        COND_SLACKS_FACT_solve(ss, qp, qp_sol, arg, ws);
    } else if (nb[ss] + ng[ss] > 0) {
        daxpy(nb[ss] + ng[ss], 1.0, Gamma + ss, nb[ss] + ng[ss], Gamma + ss, 0, tmp_nbgM + 0, 0);
        daxpy(nb[ss] + ng[ss], -1.0, gamma + ss, nb[ss] + ng[ss], gamma + ss, 0, tmp_nbgM + 1, 0);
    }
    if (nb[ss] > 0) {
        for (ii = 0; ii < nb[ss]; ii++) {
            tmp = VECEL(tmp_nbgM + 0, ii);
            tmp = tmp >= 0.0 ? tmp : 0.0;
            tmp = sqrt(tmp);
            MATEL(lq0, idxb[ss][ii], nu[ss] + nx[ss] + idxb[ss][ii]) = tmp > 0.0 ? tmp : 0.0;
        }
        dvecad_sp(nb[ss], 1.0, tmp_nbgM + 1, 0, idxb[ss], dux + ss, 0);
    }
    if (ng[ss] > 0) {
        for (ii = 0; ii < ng[ss]; ii++) {
            tmp = VECEL(tmp_nbgM + 0, nb[ss] + ii);
            tmp = tmp >= 0.0 ? tmp : 0.0;
            tmp = sqrt(tmp);
            VECEL(tmp_nbgM + 0, nb[ss] + ii) = tmp;
        }
        dgemm_nd(nu[ss] + nx[ss], ng[ss], 1.0, DCt + ss, 0, 0, tmp_nbgM + 0, nb[ss], 0.0, lq0, 0, 2 * nu[ss] + 2 * nx[ss], lq0, 0, 2 * nu[ss] + 2 * nx[ss]);
        dgemv_n(nu[ss] + nx[ss], ng[ss], 1.0, DCt + ss, 0, 0, tmp_nbgM + 1, nb[ss], 1.0, dux + ss, 0, dux + ss, 0);
    }

    if (ws->use_hess_fact[ss] == 0) {
        dtrcp_l(nu[ss] + nx[ss], RSQrq + ss, 0, 0, Lh + ss, 0, 0);
        ddiare(nu[ss] + nx[ss], arg->reg_prim, Lh + ss, 0, 0);
        dpotrf_l(nu[ss] + nx[ss], Lh + ss, 0, 0, Lh + ss, 0, 0);
        ws->use_hess_fact[ss] = 1;
    }

#if defined(LA_HIGH_PERFORMANCE) | defined(LA_REFERENCE)
    dtrcp_l(nu[ss] + nx[ss], Lh + ss, 0, 0, L + ss, 0, 0);
    //	dgelqf_pd(nu[ss]+nx[ss], 2*nu[ss]+2*nx[ss]+nx[ss+1]+ng[ss], lq0, 0, 0, lq0, 0, 0, lq_work0);
    dgelqf_pd_lla(nu[ss] + nx[ss], nx[ss + 1] + ng[ss], L + ss, 0, 0, lq0, 0, nu[ss] + nx[ss], lq0, 0, 2 * nu[ss] + 2 * nx[ss], lq_work0);  // TODO reduce lq1 size !!!
#else
    dtrcp_l(nu[ss] + nx[ss], Lh + ss, 0, 0, lq0, 0, 0);
    dgelqf(nu[ss] + nx[ss], 2 * nu[ss] + 2 * nx[ss] + nx[ss + 1] + ng[ss], lq0, 0, 0, lq0, 0, 0, lq_work0);
    dtrcp_l(nu[ss] + nx[ss], lq0, 0, 0, L + ss, 0, 0);
    for (ii = 0; ii < nu[ss] + nx[ss]; ii++)
        if (MATEL(L + ss, ii, ii) < 0)
            dcolsc(nu[ss] + nx[ss] - ii, -1.0, L + ss, ii, ii);
#endif

    dtrsv_lnn(nu[ss] + nx[ss], L + ss, 0, 0, dux + ss, 0, dux + ss, 0);

    // forward substitution

    // first stage
    ss = 0;
    dveccp(nu[ss] + nx[ss], dux + ss, 0, l + ss, 0);
    dveccp(nx[ss + 1], dux + ss + 1, nu[ss + 1], dpi + ss, 0);
    dveccp(nx[ss + 1], dux + ss + 1, nu[ss + 1], l + ss + 1, nu[ss + 1]);
    dvecsc(nu[ss] + nx[ss], -1.0, dux + ss, 0);
    dtrsv_ltn(nu[ss] + nx[ss], L + ss, 0, 0, dux + ss, 0, dux + ss, 0);
    dgemv_t(nu[ss] + nx[ss], nx[ss + 1], 1.0, BAbt + ss, 0, 0, dux + ss, 0, 1.0, res_b + ss, 0, dux + ss + 1, nu[ss + 1]);
    dveccp(nx[ss + 1], dux + ss + 1, nu[ss + 1], tmp_nuxM, 0);
    dtrmv_ltn(nx[ss + 1], L + ss + 1, nu[ss + 1], nu[ss + 1], tmp_nuxM, 0, tmp_nuxM, 0);
    dtrmv_lnn(nx[ss + 1], L + ss + 1, nu[ss + 1], nu[ss + 1], tmp_nuxM, 0, tmp_nuxM, 0);
    daxpy(nx[ss + 1], 1.0, tmp_nuxM, 0, dpi + ss, 0, dpi + ss, 0);

    // middle stages
    for (ss = 1; ss < N; ss++) {
        dveccp(nu[ss], dux + ss, 0, l + ss, 0);
        dveccp(nx[ss + 1], dux + ss + 1, nu[ss + 1], dpi + ss, 0);
        dveccp(nx[ss + 1], dux + ss + 1, nu[ss + 1], l + ss + 1, nu[ss + 1]);
        dvecsc(nu[ss], -1.0, dux + ss, 0);
        dtrsv_ltn_mn(nu[ss] + nx[ss], nu[ss], L + ss, 0, 0, dux + ss, 0, dux + ss, 0);
        dgemv_t(nu[ss] + nx[ss], nx[ss + 1], 1.0, BAbt + ss, 0, 0, dux + ss, 0, 1.0, res_b + ss, 0, dux + ss + 1, nu[ss + 1]);
        dveccp(nx[ss + 1], dux + ss + 1, nu[ss + 1], tmp_nuxM, 0);
        dtrmv_ltn(nx[ss + 1], L + ss + 1, nu[ss + 1], nu[ss + 1], tmp_nuxM, 0, tmp_nuxM, 0);
        dtrmv_lnn(nx[ss + 1], L + ss + 1, nu[ss + 1], nu[ss + 1], tmp_nuxM, 0, tmp_nuxM, 0);
        daxpy(nx[ss + 1], 1.0, tmp_nuxM, 0, dpi + ss, 0, dpi + ss, 0);
    }

    ss = N;
    dvecsc(nu[ss], -1.0, dux + ss, 0);
    dveccp(nu[ss], dux + ss, 0, l + ss, 0);
    dtrsv_ltn_mn(nu[ss] + nx[ss], nu[ss], L + ss, 0, 0, dux + ss, 0, dux + ss, 0);


    for (ss = 0; ss <= N; ss++)
        dvecex_sp(nb[ss], 1.0, idxb[ss], dux + ss, 0, dt + ss, 0);
    for (ss = 0; ss <= N; ss++)
        dgemv_t(nu[ss] + nx[ss], ng[ss], 1.0, DCt + ss, 0, 0, dux + ss, 0, 0.0, dt + ss, nb[ss], dt + ss, nb[ss]);

    for (ss = 0; ss <= N; ss++) {
        dveccp(nb[ss] + ng[ss], dt + ss, 0, dt + ss, nb[ss] + ng[ss]);
        dvecsc(nb[ss] + ng[ss], -1.0, dt + ss, nb[ss] + ng[ss]);
    }

    for (ss = 0; ss <= N; ss++) {
        if (ns[ss] > 0)
            EXPAND_SLACKS(ss, qp, qp_sol, ws);
    }

    d_compute_lam_t_qp(res_d[0].pa, res_m[0].pa, dlam[0].pa, dt[0].pa, cws);
}


// backward Riccati recursion
void d_ocp_qp_solve_kkt_step(struct d_ocp_qp* qp, struct d_ocp_qp_sol* qp_sol, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws) {

    int N = qp->dim->N;
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;
    int* ns = qp->dim->ns;

    struct mat* BAbt = qp->BAbt;
    //	struct mat *RSQrq = qp->RSQrq;
    struct mat* DCt = qp->DCt;
    struct vec* res_g = qp->rqz;
    struct vec* res_b = qp->b;
    struct vec* res_d = qp->d;
    struct vec* res_m = qp->m;
    int** idxb = qp->idxb;

    struct vec* dux = qp_sol->ux;
    struct vec* dpi = qp_sol->pi;
    struct vec* dlam = qp_sol->lam;
    struct vec* dt = qp_sol->t;

    struct mat* L = ws->L;
    struct vec* l = ws->l;
    struct vec* gamma = ws->gamma;
    struct vec* Pb = ws->Pb;
    struct vec* tmp_nuxM = ws->tmp_nuxM;
    struct vec* tmp_nbgM = ws->tmp_nbgM;

    //
    int ss, nn, ii;

    struct d_core_qp_ipm_workspace* cws = ws->core_workspace;

    // printf("\nin solve\n");
    d_compute_gamma_qp(res_d[0].pa, res_m[0].pa, cws);

    if (ws->square_root_alg) {
        ws->valid_ric_p = 1;

        // backward substitution

        // last stage
        ss = N;
        // blasfeo_print_exp_tran_dvec(2*nb[ss]+2*ng[ss], gamma+ss, 0);
        dveccp(nu[ss] + nx[ss], res_g + ss, 0, dux + ss, 0);
        // blasfeo_print_exp_tran_dvec(nu[ss]+nx[ss], dux+ss, 0);
        if (ns[ss] > 0) {
            COND_SLACKS_solve(ss, qp, qp_sol, ws);
        } else if (nb[ss] + ng[ss] > 0) {
            daxpy(nb[ss] + ng[ss], -1.0, gamma + ss, nb[ss] + ng[ss], gamma + ss, 0, tmp_nbgM + 1, 0);
        }
        if (nb[ss] > 0) {
            dvecad_sp(nb[ss], 1.0, tmp_nbgM + 1, 0, idxb[ss], dux + ss, 0);
        }
        // blasfeo_print_exp_tran_dvec(nu[ss]+nx[ss], dux+ss, 0);
        if (ng[ss] > 0) {
            dgemv_n(nu[ss] + nx[ss], ng[ss], 1.0, DCt + ss, 0, 0, tmp_nbgM + 1, nb[ss], 1.0, dux + ss, 0, dux + ss, 0);
        }
        // blasfeo_print_exp_tran_dvec(nu[ss]+nx[ss], dux+ss, 0);
        dtrsv_lnn_mn(nu[ss] + nx[ss], nu[ss], L + ss, 0, 0, dux + ss, 0, dux + ss, 0);
        // blasfeo_print_exp_tran_dvec(nu[ss]+nx[ss], dux+ss, 0);

        // middle stages
        for (nn = 0; nn < N - 1; nn++) {
            ss = N - nn - 1;
            dveccp(nu[ss] + nx[ss], res_g + ss, 0, dux + ss, 0);
            if (ns[ss] > 0) {
                COND_SLACKS_solve(ss, qp, qp_sol, ws);
            } else if (nb[ss] + ng[ss] > 0) {
                daxpy(nb[ss] + ng[ss], -1.0, gamma + ss, nb[ss] + ng[ss], gamma + ss, 0, tmp_nbgM + 1, 0);
            }
            if (nb[ss] > 0) {
                dvecad_sp(nb[ss], 1.0, tmp_nbgM + 1, 0, idxb[ss], dux + ss, 0);
            }
            if (ng[ss] > 0) {
                dgemv_n(nu[ss] + nx[ss], ng[ss], 1.0, DCt + ss, 0, 0, tmp_nbgM + 1, nb[ss], 1.0, dux + ss, 0, dux + ss, 0);
            }
            if (ws->use_Pb) {
                daxpy(nx[ss + 1], 1.0, dux + ss + 1, nu[ss + 1], Pb + ss, 0, tmp_nuxM, 0);
            } else {
                dtrmv_ltn(nx[ss + 1], L + ss + 1, nu[ss + 1], nu[ss + 1], res_b + ss, 0, tmp_nuxM, 0);
                dtrmv_lnn(nx[ss + 1], L + ss + 1, nu[ss + 1], nu[ss + 1], tmp_nuxM, 0, tmp_nuxM, 0);
                daxpy(nx[ss + 1], 1.0, dux + ss + 1, nu[ss + 1], tmp_nuxM, 0, tmp_nuxM, 0);
            }
            dgemv_n(nu[ss] + nx[ss], nx[ss + 1], 1.0, BAbt + ss, 0, 0, tmp_nuxM, 0, 1.0, dux + ss, 0, dux + ss, 0);
            dtrsv_lnn_mn(nu[ss] + nx[ss], nu[ss], L + ss, 0, 0, dux + ss, 0, dux + ss, 0);
        }

        // first stage
        nn = N - 1;
        ss = N - nn - 1;
        dveccp(nu[ss] + nx[ss], res_g + ss, 0, dux + ss, 0);
        if (ns[ss] > 0) {
            COND_SLACKS_solve(ss, qp, qp_sol, ws);
        } else if (nb[ss] + ng[ss] > 0) {
            daxpy(nb[ss] + ng[ss], -1.0, gamma + ss, nb[ss] + ng[ss], gamma + ss, 0, tmp_nbgM + 1, 0);
        }
        if (nb[ss] > 0) {
            dvecad_sp(nb[ss], 1.0, tmp_nbgM + 1, 0, idxb[ss], dux + ss, 0);
        }
        if (ng[ss] > 0) {
            dgemv_n(nu[ss] + nx[ss], ng[ss], 1.0, DCt + ss, 0, 0, tmp_nbgM + 1, nb[ss], 1.0, dux + ss, 0, dux + ss, 0);
        }
        if (ws->use_Pb) {
            daxpy(nx[ss + 1], 1.0, dux + ss + 1, nu[ss + 1], Pb + ss, 0, tmp_nuxM, 0);
        } else {
            dtrmv_ltn(nx[ss + 1], L + ss + 1, nu[ss + 1], nu[ss + 1], res_b + ss, 0, tmp_nuxM, 0);
            dtrmv_lnn(nx[ss + 1], L + ss + 1, nu[ss + 1], nu[ss + 1], tmp_nuxM, 0, tmp_nuxM, 0);
            daxpy(nx[ss + 1], 1.0, dux + ss + 1, nu[ss + 1], tmp_nuxM, 0, tmp_nuxM, 0);
        }
        dgemv_n(nu[ss] + nx[ss], nx[ss + 1], 1.0, BAbt + ss, 0, 0, tmp_nuxM, 0, 1.0, dux + ss, 0, dux + ss, 0);
        dtrsv_lnn(nu[ss] + nx[ss], L + ss, 0, 0, dux + ss, 0, dux + ss, 0);

        // forward substitution

        // first stage
        ss = 0;
        if (arg->comp_dual_sol_eq) {
            dveccp(nx[ss + 1], dux + ss + 1, nu[ss + 1], dpi + ss, 0);
        }
        dveccp(nu[ss] + nx[ss], dux + ss, 0, l + ss, 0);
        dveccp(nx[ss + 1], dux + ss + 1, nu[ss + 1], l + ss + 1, nu[ss + 1]);
        dvecsc(nu[ss] + nx[ss], -1.0, dux + ss, 0);
        dtrsv_ltn(nu[ss] + nx[ss], L + ss, 0, 0, dux + ss, 0, dux + ss, 0);
        dgemv_t(nu[ss] + nx[ss], nx[ss + 1], 1.0, BAbt + ss, 0, 0, dux + ss, 0, 1.0, res_b + ss, 0, dux + ss + 1, nu[ss + 1]);
        if (arg->comp_dual_sol_eq) {
            dtrmv_ltn(nx[ss + 1], L + ss + 1, nu[ss + 1], nu[ss + 1], dux + ss + 1, nu[ss + 1], tmp_nuxM, 0);
            dtrmv_lnn(nx[ss + 1], L + ss + 1, nu[ss + 1], nu[ss + 1], tmp_nuxM, 0, tmp_nuxM, 0);
            daxpy(nx[ss + 1], 1.0, tmp_nuxM, 0, dpi + ss, 0, dpi + ss, 0);
        }

        // middle stages
        for (ss = 1; ss < N; ss++) {
            if (arg->comp_dual_sol_eq) {
                dveccp(nx[ss + 1], dux + ss + 1, nu[ss + 1], dpi + ss, 0);
            }
            dveccp(nu[ss], dux + ss, 0, l + ss, 0);
            dveccp(nx[ss + 1], dux + ss + 1, nu[ss + 1], l + ss + 1, nu[ss + 1]);
            dvecsc(nu[ss], -1.0, dux + ss, 0);
            dtrsv_ltn_mn(nu[ss] + nx[ss], nu[ss], L + ss, 0, 0, dux + ss, 0, dux + ss, 0);
            dgemv_t(nu[ss] + nx[ss], nx[ss + 1], 1.0, BAbt + ss, 0, 0, dux + ss, 0, 1.0, res_b + ss, 0, dux + ss + 1, nu[ss + 1]);
            if (arg->comp_dual_sol_eq) {
                dtrmv_ltn(nx[ss + 1], L + ss + 1, nu[ss + 1], nu[ss + 1], dux + ss + 1, nu[ss + 1], tmp_nuxM, 0);
                dtrmv_lnn(nx[ss + 1], L + ss + 1, nu[ss + 1], nu[ss + 1], tmp_nuxM, 0, tmp_nuxM, 0);
                daxpy(nx[ss + 1], 1.0, tmp_nuxM, 0, dpi + ss, 0, dpi + ss, 0);
            }
        }

        ss = N;
        dveccp(nu[ss], dux + ss, 0, l + ss, 0);
        dvecsc(nu[ss], -1.0, dux + ss, 0);
        dtrsv_ltn_mn(nu[ss] + nx[ss], nu[ss], L + ss, 0, 0, dux + ss, 0, dux + ss, 0);

    } else  // classical algirthm
    {
        ws->valid_ric_p = 1;

        struct mat* P = ws->P;

        // backward substitution

        // last stage
        ss = N;
        // blasfeo_print_exp_tran_dvec(2*nb[ss]+2*ng[ss], gamma+ss, 0);
        dveccp(nu[ss] + nx[ss], res_g + ss, 0, dux + ss, 0);
        // blasfeo_print_exp_tran_dvec(nu[ss]+nx[ss], dux+ss, 0);
        if (ns[ss] > 0) {
            COND_SLACKS_solve(ss, qp, qp_sol, ws);
        } else if (nb[ss] + ng[ss] > 0) {
            daxpy(nb[ss] + ng[ss], -1.0, gamma + ss, nb[ss] + ng[ss], gamma + ss, 0, tmp_nbgM + 1, 0);
        }
        if (nb[ss] > 0) {
            dvecad_sp(nb[ss], 1.0, tmp_nbgM + 1, 0, idxb[ss], dux + ss, 0);
        }
        // blasfeo_print_exp_tran_dvec(nu[ss]+nx[ss], dux+ss, 0);
        if (ng[ss] > 0) {
            dgemv_n(nu[ss] + nx[ss], ng[ss], 1.0, DCt + ss, 0, 0, tmp_nbgM + 1, nb[ss], 1.0, dux + ss, 0, dux + ss, 0);
        }
        // blasfeo_print_dmat(nu[ss]+nx[ss], nu[ss], L+ss, 0, 0);
        // blasfeo_print_exp_tran_dvec(nu[ss]+nx[ss], dux+ss, 0);
        dtrsv_lnn_mn(nu[ss] + nx[ss], nu[ss], L + ss, 0, 0, dux + ss, 0, dux + ss, 0);
        // blasfeo_print_exp_tran_dvec(nu[ss]+nx[ss], dux+ss, 0);

        // middle stages
        for (nn = 0; nn < N - 1; nn++) {
            ss = N - nn - 1;
            dveccp(nu[ss] + nx[ss], res_g + ss, 0, dux + ss, 0);
            if (ns[ss] > 0) {
                COND_SLACKS_solve(ss, qp, qp_sol, ws);
            } else if (nb[ss] + ng[ss] > 0) {
                daxpy(nb[ss] + ng[ss], -1.0, gamma + ss, nb[ss] + ng[ss], gamma + ss, 0, tmp_nbgM + 1, 0);
            }
            if (nb[ss] > 0) {
                dvecad_sp(nb[ss], 1.0, tmp_nbgM + 1, 0, idxb[ss], dux + ss, 0);
            }
            if (ng[ss] > 0) {
                dgemv_n(nu[ss] + nx[ss], ng[ss], 1.0, DCt + ss, 0, 0, tmp_nbgM + 1, nb[ss], 1.0, dux + ss, 0, dux + ss, 0);
            }
            if (ws->use_Pb) {
                daxpy(nx[ss + 1], 1.0, dux + ss + 1, nu[ss + 1], Pb + ss, 0, tmp_nuxM, 0);
            } else {
                dgemv_n(nx[ss + 1], nx[ss + 1], 1.0, P + ss + 1, 0, 0, res_b + ss, 0, 1.0, dux + ss + 1, nu[ss + 1], tmp_nuxM, 0);
            }
            dgemv_n(nu[ss] + nx[ss], nx[ss + 1], 1.0, BAbt + ss, 0, 0, tmp_nuxM, 0, 1.0, dux + ss, 0, dux + ss, 0);
            // blasfeo_print_dmat(nu[ss]+nx[ss], nu[ss], L+ss, 0, 0);
            // blasfeo_print_exp_tran_dvec(nu[ss]+nx[ss], dux+ss, 0);
            dtrsv_lnn_mn(nu[ss] + nx[ss], nu[ss], L + ss, 0, 0, dux + ss, 0, dux + ss, 0);
            // blasfeo_print_exp_tran_dvec(nu[ss]+nx[ss], dux+ss, 0);
            // exit(1);
        }

        // first stage
        nn = N - 1;
        ss = N - nn - 1;
        dveccp(nu[ss] + nx[ss], res_g + ss, 0, dux + ss, 0);
        if (ns[ss] > 0) {
            COND_SLACKS_solve(ss, qp, qp_sol, ws);
        } else if (nb[ss] + ng[ss] > 0) {
            daxpy(nb[ss] + ng[ss], -1.0, gamma + ss, nb[ss] + ng[ss], gamma + ss, 0, tmp_nbgM + 1, 0);
        }
        if (nb[ss] > 0) {
            dvecad_sp(nb[ss], 1.0, tmp_nbgM + 1, 0, idxb[ss], dux + ss, 0);
        }
        if (ng[ss] > 0) {
            dgemv_n(nu[ss] + nx[ss], ng[ss], 1.0, DCt + ss, 0, 0, tmp_nbgM + 1, nb[ss], 1.0, dux + ss, 0, dux + ss, 0);
        }
        if (ws->use_Pb) {
            daxpy(nx[ss + 1], 1.0, dux + ss + 1, nu[ss + 1], Pb + ss, 0, tmp_nuxM, 0);
        } else {
            dgemv_n(nx[ss + 1], nx[ss + 1], 1.0, P + ss + 1, 0, 0, res_b + ss, 0, 1.0, dux + ss + 1, nu[ss + 1], tmp_nuxM, 0);
        }
        dgemv_n(nu[ss] + nx[ss], nx[ss + 1], 1.0, BAbt + ss, 0, 0, tmp_nuxM, 0, 1.0, dux + ss, 0, dux + ss, 0);
        dtrsv_lnn(nu[ss] + nx[ss], L + ss, 0, 0, dux + ss, 0, dux + ss, 0);

        // forward substitution

        // first stage
        ss = 0;
        if (arg->comp_dual_sol_eq) {
            dveccp(nx[ss + 1], dux + ss + 1, nu[ss + 1], dpi + ss, 0);
        }
        dveccp(nu[ss] + nx[ss], dux + ss, 0, l + ss, 0);
        dveccp(nx[ss + 1], dux + ss + 1, nu[ss + 1], l + ss + 1, nu[ss + 1]);
        dvecsc(nu[ss] + nx[ss], -1.0, dux + ss, 0);
        dtrsv_ltn(nu[ss] + nx[ss], L + ss, 0, 0, dux + ss, 0, dux + ss, 0);
        dgemv_t(nu[ss] + nx[ss], nx[ss + 1], 1.0, BAbt + ss, 0, 0, dux + ss, 0, 1.0, res_b + ss, 0, dux + ss + 1, nu[ss + 1]);
        if (arg->comp_dual_sol_eq) {
            dgemv_n(nx[ss + 1], nx[ss + 1], 1.0, P + ss + 1, 0, 0, dux + ss + 1, nu[ss + 1], 1.0, dpi + ss, 0, dpi + ss, 0);
        }

        // middle stages
        for (ss = 1; ss < N; ss++) {
            if (arg->comp_dual_sol_eq) {
                dveccp(nx[ss + 1], dux + ss + 1, nu[ss + 1], dpi + ss, 0);
            }
            dveccp(nu[ss], dux + ss, 0, l + ss, 0);
            dveccp(nx[ss + 1], dux + ss + 1, nu[ss + 1], l + ss + 1, nu[ss + 1]);
            dvecsc(nu[ss], -1.0, dux + ss, 0);
            dtrsv_ltn_mn(nu[ss] + nx[ss], nu[ss], L + ss, 0, 0, dux + ss, 0, dux + ss, 0);
            dgemv_t(nu[ss] + nx[ss], nx[ss + 1], 1.0, BAbt + ss, 0, 0, dux + ss, 0, 1.0, res_b + ss, 0, dux + ss + 1, nu[ss + 1]);
            if (arg->comp_dual_sol_eq) {
                dgemv_n(nx[ss + 1], nx[ss + 1], 1.0, P + ss + 1, 0, 0, dux + ss + 1, nu[ss + 1], 1.0, dpi + ss, 0, dpi + ss, 0);
            }
        }

        ss = N;
        dveccp(nu[ss], dux + ss, 0, l + ss, 0);
        dvecsc(nu[ss], -1.0, dux + ss, 0);
        dtrsv_ltn_mn(nu[ss] + nx[ss], nu[ss], L + ss, 0, 0, dux + ss, 0, dux + ss, 0);
    }

    for (ss = 0; ss <= N; ss++)
        dvecex_sp(nb[ss], 1.0, idxb[ss], dux + ss, 0, dt + ss, 0);
    for (ss = 0; ss <= N; ss++)
        dgemv_t(nu[ss] + nx[ss], ng[ss], 1.0, DCt + ss, 0, 0, dux + ss, 0, 0.0, dt + ss, nb[ss], dt + ss, nb[ss]);

    for (ss = 0; ss <= N; ss++) {
        dveccp(nb[ss] + ng[ss], dt + ss, 0, dt + ss, nb[ss] + ng[ss]);
        dvecsc(nb[ss] + ng[ss], -1.0, dt + ss, nb[ss] + ng[ss]);
    }

    for (ss = 0; ss <= N; ss++) {
        if (ns[ss] > 0)
            EXPAND_SLACKS(ss, qp, qp_sol, ws);
    }

    d_compute_lam_t_qp(res_d[0].pa, res_m[0].pa, dlam[0].pa, dt[0].pa, cws);
}
