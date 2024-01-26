#include <stdio.h>
#include <stdlib.h>

#include "tinyhpipm/common.h"
#include "tinyhpipm/ocp/d_ocp_qcqp.h"
#include "tinyhpipm/ocp/d_ocp_qcqp_dim.h"
#include "tinyhpipm/ocp/d_ocp_qcqp_ipm.h"
#include "tinyhpipm/ocp/d_ocp_qcqp_sol.h"

/************************************************
 * import dims data
 ************************************************/
extern int N;
extern int* nx;
extern int* nu;
extern int* nbu;
extern int* nbx;
extern int* ng;
extern int* nq;
extern int* nsbx;
extern int* nsbu;
extern int* nsg;
extern int* nsq;
extern double** hA;
extern double** hB;
extern double** hb;
extern double** hQ;
extern double** hR;
extern double** hS;
extern double** hq;
extern double** hr;
extern int** hidxbx;
extern double** hlbx;
extern double** hubx;
extern int** hidxbu;
extern double** hlbu;
extern double** hubu;
extern double** hC;
extern double** hD;
extern double** hlg;
extern double** hug;
extern double** hQq;
extern double** hRq;
extern double** hSq;
extern double** hqq;
extern double** hrq;
extern double** huq;
extern double** hZl;
extern double** hZu;
extern double** hzl;
extern double** hzu;
extern int** hidxs;
extern double** hlls;
extern double** hlus;

/************************************************
 * import solver options data
 ************************************************/
extern int mode;
extern int iter_max;
extern double alpha_min;
extern double mu0;
extern double tol_stat;
extern double tol_eq;
extern double tol_ineq;
extern double tol_comp;
extern double reg_prim;
extern int warm_start;
extern int pred_corr;
extern int ric_alg;
extern int split_step;

int main(int argc, char** argv) {
    int ii;  // stage index reused between initialization loops
    int hpipm_status;  // return flag from HPIPM functions
    int rep, nrep = 50, nwarmup = 5;  // repetitions of each timing experiment
                                      // (partial conditioning, solving, expanding)
    struct tinyhpipm_timer timer;  // tic-toc style timer used to time HPIPM functions

    /********************************************************
     * QCQP initialization **********************************
     ********************************************************/
    // create qcqp_dim
    void* dim_mem = malloc(d_ocp_qcqp_dim_memsize(N));
    struct d_ocp_qcqp_dim qcqp_dim;
    d_ocp_qcqp_dim_create(N, &qcqp_dim, dim_mem);
    for (ii = 0; ii <= N; ii++) {
        d_ocp_qcqp_dim_set_nx(ii, nx[ii], &qcqp_dim);
        d_ocp_qcqp_dim_set_nu(ii, nu[ii], &qcqp_dim);
        d_ocp_qcqp_dim_set_nbx(ii, nbx[ii], &qcqp_dim);
        d_ocp_qcqp_dim_set_nbu(ii, nbu[ii], &qcqp_dim);
        d_ocp_qcqp_dim_set_ng(ii, ng[ii], &qcqp_dim);
        d_ocp_qcqp_dim_set_nq(ii, nq[ii], &qcqp_dim);
        d_ocp_qcqp_dim_set_nsbx(ii, nsbx[ii], &qcqp_dim);
        d_ocp_qcqp_dim_set_nsbu(ii, nsbu[ii], &qcqp_dim);
        d_ocp_qcqp_dim_set_nsg(ii, nsg[ii], &qcqp_dim);
        d_ocp_qcqp_dim_set_nsq(ii, nsq[ii], &qcqp_dim);
    }

    // create qcqp
    void* qcqp_mem = malloc(d_ocp_qcqp_memsize(&qcqp_dim));
    struct d_ocp_qcqp qcqp;
    d_ocp_qcqp_create(&qcqp_dim, &qcqp, qcqp_mem);
    for (ii = 0; ii <= N; ii++) {
        // dynamics
        if (ii != N) {
            d_ocp_qcqp_set_A(ii, hA[ii], &qcqp);
            d_ocp_qcqp_set_B(ii, hB[ii], &qcqp);
            d_ocp_qcqp_set_b(ii, hb[ii], &qcqp);
        }
        // cost
        d_ocp_qcqp_set_Q(ii, hQ[ii], &qcqp);
        d_ocp_qcqp_set_S(ii, hS[ii], &qcqp);
        d_ocp_qcqp_set_R(ii, hR[ii], &qcqp);
        d_ocp_qcqp_set_q(ii, hq[ii], &qcqp);
        d_ocp_qcqp_set_r(ii, hr[ii], &qcqp);
        d_ocp_qcqp_set_Zl(ii, hZl[ii], &qcqp);
        d_ocp_qcqp_set_Zu(ii, hZu[ii], &qcqp);
        d_ocp_qcqp_set_zl(ii, hzl[ii], &qcqp);
        d_ocp_qcqp_set_zu(ii, hzu[ii], &qcqp);
        // constraints
        d_ocp_qcqp_set_idxbx(ii, hidxbx[ii], &qcqp);
        d_ocp_qcqp_set_lbx(ii, hlbx[ii], &qcqp);
        d_ocp_qcqp_set_ubx(ii, hubx[ii], &qcqp);
        d_ocp_qcqp_set_idxbu(ii, hidxbu[ii], &qcqp);
        d_ocp_qcqp_set_lbu(ii, hlbu[ii], &qcqp);
        d_ocp_qcqp_set_ubu(ii, hubu[ii], &qcqp);
        d_ocp_qcqp_set_C(ii, hC[ii], &qcqp);
        d_ocp_qcqp_set_D(ii, hD[ii], &qcqp);
        d_ocp_qcqp_set_lg(ii, hlg[ii], &qcqp);
        d_ocp_qcqp_set_ug(ii, hug[ii], &qcqp);
        d_ocp_qcqp_set_Rq(ii, hRq[ii], &qcqp);
        d_ocp_qcqp_set_Sq(ii, hSq[ii], &qcqp);
        d_ocp_qcqp_set_Qq(ii, hQq[ii], &qcqp);
        d_ocp_qcqp_set_rq(ii, hrq[ii], &qcqp);
        d_ocp_qcqp_set_qq(ii, hqq[ii], &qcqp);
        d_ocp_qcqp_set_uq(ii, huq[ii], &qcqp);
        d_ocp_qcqp_set_idxs(ii, hidxs[ii], &qcqp);
        d_ocp_qcqp_set_lls(ii, hlls[ii], &qcqp);
        d_ocp_qcqp_set_lus(ii, hlus[ii], &qcqp);
    }

    // create qcqp_sol
    void* qcqp_sol_mem = malloc(d_ocp_qcqp_sol_memsize(&qcqp_dim));
    struct d_ocp_qcqp_sol qcqp_sol;
    d_ocp_qcqp_sol_create(&qcqp_dim, &qcqp_sol, qcqp_sol_mem);

    // create qcqp_ipm_arg
    void* ipm_arg_mem = malloc(d_ocp_qcqp_ipm_arg_memsize(&qcqp_dim));
    struct d_ocp_qcqp_ipm_arg qcqp_ipm_arg;
    d_ocp_qcqp_ipm_arg_create(&qcqp_dim, &qcqp_ipm_arg, ipm_arg_mem);

    d_ocp_qcqp_ipm_arg_set_default(mode, &qcqp_ipm_arg);

    d_ocp_qcqp_ipm_arg_set_mu0(&mu0, &qcqp_ipm_arg);
    d_ocp_qcqp_ipm_arg_set_iter_max(&iter_max, &qcqp_ipm_arg);
    d_ocp_qcqp_ipm_arg_set_tol_stat(&tol_stat, &qcqp_ipm_arg);
    d_ocp_qcqp_ipm_arg_set_tol_eq(&tol_eq, &qcqp_ipm_arg);
    d_ocp_qcqp_ipm_arg_set_tol_ineq(&tol_ineq, &qcqp_ipm_arg);
    d_ocp_qcqp_ipm_arg_set_tol_comp(&tol_comp, &qcqp_ipm_arg);
    d_ocp_qcqp_ipm_arg_set_reg_prim(&reg_prim, &qcqp_ipm_arg);
    d_ocp_qcqp_ipm_arg_set_warm_start(&warm_start, &qcqp_ipm_arg);
    d_ocp_qcqp_ipm_arg_set_pred_corr(&pred_corr, &qcqp_ipm_arg);
    d_ocp_qcqp_ipm_arg_set_ric_alg(&ric_alg, &qcqp_ipm_arg);
    d_ocp_qcqp_ipm_arg_set_split_step(&split_step, &qcqp_ipm_arg);

    /*****************************************************************
     * Solve QCQP ****************************************************
     *****************************************************************/

    // create qcqp_ipm_ws (based on the partially condensed problem)
    hpipm_size_t ipm_size = d_ocp_qcqp_ipm_ws_memsize(&qcqp_dim, &qcqp_ipm_arg);
    void* ipm_mem = malloc(ipm_size);
    struct d_ocp_qcqp_ipm_ws qcqp_ipm_ws;
    d_ocp_qcqp_ipm_ws_create(&qcqp_dim, &qcqp_ipm_arg, &qcqp_ipm_ws, ipm_mem);

    // perform solving
    hpipm_timer timer2;

    double xinit[16] = {0., 0., 0.28, 0., 0., 0., 0., 0.,
                        -0.001, 0., 0., 0., 0.0, 0.0, 0.0, 0.0};
    double uinit[4] = {0.0, 0.0, 0.0, 0.0};
    for (rep = 0; rep < nrep; rep++) {
        if (rep == nwarmup) {
            hpipm_tic(&timer);
        }
        // reset initial guess
        for (int ii = 0; ii < N; ii++) {
            d_ocp_qcqp_sol_set_x(ii, xinit, &qcqp_sol);
            d_ocp_qcqp_sol_set_u(ii, uinit, &qcqp_sol);
        }
        d_ocp_qcqp_sol_set_x(N, xinit, &qcqp_sol);
        hpipm_tic(&timer2);
        d_ocp_qcqp_ipm_solve(&qcqp, &qcqp_sol, &qcqp_ipm_arg, &qcqp_ipm_ws);
        printf("ipm time = %.3f ms\n", 1000 * hpipm_toc(&timer2));
        d_ocp_qcqp_ipm_get_status(&qcqp_ipm_ws, &hpipm_status);
    }
    double time_ipm = hpipm_toc(&timer) / nrep;

    /*****************************************************************
     * Print solution and info
     *****************************************************************/

#ifndef NDEBUG
    printf("\nHPIPM returned with flag %i ", hpipm_status);
#endif
    if (hpipm_status == 0) {
#ifndef NDEBUG
        printf("-> QP solved!\n");

        // print solution
        // u
        int nu_max = nu[0];
        for (ii = 1; ii <= N; ii++) {
            if (nu[ii] > nu_max) {
                nu_max = nu[ii];
            }
        }
        double* u = malloc(nu_max * sizeof(double));
        printf("\nu = \n");
        for (ii = 0; ii <= N; ii++) {
            d_ocp_qcqp_sol_get_u(ii, &qcqp_sol, u);
            d_print_mat(1, nu[ii], u, 1);
        }
        // x
        int nx_max = nx[0];
        for (ii = 1; ii <= N; ii++) {
            if (nx[ii] > nx_max) {
                nx_max = nx[ii];
            }
        }
        double* x = malloc(nx_max * sizeof(double));
        printf("\nx = \n");
        for (ii = 0; ii <= N; ii++) {
            d_ocp_qcqp_sol_get_x(ii, &qcqp_sol, x);
            d_print_mat(1, nx[ii], x, 1);
        }
        free(u);
        free(x);
#endif
    } else {
        if (hpipm_status == 1) {
            printf("-> Solver failed! Maximum number of iterations reached\n");
        } else if (hpipm_status == 2) {
            printf("-> Solver failed! Minimum step lenght reached\n");
        } else if (hpipm_status == 3) {
            printf("-> Solver failed! NaN in computations\n");
        } else if (hpipm_status == 4) {
            printf("-> Solver failed! Inconsistent inequality constraints\n");
        } else {
            printf("-> Solver failed! Unknown return flag\n");
        }
        // print ipm statistics
        int iter;
        d_ocp_qcqp_ipm_get_iter(&qcqp_ipm_ws, &iter);
        double res_stat;
        d_ocp_qcqp_ipm_get_max_res_stat(&qcqp_ipm_ws, &res_stat);
        double res_eq;
        d_ocp_qcqp_ipm_get_max_res_eq(&qcqp_ipm_ws, &res_eq);
        double res_ineq;
        d_ocp_qcqp_ipm_get_max_res_ineq(&qcqp_ipm_ws, &res_ineq);
        double res_comp;
        d_ocp_qcqp_ipm_get_max_res_comp(&qcqp_ipm_ws, &res_comp);
        double* stat;
        d_ocp_qcqp_ipm_get_stat(&qcqp_ipm_ws, &stat);
        int stat_m;
        d_ocp_qcqp_ipm_get_stat_m(&qcqp_ipm_ws, &stat_m);

        printf(
                "\nipm residuals max: res_g = %e, res_b = %e, res_d = %e, res_m = %e\n",
                res_stat, res_eq, res_ineq, res_comp);
        printf("ipm iter = %d\n", iter);
        printf("\niteration stats:"
               "\nalpha_aff\tmu_aff\tsigma\talpha_prim\talpha_dual\tmu\t"
               "res_stat\tres_eq\tres_ineq\tres_comp\tobj\tlq_fact\t"
               "itref_pred\titref_corr\tlin_res_stat\tlin_res_eq\t"
               "lin_res_ineq\tlin_res_comp\n");
        d_print_exp_tran_mat(stat_m, iter + 1, stat, stat_m);
    }
    // print timing info
    printf("\nAverage solution time over %i runs: %.3f ms\n", nrep,
           1000 * time_ipm);
    printf("ocp ipm time = %.3f ms\n", 1000 * time_ipm);

    // free memory
    free(dim_mem);
    free(qcqp_mem);
    free(qcqp_sol_mem);
    free(ipm_arg_mem);
    free(ipm_mem);

    return 0;
}
