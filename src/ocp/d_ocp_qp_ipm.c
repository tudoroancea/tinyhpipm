#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "tinyhpipm/blas.h"
#include "tinyhpipm/common.h"
#include "tinyhpipm/ipm_core/d_core_qp_ipm.h"
#include "tinyhpipm/ipm_core/d_core_qp_ipm_aux.h"
#include "tinyhpipm/ocp/d_ocp_qp.h"
#include "tinyhpipm/ocp/d_ocp_qp_dim.h"
#include "tinyhpipm/ocp/d_ocp_qp_ipm.h"
#include "tinyhpipm/ocp/d_ocp_qp_kkt.h"
#include "tinyhpipm/ocp/d_ocp_qp_res.h"
#include "tinyhpipm/ocp/d_ocp_qp_sol.h"
#include "tinyhpipm/ocp/d_ocp_qp_utils.h"


hpipm_size_t d_ocp_qp_ipm_arg_strsize() {
    return sizeof(struct d_ocp_qp_ipm_arg);
}


hpipm_size_t d_ocp_qp_ipm_arg_memsize(struct d_ocp_qp_dim* dim) {

    return 0;
}


void d_ocp_qp_ipm_arg_create(struct d_ocp_qp_dim* dim, struct d_ocp_qp_ipm_arg* arg, void* mem) {

    // zero memory (to avoid corrupted memory like e.g. NaN)
    //	hpipm_size_t memsize = d_ocp_qp_ipm_arg_memsize(dim);
    //	hpipm_zero_memset(memsize, mem);

    arg->memsize = 0;
}


void d_ocp_qp_ipm_arg_set_default(enum tinyhpipm_mode mode, struct d_ocp_qp_ipm_arg* arg) {

    double mu0, alpha_min, res_g_max, res_b_max, res_d_max, res_m_max, reg_prim, lam_min, t_min, tau_min;
    int iter_max, stat_max, pred_corr, cond_pred_corr, itref_pred_max, itref_corr_max, lq_fact, warm_start, abs_form, comp_res_exit, comp_res_pred, square_root_alg, comp_dual_sol_eq, split_step, var_init_scheme, t_lam_min;

    if (mode == SPEED_ABS) {
        mu0 = 1e1;
        alpha_min = 1e-12;
        res_g_max = 1e0;  // not used
        res_b_max = 1e0;  // not used
        res_d_max = 1e0;  // not used
        res_m_max = 1e-8;
        iter_max = 15;
        stat_max = 15;
        pred_corr = 1;
        cond_pred_corr = 1;
        itref_pred_max = 0;  // not used
        itref_corr_max = 0;  // not used
        reg_prim = 1e-15;
        square_root_alg = 1;
        lq_fact = 0;  // not used
        lam_min = 1e-16;
        t_min = 1e-16;
        tau_min = 1e-16;
        warm_start = 0;
        abs_form = 1;
        comp_dual_sol_eq = 0;
        comp_res_exit = 0;
        comp_res_pred = 0;
        split_step = 1;
        var_init_scheme = 0;
        t_lam_min = 2;
    } else if (mode == SPEED) {
        mu0 = 1e1;
        alpha_min = 1e-12;
        res_g_max = 1e-6;
        res_b_max = 1e-8;
        res_d_max = 1e-8;
        res_m_max = 1e-8;
        iter_max = 15;
        stat_max = 15;
        pred_corr = 1;
        cond_pred_corr = 1;
        itref_pred_max = 0;
        itref_corr_max = 0;
        reg_prim = 1e-15;
        square_root_alg = 1;
        lq_fact = 0;
        lam_min = 1e-16;
        t_min = 1e-16;
        tau_min = 1e-16;
        warm_start = 0;
        abs_form = 0;
        comp_dual_sol_eq = 1;
        comp_res_exit = 1;
        comp_res_pred = 1;
        split_step = 1;
        var_init_scheme = 0;
        t_lam_min = 2;
    } else if (mode == BALANCE) {
        mu0 = 1e1;
        alpha_min = 1e-12;
        res_g_max = 1e-6;
        res_b_max = 1e-8;
        res_d_max = 1e-8;
        res_m_max = 1e-8;
        iter_max = 30;
        stat_max = 30;
        pred_corr = 1;
        cond_pred_corr = 1;
        itref_pred_max = 0;
        itref_corr_max = 2;
        reg_prim = 1e-15;
        square_root_alg = 1;
        lq_fact = 1;
        lam_min = 1e-16;
        t_min = 1e-16;
        tau_min = 1e-16;
        warm_start = 0;
        abs_form = 0;
        comp_dual_sol_eq = 1;
        comp_res_exit = 1;
        comp_res_pred = 1;
        split_step = 0;
        var_init_scheme = 0;
        t_lam_min = 2;
    } else if (mode == ROBUST) {
        mu0 = 1e2;
        alpha_min = 1e-12;
        res_g_max = 1e-6;
        res_b_max = 1e-8;
        res_d_max = 1e-8;
        res_m_max = 1e-8;
        iter_max = 100;
        stat_max = 100;
        pred_corr = 1;
        cond_pred_corr = 1;
        itref_pred_max = 0;
        itref_corr_max = 4;
        reg_prim = 1e-15;
        square_root_alg = 1;
        lq_fact = 2;
        lam_min = 1e-16;
        t_min = 1e-16;
        tau_min = 1e-16;
        warm_start = 0;
        abs_form = 0;
        comp_dual_sol_eq = 1;
        comp_res_exit = 1;
        comp_res_pred = 1;
        split_step = 0;
        var_init_scheme = 0;
        t_lam_min = 2;
    } else {
        printf("\nerror: d_ocp_qp_ipm_arg_set_default: wrong set default mode\n");
        exit(1);
    }

    // use individual setters when available
    d_ocp_qp_ipm_arg_set_mu0(&mu0, arg);
    d_ocp_qp_ipm_arg_set_alpha_min(&alpha_min, arg);
    d_ocp_qp_ipm_arg_set_tol_stat(&res_g_max, arg);
    d_ocp_qp_ipm_arg_set_tol_eq(&res_b_max, arg);
    d_ocp_qp_ipm_arg_set_tol_ineq(&res_d_max, arg);
    d_ocp_qp_ipm_arg_set_tol_comp(&res_m_max, arg);
    d_ocp_qp_ipm_arg_set_iter_max(&iter_max, arg);
    arg->stat_max = stat_max;
    d_ocp_qp_ipm_arg_set_pred_corr(&pred_corr, arg);
    d_ocp_qp_ipm_arg_set_cond_pred_corr(&cond_pred_corr, arg);
    d_ocp_qp_ipm_arg_set_ric_alg(&square_root_alg, arg);
    arg->itref_pred_max = itref_pred_max;
    arg->itref_corr_max = itref_corr_max;
    d_ocp_qp_ipm_arg_set_reg_prim(&reg_prim, arg);
    arg->lq_fact = lq_fact;
    d_ocp_qp_ipm_arg_set_lam_min(&lam_min, arg);
    d_ocp_qp_ipm_arg_set_t_min(&t_min, arg);
    d_ocp_qp_ipm_arg_set_tau_min(&tau_min, arg);
    d_ocp_qp_ipm_arg_set_warm_start(&warm_start, arg);
    arg->abs_form = abs_form;
    d_ocp_qp_ipm_arg_set_comp_dual_sol_eq(&comp_dual_sol_eq, arg);
    d_ocp_qp_ipm_arg_set_comp_res_pred(&comp_res_pred, arg);
    d_ocp_qp_ipm_arg_set_comp_res_exit(&comp_res_pred, arg);
    d_ocp_qp_ipm_arg_set_split_step(&split_step, arg);
    d_ocp_qp_ipm_arg_set_var_init_scheme(&var_init_scheme, arg);
    d_ocp_qp_ipm_arg_set_t_lam_min(&t_lam_min, arg);
    arg->mode = mode;
}


void d_ocp_qp_ipm_arg_set(char* field, void* value, struct d_ocp_qp_ipm_arg* arg) {
    if (hpipm_strcmp(field, "iter_max")) {
        d_ocp_qp_ipm_arg_set_iter_max(value, arg);
    } else if (hpipm_strcmp(field, "alpha_min")) {
        d_ocp_qp_ipm_arg_set_alpha_min(value, arg);
    } else if (hpipm_strcmp(field, "mu0")) {
        d_ocp_qp_ipm_arg_set_mu0(value, arg);
    } else if (hpipm_strcmp(field, "tol_stat")) {
        d_ocp_qp_ipm_arg_set_tol_stat(value, arg);
    } else if (hpipm_strcmp(field, "tol_eq")) {
        d_ocp_qp_ipm_arg_set_tol_eq(value, arg);
    } else if (hpipm_strcmp(field, "tol_ineq")) {
        d_ocp_qp_ipm_arg_set_tol_ineq(value, arg);
    } else if (hpipm_strcmp(field, "tol_comp")) {
        d_ocp_qp_ipm_arg_set_tol_comp(value, arg);
    } else if (hpipm_strcmp(field, "reg_prim")) {
        d_ocp_qp_ipm_arg_set_reg_prim(value, arg);
    } else if (hpipm_strcmp(field, "warm_start")) {
        d_ocp_qp_ipm_arg_set_warm_start(value, arg);
    } else if (hpipm_strcmp(field, "pred_corr")) {
        d_ocp_qp_ipm_arg_set_pred_corr(value, arg);
    } else if (hpipm_strcmp(field, "cond_pred_corr")) {
        d_ocp_qp_ipm_arg_set_cond_pred_corr(value, arg);
    } else if (hpipm_strcmp(field, "ric_alg")) {
        d_ocp_qp_ipm_arg_set_ric_alg(value, arg);
    } else if (hpipm_strcmp(field, "comp_dual_sol_eq")) {
        d_ocp_qp_ipm_arg_set_comp_dual_sol_eq(value, arg);
    } else if (hpipm_strcmp(field, "comp_res_exit")) {
        d_ocp_qp_ipm_arg_set_comp_res_exit(value, arg);
    } else if (hpipm_strcmp(field, "comp_res_pred")) {
        d_ocp_qp_ipm_arg_set_comp_res_pred(value, arg);
    } else if (hpipm_strcmp(field, "lam_min")) {
        d_ocp_qp_ipm_arg_set_lam_min(value, arg);
    } else if (hpipm_strcmp(field, "t_min")) {
        d_ocp_qp_ipm_arg_set_t_min(value, arg);
    } else if (hpipm_strcmp(field, "tau_min")) {
        d_ocp_qp_ipm_arg_set_tau_min(value, arg);
    } else if (hpipm_strcmp(field, "split_step")) {
        d_ocp_qp_ipm_arg_set_split_step(value, arg);
    } else if (hpipm_strcmp(field, "var_init_scheme")) {
        d_ocp_qp_ipm_arg_set_var_init_scheme(value, arg);
    } else if (hpipm_strcmp(field, "t_lam_min")) {
        d_ocp_qp_ipm_arg_set_t_lam_min(value, arg);
    } else {
        printf("error: d_ocp_qp_ipm_arg_set: wrong field %s\n", field);
        exit(1);
    }
}


void d_ocp_qp_ipm_arg_set_iter_max(int* iter_max, struct d_ocp_qp_ipm_arg* arg) {
    arg->iter_max = *iter_max;
}


void d_ocp_qp_ipm_arg_set_alpha_min(double* alpha_min, struct d_ocp_qp_ipm_arg* arg) {
    arg->alpha_min = *alpha_min;
}


void d_ocp_qp_ipm_arg_set_mu0(double* mu0, struct d_ocp_qp_ipm_arg* arg) {
    arg->mu0 = *mu0;
}


void d_ocp_qp_ipm_arg_set_tol_stat(double* tol_stat, struct d_ocp_qp_ipm_arg* arg) {
    arg->res_g_max = *tol_stat;
}


void d_ocp_qp_ipm_arg_set_tol_eq(double* tol_eq, struct d_ocp_qp_ipm_arg* arg) {
    arg->res_b_max = *tol_eq;
}


void d_ocp_qp_ipm_arg_set_tol_ineq(double* tol_ineq, struct d_ocp_qp_ipm_arg* arg) {
    arg->res_d_max = *tol_ineq;
}


void d_ocp_qp_ipm_arg_set_tol_comp(double* tol_comp, struct d_ocp_qp_ipm_arg* arg) {
    arg->res_m_max = *tol_comp;
}


void d_ocp_qp_ipm_arg_set_reg_prim(double* reg, struct d_ocp_qp_ipm_arg* arg) {
    arg->reg_prim = *reg;
}


void d_ocp_qp_ipm_arg_set_warm_start(int* warm_start, struct d_ocp_qp_ipm_arg* arg) {
    arg->warm_start = *warm_start;
}


void d_ocp_qp_ipm_arg_set_pred_corr(int* pred_corr, struct d_ocp_qp_ipm_arg* arg) {
    arg->pred_corr = *pred_corr;
}


void d_ocp_qp_ipm_arg_set_cond_pred_corr(int* value, struct d_ocp_qp_ipm_arg* arg) {
    arg->cond_pred_corr = *value;
}


void d_ocp_qp_ipm_arg_set_ric_alg(int* ric_alg, struct d_ocp_qp_ipm_arg* arg) {
    arg->square_root_alg = *ric_alg;
}


void d_ocp_qp_ipm_arg_set_comp_dual_sol_eq(int* value, struct d_ocp_qp_ipm_arg* arg) {
    arg->comp_dual_sol_eq = *value;
}


void d_ocp_qp_ipm_arg_set_comp_res_exit(int* comp_res_exit, struct d_ocp_qp_ipm_arg* arg) {
    arg->comp_res_exit = *comp_res_exit;
    if (*comp_res_exit != 0)
        arg->comp_dual_sol_eq = 1;
}


void d_ocp_qp_ipm_arg_set_comp_res_pred(int* comp_res_pred, struct d_ocp_qp_ipm_arg* arg) {
    arg->comp_res_pred = *comp_res_pred;
}


void d_ocp_qp_ipm_arg_set_lam_min(double* value, struct d_ocp_qp_ipm_arg* arg) {
    arg->lam_min = *value;
}


void d_ocp_qp_ipm_arg_set_t_min(double* value, struct d_ocp_qp_ipm_arg* arg) {
    arg->t_min = *value;
}


void d_ocp_qp_ipm_arg_set_tau_min(double* value, struct d_ocp_qp_ipm_arg* arg) {
    arg->tau_min = *value;
}


void d_ocp_qp_ipm_arg_set_split_step(int* value, struct d_ocp_qp_ipm_arg* arg) {
    arg->split_step = *value;
}


void d_ocp_qp_ipm_arg_set_var_init_scheme(int* value, struct d_ocp_qp_ipm_arg* arg) {
    arg->var_init_scheme = *value;
}


void d_ocp_qp_ipm_arg_set_t_lam_min(int* value, struct d_ocp_qp_ipm_arg* arg) {
    arg->t_lam_min = *value;
}


hpipm_size_t d_ocp_qp_ipm_ws_strsize() {
    return sizeof(struct d_ocp_qp_ipm_ws);
}


hpipm_size_t d_ocp_qp_ipm_ws_memsize(struct d_ocp_qp_dim* dim, struct d_ocp_qp_ipm_arg* arg) {

    // stat_max is at least as big as iter_max
    if (arg->iter_max > arg->stat_max)
        arg->stat_max = arg->iter_max;

    // loop index
    int ii;

    // extract ocp qp size
    int N = dim->N;
    int* nx = dim->nx;
    int* nu = dim->nu;
    int* nb = dim->nb;
    int* ng = dim->ng;
    int* ns = dim->ns;

    // compute core qp size and max size
    int nvt = 0;
    int net = 0;
    int nct = 0;
    int nxM = 0;
    int nuM = 0;
    int nbM = 0;
    int ngM = 0;
    int nsM = 0;
    for (ii = 0; ii < N; ii++) {
        nvt += nx[ii] + nu[ii] + 2 * ns[ii];
        net += nx[ii + 1];
        nct += 2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii];
        nxM = nx[ii] > nxM ? nx[ii] : nxM;
        nuM = nu[ii] > nuM ? nu[ii] : nuM;
        nbM = nb[ii] > nbM ? nb[ii] : nbM;
        ngM = ng[ii] > ngM ? ng[ii] : ngM;
        nsM = ns[ii] > nsM ? ns[ii] : nsM;
    }
    nvt += nx[ii] + nu[ii] + 2 * ns[ii];
    nct += 2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii];
    nxM = nx[ii] > nxM ? nx[ii] : nxM;
    nuM = nu[ii] > nuM ? nu[ii] : nuM;
    nbM = nb[ii] > nbM ? nb[ii] : nbM;
    ngM = ng[ii] > ngM ? ng[ii] : ngM;
    nsM = ns[ii] > nsM ? ns[ii] : nsM;

    hpipm_size_t size = 0;

    size += 1 * sizeof(struct d_core_qp_ipm_workspace);
    size += 1 * d_memsize_core_qp_ipm(nvt, net, nct);

    size += 1 * sizeof(struct d_ocp_qp_res_ws);  // res_workspace

    size += 2 * sizeof(struct d_ocp_qp);  // qp_step qp_itref

    size += 2 * sizeof(struct d_ocp_qp_sol);  // sol_step sol_itref
    size += 1 * d_ocp_qp_sol_memsize(dim);  // sol_itref

    size += 2 * sizeof(struct d_ocp_qp_res);  // res res_itref
    size += 1 * d_ocp_qp_res_memsize(dim);  // res_itref

    size += 10 * (N + 1) * sizeof(struct vec);  // res_g res_d res_m Gamma gamma Zs_inv sol_step(v,lam,t)  tmp_m
    size += 3 * N * sizeof(struct vec);  // res_b Pb sol_step(pi)
    size += 9 * sizeof(struct vec);  // tmp_nuxM (4+2)*tmp_nbgM (1+1)*tmp_nsM
    size += 1 * (N + 1) * sizeof(struct vec);  // l

    size += 1 * (N + 1) * sizeof(struct mat);  // L
    if (!arg->square_root_alg) {
        size += 1 * (N + 1) * sizeof(struct mat);  // P
    }
    size += 3 * sizeof(struct mat);  // tmp_nxM_nxM Ls
    if (arg->lq_fact > 0) {
        size += 1 * (N + 1) * sizeof(struct mat);  // Lh
    }
    size += 2 * sizeof(struct mat);  // AL
    if (arg->lq_fact > 0) {
        size += 1 * sizeof(struct mat);  // lq0
    }

    size += 1 * memsize_vec(nuM + nxM);  // tmp_nuxM
    size += 4 * memsize_vec(nbM + ngM);  // tmp_nbgM
    size += 1 * memsize_vec(nsM);  // tmp_nsM
    for (ii = 0; ii < N; ii++) size += 1 * memsize_vec(nx[ii + 1]);  // Pb
    for (ii = 0; ii <= N; ii++) size += 1 * memsize_vec(2 * ns[ii]);  // Zs_inv
    for (ii = 0; ii <= N; ii++) size += 1 * memsize_vec(nu[ii] + nx[ii]);  // l

    for (ii = 0; ii <= N; ii++) size += 1 * memsize_mat(nu[ii] + nx[ii] + 1, nu[ii] + nx[ii]);  // L
    if (!arg->square_root_alg) {
        for (ii = 0; ii <= N; ii++) size += 1 * memsize_mat(nx[ii] + 1, nx[ii]);  // P
    }
    size += 2 * memsize_mat(nxM, nxM);  // tmp_nxM_nxM
    size += 1 * memsize_mat(nxM + 1, nuM);  // Ls
    if (arg->lq_fact > 0) {
        for (ii = 0; ii <= N; ii++) size += 1 * memsize_mat(nu[ii] + nx[ii] + 1, nu[ii] + nx[ii]);  // Lh
    }
    size += 2 * memsize_mat(nuM + nxM + 1, nxM + ngM);  // AL
    if (arg->lq_fact > 0) {
        size += 1 * memsize_mat(nuM + nxM, 2 * nuM + 3 * nxM + ngM);  // lq0
    }
    size += 1 * memsize_vec(nct);  // tmp_m

    if (arg->lq_fact > 0) {
        size += 1 * dgelqf_worksize(nuM + nxM, 2 * nuM + 3 * nxM + ngM);  // lq_work0
    }

    int stat_m = 18;
    size += stat_m * (1 + arg->stat_max) * sizeof(double);  // stat

    size += (N + 1) * sizeof(int);  // use_hess_fact

    size += 1 * 64;  // align once to typical cache line size
    size += 1 * 8;  // align once to 8-byte boundary
    size = (size + 63) / 64 * 64;  // make multiple of typical cache line size

    return size;
}


void d_ocp_qp_ipm_ws_create(struct d_ocp_qp_dim* dim, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* workspace, void* mem) {

    // loop index
    int ii;

    // zero memory (to avoid corrupted memory like e.g. NaN)
    hpipm_size_t memsize = d_ocp_qp_ipm_ws_memsize(dim, arg);
    hpipm_zero_memset(memsize, mem);

    // extract ocp qp size
    int N = dim->N;
    int* nx = dim->nx;
    int* nu = dim->nu;
    int* nb = dim->nb;
    int* ng = dim->ng;
    int* ns = dim->ns;


    // compute core qp size and max size
    int nvt = 0;
    int net = 0;
    int nct = 0;
    int nxM = 0;
    int nuM = 0;
    int nbM = 0;
    int ngM = 0;
    int nsM = 0;
    for (ii = 0; ii < N; ii++) {
        nvt += nx[ii] + nu[ii] + 2 * ns[ii];
        net += nx[ii + 1];
        nct += 2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii];
        nxM = nx[ii] > nxM ? nx[ii] : nxM;
        nuM = nu[ii] > nuM ? nu[ii] : nuM;
        nbM = nb[ii] > nbM ? nb[ii] : nbM;
        ngM = ng[ii] > ngM ? ng[ii] : ngM;
        nsM = ns[ii] > nsM ? ns[ii] : nsM;
    }
    nvt += nx[ii] + nu[ii] + 2 * ns[ii];
    nct += 2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii];
    nxM = nx[ii] > nxM ? nx[ii] : nxM;
    nuM = nu[ii] > nuM ? nu[ii] : nuM;
    nbM = nb[ii] > nbM ? nb[ii] : nbM;
    ngM = ng[ii] > ngM ? ng[ii] : ngM;
    nsM = ns[ii] > nsM ? ns[ii] : nsM;


    // core struct
    struct d_core_qp_ipm_workspace* sr_ptr = mem;

    // core workspace
    workspace->core_workspace = sr_ptr;
    sr_ptr += 1;
    struct d_core_qp_ipm_workspace* cws = workspace->core_workspace;


    // res struct
    struct d_ocp_qp_res* res_ptr = (struct d_ocp_qp_res*) sr_ptr;
    workspace->res = res_ptr;
    res_ptr += 1;
    workspace->res_itref = res_ptr;
    res_ptr += 1;


    // res workspace struct
    struct d_ocp_qp_res_ws* res_ws_ptr = (struct d_ocp_qp_res_ws*) res_ptr;
    workspace->res_workspace = res_ws_ptr;
    res_ws_ptr += 1;


    // qp sol struct
    struct d_ocp_qp_sol* qp_sol_ptr = (struct d_ocp_qp_sol*) res_ws_ptr;

    workspace->sol_step = qp_sol_ptr;
    qp_sol_ptr += 1;
    workspace->sol_itref = qp_sol_ptr;
    qp_sol_ptr += 1;


    // qp struct
    struct d_ocp_qp* qp_ptr = (struct d_ocp_qp*) qp_sol_ptr;

    workspace->qp_step = qp_ptr;
    qp_ptr += 1;
    workspace->qp_itref = qp_ptr;
    qp_ptr += 1;


    // matrix struct
    struct mat* sm_ptr = (struct mat*) qp_ptr;

    workspace->L = sm_ptr;
    sm_ptr += N + 1;
    if (!arg->square_root_alg) {
        workspace->P = sm_ptr;
        sm_ptr += N + 1;
    }
    workspace->tmp_nxM_nxM = sm_ptr;
    sm_ptr += 2;
    workspace->Ls = sm_ptr;
    sm_ptr += 1;
    if (arg->lq_fact > 0) {
        workspace->Lh = sm_ptr;
        sm_ptr += N + 1;
    }
    workspace->AL = sm_ptr;
    sm_ptr += 2;
    if (arg->lq_fact > 0) {
        workspace->lq0 = sm_ptr;
        sm_ptr += 1;
    }


    // vector struct
    struct vec* sv_ptr = (struct vec*) sm_ptr;

    workspace->l = sv_ptr;
    sv_ptr += N + 1;
    workspace->sol_step->ux = sv_ptr;
    sv_ptr += N + 1;
    workspace->sol_step->pi = sv_ptr;
    sv_ptr += N;
    workspace->sol_step->lam = sv_ptr;
    sv_ptr += N + 1;
    workspace->sol_step->t = sv_ptr;
    sv_ptr += N + 1;
    workspace->res->res_g = sv_ptr;
    sv_ptr += N + 1;
    workspace->res->res_b = sv_ptr;
    sv_ptr += N;
    workspace->res->res_d = sv_ptr;
    sv_ptr += N + 1;
    workspace->res->res_m = sv_ptr;
    sv_ptr += N + 1;
    workspace->Gamma = sv_ptr;
    sv_ptr += N + 1;
    workspace->gamma = sv_ptr;
    sv_ptr += N + 1;
    workspace->Pb = sv_ptr;
    sv_ptr += N;
    workspace->Zs_inv = sv_ptr;
    sv_ptr += N + 1;
    workspace->tmp_nuxM = sv_ptr;
    sv_ptr += 1;
    workspace->tmp_nbgM = sv_ptr;
    sv_ptr += 4;
    workspace->res_workspace->tmp_nbgM = sv_ptr;
    sv_ptr += 2;
    workspace->tmp_nsM = sv_ptr;
    sv_ptr += 1;
    workspace->res_workspace->tmp_nsM = sv_ptr;
    sv_ptr += 1;
    workspace->tmp_m = sv_ptr;
    sv_ptr += N + 1;


    // align to 8-byte boundary
    hpipm_size_t s_ptr = (hpipm_size_t) sv_ptr;
    s_ptr = (s_ptr + 7) / 8 * 8;


    // double/float stuff
    double* d_ptr = (double*) s_ptr;

    workspace->stat = d_ptr;
    int stat_m = 18;
    d_ptr += stat_m * (1 + arg->stat_max);


    // int stuff
    int* i_ptr = (int*) d_ptr;

    workspace->use_hess_fact = i_ptr;
    i_ptr += N + 1;


    // align to typical cache line size
    s_ptr = (hpipm_size_t) i_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;


    // void stuf
    char* c_ptr = (char*) s_ptr;

    d_ocp_qp_sol_create(dim, workspace->sol_itref, c_ptr);
    c_ptr += workspace->sol_itref->memsize;

    d_ocp_qp_res_create(dim, workspace->res_itref, c_ptr);
    c_ptr += workspace->res_itref->memsize;

    for (ii = 0; ii <= N; ii++) {
        create_mat(nu[ii] + nx[ii] + 1, nu[ii] + nx[ii], workspace->L + ii, c_ptr);
        c_ptr += (workspace->L + ii)->memsize;
    }
    if (!arg->square_root_alg) {
        for (ii = 0; ii <= N; ii++) {
            create_mat(nx[ii] + 1, nx[ii], workspace->P + ii, c_ptr);
            c_ptr += (workspace->P + ii)->memsize;
        }
    }
    create_mat(nxM, nxM, workspace->tmp_nxM_nxM + 0, c_ptr);
    c_ptr += (workspace->tmp_nxM_nxM + 0)->memsize;
    create_mat(nxM, nxM, workspace->tmp_nxM_nxM + 1, c_ptr);
    c_ptr += (workspace->tmp_nxM_nxM + 1)->memsize;
    create_mat(nxM + 1, nuM, workspace->Ls, c_ptr);
    c_ptr += (workspace->Ls)->memsize;

    if (arg->lq_fact > 0) {
        for (ii = 0; ii <= N; ii++) {
            create_mat(nu[ii] + nx[ii] + 1, nu[ii] + nx[ii], workspace->Lh + ii, c_ptr);
            c_ptr += (workspace->Lh + ii)->memsize;
        }
    }

    create_mat(nuM + nxM + 1, nxM + ngM, workspace->AL + 0, c_ptr);
    c_ptr += (workspace->AL + 0)->memsize;

    create_mat(nuM + nxM + 1, nxM + ngM, workspace->AL + 1, c_ptr);
    c_ptr += (workspace->AL + 1)->memsize;

    if (arg->lq_fact > 0) {
        create_mat(nuM + nxM, 2 * nuM + 3 * nxM + ngM, workspace->lq0, c_ptr);
        c_ptr += (workspace->lq0)->memsize;
    }

    for (ii = 0; ii <= N; ii++) {
        create_vec(nu[ii] + nx[ii], workspace->l + ii, c_ptr);
        c_ptr += (workspace->l + ii)->memsize;
    }

    for (ii = 0; ii < N; ii++) {
        create_vec(nx[ii + 1], workspace->Pb + ii, c_ptr);
        c_ptr += (workspace->Pb + ii)->memsize;
    }

    for (ii = 0; ii < N + 1; ii++) {
        create_vec(2 * ns[ii], workspace->Zs_inv + ii, c_ptr);
        c_ptr += (workspace->Zs_inv + ii)->memsize;
    }

    create_vec(nuM + nxM, workspace->tmp_nuxM, c_ptr);
    c_ptr += workspace->tmp_nuxM->memsize;

    create_vec(nbM + ngM, workspace->tmp_nbgM + 0, c_ptr);
    create_vec(nbM + ngM, workspace->res_workspace->tmp_nbgM + 0, c_ptr);
    c_ptr += (workspace->tmp_nbgM + 0)->memsize;

    create_vec(nbM + ngM, workspace->tmp_nbgM + 1, c_ptr);
    create_vec(nbM + ngM, workspace->res_workspace->tmp_nbgM + 1, c_ptr);
    c_ptr += (workspace->tmp_nbgM + 1)->memsize;

    create_vec(nbM + ngM, workspace->tmp_nbgM + 2, c_ptr);
    c_ptr += (workspace->tmp_nbgM + 2)->memsize;

    create_vec(nbM + ngM, workspace->tmp_nbgM + 3, c_ptr);
    c_ptr += (workspace->tmp_nbgM + 3)->memsize;

    create_vec(nsM, workspace->tmp_nsM + 0, c_ptr);
    create_vec(nsM, workspace->res_workspace->tmp_nsM + 0, c_ptr);
    c_ptr += (workspace->tmp_nsM + 0)->memsize;

    for (ii = 0; ii <= N; ii++) {
        create_vec(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], workspace->tmp_m + ii, c_ptr);
        c_ptr += (2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii]) * sizeof(double);
    }

    d_create_core_qp_ipm(nvt, net, nct, cws, c_ptr);
    c_ptr += workspace->core_workspace->memsize;

    if (arg->lq_fact > 0) {
        workspace->lq_work0 = c_ptr;
        c_ptr += dgelqf_worksize(nuM + nxM, 2 * nuM + 3 * nxM + ngM);
    }


    // alias members of workspace and core_workspace
    //
    c_ptr = (char*) cws->dv;
    for (ii = 0; ii <= N; ii++) {
        create_vec(nu[ii] + nx[ii] + 2 * ns[ii], workspace->sol_step->ux + ii, c_ptr);
        c_ptr += (nu[ii] + nx[ii]) * sizeof(double);
        c_ptr += ns[ii] * sizeof(double);
        c_ptr += ns[ii] * sizeof(double);
    }
    //
    c_ptr = (char*) cws->dpi;
    for (ii = 0; ii < N; ii++) {
        create_vec(nx[ii + 1], workspace->sol_step->pi + ii, c_ptr);
        c_ptr += (nx[ii + 1]) * sizeof(double);
    }
    //
    c_ptr = (char*) cws->dlam;
    for (ii = 0; ii <= N; ii++) {
        create_vec(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], workspace->sol_step->lam + ii, c_ptr);
        c_ptr += nb[ii] * sizeof(double);
        c_ptr += ng[ii] * sizeof(double);
        c_ptr += nb[ii] * sizeof(double);
        c_ptr += ng[ii] * sizeof(double);
        c_ptr += ns[ii] * sizeof(double);
        c_ptr += ns[ii] * sizeof(double);
    }
    //
    c_ptr = (char*) cws->dt;
    for (ii = 0; ii <= N; ii++) {
        create_vec(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], workspace->sol_step->t + ii, c_ptr);
        c_ptr += nb[ii] * sizeof(double);
        c_ptr += ng[ii] * sizeof(double);
        c_ptr += nb[ii] * sizeof(double);
        c_ptr += ng[ii] * sizeof(double);
        c_ptr += ns[ii] * sizeof(double);
        c_ptr += ns[ii] * sizeof(double);
    }
    //
    c_ptr = (char*) cws->res_g;
    for (ii = 0; ii <= N; ii++) {
        create_vec(nu[ii] + nx[ii] + 2 * ns[ii], workspace->res->res_g + ii, c_ptr);
        c_ptr += (nu[ii] + nx[ii]) * sizeof(double);
        c_ptr += ns[ii] * sizeof(double);
        c_ptr += ns[ii] * sizeof(double);
    }
    //
    c_ptr = (char*) cws->res_b;
    for (ii = 0; ii < N; ii++) {
        create_vec(nx[ii + 1], workspace->res->res_b + ii, c_ptr);
        c_ptr += (nx[ii + 1]) * sizeof(double);
    }
    //
    c_ptr = (char*) cws->res_d;
    for (ii = 0; ii <= N; ii++) {
        create_vec(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], workspace->res->res_d + ii, c_ptr);
        c_ptr += nb[ii] * sizeof(double);
        c_ptr += ng[ii] * sizeof(double);
        c_ptr += nb[ii] * sizeof(double);
        c_ptr += ng[ii] * sizeof(double);
        c_ptr += ns[ii] * sizeof(double);
        c_ptr += ns[ii] * sizeof(double);
    }
    //
    c_ptr = (char*) cws->res_m;
    for (ii = 0; ii <= N; ii++) {
        create_vec(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], workspace->res->res_m + ii, c_ptr);
        c_ptr += nb[ii] * sizeof(double);
        c_ptr += ng[ii] * sizeof(double);
        c_ptr += nb[ii] * sizeof(double);
        c_ptr += ng[ii] * sizeof(double);
        c_ptr += ns[ii] * sizeof(double);
        c_ptr += ns[ii] * sizeof(double);
    }
    //
    c_ptr = (char*) cws->Gamma;
    for (ii = 0; ii <= N; ii++) {
        create_vec(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], workspace->Gamma + ii, c_ptr);
        c_ptr += nb[ii] * sizeof(double);
        c_ptr += ng[ii] * sizeof(double);
        c_ptr += nb[ii] * sizeof(double);
        c_ptr += ng[ii] * sizeof(double);
        c_ptr += ns[ii] * sizeof(double);
        c_ptr += ns[ii] * sizeof(double);
    }
    //
    c_ptr = (char*) cws->gamma;
    for (ii = 0; ii <= N; ii++) {
        create_vec(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], workspace->gamma + ii, c_ptr);
        c_ptr += nb[ii] * sizeof(double);
        c_ptr += ng[ii] * sizeof(double);
        c_ptr += nb[ii] * sizeof(double);
        c_ptr += ng[ii] * sizeof(double);
        c_ptr += ns[ii] * sizeof(double);
        c_ptr += ns[ii] * sizeof(double);
    }


    workspace->res->dim = dim;

    workspace->stat_max = arg->stat_max;
    workspace->stat_m = stat_m;

    for (ii = 0; ii <= N; ii++)
        workspace->use_hess_fact[ii] = 0;

    workspace->use_Pb = 0;

    workspace->valid_ric_vec = 0;

    // cache stuff
    workspace->dim = dim;
    workspace->square_root_alg = arg->square_root_alg;
    workspace->lq_fact = arg->lq_fact;

    workspace->memsize = memsize;  // d_ocp_qp_ipm_ws_memsize(dim, arg);


#if defined(RUNTIME_CHECKS)
    if (c_ptr > ((char*) mem) + workspace->memsize) {
        printf("\nCreate_ocp_qp_ipm: outside memory bounds!\n\n");
        exit(1);
    }
#endif
}


void d_ocp_qp_ipm_get(char* field, struct d_ocp_qp_ipm_ws* ws, void* value) {
    if (hpipm_strcmp(field, "status")) {
        d_ocp_qp_ipm_get_status(ws, value);
    } else if (hpipm_strcmp(field, "iter")) {
        d_ocp_qp_ipm_get_iter(ws, value);
    } else if (hpipm_strcmp(field, "max_res_stat")) {
        d_ocp_qp_ipm_get_max_res_stat(ws, value);
    } else if (hpipm_strcmp(field, "max_res_eq")) {
        d_ocp_qp_ipm_get_max_res_eq(ws, value);
    } else if (hpipm_strcmp(field, "max_res_ineq")) {
        d_ocp_qp_ipm_get_max_res_ineq(ws, value);
    } else if (hpipm_strcmp(field, "max_res_comp")) {
        d_ocp_qp_ipm_get_max_res_comp(ws, value);
    } else if (hpipm_strcmp(field, "obj")) {
        d_ocp_qp_ipm_get_obj(ws, value);
    } else if (hpipm_strcmp(field, "stat")) {
        d_ocp_qp_ipm_get_stat(ws, value);
    } else if (hpipm_strcmp(field, "stat_m")) {
        d_ocp_qp_ipm_get_stat_m(ws, value);
    } else {
        printf("error: d_ocp_qp_ipm_get: wrong field %s\n", field);
        exit(1);
    }
}


void d_ocp_qp_ipm_get_status(struct d_ocp_qp_ipm_ws* ws, int* status) {
    *status = ws->status;
}


void d_ocp_qp_ipm_get_iter(struct d_ocp_qp_ipm_ws* ws, int* iter) {
    *iter = ws->iter;
}


void d_ocp_qp_ipm_get_max_res_stat(struct d_ocp_qp_ipm_ws* ws, double* res_stat) {
    *res_stat = ws->res->res_max[0];
}


void d_ocp_qp_ipm_get_max_res_eq(struct d_ocp_qp_ipm_ws* ws, double* res_eq) {
    *res_eq = ws->res->res_max[1];
}


void d_ocp_qp_ipm_get_max_res_ineq(struct d_ocp_qp_ipm_ws* ws, double* res_ineq) {
    *res_ineq = ws->res->res_max[2];
}


void d_ocp_qp_ipm_get_max_res_comp(struct d_ocp_qp_ipm_ws* ws, double* res_comp) {
    *res_comp = ws->res->res_max[3];
}


void d_ocp_qp_ipm_get_obj(struct d_ocp_qp_ipm_ws* ws, double* obj) {
    *obj = ws->res->obj;
}


void d_ocp_qp_ipm_get_stat(struct d_ocp_qp_ipm_ws* ws, double** stat) {
    *stat = ws->stat;
}


void d_ocp_qp_ipm_get_stat_m(struct d_ocp_qp_ipm_ws* ws, int* stat_m) {
    *stat_m = ws->stat_m;
}


void d_ocp_qp_ipm_get_ric_Lr(struct d_ocp_qp* qp, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws, int stage, double* Lr) {
    int* nu = ws->dim->nu;

    int nu0 = nu[stage];

    unpack_mat(nu0, nu0, ws->L + stage, 0, 0, Lr, nu0);
}


void d_ocp_qp_ipm_get_ric_Ls(struct d_ocp_qp* qp, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws, int stage, double* Ls) {
    int* nu = ws->dim->nu;
    int* nx = ws->dim->nx;

    int nu0 = nu[stage];
    int nx0 = nx[stage];

    unpack_mat(nx0, nu0, ws->L + stage, nu0, 0, Ls, nx0);
}


void d_ocp_qp_ipm_get_ric_P(struct d_ocp_qp* qp, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws, int stage, double* P) {
    int* nu = ws->dim->nu;
    int* nx = ws->dim->nx;

    int nu0 = nu[stage];
    int nx0 = nx[stage];

    if (ws->square_root_alg == 1 | stage == 0) {
        dgese(nx0, nx0, 0.0, ws->tmp_nxM_nxM + 0, 0, 0);
        dtrcp_l(nx0, ws->L + stage, nu0, nu0, ws->tmp_nxM_nxM + 0, 0, 0);
        dsyrk_ln(nx0, nx0, 1.0, ws->tmp_nxM_nxM + 0, 0, 0, ws->tmp_nxM_nxM + 0, 0, 0, 0.0, ws->tmp_nxM_nxM + 1, 0, 0, ws->tmp_nxM_nxM + 1, 0, 0);  // TODO lauum
        dtrtr_l(nx0, ws->tmp_nxM_nxM + 1, 0, 0, ws->tmp_nxM_nxM + 1, 0, 0);
        unpack_mat(nx0, nx0, ws->tmp_nxM_nxM + 1, 0, 0, P, nx0);
    } else {
        unpack_mat(nx0, nx0, ws->P + stage, 0, 0, P, nx0);
    }
}


//
void d_ocp_qp_ipm_get_ric_lr(struct d_ocp_qp* qp, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws, int stage, double* lr) {
    int* nu = ws->dim->nu;
    int* nx = ws->dim->nx;

    int nu0 = nu[stage];
    int nx0 = nx[stage];

    if (ws->valid_ric_vec == 0) {
        int ii;

        struct d_core_qp_ipm_workspace* cws = ws->core_workspace;

        // arg to core workspace
        cws->lam_min = arg->lam_min;
        cws->t_min = arg->t_min;
        cws->t_min_inv = arg->t_min > 0 ? 1.0 / arg->t_min : 1e30;
        cws->t_lam_min = arg->t_lam_min;

        // alias qp vectors into qp_sol
        cws->v = ws->sol_itref->ux->pa;
        cws->pi = ws->sol_itref->pi->pa;
        cws->lam = ws->sol_itref->lam->pa;
        cws->t = ws->sol_itref->t->pa;

        // load sol from bkp
        for (ii = 0; ii < cws->nv; ii++)
            cws->v[ii] = cws->v_bkp[ii];
        for (ii = 0; ii < cws->ne; ii++)
            cws->pi[ii] = cws->pi_bkp[ii];
        for (ii = 0; ii < cws->nc; ii++)
            cws->lam[ii] = cws->lam_bkp[ii];
        for (ii = 0; ii < cws->nc; ii++)
            cws->t[ii] = cws->t_bkp[ii];

        ws->use_Pb = 0;
        d_ocp_qp_solve_kkt_step(qp, ws->sol_itref, arg, ws);

        ws->valid_ric_vec = 1;
    }

    //	unpack_mat(1, nu0, ws->L+stage, nu0+nx0, 0, lr, 1);
    unpack_vec(nu0, ws->l + stage, 0, lr, 1);
}


//
void d_ocp_qp_ipm_get_ric_p(struct d_ocp_qp* qp, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws, int stage, double* p) {
    int* nu = ws->dim->nu;
    int* nx = ws->dim->nx;

    int nu0 = nu[stage];
    int nx0 = nx[stage];

    if (ws->valid_ric_vec == 0) {
        int ii;

        struct d_core_qp_ipm_workspace* cws = ws->core_workspace;

        // arg to core workspace
        cws->lam_min = arg->lam_min;
        cws->t_min = arg->t_min;
        cws->t_min_inv = arg->t_min > 0 ? 1.0 / arg->t_min : 1e30;
        cws->t_lam_min = arg->t_lam_min;

        // alias qp vectors into qp_sol
        cws->v = ws->sol_itref->ux->pa;
        cws->pi = ws->sol_itref->pi->pa;
        cws->lam = ws->sol_itref->lam->pa;
        cws->t = ws->sol_itref->t->pa;

        // load sol from bkp
        for (ii = 0; ii < cws->nv; ii++)
            cws->v[ii] = cws->v_bkp[ii];
        for (ii = 0; ii < cws->ne; ii++)
            cws->pi[ii] = cws->pi_bkp[ii];
        for (ii = 0; ii < cws->nc; ii++)
            cws->lam[ii] = cws->lam_bkp[ii];
        for (ii = 0; ii < cws->nc; ii++)
            cws->t[ii] = cws->t_bkp[ii];

        ws->use_Pb = 0;
        d_ocp_qp_solve_kkt_step(qp, ws->sol_itref, arg, ws);

        ws->valid_ric_vec = 1;
    }

    //	if(ws->square_root_alg)
    if (ws->valid_ric_p == 0 | stage == 0) {
        //		drowex(nx0, 1.0, ws->L+stage, nu0+nx0, nu0, ws->tmp_nuxM, 0);
        //		dtrmv_lnn(nx0, ws->L+stage, nu0, nu0, ws->tmp_nuxM, 0, ws->tmp_nuxM, 0);
        dtrmv_lnn(nx0, ws->L + stage, nu0, nu0, ws->l + stage, nu0, ws->tmp_nuxM, 0);
        unpack_vec(nx0, ws->tmp_nuxM, 0, p, 1);
    } else {
        //		unpack_mat(1, nx0, ws->P+stage, nx0, 0, p, 1);
        unpack_vec(nx0, ws->l + stage, nu0, p, 1);
    }
}


// void d_ocp_qp_ipm_get_ric_K(struct d_ocp_qp* qp, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws, int stage, double* K) {
//     int* nu = ws->dim->nu;
//     int* nx = ws->dim->nx;

//     int nu0 = nu[stage];
//     int nx0 = nx[stage];

//     // XXX when implemented in HP, better copy-to-align first
//     dtrsm_rlnn(nx0, nu0, -1.0, ws->L + stage, 0, 0, ws->L + stage, nu0, 0, ws->Ls, 0, 0);

//     unpack_tran_mat(nx0, nu0, ws->Ls, 0, 0, K, nu0);
// }


// void d_ocp_qp_ipm_get_ric_k(struct d_ocp_qp* qp, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws, int stage, double* k) {
//     int* nu = ws->dim->nu;
//     int* nx = ws->dim->nx;

//     int nu0 = nu[stage];
//     int nx0 = nx[stage];

//     if (ws->valid_ric_vec == 0) {
//         int ii;

//         struct d_core_qp_ipm_workspace* cws = ws->core_workspace;

//         // arg to core workspace
//         cws->lam_min = arg->lam_min;
//         cws->t_min = arg->t_min;
//         cws->t_min_inv = arg->t_min > 0 ? 1.0 / arg->t_min : 1e30;
//         cws->t_lam_min = arg->t_lam_min;

//         // alias qp vectors into qp_sol
//         cws->v = ws->sol_itref->ux->pa;
//         cws->pi = ws->sol_itref->pi->pa;
//         cws->lam = ws->sol_itref->lam->pa;
//         cws->t = ws->sol_itref->t->pa;

//         // load sol from bkp
//         for (ii = 0; ii < cws->nv; ii++)
//             cws->v[ii] = cws->v_bkp[ii];
//         for (ii = 0; ii < cws->ne; ii++)
//             cws->pi[ii] = cws->pi_bkp[ii];
//         for (ii = 0; ii < cws->nc; ii++)
//             cws->lam[ii] = cws->lam_bkp[ii];
//         for (ii = 0; ii < cws->nc; ii++)
//             cws->t[ii] = cws->t_bkp[ii];

//         ws->use_Pb = 0;
//         d_ocp_qp_solve_kkt_step(qp, ws->sol_itref, arg, ws);

//         ws->valid_ric_vec = 1;
//     }

//     dtrsv_ltn(nu0, ws->L + stage, 0, 0, ws->l + stage, 0, ws->tmp_nuxM, 0);
//     dvecsc(nu0, -1.0, ws->tmp_nuxM, 0);
//     unpack_vec(nu0, ws->tmp_nuxM, 0, k, 1);
// }


// with warm_start==2 also init dual variables
void d_ocp_qp_init_var(struct d_ocp_qp* qp, struct d_ocp_qp_sol* qp_sol, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws) {

    //	struct d_core_qp_ipm_workspace *cws = ws->core_workspace;

    // loop index
    int ii, jj;

    //
    int N = qp->dim->N;
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;
    int* ns = qp->dim->ns;

    double mu0 = arg->mu0;

    //
    double *ux, *s, *pi, *d_lb, *d_ub, *d_lg, *d_ug, *d_ls, *lam_lb, *lam_ub, *lam_lg, *lam_ug, *lam_ls, *t_lb, *t_ub, *t_lg, *t_ug, *t_ls;
    int *idxb, *idxs_rev;
    int idx;

    double thr0 = 1e-1;


    // primal and dual variables
    if (arg->warm_start == 2) {

        thr0 = 1e-1;

        for (ii = 0; ii <= N; ii++) {
            lam_lb = qp_sol->lam[ii].pa + 0;
            t_lb = qp_sol->t[ii].pa + 0;

            for (jj = 0; jj < 2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii]; jj++) {
                if (lam_lb[jj] < thr0)
                    lam_lb[jj] = thr0;
                if (t_lb[jj] < thr0)
                    t_lb[jj] = thr0;
            }
        }
    }


    // ux
    if (arg->warm_start == 0) {
        // cold start
        for (ii = 0; ii <= N; ii++) {
            ux = qp_sol->ux[ii].pa;
            for (jj = 0; jj < nu[ii] + nx[ii] + 2 * ns[ii]; jj++) {
                ux[jj] = 0.0;
            }
        }
    } else {
        // warm start (keep u and x in solution)
        for (ii = 0; ii <= N; ii++) {
            ux = qp_sol->ux[ii].pa;
            for (jj = nu[ii] + nx[ii]; jj < nu[ii] + nx[ii] + 2 * ns[ii]; jj++) {
                ux[jj] = 0.0;
            }
        }
    }

    // pi
    for (ii = 0; ii < N; ii++) {
        pi = qp_sol->pi[ii].pa;
        for (jj = 0; jj < nx[ii + 1]; jj++) {
            pi[jj] = 0.0;
        }
    }

    if (arg->var_init_scheme == 0)  // safest scheme, no tailored init for soft constr
    {

        // box constraints
        for (ii = 0; ii <= N; ii++) {
            ux = qp_sol->ux[ii].pa;
            d_lb = qp->d[ii].pa + 0;
            d_ub = qp->d[ii].pa + nb[ii] + ng[ii];
            lam_lb = qp_sol->lam[ii].pa + 0;
            lam_ub = qp_sol->lam[ii].pa + nb[ii] + ng[ii];
            t_lb = qp_sol->t[ii].pa + 0;
            t_ub = qp_sol->t[ii].pa + nb[ii] + ng[ii];
            idxb = qp->idxb[ii];
            for (jj = 0; jj < nb[ii]; jj++) {
#if 1
                t_lb[jj] = -d_lb[jj] + ux[idxb[jj]];
                t_ub[jj] = -d_ub[jj] - ux[idxb[jj]];
                //			printf("\n%d %f %f\n", jj, t_lb[jj], t_ub[jj]);
                if (t_lb[jj] < thr0) {
                    if (t_ub[jj] < thr0) {
                        //					ux[idxb[jj]] = 0.5*(d_lb[jj] + d_ub[jj]);
                        ux[idxb[jj]] = 0.5 * (d_lb[jj] - d_ub[jj]);
                        t_lb[jj] = thr0;
                        t_ub[jj] = thr0;
                    } else {
                        t_lb[jj] = thr0;
                        ux[idxb[jj]] = d_lb[jj] + thr0;
                    }
                } else if (t_ub[jj] < thr0) {
                    t_ub[jj] = thr0;
                    ux[idxb[jj]] = -d_ub[jj] - thr0;
                }
#else
                t_lb[jj] = 1.0;
                t_ub[jj] = 1.0;
#endif
                lam_lb[jj] = mu0 / t_lb[jj];
                lam_ub[jj] = mu0 / t_ub[jj];
            }
            //		print_tran_vec(nb[ii], qp->d+ii, 0);
            //		print_tran_vec(nb[ii], qp->d+ii, nb[ii]+ng[ii]);
            //		print_tran_vec(nu[ii]+nx[ii], qp_sol->ux+ii, 0);
            //		print_tran_vec(nb[ii], qp_sol->t+ii, 0);
            //		print_tran_vec(nb[ii], qp_sol->t+ii, nb[ii]+ng[ii]);
            //		exit(1);
        }

        // general constraints
        for (ii = 0; ii <= N; ii++) {
            t_lg = qp_sol->t[ii].pa + nb[ii];
            t_ug = qp_sol->t[ii].pa + 2 * nb[ii] + ng[ii];
            lam_lg = qp_sol->lam[ii].pa + nb[ii];
            lam_ug = qp_sol->lam[ii].pa + 2 * nb[ii] + ng[ii];
            d_lg = qp->d[ii].pa + nb[ii];
            d_ug = qp->d[ii].pa + 2 * nb[ii] + ng[ii];
            ux = qp_sol->ux[ii].pa;
            dgemv_t(nu[ii] + nx[ii], ng[ii], 1.0, qp->DCt + ii, 0, 0, qp_sol->ux + ii, 0, 0.0, qp_sol->t + ii, nb[ii], qp_sol->t + ii, nb[ii]);
            for (jj = 0; jj < ng[ii]; jj++) {
#if 1
                t_ug[jj] = -t_lg[jj];
                t_lg[jj] -= d_lg[jj];
                t_ug[jj] -= d_ug[jj];
                //			t_lg[jj] = fmax(thr0, t_lg[jj]);
                //			t_ug[jj] = fmax(thr0, t_ug[jj]);
                t_lg[jj] = thr0 > t_lg[jj] ? thr0 : t_lg[jj];
                t_ug[jj] = thr0 > t_ug[jj] ? thr0 : t_ug[jj];
#else
                t_lg[jj] = 1.0;
                t_ug[jj] = 1.0;
#endif
                lam_lg[jj] = mu0 / t_lg[jj];
                lam_ug[jj] = mu0 / t_ug[jj];
            }
        }

        // soft constraints
        for (ii = 0; ii <= N; ii++) {
            lam_lb = qp_sol->lam[ii].pa + 2 * nb[ii] + 2 * ng[ii];
            lam_ub = qp_sol->lam[ii].pa + 2 * nb[ii] + 2 * ng[ii] + ns[ii];
            t_lb = qp_sol->t[ii].pa + 2 * nb[ii] + 2 * ng[ii];
            t_ub = qp_sol->t[ii].pa + 2 * nb[ii] + 2 * ng[ii] + ns[ii];
            for (jj = 0; jj < ns[ii]; jj++) {
                t_lb[jj] = 1.0;  // thr0;
                t_ub[jj] = 1.0;  // thr0;
                //			t_lb[jj] = sqrt(mu0); // thr0;
                //			t_ub[jj] = sqrt(mu0); // thr0;
                lam_lb[jj] = mu0 / t_lb[jj];
                lam_ub[jj] = mu0 / t_ub[jj];
            }
        }

    } else  // alternative scheme for soft constr
    {

        for (ii = 0; ii <= N; ii++) {

            //		printf("\nii = %d\n", ii);

            ux = qp_sol->ux[ii].pa;
            s = qp_sol->ux[ii].pa + nu[ii] + nx[ii];
            d_lb = qp->d[ii].pa + 0;
            d_ub = qp->d[ii].pa + nb[ii] + ng[ii];
            d_lg = qp->d[ii].pa + nb[ii];
            d_ug = qp->d[ii].pa + 2 * nb[ii] + ng[ii];
            d_ls = qp->d[ii].pa + 2 * nb[ii] + 2 * ng[ii];
            lam_lb = qp_sol->lam[ii].pa + 0;
            lam_ub = qp_sol->lam[ii].pa + nb[ii] + ng[ii];
            lam_lg = qp_sol->lam[ii].pa + nb[ii];
            lam_ug = qp_sol->lam[ii].pa + 2 * nb[ii] + ng[ii];
            lam_ls = qp_sol->lam[ii].pa + 2 * nb[ii] + 2 * ng[ii];
            t_lb = qp_sol->t[ii].pa + 0;
            t_ub = qp_sol->t[ii].pa + nb[ii] + ng[ii];
            t_lg = qp_sol->t[ii].pa + nb[ii];
            t_ug = qp_sol->t[ii].pa + 2 * nb[ii] + ng[ii];
            t_ls = qp_sol->t[ii].pa + 2 * nb[ii] + 2 * ng[ii];
            idxb = qp->idxb[ii];
            idxs_rev = qp->idxs_rev[ii];

            // lower bound on slacks
            daxpy(2 * ns[ii], -1.0, qp->d + ii, 2 * nb[ii] + 2 * ng[ii], qp_sol->ux + ii, nu[ii] + nx[ii], qp_sol->t + ii, 2 * nb[ii] + 2 * ng[ii]);
            for (jj = 0; jj < 2 * ns[ii]; jj++) {
#if 1
                if (t_ls[jj] < thr0) {
                    t_ls[jj] = thr0;  // 1.0;
                    s[jj] = d_ls[jj] + t_ls[jj];
                }
#else
                t_ls[jj] = 1.0;
                //			t_ls[jj] = sqrt(mu0);
#endif
            }
            //		print_tran_vec(2*ns[ii], qp_sol->ux+ii, nu[ii]+nx[ii]);
            //		print_tran_vec(2*ns[ii], qp_sol->t+ii, 2*nb[ii]+2*ng[ii]);

            // upper and lower bounds on inputs and states
            dvecex_sp(nb[ii], 1.0, qp->idxb[ii], qp_sol->ux + ii, 0, qp_sol->t + ii, 0);
            dveccpsc(nb[ii], -1.0, qp_sol->t + ii, 0, qp_sol->t + ii, nb[ii] + ng[ii]);
            for (jj = 0; jj < nb[ii]; jj++) {
                idx = idxs_rev[jj];
                if (idx != -1) {
                    // softed bound
                    t_lb[jj] += s[idx];
                    t_ub[jj] += s[ns[ii] + idx];
                }
            }
            daxpy(nb[ii], -1.0, qp->d + ii, 0, qp_sol->t + ii, 0, qp_sol->t + ii, 0);
            daxpy(nb[ii], -1.0, qp->d + ii, nb[ii] + ng[ii], qp_sol->t + ii, nb[ii] + ng[ii], qp_sol->t + ii, nb[ii] + ng[ii]);
            //		print_tran_vec(nb[ii], qp_sol->t+ii, 0);
            //		print_tran_vec(nb[ii], qp_sol->t+ii, nb[ii]+ng[ii]);
            for (jj = 0; jj < nb[ii]; jj++) {
#if 1
                if (t_lb[jj] < thr0) {
                    if (t_ub[jj] < thr0) {
                        //					ux[idxb[jj]] = 0.5*(d_lb[jj] + d_ub[jj]);
                        ux[idxb[jj]] = 0.5 * (d_lb[jj] - d_ub[jj]);
                        t_lb[jj] = thr0;
                        t_ub[jj] = thr0;
                    } else {
                        t_lb[jj] = thr0;
                        ux[idxb[jj]] = d_lb[jj] + thr0;
                    }
                } else if (t_ub[jj] < thr0) {
                    t_ub[jj] = thr0;
                    ux[idxb[jj]] = -d_ub[jj] - thr0;
                }
#else
                t_lb[jj] = 1.0;
                t_ub[jj] = 1.0;
#endif
            }
            //		print_tran_vec(nu[ii]+nx[ii], qp_sol->ux+ii, 0);
            //		print_tran_vec(nb[ii], qp_sol->t+ii, 0);
            //		print_tran_vec(nb[ii], qp_sol->t+ii, nb[ii]+ng[ii]);

            // upper and lower general constaints
            dgemv_t(nu[ii] + nx[ii], ng[ii], 1.0, qp->DCt + ii, 0, 0, qp_sol->ux + ii, 0, 0.0, qp_sol->t + ii, nb[ii], qp_sol->t + ii, nb[ii]);
            dveccpsc(ng[ii], -1.0, qp_sol->t + ii, nb[ii], qp_sol->t + ii, 2 * nb[ii] + ng[ii]);
            //		print_tran_vec(ng[ii], qp_sol->t+ii, nb[ii]);
            //		print_tran_vec(ng[ii], qp_sol->t+ii, 2*nb[ii]+ng[ii]);
            for (jj = 0; jj < ng[ii]; jj++) {
                idx = idxs_rev[nb[ii] + jj];
                if (idx != -1) {
                    // softed general constraint
                    t_lb[nb[ii] + jj] += s[idx];
                    t_ub[nb[ii] + jj] += s[ns[ii] + idx];
                }
            }
            //		print_tran_vec(ng[ii], qp_sol->t+ii, nb[ii]);
            //		print_tran_vec(ng[ii], qp_sol->t+ii, 2*nb[ii]+ng[ii]);
            daxpy(ng[ii], -1.0, qp->d + ii, nb[ii], qp_sol->t + ii, nb[ii], qp_sol->t + ii, nb[ii]);
            daxpy(ng[ii], -1.0, qp->d + ii, 2 * nb[ii] + ng[ii], qp_sol->t + ii, 2 * nb[ii] + ng[ii], qp_sol->t + ii, 2 * nb[ii] + ng[ii]);
            for (jj = 0; jj < ng[ii]; jj++) {
#if 1
                t_lg[jj] = thr0 > t_lg[jj] ? thr0 : t_lg[jj];
                t_ug[jj] = thr0 > t_ug[jj] ? thr0 : t_ug[jj];
#else
                t_lg[jj] = 1.0;
                t_ug[jj] = 1.0;
#endif
            }
            //		print_tran_vec(ng[ii], qp_sol->t+ii, nb[ii]);
            //		print_tran_vec(ng[ii], qp_sol->t+ii, 2*nb[ii]+ng[ii]);

            // multipliers
            for (jj = 0; jj < 2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii]; jj++)
                lam_lb[jj] = mu0 / t_lb[jj];
        }
    }
}


void d_ocp_qp_ipm_abs_step(int kk, struct d_ocp_qp* qp, struct d_ocp_qp_sol* qp_sol, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws) {

    struct d_core_qp_ipm_workspace* cws = ws->core_workspace;

    int N = qp->dim->N;
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;
    int* ns = qp->dim->ns;

    int ii;

    double tmp;
    double mu_aff0;  //, mu;

    double* stat = ws->stat;
    int stat_m = ws->stat_m;

    dvecsc(cws->nc, -1.0, ws->tmp_m, 0);

    d_backup_res_m(cws);

    // tau_min as barrier parameter for affine step
    d_compute_tau_min_qp(cws);

    // fact solve
    // d_ocp_qp_print(ws->qp_step->dim, ws->qp_step);
    // exit(1);
    d_ocp_qp_fact_solve_kkt_step(ws->qp_step, ws->sol_step, arg, ws);
    // print_tran_vec(cws->nv, ws->sol_step->ux, 0);
    // d_ocp_qp_sol_print(ws->qp_step->dim, ws->sol_step);
    // exit(1);

    // compute step
    daxpy(cws->nv, -1.0, qp_sol->ux, 0, ws->sol_step->ux, 0, ws->sol_step->ux, 0);
    daxpy(cws->ne, -1.0, qp_sol->pi, 0, ws->sol_step->pi, 0, ws->sol_step->pi, 0);
    daxpy(cws->nc, -1.0, qp_sol->lam, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
    daxpy(cws->nc, -1.0, qp_sol->t, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);
    if (ws->mask_constr) {
        // mask out disregarded constraints
        for (ii = 0; ii <= N; ii++)
            dvecmul(2 * ns[ii], qp->d_mask + ii, 2 * nb[ii] + 2 * ng[ii], ws->sol_step->ux + ii, nu[ii] + nx[ii], ws->sol_step->ux + ii, nu[ii] + nx[ii]);
        dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);
        dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
    }

    // alpha
    d_compute_alpha_qp(cws);
    if (kk + 1 < ws->stat_max)
        stat[stat_m * (kk + 1) + 0] = cws->alpha;

    // Mehrotra's predictor-corrector
    if (arg->pred_corr == 1) {
        // mu_aff
        d_compute_mu_aff_qp(cws);
        if (kk + 1 < ws->stat_max)
            stat[stat_m * (kk + 1) + 1] = cws->mu_aff;

        tmp = cws->mu_aff / cws->mu;
        cws->sigma = tmp * tmp * tmp;
        if (kk + 1 < ws->stat_max)
            stat[stat_m * (kk + 1) + 2] = cws->sigma;

        d_compute_centering_correction_qp(cws);
        if (ws->mask_constr) {
            // mask out disregarded constraints
            dvecmul(cws->nc, qp->d_mask, 0, ws->res->res_m, 0, ws->res->res_m, 0);
        }

        // fact and solve kkt
        ws->use_Pb = 1;
        d_ocp_qp_solve_kkt_step(ws->qp_step, ws->sol_step, arg, ws);

        // compute step
        daxpy(cws->nv, -1.0, qp_sol->ux, 0, ws->sol_step->ux, 0, ws->sol_step->ux, 0);
        daxpy(cws->ne, -1.0, qp_sol->pi, 0, ws->sol_step->pi, 0, ws->sol_step->pi, 0);
        daxpy(cws->nc, -1.0, qp_sol->lam, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
        daxpy(cws->nc, -1.0, qp_sol->t, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);
        if (ws->mask_constr) {
            // mask out disregarded constraints
            for (ii = 0; ii <= N; ii++)
                dvecmul(2 * ns[ii], qp->d_mask + ii, 2 * nb[ii] + 2 * ng[ii], ws->sol_step->ux + ii, nu[ii] + nx[ii], ws->sol_step->ux + ii, nu[ii] + nx[ii]);
            dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);
            dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
        }

        // alpha
        d_compute_alpha_qp(cws);
        if (kk + 1 < ws->stat_max) {
            stat[stat_m * (kk + 1) + 3] = cws->alpha_prim;
            stat[stat_m * (kk + 1) + 4] = cws->alpha_dual;
        }

        // conditional Mehrotra's predictor-corrector
        if (arg->cond_pred_corr == 1) {

            // save mu_aff (from prediction sol_step)
            mu_aff0 = cws->mu_aff;

            // compute mu for predictor-corrector-centering
            d_compute_mu_aff_qp(cws);

            //			if(cws->mu_aff > 2.0*cws->mu)
            if (cws->mu_aff > 2.0 * mu_aff0) {

                // centering direction
                d_compute_centering_qp(cws);
                if (ws->mask_constr) {
                    // mask out disregarded constraints
                    dvecmul(cws->nc, qp->d_mask, 0, ws->res->res_m, 0, ws->res->res_m, 0);
                }

                // solve kkt
                d_ocp_qp_solve_kkt_step(ws->qp_step, ws->sol_step, arg, ws);
                // compute step
                daxpy(cws->nv, -1.0, qp_sol->ux, 0, ws->sol_step->ux, 0, ws->sol_step->ux, 0);
                daxpy(cws->ne, -1.0, qp_sol->pi, 0, ws->sol_step->pi, 0, ws->sol_step->pi, 0);
                daxpy(cws->nc, -1.0, qp_sol->lam, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
                daxpy(cws->nc, -1.0, qp_sol->t, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);
                if (ws->mask_constr) {
                    // mask out disregarded constraints
                    for (ii = 0; ii <= N; ii++)
                        dvecmul(2 * ns[ii], qp->d_mask + ii, 2 * nb[ii] + 2 * ng[ii], ws->sol_step->ux + ii, nu[ii] + nx[ii], ws->sol_step->ux + ii, nu[ii] + nx[ii]);
                    dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);
                    dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
                }

                // alpha
                d_compute_alpha_qp(cws);
                if (kk + 1 < ws->stat_max) {
                    stat[stat_m * (kk + 1) + 3] = cws->alpha_prim;
                    stat[stat_m * (kk + 1) + 4] = cws->alpha_dual;
                }
            }
        }
    }

    //
    d_update_var_qp(cws);
    if (ws->mask_constr) {
        // mask out disregarded constraints
        dvecmul(cws->nc, qp->d_mask, 0, qp_sol->lam, 0, qp_sol->lam, 0);
    }
}


void d_ocp_qp_ipm_delta_step(int kk, struct d_ocp_qp* qp, struct d_ocp_qp_sol* qp_sol, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws) {

    // d_ocp_qp_print(qp->dim, qp);
    // d_ocp_qp_sol_print(qp->dim, qp_sol);
    // exit(1);
    //  dim
    int N = qp->dim->N;
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;
    int* ns = qp->dim->ns;

    struct d_core_qp_ipm_workspace* cws = ws->core_workspace;

    int itref0 = 0, itref1 = 0, iter_ref_step;
    int ii;
    double tmp;
    double mu_aff0, mu;

    double* stat = ws->stat;
    int stat_m = ws->stat_m;

    double itref_qp_norm[4] = {0, 0, 0, 0};
    double itref_qp_norm0[4] = {0, 0, 0, 0};

    int force_lq = 0;

    // step body

    d_backup_res_m(cws);

    // tau_min as barrier parameter for affine step
    d_compute_tau_min_qp(cws);

    // fact and solve kkt
    if (ws->lq_fact == 0) {

        // syrk+cholesky
        d_ocp_qp_fact_solve_kkt_step(ws->qp_step, ws->sol_step, arg, ws);
        if (ws->mask_constr) {
            // mask out disregarded constraints
            for (ii = 0; ii <= N; ii++)
                dvecmul(2 * ns[ii], qp->d_mask + ii, 2 * nb[ii] + 2 * ng[ii], ws->sol_step->ux + ii, nu[ii] + nx[ii], ws->sol_step->ux + ii, nu[ii] + nx[ii]);
            dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);
            dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
        }
        if (kk + 1 < ws->stat_max)
            stat[stat_m * (kk + 1) + 11] = 0;
    } else if (ws->lq_fact == 1 & force_lq == 0) {

        // syrk+chol, switch to lq when needed
        d_ocp_qp_fact_solve_kkt_step(ws->qp_step, ws->sol_step, arg, ws);
        if (ws->mask_constr) {
            // mask out disregarded constraints
            for (ii = 0; ii <= N; ii++)
                dvecmul(2 * ns[ii], qp->d_mask + ii, 2 * nb[ii] + 2 * ng[ii], ws->sol_step->ux + ii, nu[ii] + nx[ii], ws->sol_step->ux + ii, nu[ii] + nx[ii]);
            dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);
            dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
        }

        // compute res of linear system
        d_ocp_qp_res_compute_lin(ws->qp_step, qp_sol, ws->sol_step, ws->res_itref, ws->res_workspace);
        if (ws->mask_constr) {
            // mask out disregarded constraints
            for (ii = 0; ii <= N; ii++)
                dvecmul(2 * ns[ii], qp->d_mask + ii, 2 * nb[ii] + 2 * ng[ii], ws->res_itref->res_g + ii, nu[ii] + nx[ii], ws->res_itref->res_g + ii, nu[ii] + nx[ii]);
            dvecmul(cws->nc, qp->d_mask, 0, ws->res_itref->res_d, 0, ws->res_itref->res_d, 0);
            dvecmul(cws->nc, qp->d_mask, 0, ws->res_itref->res_m, 0, ws->res_itref->res_m, 0);
        }
        d_ocp_qp_res_compute_inf_norm(ws->res_itref);
        itref_qp_norm[0] = ws->res_itref->res_max[0];
        itref_qp_norm[1] = ws->res_itref->res_max[1];
        itref_qp_norm[2] = ws->res_itref->res_max[2];
        itref_qp_norm[3] = ws->res_itref->res_max[3];

        if (kk + 1 < ws->stat_max)
            stat[stat_m * (kk + 1) + 11] = 0;

        // inaccurate factorization: switch to lq
        if (

                ((itref_qp_norm[0] == 0.0) & isnan(VECEL(ws->res_itref->res_g + 0, 0))) |
                // #else
                //                 ((itref_qp_norm[0] == 0.0) & (VECEL(ws->res_itref->res_g + 0, 0) != VECEL(ws->res_itref->res_g + 0, 0))) |

                (itref_qp_norm[0] > 1e-5) |
                (itref_qp_norm[1] > 1e-5) |
                (itref_qp_norm[2] > 1e-5) |
                (itref_qp_norm[3] > 1e-5)) {

            // refactorize using lq
            d_ocp_qp_fact_lq_solve_kkt_step(ws->qp_step, ws->sol_step, arg, ws);
            if (ws->mask_constr) {
                // mask out disregarded constraints
                for (ii = 0; ii <= N; ii++)
                    dvecmul(2 * ns[ii], qp->d_mask + ii, 2 * nb[ii] + 2 * ng[ii], ws->sol_step->ux + ii, nu[ii] + nx[ii], ws->sol_step->ux + ii, nu[ii] + nx[ii]);
                dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);
                dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
            }

            // switch to lq
            force_lq = 1;

            if (kk + 1 < ws->stat_max)
                stat[stat_m * (kk + 1) + 11] = 1;
        }

    } else  // ws->lq_fact==2
    {

        d_ocp_qp_fact_lq_solve_kkt_step(ws->qp_step, ws->sol_step, arg, ws);
        if (ws->mask_constr) {
            // mask out disregarded constraints
            for (ii = 0; ii <= N; ii++)
                dvecmul(2 * ns[ii], qp->d_mask + ii, 2 * nb[ii] + 2 * ng[ii], ws->sol_step->ux + ii, nu[ii] + nx[ii], ws->sol_step->ux + ii, nu[ii] + nx[ii]);
            dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);
            dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
        }
        if (kk + 1 < ws->stat_max)
            stat[stat_m * (kk + 1) + 11] = 1;
    }

    // iterative refinement on prediction step
    if (arg->itref_pred_max == 0) {
        if (kk + 1 < ws->stat_max) {
            stat[stat_m * (kk + 1) + 14] = 0.0;
            stat[stat_m * (kk + 1) + 15] = 0.0;
            stat[stat_m * (kk + 1) + 16] = 0.0;
            stat[stat_m * (kk + 1) + 17] = 0.0;
        }
    } else {
        for (itref0 = 0; itref0 < arg->itref_pred_max; itref0++) {

            d_ocp_qp_res_compute_lin(ws->qp_step, qp_sol, ws->sol_step, ws->res_itref, ws->res_workspace);
            if (ws->mask_constr) {
                // mask out disregarded constraints
                for (ii = 0; ii <= N; ii++)
                    dvecmul(2 * ns[ii], qp->d_mask + ii, 2 * nb[ii] + 2 * ng[ii], ws->res_itref->res_g + ii, nu[ii] + nx[ii], ws->res_itref->res_g + ii, nu[ii] + nx[ii]);
                dvecmul(cws->nc, qp->d_mask, 0, ws->res_itref->res_d, 0, ws->res_itref->res_d, 0);
                dvecmul(cws->nc, qp->d_mask, 0, ws->res_itref->res_m, 0, ws->res_itref->res_m, 0);
            }
            d_ocp_qp_res_compute_inf_norm(ws->res_itref);
            itref_qp_norm[0] = ws->res_itref->res_max[0];
            itref_qp_norm[1] = ws->res_itref->res_max[1];
            itref_qp_norm[2] = ws->res_itref->res_max[2];
            itref_qp_norm[3] = ws->res_itref->res_max[3];
            if (kk + 1 < ws->stat_max) {
                stat[stat_m * (kk + 1) + 14] = itref_qp_norm[0];
                stat[stat_m * (kk + 1) + 15] = itref_qp_norm[1];
                stat[stat_m * (kk + 1) + 16] = itref_qp_norm[2];
                stat[stat_m * (kk + 1) + 17] = itref_qp_norm[3];
            }

            if (itref0 == 0) {
                itref_qp_norm0[0] = itref_qp_norm[0];
                itref_qp_norm0[1] = itref_qp_norm[1];
                itref_qp_norm0[2] = itref_qp_norm[2];
                itref_qp_norm0[3] = itref_qp_norm[3];
            }

            if (
                    (itref_qp_norm[0] < 1e0 * arg->res_g_max | itref_qp_norm[0] < 1e-3 * ws->res->res_max[0]) &
                    (itref_qp_norm[1] < 1e0 * arg->res_b_max | itref_qp_norm[1] < 1e-3 * ws->res->res_max[1]) &
                    (itref_qp_norm[2] < 1e0 * arg->res_d_max | itref_qp_norm[2] < 1e-3 * ws->res->res_max[2]) &
                    (itref_qp_norm[3] < 1e0 * arg->res_m_max | itref_qp_norm[3] < 1e-3 * ws->res->res_max[3])) {
                break;
            }

            ws->use_Pb = 0;
            d_ocp_qp_solve_kkt_step(ws->qp_itref, ws->sol_itref, arg, ws);
            if (ws->mask_constr) {
                // mask out disregarded constraints
                for (ii = 0; ii <= N; ii++)
                    dvecmul(2 * ns[ii], qp->d_mask + ii, 2 * nb[ii] + 2 * ng[ii], ws->sol_step->ux + ii, nu[ii] + nx[ii], ws->sol_step->ux + ii, nu[ii] + nx[ii]);
                dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);
                dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
            }

            for (ii = 0; ii <= N; ii++)
                daxpy(nu[ii] + nx[ii] + 2 * ns[ii], 1.0, ws->sol_itref->ux + ii, 0, ws->sol_step->ux + ii, 0, ws->sol_step->ux + ii, 0);
            for (ii = 0; ii < N; ii++)
                daxpy(nx[ii + 1], 1.0, ws->sol_itref->pi + ii, 0, ws->sol_step->pi + ii, 0, ws->sol_step->pi + ii, 0);
            for (ii = 0; ii <= N; ii++)
                daxpy(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], 1.0, ws->sol_itref->lam + ii, 0, ws->sol_step->lam + ii, 0, ws->sol_step->lam + ii, 0);
            for (ii = 0; ii <= N; ii++)
                daxpy(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], 1.0, ws->sol_itref->t + ii, 0, ws->sol_step->t + ii, 0, ws->sol_step->t + ii, 0);
        }
        if (itref0 == arg->itref_pred_max) {
            d_ocp_qp_res_compute_lin(ws->qp_step, qp_sol, ws->sol_step, ws->res_itref, ws->res_workspace);
            if (ws->mask_constr) {
                // mask out disregarded constraints
                for (ii = 0; ii <= N; ii++)
                    dvecmul(2 * ns[ii], qp->d_mask + ii, 2 * nb[ii] + 2 * ng[ii], ws->res_itref->res_g + ii, nu[ii] + nx[ii], ws->res_itref->res_g + ii, nu[ii] + nx[ii]);
                dvecmul(cws->nc, qp->d_mask, 0, ws->res_itref->res_d, 0, ws->res_itref->res_d, 0);
                dvecmul(cws->nc, qp->d_mask, 0, ws->res_itref->res_m, 0, ws->res_itref->res_m, 0);
            }
            d_ocp_qp_res_compute_inf_norm(ws->res_itref);
            itref_qp_norm[0] = ws->res_itref->res_max[0];
            itref_qp_norm[1] = ws->res_itref->res_max[1];
            itref_qp_norm[2] = ws->res_itref->res_max[2];
            itref_qp_norm[3] = ws->res_itref->res_max[3];
            if (kk + 1 < ws->stat_max) {
                stat[stat_m * (kk + 1) + 14] = itref_qp_norm[0];
                stat[stat_m * (kk + 1) + 15] = itref_qp_norm[1];
                stat[stat_m * (kk + 1) + 16] = itref_qp_norm[2];
                stat[stat_m * (kk + 1) + 17] = itref_qp_norm[3];
            }
        }
    }

    if (kk + 1 < ws->stat_max)
        stat[stat_m * (kk + 1) + 12] = itref0;

    // alpha
    d_compute_alpha_qp(cws);
    if (kk + 1 < ws->stat_max)
        stat[stat_m * (kk + 1) + 0] = cws->alpha;

    // Mehrotra's predictor-corrector
    if (arg->pred_corr == 1) {
        // mu_aff
        d_compute_mu_aff_qp(cws);
        if (kk + 1 < ws->stat_max)
            stat[stat_m * (kk + 1) + 1] = cws->mu_aff;

        tmp = cws->mu_aff / cws->mu;
        cws->sigma = tmp * tmp * tmp;
        if (kk + 1 < ws->stat_max)
            stat[stat_m * (kk + 1) + 2] = cws->sigma;

        d_compute_centering_correction_qp(cws);
        if (ws->mask_constr) {
            // mask out disregarded constraints
            dvecmul(cws->nc, qp->d_mask, 0, ws->res->res_m, 0, ws->res->res_m, 0);
        }

        // fact and solve kkt
        ws->use_Pb = 1;
        d_ocp_qp_solve_kkt_step(ws->qp_step, ws->sol_step, arg, ws);
        if (ws->mask_constr) {
            // mask out disregarded constraints
            for (ii = 0; ii <= N; ii++)
                dvecmul(2 * ns[ii], qp->d_mask + ii, 2 * nb[ii] + 2 * ng[ii], ws->sol_step->ux + ii, nu[ii] + nx[ii], ws->sol_step->ux + ii, nu[ii] + nx[ii]);
            dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);
            dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
        }

        // alpha
        d_compute_alpha_qp(cws);
        if (kk + 1 < ws->stat_max) {
            stat[stat_m * (kk + 1) + 3] = cws->alpha_prim;
            stat[stat_m * (kk + 1) + 4] = cws->alpha_dual;
        }

        // conditional Mehrotra's predictor-corrector
        if (arg->cond_pred_corr == 1) {

            // save mu_aff (from prediction step)
            mu_aff0 = cws->mu_aff;

            // compute mu for predictor-corrector-centering
            d_compute_mu_aff_qp(cws);

            //				if(cws->mu_aff > 2.0*cws->mu)
            if (cws->mu_aff > 2.0 * mu_aff0) {

                // centering direction
                d_compute_centering_qp(cws);
                if (ws->mask_constr) {
                    // mask out disregarded constraints
                    dvecmul(cws->nc, qp->d_mask, 0, ws->res->res_m, 0, ws->res->res_m, 0);
                }

                // solve kkt
                ws->use_Pb = 1;
                d_ocp_qp_solve_kkt_step(ws->qp_step, ws->sol_step, arg, ws);
                if (ws->mask_constr) {
                    // mask out disregarded constraints
                    for (ii = 0; ii <= N; ii++)
                        dvecmul(2 * ns[ii], qp->d_mask + ii, 2 * nb[ii] + 2 * ng[ii], ws->sol_step->ux + ii, nu[ii] + nx[ii], ws->sol_step->ux + ii, nu[ii] + nx[ii]);
                    dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);
                    dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
                }

                // alpha
                d_compute_alpha_qp(cws);
                if (kk + 1 < ws->stat_max) {
                    stat[stat_m * (kk + 1) + 3] = cws->alpha_prim;
                    stat[stat_m * (kk + 1) + 4] = cws->alpha_dual;
                }
            }
        }

        iter_ref_step = 0;
        if (arg->itref_corr_max > 0) {
            for (itref1 = 0; itref1 < arg->itref_corr_max; itref1++) {

                d_ocp_qp_res_compute_lin(ws->qp_step, qp_sol, ws->sol_step, ws->res_itref, ws->res_workspace);
                if (ws->mask_constr) {
                    // mask out disregarded constraints
                    for (ii = 0; ii <= N; ii++)
                        dvecmul(2 * ns[ii], qp->d_mask + ii, 2 * nb[ii] + 2 * ng[ii], ws->res_itref->res_g + ii, nu[ii] + nx[ii], ws->res_itref->res_g + ii, nu[ii] + nx[ii]);
                    dvecmul(cws->nc, qp->d_mask, 0, ws->res_itref->res_d, 0, ws->res_itref->res_d, 0);
                    dvecmul(cws->nc, qp->d_mask, 0, ws->res_itref->res_m, 0, ws->res_itref->res_m, 0);
                }
                d_ocp_qp_res_compute_inf_norm(ws->res_itref);
                itref_qp_norm[0] = ws->res_itref->res_max[0];
                itref_qp_norm[1] = ws->res_itref->res_max[1];
                itref_qp_norm[2] = ws->res_itref->res_max[2];
                itref_qp_norm[3] = ws->res_itref->res_max[3];
                if (kk + 1 < ws->stat_max) {
                    stat[stat_m * (kk + 1) + 14] = itref_qp_norm[0];
                    stat[stat_m * (kk + 1) + 15] = itref_qp_norm[1];
                    stat[stat_m * (kk + 1) + 16] = itref_qp_norm[2];
                    stat[stat_m * (kk + 1) + 17] = itref_qp_norm[3];
                }

                if (itref1 == 0) {
                    itref_qp_norm0[0] = itref_qp_norm[0];
                    itref_qp_norm0[1] = itref_qp_norm[1];
                    itref_qp_norm0[2] = itref_qp_norm[2];
                    itref_qp_norm0[3] = itref_qp_norm[3];
                }

                if (
                        (itref_qp_norm[0] < 1e0 * arg->res_g_max | itref_qp_norm[0] < 1e-3 * ws->res->res_max[0]) &
                        (itref_qp_norm[1] < 1e0 * arg->res_b_max | itref_qp_norm[1] < 1e-3 * ws->res->res_max[1]) &
                        (itref_qp_norm[2] < 1e0 * arg->res_d_max | itref_qp_norm[2] < 1e-3 * ws->res->res_max[2]) &
                        (itref_qp_norm[3] < 1e0 * arg->res_m_max | itref_qp_norm[3] < 1e-3 * ws->res->res_max[3])) {
                    break;
                }

                ws->use_Pb = 0;
                d_ocp_qp_solve_kkt_step(ws->qp_itref, ws->sol_itref, arg, ws);
                if (ws->mask_constr) {
                    // mask out disregarded constraints
                    for (ii = 0; ii <= N; ii++)
                        dvecmul(2 * ns[ii], qp->d_mask + ii, 2 * nb[ii] + 2 * ng[ii], ws->sol_step->ux + ii, nu[ii] + nx[ii], ws->sol_step->ux + ii, nu[ii] + nx[ii]);
                    dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);
                    dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
                }
                iter_ref_step = 1;

                for (ii = 0; ii <= N; ii++)
                    daxpy(nu[ii] + nx[ii] + 2 * ns[ii], 1.0, ws->sol_itref->ux + ii, 0, ws->sol_step->ux + ii, 0, ws->sol_step->ux + ii, 0);
                for (ii = 0; ii < N; ii++)
                    daxpy(nx[ii + 1], 1.0, ws->sol_itref->pi + ii, 0, ws->sol_step->pi + ii, 0, ws->sol_step->pi + ii, 0);
                for (ii = 0; ii <= N; ii++)
                    daxpy(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], 1.0, ws->sol_itref->lam + ii, 0, ws->sol_step->lam + ii, 0, ws->sol_step->lam + ii, 0);
                for (ii = 0; ii <= N; ii++)
                    daxpy(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], 1.0, ws->sol_itref->t + ii, 0, ws->sol_step->t + ii, 0, ws->sol_step->t + ii, 0);
            }
            if (itref1 == arg->itref_corr_max) {
                d_ocp_qp_res_compute_lin(ws->qp_step, qp_sol, ws->sol_step, ws->res_itref, ws->res_workspace);
                if (ws->mask_constr) {
                    // mask out disregarded constraints
                    for (ii = 0; ii <= N; ii++)
                        dvecmul(2 * ns[ii], qp->d_mask + ii, 2 * nb[ii] + 2 * ng[ii], ws->res_itref->res_g + ii, nu[ii] + nx[ii], ws->res_itref->res_g + ii, nu[ii] + nx[ii]);
                    dvecmul(cws->nc, qp->d_mask, 0, ws->res_itref->res_d, 0, ws->res_itref->res_d, 0);
                    dvecmul(cws->nc, qp->d_mask, 0, ws->res_itref->res_m, 0, ws->res_itref->res_m, 0);
                }
                d_ocp_qp_res_compute_inf_norm(ws->res_itref);
                itref_qp_norm[0] = ws->res_itref->res_max[0];
                itref_qp_norm[1] = ws->res_itref->res_max[1];
                itref_qp_norm[2] = ws->res_itref->res_max[2];
                itref_qp_norm[3] = ws->res_itref->res_max[3];
                if (kk + 1 < ws->stat_max) {
                    stat[stat_m * (kk + 1) + 14] = itref_qp_norm[0];
                    stat[stat_m * (kk + 1) + 15] = itref_qp_norm[1];
                    stat[stat_m * (kk + 1) + 16] = itref_qp_norm[2];
                    stat[stat_m * (kk + 1) + 17] = itref_qp_norm[3];
                }
            }
        }

        if (iter_ref_step) {
            // alpha
            d_compute_alpha_qp(cws);
            if (kk + 1 < ws->stat_max) {
                stat[stat_m * (kk + 1) + 3] = cws->alpha_prim;
                stat[stat_m * (kk + 1) + 4] = cws->alpha_dual;
            }
        }
    }
    if (arg->itref_corr_max == 0) {
        if (kk + 1 < ws->stat_max) {
            stat[stat_m * (kk + 1) + 14] = 0.0;
            stat[stat_m * (kk + 1) + 15] = 0.0;
            stat[stat_m * (kk + 1) + 16] = 0.0;
            stat[stat_m * (kk + 1) + 17] = 0.0;
        }
    }
    if (kk + 1 < ws->stat_max)
        stat[stat_m * (kk + 1) + 13] = itref1;

    // TODO step length computation

    //
    d_update_var_qp(cws);
    if (ws->mask_constr) {
        // mask out disregarded constraints
        dvecmul(cws->nc, qp->d_mask, 0, qp_sol->lam, 0, qp_sol->lam, 0);
    }
}


void d_ocp_qp_ipm_solve(struct d_ocp_qp* qp, struct d_ocp_qp_sol* qp_sol, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws) {

#if 0
	d_ocp_qp_dim_print(qp->dim);
	d_ocp_qp_print(qp->dim, qp);
#endif

    // dim
    int N = qp->dim->N;
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;
    int* ns = qp->dim->ns;

    struct d_core_qp_ipm_workspace* cws = ws->core_workspace;

    int kk, ii;
    double mu;

    double* stat = ws->stat;
    int stat_m = ws->stat_m;
    double tau_min = arg->tau_min;

    // arg to core workspace
    cws->lam_min = arg->lam_min;
    cws->t_min = arg->t_min;
    cws->t_min_inv = arg->t_min > 0 ? 1.0 / arg->t_min : 1e30;
    cws->tau_min = arg->tau_min;
    cws->split_step = arg->split_step;
    cws->t_lam_min = arg->t_lam_min;

    // alias qp vectors into qp_sol
    cws->v = qp_sol->ux->pa;
    cws->pi = qp_sol->pi->pa;
    cws->lam = qp_sol->lam->pa;
    cws->t = qp_sol->t->pa;

    // alias members of qp_step
    ws->qp_step->dim = qp->dim;
    ws->qp_step->RSQrq = qp->RSQrq;
    ws->qp_step->BAbt = qp->BAbt;
    ws->qp_step->DCt = qp->DCt;
    ws->qp_step->Z = qp->Z;
    ws->qp_step->idxb = qp->idxb;
    ws->qp_step->idxs_rev = qp->idxs_rev;
    ws->qp_step->rqz = ws->res->res_g;
    ws->qp_step->b = ws->res->res_b;
    ws->qp_step->d = ws->res->res_d;
    ws->qp_step->m = ws->res->res_m;

    // alias members of qp_itref
    ws->qp_itref->dim = qp->dim;
    ws->qp_itref->RSQrq = qp->RSQrq;
    ws->qp_itref->BAbt = qp->BAbt;
    ws->qp_itref->DCt = qp->DCt;
    ws->qp_itref->Z = qp->Z;
    ws->qp_itref->idxb = qp->idxb;
    ws->qp_itref->idxs_rev = qp->idxs_rev;
    ws->qp_itref->rqz = ws->res_itref->res_g;
    ws->qp_itref->b = ws->res_itref->res_b;
    ws->qp_itref->d = ws->res_itref->res_d;
    ws->qp_itref->m = ws->res_itref->res_m;

    double* qp_res_max = ws->res->res_max;
    qp_res_max[0] = 0;
    qp_res_max[1] = 0;
    qp_res_max[2] = 0;
    qp_res_max[3] = 0;

    ws->valid_ric_vec = 0;


    // detect constr mask
    int mask_unconstr;
    int nc_mask = 0;
    for (ii = 0; ii < cws->nc; ii++) {
        if (qp->d_mask->pa[ii] != 0.0)
            nc_mask++;
    }
    if (nc_mask < cws->nc) {
        ws->mask_constr = 1;
    } else {
        ws->mask_constr = 0;
    }
    if (nc_mask == 0) {
        mask_unconstr = 1;
        cws->nc_mask = 0;
        cws->nc_mask_inv = 0.0;
    } else {
        mask_unconstr = 0;
        cws->nc_mask = nc_mask;
        cws->nc_mask_inv = 1.0 / nc_mask;
    }


    // no constraints
    if (cws->nc == 0 | mask_unconstr == 1) {
        ws->valid_ric_vec = 1;
        d_ocp_qp_fact_solve_kkt_unconstr(qp, qp_sol, arg, ws);
        if (arg->comp_res_exit) {
            d_ocp_qp_res_compute(qp, qp_sol, ws->res, ws->res_workspace);
            // XXX no constraints, so no mask
            d_ocp_qp_res_compute_inf_norm(ws->res);
            if (0 < ws->stat_max) {
                stat[6] = qp_res_max[0];
                stat[7] = qp_res_max[1];
                stat[8] = qp_res_max[2];
                stat[9] = qp_res_max[3];
                stat[10] = ws->res->obj;
            }
            cws->mu = ws->res->res_mu;
        }
        ws->iter = 0;

        if (isnan(VECEL(qp_sol->ux + 0, 0))) {
            // NaN in the solution
            ws->status = NAN_SOL;
        }
        // #else
        //         if (VECEL(qp_sol->ux + 0, 0) != VECEL(qp_sol->ux + 0, 0)) {
        //             // NaN in the solution
        //             ws->status = NAN_sol;
        //         }

        else {
            // normal return
            ws->status = SUCCESS;
        }
    }


    // init solver
    d_ocp_qp_init_var(qp, qp_sol, arg, ws);
    if (ws->mask_constr) {
        // mask out disregarded constraints
        dvecmul(cws->nc, qp->d_mask, 0, qp_sol->lam, 0, qp_sol->lam, 0);
    }

    cws->alpha = 1.0;


    // absolute IPM formulation

    if (arg->abs_form) {

        // alias members of qp_step
        ws->qp_step->dim = qp->dim;
        ws->qp_step->RSQrq = qp->RSQrq;
        ws->qp_step->BAbt = qp->BAbt;
        ws->qp_step->DCt = qp->DCt;
        ws->qp_step->Z = qp->Z;
        ws->qp_step->idxb = qp->idxb;
        ws->qp_step->idxe = qp->idxe;
        ws->qp_step->idxs_rev = qp->idxs_rev;
        ws->qp_step->rqz = qp->rqz;
        ws->qp_step->b = qp->b;
        ws->qp_step->d = qp->d;
        ws->qp_step->d_mask = qp->d_mask;
        ws->qp_step->m = ws->tmp_m;

        // d_ocp_qp_dim_print(ws->qp_step->dim);
        // d_ocp_qp_print(ws->qp_step->dim, ws->qp_step);
        // exit(1);
        //  alias core workspace
        cws->res_m = ws->qp_step->m->pa;
        cws->res_m_bkp = ws->qp_step->m->pa;  // TODO remove (as in dense qp) ???

        ws->valid_ric_vec = 1;


        mu = dvecmuldot(cws->nc, qp_sol->lam, 0, qp_sol->t, 0, ws->tmp_m, 0);
        mu /= cws->nc;
        cws->mu = mu;

        // IPM loop (absolute formulation)
        for (kk = 0;
             kk<arg->iter_max &
                cws->alpha>
                     arg->alpha_min &
             fabs(mu - tau_min) > arg->res_m_max;
             kk++) {

            // compute delta step
            d_ocp_qp_ipm_abs_step(kk, qp, qp_sol, arg, ws);

            // compute mu
            mu = dvecmuldot(cws->nc, qp_sol->lam, 0, qp_sol->t, 0, ws->tmp_m, 0);
            mu /= cws->nc;
            cws->mu = mu;
            if (kk + 1 < ws->stat_max)
                stat[stat_m * (kk + 1) + 5] = mu;
        }

        if (arg->comp_res_exit & arg->comp_dual_sol_eq) {
            // compute residuals
            d_ocp_qp_res_compute(qp, qp_sol, ws->res, ws->res_workspace);
            if (ws->mask_constr) {
                // mask out disregarded constraints
                for (ii = 0; ii <= N; ii++)
                    dvecmul(2 * ns[ii], qp->d_mask + ii, 2 * nb[ii] + 2 * ng[ii], ws->res->res_g + ii, nu[ii] + nx[ii], ws->res->res_g + ii, nu[ii] + nx[ii]);
                dvecmul(cws->nc, qp->d_mask, 0, ws->res->res_d, 0, ws->res->res_d, 0);
                dvecmul(cws->nc, qp->d_mask, 0, ws->res->res_m, 0, ws->res->res_m, 0);
            }
            d_ocp_qp_res_compute_inf_norm(ws->res);
            // save infinity norm of residuals
            // XXX it is already kk+1
            if (kk < ws->stat_max) {
                stat[stat_m * (kk + 0) + 6] = qp_res_max[0];
                stat[stat_m * (kk + 0) + 7] = qp_res_max[1];
                stat[stat_m * (kk + 0) + 8] = qp_res_max[2];
                stat[stat_m * (kk + 0) + 9] = qp_res_max[3];
                stat[stat_m * (kk + 0) + 10] = ws->res->obj;
            }
        }

        goto set_status;
    }


    // compute residuals
    d_ocp_qp_res_compute(qp, qp_sol, ws->res, ws->res_workspace);
    if (ws->mask_constr) {
        // mask out disregarded constraints
        for (ii = 0; ii <= N; ii++)
            dvecmul(2 * ns[ii], qp->d_mask + ii, 2 * nb[ii] + 2 * ng[ii], ws->res->res_g + ii, nu[ii] + nx[ii], ws->res->res_g + ii, nu[ii] + nx[ii]);
        dvecmul(cws->nc, qp->d_mask, 0, ws->res->res_d, 0, ws->res->res_d, 0);
        dvecmul(cws->nc, qp->d_mask, 0, ws->res->res_m, 0, ws->res->res_m, 0);
    }
    cws->mu = ws->res->res_mu;
    d_ocp_qp_res_compute_inf_norm(ws->res);
    // save infinity norm of residuals
    if (0 < ws->stat_max) {
        stat[stat_m * (0) + 6] = qp_res_max[0];
        stat[stat_m * (0) + 7] = qp_res_max[1];
        stat[stat_m * (0) + 8] = qp_res_max[2];
        stat[stat_m * (0) + 9] = qp_res_max[3];
        stat[stat_m * (0) + 10] = ws->res->obj;
    }


    // relative (delta) IPM formulation

    // IPM loop
    for (kk = 0; kk<arg->iter_max &
                    cws->alpha>
                         arg->alpha_min &
         (qp_res_max[0] > arg->res_g_max |
          qp_res_max[1] > arg->res_b_max |
          qp_res_max[2] > arg->res_d_max |
          fabs(qp_res_max[3] - tau_min) > arg->res_m_max);
         kk++) {

        // compute delta step
        d_ocp_qp_ipm_delta_step(kk, qp, qp_sol, arg, ws);

        // compute residuals
        d_ocp_qp_res_compute(qp, qp_sol, ws->res, ws->res_workspace);
        if (ws->mask_constr) {
            // mask out disregarded constraints
            for (ii = 0; ii <= N; ii++)
                dvecmul(2 * ns[ii], qp->d_mask + ii, 2 * nb[ii] + 2 * ng[ii], ws->res->res_g + ii, nu[ii] + nx[ii], ws->res->res_g + ii, nu[ii] + nx[ii]);
            dvecmul(cws->nc, qp->d_mask, 0, ws->res->res_d, 0, ws->res->res_d, 0);
            dvecmul(cws->nc, qp->d_mask, 0, ws->res->res_m, 0, ws->res->res_m, 0);
        }
        cws->mu = ws->res->res_mu;
        d_ocp_qp_res_compute_inf_norm(ws->res);
        // save infinity norm of residuals
        if (kk + 1 < ws->stat_max) {
            stat[stat_m * (kk + 1) + 5] = ws->res->res_mu;
            stat[stat_m * (kk + 1) + 6] = qp_res_max[0];
            stat[stat_m * (kk + 1) + 7] = qp_res_max[1];
            stat[stat_m * (kk + 1) + 8] = qp_res_max[2];
            stat[stat_m * (kk + 1) + 9] = qp_res_max[3];
            stat[stat_m * (kk + 1) + 10] = ws->res->obj;
        }
    }

set_status:

    // save info before return
    ws->iter = kk;

    if (kk == arg->iter_max) {
        // max iteration number reached
        ws->status = MAX_ITER;
    } else if (cws->alpha <= arg->alpha_min) {
        // min step lenght
        ws->status = MIN_STEP;
    }

    else if (isnan(cws->mu)) {
        // NaN in the solution
        ws->status = NAN_SOL;
    }
    // #else
    //     else if ((cws->mu != cws->mu)) {
    //         // NaN in the solution
    //         ws->status = NAN_sol;
    //     }

    else {
        // normal return
        ws->status = SUCCESS;
    }

call_return:

    // TODO compute objective

    // return
}


void d_ocp_qp_ipm_predict(struct d_ocp_qp* qp, struct d_ocp_qp_sol* qp_sol, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws) {

#if 0
	d_ocp_qp_dim_print(qp->dim);
	d_ocp_qp_print(qp->dim, qp);
#endif

    int ii;

    struct d_core_qp_ipm_workspace* cws = ws->core_workspace;

    // arg to core workspace
    cws->lam_min = arg->lam_min;
    cws->t_min = arg->t_min;
    cws->t_min_inv = arg->t_min > 0 ? 1.0 / arg->t_min : 1e30;
    cws->t_lam_min = arg->t_lam_min;

    // alias qp vectors into qp_sol
    cws->v = qp_sol->ux->pa;
    cws->pi = qp_sol->pi->pa;
    cws->lam = qp_sol->lam->pa;
    cws->t = qp_sol->t->pa;

    // load sol from bkp
    for (ii = 0; ii < cws->nv; ii++)
        cws->v[ii] = cws->v_bkp[ii];
    for (ii = 0; ii < cws->ne; ii++)
        cws->pi[ii] = cws->pi_bkp[ii];
    for (ii = 0; ii < cws->nc; ii++)
        cws->lam[ii] = cws->lam_bkp[ii];
    for (ii = 0; ii < cws->nc; ii++)
        cws->t[ii] = cws->t_bkp[ii];

    if (arg->abs_form) {
        // solve kkt
        ws->use_Pb = 0;
        d_ocp_qp_solve_kkt_step(qp, qp_sol, arg, ws);

        if (arg->comp_res_exit & arg->comp_dual_sol_eq) {
            // blasfeo alias for residuals
            struct vec str_res_g;
            struct vec str_res_b;
            struct vec str_res_d;
            struct vec str_res_m;
            str_res_g.m = cws->nv;
            str_res_b.m = cws->ne;
            str_res_d.m = cws->nc;
            str_res_m.m = cws->nc;
            str_res_g.pa = cws->res_g;
            str_res_b.pa = cws->res_b;
            str_res_d.pa = cws->res_d;
            str_res_m.pa = cws->res_m;

            double* qp_res_max = ws->res->res_max;
            qp_res_max[0] = 0;
            qp_res_max[1] = 0;
            qp_res_max[2] = 0;
            qp_res_max[3] = 0;

            // compute residuals
            d_ocp_qp_res_compute(qp, qp_sol, ws->res, ws->res_workspace);

            // TODO mask out disregarded constraints ???

            // compute infinity norm of residuals
            dvecnrm_inf(cws->nv, &str_res_g, 0, &qp_res_max[0]);
            dvecnrm_inf(cws->ne, &str_res_b, 0, &qp_res_max[1]);
            dvecnrm_inf(cws->nc, &str_res_d, 0, &qp_res_max[2]);
            dvecnrm_inf(cws->nc, &str_res_m, 0, &qp_res_max[3]);
        }
    }

    // alias members of qp_step
    ws->qp_step->dim = qp->dim;
    ws->qp_step->RSQrq = qp->RSQrq;
    ws->qp_step->BAbt = qp->BAbt;
    ws->qp_step->DCt = qp->DCt;
    ws->qp_step->Z = qp->Z;
    ws->qp_step->idxb = qp->idxb;
    ws->qp_step->idxs_rev = qp->idxs_rev;
    ws->qp_step->rqz = ws->res->res_g;
    ws->qp_step->b = ws->res->res_b;
    ws->qp_step->d = ws->res->res_d;
    ws->qp_step->m = ws->res->res_m;

    // TODO ?

    // blasfeo alias for residuals
    struct vec str_res_g;
    struct vec str_res_b;
    struct vec str_res_d;
    struct vec str_res_m;
    str_res_g.m = cws->nv;
    str_res_b.m = cws->ne;
    str_res_d.m = cws->nc;
    str_res_m.m = cws->nc;
    str_res_g.pa = cws->res_g;
    str_res_b.pa = cws->res_b;
    str_res_d.pa = cws->res_d;
    str_res_m.pa = cws->res_m;

    double* qp_res_max = ws->res->res_max;
    qp_res_max[0] = 0;
    qp_res_max[1] = 0;
    qp_res_max[2] = 0;
    qp_res_max[3] = 0;

    // TODO robust formulation !!!!!

    // compute residuals
    d_ocp_qp_res_compute(qp, qp_sol, ws->res, ws->res_workspace);

    // compute infinity norm of residuals
    dvecnrm_inf(cws->nv, &str_res_g, 0, &qp_res_max[0]);
    dvecnrm_inf(cws->ne, &str_res_b, 0, &qp_res_max[1]);
    dvecnrm_inf(cws->nc, &str_res_d, 0, &qp_res_max[2]);
    dvecnrm_inf(cws->nc, &str_res_m, 0, &qp_res_max[3]);

    // printf("\npredict\t%e\t%e\t%e\t%e\n", qp_res_max[0], qp_res_max[1], qp_res_max[2], qp_res_max[3]);

    // solve kkt
    ws->use_Pb = 0;
    d_ocp_qp_solve_kkt_step(ws->qp_step, ws->sol_step, arg, ws);

    // alpha TODO fix alpha=1 !!!!!
    //	d_compute_alpha_qp(cws);
    cws->alpha = 1.0;

    //
    d_update_var_qp(cws);

    if (arg->comp_res_pred) {
        // compute residuals in exit
        d_ocp_qp_res_compute(qp, qp_sol, ws->res, ws->res_workspace);

        // TODO mask out disregarded constraints ???

        // compute infinity norm of residuals
        dvecnrm_inf(cws->nv, &str_res_g, 0, &qp_res_max[0]);
        dvecnrm_inf(cws->ne, &str_res_b, 0, &qp_res_max[1]);
        dvecnrm_inf(cws->nc, &str_res_d, 0, &qp_res_max[2]);
        dvecnrm_inf(cws->nc, &str_res_m, 0, &qp_res_max[3]);
    }

    // printf("\npredict\t%e\t%e\t%e\t%e\n", qp_res_max[0], qp_res_max[1], qp_res_max[2], qp_res_max[3]);

    // TODO

    // do not change status
}


void d_ocp_qp_ipm_sens(struct d_ocp_qp* qp, struct d_ocp_qp_sol* qp_sol, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws) {

#if 0
	d_ocp_qp_dim_print(qp->dim);
	d_ocp_qp_print(qp->dim, qp);
#endif

    int ii;

    struct d_core_qp_ipm_workspace* cws = ws->core_workspace;

    // arg to core workspace
    cws->lam_min = arg->lam_min;
    cws->t_min = arg->t_min;
    cws->t_min_inv = arg->t_min > 0 ? 1.0 / arg->t_min : 1e30;
    cws->t_lam_min = arg->t_lam_min;

    // alias qp vectors into qp_sol
    cws->v = qp_sol->ux->pa;
    cws->pi = qp_sol->pi->pa;
    cws->lam = qp_sol->lam->pa;
    cws->t = qp_sol->t->pa;

#if 0
	// alias members of qp_step
	ws->qp_step->dim = qp->dim;
	ws->qp_step->RSQrq = qp->RSQrq;
	ws->qp_step->BAbt = qp->BAbt;
	ws->qp_step->DCt = qp->DCt;
	ws->qp_step->Z = qp->Z;
	ws->qp_step->idxb = qp->idxb;
	ws->qp_step->idxs_rev = qp->idxs_rev;
	ws->qp_step->rqz = res->res_g;
	ws->qp_step->b = res->res_b;
	ws->qp_step->d = res->res_d;
	ws->qp_step->m = res->res_m;
#endif

    // TODO ?

#if 0
	// blasfeo alias for residuals
	struct vec str_res_g;
	struct vec str_res_b;
	struct vec str_res_d;
	struct vec str_res_m;
	str_res_g.m = cws->nv;
	str_res_b.m = cws->ne;
	str_res_d.m = cws->nc;
	str_res_m.m = cws->nc;
	str_res_g.pa = cws->res_g;
	str_res_b.pa = cws->res_b;
	str_res_d.pa = cws->res_d;
	str_res_m.pa = cws->res_m;

	double *qp_res_max = ws->qp_res_max;
	qp_res_max[0] = 0;
	qp_res_max[1] = 0;
	qp_res_max[2] = 0;
	qp_res_max[3] = 0;
#endif

#if 0
// compute residuals
d_ocp_qp_res_compute(qp, qp_sol, ws->res, ws->res_workspace);

// compute infinity norm of residuals
dvecnrm_inf(cws->nv, &str_res_g, 0, &qp_res_max[0]);
dvecnrm_inf(cws->ne, &str_res_b, 0, &qp_res_max[1]);
dvecnrm_inf(cws->nc, &str_res_d, 0, &qp_res_max[2]);
dvecnrm_inf(cws->nc, &str_res_m, 0, &qp_res_max[3]);

printf("\npredict\t%e\t%e\t%e\t%e\n", qp_res_max[0], qp_res_max[1], qp_res_max[2], qp_res_max[3]);
#endif

    // load sol from bkp
    for (ii = 0; ii < cws->nv; ii++)
        cws->v[ii] = cws->v_bkp[ii];
    for (ii = 0; ii < cws->ne; ii++)
        cws->pi[ii] = cws->pi_bkp[ii];
    for (ii = 0; ii < cws->nc; ii++)
        cws->lam[ii] = cws->lam_bkp[ii];
    for (ii = 0; ii < cws->nc; ii++)
        cws->t[ii] = cws->t_bkp[ii];

    // printf("\npredict\t%e\t%e\t%e\t%e\n", qp_res_max[0], qp_res_max[1], qp_res_max[2], qp_res_max[3]);

    // solve kkt
    ws->use_Pb = 0;
    //	d_ocp_qp_solve_kkt_step(ws->qp_step, ws->sol_step, arg, ws);
    d_ocp_qp_solve_kkt_step(qp, qp_sol, arg, ws);

#if 0
	// alpha TODO fix alpha=1 !!!!!
//	d_compute_alpha_qp(cws);
	cws->alpha = 1.0;

	//
	d_update_var_qp(cws);

	if(arg->comp_res_pred)
		{
		// compute residuals in exit
		d_ocp_qp_res_compute(qp, qp_sol, ws->res, ws->res_workspace);

		// compute infinity norm of residuals
		dvecnrm_inf(cws->nv, &str_res_g, 0, &qp_res_max[0]);
		dvecnrm_inf(cws->ne, &str_res_b, 0, &qp_res_max[1]);
		dvecnrm_inf(cws->nc, &str_res_d, 0, &qp_res_max[2]);
		dvecnrm_inf(cws->nc, &str_res_m, 0, &qp_res_max[3]);
		}
#endif

    // printf("\npredict\t%e\t%e\t%e\t%e\n", qp_res_max[0], qp_res_max[1], qp_res_max[2], qp_res_max[3]);

    // TODO

    // do not change status
}
