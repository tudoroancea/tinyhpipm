#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "tinyhpipm/blas.h"
#include "tinyhpipm/common.h"
#include "tinyhpipm/dense/d_dense_qp.h"
#include "tinyhpipm/dense/d_dense_qp_dim.h"
#include "tinyhpipm/dense/d_dense_qp_ipm.h"
#include "tinyhpipm/dense/d_dense_qp_kkt.h"
#include "tinyhpipm/dense/d_dense_qp_res.h"
#include "tinyhpipm/dense/d_dense_qp_sol.h"
#include "tinyhpipm/dense/d_dense_qp_utils.h"
#include "tinyhpipm/ipm_core/d_core_qp_ipm.h"
#include "tinyhpipm/ipm_core/d_core_qp_ipm_aux.h"

hpipm_size_t d_dense_qp_ipm_arg_strsize() {
    return sizeof(struct d_dense_qp_ipm_arg);
}


hpipm_size_t d_dense_qp_ipm_arg_memsize(struct d_dense_qp_dim* dim) {
    return 0;
}


void d_dense_qp_ipm_arg_create(struct d_dense_qp_dim* dim, struct d_dense_qp_ipm_arg* arg, void* mem) {

    // zero memory (to avoid corrupted memory like e.g. NaN)
    //	hpipm_size_t memsize = d_dense_qp_ipm_arg_memsize(dim);
    //	hpipm_zero_memset(memsize, mem);

    arg->memsize = 0;
}


void d_dense_qp_ipm_arg_set_default(enum tinyhpipm_mode mode, struct d_dense_qp_ipm_arg* arg) {

    double mu0, alpha_min, res_g, res_b, res_d, res_m, reg_prim, reg_dual, lam_min, t_min, tau_min;
    int iter_max, stat_max, pred_corr, cond_pred_corr, itref_pred_max, itref_corr_max, lq_fact, scale, warm_start, abs_form, comp_res_exit, comp_res_pred, kkt_fact_alg, remove_lin_dep_eq, compute_obj, split_step, t_lam_min;

    if (mode == SPEED_ABS) {
        mu0 = 1e1;
        alpha_min = 1e-12;
        res_g = 1e0;
        res_b = 1e0;
        res_d = 1e0;
        res_m = 1e-8;
        iter_max = 15;
        stat_max = 15;
        pred_corr = 1;
        cond_pred_corr = 1;
        itref_pred_max = 0;
        itref_corr_max = 0;
        reg_prim = 1e-15;
        reg_dual = 1e-15;
        lq_fact = 0;
        scale = 0;
        lam_min = 1e-16;
        t_min = 1e-16;
        tau_min = 1e-16;
        warm_start = 0;
        abs_form = 1;
        comp_res_exit = 0;
        comp_res_pred = 0;
        kkt_fact_alg = 1;
        remove_lin_dep_eq = 0;
        compute_obj = 0;
        split_step = 1;
        t_lam_min = 2;
    } else if (mode == SPEED) {
        mu0 = 1e1;
        alpha_min = 1e-12;
        res_g = 1e-6;
        res_b = 1e-8;
        res_d = 1e-8;
        res_m = 1e-8;
        iter_max = 15;
        stat_max = 15;
        pred_corr = 1;
        cond_pred_corr = 1;
        itref_pred_max = 0;
        itref_corr_max = 0;
        reg_prim = 1e-15;
        reg_dual = 1e-15;
        lq_fact = 0;
        scale = 0;
        lam_min = 1e-16;
        t_min = 1e-16;
        tau_min = 1e-16;
        warm_start = 0;
        abs_form = 0;
        comp_res_exit = 1;
        comp_res_pred = 0;
        kkt_fact_alg = 1;
        remove_lin_dep_eq = 0;
        compute_obj = 0;
        split_step = 1;
        t_lam_min = 2;
    } else if (mode == BALANCE) {
        mu0 = 1e1;
        alpha_min = 1e-12;
        res_g = 1e-6;
        res_b = 1e-8;
        res_d = 1e-8;
        res_m = 1e-8;
        iter_max = 30;
        stat_max = 30;
        pred_corr = 1;
        cond_pred_corr = 1;
        itref_pred_max = 0;
        itref_corr_max = 2;
        reg_prim = 1e-15;
        reg_dual = 1e-15;
        lq_fact = 1;
        scale = 0;
        lam_min = 1e-16;
        t_min = 1e-16;
        tau_min = 1e-16;
        warm_start = 0;
        abs_form = 0;
        comp_res_exit = 1;
        comp_res_pred = 0;
        kkt_fact_alg = 1;
        remove_lin_dep_eq = 0;
        compute_obj = 0;
        split_step = 0;
        t_lam_min = 2;
    } else if (mode == ROBUST) {
        mu0 = 1e2;
        alpha_min = 1e-12;
        res_g = 1e-6;
        res_b = 1e-8;
        res_d = 1e-8;
        res_m = 1e-8;
        iter_max = 100;
        stat_max = 100;
        pred_corr = 1;
        cond_pred_corr = 1;
        itref_pred_max = 0;
        itref_corr_max = 4;
        reg_prim = 1e-15;
        reg_dual = 1e-15;
        lq_fact = 2;
        scale = 0;
        lam_min = 1e-16;
        t_min = 1e-16;
        tau_min = 1e-16;
        warm_start = 0;
        abs_form = 0;
        comp_res_exit = 1;
        comp_res_pred = 0;
        kkt_fact_alg = 1;
        remove_lin_dep_eq = 0;
        compute_obj = 0;
        split_step = 0;
        t_lam_min = 2;
    } else {
        printf("\nerror: d_dense_qp_ipm_arg_set_default: wrong set default mode\n");
        exit(1);
    }

    // use individual setters when available
    d_dense_qp_ipm_arg_set_mu0(&mu0, arg);
    d_dense_qp_ipm_arg_set_alpha_min(&alpha_min, arg);
    d_dense_qp_ipm_arg_set_tol_stat(&res_g, arg);
    d_dense_qp_ipm_arg_set_tol_eq(&res_b, arg);
    d_dense_qp_ipm_arg_set_tol_ineq(&res_d, arg);
    d_dense_qp_ipm_arg_set_tol_comp(&res_m, arg);
    d_dense_qp_ipm_arg_set_iter_max(&iter_max, arg);
    arg->stat_max = stat_max;
    d_dense_qp_ipm_arg_set_pred_corr(&pred_corr, arg);
    d_dense_qp_ipm_arg_set_cond_pred_corr(&cond_pred_corr, arg);
    arg->itref_pred_max = itref_pred_max;
    arg->itref_corr_max = itref_corr_max;
    d_dense_qp_ipm_arg_set_reg_prim(&reg_prim, arg);
    d_dense_qp_ipm_arg_set_reg_dual(&reg_prim, arg);
    arg->lq_fact = lq_fact;
    arg->scale = scale;
    d_dense_qp_ipm_arg_set_lam_min(&lam_min, arg);
    d_dense_qp_ipm_arg_set_t_min(&t_min, arg);
    d_dense_qp_ipm_arg_set_tau_min(&tau_min, arg);
    d_dense_qp_ipm_arg_set_warm_start(&warm_start, arg);
    arg->abs_form = abs_form;
    d_dense_qp_ipm_arg_set_comp_res_exit(&comp_res_exit, arg);
    d_dense_qp_ipm_arg_set_comp_res_pred(&comp_res_pred, arg);
    d_dense_qp_ipm_arg_set_kkt_fact_alg(&kkt_fact_alg, arg);
    arg->mode = mode;
    d_dense_qp_ipm_arg_set_remove_lin_dep_eq(&remove_lin_dep_eq, arg);
    d_dense_qp_ipm_arg_set_compute_obj(&compute_obj, arg);
    d_dense_qp_ipm_arg_set_split_step(&split_step, arg);
    d_dense_qp_ipm_arg_set_t_lam_min(&t_lam_min, arg);
}


void d_dense_qp_ipm_arg_set(char* field, void* value, struct d_dense_qp_ipm_arg* arg) {
    if (hpipm_strcmp(field, "iter_max")) {
        d_dense_qp_ipm_arg_set_iter_max(value, arg);
    } else if (hpipm_strcmp(field, "alpha_min")) {
        d_dense_qp_ipm_arg_set_alpha_min(value, arg);
    } else if (hpipm_strcmp(field, "mu0")) {
        d_dense_qp_ipm_arg_set_mu0(value, arg);
    } else if (hpipm_strcmp(field, "tol_stat")) {
        d_dense_qp_ipm_arg_set_tol_stat(value, arg);
    } else if (hpipm_strcmp(field, "tol_eq")) {
        d_dense_qp_ipm_arg_set_tol_eq(value, arg);
    } else if (hpipm_strcmp(field, "tol_ineq")) {
        d_dense_qp_ipm_arg_set_tol_ineq(value, arg);
    } else if (hpipm_strcmp(field, "tol_comp")) {
        d_dense_qp_ipm_arg_set_tol_comp(value, arg);
    } else if (hpipm_strcmp(field, "reg_prim")) {
        d_dense_qp_ipm_arg_set_reg_prim(value, arg);
    } else if (hpipm_strcmp(field, "reg_dual")) {
        d_dense_qp_ipm_arg_set_reg_dual(value, arg);
    } else if (hpipm_strcmp(field, "warm_start")) {
        d_dense_qp_ipm_arg_set_warm_start(value, arg);
    } else if (hpipm_strcmp(field, "pred_corr")) {
        d_dense_qp_ipm_arg_set_pred_corr(value, arg);
    } else if (hpipm_strcmp(field, "cond_pred_corr")) {
        d_dense_qp_ipm_arg_set_cond_pred_corr(value, arg);
    } else if (hpipm_strcmp(field, "comp_res_exit")) {
        d_dense_qp_ipm_arg_set_comp_res_exit(value, arg);
    } else if (hpipm_strcmp(field, "comp_res_pred")) {
        d_dense_qp_ipm_arg_set_comp_res_pred(value, arg);
    } else if (hpipm_strcmp(field, "lam_min")) {
        d_dense_qp_ipm_arg_set_lam_min(value, arg);
    } else if (hpipm_strcmp(field, "t_min")) {
        d_dense_qp_ipm_arg_set_t_min(value, arg);
    } else if (hpipm_strcmp(field, "tau_min")) {
        d_dense_qp_ipm_arg_set_tau_min(value, arg);
    } else if (hpipm_strcmp(field, "kkt_fact_alg")) {
        d_dense_qp_ipm_arg_set_kkt_fact_alg(value, arg);
    } else if (hpipm_strcmp(field, "remove_lin_dep_eq")) {
        d_dense_qp_ipm_arg_set_remove_lin_dep_eq(value, arg);
    } else if (hpipm_strcmp(field, "compute_obj")) {
        d_dense_qp_ipm_arg_set_compute_obj(value, arg);
    } else if (hpipm_strcmp(field, "split_step")) {
        d_dense_qp_ipm_arg_set_split_step(value, arg);
    } else if (hpipm_strcmp(field, "t_lam_min")) {
        d_dense_qp_ipm_arg_set_t_lam_min(value, arg);
    } else {
        printf("error: d_dense_qp_ipm_arg_set: wrong field %s\n", field);
        exit(1);
    }
}


void d_dense_qp_ipm_arg_set_iter_max(int* iter_max, struct d_dense_qp_ipm_arg* arg) {
    arg->iter_max = *iter_max;
}


void d_dense_qp_ipm_arg_set_alpha_min(double* alpha_min, struct d_dense_qp_ipm_arg* arg) {
    arg->alpha_min = *alpha_min;
}


void d_dense_qp_ipm_arg_set_mu0(double* mu0, struct d_dense_qp_ipm_arg* arg) {
    arg->mu0 = *mu0;
}


void d_dense_qp_ipm_arg_set_tol_stat(double* tol_stat, struct d_dense_qp_ipm_arg* arg) {
    arg->res_g_max = *tol_stat;
}


void d_dense_qp_ipm_arg_set_tol_eq(double* tol_eq, struct d_dense_qp_ipm_arg* arg) {
    arg->res_b_max = *tol_eq;
}


void d_dense_qp_ipm_arg_set_tol_ineq(double* tol_ineq, struct d_dense_qp_ipm_arg* arg) {
    arg->res_d_max = *tol_ineq;
}


void d_dense_qp_ipm_arg_set_tol_comp(double* tol_comp, struct d_dense_qp_ipm_arg* arg) {
    arg->res_m_max = *tol_comp;
}


void d_dense_qp_ipm_arg_set_reg_prim(double* reg, struct d_dense_qp_ipm_arg* arg) {
    arg->reg_prim = *reg;
}


void d_dense_qp_ipm_arg_set_reg_dual(double* reg, struct d_dense_qp_ipm_arg* arg) {
    arg->reg_dual = *reg;
}


void d_dense_qp_ipm_arg_set_warm_start(int* warm_start, struct d_dense_qp_ipm_arg* arg) {
    arg->warm_start = *warm_start;
}


void d_dense_qp_ipm_arg_set_pred_corr(int* pred_corr, struct d_dense_qp_ipm_arg* arg) {
    arg->pred_corr = *pred_corr;
}


void d_dense_qp_ipm_arg_set_cond_pred_corr(int* cond_pred_corr, struct d_dense_qp_ipm_arg* arg) {
    arg->cond_pred_corr = *cond_pred_corr;
}


void d_dense_qp_ipm_arg_set_comp_res_pred(int* comp_res_pred, struct d_dense_qp_ipm_arg* arg) {
    arg->comp_res_pred = *comp_res_pred;
}


void d_dense_qp_ipm_arg_set_comp_res_exit(int* comp_res_exit, struct d_dense_qp_ipm_arg* arg) {
    arg->comp_res_exit = *comp_res_exit;
}


void d_dense_qp_ipm_arg_set_lam_min(double* value, struct d_dense_qp_ipm_arg* arg) {
    arg->lam_min = *value;
}


void d_dense_qp_ipm_arg_set_t_min(double* value, struct d_dense_qp_ipm_arg* arg) {
    arg->t_min = *value;
}


void d_dense_qp_ipm_arg_set_tau_min(double* value, struct d_dense_qp_ipm_arg* arg) {
    arg->tau_min = *value;
}


void d_dense_qp_ipm_arg_set_kkt_fact_alg(int* value, struct d_dense_qp_ipm_arg* arg) {
    arg->kkt_fact_alg = *value;
}


void d_dense_qp_ipm_arg_set_remove_lin_dep_eq(int* value, struct d_dense_qp_ipm_arg* arg) {
    arg->remove_lin_dep_eq = *value;
}


void d_dense_qp_ipm_arg_set_compute_obj(int* value, struct d_dense_qp_ipm_arg* arg) {
    arg->compute_obj = *value;
}


void d_dense_qp_ipm_arg_set_split_step(int* value, struct d_dense_qp_ipm_arg* arg) {
    arg->split_step = *value;
}


void d_dense_qp_ipm_arg_set_t_lam_min(int* value, struct d_dense_qp_ipm_arg* arg) {
    arg->t_lam_min = *value;
}


hpipm_size_t d_dense_qp_ipm_ws_strsize() {
    return sizeof(struct d_dense_qp_ipm_ws);
}


hpipm_size_t d_dense_qp_ipm_ws_memsize(struct d_dense_qp_dim* dim, struct d_dense_qp_ipm_arg* arg) {

    int nv = dim->nv;
    int ne = dim->ne;
    int nb = dim->nb;
    int ng = dim->ng;
    int ns = dim->ns;

    hpipm_size_t size = 0;

    size += 1 * sizeof(struct d_core_qp_ipm_workspace);
    size += 1 * d_memsize_core_qp_ipm(nv + 2 * ns, ne, 2 * nb + 2 * ng + 2 * ns);

    size += 1 * sizeof(struct d_dense_qp_res_ws);  // res_ws

    size += 2 * sizeof(struct d_dense_qp);  // qp_step qp_itref

    size += 2 * sizeof(struct d_dense_qp_sol);  // sol_step sol_itref
    size += 1 * d_dense_qp_sol_memsize(dim);  // sol_itref

    size += 3 * sizeof(struct d_dense_qp_res);  // res res_itref res_step
    size += 2 * d_dense_qp_res_memsize(dim);  // res_itref res_step

    size += 27 * sizeof(struct vec);  // sol_step(v,pi,lam,t) res_g res_b res_d res_m lv (4+2)*tmp_nbg (1+1)*tmp_ns Gamma gamma Zs_inv sv se tmp_m tmp_nv tmp_2ns tmp_nv2ns
    size += 5 * sizeof(struct mat);  // 2*Lv AL Le Ctx
    if (arg->lq_fact > 0) {
        size += 2 * sizeof(struct mat);  // lq0 lq1
    }
    if (arg->kkt_fact_alg == 0) {
        size += 5 * sizeof(struct mat);  // A_LQ A_Q Zt ZtH ZthZ
        size += 3 * sizeof(struct vec);  // xy Yxy xz
    } else {
        //		TODO
    }
    if (arg->remove_lin_dep_eq != 0) {
        size += 2 * sizeof(struct mat);  // A_li Ab_LU
        size += 1 * sizeof(struct vec);  // b_li
    }

    size += 4 * memsize_vec(nb + ng);  // 4*tmp_nbg
    size += 1 * memsize_vec(ns);  // tmp_ns
    size += 2 * memsize_vec(nv);  // lv sv
    size += 1 * memsize_vec(ne);  // se
    size += 1 * memsize_vec(2 * ns);  // Zs_inv
    size += 2 * memsize_mat(nv + 1, nv);  // Lv
    size += 1 * memsize_mat(ne, nv);  // AL
    size += 1 * memsize_mat(ne, ne);  // Le
    size += 1 * memsize_mat(nv + 1, ng);  // Ctx
    size += 1 * memsize_vec(2 * nb + 2 * ng + 2 * ns);  // tmp_m
    size += ne > 0 ? 1 * dgelqf_worksize(ne, nv) : 0;  // lq_work_null
    size += 1 * memsize_vec(nv);  // tmp_nv
    size += 1 * memsize_vec(2 * ns);  // tmp_2ns
    size += 2 * memsize_vec(nv + 2 * ns);  // tmp_nv2ns
    if (arg->lq_fact > 0) {
        size += 1 * memsize_mat(ne, ne + nv);  // lq0
        size += 1 * memsize_mat(nv, nv + nv + ng);  // lq1
    }
    if (arg->kkt_fact_alg == 0) {
        size += 1 * memsize_mat(ne, nv);  // A_LQ
        size += 1 * memsize_mat(nv, nv);  // A_Q
        size += 1 * memsize_vec(nv);  // Yxy
        size += 1 * memsize_vec(ne);  // xy
        if (arg->remove_lin_dep_eq != 0) {
            size += 2 * memsize_mat(nv, nv);  // Zt ZtH
            size += 1 * memsize_mat(nv, nv);  // ZtHZ
            size += 1 * memsize_vec(nv);  // xz
        } else {
            size += 2 * memsize_mat(nv - ne, nv);  // Zt ZtH
            size += 1 * memsize_mat(nv - ne, nv - ne);  // ZtHZ
            size += 1 * memsize_vec(nv - ne);  // xz
        }
    } else {
        // TODO
    }
    if (arg->remove_lin_dep_eq != 0) {
        size += 1 * memsize_mat(ne, nv);  // A_li
        size += 1 * memsize_vec(ne);  // b_li
        size += 1 * memsize_mat(ne, nv + 1);  // Ab_LU
    }

    if (arg->lq_fact > 0) {
        size += ne > 0 ? 1 * dgelqf_worksize(ne, nv) : 0;  // lq_work0
        size += 1 * dgelqf_worksize(nv, nv + nv + ng);  // lq_work1
    }
    if (arg->kkt_fact_alg == 0) {
        size += 1 * dorglq_worksize(nv, nv, ne);  // orglq_work_null
    } else {
        // TODO
    }

    // cache value of stat_max used for memory allocation
    if (arg->stat_max < arg->iter_max)
        arg->stat_max = arg->iter_max;

    int stat_m = 18;
    size += stat_m * (1 + arg->stat_max) * sizeof(double);  // stat

    size += nv * sizeof(int);  // ipiv_v
    size += 2 * ne * sizeof(int);  // ipiv_e ipiv_e1

    size += 1 * 64;  // align once to typical cache line size
    size += 1 * 8;  // align once to 8-byte boundary
    size = (size + 63) / 64 * 64;  // make multiple of typical cache line size

    return size;
}


void d_dense_qp_ipm_ws_create(struct d_dense_qp_dim* dim, struct d_dense_qp_ipm_arg* arg, struct d_dense_qp_ipm_ws* workspace, void* mem) {

    // zero memory (to avoid corrupted memory like e.g. NaN)
    hpipm_size_t memsize = d_dense_qp_ipm_ws_memsize(dim, arg);
    hpipm_zero_memset(memsize, mem);

    int nv = dim->nv;
    int ne = dim->ne;
    int nb = dim->nb;
    int ng = dim->ng;
    int ns = dim->ns;


    // core struct
    struct d_core_qp_ipm_workspace* sr_ptr = mem;

    // core workspace
    workspace->core_workspace = sr_ptr;
    sr_ptr += 1;
    struct d_core_qp_ipm_workspace* cws = workspace->core_workspace;


    // res struct
    struct d_dense_qp_res* res_ptr = (struct d_dense_qp_res*) sr_ptr;

    workspace->res = res_ptr;
    workspace->res->dim = dim;
    res_ptr += 1;
    workspace->res_itref = res_ptr;
    res_ptr += 1;
    workspace->res_step = res_ptr;
    res_ptr += 1;


    // res workspace struct
    struct d_dense_qp_res_ws* res_ws_ptr = (struct d_dense_qp_res_ws*) res_ptr;

    workspace->res_ws = res_ws_ptr;
    res_ws_ptr += 1;


    // qp sol struct
    struct d_dense_qp_sol* qp_sol_ptr = (struct d_dense_qp_sol*) res_ws_ptr;

    workspace->sol_step = qp_sol_ptr;
    qp_sol_ptr += 1;
    workspace->sol_itref = qp_sol_ptr;
    qp_sol_ptr += 1;


    // qp struct
    struct d_dense_qp* qp_ptr = (struct d_dense_qp*) qp_sol_ptr;

    workspace->qp_step = qp_ptr;
    qp_ptr += 1;
    workspace->qp_itref = qp_ptr;
    qp_ptr += 1;


    // matrix struct
    struct mat* sm_ptr = (struct mat*) qp_ptr;

    workspace->Lv = sm_ptr;
    sm_ptr += 2;
    workspace->AL = sm_ptr;
    sm_ptr += 1;
    workspace->Le = sm_ptr;
    sm_ptr += 1;
    workspace->Ctx = sm_ptr;
    sm_ptr += 1;
    if (arg->lq_fact > 0) {
        workspace->lq0 = sm_ptr;
        sm_ptr += 1;
        workspace->lq1 = sm_ptr;
        sm_ptr += 1;
    }
    if (arg->kkt_fact_alg == 0) {
        workspace->A_LQ = sm_ptr;
        sm_ptr += 1;
        workspace->A_Q = sm_ptr;
        sm_ptr += 1;
        workspace->Zt = sm_ptr;
        sm_ptr += 1;
        workspace->ZtH = sm_ptr;
        sm_ptr += 1;
        workspace->ZtHZ = sm_ptr;
        sm_ptr += 1;
    } else {
        // TODO
    }
    if (arg->remove_lin_dep_eq != 0) {
        workspace->A_li = sm_ptr;
        sm_ptr += 1;
        workspace->Ab_LU = sm_ptr;
        sm_ptr += 1;
    }


    // vector struct
    struct vec* sv_ptr = (struct vec*) sm_ptr;

    workspace->sol_step->v = sv_ptr;
    sv_ptr += 1;
    workspace->sol_step->pi = sv_ptr;
    sv_ptr += 1;
    workspace->sol_step->lam = sv_ptr;
    sv_ptr += 1;
    workspace->sol_step->t = sv_ptr;
    sv_ptr += 1;
    workspace->res->res_g = sv_ptr;
    sv_ptr += 1;
    workspace->res->res_b = sv_ptr;
    sv_ptr += 1;
    workspace->res->res_d = sv_ptr;
    sv_ptr += 1;
    workspace->res->res_m = sv_ptr;
    sv_ptr += 1;
    workspace->Gamma = sv_ptr;
    sv_ptr += 1;
    workspace->gamma = sv_ptr;
    sv_ptr += 1;
    workspace->Zs_inv = sv_ptr;
    sv_ptr += 1;
    workspace->lv = sv_ptr;
    sv_ptr += 1;
    workspace->sv = sv_ptr;
    sv_ptr += 1;
    workspace->se = sv_ptr;
    sv_ptr += 1;
    workspace->tmp_nbg = sv_ptr;
    sv_ptr += 4;
    workspace->res_ws->tmp_nbg = sv_ptr;
    sv_ptr += 2;
    workspace->tmp_ns = sv_ptr;
    sv_ptr += 1;
    workspace->res_ws->tmp_ns = sv_ptr;
    sv_ptr += 1;
    workspace->tmp_m = sv_ptr;
    sv_ptr += 1;
    workspace->tmp_nv = sv_ptr;
    sv_ptr += 1;
    workspace->tmp_2ns = sv_ptr;
    sv_ptr += 1;
    workspace->tmp_nv2ns = sv_ptr;
    sv_ptr += 2;
    if (arg->kkt_fact_alg == 0) {
        workspace->xy = sv_ptr;
        sv_ptr += 1;
        workspace->Yxy = sv_ptr;
        sv_ptr += 1;
        workspace->xz = sv_ptr;
        sv_ptr += 1;
    } else {
        // TODO
    }
    if (arg->remove_lin_dep_eq != 0) {
        workspace->b_li = sv_ptr;
        sv_ptr += 1;
    }


    // align to 8-byte boundary
    hpipm_size_t s_ptr = (hpipm_size_t) sv_ptr;
    s_ptr = (s_ptr + 7) / 8 * 8;


    // double/float stuff
    double* d_ptr = (double*) s_ptr;

    workspace->stat = d_ptr;
    int stat_m = 18;
    d_ptr += stat_m * (1 + arg->stat_max);


    // int suff
    int* i_ptr = (int*) d_ptr;

    workspace->ipiv_v = i_ptr;
    i_ptr += nv;

    workspace->ipiv_e = i_ptr;
    i_ptr += ne;

    workspace->ipiv_e1 = i_ptr;
    i_ptr += ne;


    // align to typical cache line size
    s_ptr = (hpipm_size_t) i_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;


    // void stuf
    char* c_ptr = (char*) s_ptr;

    d_dense_qp_sol_create(dim, workspace->sol_itref, c_ptr);
    c_ptr += workspace->sol_itref->memsize;

    d_dense_qp_res_create(dim, workspace->res_itref, c_ptr);
    c_ptr += workspace->res_itref->memsize;

    d_dense_qp_res_create(dim, workspace->res_step, c_ptr);
    c_ptr += workspace->res_step->memsize;

    create_mat(nv + 1, nv, workspace->Lv, c_ptr);
    c_ptr += workspace->Lv->memsize;

    create_mat(nv + 1, nv, workspace->Lv + 1, c_ptr);
    c_ptr += workspace->Lv[1].memsize;

    create_mat(ne, nv, workspace->AL, c_ptr);
    c_ptr += workspace->AL->memsize;

    create_mat(ne, ne, workspace->Le, c_ptr);
    c_ptr += workspace->Le->memsize;

    create_mat(nv + 1, ng, workspace->Ctx, c_ptr);
    c_ptr += workspace->Ctx->memsize;

    if (arg->lq_fact > 0) {
        create_mat(ne, ne + nv, workspace->lq0, c_ptr);
        c_ptr += workspace->lq0->memsize;

        create_mat(nv, nv + nv + ng, workspace->lq1, c_ptr);
        c_ptr += workspace->lq1->memsize;
    }

    if (arg->kkt_fact_alg == 0) {
        create_mat(ne, nv, workspace->A_LQ, c_ptr);
        c_ptr += workspace->A_LQ->memsize;

        create_mat(nv, nv, workspace->A_Q, c_ptr);
        c_ptr += workspace->A_Q->memsize;

        if (arg->remove_lin_dep_eq != 0) {
            create_mat(nv, nv, workspace->Zt, c_ptr);
            c_ptr += workspace->Zt->memsize;

            create_mat(nv, nv, workspace->ZtH, c_ptr);
            c_ptr += workspace->ZtH->memsize;

            create_mat(nv, nv, workspace->ZtHZ, c_ptr);
            c_ptr += workspace->ZtHZ->memsize;
        } else {
            create_mat(nv - ne, nv, workspace->Zt, c_ptr);
            c_ptr += workspace->Zt->memsize;

            create_mat(nv - ne, nv, workspace->ZtH, c_ptr);
            c_ptr += workspace->ZtH->memsize;

            create_mat(nv - ne, nv - ne, workspace->ZtHZ, c_ptr);
            c_ptr += workspace->ZtHZ->memsize;
        }
    } else {
        // TODO
    }
    if (arg->remove_lin_dep_eq != 0) {
        create_mat(ne, nv, workspace->A_li, c_ptr);
        c_ptr += workspace->A_li->memsize;

        create_mat(ne, nv + 1, workspace->Ab_LU, c_ptr);
        c_ptr += workspace->Ab_LU->memsize;
    }

    create_vec(nv, workspace->lv, c_ptr);
    c_ptr += workspace->lv->memsize;

    create_vec(nv, workspace->sv, c_ptr);
    c_ptr += workspace->sv->memsize;

    create_vec(ne, workspace->se, c_ptr);
    c_ptr += workspace->se->memsize;

    create_vec(2 * ns, workspace->Zs_inv, c_ptr);
    c_ptr += workspace->Zs_inv->memsize;

    create_vec(nb + ng, workspace->tmp_nbg + 0, c_ptr);
    create_vec(nb + ng, workspace->res_ws->tmp_nbg + 0, c_ptr);
    c_ptr += (workspace->tmp_nbg + 0)->memsize;

    create_vec(nb + ng, workspace->tmp_nbg + 1, c_ptr);
    create_vec(nb + ng, workspace->res_ws->tmp_nbg + 1, c_ptr);
    c_ptr += (workspace->tmp_nbg + 1)->memsize;

    create_vec(nb + ng, workspace->tmp_nbg + 2, c_ptr);
    c_ptr += (workspace->tmp_nbg + 2)->memsize;

    create_vec(nb + ng, workspace->tmp_nbg + 3, c_ptr);
    c_ptr += (workspace->tmp_nbg + 3)->memsize;

    create_vec(ns, workspace->tmp_ns + 0, c_ptr);
    create_vec(ns, workspace->res_ws->tmp_ns + 0, c_ptr);
    c_ptr += (workspace->tmp_ns + 0)->memsize;

    create_vec(2 * nb + 2 * ng + 2 * ns, workspace->tmp_m, c_ptr);
    c_ptr += (workspace->tmp_m)->memsize;

    create_vec(nv, workspace->tmp_nv, c_ptr);
    c_ptr += workspace->tmp_nv->memsize;

    create_vec(2 * ns, workspace->tmp_2ns, c_ptr);
    c_ptr += workspace->tmp_2ns->memsize;

    create_vec(nv + 2 * ns, workspace->tmp_nv2ns + 0, c_ptr);
    c_ptr += (workspace->tmp_nv2ns + 0)->memsize;

    create_vec(nv + 2 * ns, workspace->tmp_nv2ns + 1, c_ptr);
    c_ptr += (workspace->tmp_nv2ns + 1)->memsize;

    if (arg->kkt_fact_alg == 0) {
        create_vec(ne, workspace->xy, c_ptr);
        c_ptr += workspace->xy->memsize;

        create_vec(nv, workspace->Yxy, c_ptr);
        c_ptr += workspace->Yxy->memsize;

        if (arg->remove_lin_dep_eq != 0) {
            create_vec(nv, workspace->xz, c_ptr);
            c_ptr += workspace->xz->memsize;
        } else {
            create_vec(nv - ne, workspace->xz, c_ptr);
            c_ptr += workspace->xz->memsize;
        }

    } else {
        // TODO
    }
    if (arg->remove_lin_dep_eq != 0) {
        create_vec(ne, workspace->b_li, c_ptr);
        c_ptr += workspace->b_li->memsize;
    }

    d_create_core_qp_ipm(nv + 2 * ns, ne, 2 * nb + 2 * ng + 2 * ns, cws, c_ptr);
    c_ptr += workspace->core_workspace->memsize;

    workspace->lq_work_null = c_ptr;
    c_ptr += ne > 0 ? dgelqf_worksize(ne, nv) : 0;

    if (arg->lq_fact > 0) {
        workspace->lq_work0 = c_ptr;
        c_ptr += ne > 0 ? dgelqf_worksize(ne, nv) : 0;

        workspace->lq_work1 = c_ptr;
        c_ptr += dgelqf_worksize(nv, nv + nv + ng);
    }

    if (arg->kkt_fact_alg == 0) {
        workspace->orglq_work_null = c_ptr;
        c_ptr += dorglq_worksize(nv, nv, ne);
    } else {
        // TODO
    }

    // alias members of workspace and core_workspace
    //
    create_vec(nv + 2 * ns, workspace->sol_step->v, cws->dv);
    //
    create_vec(ne, workspace->sol_step->pi, cws->dpi);
    //
    create_vec(2 * nb + 2 * ng + 2 * ns, workspace->sol_step->lam, cws->dlam);
    //
    create_vec(2 * nb + 2 * ng + 2 * ns, workspace->sol_step->t, cws->dt);
    //
    create_vec(nv + 2 * ns, workspace->res->res_g, cws->res_g);
    //
    create_vec(ne, workspace->res->res_b, cws->res_b);
    //
    create_vec(2 * nb + 2 * ng + 2 * ns, workspace->res->res_d, cws->res_d);
    //
    create_vec(2 * nb + 2 * ng + 2 * ns, workspace->res->res_m, cws->res_m);
    //
    create_vec(2 * nb + 2 * ng + 2 * ns, workspace->Gamma, cws->Gamma);
    //
    create_vec(2 * nb + 2 * ng + 2 * ns, workspace->gamma, cws->gamma);

    //
    workspace->sol_step->dim = dim;

    workspace->stat_max = arg->stat_max;
    workspace->stat_m = stat_m;

    //
    workspace->use_hess_fact = 0;
    workspace->use_A_fact = 0;

    // cache stuff
    workspace->lq_fact = arg->lq_fact;

    //
    workspace->memsize = memsize;  // d_dense_qp_ipm_ws_memsize(dim, arg);


#if defined(RUNTIME_CHECKS)
    if (c_ptr > ((char*) mem) + workspace->memsize) {
        printf("\nCreate_dense_qp_ipm: outside memory bounds!\n\n");
        exit(1);
    }
#endif
}


void d_dense_qp_ipm_get(char* field, struct d_dense_qp_ipm_ws* ws, void* value) {
    if (hpipm_strcmp(field, "status")) {
        d_dense_qp_ipm_get_status(ws, value);
    } else if (hpipm_strcmp(field, "iter")) {
        d_dense_qp_ipm_get_iter(ws, value);
    } else if (hpipm_strcmp(field, "max_res_stat")) {
        d_dense_qp_ipm_get_max_res_stat(ws, value);
    } else if (hpipm_strcmp(field, "max_res_eq")) {
        d_dense_qp_ipm_get_max_res_eq(ws, value);
    } else if (hpipm_strcmp(field, "max_res_ineq")) {
        d_dense_qp_ipm_get_max_res_ineq(ws, value);
    } else if (hpipm_strcmp(field, "max_res_comp")) {
        d_dense_qp_ipm_get_max_res_comp(ws, value);
    } else if (hpipm_strcmp(field, "obj")) {
        d_dense_qp_ipm_get_obj(ws, value);
    } else if (hpipm_strcmp(field, "stat")) {
        d_dense_qp_ipm_get_stat(ws, value);
    } else if (hpipm_strcmp(field, "stat_m")) {
        d_dense_qp_ipm_get_stat_m(ws, value);
    } else {
        printf("error: d_dense_qp_ipm_get: wrong field %s\n", field);
        exit(1);
    }
}


void d_dense_qp_ipm_get_status(struct d_dense_qp_ipm_ws* ws, int* status) {
    *status = ws->status;
}


void d_dense_qp_ipm_get_iter(struct d_dense_qp_ipm_ws* ws, int* iter) {
    *iter = ws->iter;
}


void d_dense_qp_ipm_get_max_res_stat(struct d_dense_qp_ipm_ws* ws, double* res_stat) {
    *res_stat = ws->res->res_max[0];
}


void d_dense_qp_ipm_get_max_res_eq(struct d_dense_qp_ipm_ws* ws, double* res_eq) {
    *res_eq = ws->res->res_max[1];
}


void d_dense_qp_ipm_get_max_res_ineq(struct d_dense_qp_ipm_ws* ws, double* res_ineq) {
    *res_ineq = ws->res->res_max[2];
}


void d_dense_qp_ipm_get_max_res_comp(struct d_dense_qp_ipm_ws* ws, double* res_comp) {
    *res_comp = ws->res->res_max[3];
}


void d_dense_qp_ipm_get_obj(struct d_dense_qp_ipm_ws* ws, double* obj) {
    *obj = ws->res->obj;
}


void d_dense_qp_ipm_get_stat(struct d_dense_qp_ipm_ws* ws, double** stat) {
    *stat = ws->stat;
}


void d_dense_qp_ipm_get_stat_m(struct d_dense_qp_ipm_ws* ws, int* stat_m) {
    *stat_m = ws->stat_m;
}


void d_dense_qp_init_var(struct d_dense_qp* qp, struct d_dense_qp_sol* qp_sol, struct d_dense_qp_ipm_arg* arg, struct d_dense_qp_ipm_ws* ws) {

    //	struct d_core_qp_ipm_workspace *cws = ws->core_workspace;

    // extract cws members
    int nv = qp->dim->nv;
    int ne = qp->dim->ne;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;
    int ns = qp->dim->ns;

    double* d = qp->d->pa;
    int* idxb = qp->idxb;

    double* v = qp_sol->v->pa;
    double* pi = qp_sol->pi->pa;
    double* lam = qp_sol->lam->pa;
    double* t = qp_sol->t->pa;

    double mu0 = arg->mu0;

    // local variables
    int ii;
    int idxb0;
    double thr0 = 0.5;


    // primal and dual variables
    if (arg->warm_start == 2) {

        thr0 = 1e-1;

        for (ii = 0; ii < 2 * nb + 2 * ng + 2 * ns; ii++) {
            if (lam[ii] < thr0)
                lam[ii] = thr0;
            if (t[ii] < thr0)
                t[ii] = thr0;
        }
    }


    // primal variables
    if (arg->warm_start == 0) {
        // cold start
        for (ii = 0; ii < nv + 2 * ns; ii++) {
            v[ii] = 0.0;
        }
    }

    // equality constraints
    for (ii = 0; ii < ne; ii++) {
        pi[ii] = 0.0;
    }

    // box constraints
    for (ii = 0; ii < nb; ii++) {
#if 1
        idxb0 = idxb[ii];
        t[0 + ii] = -d[0 + ii] + v[idxb0];
        t[nb + ng + ii] = -d[nb + ng + ii] - v[idxb0];
        if (t[0 + ii] < thr0) {
            if (t[nb + ng + ii] < thr0) {
                v[idxb0] = 0.5 * (d[0 + ii] + d[nb + ng + ii]);
                t[0 + ii] = thr0;
                t[nb + ng + ii] = thr0;
            } else {
                t[0 + ii] = thr0;
                v[idxb0] = d[0 + ii] + thr0;
            }
        } else if (t[nb + ng + ii] < thr0) {
            t[nb + ng + ii] = thr0;
            v[idxb0] = -d[nb + ng + ii] - thr0;
        }
#else
        t[0 + ii] = 1.0;
        t[nb + ng + ii] = 1.0;
#endif
        lam[0 + ii] = mu0 / t[0 + ii];
        lam[nb + ng + ii] = mu0 / t[nb + ng + ii];
    }

    // general constraints
    if (ng > 0) {
        dgemv_t(nv, ng, 1.0, qp->Ct, 0, 0, qp_sol->v, 0, 0.0, qp_sol->t, nb, qp_sol->t, nb);
        for (ii = 0; ii < ng; ii++) {
#if 1
            t[2 * nb + ng + ii] = t[nb + ii];
            t[nb + ii] -= d[nb + ii];
            t[2 * nb + ng + ii] -= d[2 * nb + ng + ii];
            //		t[nb+ii]      = fmax( thr0, t[nb+ii] );
            //		t[2*nb+ng+ii] = fmax( thr0, t[2*nb+ng+ii] );
            t[nb + ii] = thr0 > t[nb + ii] ? thr0 : t[nb + ii];
            t[2 * nb + ng + ii] = thr0 > t[2 * nb + ng + ii] ? thr0 : t[2 * nb + ng + ii];
#else
            t[nb + ii] = 1.0;
            t[2 * nb + ng + ii] = 1.0;
#endif
            lam[nb + ii] = mu0 / t[nb + ii];
            lam[2 * nb + ng + ii] = mu0 / t[2 * nb + ng + ii];
        }
    }

    // soft constraints
    for (ii = 0; ii < ns; ii++) {
        t[2 * nb + 2 * ng + ii] = 1.0;  // thr0;
        t[2 * nb + 2 * ng + ns + ii] = 1.0;  // thr0;
        lam[2 * nb + 2 * ng + ii] = mu0 / t[2 * nb + 2 * ng + ii];
        lam[2 * nb + 2 * ng + ns + ii] = mu0 / t[2 * nb + 2 * ng + ns + ii];
    }
}


void d_dense_qp_ipm_abs_step(int kk, struct d_dense_qp* qp, struct d_dense_qp_sol* qp_sol, struct d_dense_qp_ipm_arg* arg, struct d_dense_qp_ipm_ws* ws) {

    struct d_core_qp_ipm_workspace* cws = ws->core_workspace;

    int nv = qp->dim->nv;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;
    int ns = qp->dim->ns;

    double tmp;
    double mu_aff0;  //, mu;

    dvecsc(cws->nc, -1.0, ws->tmp_m, 0);

    d_backup_res_m(cws);

    // tau_min as barrier parameter for affine step
    d_compute_tau_min_qp(cws);

    // fact solve
    d_fact_solve_kkt_step_dense_qp(ws->qp_step, ws->sol_step, arg, ws);
    // compute step
    daxpy(cws->nv, -1.0, qp_sol->v, 0, ws->sol_step->v, 0, ws->sol_step->v, 0);
    daxpy(cws->ne, -1.0, qp_sol->pi, 0, ws->sol_step->pi, 0, ws->sol_step->pi, 0);
    daxpy(cws->nc, -1.0, qp_sol->lam, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
    daxpy(cws->nc, -1.0, qp_sol->t, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);
    if (ws->mask_constr) {
        // mask out disregarded constraints
        dvecmul(2 * ns, qp->d_mask, 2 * nb + 2 * ng, ws->sol_step->v, nv, ws->sol_step->v, nv);
        dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);
        dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
    }

    // alpha
    d_compute_alpha_qp(cws);
    if (kk + 1 < ws->stat_max)
        ws->stat[ws->stat_m * (kk + 1) + 0] = cws->alpha;

    // Mehrotra's predictor-corrector
    if (arg->pred_corr == 1) {
        // mu_aff
        d_compute_mu_aff_qp(cws);
        if (kk + 1 < ws->stat_max)
            ws->stat[ws->stat_m * (kk + 1) + 1] = cws->mu_aff;

        tmp = cws->mu_aff / cws->mu;
        cws->sigma = tmp * tmp * tmp;
        if (kk + 1 < ws->stat_max)
            ws->stat[ws->stat_m * (kk + 1) + 2] = cws->sigma;

        d_compute_centering_correction_qp(cws);
        if (ws->mask_constr) {
            // mask out disregarded constraints
            dvecmul(cws->nc, qp->d_mask, 0, ws->res->res_m, 0, ws->res->res_m, 0);
        }

        // fact and solve kkt
        d_solve_kkt_step_dense_qp(ws->qp_step, ws->sol_step, arg, ws);
        // compute step
        daxpy(cws->nv, -1.0, qp_sol->v, 0, ws->sol_step->v, 0, ws->sol_step->v, 0);
        daxpy(cws->ne, -1.0, qp_sol->pi, 0, ws->sol_step->pi, 0, ws->sol_step->pi, 0);
        daxpy(cws->nc, -1.0, qp_sol->lam, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
        daxpy(cws->nc, -1.0, qp_sol->t, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);
        if (ws->mask_constr) {
            // mask out disregarded constraints
            dvecmul(2 * ns, qp->d_mask, 2 * nb + 2 * ng, ws->sol_step->v, nv, ws->sol_step->v, nv);
            dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);
            dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
        }

        // alpha
        d_compute_alpha_qp(cws);
        if (kk + 1 < ws->stat_max) {
            ws->stat[ws->stat_m * (kk + 1) + 3] = cws->alpha_prim;
            ws->stat[ws->stat_m * (kk + 1) + 4] = cws->alpha_dual;
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
                d_solve_kkt_step_dense_qp(ws->qp_step, ws->sol_step, arg, ws);
                // compute step
                daxpy(cws->nv, -1.0, qp_sol->v, 0, ws->sol_step->v, 0, ws->sol_step->v, 0);
                daxpy(cws->ne, -1.0, qp_sol->pi, 0, ws->sol_step->pi, 0, ws->sol_step->pi, 0);
                daxpy(cws->nc, -1.0, qp_sol->lam, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
                daxpy(cws->nc, -1.0, qp_sol->t, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);
                if (ws->mask_constr) {
                    // mask out disregarded constraints
                    dvecmul(2 * ns, qp->d_mask, 2 * nb + 2 * ng, ws->sol_step->v, nv, ws->sol_step->v, nv);
                    dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);
                    dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
                }

                // alpha
                d_compute_alpha_qp(cws);
                if (kk + 1 < ws->stat_max) {
                    ws->stat[ws->stat_m * (kk + 1) + 3] = cws->alpha_prim;
                    ws->stat[ws->stat_m * (kk + 1) + 4] = cws->alpha_dual;
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


void d_dense_qp_ipm_delta_step(int kk, struct d_dense_qp* qp, struct d_dense_qp_sol* qp_sol, struct d_dense_qp_ipm_arg* arg, struct d_dense_qp_ipm_ws* ws) {

    // dim
    int nv = qp->dim->nv;
    int ne = qp->dim->ne;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;
    int ns = qp->dim->ns;

    struct d_core_qp_ipm_workspace* cws = ws->core_workspace;

    int itref0 = 0, itref1 = 0, iter_ref_step;
    double tmp;
    double mu_aff0;  //, mu;

    double itref_qp_norm[4] = {0, 0, 0, 0};
    double itref_qp_norm0[4] = {0, 0, 0, 0};
    //	int ndp0, ndp1;

    double* qp_res_max = ws->res->res_max;

    int force_lq = 0;

    ws->scale = arg->scale;

    // step body

    d_backup_res_m(cws);

    // tau_min as barrier parameter for affine step
    d_compute_tau_min_qp(cws);

    // fact and solve kkt
    if (ws->lq_fact == 0) {
        // syrk+cholesky
        d_fact_solve_kkt_step_dense_qp(ws->qp_step, ws->sol_step, arg, ws);
        if (ws->mask_constr) {
            // mask out disregarded constraints
            dvecmul(2 * ns, qp->d_mask, 2 * nb + 2 * ng, ws->sol_step->v, nv, ws->sol_step->v, nv);
            dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);
            dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
        }
        if (kk + 1 < ws->stat_max)
            ws->stat[ws->stat_m * (kk + 1) + 11] = 0;
    } else if (ws->lq_fact == 1 & force_lq == 0) {
        // syrk+chol, switch to lq when needed
        d_fact_solve_kkt_step_dense_qp(ws->qp_step, ws->sol_step, arg, ws);
        if (ws->mask_constr) {
            // mask out disregarded constraints
            dvecmul(2 * ns, qp->d_mask, 2 * nb + 2 * ng, ws->sol_step->v, nv, ws->sol_step->v, nv);
            dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);
            dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
        }

        // compute res of linear system
        d_dense_qp_res_compute_lin(ws->qp_step, qp_sol, ws->sol_step, ws->res_itref, ws->res_ws);
        // TODO mask
        if (ws->mask_constr) {
            // mask out disregarded constraints
            dvecmul(2 * ns, qp->d_mask, 2 * nb + 2 * ng, ws->res_itref->res_g, nv, ws->res_itref->res_g, nv);
            dvecmul(cws->nc, qp->d_mask, 0, ws->res_itref->res_d, 0, ws->res_itref->res_d, 0);
            dvecmul(cws->nc, qp->d_mask, 0, ws->res_itref->res_m, 0, ws->res_itref->res_m, 0);
        }
        d_dense_qp_res_compute_inf_norm(ws->res_itref);
        itref_qp_norm[0] = ws->res_itref->res_max[0];
        itref_qp_norm[1] = ws->res_itref->res_max[1];
        itref_qp_norm[2] = ws->res_itref->res_max[2];
        itref_qp_norm[3] = ws->res_itref->res_max[3];

        if (kk + 1 < ws->stat_max)
            ws->stat[ws->stat_m * (kk + 1) + 11] = 0;

        // printf("\n%e\t%e\t%e\t%e\n", itref_qp_norm[0], itref_qp_norm[1], itref_qp_norm[2], itref_qp_norm[3]);

        // inaccurate factorization: switch to lq
        if (

                ((itref_qp_norm[0] == 0.0) & isnan(VECEL(ws->res_itref->res_g, 0))) |
                // #else
                //                 (itref_qp_norm[0] == 0.0 & VECEL(ws->res_itref->res_g, 0) != VECEL(ws->res_itref->res_g, 0)) |

                (itref_qp_norm[0] > 1e-5) |
                (itref_qp_norm[1] > 1e-5) |
                (itref_qp_norm[2] > 1e-5) |
                (itref_qp_norm[3] > 1e-5)) {

            // refactorize using lq
            d_fact_lq_solve_kkt_step_dense_qp(ws->qp_step, ws->sol_step, arg, ws);
            if (ws->mask_constr) {
                // mask out disregarded constraints
                dvecmul(2 * ns, qp->d_mask, 2 * nb + 2 * ng, ws->sol_step->v, nv, ws->sol_step->v, nv);
                dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);
                dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
            }

            // switch to lq
            force_lq = 1;

            if (kk + 1 < ws->stat_max)
                ws->stat[ws->stat_m * (kk + 1) + 11] = 1;
        }
    } else  // ws->lq_fact==2
    {
        // lq
        d_fact_lq_solve_kkt_step_dense_qp(ws->qp_step, ws->sol_step, arg, ws);
        if (ws->mask_constr) {
            // mask out disregarded constraints
            dvecmul(2 * ns, qp->d_mask, 2 * nb + 2 * ng, ws->sol_step->v, nv, ws->sol_step->v, nv);
            dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);
            dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
        }
        if (kk + 1 < ws->stat_max)
            ws->stat[ws->stat_m * (kk + 1) + 10] = 1;
    }

    // iterative refinement on prediction step
    if (arg->itref_pred_max == 0) {
        if (kk + 1 < ws->stat_max) {
            ws->stat[ws->stat_m * (kk + 1) + 14] = 0.0;
            ws->stat[ws->stat_m * (kk + 1) + 15] = 0.0;
            ws->stat[ws->stat_m * (kk + 1) + 16] = 0.0;
            ws->stat[ws->stat_m * (kk + 1) + 17] = 0.0;
        }
    } else {
        for (itref0 = 0; itref0 < arg->itref_pred_max; itref0++) {

            d_dense_qp_res_compute_lin(ws->qp_step, qp_sol, ws->sol_step, ws->res_itref, ws->res_ws);
            if (ws->mask_constr) {
                // mask out disregarded constraints
                dvecmul(2 * ns, qp->d_mask, 2 * nb + 2 * ng, ws->res_itref->res_g, nv, ws->res_itref->res_g, nv);
                dvecmul(cws->nc, qp->d_mask, 0, ws->res_itref->res_d, 0, ws->res_itref->res_d, 0);
                dvecmul(cws->nc, qp->d_mask, 0, ws->res_itref->res_m, 0, ws->res_itref->res_m, 0);
            }
            d_dense_qp_res_compute_inf_norm(ws->res_itref);
            itref_qp_norm[0] = ws->res_itref->res_max[0];
            itref_qp_norm[1] = ws->res_itref->res_max[1];
            itref_qp_norm[2] = ws->res_itref->res_max[2];
            itref_qp_norm[3] = ws->res_itref->res_max[3];
            if (kk + 1 < ws->stat_max) {
                ws->stat[ws->stat_m * (kk + 1) + 14] = itref_qp_norm[0];
                ws->stat[ws->stat_m * (kk + 1) + 15] = itref_qp_norm[1];
                ws->stat[ws->stat_m * (kk + 1) + 16] = itref_qp_norm[2];
                ws->stat[ws->stat_m * (kk + 1) + 17] = itref_qp_norm[3];
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

            d_solve_kkt_step_dense_qp(ws->qp_itref, ws->sol_itref, arg, ws);
            if (ws->mask_constr) {
                // mask out disregarded constraints
                dvecmul(2 * ns, qp->d_mask, 2 * nb + 2 * ng, ws->sol_step->v, nv, ws->sol_step->v, nv);
                dvecmul(cws->nc, qp->d_mask, 0, ws->sol_itref->t, 0, ws->sol_itref->t, 0);
                dvecmul(cws->nc, qp->d_mask, 0, ws->sol_itref->lam, 0, ws->sol_itref->lam, 0);
            }

            daxpy(nv + 2 * ns, 1.0, ws->sol_itref->v, 0, ws->sol_step->v, 0, ws->sol_step->v, 0);
            daxpy(ne, 1.0, ws->sol_itref->pi, 0, ws->sol_step->pi, 0, ws->sol_step->pi, 0);
            daxpy(2 * nb + 2 * ng + 2 * ns, 1.0, ws->sol_itref->lam, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
            daxpy(2 * nb + 2 * ng + 2 * ns, 1.0, ws->sol_itref->t, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);
        }
        if (itref0 == arg->itref_pred_max) {
            d_dense_qp_res_compute_lin(ws->qp_step, qp_sol, ws->sol_step, ws->res_itref, ws->res_ws);
            if (ws->mask_constr) {
                // mask out disregarded constraints
                dvecmul(2 * ns, qp->d_mask, 2 * nb + 2 * ng, ws->res_itref->res_g, nv, ws->res_itref->res_g, nv);
                dvecmul(cws->nc, qp->d_mask, 0, ws->res_itref->res_d, 0, ws->res_itref->res_d, 0);
                dvecmul(cws->nc, qp->d_mask, 0, ws->res_itref->res_m, 0, ws->res_itref->res_m, 0);
            }
            d_dense_qp_res_compute_inf_norm(ws->res_itref);
            itref_qp_norm[0] = ws->res_itref->res_max[0];
            itref_qp_norm[1] = ws->res_itref->res_max[1];
            itref_qp_norm[2] = ws->res_itref->res_max[2];
            itref_qp_norm[3] = ws->res_itref->res_max[3];
            if (kk + 1 < ws->stat_max) {
                ws->stat[ws->stat_m * (kk + 1) + 14] = itref_qp_norm[0];
                ws->stat[ws->stat_m * (kk + 1) + 15] = itref_qp_norm[1];
                ws->stat[ws->stat_m * (kk + 1) + 16] = itref_qp_norm[2];
                ws->stat[ws->stat_m * (kk + 1) + 17] = itref_qp_norm[3];
            }
        }
    }

    if (kk + 1 < ws->stat_max)
        ws->stat[ws->stat_m * (kk + 1) + 12] = itref0;

#if 0
	ndp0 = 0;
	for(ii=0; ii<qp->dim->nv; ii++)
		{
		if(ws->Lv->dA[ii]<=0)
			{
			ndp0 = ii;
			break;
			}
		}
	ndp1 = 0;
	for(ii=0; ii<qp->dim->ne; ii++)
		{
		if(ws->Le->dA[ii]<=0)
			{
			ndp1 = ii;
			break;
			}
		}
#endif

    // alpha
    d_compute_alpha_qp(cws);
    if (kk + 1 < ws->stat_max)
        ws->stat[ws->stat_m * (kk + 1) + 0] = cws->alpha;

    // Mehrotra's predictor-corrector
    if (arg->pred_corr == 1) {
        // mu_aff
        d_compute_mu_aff_qp(cws);
        if (kk + 1 < ws->stat_max)
            ws->stat[ws->stat_m * (kk + 1) + 1] = cws->mu_aff;

        // compute centering parameter
        tmp = cws->mu_aff / cws->mu;
        cws->sigma = tmp * tmp * tmp;
        if (kk + 1 < ws->stat_max)
            ws->stat[ws->stat_m * (kk + 1) + 2] = cws->sigma;

        d_compute_centering_correction_qp(cws);
        if (ws->mask_constr) {
            // mask out disregarded constraints
            dvecmul(cws->nc, qp->d_mask, 0, ws->res->res_m, 0, ws->res->res_m, 0);
        }

        // solve kkt
        d_solve_kkt_step_dense_qp(ws->qp_step, ws->sol_step, arg, ws);
        if (ws->mask_constr) {
            // mask out disregarded constraints
            dvecmul(2 * ns, qp->d_mask, 2 * nb + 2 * ng, ws->sol_step->v, nv, ws->sol_step->v, nv);
            dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);
            dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
        }

        // alpha
        d_compute_alpha_qp(cws);
        if (kk + 1 < ws->stat_max) {
            ws->stat[ws->stat_m * (kk + 1) + 3] = cws->alpha_prim;
            ws->stat[ws->stat_m * (kk + 1) + 4] = cws->alpha_dual;
        }

        // conditional Mehrotra's predictor-corrector
        if (arg->cond_pred_corr == 1) {

            // save mu_aff (from prediction sol_step)
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
                d_solve_kkt_step_dense_qp(ws->qp_step, ws->sol_step, arg, ws);
                if (ws->mask_constr) {
                    // mask out disregarded constraints
                    dvecmul(2 * ns, qp->d_mask, 2 * nb + 2 * ng, ws->sol_step->v, nv, ws->sol_step->v, nv);
                    dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);
                    dvecmul(cws->nc, qp->d_mask, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
                }

                // alpha
                d_compute_alpha_qp(cws);
                if (kk + 1 < ws->stat_max) {
                    ws->stat[ws->stat_m * (kk + 1) + 3] = cws->alpha_prim;
                    ws->stat[ws->stat_m * (kk + 1) + 4] = cws->alpha_dual;
                }
            }
        }

        iter_ref_step = 0;
        if (arg->itref_corr_max > 0) {
            for (itref1 = 0; itref1 < arg->itref_corr_max; itref1++) {

                d_dense_qp_res_compute_lin(ws->qp_step, qp_sol, ws->sol_step, ws->res_itref, ws->res_ws);
                if (ws->mask_constr) {
                    // mask out disregarded constraints
                    dvecmul(2 * ns, qp->d_mask, 2 * nb + 2 * ng, ws->res_itref->res_g, nv, ws->res_itref->res_g, nv);
                    dvecmul(cws->nc, qp->d_mask, 0, ws->res_itref->res_d, 0, ws->res_itref->res_d, 0);
                    dvecmul(cws->nc, qp->d_mask, 0, ws->res_itref->res_m, 0, ws->res_itref->res_m, 0);
                }
                d_dense_qp_res_compute_inf_norm(ws->res_itref);
                itref_qp_norm[0] = ws->res_itref->res_max[0];
                itref_qp_norm[1] = ws->res_itref->res_max[1];
                itref_qp_norm[2] = ws->res_itref->res_max[2];
                itref_qp_norm[3] = ws->res_itref->res_max[3];
                if (kk + 1 < ws->stat_max) {
                    ws->stat[ws->stat_m * (kk + 1) + 14] = itref_qp_norm[0];
                    ws->stat[ws->stat_m * (kk + 1) + 15] = itref_qp_norm[1];
                    ws->stat[ws->stat_m * (kk + 1) + 16] = itref_qp_norm[2];
                    ws->stat[ws->stat_m * (kk + 1) + 17] = itref_qp_norm[3];
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

                d_solve_kkt_step_dense_qp(ws->qp_itref, ws->sol_itref, arg, ws);
                if (ws->mask_constr) {
                    // mask out disregarded constraints
                    dvecmul(2 * ns, qp->d_mask, 2 * nb + 2 * ng, ws->sol_step->v, nv, ws->sol_step->v, nv);
                    dvecmul(cws->nc, qp->d_mask, 0, ws->sol_itref->t, 0, ws->sol_itref->t, 0);
                    dvecmul(cws->nc, qp->d_mask, 0, ws->sol_itref->lam, 0, ws->sol_itref->lam, 0);
                }
                iter_ref_step = 1;

                daxpy(nv + 2 * ns, 1.0, ws->sol_itref->v, 0, ws->sol_step->v, 0, ws->sol_step->v, 0);
                daxpy(ne, 1.0, ws->sol_itref->pi, 0, ws->sol_step->pi, 0, ws->sol_step->pi, 0);
                daxpy(2 * nb + 2 * ng + 2 * ns, 1.0, ws->sol_itref->lam, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
                daxpy(2 * nb + 2 * ng + 2 * ns, 1.0, ws->sol_itref->t, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);
            }
            if (itref1 == arg->itref_corr_max) {
                d_dense_qp_res_compute_lin(ws->qp_step, qp_sol, ws->sol_step, ws->res_itref, ws->res_ws);
                if (ws->mask_constr) {
                    // mask out disregarded constraints
                    dvecmul(2 * ns, qp->d_mask, 2 * nb + 2 * ng, ws->res_itref->res_g, nv, ws->res_itref->res_g, nv);
                    dvecmul(cws->nc, qp->d_mask, 0, ws->res_itref->res_d, 0, ws->res_itref->res_d, 0);
                    dvecmul(cws->nc, qp->d_mask, 0, ws->res_itref->res_m, 0, ws->res_itref->res_m, 0);
                }
                d_dense_qp_res_compute_inf_norm(ws->res_itref);
                itref_qp_norm[0] = ws->res_itref->res_max[0];
                itref_qp_norm[1] = ws->res_itref->res_max[1];
                itref_qp_norm[2] = ws->res_itref->res_max[2];
                itref_qp_norm[3] = ws->res_itref->res_max[3];
                if (kk + 1 < ws->stat_max) {
                    ws->stat[ws->stat_m * (kk + 1) + 14] = itref_qp_norm[0];
                    ws->stat[ws->stat_m * (kk + 1) + 15] = itref_qp_norm[1];
                    ws->stat[ws->stat_m * (kk + 1) + 16] = itref_qp_norm[2];
                    ws->stat[ws->stat_m * (kk + 1) + 17] = itref_qp_norm[3];
                }
            }
        }

        if (iter_ref_step) {
            // alpha
            d_compute_alpha_qp(cws);
            if (kk + 1 < ws->stat_max) {
                ws->stat[ws->stat_m * (kk + 1) + 3] = cws->alpha_prim;
                ws->stat[ws->stat_m * (kk + 1) + 4] = cws->alpha_dual;
            }
        }
    }
    if (arg->itref_corr_max == 0) {
        if (kk + 1 < ws->stat_max) {
            ws->stat[ws->stat_m * (kk + 1) + 14] = 0.0;
            ws->stat[ws->stat_m * (kk + 1) + 15] = 0.0;
            ws->stat[ws->stat_m * (kk + 1) + 16] = 0.0;
            ws->stat[ws->stat_m * (kk + 1) + 17] = 0.0;
        }
    }
    if (kk + 1 < ws->stat_max)
        ws->stat[ws->stat_m * (kk + 1) + 13] = itref1;

    // TODO check for step length computation
    if (1) {
        // compute step minimizing phi
        d_dense_qp_compute_step_length(qp, qp_sol, arg, ws);

        // TODO put in new stat col ???
        if (kk + 1 < ws->stat_max) {
            ws->stat[ws->stat_m * (kk + 1) + 3] = cws->alpha_prim;
            ws->stat[ws->stat_m * (kk + 1) + 4] = cws->alpha_dual;
        }
    }

    //
    d_update_var_qp(cws);
    if (ws->mask_constr) {
        // mask out disregarded constraints
        dvecmul(cws->nc, qp->d_mask, 0, qp_sol->lam, 0, qp_sol->lam, 0);
    }
}


void d_dense_qp_ipm_solve(struct d_dense_qp* qp, struct d_dense_qp_sol* qp_sol, struct d_dense_qp_ipm_arg* arg, struct d_dense_qp_ipm_ws* ws) {
    // dim
    int nv = qp->dim->nv;
    int ne = qp->dim->ne;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;
    int ns = qp->dim->ns;

    struct d_core_qp_ipm_workspace* cws = ws->core_workspace;

    if (arg->remove_lin_dep_eq) {
        d_dense_qp_remove_lin_dep_eq(qp, arg, ws);
        if (ws->status == INCONS_EQ) {
            ws->iter = 0;
            goto call_return;
        }
    }

    int kk, ii;
    double mu;

    double* stat = ws->stat;
    int stat_m = ws->stat_m;
    double tau_min = arg->tau_min;

    // arg to core workspace
    cws->lam_min = arg->lam_min;
    cws->t_min = arg->t_min;
    cws->t_min_inv = arg->t_min > 0.0 ? 1.0 / arg->t_min : 1e30;
    cws->tau_min = arg->tau_min;
    cws->split_step = arg->split_step;
    cws->t_lam_min = arg->t_lam_min;

    // alias qp vectors into qp_sol
    cws->v = qp_sol->v->pa;
    cws->pi = qp_sol->pi->pa;
    cws->lam = qp_sol->lam->pa;
    cws->t = qp_sol->t->pa;

    // alias members of qp_step
    ws->qp_step->dim = qp->dim;
    ws->qp_step->Hv = qp->Hv;
    ws->qp_step->A = qp->A;
    ws->qp_step->Ct = qp->Ct;
    ws->qp_step->Z = qp->Z;
    ws->qp_step->idxb = qp->idxb;
    ws->qp_step->idxs_rev = qp->idxs_rev;
    ws->qp_step->gz = ws->res->res_g;
    ws->qp_step->b = ws->res->res_b;
    ws->qp_step->d = ws->res->res_d;
    ws->qp_step->m = ws->res->res_m;
    ws->qp_step->d_mask = qp->d_mask;  // XXX

    // alias members of qp_itref
    ws->qp_itref->dim = qp->dim;
    ws->qp_itref->Hv = qp->Hv;
    ws->qp_itref->A = qp->A;
    ws->qp_itref->Ct = qp->Ct;
    ws->qp_itref->Z = qp->Z;
    ws->qp_itref->idxb = qp->idxb;
    ws->qp_itref->idxs_rev = qp->idxs_rev;
    ws->qp_itref->gz = ws->res_itref->res_g;
    ws->qp_itref->b = ws->res_itref->res_b;
    ws->qp_itref->d = ws->res_itref->res_d;
    ws->qp_itref->m = ws->res_itref->res_m;
    ws->qp_itref->d_mask = qp->d_mask;  // XXX

    double* qp_res_max = ws->res->res_max;
    qp_res_max[0] = 0;
    qp_res_max[1] = 0;
    qp_res_max[2] = 0;
    qp_res_max[3] = 0;

    ws->use_hess_fact = 0;
    ws->use_A_fact = 0;

    //	for(ii=0; ii<cws->nc; ii++)
    //		{
    //		qp->m[0].pa[ii] = 1e-2;
    //		}
    // d_dense_qp_print(qp->dim, qp);
    // exit(1);

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
        d_fact_solve_kkt_unconstr_dense_qp(qp, qp_sol, arg, ws);
        if (arg->comp_res_exit) {
            // compute residuals
            d_dense_qp_res_compute(qp, qp_sol, ws->res, ws->res_ws);
            // XXX no constraints, so no mask
            d_dense_qp_res_compute_inf_norm(ws->res);
            // save infinity norm of residuals
            if (0 < ws->stat_max) {
                stat[6] = qp_res_max[0];
                stat[7] = qp_res_max[1];
                stat[8] = qp_res_max[2];
                stat[9] = qp_res_max[3];
                stat[10] = ws->res->obj;
            }
            cws->mu = ws->res->res_mu;
        }
        // save info before return
        ws->iter = 0;

        if (isnan(VECEL(qp_sol->v, 0))) {
            // NaN in the solution
            ws->status = NAN_SOL;
        }
        // #else
        //         if (VECEL(qp_sol->v, 0) != VECEL(qp_sol->v, 0)) {
        //             // NaN in the solution
        //             ws->status = NAN_sol;
        //         }

        else {
            // normal return
            ws->status = SUCCESS;
        }
        goto call_return;
    }


    // init solver
    d_dense_qp_init_var(qp, qp_sol, arg, ws);
    if (ws->mask_constr) {
        // mask out disregarded constraints
        dvecmul(cws->nc, qp->d_mask, 0, qp_sol->lam, 0, qp_sol->lam, 0);
    }

    cws->alpha = 1.0;


    // absolute IPM formulation
    if (arg->abs_form) {

        // alias members of qp_step
        ws->qp_step->dim = qp->dim;
        ws->qp_step->Hv = qp->Hv;
        ws->qp_step->A = qp->A;
        ws->qp_step->Ct = qp->Ct;
        ws->qp_step->Z = qp->Z;
        ws->qp_step->idxb = qp->idxb;
        ws->qp_step->idxs_rev = qp->idxs_rev;
        ws->qp_step->gz = qp->gz;
        ws->qp_step->b = qp->b;
        ws->qp_step->d = qp->d;
        ws->qp_step->m = ws->tmp_m;
        ws->qp_step->d_mask = qp->d_mask;  // XXX

        // alias core workspace
        cws->res_m = ws->qp_step->m->pa;
        //		cws->res_m_bkp = ws->qp_step->m->pa;

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
            d_dense_qp_ipm_abs_step(kk, qp, qp_sol, arg, ws);

            // compute mu
            mu = dvecmuldot(cws->nc, qp_sol->lam, 0, qp_sol->t, 0, ws->tmp_m, 0);
            mu /= cws->nc;
            cws->mu = mu;
            if (kk + 1 < ws->stat_max)
                stat[stat_m * (kk + 1) + 5] = mu;
        }

        if (arg->comp_res_exit) {
            // compute residuals
            d_dense_qp_res_compute(qp, qp_sol, ws->res, ws->res_ws);
            if (ws->mask_constr) {
                // mask out disregarded constraints
                dvecmul(2 * ns, qp->d_mask, 2 * nb + 2 * ng, ws->res->res_g, nv, ws->res->res_g, nv);
                dvecmul(cws->nc, qp->d_mask, 0, ws->res->res_d, 0, ws->res->res_d, 0);
                dvecmul(cws->nc, qp->d_mask, 0, ws->res->res_m, 0, ws->res->res_m, 0);
            }
            d_dense_qp_res_compute_inf_norm(ws->res);
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
    d_dense_qp_res_compute(qp, qp_sol, ws->res, ws->res_ws);
    if (ws->mask_constr) {
        // mask out disregarded constraints
        dvecmul(2 * ns, qp->d_mask, 2 * nb + 2 * ng, ws->res->res_g, nv, ws->res->res_g, nv);
        dvecmul(cws->nc, qp->d_mask, 0, ws->res->res_d, 0, ws->res->res_d, 0);
        dvecmul(cws->nc, qp->d_mask, 0, ws->res->res_m, 0, ws->res->res_m, 0);
    }
    d_dense_qp_res_compute_inf_norm(ws->res);
    cws->mu = ws->res->res_mu;
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
    for (kk = 0;
         kk<arg->iter_max &
            cws->alpha>
                 arg->alpha_min &
         (qp_res_max[0] > arg->res_g_max |
          qp_res_max[1] > arg->res_b_max |
          qp_res_max[2] > arg->res_d_max |
          fabs(qp_res_max[3] - tau_min) > arg->res_m_max);
         kk++) {

        // compute delta step
        d_dense_qp_ipm_delta_step(kk, qp, qp_sol, arg, ws);

        // compute residuals
        d_dense_qp_res_compute(qp, qp_sol, ws->res, ws->res_ws);
        if (ws->mask_constr) {
            // mask out disregarded constraints
            dvecmul(2 * ns, qp->d_mask, 2 * nb + 2 * ng, ws->res->res_g, nv, ws->res->res_g, nv);
            dvecmul(cws->nc, qp->d_mask, 0, ws->res->res_d, 0, ws->res->res_d, 0);
            dvecmul(cws->nc, qp->d_mask, 0, ws->res->res_m, 0, ws->res->res_m, 0);
        }
        d_dense_qp_res_compute_inf_norm(ws->res);
        cws->mu = ws->res->res_mu;
        // save infinity norm of residuals
        if (kk + 1 < ws->stat_max) {
            stat[ws->stat_m * (kk + 1) + 5] = ws->res->res_mu;
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
    //     else if (cws->mu != cws->mu) {
    //         // NaN in the solution
    //         ws->status = NAN_sol;
    //     }

    else {
        // normal return
        ws->status = SUCCESS;
    }

call_return:

    // TODO recover original pi

    if (arg->remove_lin_dep_eq) {
        d_dense_qp_restore_lin_dep_eq(qp, arg, ws);
    }

    // compute obj
    if (arg->compute_obj) {
        d_dense_qp_compute_obj(qp, qp_sol, arg, ws);
    }

#if 0
	d_dense_qp_sol_print(qp->dim, qp_sol);
exit(1);
#endif

    // return
}


void d_dense_qp_ipm_predict(struct d_dense_qp* qp, struct d_dense_qp_sol* qp_sol, struct d_dense_qp_ipm_arg* arg, struct d_dense_qp_ipm_ws* ws) {

#if 0
	d_dense_qp_dim_print(qp->dim);
	d_dense_qp_print(qp->dim, qp);
#endif

    int ii;

    struct d_core_qp_ipm_workspace* cws = ws->core_workspace;

    // arg to core workspace
    cws->lam_min = arg->lam_min;
    cws->t_min = arg->t_min;
    cws->t_min_inv = arg->t_min > 0.0 ? 1.0 / arg->t_min : 1e30;
    cws->tau_min = arg->tau_min;
    cws->t_lam_min = arg->t_lam_min;

    // alias qp vectors into qp_sol
    cws->v = qp_sol->v->pa;
    cws->pi = qp_sol->pi->pa;
    cws->lam = qp_sol->lam->pa;
    cws->t = qp_sol->t->pa;

    // alias members of qp_step
    ws->qp_step->dim = qp->dim;
    ws->qp_step->Hv = qp->Hv;
    ws->qp_step->A = qp->A;
    ws->qp_step->Ct = qp->Ct;
    ws->qp_step->Z = qp->Z;
    ws->qp_step->idxb = qp->idxb;
    ws->qp_step->idxs_rev = qp->idxs_rev;
    ws->qp_step->gz = ws->res->res_g;
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

#if 0
// compute residuals
d_dense_qp_res_compute(qp, qp_sol, ws->res, ws->res_ws);

printf("\npredict\t%e\t%e\t%e\t%e\n", qp_res[0], qp_res[1], qp_res[2], qp_res[3]);
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

    // TODO absolute formulation !!!!!

    // TODO robust formulation !!!!!

    // compute residuals
    d_dense_qp_res_compute(qp, qp_sol, ws->res, ws->res_ws);

    // printf("\npredict\t%e\t%e\t%e\t%e\n", qp_res[0], qp_res[1], qp_res[2], qp_res[3]);

    // solve kkt
    d_solve_kkt_step_dense_qp(ws->qp_step, ws->sol_step, arg, ws);

    // alpha TODO fix alpha=1 !!!!!
    //	d_compute_alpha_qp(cws->dlam, cws->dt, cws);
    cws->alpha = 1.0;

    //
    d_update_var_qp(cws);

    if (arg->comp_res_pred) {
        // compute residuals in exit
        d_dense_qp_res_compute(qp, qp_sol, ws->res, ws->res_ws);
    }

    // printf("\npredict\t%e\t%e\t%e\t%e\n", qp_res[0], qp_res[1], qp_res[2], qp_res[3]);

    // TODO

    // do not change status
}


void d_dense_qp_ipm_sens(struct d_dense_qp* qp, struct d_dense_qp_sol* qp_sol, struct d_dense_qp_ipm_arg* arg, struct d_dense_qp_ipm_ws* ws) {

#if 0
	d_dense_qp_dim_print(qp->dim);
	d_dense_qp_print(qp->dim, qp);
#endif

    int ii;

    struct d_core_qp_ipm_workspace* cws = ws->core_workspace;

    // arg to core workspace
    cws->lam_min = arg->lam_min;
    cws->t_min = arg->t_min;
    cws->t_min_inv = arg->t_min > 0.0 ? 1.0 / arg->t_min : 1e30;
    cws->tau_min = arg->tau_min;
    cws->t_lam_min = arg->t_lam_min;

    // alias qp vectors into qp_sol
    cws->v = qp_sol->v->pa;
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

    // solve kkt
    d_solve_kkt_step_dense_qp(qp, qp_sol, arg, ws);

#if 0
	// alpha TODO fix alpha=1 !!!!!
//	d_compute_alpha_qp(cws->dlam, cws->dt, cws);
	cws->alpha = 1.0;

	//
	d_update_var_qp(cws);

	if(arg->comp_res_pred)
		{
		// compute residuals in exit
		d_dense_qp_res_compute(qp, qp_sol, ws->res, ws->res_ws);
		}
#endif

    // printf("\npredict\t%e\t%e\t%e\t%e\n", qp_res[0], qp_res[1], qp_res[2], qp_res[3]);

    // TODO

    // do not change status
}


void d_dense_qp_compute_step_length(struct d_dense_qp* qp, struct d_dense_qp_sol* qp_sol, struct d_dense_qp_ipm_arg* arg, struct d_dense_qp_ipm_ws* ws) {


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

    struct d_core_qp_ipm_workspace* cws = ws->core_workspace;
    struct d_dense_qp_sol* qp_sol_step = ws->sol_step;
    struct d_dense_qp_res* res = ws->res;
    struct d_dense_qp_res* res_step = ws->res_step;

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

    struct vec* dv = qp_sol_step->v;
    struct vec* dpi = qp_sol_step->pi;
    struct vec* dlam = qp_sol_step->lam;
    struct vec* dt = qp_sol_step->t;

    struct vec* res_g = res->res_g;
    struct vec* res_b = res->res_b;
    struct vec* res_d = res->res_d;
    struct vec* res_m = res->res_m;

    struct vec* step_res_g_p = res_step->res_g;
    struct vec* step_res_g_d = ws->tmp_nv2ns;
    struct vec* step_res_b = res_step->res_b;
    struct vec* step_res_d = res_step->res_d;
    struct vec* step_res_m = res_step->res_m;

    struct vec* tmp_nbg = ws->res_ws->tmp_nbg;
    struct vec* tmp_ns = ws->res_ws->tmp_ns;
    struct vec* tmp_nv2ns = ws->tmp_nv2ns + 1;

    double mu, tmp;
    double phi_H, phi_g, phi_rho, phi_obj0, phi_obj1;

    double alpha_stat, alpha0, alpha1, alpha_p, alpha_d;

    int ii, idx;

    // compute (strictly) linear (i.e. no constr) part of res wrt step

    // res g
    dsymv_l(nv, 1.0, Hg, 0, 0, dv, 0, 0.0, step_res_g_p, 0, step_res_g_p, 0);

    if (nb + ng > 0) {
        daxpy(nb + ng, -1.0, dlam, 0, dlam, nb + ng, tmp_nbg + 0, 0);
        //		daxpy(nb+ng,  1.0, d, 0, t, 0, res_d, 0);
        //		daxpy(nb+ng,  1.0, d, nb+ng, t, nb+ng, res_d, nb+ng);
        //		daxpy(2*nb+2*ng,  1.0, d, 0, t, 0, res_d, 0);
        dveccp(2 * nb + 2 * ng, dt, 0, step_res_d, 0);
        // box
        if (nb > 0) {
            dvecse(nv, 0.0, step_res_g_d, 0);
            dvecad_sp(nb, 1.0, tmp_nbg + 0, 0, idxb, step_res_g_d, 0);
            dvecex_sp(nb, 1.0, idxb, dv, 0, tmp_nbg + 1, 0);
        }
        // general
        if (ng > 0) {
            dgemv_nt(nv, ng, 1.0, 1.0, Ct, 0, 0, tmp_nbg + 0, nb, dv, 0, 1.0, 0.0, step_res_g_d, 0, tmp_nbg + 1, nb, step_res_g_d, 0, tmp_nbg + 1, nb);
        }
        daxpy(nb + ng, -1.0, tmp_nbg + 1, 0, step_res_d, 0, step_res_d, 0);
        daxpy(nb + ng, 1.0, tmp_nbg + 1, 0, step_res_d, nb + ng, step_res_d, nb + ng);
    }
    if (ns > 0) {
        // res_g
        dgemv_d(2 * ns, 1.0, Z, 0, dv, nv, 0.0, step_res_g_p, nv, step_res_g_p, nv);
        //		daxpy(2*ns, -1.0, dlam, 2*nb+2*ng, step_res_g, nv, step_res_g, nv);
        dveccpsc(2 * ns, -1.0, dlam, 2 * nb + 2 * ng, step_res_g_d, nv);
        for (ii = 0; ii < nb + ng; ii++) {
            idx = idxs_rev[ii];
            if (idx != -1) {
                VECEL(step_res_g_d, nv + idx) -= VECEL(dlam, ii);
                VECEL(step_res_g_d, nv + ns + idx) -= VECEL(dlam, nb + ng + ii);
                // res_d
                VECEL(step_res_d, ii) -= VECEL(dv + ii, nv + idx);
                VECEL(step_res_d, nb + ng + ii) -= VECEL(dv + ii, nv + ns + idx);
            }
        }
        // res_d
        daxpy(2 * ns, -1.0, dv, nv, dt, 2 * nb + 2 * ng, step_res_d, 2 * nb + 2 * ng);
        //		daxpy(2*ns, 1.0, d, 2*nb+2*ng, res_d, 2*nb+2*ng, res_d, 2*nb+2*ng);
    }

    // res b, res g
    if (ne > 0)
        dgemv_nt(ne, nv, -1.0, -1.0, A, 0, 0, dv, 0, dpi, 0, 0.0, 1.0, step_res_b, 0, step_res_g_d, 0, step_res_b, 0, step_res_g_d, 0);

        // compute alpha QP

#if 1

    alpha_p = cws->alpha_prim;
    alpha_d = cws->alpha_dual;

    if (alpha_p < alpha_d) {
        // fix alpha_p
        daxpy(nv + 2 * ns, alpha_p, step_res_g_p, 0, res_g, 0, tmp_nv2ns, 0);
        // TODO
        phi_H = 0.0;
        phi_H += ddot(nv + 2 * ns, step_res_g_d, 0, step_res_g_d, 0);
        // TODO
        phi_g = 0.0;
        phi_g += 2 * ddot(nv + 2 * ns, step_res_g_d, 0, tmp_nv2ns, 0);
        phi_g += ddot(2 * nb + 2 * ng + 2 * ns, dlam, 0, t, 0);
        phi_g += alpha_p * ddot(2 * nb + 2 * ng + 2 * ns, dlam, 0, dt, 0);
        // TODO
        phi_rho = 0.0;
        phi_rho += ddot(nv + 2 * ns, tmp_nv2ns, 0, tmp_nv2ns, 0);
        phi_rho += alpha_p * alpha_p * ddot(ne, step_res_b, 0, step_res_b, 0);
        phi_rho += alpha_p * 2 * ddot(ne, step_res_b, 0, res_b, 0);
        phi_rho += ddot(ne, res_b, 0, res_b, 0);
        phi_rho += alpha_p * alpha_p * ddot(2 * nb + 2 * ng + 2 * ns, step_res_d, 0, step_res_d, 0);
        phi_rho += alpha_p * 2 * ddot(2 * nb + 2 * ng + 2 * ns, step_res_d, 0, res_d, 0);
        phi_rho += ddot(2 * nb + 2 * ng + 2 * ns, res_d, 0, res_d, 0);
        phi_rho += ddot(2 * nb + 2 * ng + 2 * ns, lam, 0, t, 0);
        phi_rho += alpha_p * ddot(2 * nb + 2 * ng + 2 * ns, lam, 0, dt, 0);
        // TODO
        // unconstr minimizer
        alpha_stat = -0.5 * phi_g / phi_H;
        // clipping
        alpha1 = alpha_stat > alpha_p ? alpha_stat : alpha_p;
        alpha1 = alpha1 < alpha_d ? alpha1 : alpha_d;
        // obj
        phi_obj1 = phi_H * alpha1 * alpha1 + phi_g * alpha1 + phi_rho;
    } else {
        // fix alpha_d
        daxpy(nv + 2 * ns, alpha_d, step_res_g_d, 0, res_g, 0, tmp_nv2ns, 0);
        // TODO
        phi_H = 0.0;
        phi_H += ddot(nv + 2 * ns, step_res_g_p, 0, step_res_g_p, 0);
        phi_H += ddot(ne, step_res_b, 0, step_res_b, 0);
        phi_H += ddot(2 * nb + 2 * ng + 2 * ns, step_res_d, 0, step_res_d, 0);
        // TODO
        phi_g = 0.0;
        phi_g += 2 * ddot(nv + 2 * ns, step_res_g_p, 0, tmp_nv2ns, 0);
        phi_g += 2 * ddot(ne, step_res_b, 0, res_b, 0);
        phi_g += 2 * ddot(2 * nb + 2 * ng + 2 * ns, step_res_d, 0, res_d, 0);
        phi_g += ddot(2 * nb + 2 * ng + 2 * ns, lam, 0, dt, 0);
        phi_g += alpha_d * ddot(2 * nb + 2 * ng + 2 * ns, dlam, 0, dt, 0);
        // TODO
        phi_rho = 0.0;
        phi_rho += ddot(nv + 2 * ns, tmp_nv2ns, 0, tmp_nv2ns, 0);
        phi_rho += ddot(ne, res_b, 0, res_b, 0);
        phi_rho += ddot(2 * nb + 2 * ng + 2 * ns, res_d, 0, res_d, 0);
        phi_rho += ddot(2 * nb + 2 * ng + 2 * ns, lam, 0, t, 0);
        phi_rho += alpha_d * ddot(2 * nb + 2 * ng + 2 * ns, dlam, 0, t, 0);
        // TODO
        // unconstr minimizer
        alpha_stat = -0.5 * phi_g / phi_H;
        // clipping
        alpha1 = alpha_stat > alpha_d ? alpha_stat : alpha_d;
        alpha1 = alpha1 < alpha_p ? alpha1 : alpha_p;
        // obj
        phi_obj1 = phi_H * alpha1 * alpha1 + phi_g * alpha1 + phi_rho;
    }

#endif

    daxpy(nv + 2 * ns, 1.0, step_res_g_d, 0, step_res_g_p, 0, step_res_g_p, 0);

    phi_H = 0.0;
    phi_H += ddot(nv + 2 * ns, step_res_g_p, 0, step_res_g_p, 0);
    phi_H += ddot(ne, step_res_b, 0, step_res_b, 0);
    phi_H += ddot(2 * nb + 2 * ng + 2 * ns, step_res_d, 0, step_res_d, 0);
    phi_H += ddot(2 * nb + 2 * ng + 2 * ns, dlam, 0, dt, 0);

    phi_g = 0.0;
    phi_g += 2 * ddot(nv + 2 * ns, step_res_g_p, 0, res_g, 0);
    phi_g += 2 * ddot(ne, step_res_b, 0, res_b, 0);
    phi_g += 2 * ddot(2 * nb + 2 * ng + 2 * ns, step_res_d, 0, res_d, 0);
    phi_g += ddot(2 * nb + 2 * ng + 2 * ns, dlam, 0, t, 0);
    phi_g += ddot(2 * nb + 2 * ng + 2 * ns, lam, 0, dt, 0);

    phi_rho = 0.0;
    phi_rho += ddot(nv + 2 * ns, res_g, 0, res_g, 0);
    phi_rho += ddot(ne, res_b, 0, res_b, 0);
    phi_rho += ddot(2 * nb + 2 * ng + 2 * ns, res_d, 0, res_d, 0);
    phi_rho += ddot(2 * nb + 2 * ng + 2 * ns, lam, 0, t, 0);

    alpha_stat = -0.5 * phi_g / phi_H;

    double step_min = 1e-3;  // minimum step length

    if (phi_H > 0)  // alpha_stat is a minimum
    {
        alpha0 = alpha_stat > step_min ? alpha_stat : step_min;
        alpha0 = alpha0 < cws->alpha ? alpha0 : cws->alpha;
    } else  // alpha_stat is a maximum
    {
        alpha0 = cws->alpha;
        //		if(alpha_stat<=0.0)
        //			{
        //			printf("\nmax < 0\n");
        //			alpha0 = cws->alpha;
        //			}
        //		else if(alpha_stat>=cws->alpha)
        //			{
        //			alpha0 = 0.0;
        //			}
        //		else // 0 < alpha_stat < alpha
        //			{
        //			phi_obj0 = phi_rho;
        //			phi_obj1 = phi_H*cws->alpha*cws->alpha + phi_g*cws->alpha + phi_rho;
        //			printf("\nphi obj %e %e\n", phi_obj0, phi_obj1);
        //			alpha0 = phi_obj0<phi_obj1 ? 0.0 : cws->alpha;
        //			alpha0 = cws->alpha;
        //			}
    }
    // obj
    phi_obj0 = phi_H * alpha0 * alpha0 + phi_g * alpha0 + phi_rho;

    cws->alpha = alpha0;

    if (phi_obj0 <= phi_obj1) {
        cws->alpha_prim = alpha0;
        cws->alpha_dual = alpha0;
    } else {
        if (alpha_p < alpha_d) {
            //			cws->alpha_prim = alpha_p;
            cws->alpha_dual = alpha1;
        } else {
            // TODO
            cws->alpha_prim = alpha1;
            //			cws->alpha_dual = alpha_d;
        }
    }

    //	printf("\nobj %e %e alpha %e %e %e %e\n", phi_obj0, phi_obj1, alpha0, alpha_p, alpha_d, alpha1);

    // res_m res_mu
    //	mu = dvecmuldot(nct, lam, 0, t, 0, res_m, 0);
    //	daxpy(nct, -1.0, m, 0, res_m, 0, res_m, 0);
    // TODO use nc_mask_inv from cws if available !!!!!
    //	res->res_mu = mu*nct_inv;
}
