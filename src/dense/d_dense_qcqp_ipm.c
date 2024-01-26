#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "tinyhpipm/blas.h"
#include "tinyhpipm/dense/d_dense_qcqp.h"
#include "tinyhpipm/dense/d_dense_qcqp_dim.h"
#include "tinyhpipm/dense/d_dense_qcqp_ipm.h"
#include "tinyhpipm/dense/d_dense_qcqp_res.h"
#include "tinyhpipm/dense/d_dense_qcqp_sol.h"
#include "tinyhpipm/dense/d_dense_qcqp_utils.h"
#include "tinyhpipm/dense/d_dense_qp.h"
#include "tinyhpipm/dense/d_dense_qp_dim.h"
#include "tinyhpipm/dense/d_dense_qp_ipm.h"
#include "tinyhpipm/dense/d_dense_qp_kkt.h"
#include "tinyhpipm/dense/d_dense_qp_res.h"
#include "tinyhpipm/dense/d_dense_qp_sol.h"
#include "tinyhpipm/dense/d_dense_qp_utils.h"
#include "tinyhpipm/ipm_core/d_core_qp_ipm.h"
#include "tinyhpipm/ipm_core/d_core_qp_ipm_aux.h"


hpipm_size_t d_dense_qcqp_ipm_arg_strsize() {
    return sizeof(struct d_dense_qcqp_ipm_arg);
}


hpipm_size_t d_dense_qcqp_ipm_arg_memsize(struct d_dense_qcqp_dim* dim) {

    hpipm_size_t size = 0;

    size += 1 * sizeof(struct d_dense_qp_ipm_arg);
    size += 1 * d_dense_qp_ipm_arg_memsize(dim->qp_dim);

    size = (size + 63) / 64 * 64;  // make multiple of typical cache line size
    size += 1 * 64;  // align once to typical cache line size

    return size;
}


void d_dense_qcqp_ipm_arg_create(struct d_dense_qcqp_dim* dim, struct d_dense_qcqp_ipm_arg* arg, void* mem) {
    // loop index
    int ii;

    // zero memory (to avoid corrupted memory like e.g. NaN)
    hpipm_size_t memsize = d_dense_qcqp_ipm_arg_memsize(dim);
    hpipm_zero_memset(memsize, mem);

    // qp_dim struct
    struct d_dense_qp_ipm_arg* arg_ptr = mem;

    arg->qp_arg = arg_ptr;
    arg_ptr += 1;

    // align to typical cache line size
    hpipm_size_t s_ptr = (hpipm_size_t) arg_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;

    // void
    char* c_ptr = (char*) s_ptr;

    d_dense_qp_ipm_arg_create(dim->qp_dim, arg->qp_arg, c_ptr);
    c_ptr += arg->qp_arg->memsize;


    arg->memsize = d_dense_qcqp_ipm_arg_memsize(dim);

#if defined(RUNTIME_CHECKS)
    if (c_ptr > ((char*) mem) + arg->memsize) {
        printf("\nerror: d_dense_qcqp_ipm_arg_create: outside memory bounds!\n\n");
        exit(1);
    }
#endif
}


void d_dense_qcqp_ipm_arg_set_default(enum tinyhpipm_mode mode, struct d_dense_qcqp_ipm_arg* arg) {

    d_dense_qp_ipm_arg_set_default(mode, arg->qp_arg);

    double mu0, alpha_min, res_g, res_b, res_d, res_m, reg_prim, reg_dual, lam_min, t_min;
    int iter_max, stat_max, pred_corr, cond_pred_corr, itref_pred_max, itref_corr_max, lq_fact, scale, warm_start, abs_form, comp_res_exit, comp_res_pred, split_step, t_lam_min;

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
        warm_start = 0;
        abs_form = 1;
        comp_res_exit = 0;
        comp_res_pred = 0;
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
        warm_start = 0;
        abs_form = 0;
        comp_res_exit = 1;
        comp_res_pred = 0;
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
        warm_start = 0;
        abs_form = 0;
        comp_res_exit = 1;
        comp_res_pred = 0;
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
        warm_start = 0;
        abs_form = 0;
        comp_res_exit = 1;
        comp_res_pred = 0;
        split_step = 0;
        t_lam_min = 2;
    } else {
        printf("\nerror: d_dense_qp_ipm_arg_set_default: wrong set default mode\n");
        exit(1);
    }

    // use individual setters when available
    d_dense_qcqp_ipm_arg_set_mu0(&mu0, arg);
    d_dense_qcqp_ipm_arg_set_alpha_min(&alpha_min, arg);
    d_dense_qcqp_ipm_arg_set_tol_stat(&res_g, arg);
    d_dense_qcqp_ipm_arg_set_tol_eq(&res_b, arg);
    d_dense_qcqp_ipm_arg_set_tol_ineq(&res_d, arg);
    d_dense_qcqp_ipm_arg_set_tol_comp(&res_m, arg);
    d_dense_qcqp_ipm_arg_set_iter_max(&iter_max, arg);
    arg->stat_max = stat_max;
    d_dense_qcqp_ipm_arg_set_pred_corr(&pred_corr, arg);
    d_dense_qcqp_ipm_arg_set_cond_pred_corr(&cond_pred_corr, arg);
    arg->itref_pred_max = itref_pred_max;
    arg->itref_corr_max = itref_corr_max;
    d_dense_qcqp_ipm_arg_set_reg_prim(&reg_prim, arg);
    d_dense_qcqp_ipm_arg_set_reg_dual(&reg_prim, arg);
    arg->lq_fact = lq_fact;
    arg->scale = scale;
    arg->lam_min = lam_min;
    arg->t_min = t_min;
    d_dense_qcqp_ipm_arg_set_lam_min(&lam_min, arg);
    d_dense_qcqp_ipm_arg_set_t_min(&t_min, arg);
    d_dense_qcqp_ipm_arg_set_warm_start(&warm_start, arg);
    arg->abs_form = abs_form;
    d_dense_qcqp_ipm_arg_set_comp_res_pred(&comp_res_pred, arg);
    d_dense_qcqp_ipm_arg_set_comp_res_exit(&comp_res_pred, arg);
    d_dense_qcqp_ipm_arg_set_split_step(&split_step, arg);
    d_dense_qcqp_ipm_arg_set_t_lam_min(&t_lam_min, arg);
    arg->mode = mode;
}


void d_dense_qcqp_ipm_arg_set(char* field, void* value, struct d_dense_qcqp_ipm_arg* arg) {
    if (hpipm_strcmp(field, "iter_max")) {
        d_dense_qcqp_ipm_arg_set_iter_max(value, arg);
    } else if (hpipm_strcmp(field, "alpha_min")) {
        d_dense_qcqp_ipm_arg_set_alpha_min(value, arg);
    } else if (hpipm_strcmp(field, "mu0")) {
        d_dense_qcqp_ipm_arg_set_mu0(value, arg);
    } else if (hpipm_strcmp(field, "tol_stat")) {
        d_dense_qcqp_ipm_arg_set_tol_stat(value, arg);
    } else if (hpipm_strcmp(field, "tol_eq")) {
        d_dense_qcqp_ipm_arg_set_tol_eq(value, arg);
    } else if (hpipm_strcmp(field, "tol_ineq")) {
        d_dense_qcqp_ipm_arg_set_tol_ineq(value, arg);
    } else if (hpipm_strcmp(field, "tol_comp")) {
        d_dense_qcqp_ipm_arg_set_tol_comp(value, arg);
    } else if (hpipm_strcmp(field, "reg_prim")) {
        d_dense_qcqp_ipm_arg_set_reg_prim(value, arg);
    } else if (hpipm_strcmp(field, "reg_dual")) {
        d_dense_qcqp_ipm_arg_set_reg_dual(value, arg);
    } else if (hpipm_strcmp(field, "warm_start")) {
        d_dense_qcqp_ipm_arg_set_warm_start(value, arg);
    } else if (hpipm_strcmp(field, "pred_corr")) {
        d_dense_qcqp_ipm_arg_set_pred_corr(value, arg);
    } else if (hpipm_strcmp(field, "cond_pred_corr")) {
        d_dense_qcqp_ipm_arg_set_cond_pred_corr(value, arg);
    } else if (hpipm_strcmp(field, "comp_res_pred")) {
        d_dense_qcqp_ipm_arg_set_comp_res_pred(value, arg);
    } else if (hpipm_strcmp(field, "comp_res_exit")) {
        d_dense_qcqp_ipm_arg_set_comp_res_exit(value, arg);
    } else if (hpipm_strcmp(field, "lam_min")) {
        d_dense_qcqp_ipm_arg_set_lam_min(value, arg);
    } else if (hpipm_strcmp(field, "t_min")) {
        d_dense_qcqp_ipm_arg_set_t_min(value, arg);
    } else if (hpipm_strcmp(field, "split_step")) {
        d_dense_qcqp_ipm_arg_set_split_step(value, arg);
    } else if (hpipm_strcmp(field, "t_lam_min")) {
        d_dense_qcqp_ipm_arg_set_t_lam_min(value, arg);
    } else {
        printf("error: d_dense_qcqp_ipm_arg_set: wrong field %s\n", field);
        exit(1);
    }
}


void d_dense_qcqp_ipm_arg_set_iter_max(int* value, struct d_dense_qcqp_ipm_arg* arg) {
    arg->iter_max = *value;
    d_dense_qp_ipm_arg_set_iter_max(value, arg->qp_arg);
}


void d_dense_qcqp_ipm_arg_set_alpha_min(double* value, struct d_dense_qcqp_ipm_arg* arg) {
    arg->alpha_min = *value;
    d_dense_qp_ipm_arg_set_alpha_min(value, arg->qp_arg);
}


void d_dense_qcqp_ipm_arg_set_mu0(double* value, struct d_dense_qcqp_ipm_arg* arg) {
    arg->mu0 = *value;
    d_dense_qp_ipm_arg_set_mu0(value, arg->qp_arg);
}


void d_dense_qcqp_ipm_arg_set_tol_stat(double* value, struct d_dense_qcqp_ipm_arg* arg) {
    arg->res_g_max = *value;
    d_dense_qp_ipm_arg_set_tol_stat(value, arg->qp_arg);
}


void d_dense_qcqp_ipm_arg_set_tol_eq(double* value, struct d_dense_qcqp_ipm_arg* arg) {
    arg->res_b_max = *value;
    d_dense_qp_ipm_arg_set_tol_eq(value, arg->qp_arg);
}


void d_dense_qcqp_ipm_arg_set_tol_ineq(double* value, struct d_dense_qcqp_ipm_arg* arg) {
    arg->res_d_max = *value;
    d_dense_qp_ipm_arg_set_tol_ineq(value, arg->qp_arg);
}


void d_dense_qcqp_ipm_arg_set_tol_comp(double* value, struct d_dense_qcqp_ipm_arg* arg) {
    arg->res_m_max = *value;
    d_dense_qp_ipm_arg_set_tol_comp(value, arg->qp_arg);
}


void d_dense_qcqp_ipm_arg_set_reg_prim(double* value, struct d_dense_qcqp_ipm_arg* arg) {
    arg->reg_prim = *value;
    d_dense_qp_ipm_arg_set_reg_prim(value, arg->qp_arg);
}


void d_dense_qcqp_ipm_arg_set_reg_dual(double* value, struct d_dense_qcqp_ipm_arg* arg) {
    arg->reg_dual = *value;
    d_dense_qp_ipm_arg_set_reg_dual(value, arg->qp_arg);
}


void d_dense_qcqp_ipm_arg_set_warm_start(int* value, struct d_dense_qcqp_ipm_arg* arg) {
    arg->warm_start = *value;
    d_dense_qp_ipm_arg_set_warm_start(value, arg->qp_arg);
}


void d_dense_qcqp_ipm_arg_set_pred_corr(int* value, struct d_dense_qcqp_ipm_arg* arg) {
    arg->pred_corr = *value;
    d_dense_qp_ipm_arg_set_pred_corr(value, arg->qp_arg);
}


void d_dense_qcqp_ipm_arg_set_cond_pred_corr(int* value, struct d_dense_qcqp_ipm_arg* arg) {
    arg->cond_pred_corr = *value;
    d_dense_qp_ipm_arg_set_cond_pred_corr(value, arg->qp_arg);
}


void d_dense_qcqp_ipm_arg_set_comp_res_exit(int* value, struct d_dense_qcqp_ipm_arg* arg) {
    arg->comp_res_exit = *value;
    d_dense_qp_ipm_arg_set_comp_res_exit(value, arg->qp_arg);
}


void d_dense_qcqp_ipm_arg_set_comp_res_pred(int* value, struct d_dense_qcqp_ipm_arg* arg) {
    arg->comp_res_pred = *value;
    d_dense_qp_ipm_arg_set_comp_res_pred(value, arg->qp_arg);
}


void d_dense_qcqp_ipm_arg_set_lam_min(double* value, struct d_dense_qcqp_ipm_arg* arg) {
    arg->lam_min = *value;
    d_dense_qp_ipm_arg_set_lam_min(value, arg->qp_arg);
}


void d_dense_qcqp_ipm_arg_set_t_min(double* value, struct d_dense_qcqp_ipm_arg* arg) {
    arg->t_min = *value;
    d_dense_qp_ipm_arg_set_t_min(value, arg->qp_arg);
}


void d_dense_qcqp_ipm_arg_set_split_step(int* value, struct d_dense_qcqp_ipm_arg* arg) {
    arg->split_step = *value;
    d_dense_qp_ipm_arg_set_split_step(value, arg->qp_arg);
}


void d_dense_qcqp_ipm_arg_set_t_lam_min(int* value, struct d_dense_qcqp_ipm_arg* arg) {
    arg->t_lam_min = *value;
    d_dense_qp_ipm_arg_set_t_lam_min(value, arg->qp_arg);
}

hpipm_size_t d_dense_qcqp_ipm_ws_strsize() {
    return sizeof(struct d_dense_qcqp_ipm_ws);
}


hpipm_size_t d_dense_qcqp_ipm_ws_memsize(struct d_dense_qcqp_dim* dim, struct d_dense_qcqp_ipm_arg* arg) {

    int nv = dim->nv;

    hpipm_size_t size = 0;

    size += 1 * sizeof(struct d_dense_qp_ipm_ws);
    size += 1 * d_dense_qp_ipm_ws_memsize(dim->qp_dim, arg->qp_arg);

    size += 1 * sizeof(struct d_dense_qcqp_res_ws);  // qcqp_res_ws
    size += 1 * d_dense_qcqp_res_ws_memsize(dim);  // qcqp_res_ws

    size += 1 * sizeof(struct d_dense_qcqp_res);  // qcqp_res
    size += 1 * d_dense_qcqp_res_memsize(dim);  // qcqp_res

    size += 1 * sizeof(struct d_dense_qp);  // qp
    size += 1 * d_dense_qp_memsize(dim->qp_dim);  // qp

    size += 1 * sizeof(struct d_dense_qp_sol);  // qp_sol
    size += 1 * d_dense_qp_sol_memsize(dim->qp_dim);  // qp_sol

    size += 2 * sizeof(struct vec);  // tmp_nv
    size += 2 * memsize_vec(nv);  // tmp_nv

    size = (size + 63) / 64 * 64;  // make multiple of typical cache line size
    size += 1 * 64;  // align once to typical cache line size

    return size;
}


void d_dense_qcqp_ipm_ws_create(struct d_dense_qcqp_dim* dim, struct d_dense_qcqp_ipm_arg* arg, struct d_dense_qcqp_ipm_ws* workspace, void* mem) {

    int ii;

    // zero memory (to avoid corrupted memory like e.g. NaN)
    hpipm_size_t memsize = d_dense_qcqp_ipm_ws_memsize(dim, arg);
    hpipm_zero_memset(memsize, mem);

    int nv = dim->nv;

    char* c_ptr = mem;


    // structures
    workspace->qp_ws = (struct d_dense_qp_ipm_ws*) c_ptr;
    c_ptr += sizeof(struct d_dense_qp_ipm_ws);

    workspace->qp = (struct d_dense_qp*) c_ptr;
    c_ptr += sizeof(struct d_dense_qp);

    workspace->qp_sol = (struct d_dense_qp_sol*) c_ptr;
    c_ptr += sizeof(struct d_dense_qp_sol);

    workspace->qcqp_res_ws = (struct d_dense_qcqp_res_ws*) c_ptr;
    c_ptr += sizeof(struct d_dense_qcqp_res_ws);

    workspace->qcqp_res = (struct d_dense_qcqp_res*) c_ptr;
    c_ptr += sizeof(struct d_dense_qcqp_res);


    // vector struct
    struct vec* sv_ptr = (struct vec*) c_ptr;

    workspace->tmp_nv = sv_ptr;
    sv_ptr += 2;


    // align to typical cache line size
    hpipm_size_t s_ptr = (hpipm_size_t) sv_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;


    // memory of structures
    c_ptr = (char*) s_ptr;

    d_dense_qp_ipm_ws_create(dim->qp_dim, arg->qp_arg, workspace->qp_ws, c_ptr);
    c_ptr += workspace->qp_ws->memsize;

    d_dense_qp_create(dim->qp_dim, workspace->qp, c_ptr);
    c_ptr += workspace->qp->memsize;

    d_dense_qp_sol_create(dim->qp_dim, workspace->qp_sol, c_ptr);
    c_ptr += workspace->qp_sol->memsize;

    d_dense_qcqp_res_ws_create(dim, workspace->qcqp_res_ws, c_ptr);
    c_ptr += workspace->qcqp_res_ws->memsize;

    d_dense_qcqp_res_create(dim, workspace->qcqp_res, c_ptr);
    c_ptr += workspace->qcqp_res->memsize;

    create_vec(nv, workspace->tmp_nv + 0, c_ptr);
    c_ptr += (workspace->tmp_nv + 0)->memsize;
    create_vec(nv, workspace->tmp_nv + 1, c_ptr);
    c_ptr += (workspace->tmp_nv + 1)->memsize;


    //
    workspace->memsize = d_dense_qcqp_ipm_ws_memsize(dim, arg);


#if defined(RUNTIME_CHECKS)
    if (c_ptr > ((char*) mem) + workspace->memsize) {
        printf("\nCreate_dense_qp_ipm: outside memory bounds!\n\n");
        exit(1);
    }
#endif
}


void d_dense_qcqp_ipm_get(char* field, struct d_dense_qcqp_ipm_ws* ws, void* value) {
    if (hpipm_strcmp(field, "status")) {
        d_dense_qcqp_ipm_get_status(ws, value);
    } else if (hpipm_strcmp(field, "iter")) {
        d_dense_qcqp_ipm_get_iter(ws, value);
    } else if (hpipm_strcmp(field, "max_res_stat")) {
        d_dense_qcqp_ipm_get_max_res_stat(ws, value);
    } else if (hpipm_strcmp(field, "max_res_eq")) {
        d_dense_qcqp_ipm_get_max_res_eq(ws, value);
    } else if (hpipm_strcmp(field, "max_res_ineq")) {
        d_dense_qcqp_ipm_get_max_res_ineq(ws, value);
    } else if (hpipm_strcmp(field, "max_res_comp")) {
        d_dense_qcqp_ipm_get_max_res_comp(ws, value);
    } else if (hpipm_strcmp(field, "obj")) {
        d_dense_qcqp_ipm_get_obj(ws, value);
    } else if (hpipm_strcmp(field, "stat")) {
        d_dense_qcqp_ipm_get_stat(ws, value);
    } else if (hpipm_strcmp(field, "stat_m")) {
        d_dense_qcqp_ipm_get_stat_m(ws, value);
    } else {
        printf("error: d_dense_qcqp_ipm_get: wrong field %s\n", field);
        exit(1);
    }
}


void d_dense_qcqp_ipm_get_status(struct d_dense_qcqp_ipm_ws* ws, int* status) {
    *status = ws->status;
}


void d_dense_qcqp_ipm_get_iter(struct d_dense_qcqp_ipm_ws* ws, int* iter) {
    *iter = ws->iter;
}


void d_dense_qcqp_ipm_get_max_res_stat(struct d_dense_qcqp_ipm_ws* ws, double* res_stat) {
    *res_stat = ws->qcqp_res->res_max[0];
}


void d_dense_qcqp_ipm_get_max_res_eq(struct d_dense_qcqp_ipm_ws* ws, double* res_eq) {
    *res_eq = ws->qcqp_res->res_max[1];
}


void d_dense_qcqp_ipm_get_max_res_ineq(struct d_dense_qcqp_ipm_ws* ws, double* res_ineq) {
    *res_ineq = ws->qcqp_res->res_max[2];
}


void d_dense_qcqp_ipm_get_max_res_comp(struct d_dense_qcqp_ipm_ws* ws, double* res_comp) {
    *res_comp = ws->qcqp_res->res_max[3];
}


void d_dense_qcqp_ipm_get_obj(struct d_dense_qcqp_ipm_ws* ws, double* obj) {
    *obj = ws->qcqp_res->obj;
}


void d_dense_qcqp_ipm_get_stat(struct d_dense_qcqp_ipm_ws* ws, double** stat) {
    d_dense_qp_ipm_get_stat(ws->qp_ws, stat);
}


void d_dense_qcqp_ipm_get_stat_m(struct d_dense_qcqp_ipm_ws* ws, int* stat_m) {
    d_dense_qp_ipm_get_stat_m(ws->qp_ws, stat_m);
}


void d_dense_qcqp_init_var(struct d_dense_qcqp* qcqp, struct d_dense_qcqp_sol* qcqp_sol, struct d_dense_qcqp_ipm_arg* arg, struct d_dense_qcqp_ipm_ws* ws) {

    //	struct CORE_qcqp_ipm_workspace *cws = ws->core_workspace;

    // extract cws members
    int nv = qcqp->dim->nv;
    int ne = qcqp->dim->ne;
    int nb = qcqp->dim->nb;
    int ng = qcqp->dim->ng;
    int nq = qcqp->dim->nq;
    int ns = qcqp->dim->ns;

    double* d = qcqp->d->pa;
    int* idxb = qcqp->idxb;

    double* v = qcqp_sol->v->pa;
    double* pi = qcqp_sol->pi->pa;
    double* lam = qcqp_sol->lam->pa;
    double* t = qcqp_sol->t->pa;

    double mu0 = arg->mu0;

    // local variables
    int ii;
    int idxb0;
    double tmp;

    // TODO move to args ???
    double thr0 = 0.5;


    // primal and dual variables
    if (arg->warm_start == 2) {

        thr0 = 1e-1;

        for (ii = 0; ii < 2 * nb + 2 * ng + 2 * nq + 2 * ns; ii++) {
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
        t[nb + ng + nq + ii] = -d[nb + ng + nq + ii] - v[idxb0];
        if (t[0 + ii] < thr0) {
            if (t[nb + ng + nq + ii] < thr0) {
                v[idxb0] = 0.5 * (d[0 + ii] + d[nb + ng + nq + ii]);
                t[0 + ii] = thr0;
                t[nb + ng + nq + ii] = thr0;
            } else {
                t[0 + ii] = thr0;
                v[idxb0] = d[0 + ii] + thr0;
            }
        } else if (t[nb + ng + nq + ii] < thr0) {
            t[nb + ng + nq + ii] = thr0;
            v[idxb0] = -d[nb + ng + nq + ii] - thr0;
        }
#else
        t[0 + ii] = 1.0;
        t[nb + ng + nq + ii] = 1.0;
#endif
        lam[0 + ii] = mu0 / t[0 + ii];
        lam[nb + ng + nq + ii] = mu0 / t[nb + ng + nq + ii];
    }

    // general constraints
    if (ng > 0) {
        dgemv_t(nv, ng, 1.0, qcqp->Ct, 0, 0, qcqp_sol->v, 0, 0.0, qcqp_sol->t, nb, qcqp_sol->t, nb);
        for (ii = 0; ii < ng; ii++) {
#if 1
            t[2 * nb + ng + nq + ii] = t[nb + ii];
            t[nb + ii] -= d[nb + ii];
            t[2 * nb + ng + nq + ii] -= d[2 * nb + ng + nq + ii];
            //		t[nb+ii]      = fmax( thr0, t[nb+ii] );
            //		t[2*nb+ng+nq+ii] = fmax( thr0, t[2*nb+ng+nq+ii] );
            t[nb + ii] = thr0 > t[nb + ii] ? thr0 : t[nb + ii];
            t[2 * nb + ng + nq + ii] = thr0 > t[2 * nb + ng + nq + ii] ? thr0 : t[2 * nb + ng + nq + ii];
#else
            t[nb + ii] = 1.0;
            t[2 * nb + ng + nq + ii] = 1.0;
#endif
            lam[nb + ii] = mu0 / t[nb + ii];
            lam[2 * nb + ng + nq + ii] = mu0 / t[2 * nb + ng + nq + ii];
        }
    }

    // soft constraints
    for (ii = 0; ii < ns; ii++) {
        t[2 * nb + 2 * ng + 2 * nq + ii] = 1.0;  // thr0;
        t[2 * nb + 2 * ng + 2 * nq + ns + ii] = 1.0;  // thr0;
        lam[2 * nb + 2 * ng + 2 * nq + ii] = mu0 / t[2 * nb + 2 * ng + 2 * nq + ii];
        lam[2 * nb + 2 * ng + 2 * nq + ns + ii] = mu0 / t[2 * nb + 2 * ng + 2 * nq + ns + ii];
    }

    //  quadratic constraints
    double sqrt_mu0 = sqrt(mu0);
    sqrt_mu0 = thr0 > sqrt_mu0 ? thr0 : sqrt_mu0;
    double mu0_div_sqrt_mu0 = mu0 / sqrt_mu0;

    for (ii = 0; ii < nq; ii++) {
        // disregard lower
        lam[nb + ng + ii] = 0.0;
        t[nb + ng + ii] = 1.0;
        // upper
#if 1
        t[2 * nb + 2 * ng + nq + ii] = sqrt_mu0;
        lam[2 * nb + 2 * ng + nq + ii] = mu0_div_sqrt_mu0;
#else
        //		t[2*nb+2*ng+nq+ii] = 1.0; // thr0;
        dcolex(nv, qcqp->Ct, 0, ng + ii, ws->tmp_nv, 0);
        dsymv_l(nv, 0.5, qcqp->Hq + ii, 0, 0, qcqp_sol->v, 0, 1.0, ws->tmp_nv, 0, ws->tmp_nv, 0);
        tmp = ddot(nv, ws->tmp_nv, 0, qcqp_sol->v, 0);
        tmp = -d[2 * nb + 2 * ng + nq + ii] - tmp;
        t[2 * nb + 2 * ng + nq + ii] = thr0 > tmp ? thr0 : tmp;
        lam[2 * nb + 2 * ng + nq + ii] = mu0 / t[2 * nb + 2 * ng + nq + ii];
#endif
    }

    // TODO rewrite all the above taking some pointers to key parts, e.g. lam_lb, lam_ub, and make relative to them
}


void d_dense_qcqp_approx_qp(struct d_dense_qcqp* qcqp, struct d_dense_qcqp_sol* qcqp_sol, struct d_dense_qp* qp, struct d_dense_qcqp_ipm_ws* ws) {

    int nv = qcqp->dim->nv;
    int ne = qcqp->dim->ne;
    int nb = qcqp->dim->nb;
    int ng = qcqp->dim->ng;
    int nq = qcqp->dim->nq;
    int ns = qcqp->dim->ns;

    double tmp;

    int ii;


    dveccp(2 * nb + 2 * ng + 2 * nq + 2 * ns, qcqp->d, 0, qp->d, 0);

    dgecp(nv, nv, qcqp->Hv, 0, 0, qp->Hv, 0, 0);

    dvecse(nv, 0.0, ws->qcqp_res_ws->q_adj, 0);

    for (ii = 0; ii < nq; ii++) {
        tmp = -VECEL(qcqp_sol->lam, nb + ng + ii) + VECEL(qcqp_sol->lam, 2 * nb + 2 * ng + nq + ii);
        dgead(nv, nv, tmp, qcqp->Hq + ii, 0, 0, qp->Hv, 0, 0);

        dsymv_l(nv, 1.0, qcqp->Hq + ii, 0, 0, qcqp_sol->v, 0, 0.0, ws->tmp_nv + 0, 0, ws->tmp_nv + 0, 0);
        dcolex(nv, qcqp->Ct, 0, ng + ii, ws->tmp_nv + 1, 0);
        daxpy(nv, 1.0, ws->tmp_nv + 0, 0, ws->tmp_nv + 1, 0, ws->tmp_nv + 1, 0);
        dcolin(nv, ws->tmp_nv + 1, 0, qp->Ct, 0, ng + ii);
        daxpy(nv, tmp, ws->tmp_nv + 1, 0, ws->qcqp_res_ws->q_adj, 0, ws->qcqp_res_ws->q_adj, 0);

        dcolex(nv, qcqp->Ct, 0, ng + ii, ws->tmp_nv + 1, 0);
        daxpy(nv, 0.5, ws->tmp_nv + 0, 0, ws->tmp_nv + 1, 0, ws->tmp_nv + 1, 0);
        tmp = ddot(nv, ws->tmp_nv + 1, 0, qcqp_sol->v, 0);
        // TODO maybe swap signs?
        VECEL(qp->d, nb + ng + ii) += -tmp;
        VECEL(qp->d, 2 * nb + 2 * ng + nq + ii) += +tmp;
        VECEL(ws->qcqp_res_ws->q_fun, ii) = tmp;
    }

    dveccp(2 * nb + 2 * ng + 2 * nq + 2 * ns, qcqp->d_mask, 0, qp->d_mask, 0);

    dgecp(ne, nv, qcqp->A, 0, 0, qp->A, 0, 0);

    dgecp(nv, ng, qcqp->Ct, 0, 0, qp->Ct, 0, 0);

    dveccp(nv + 2 * ns, qcqp->gz, 0, qp->gz, 0);

    dveccp(ne, qcqp->b, 0, qp->b, 0);

    dveccp(2 * nb + 2 * ng + 2 * nq + 2 * ns, qcqp->m, 0, qp->m, 0);

    dveccp(2 * ns, qcqp->Z, 0, qp->Z, 0);

    for (ii = 0; ii < nb; ii++)
        qp->idxb[ii] = qcqp->idxb[ii];

    for (ii = 0; ii < nb + ng + nq; ii++)
        qp->idxs_rev[ii] = qcqp->idxs_rev[ii];
}


void d_dense_qcqp_update_qp(struct d_dense_qcqp* qcqp, struct d_dense_qcqp_sol* qcqp_sol, struct d_dense_qp* qp, struct d_dense_qcqp_ipm_ws* ws) {

    int nv = qcqp->dim->nv;
    int ne = qcqp->dim->ne;
    int nb = qcqp->dim->nb;
    int ng = qcqp->dim->ng;
    int nq = qcqp->dim->nq;
    int ns = qcqp->dim->ns;

    double tmp;

    int ii;

    // TODO only the 2*nq part needed !!!!!
    dveccp(2 * nb + 2 * ng + 2 * nq + 2 * ns, qcqp->d, 0, qp->d, 0);

    dgecp(nv, nv, qcqp->Hv, 0, 0, qp->Hv, 0, 0);

    dvecse(nv, 0.0, ws->qcqp_res_ws->q_adj, 0);

    for (ii = 0; ii < nq; ii++) {
        tmp = -VECEL(qcqp_sol->lam, nb + ng + ii) + VECEL(qcqp_sol->lam, 2 * nb + 2 * ng + nq + ii);
        dgead(nv, nv, tmp, qcqp->Hq + ii, 0, 0, qp->Hv, 0, 0);

        dsymv_l(nv, 1.0, qcqp->Hq + ii, 0, 0, qcqp_sol->v, 0, 0.0, ws->tmp_nv + 0, 0, ws->tmp_nv + 0, 0);
        dcolex(nv, qcqp->Ct, 0, ng + ii, ws->tmp_nv + 1, 0);
        daxpy(nv, 1.0, ws->tmp_nv + 0, 0, ws->tmp_nv + 1, 0, ws->tmp_nv + 1, 0);
        dcolin(nv, ws->tmp_nv + 1, 0, qp->Ct, 0, ng + ii);
        daxpy(nv, tmp, ws->tmp_nv + 1, 0, ws->qcqp_res_ws->q_adj, 0, ws->qcqp_res_ws->q_adj, 0);

        dcolex(nv, qcqp->Ct, 0, ng + ii, ws->tmp_nv + 1, 0);
        daxpy(nv, 0.5, ws->tmp_nv + 0, 0, ws->tmp_nv + 1, 0, ws->tmp_nv + 1, 0);
        tmp = ddot(nv, ws->tmp_nv + 1, 0, qcqp_sol->v, 0);
        // TODO maybe swap signs?
        VECEL(qp->d, nb + ng + ii) += -tmp;
        VECEL(qp->d, 2 * nb + 2 * ng + nq + ii) += +tmp;
        VECEL(ws->qcqp_res_ws->q_fun, ii) = tmp;
    }

    // TODO needed ?????
    dveccp(2 * nb + 2 * ng + 2 * nq + 2 * ns, qcqp->m, 0, qp->m, 0);
}


void d_dense_qcqp_update_qp_abs_step(struct d_dense_qcqp* qcqp, struct d_dense_qcqp_sol* qcqp_sol, struct d_dense_qp* qp, struct d_dense_qcqp_ipm_ws* ws) {

    int nv = qcqp->dim->nv;
    int ne = qcqp->dim->ne;
    int nb = qcqp->dim->nb;
    int ng = qcqp->dim->ng;
    int nq = qcqp->dim->nq;
    int ns = qcqp->dim->ns;

    double tmp;

    int ii;

    dveccp(2 * nb + 2 * ng + 2 * nq + 2 * ns, qcqp->d, 0, qp->d, 0);

    dgecp(nv, nv, qcqp->Hv, 0, 0, qp->Hv, 0, 0);

    dvecse(nv, 0.0, ws->qcqp_res_ws->q_adj, 0);

    for (ii = 0; ii < nq; ii++) {
        tmp = -VECEL(qcqp_sol->lam, nb + ng + ii) + VECEL(qcqp_sol->lam, 2 * nb + 2 * ng + nq + ii);
        dgead(nv, nv, tmp, qcqp->Hq + ii, 0, 0, qp->Hv, 0, 0);

        dsymv_l(nv, 1.0, qcqp->Hq + ii, 0, 0, qcqp_sol->v, 0, 0.0, ws->tmp_nv + 0, 0, ws->tmp_nv + 0, 0);
        dcolex(nv, qcqp->Ct, 0, ng + ii, ws->tmp_nv + 1, 0);
        daxpy(nv, 1.0, ws->tmp_nv + 0, 0, ws->tmp_nv + 1, 0, ws->tmp_nv + 1, 0);
        dcolin(nv, ws->tmp_nv + 1, 0, qp->Ct, 0, ng + ii);
        daxpy(nv, tmp, ws->tmp_nv + 1, 0, ws->qcqp_res_ws->q_adj, 0, ws->qcqp_res_ws->q_adj, 0);

        ////		daxpy(nv, 0.5, ws->tmp_nv+0, 0, qcqp->gq+ii, 0, ws->tmp_nv+1, 0);
        //		daxpy(nv, 0.5, ws->tmp_nv+0, 0, qcqp->gq+ii, 0, ws->tmp_nv+0, 0);
        //		daxpy(nv, -1.0, ws->tmp_nv+1, 0, ws->tmp_nv+0, 0, ws->tmp_nv+1, 0);
        daxpby(nv, -1.0, ws->tmp_nv + 1, 0, 0.5, ws->tmp_nv + 0, 0, ws->tmp_nv + 1, 0);
        dcolex(nv, qcqp->Ct, 0, ng + ii, ws->tmp_nv + 0, 0);
        daxpy(nv, 1.0, ws->tmp_nv + 0, 0, ws->tmp_nv + 1, 0, ws->tmp_nv + 1, 0);
        tmp = ddot(nv, ws->tmp_nv + 1, 0, qcqp_sol->v, 0);
        // TODO maybe swap signs?
        VECEL(qp->d, nb + ng + ii) += -tmp;
        VECEL(qp->d, 2 * nb + 2 * ng + nq + ii) += +tmp;
        VECEL(ws->qcqp_res_ws->q_fun, ii) = tmp;
    }

    dveccp(2 * nb + 2 * ng + 2 * nq + 2 * ns, qcqp->m, 0, qp->m, 0);
}


void d_dense_qcqp_sol_conv_qp_sol(struct d_dense_qcqp_sol* qcqp_sol, struct d_dense_qp_sol* qp_sol) {

    int nv = qcqp_sol->dim->nv;
    int ne = qcqp_sol->dim->ne;
    int nb = qcqp_sol->dim->nb;
    int ng = qcqp_sol->dim->ng;
    int nq = qcqp_sol->dim->nq;
    int ns = qcqp_sol->dim->ns;

    dveccp(nv + 2 * ns, qcqp_sol->v, 0, qp_sol->v, 0);

    dveccp(ne, qcqp_sol->pi, 0, qp_sol->pi, 0);

    dveccp(2 * nb + 2 * ng + 2 * nq + 2 * ns, qcqp_sol->lam, 0, qp_sol->lam, 0);

    dveccp(2 * nb + 2 * ng + 2 * nq + 2 * ns, qcqp_sol->t, 0, qp_sol->t, 0);
}


void d_dense_qp_sol_conv_qcqp_sol(struct d_dense_qp_sol* qp_sol, struct d_dense_qcqp_sol* qcqp_sol) {

    int nv = qcqp_sol->dim->nv;
    int ne = qcqp_sol->dim->ne;
    int nb = qcqp_sol->dim->nb;
    int ng = qcqp_sol->dim->ng;
    int nq = qcqp_sol->dim->nq;
    int ns = qcqp_sol->dim->ns;

    dveccp(nv + 2 * ns, qp_sol->v, 0, qcqp_sol->v, 0);

    dveccp(ne, qp_sol->pi, 0, qcqp_sol->pi, 0);

    dveccp(2 * nb + 2 * ng + 2 * nq + 2 * ns, qp_sol->lam, 0, qcqp_sol->lam, 0);

    dveccp(2 * nb + 2 * ng + 2 * nq + 2 * ns, qp_sol->t, 0, qcqp_sol->t, 0);
}


void d_dense_qcqp_res_conv_qp_res(struct d_dense_qcqp_res* qcqp_res, struct d_dense_qp_res* qp_res) {

    int nv = qcqp_res->dim->nv;
    int ne = qcqp_res->dim->ne;
    int nb = qcqp_res->dim->nb;
    int ng = qcqp_res->dim->ng;
    int nq = qcqp_res->dim->nq;
    int ns = qcqp_res->dim->ns;

    dveccp(nv + 2 * ns, qcqp_res->res_g, 0, qp_res->res_g, 0);

    dveccp(ne, qcqp_res->res_b, 0, qp_res->res_b, 0);

    dveccp(2 * nb + 2 * ng + 2 * nq + 2 * ns, qcqp_res->res_d, 0, qp_res->res_d, 0);

    dveccp(2 * nb + 2 * ng + 2 * nq + 2 * ns, qcqp_res->res_m, 0, qp_res->res_m, 0);

    qp_res->res_mu = qcqp_res->res_mu;
}


void d_dense_qcqp_ipm_solve(struct d_dense_qcqp* qcqp, struct d_dense_qcqp_sol* qcqp_sol, struct d_dense_qcqp_ipm_arg* qcqp_arg, struct d_dense_qcqp_ipm_ws* qcqp_ws) {

    int nv = qcqp->dim->nv;
    int ne = qcqp->dim->ne;
    int nb = qcqp->dim->nb;
    int ng = qcqp->dim->ng;
    int nq = qcqp->dim->nq;
    int ns = qcqp->dim->ns;

    // extract stuff

    struct d_dense_qp* qp = qcqp_ws->qp;
    struct d_dense_qp_sol* qp_sol = qcqp_ws->qp_sol;
    struct d_dense_qp_ipm_ws* qp_ws = qcqp_ws->qp_ws;
    struct d_dense_qp_ipm_arg* qp_arg = qcqp_arg->qp_arg;

    struct d_dense_qcqp_dim* qcqp_dim = qcqp->dim;
    struct d_dense_qcqp_res* qcqp_res = qcqp_ws->qcqp_res;
    struct d_dense_qcqp_res_ws* qcqp_res_ws = qcqp_ws->qcqp_res_ws;

    struct d_core_qp_ipm_workspace* cws = qp_ws->core_workspace;

    int kk, ii, idx;
    double mu;

    double* stat = qp_ws->stat;
    int stat_m = qp_ws->stat_m;
    int stat_max = qp_ws->stat_max;

    //	int qcqp_nv = qcqp->dim->nv + 2*qcqp->dim->ns;
    //	int qcqp_ne = qcqp->dim->ne;
    //	int qcqp_nc = 2*qcqp->dim->nb + 2*qcqp->dim->ng + qcqp->dim->nq + 2*qcqp->dim->ns;


    // qp_arg to core workspace
    cws->lam_min = qp_arg->lam_min;
    cws->t_min = qp_arg->t_min;
    cws->t_min_inv = qp_arg->t_min > 0.0 ? 1.0 / qp_arg->t_min : 1e30;
    cws->split_step = qp_arg->split_step;
    cws->t_lam_min = qp_arg->t_lam_min;

    // alias qp vectors into qp_sol
    cws->v = qp_sol->v->pa;
    cws->pi = qp_sol->pi->pa;
    cws->lam = qp_sol->lam->pa;
    cws->t = qp_sol->t->pa;

    // alias members of qp_step
    qp_ws->qp_step->dim = qp->dim;
    qp_ws->qp_step->Hv = qp->Hv;
    qp_ws->qp_step->A = qp->A;
    qp_ws->qp_step->Ct = qp->Ct;
    qp_ws->qp_step->Z = qp->Z;
    qp_ws->qp_step->idxb = qp->idxb;
    qp_ws->qp_step->idxs_rev = qp->idxs_rev;
    qp_ws->qp_step->gz = qp_ws->res->res_g;
    qp_ws->qp_step->b = qp_ws->res->res_b;
    qp_ws->qp_step->d = qp_ws->res->res_d;
    qp_ws->qp_step->m = qp_ws->res->res_m;
    qp_ws->qp_step->d_mask = qp->d_mask;

    // alias members of qp_itref
    qp_ws->qp_itref->dim = qp->dim;
    qp_ws->qp_itref->Hv = qp->Hv;
    qp_ws->qp_itref->A = qp->A;
    qp_ws->qp_itref->Ct = qp->Ct;
    qp_ws->qp_itref->Z = qp->Z;
    qp_ws->qp_itref->idxb = qp->idxb;
    qp_ws->qp_itref->idxs_rev = qp->idxs_rev;
    qp_ws->qp_itref->gz = qp_ws->res_itref->res_g;
    qp_ws->qp_itref->b = qp_ws->res_itref->res_b;
    qp_ws->qp_itref->d = qp_ws->res_itref->res_d;
    qp_ws->qp_itref->m = qp_ws->res_itref->res_m;
    qp_ws->qp_itref->d_mask = qp->d_mask;

    double* qcqp_res_max = qcqp_res->res_max;

    qp_ws->use_A_fact = 0;

    // cache q_fun & q_adj from approx/update for res
    qcqp_ws->qcqp_res_ws->use_q_fun = 1;
    qcqp_ws->qcqp_res_ws->use_q_adj = 1;


    // disregard soft constr on (disregarded) lower quard constr
    dvecse(nq, 0.0, qcqp->d_mask, nb + ng);  // TODO needed ???
    // TODO probably remove when using only idxs_rev, as the same slack may be associated with other constraints !!!!!
    for (ii = 0; ii < nq; ii++) {
        idx = qcqp->idxs_rev[nb + ng + ii];
        if (idx != -1) {
            VECEL(qcqp->d_mask, 2 * nb + 2 * ng + 2 * nq + idx) = 0.0;
        }
    }


    // initialize qcqp & qp
    d_dense_qcqp_init_var(qcqp, qcqp_sol, qcqp_arg, qcqp_ws);
    d_dense_qcqp_sol_conv_qp_sol(qcqp_sol, qp_sol);

    // approximate qcqp with a qp
    d_dense_qcqp_approx_qp(qcqp, qcqp_sol, qp, qcqp_ws);


    // detect constr mask
    int mask_unconstr;
    int nc_mask = 0;
    for (ii = 0; ii < cws->nc; ii++) {
        if (qp->d_mask->pa[ii] != 0.0)
            nc_mask++;
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
    // always mask lower quadratic constr
    qp_ws->mask_constr = 1;


    // no constraints
    if (cws->nc == 0 | mask_unconstr == 1) {
        d_fact_solve_kkt_unconstr_dense_qp(qp, qp_sol, qp_arg, qp_ws);
        d_dense_qp_sol_conv_qcqp_sol(qp_sol, qcqp_sol);
        if (qp_arg->comp_res_exit) {
            // compute residuals
            d_dense_qcqp_res_compute(qcqp, qcqp_sol, qcqp_res, qcqp_res_ws);
            // XXX no constraints, so no mask
            d_dense_qcqp_res_compute_inf_norm(qcqp_res);
            // save infinity norm of residuals
            if (0 < stat_max) {
                stat[6] = qcqp_res_max[0];
                stat[7] = qcqp_res_max[1];
                stat[8] = qcqp_res_max[2];
                stat[9] = qcqp_res_max[3];
                stat[10] = qcqp_res->obj;
            }
            cws->mu = qcqp_res->res_mu;
        }
        // save info before return
        qcqp_ws->iter = 0;

        if (isnan(VECEL(qcqp_sol->v, 0))) {
            // NaN in the solution
            qcqp_ws->status = NAN_SOL;
        }
        // #else
        //         if (VECEL(qcqp_sol->v, 0) != VECEL(qcqp_sol->v, 0)) {
        //             // NaN in the solution
        //             qcqp_ws->status = NAN_sol;
        //         }

        else {
            // normal return
            qcqp_ws->status = SUCCESS;
        }
    }


    cws->alpha = 1.0;


    // absolute IPM formulation
    if (qp_arg->abs_form) {

        // alias members of qp_step
        qp_ws->qp_step->dim = qp->dim;
        qp_ws->qp_step->Hv = qp->Hv;
        qp_ws->qp_step->A = qp->A;
        qp_ws->qp_step->Ct = qp->Ct;
        qp_ws->qp_step->Z = qp->Z;
        qp_ws->qp_step->idxb = qp->idxb;
        qp_ws->qp_step->idxs_rev = qp->idxs_rev;
        qp_ws->qp_step->gz = qp->gz;
        qp_ws->qp_step->b = qp->b;
        qp_ws->qp_step->d = qp->d;
        qp_ws->qp_step->m = qp_ws->tmp_m;

        // alias core workspace
        cws->res_m = qp_ws->qp_step->m->pa;

        // update approximation of qcqp as qp for absolute step
        d_dense_qcqp_update_qp_abs_step(qcqp, qcqp_sol, qp, qcqp_ws);

        // compute mu
        mu = dvecmuldot(cws->nc, qp_sol->lam, 0, qp_sol->t, 0, qp_ws->tmp_m, 0);
        mu /= cws->nc;
        cws->mu = mu;

        // IPM loop (absolute formulation)
        for (kk = 0;
             kk<qcqp_arg->iter_max &
                cws->alpha>
                     qcqp_arg->alpha_min &
             mu > qcqp_arg->res_m_max;
             kk++) {

            // compute delta step
            d_dense_qp_ipm_abs_step(kk, qp, qp_sol, qp_arg, qp_ws);
            // blasfeo_print_exp_tran_dvec(cws->nc, qp_sol->lam, 0);
            d_dense_qp_sol_conv_qcqp_sol(qp_sol, qcqp_sol);

            // update approximation of qcqp as qp for absolute step
            d_dense_qcqp_update_qp_abs_step(qcqp, qcqp_sol, qp, qcqp_ws);

            // compute mu
            mu = dvecmuldot(cws->nc, qcqp_sol->lam, 0, qcqp_sol->t, 0, qp_ws->tmp_m, 0);
            mu /= cws->nc;
            cws->mu = mu;
            if (kk + 1 < stat_max)
                stat[stat_m * (kk + 1) + 5] = mu;
        }

        if (qp_arg->comp_res_exit) {
            // compute residuals
            d_dense_qcqp_res_compute(qcqp, qcqp_sol, qcqp_res, qcqp_res_ws);
            if (qp_ws->mask_constr) {
                // mask out disregarded constraints
                dvecmul(2 * ns, qp->d_mask, 2 * nb + 2 * ng + 2 * nq, qcqp_res->res_g, nv, qcqp_res->res_g, nv);
                dvecmul(cws->nc, qp->d_mask, 0, qcqp_res->res_d, 0, qcqp_res->res_d, 0);
                dvecmul(cws->nc, qp->d_mask, 0, qcqp_res->res_m, 0, qcqp_res->res_m, 0);
            }
            d_dense_qcqp_res_compute_inf_norm(qcqp_res);
            // save infinity norm of residuals
            // XXX it is already kk+1
            if (kk < stat_max) {
                stat[stat_m * (kk + 0) + 6] = qcqp_res_max[0];
                stat[stat_m * (kk + 0) + 7] = qcqp_res_max[1];
                stat[stat_m * (kk + 0) + 8] = qcqp_res_max[2];
                stat[stat_m * (kk + 0) + 9] = qcqp_res_max[3];
                stat[stat_m * (kk + 0) + 10] = qcqp_res->obj;
            }
        }

        goto set_status;
    }


    // compute residuals
    d_dense_qcqp_res_compute(qcqp, qcqp_sol, qcqp_res, qcqp_res_ws);
    if (qp_ws->mask_constr) {
        // mask out disregarded constraints
        dvecmul(2 * ns, qp->d_mask, 2 * nb + 2 * ng + 2 * nq, qcqp_res->res_g, nv, qcqp_res->res_g, nv);
        dvecmul(cws->nc, qp->d_mask, 0, qcqp_res->res_d, 0, qcqp_res->res_d, 0);
        dvecmul(cws->nc, qp->d_mask, 0, qcqp_res->res_m, 0, qcqp_res->res_m, 0);
    }
    d_dense_qcqp_res_compute_inf_norm(qcqp_res);
    d_dense_qcqp_res_conv_qp_res(qcqp_res, qp_ws->res);
    cws->mu = qcqp_res->res_mu;
    // save infinity norm of residuals
    if (0 < stat_max) {
        stat[stat_m * (0) + 6] = qcqp_res_max[0];
        stat[stat_m * (0) + 7] = qcqp_res_max[1];
        stat[stat_m * (0) + 8] = qcqp_res_max[2];
        stat[stat_m * (0) + 9] = qcqp_res_max[3];
        stat[stat_m * (0) + 10] = qcqp_res->obj;
    }


    // relative (delta) IPM formulation
    for (kk = 0;
         kk<qcqp_arg->iter_max &
            cws->alpha>
                 qcqp_arg->alpha_min &
         (qcqp_res_max[0] > qcqp_arg->res_g_max |
          qcqp_res_max[1] > qcqp_arg->res_b_max |
          qcqp_res_max[2] > qcqp_arg->res_d_max |
          qcqp_res_max[3] > qcqp_arg->res_m_max);
         kk++) {

        // hessian is updated with quad constr: can not reuse hessian factorization !!!
        qp_ws->use_hess_fact = 0;

        // compute delta step
        d_dense_qp_ipm_delta_step(kk, qp, qp_sol, qp_arg, qp_ws);
        // blasfeo_print_exp_tran_dvec(cws->nc, qp_sol->lam, 0);
        d_dense_qp_sol_conv_qcqp_sol(qp_sol, qcqp_sol);
        // XXX maybe not needed
        if (qp_ws->mask_constr) {
            // mask out disregarded constraints
            dvecmul(cws->nc, qp->d_mask, 0, qcqp_sol->lam, 0, qcqp_sol->lam, 0);
        }

        // update approximation of qcqp as qp
        d_dense_qcqp_update_qp(qcqp, qcqp_sol, qp, qcqp_ws);

        // compute residuals
        d_dense_qcqp_res_compute(qcqp, qcqp_sol, qcqp_res, qcqp_res_ws);
        if (qp_ws->mask_constr) {
            // mask out disregarded constraints
            dvecmul(2 * ns, qp->d_mask, 2 * nb + 2 * ng + 2 * nq, qcqp_res->res_g, nv, qcqp_res->res_g, nv);
            dvecmul(cws->nc, qp->d_mask, 0, qcqp_res->res_d, 0, qcqp_res->res_d, 0);
            dvecmul(cws->nc, qp->d_mask, 0, qcqp_res->res_m, 0, qcqp_res->res_m, 0);
        }
        d_dense_qcqp_res_compute_inf_norm(qcqp_res);
        d_dense_qcqp_res_conv_qp_res(qcqp_res, qp_ws->res);
        cws->mu = qcqp_res->res_mu;
        // save infinity norm of residuals
        if (kk + 1 < stat_max) {
            stat[stat_m * (kk + 1) + 5] = qcqp_res->res_mu;
            stat[stat_m * (kk + 1) + 6] = qcqp_res_max[0];
            stat[stat_m * (kk + 1) + 7] = qcqp_res_max[1];
            stat[stat_m * (kk + 1) + 8] = qcqp_res_max[2];
            stat[stat_m * (kk + 1) + 9] = qcqp_res_max[3];
            stat[stat_m * (kk + 1) + 10] = qcqp_res->obj;
        }
    }

set_status:

    // save info before return
    qcqp_ws->iter = kk;

    if (kk == qcqp_arg->iter_max) {
        // max iteration number reached
        qcqp_ws->status = MAX_ITER;
    } else if (cws->alpha <= qcqp_arg->alpha_min) {
        // min step lenght
        qcqp_ws->status = MIN_STEP;
    }

    else if (isnan(cws->mu)) {
        // NaN in the solution
        qcqp_ws->status = NAN_SOL;
    }
    // #else
    //     else if (cws->mu != cws->mu) {
    //         // NaN in the solution
    //         qcqp_ws->status = NAN_sol;
    //     }

    else {
        // normal return
        qcqp_ws->status = SUCCESS;
    }
}
