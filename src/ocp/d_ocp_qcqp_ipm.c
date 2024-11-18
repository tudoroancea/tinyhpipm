#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "tinyhpipm/blas.h"
#include "tinyhpipm/common.h"
#include "tinyhpipm/ipm_core/d_core_qp_ipm.h"
#include "tinyhpipm/ipm_core/d_core_qp_ipm_aux.h"
#include "tinyhpipm/ocp/d_ocp_qcqp.h"
#include "tinyhpipm/ocp/d_ocp_qcqp_dim.h"
#include "tinyhpipm/ocp/d_ocp_qcqp_ipm.h"
#include "tinyhpipm/ocp/d_ocp_qcqp_res.h"
#include "tinyhpipm/ocp/d_ocp_qcqp_sol.h"
#include "tinyhpipm/ocp/d_ocp_qp.h"
#include "tinyhpipm/ocp/d_ocp_qp_dim.h"
#include "tinyhpipm/ocp/d_ocp_qp_ipm.h"
#include "tinyhpipm/ocp/d_ocp_qp_kkt.h"
#include "tinyhpipm/ocp/d_ocp_qp_res.h"
#include "tinyhpipm/ocp/d_ocp_qp_sol.h"
#include "tinyhpipm/ocp/d_ocp_qp_utils.h"
// #include <hpipm_d_ocp_qcqp_utils.h>


hpipm_size_t d_ocp_qcqp_ipm_arg_strsize() {
    return sizeof(struct d_ocp_qcqp_ipm_arg);
}


hpipm_size_t d_ocp_qcqp_ipm_arg_memsize(struct d_ocp_qcqp_dim* dim) {

    hpipm_size_t size = 0;

    size += 1 * sizeof(struct d_ocp_qp_ipm_arg);
    size += 1 * d_ocp_qp_ipm_arg_memsize(dim->qp_dim);

    size = (size + 63) / 64 * 64;  // make multiple of typical cache line size
    size += 1 * 64;  // align once to typical cache line size

    return size;
}


void d_ocp_qcqp_ipm_arg_create(struct d_ocp_qcqp_dim* dim, struct d_ocp_qcqp_ipm_arg* arg, void* mem) {

    // loop index
    int ii;

    // zero memory (to avoid corrupted memory like e.g. NaN)
    hpipm_size_t memsize = d_ocp_qcqp_ipm_arg_memsize(dim);
    hpipm_zero_memset(memsize, mem);

    // qp_dim struct
    struct d_ocp_qp_ipm_arg* arg_ptr = mem;

    arg->qp_arg = arg_ptr;
    arg_ptr += 1;

    // align to typical cache line size
    hpipm_size_t s_ptr = (hpipm_size_t) arg_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;

    // void
    char* c_ptr = (char*) s_ptr;

    d_ocp_qp_ipm_arg_create(dim->qp_dim, arg->qp_arg, c_ptr);
    c_ptr += arg->qp_arg->memsize;


    arg->memsize = memsize;  // d_ocp_qcqp_ipm_arg_memsize(dim);

#if defined(RUNTIME_CHECKS)
    if (c_ptr > ((char*) mem) + arg->memsize) {
        printf("\nerror: d_ocp_qcqp_ipm_arg_create: outside memory bounds!\n\n");
        exit(1);
    }
#endif
}


void d_ocp_qcqp_ipm_arg_set_default(enum tinyhpipm_mode mode, struct d_ocp_qcqp_ipm_arg* arg) {

    d_ocp_qp_ipm_arg_set_default(mode, arg->qp_arg);

    double mu0, alpha_min, res_g_max, res_b_max, res_d_max, res_m_max, reg_prim, lam_min, t_min;
    int iter_max, stat_max, pred_corr, cond_pred_corr, itref_pred_max, itref_corr_max, lq_fact, warm_start, abs_form, comp_res_exit, comp_res_pred, square_root_alg, comp_dual_sol_eq, split_step, t_lam_min;

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
        warm_start = 0;
        abs_form = 1;
        comp_dual_sol_eq = 0;
        comp_res_exit = 0;
        comp_res_pred = 0;
        split_step = 1;
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
        warm_start = 0;
        abs_form = 0;
        comp_dual_sol_eq = 1;
        comp_res_exit = 1;
        comp_res_pred = 1;
        split_step = 1;
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
        warm_start = 0;
        abs_form = 0;
        comp_dual_sol_eq = 1;
        comp_res_exit = 1;
        comp_res_pred = 1;
        split_step = 0;
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
        warm_start = 0;
        abs_form = 0;
        comp_dual_sol_eq = 1;
        comp_res_exit = 1;
        comp_res_pred = 1;
        split_step = 0;
        t_lam_min = 2;
    } else {
        printf("\nerror: d_ocp_qcqp_ipm_arg_set_default: wrong set default mode\n");
        exit(1);
    }

    // use individual setters when available
    d_ocp_qcqp_ipm_arg_set_mu0(&mu0, arg);
    d_ocp_qcqp_ipm_arg_set_alpha_min(&alpha_min, arg);
    d_ocp_qcqp_ipm_arg_set_tol_stat(&res_g_max, arg);
    d_ocp_qcqp_ipm_arg_set_tol_eq(&res_b_max, arg);
    d_ocp_qcqp_ipm_arg_set_tol_ineq(&res_d_max, arg);
    d_ocp_qcqp_ipm_arg_set_tol_comp(&res_m_max, arg);
    d_ocp_qcqp_ipm_arg_set_iter_max(&iter_max, arg);
    arg->stat_max = stat_max;
    d_ocp_qcqp_ipm_arg_set_pred_corr(&pred_corr, arg);
    d_ocp_qcqp_ipm_arg_set_cond_pred_corr(&cond_pred_corr, arg);
    d_ocp_qcqp_ipm_arg_set_ric_alg(&square_root_alg, arg);
    arg->itref_pred_max = itref_pred_max;
    arg->itref_corr_max = itref_corr_max;
    d_ocp_qcqp_ipm_arg_set_reg_prim(&reg_prim, arg);
    arg->lq_fact = lq_fact;
    d_ocp_qcqp_ipm_arg_set_lam_min(&lam_min, arg);
    d_ocp_qcqp_ipm_arg_set_t_min(&t_min, arg);
    d_ocp_qcqp_ipm_arg_set_warm_start(&warm_start, arg);
    arg->abs_form = abs_form;
    d_ocp_qcqp_ipm_arg_set_comp_res_pred(&comp_res_pred, arg);
    d_ocp_qcqp_ipm_arg_set_comp_res_exit(&comp_res_pred, arg);
    d_ocp_qcqp_ipm_arg_set_split_step(&split_step, arg);
    d_ocp_qcqp_ipm_arg_set_t_lam_min(&t_lam_min, arg);
    arg->mode = mode;
}


void d_ocp_qcqp_ipm_arg_set(char* field, void* value, struct d_ocp_qcqp_ipm_arg* arg) {
    if (hpipm_strcmp(field, "iter_max")) {
        d_ocp_qcqp_ipm_arg_set_iter_max(value, arg);
    } else if (hpipm_strcmp(field, "alpha_min")) {
        d_ocp_qcqp_ipm_arg_set_alpha_min(value, arg);
    } else if (hpipm_strcmp(field, "mu0")) {
        d_ocp_qcqp_ipm_arg_set_mu0(value, arg);
    } else if (hpipm_strcmp(field, "tol_stat")) {
        d_ocp_qcqp_ipm_arg_set_tol_stat(value, arg);
    } else if (hpipm_strcmp(field, "tol_eq")) {
        d_ocp_qcqp_ipm_arg_set_tol_eq(value, arg);
    } else if (hpipm_strcmp(field, "tol_ineq")) {
        d_ocp_qcqp_ipm_arg_set_tol_ineq(value, arg);
    } else if (hpipm_strcmp(field, "tol_comp")) {
        d_ocp_qcqp_ipm_arg_set_tol_comp(value, arg);
    } else if (hpipm_strcmp(field, "reg_prim")) {
        d_ocp_qcqp_ipm_arg_set_reg_prim(value, arg);
    } else if (hpipm_strcmp(field, "warm_start")) {
        d_ocp_qcqp_ipm_arg_set_warm_start(value, arg);
    } else if (hpipm_strcmp(field, "pred_corr")) {
        d_ocp_qcqp_ipm_arg_set_pred_corr(value, arg);
    } else if (hpipm_strcmp(field, "cond_pred_corr")) {
        d_ocp_qcqp_ipm_arg_set_cond_pred_corr(value, arg);
    } else if (hpipm_strcmp(field, "ric_alg")) {
        d_ocp_qcqp_ipm_arg_set_ric_alg(value, arg);
    } else if (hpipm_strcmp(field, "comp_res_exit")) {
        d_ocp_qcqp_ipm_arg_set_comp_res_exit(value, arg);
    } else if (hpipm_strcmp(field, "comp_res_pred")) {
        d_ocp_qcqp_ipm_arg_set_comp_res_pred(value, arg);
    } else if (hpipm_strcmp(field, "lam_min")) {
        d_ocp_qcqp_ipm_arg_set_lam_min(value, arg);
    } else if (hpipm_strcmp(field, "t_min")) {
        d_ocp_qcqp_ipm_arg_set_t_min(value, arg);
    } else if (hpipm_strcmp(field, "split_step")) {
        d_ocp_qcqp_ipm_arg_set_split_step(value, arg);
    } else if (hpipm_strcmp(field, "t_lam_min")) {
        d_ocp_qcqp_ipm_arg_set_t_lam_min(value, arg);
    } else {
        printf("error: d_ocp_qcqp_ipm_arg_set: wrong field %s\n", field);
        exit(1);
    }
}


void d_ocp_qcqp_ipm_arg_set_iter_max(int* value, struct d_ocp_qcqp_ipm_arg* arg) {
    arg->iter_max = *value;
    d_ocp_qp_ipm_arg_set_iter_max(value, arg->qp_arg);
}


void d_ocp_qcqp_ipm_arg_set_alpha_min(double* value, struct d_ocp_qcqp_ipm_arg* arg) {
    arg->alpha_min = *value;
    d_ocp_qp_ipm_arg_set_alpha_min(value, arg->qp_arg);
}


void d_ocp_qcqp_ipm_arg_set_mu0(double* value, struct d_ocp_qcqp_ipm_arg* arg) {
    arg->mu0 = *value;
    d_ocp_qp_ipm_arg_set_mu0(value, arg->qp_arg);
}


void d_ocp_qcqp_ipm_arg_set_tol_stat(double* value, struct d_ocp_qcqp_ipm_arg* arg) {
    arg->res_g_max = *value;
    d_ocp_qp_ipm_arg_set_tol_stat(value, arg->qp_arg);
}


void d_ocp_qcqp_ipm_arg_set_tol_eq(double* value, struct d_ocp_qcqp_ipm_arg* arg) {
    arg->res_b_max = *value;
    d_ocp_qp_ipm_arg_set_tol_stat(value, arg->qp_arg);
}


void d_ocp_qcqp_ipm_arg_set_tol_ineq(double* value, struct d_ocp_qcqp_ipm_arg* arg) {
    arg->res_d_max = *value;
    d_ocp_qp_ipm_arg_set_tol_stat(value, arg->qp_arg);
}


void d_ocp_qcqp_ipm_arg_set_tol_comp(double* value, struct d_ocp_qcqp_ipm_arg* arg) {
    arg->res_m_max = *value;
    d_ocp_qp_ipm_arg_set_tol_comp(value, arg->qp_arg);
}


void d_ocp_qcqp_ipm_arg_set_reg_prim(double* value, struct d_ocp_qcqp_ipm_arg* arg) {
    arg->reg_prim = *value;
    d_ocp_qp_ipm_arg_set_reg_prim(value, arg->qp_arg);
}


void d_ocp_qcqp_ipm_arg_set_warm_start(int* value, struct d_ocp_qcqp_ipm_arg* arg) {
    arg->warm_start = *value;
    d_ocp_qp_ipm_arg_set_warm_start(value, arg->qp_arg);
}


void d_ocp_qcqp_ipm_arg_set_pred_corr(int* value, struct d_ocp_qcqp_ipm_arg* arg) {
    arg->pred_corr = *value;
    d_ocp_qp_ipm_arg_set_pred_corr(value, arg->qp_arg);
}


void d_ocp_qcqp_ipm_arg_set_cond_pred_corr(int* value, struct d_ocp_qcqp_ipm_arg* arg) {
    arg->cond_pred_corr = *value;
    d_ocp_qp_ipm_arg_set_cond_pred_corr(value, arg->qp_arg);
}


void d_ocp_qcqp_ipm_arg_set_ric_alg(int* value, struct d_ocp_qcqp_ipm_arg* arg) {
    arg->square_root_alg = *value;
    d_ocp_qp_ipm_arg_set_ric_alg(value, arg->qp_arg);
}


void d_ocp_qcqp_ipm_arg_set_comp_res_exit(int* value, struct d_ocp_qcqp_ipm_arg* arg) {
    arg->comp_res_exit = *value;
    if (*value != 0)
        arg->comp_dual_sol_eq = 1;
    d_ocp_qp_ipm_arg_set_comp_res_exit(value, arg->qp_arg);
}


void d_ocp_qcqp_ipm_arg_set_comp_res_pred(int* value, struct d_ocp_qcqp_ipm_arg* arg) {
    arg->comp_res_pred = *value;
    d_ocp_qp_ipm_arg_set_comp_res_pred(value, arg->qp_arg);
}


void d_ocp_qcqp_ipm_arg_set_lam_min(double* value, struct d_ocp_qcqp_ipm_arg* arg) {
    arg->lam_min = *value;
    d_ocp_qp_ipm_arg_set_lam_min(value, arg->qp_arg);
}


void d_ocp_qcqp_ipm_arg_set_t_min(double* value, struct d_ocp_qcqp_ipm_arg* arg) {
    arg->t_min = *value;
    d_ocp_qp_ipm_arg_set_t_min(value, arg->qp_arg);
}


void d_ocp_qcqp_ipm_arg_set_split_step(int* value, struct d_ocp_qcqp_ipm_arg* arg) {
    arg->split_step = *value;
    d_ocp_qp_ipm_arg_set_split_step(value, arg->qp_arg);
}


void d_ocp_qcqp_ipm_arg_set_t_lam_min(int* value, struct d_ocp_qcqp_ipm_arg* arg) {
    arg->t_lam_min = *value;
    d_ocp_qp_ipm_arg_set_t_lam_min(value, arg->qp_arg);
}


hpipm_size_t d_ocp_qcqp_ipm_ws_strsize() {
    return sizeof(struct d_ocp_qcqp_ipm_ws);
}


hpipm_size_t d_ocp_qcqp_ipm_ws_memsize(struct d_ocp_qcqp_dim* dim, struct d_ocp_qcqp_ipm_arg* arg) {

    int N = dim->N;
    int* nu = dim->nu;
    int* nx = dim->nx;

    int ii;

    int nuM = 0;
    int nxM = 0;
    for (ii = 0; ii <= N; ii++) {
        nuM = nu[ii] > nuM ? nu[ii] : nuM;
        nxM = nx[ii] > nxM ? nx[ii] : nxM;
    }

    hpipm_size_t size = 0;

    size += 1 * sizeof(struct d_ocp_qp_ipm_ws);
    size += 1 * d_ocp_qp_ipm_ws_memsize(dim->qp_dim, arg->qp_arg);

    size += 1 * sizeof(struct d_ocp_qcqp_res_ws);  // qcqp_res_ws
    size += 1 * d_ocp_qcqp_res_ws_memsize(dim);  // qcqp_res_ws

    size += 1 * sizeof(struct d_ocp_qcqp_res);  // qcqp_res
    size += 1 * d_ocp_qcqp_res_memsize(dim);  // qcqp_res

    size += 1 * sizeof(struct d_ocp_qp);  // qp
    size += 1 * d_ocp_qp_memsize(dim->qp_dim);  // qp

    size += 1 * sizeof(struct d_ocp_qp_sol);  // qp_sol
    size += 1 * d_ocp_qp_sol_memsize(dim->qp_dim);  // qp_sol

    size += 2 * sizeof(struct vec);  // tmp_nuxM
    size += 2 * memsize_vec(nuM + nxM);  // tmp_nuxM

    size = (size + 63) / 64 * 64;  // make multiple of typical cache line size
    size += 1 * 64;  // align once to typical cache line size

    return size;
}


void d_ocp_qcqp_ipm_ws_create(struct d_ocp_qcqp_dim* dim, struct d_ocp_qcqp_ipm_arg* arg, struct d_ocp_qcqp_ipm_ws* workspace, void* mem) {

    // loop index
    int ii;

    // zero memory (to avoid corrupted memory like e.g. NaN)
    hpipm_size_t memsize = d_ocp_qcqp_ipm_ws_memsize(dim, arg);
    hpipm_zero_memset(memsize, mem);

    int N = dim->N;
    int* nu = dim->nu;
    int* nx = dim->nx;

    int nuM = 0;
    int nxM = 0;
    for (ii = 0; ii <= N; ii++) {
        nuM = nu[ii] > nuM ? nu[ii] : nuM;
        nxM = nx[ii] > nxM ? nx[ii] : nxM;
    }


    char* c_ptr = mem;


    // structures
    workspace->qp_ws = (struct d_ocp_qp_ipm_ws*) c_ptr;
    c_ptr += sizeof(struct d_ocp_qp_ipm_ws);

    workspace->qp = (struct d_ocp_qp*) c_ptr;
    c_ptr += sizeof(struct d_ocp_qp);

    workspace->qp_sol = (struct d_ocp_qp_sol*) c_ptr;
    c_ptr += sizeof(struct d_ocp_qp_sol);

    workspace->qcqp_res_ws = (struct d_ocp_qcqp_res_ws*) c_ptr;
    c_ptr += sizeof(struct d_ocp_qcqp_res_ws);

    workspace->qcqp_res = (struct d_ocp_qcqp_res*) c_ptr;
    c_ptr += sizeof(struct d_ocp_qcqp_res);


    // vector struct
    struct vec* sv_ptr = (struct vec*) c_ptr;

    workspace->tmp_nuxM = sv_ptr;
    sv_ptr += 2;


    // align to typical cache line size
    hpipm_size_t s_ptr = (hpipm_size_t) sv_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;


    // memory of structures
    c_ptr = (char*) s_ptr;

    d_ocp_qp_ipm_ws_create(dim->qp_dim, arg->qp_arg, workspace->qp_ws, c_ptr);
    c_ptr += workspace->qp_ws->memsize;

    d_ocp_qp_create(dim->qp_dim, workspace->qp, c_ptr);
    c_ptr += workspace->qp->memsize;

    d_ocp_qp_sol_create(dim->qp_dim, workspace->qp_sol, c_ptr);
    c_ptr += workspace->qp_sol->memsize;

    d_ocp_qcqp_res_ws_create(dim, workspace->qcqp_res_ws, c_ptr);
    c_ptr += workspace->qcqp_res_ws->memsize;

    d_ocp_qcqp_res_create(dim, workspace->qcqp_res, c_ptr);
    c_ptr += workspace->qcqp_res->memsize;

    create_vec(nuM + nxM, workspace->tmp_nuxM + 0, c_ptr);
    c_ptr += (workspace->tmp_nuxM + 0)->memsize;
    create_vec(nuM + nxM, workspace->tmp_nuxM + 1, c_ptr);
    c_ptr += (workspace->tmp_nuxM + 1)->memsize;


    //
    workspace->memsize = memsize;  // d_ocp_qcqp_ipm_ws_memsize(dim, arg);


#if defined(RUNTIME_CHECKS)
    if (c_ptr > ((char*) mem) + workspace->memsize) {
        printf("\nCreate_dense_qp_ipm: outside memory bounds!\n\n");
        exit(1);
    }
#endif
}


// TODO
void d_ocp_qcqp_ipm_get(char* field, struct d_ocp_qcqp_ipm_ws* ws, void* value) {
    if (hpipm_strcmp(field, "status")) {
        d_ocp_qcqp_ipm_get_status(ws, value);
    } else if (hpipm_strcmp(field, "iter")) {
        d_ocp_qcqp_ipm_get_iter(ws, value);
    } else if (hpipm_strcmp(field, "max_res_stat")) {
        d_ocp_qcqp_ipm_get_max_res_stat(ws, value);
    } else if (hpipm_strcmp(field, "max_res_eq")) {
        d_ocp_qcqp_ipm_get_max_res_eq(ws, value);
    } else if (hpipm_strcmp(field, "max_res_ineq")) {
        d_ocp_qcqp_ipm_get_max_res_ineq(ws, value);
    } else if (hpipm_strcmp(field, "max_res_comp")) {
        d_ocp_qcqp_ipm_get_max_res_comp(ws, value);
    } else if (hpipm_strcmp(field, "obj")) {
        d_ocp_qcqp_ipm_get_obj(ws, value);
    } else if (hpipm_strcmp(field, "stat")) {
        d_ocp_qcqp_ipm_get_stat(ws, value);
    } else if (hpipm_strcmp(field, "stat_m")) {
        d_ocp_qcqp_ipm_get_stat_m(ws, value);
    } else {
        printf("error: d_ocp_qcqp_ipm_get: wrong field %s\n", field);
        exit(1);
    }
}


void d_ocp_qcqp_ipm_get_status(struct d_ocp_qcqp_ipm_ws* ws, int* status) {
    *status = ws->status;
}


void d_ocp_qcqp_ipm_get_iter(struct d_ocp_qcqp_ipm_ws* ws, int* iter) {
    *iter = ws->iter;
}


void d_ocp_qcqp_ipm_get_max_res_stat(struct d_ocp_qcqp_ipm_ws* ws, double* res_stat) {
    *res_stat = ws->qcqp_res->res_max[0];
}


void d_ocp_qcqp_ipm_get_max_res_eq(struct d_ocp_qcqp_ipm_ws* ws, double* res_eq) {
    *res_eq = ws->qcqp_res->res_max[1];
}


void d_ocp_qcqp_ipm_get_max_res_ineq(struct d_ocp_qcqp_ipm_ws* ws, double* res_ineq) {
    *res_ineq = ws->qcqp_res->res_max[2];
}


void d_ocp_qcqp_ipm_get_max_res_comp(struct d_ocp_qcqp_ipm_ws* ws, double* res_comp) {
    *res_comp = ws->qcqp_res->res_max[3];
}


void d_ocp_qcqp_ipm_get_obj(struct d_ocp_qcqp_ipm_ws* ws, double* obj) {
    *obj = ws->qcqp_res->obj;
}


void d_ocp_qcqp_ipm_get_stat(struct d_ocp_qcqp_ipm_ws* ws, double** stat) {
    d_ocp_qp_ipm_get_stat(ws->qp_ws, stat);
}


void d_ocp_qcqp_ipm_get_stat_m(struct d_ocp_qcqp_ipm_ws* ws, int* stat_m) {
    d_ocp_qp_ipm_get_stat_m(ws->qp_ws, stat_m);
}


// with warm_start==2 also init dual variables
void d_ocp_qcqp_init_var(struct d_ocp_qcqp* qp, struct d_ocp_qcqp_sol* qp_sol, struct d_ocp_qcqp_ipm_arg* arg, struct d_ocp_qcqp_ipm_ws* ws) {

    //	struct CORE_qcqp_ipm_workspace *cws = ws->core_workspace;

    // loop index
    int ii, jj;

    //
    int N = qp->dim->N;
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;
    int* ns = qp->dim->ns;

    double mu0 = arg->mu0;
    double tmp;

    //
    double *ux, *s, *pi, *d_lb, *d_ub, *d_lg, *d_ug, *d_ls, *d_uq, *lam_lb, *lam_ub, *lam_lg, *lam_ug, *lam_ls, *lam_lq, *lam_uq, *t_lb, *t_ub, *t_lg, *t_ug, *t_ls, *t_lq, *t_uq;
    int *idxb, *idxs_rev;
    int idx;

    double thr0 = 1e-1;


    // primal and dual variables
    if (arg->warm_start == 2) {

        thr0 = 1e-1;

        for (ii = 0; ii <= N; ii++) {
            lam_lb = qp_sol->lam[ii].pa + 0;
            t_lb = qp_sol->t[ii].pa + 0;

            for (jj = 0; jj < 2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii]; jj++) {
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

#if 1  // old version


    // box constraints
    for (ii = 0; ii <= N; ii++) {
        ux = qp_sol->ux[ii].pa;
        d_lb = qp->d[ii].pa + 0;
        d_ub = qp->d[ii].pa + nb[ii] + ng[ii] + nq[ii];
        lam_lb = qp_sol->lam[ii].pa + 0;
        lam_ub = qp_sol->lam[ii].pa + nb[ii] + ng[ii] + nq[ii];
        t_lb = qp_sol->t[ii].pa + 0;
        t_ub = qp_sol->t[ii].pa + nb[ii] + ng[ii] + nq[ii];
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
        //		blasfeo_print_tran_dvec(nb[ii], qp->d+ii, 0);
        //		blasfeo_print_tran_dvec(nb[ii], qp->d+ii, nb[ii]+ng[ii]);
        //		blasfeo_print_tran_dvec(nu[ii]+nx[ii], qp_sol->ux+ii, 0);
        //		blasfeo_print_tran_dvec(nb[ii], qp_sol->t+ii, 0);
        //		blasfeo_print_tran_dvec(nb[ii], qp_sol->t+ii, nb[ii]+ng[ii]);
        //		exit(1);
    }
    // general constraints
    for (ii = 0; ii <= N; ii++) {
        t_lg = qp_sol->t[ii].pa + nb[ii];
        t_ug = qp_sol->t[ii].pa + 2 * nb[ii] + ng[ii] + nq[ii];
        lam_lg = qp_sol->lam[ii].pa + nb[ii];
        lam_ug = qp_sol->lam[ii].pa + 2 * nb[ii] + ng[ii] + nq[ii];
        d_lg = qp->d[ii].pa + nb[ii];
        d_ug = qp->d[ii].pa + 2 * nb[ii] + ng[ii] + nq[ii];
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
        lam_ls = qp_sol->lam[ii].pa + 2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii];
        t_ls = qp_sol->t[ii].pa + 2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii];
        for (jj = 0; jj < 2 * ns[ii]; jj++) {
            t_ls[jj] = 1.0;  // thr0;
            //			t_ls[jj] = sqrt(mu0); // thr0;
            lam_ls[jj] = mu0 / t_ls[jj];
        }
    }

    //  quadratic constraints
    double sqrt_mu0 = sqrt(mu0);
    sqrt_mu0 = thr0 > sqrt_mu0 ? thr0 : sqrt_mu0;
    double mu0_div_sqrt_mu0 = mu0 / sqrt_mu0;

    for (ii = 0; ii <= N; ii++) {
        lam_lq = qp_sol->lam[ii].pa + nb[ii] + ng[ii];
        t_lq = qp_sol->t[ii].pa + nb[ii] + ng[ii];
        lam_uq = qp_sol->lam[ii].pa + 2 * nb[ii] + 2 * ng[ii] + nq[ii];
        t_uq = qp_sol->t[ii].pa + 2 * nb[ii] + 2 * ng[ii] + nq[ii];
        d_uq = qp->d[ii].pa + 2 * nb[ii] + 2 * ng[ii] + nq[ii];
        for (jj = 0; jj < nq[ii]; jj++) {
            // disregard lower
            lam_lq[jj] = 0.0;
            t_lq[jj] = 1.0;
            // upper
#if 1
            t_uq[jj] = sqrt_mu0;
            lam_uq[jj] = mu0_div_sqrt_mu0;
#else
            //		t[2*nb+2*ng+nq+jj] = 1.0; // thr0;
            dcolex(nu[ii] + nx[ii], qp->DCt + ii, 0, ng[ii] + jj, ws->tmp_nuxM, 0);
            dsymv_l(nu[ii] + nx[ii], 0.5, &qp->Hq[ii][jj], 0, 0, qp_sol->ux + ii, 0, 1.0, ws->tmp_nuxM, 0, ws->tmp_nuxM, 0);
            tmp = ddot(nu[ii] + nx[ii], ws->tmp_nuxM, 0, qp_sol->ux + ii, 0);
            tmp = -d_uq[jj] - tmp;
            t_uq[jj] = thr0 > tmp ? thr0 : tmp;
            lam_uq[jj] = mu0 / t_uq[jj];
#endif
        }
    }


#else  // new version


    for (ii = 0; ii <= N; ii++) {

        //		printf("\nii = %d\n", ii);

        ux = qp_sol->ux[ii].pa;
        s = qp_sol->ux[ii].pa + nu[ii] + nx[ii];
        d_lb = qp->d[ii].pa + 0;
        d_ub = qp->d[ii].pa + nb[ii] + ng[ii] + nq[ii];
        d_lg = qp->d[ii].pa + nb[ii];
        d_ug = qp->d[ii].pa + 2 * nb[ii] + ng[ii] + nq[ii];
        d_ls = qp->d[ii].pa + 2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii];
        lam_lb = qp_sol->lam[ii].pa + 0;
        lam_ub = qp_sol->lam[ii].pa + nb[ii] + ng[ii] + nq[ii];
        lam_lg = qp_sol->lam[ii].pa + nb[ii];
        lam_ug = qp_sol->lam[ii].pa + 2 * nb[ii] + ng[ii] + nq[ii];
        lam_ls = qp_sol->lam[ii].pa + 2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii];
        t_lb = qp_sol->t[ii].pa + 0;
        t_ub = qp_sol->t[ii].pa + nb[ii] + ng[ii] + nq[ii];
        t_lg = qp_sol->t[ii].pa + nb[ii];
        t_ug = qp_sol->t[ii].pa + 2 * nb[ii] + ng[ii] + nq[ii];
        t_ls = qp_sol->t[ii].pa + 2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii];
        idxb = qp->idxb[ii];
        idxs_rev = qp->idxs_rev[ii];

        // lower bound on slacks
        daxpy(2 * ns[ii], -1.0, qp->d + ii, 2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii], qp_sol->ux + ii, nu[ii] + nx[ii], qp_sol->t + ii, 2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii]);
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
        //		blasfeo_print_tran_dvec(2*ns[ii], qp_sol->ux+ii, nu[ii]+nx[ii]);
        //		blasfeo_print_tran_dvec(2*ns[ii], qp_sol->t+ii, 2*nb[ii]+2*ng[ii]+2*nq[ii]);

        // upper and lower bounds on inputs and states
        dvecex_sp(nb[ii], 1.0, qp->idxb[ii], qp_sol->ux + ii, 0, qp_sol->t + ii, 0);
        dveccpsc(nb[ii], -1.0, qp_sol->t + ii, 0, qp_sol->t + ii, nb[ii] + ng[ii] + nq[ii]);
        for (jj = 0; jj < nb[ii]; jj++) {
            idx = idxs_rev[jj];
            if (idx != -1) {
                // softed bound
                t_lb[jj] += s[idx];
                t_ub[jj] += s[ns[ii] + idx];
            }
        }
        daxpy(nb[ii], -1.0, qp->d + ii, 0, qp_sol->t + ii, 0, qp_sol->t + ii, 0);
        daxpy(nb[ii], -1.0, qp->d + ii, nb[ii] + ng[ii] + nq[ii], qp_sol->t + ii, nb[ii] + ng[ii] + nq[ii], qp_sol->t + ii, nb[ii] + ng[ii] + nq[ii]);
        //		blasfeo_print_tran_dvec(nb[ii], qp_sol->t+ii, 0);
        //		blasfeo_print_tran_dvec(nb[ii], qp_sol->t+ii, nb[ii]+ng[ii]+nq[ii]);
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
        //		blasfeo_print_tran_dvec(nu[ii]+nx[ii], qp_sol->ux+ii, 0);
        //		blasfeo_print_tran_dvec(nb[ii], qp_sol->t+ii, 0);
        //		blasfeo_print_tran_dvec(nb[ii], qp_sol->t+ii, nb[ii]+ng[ii]+nq[ii]);

        // upper and lower general constaints
        dgemv_t(nu[ii] + nx[ii], ng[ii], 1.0, qp->DCt + ii, 0, 0, qp_sol->ux + ii, 0, 0.0, qp_sol->t + ii, nb[ii], qp_sol->t + ii, nb[ii]);
        dveccpsc(ng[ii], -1.0, qp_sol->t + ii, nb[ii], qp_sol->t + ii, 2 * nb[ii] + ng[ii] + nq[ii]);
        //		blasfeo_print_tran_dvec(ng[ii], qp_sol->t+ii, nb[ii]);
        //		blasfeo_print_tran_dvec(ng[ii], qp_sol->t+ii, 2*nb[ii]+ng[ii]+nq[ii]);
        for (jj = 0; jj < ng[ii]; jj++) {
            idx = idxs_rev[nb[ii] + jj];
            if (idx != -1) {
                // softed general constraint
                t_lb[nb[ii] + jj] += s[idx];
                t_ub[nb[ii] + jj] += s[ns[ii] + idx];
            }
        }
        //		blasfeo_print_tran_dvec(ng[ii], qp_sol->t+ii, nb[ii]);
        //		blasfeo_print_tran_dvec(ng[ii], qp_sol->t+ii, 2*nb[ii]+ng[ii]+nq[ii]);
        daxpy(ng[ii], -1.0, qp->d + ii, nb[ii], qp_sol->t + ii, nb[ii], qp_sol->t + ii, nb[ii]);
        daxpy(ng[ii], -1.0, qp->d + ii, 2 * nb[ii] + ng[ii] + nq[ii], qp_sol->t + ii, 2 * nb[ii] + ng[ii] + nq[ii], qp_sol->t + ii, 2 * nb[ii] + ng[ii] + nq[ii]);
        for (jj = 0; jj < ng[ii]; jj++) {
#if 1
            t_lg[jj] = thr0 > t_lg[jj] ? thr0 : t_lg[jj];
            t_ug[jj] = thr0 > t_ug[jj] ? thr0 : t_ug[jj];
#else
            t_lg[jj] = 1.0;
            t_ug[jj] = 1.0;
#endif
        }
        //		blasfeo_print_tran_dvec(ng[ii], qp_sol->t+ii, nb[ii]);
        //		blasfeo_print_tran_dvec(ng[ii], qp_sol->t+ii, 2*nb[ii]+ng[ii]+nq[ii]);

        // multipliers
        for (jj = 0; jj < 2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii]; jj++)
            lam_lb[jj] = mu0 / t_lb[jj];
    }

    // TODO nq


#endif  // new version
}


void d_ocp_qcqp_approx_qp(struct d_ocp_qcqp* qcqp, struct d_ocp_qcqp_sol* qcqp_sol, struct d_ocp_qp* qp, struct d_ocp_qcqp_ipm_ws* ws) {

    int N = qcqp->dim->N;
    int* nu = qcqp->dim->nu;
    int* nx = qcqp->dim->nx;
    int* nb = qcqp->dim->nb;
    int* ng = qcqp->dim->ng;
    int* nq = qcqp->dim->nq;
    int* ns = qcqp->dim->ns;

    double tmp;

    int ii, jj;


    for (ii = 0; ii <= N; ii++) {

        dveccp(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], qcqp->d + ii, 0, qp->d + ii, 0);

        dgecp(nu[ii] + nx[ii] + 1, nu[ii] + nx[ii], qcqp->RSQrq + ii, 0, 0, qp->RSQrq + ii, 0, 0);

        dvecse(nu[ii] + nx[ii], 0.0, ws->qcqp_res_ws->q_adj + ii, 0);

        for (jj = 0; jj < nq[ii]; jj++) {
            tmp = -VECEL(qcqp_sol->lam + ii, nb[ii] + ng[ii] + jj) + VECEL(qcqp_sol->lam + ii, 2 * nb[ii] + 2 * ng[ii] + nq[ii] + jj);
            dgead(nu[ii] + nx[ii], nu[ii] + nx[ii], tmp, &qcqp->Hq[ii][jj], 0, 0, qp->RSQrq + ii, 0, 0);

            dsymv_l(nu[ii] + nx[ii], 1.0, &qcqp->Hq[ii][jj], 0, 0, qcqp_sol->ux + ii, 0, 0.0, ws->tmp_nuxM + 0, 0, ws->tmp_nuxM + 0, 0);
            dcolex(nu[ii] + nx[ii], qcqp->DCt + ii, 0, ng[ii] + jj, ws->tmp_nuxM + 1, 0);
            daxpy(nu[ii] + nx[ii], 1.0, ws->tmp_nuxM + 0, 0, ws->tmp_nuxM + 1, 0, ws->tmp_nuxM + 1, 0);
            dcolin(nu[ii] + nx[ii], ws->tmp_nuxM + 1, 0, qp->DCt + ii, 0, ng[ii] + jj);
            daxpy(nu[ii] + nx[ii], tmp, ws->tmp_nuxM + 1, 0, ws->qcqp_res_ws->q_adj + ii, 0, ws->qcqp_res_ws->q_adj + ii, 0);

            dcolex(nu[ii] + nx[ii], qcqp->DCt + ii, 0, ng[ii] + jj, ws->tmp_nuxM + 1, 0);
            daxpy(nu[ii] + nx[ii], 0.5, ws->tmp_nuxM + 0, 0, ws->tmp_nuxM + 1, 0, ws->tmp_nuxM + 1, 0);
            tmp = ddot(nu[ii] + nx[ii], ws->tmp_nuxM + 1, 0, qcqp_sol->ux + ii, 0);
            // TODO maybe swap signs?
            VECEL(qp->d + ii, nb[ii] + ng[ii] + jj) += -tmp;
            VECEL(qp->d + ii, 2 * nb[ii] + 2 * ng[ii] + nq[ii] + jj) += +tmp;
            VECEL(ws->qcqp_res_ws->q_fun + ii, jj) = tmp;
        }

        dveccp(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], qcqp->d_mask + ii, 0, qp->d_mask + ii, 0);

        dgecp(nu[ii] + nx[ii], ng[ii], qcqp->DCt + ii, 0, 0, qp->DCt + ii, 0, 0);

        dveccp(nu[ii] + nx[ii] + 2 * ns[ii], qcqp->rqz + ii, 0, qp->rqz, 0);

        dveccp(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], qcqp->m + ii, 0, qp->m + ii, 0);

        dveccp(2 * ns[ii], qcqp->Z + ii, 0, qp->Z + ii, 0);

        for (jj = 0; jj < nb[ii]; jj++)
            qp->idxb[ii][jj] = qcqp->idxb[ii][jj];

        for (jj = 0; jj < nb[ii] + ng[ii] + nq[ii]; jj++)
            qp->idxs_rev[ii][jj] = qcqp->idxs_rev[ii][jj];
    }

    for (ii = 0; ii < N; ii++) {
        dgecp(nu[ii] + nx[ii] + 1, nx[ii + 1], qcqp->BAbt + ii, 0, 0, qp->BAbt + ii, 0, 0);

        dveccp(nx[ii + 1], qcqp->b + ii, 0, qp->b + ii, 0);
    }
}


void d_ocp_qcqp_update_qp(struct d_ocp_qcqp* qcqp, struct d_ocp_qcqp_sol* qcqp_sol, struct d_ocp_qp* qp, struct d_ocp_qcqp_ipm_ws* ws) {

    int N = qcqp->dim->N;
    int* nu = qcqp->dim->nu;
    int* nx = qcqp->dim->nx;
    int* nb = qcqp->dim->nb;
    int* ng = qcqp->dim->ng;
    int* nq = qcqp->dim->nq;
    int* ns = qcqp->dim->ns;

    double tmp;

    int ii, jj, idx;


    for (ii = 0; ii <= N; ii++) {

        // TODO only the 2*nq part needed !!!!!
        dveccp(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], qcqp->d + ii, 0, qp->d + ii, 0);

        dgecp(nu[ii] + nx[ii] + 1, nu[ii] + nx[ii], qcqp->RSQrq + ii, 0, 0, qp->RSQrq + ii, 0, 0);

        dvecse(nu[ii] + nx[ii], 0.0, ws->qcqp_res_ws->q_adj + ii, 0);

        for (jj = 0; jj < nq[ii]; jj++) {
            tmp = -VECEL(qcqp_sol->lam + ii, nb[ii] + ng[ii] + jj) + VECEL(qcqp_sol->lam + ii, 2 * nb[ii] + 2 * ng[ii] + nq[ii] + jj);
            dgead(nu[ii] + nx[ii], nu[ii] + nx[ii], tmp, &qcqp->Hq[ii][jj], 0, 0, qp->RSQrq + ii, 0, 0);

            dsymv_l(nu[ii] + nx[ii], 1.0, &qcqp->Hq[ii][jj], 0, 0, qcqp_sol->ux + ii, 0, 0.0, ws->tmp_nuxM + 0, 0, ws->tmp_nuxM + 0, 0);
            dcolex(nu[ii] + nx[ii], qcqp->DCt + ii, 0, ng[ii] + jj, ws->tmp_nuxM + 1, 0);
            daxpy(nu[ii] + nx[ii], 1.0, ws->tmp_nuxM + 0, 0, ws->tmp_nuxM + 1, 0, ws->tmp_nuxM + 1, 0);
            dcolin(nu[ii] + nx[ii], ws->tmp_nuxM + 1, 0, qp->DCt + ii, 0, ng[ii] + jj);
            daxpy(nu[ii] + nx[ii], tmp, ws->tmp_nuxM + 1, 0, ws->qcqp_res_ws->q_adj + ii, 0, ws->qcqp_res_ws->q_adj + ii, 0);

            dcolex(nu[ii] + nx[ii], qcqp->DCt + ii, 0, ng[ii] + jj, ws->tmp_nuxM + 1, 0);
            daxpy(nu[ii] + nx[ii], 0.5, ws->tmp_nuxM + 0, 0, ws->tmp_nuxM + 1, 0, ws->tmp_nuxM + 1, 0);
            tmp = ddot(nu[ii] + nx[ii], ws->tmp_nuxM + 1, 0, qcqp_sol->ux + ii, 0);
            // TODO maybe swap signs?
            VECEL(qp->d + ii, nb[ii] + ng[ii] + jj) += -tmp;
            VECEL(qp->d + ii, 2 * nb[ii] + 2 * ng[ii] + nq[ii] + jj) += +tmp;
            VECEL(ws->qcqp_res_ws->q_fun + ii, jj) = tmp;
        }

        // TODO needed ?????
        dveccp(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], qcqp->m + ii, 0, qp->m + ii, 0);
    }
}


void d_ocp_qcqp_update_qp_abs_step(struct d_ocp_qcqp* qcqp, struct d_ocp_qcqp_sol* qcqp_sol, struct d_ocp_qp* qp, struct d_ocp_qcqp_ipm_ws* ws) {

    int N = qcqp->dim->N;
    int* nu = qcqp->dim->nu;
    int* nx = qcqp->dim->nx;
    int* nb = qcqp->dim->nb;
    int* ng = qcqp->dim->ng;
    int* nq = qcqp->dim->nq;
    int* ns = qcqp->dim->ns;

    double tmp;

    int ii, jj, idx;


    for (ii = 0; ii <= N; ii++) {

        // TODO only the 2*nq part needed !!!!!
        dveccp(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], qcqp->d + ii, 0, qp->d + ii, 0);

        dgecp(nu[ii] + nx[ii] + 1, nu[ii] + ns[ii], qcqp->RSQrq + ii, 0, 0, qp->RSQrq + ii, 0, 0);

        dvecse(nu[ii] + nx[ii], 0.0, ws->qcqp_res_ws->q_adj + ii, 0);

        for (jj = 0; jj < nq[ii]; jj++) {
            //			tmp = - VECEL(qcqp_sol->lam+ii, nb[ii]+ng[ii]+jj) + VECEL(qcqp_sol->lam+ii, 2*nb[ii]+2*ng[ii]+nq[ii]+jj);
            tmp = +VECEL(qcqp_sol->lam + ii, 2 * nb[ii] + 2 * ng[ii] + nq[ii] + jj);
            dgead(nu[ii] + nx[ii], nu[ii] + nx[ii], tmp, &qcqp->Hq[ii][jj], 0, 0, qp->RSQrq + ii, 0, 0);

            dsymv_l(nu[ii] + nx[ii], 1.0, &qcqp->Hq[ii][jj], 0, 0, qcqp_sol->ux + ii, 0, 0.0, ws->tmp_nuxM + 0, 0, ws->tmp_nuxM + 0, 0);
            dcolex(nu[ii] + nx[ii], qcqp->DCt + ii, 0, ng[ii] + jj, ws->tmp_nuxM + 1, 0);
            daxpy(nu[ii] + nx[ii], 1.0, ws->tmp_nuxM + 0, 0, ws->tmp_nuxM + 1, 0, ws->tmp_nuxM + 1, 0);
            dcolin(nu[ii] + nx[ii], ws->tmp_nuxM + 1, 0, qp->DCt + ii, 0, ng[ii] + jj);
            daxpy(nu[ii] + nx[ii], tmp, ws->tmp_nuxM + 1, 0, ws->qcqp_res_ws->q_adj + ii, 0, ws->qcqp_res_ws->q_adj + ii, 0);

            //	//		daxpy(nu[ii]+nx[ii], 0.5, ws->tmp_nuxM+0, 0, &qcqp->gq[ii][jj], 0, ws->tmp_nuxM+1, 0);
            //			daxpy(nu[ii]+nx[ii], 0.5, ws->tmp_nuxM+0, 0, &qcqp->gq[ii][jj], 0, ws->tmp_nuxM+0, 0);
            //			daxpy(nu[ii]+nx[ii], -1.0, ws->tmp_nuxM+1, 0, ws->tmp_nuxM+0, 0, ws->tmp_nuxM+1, 0);
            daxpby(nu[ii] + nx[ii], -1.0, ws->tmp_nuxM + 1, 0, 0.5, ws->tmp_nuxM + 0, 0, ws->tmp_nuxM + 1, 0);
            dcolex(nu[ii] + nx[ii], qcqp->DCt + ii, 0, ng[ii] + jj, ws->tmp_nuxM + 0, 0);
            daxpy(nu[ii] + nx[ii], 1.0, ws->tmp_nuxM + 0, 0, ws->tmp_nuxM + 1, 0, ws->tmp_nuxM + 1, 0);
            tmp = ddot(nu[ii] + nx[ii], ws->tmp_nuxM + 1, 0, qcqp_sol->ux + ii, 0);
            // TODO maybe swap signs?
            VECEL(qp->d + ii, nb[ii] + ng[ii] + jj) += -tmp;
            VECEL(qp->d + ii, 2 * nb[ii] + 2 * ng[ii] + nq[ii] + jj) += +tmp;
            VECEL(ws->qcqp_res_ws->q_fun + ii, jj) = tmp;
        }

        // TODO needed ?????
        dveccp(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], qcqp->m + ii, 0, qp->m + ii, 0);
    }
}


void d_ocp_qcqp_sol_conv_qp_sol(struct d_ocp_qcqp_sol* qcqp_sol, struct d_ocp_qp_sol* qp_sol) {

    int N = qcqp_sol->dim->N;
    int* nu = qcqp_sol->dim->nu;
    int* nx = qcqp_sol->dim->nx;
    int* nb = qcqp_sol->dim->nb;
    int* ng = qcqp_sol->dim->ng;
    int* nq = qcqp_sol->dim->nq;
    int* ns = qcqp_sol->dim->ns;

    int ii;

    for (ii = 0; ii <= N; ii++)
        dveccp(nu[ii] + nx[ii] + 2 * ns[ii], qcqp_sol->ux + ii, 0, qp_sol->ux + ii, 0);

    for (ii = 0; ii < N; ii++)
        dveccp(nx[ii + 1], qcqp_sol->pi + ii, 0, qp_sol->pi + ii, 0);

    for (ii = 0; ii <= N; ii++)
        dveccp(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], qcqp_sol->lam + ii, 0, qp_sol->lam + ii, 0);

    for (ii = 0; ii <= N; ii++)
        dveccp(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], qcqp_sol->t + ii, 0, qp_sol->t + ii, 0);
}


void d_ocp_qp_sol_conv_qcqp_sol(struct d_ocp_qp_sol* qp_sol, struct d_ocp_qcqp_sol* qcqp_sol) {

    int N = qcqp_sol->dim->N;
    int* nu = qcqp_sol->dim->nu;
    int* nx = qcqp_sol->dim->nx;
    int* nb = qcqp_sol->dim->nb;
    int* ng = qcqp_sol->dim->ng;
    int* nq = qcqp_sol->dim->nq;
    int* ns = qcqp_sol->dim->ns;

    int ii;

    for (ii = 0; ii <= N; ii++)
        dveccp(nu[ii] + nx[ii] + 2 * ns[ii], qp_sol->ux + ii, 0, qcqp_sol->ux + ii, 0);

    for (ii = 0; ii < N; ii++)
        dveccp(nx[ii + 1], qp_sol->pi + ii, 0, qcqp_sol->pi + ii, 0);

    for (ii = 0; ii <= N; ii++)
        dveccp(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], qp_sol->lam + ii, 0, qcqp_sol->lam + ii, 0);

    for (ii = 0; ii <= N; ii++)
        dveccp(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], qp_sol->t + ii, 0, qcqp_sol->t + ii, 0);
}


void d_ocp_qcqp_res_conv_qp_res(struct d_ocp_qcqp_res* qcqp_res, struct d_ocp_qp_res* qp_res) {

    int N = qcqp_res->dim->N;
    int* nu = qcqp_res->dim->nu;
    int* nx = qcqp_res->dim->nx;
    int* nb = qcqp_res->dim->nb;
    int* ng = qcqp_res->dim->ng;
    int* nq = qcqp_res->dim->nq;
    int* ns = qcqp_res->dim->ns;

    int ii;

    for (ii = 0; ii <= N; ii++)
        dveccp(nu[ii] + nx[ii] + 2 * ns[ii], qcqp_res->res_g + ii, 0, qp_res->res_g + ii, 0);

    for (ii = 0; ii < N; ii++)
        dveccp(nx[ii + 1], qcqp_res->res_b + ii, 0, qp_res->res_b + ii, 0);

    for (ii = 0; ii <= N; ii++)
        dveccp(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], qcqp_res->res_d + ii, 0, qp_res->res_d + ii, 0);

    for (ii = 0; ii <= N; ii++)
        dveccp(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], qcqp_res->res_m + ii, 0, qp_res->res_m + ii, 0);

    qp_res->res_mu = qcqp_res->res_mu;
}


void d_ocp_qcqp_ipm_solve(struct d_ocp_qcqp* qcqp, struct d_ocp_qcqp_sol* qcqp_sol, struct d_ocp_qcqp_ipm_arg* qcqp_arg, struct d_ocp_qcqp_ipm_ws* qcqp_ws) {

    // dim
    int N = qcqp->dim->N;
    int* nx = qcqp->dim->nx;
    int* nu = qcqp->dim->nu;
    int* nb = qcqp->dim->nb;
    int* ng = qcqp->dim->ng;
    int* nq = qcqp->dim->nq;
    int* ns = qcqp->dim->ns;

    // extract stuff

    struct d_ocp_qp* qp = qcqp_ws->qp;
    struct d_ocp_qp_sol* qp_sol = qcqp_ws->qp_sol;
    struct d_ocp_qp_ipm_ws* qp_ws = qcqp_ws->qp_ws;
    struct d_ocp_qp_ipm_arg* qp_arg = qcqp_arg->qp_arg;

    struct d_ocp_qcqp_dim* qcqp_dim = qcqp->dim;
    struct d_ocp_qcqp_res* qcqp_res = qcqp_ws->qcqp_res;
    struct d_ocp_qcqp_res_ws* qcqp_res_ws = qcqp_ws->qcqp_res_ws;

    struct d_core_qp_ipm_workspace* cws = qp_ws->core_workspace;

    int kk, ii, jj, idx;
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
    cws->t_min_inv = qp_arg->t_min > 0 ? 1.0 / qp_arg->t_min : 1e30;
    cws->split_step = qp_arg->split_step;
    cws->t_lam_min = qp_arg->t_lam_min;

    // alias qp vectors into qp_sol
    cws->v = qp_sol->ux->pa;
    cws->pi = qp_sol->pi->pa;
    cws->lam = qp_sol->lam->pa;
    cws->t = qp_sol->t->pa;

    // alias members of qp_step
    qp_ws->qp_step->dim = qp->dim;
    qp_ws->qp_step->RSQrq = qp->RSQrq;
    qp_ws->qp_step->BAbt = qp->BAbt;
    qp_ws->qp_step->DCt = qp->DCt;
    qp_ws->qp_step->Z = qp->Z;
    qp_ws->qp_step->idxb = qp->idxb;
    qp_ws->qp_step->idxs_rev = qp->idxs_rev;
    qp_ws->qp_step->rqz = qp_ws->res->res_g;
    qp_ws->qp_step->b = qp_ws->res->res_b;
    qp_ws->qp_step->d = qp_ws->res->res_d;
    qp_ws->qp_step->m = qp_ws->res->res_m;
    qp_ws->qp_step->d_mask = qp->d_mask;

    // alias members of qp_itref
    qp_ws->qp_itref->dim = qp->dim;
    qp_ws->qp_itref->RSQrq = qp->RSQrq;
    qp_ws->qp_itref->BAbt = qp->BAbt;
    qp_ws->qp_itref->DCt = qp->DCt;
    qp_ws->qp_itref->Z = qp->Z;
    qp_ws->qp_itref->idxb = qp->idxb;
    qp_ws->qp_itref->idxs_rev = qp->idxs_rev;
    qp_ws->qp_itref->rqz = qp_ws->res_itref->res_g;
    qp_ws->qp_itref->b = qp_ws->res_itref->res_b;
    qp_ws->qp_itref->d = qp_ws->res_itref->res_d;
    qp_ws->qp_itref->m = qp_ws->res_itref->res_m;
    qp_ws->qp_itref->d_mask = qp->d_mask;

    double* qcqp_res_max = qcqp_res->res_max;


    // cache q_fun & q_adj from approx/update for res
    qcqp_ws->qcqp_res_ws->use_q_fun = 1;
    qcqp_ws->qcqp_res_ws->use_q_adj = 1;


    // disregard soft constr on (disregarded) lower quard constr
    for (ii = 0; ii <= N; ii++) {
        dvecse(nq[ii], 0.0, qcqp->d_mask + ii, nb[ii] + ng[ii]);  // TODO needed ???
        // TODO probably remove when using only idxs_rev, as the same slack may be associated with other constraints !!!!!
        for (jj = 0; jj < nq[ii]; jj++) {
            idx = qcqp->idxs_rev[ii][nb[ii] + ng[ii] + jj];
            if (idx != -1) {
                VECEL(qcqp->d_mask + ii, 2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + idx) = 0.0;
            }
        }
    }


    // initialize qcqp & qp
    d_ocp_qcqp_init_var(qcqp, qcqp_sol, qcqp_arg, qcqp_ws);
    // mask out disregarded constraints
    dvecmul(cws->nc, qcqp->d_mask, 0, qcqp_sol->lam, 0, qcqp_sol->lam, 0);
    d_ocp_qcqp_sol_conv_qp_sol(qcqp_sol, qp_sol);

    // approximate qcqp with a qp
    d_ocp_qcqp_approx_qp(qcqp, qcqp_sol, qp, qcqp_ws);


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
        d_ocp_qp_fact_solve_kkt_unconstr(qp, qp_sol, qp_arg, qp_ws);
        d_ocp_qp_sol_conv_qcqp_sol(qp_sol, qcqp_sol);
        if (qp_arg->comp_res_exit) {
            // compute residuals
            d_ocp_qcqp_res_compute(qcqp, qcqp_sol, qcqp_res, qcqp_res_ws);
            // XXX no constraints, so no mask
            d_ocp_qcqp_res_compute_inf_norm(qcqp_res);
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

        if (isnan(VECEL(qcqp_sol->ux + 0, 0))) {
            // NaN in the solution
            qcqp_ws->status = NAN_SOL;
        }
        // #else
        //         if (VECEL(qcqp_sol->ux + 0, 0) != VECEL(qcqp_sol->ux + 0, 0)) {
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
        qp_ws->qp_step->RSQrq = qp->RSQrq;
        qp_ws->qp_step->BAbt = qp->BAbt;
        qp_ws->qp_step->DCt = qp->DCt;
        qp_ws->qp_step->Z = qp->Z;
        qp_ws->qp_step->idxb = qp->idxb;
        qp_ws->qp_step->idxe = qp->idxe;
        qp_ws->qp_step->idxs_rev = qp->idxs_rev;
        qp_ws->qp_step->rqz = qp->rqz;
        qp_ws->qp_step->b = qp->b;
        qp_ws->qp_step->d = qp->d;
        qp_ws->qp_step->d_mask = qp->d_mask;
        qp_ws->qp_step->m = qp_ws->tmp_m;

        // alias core workspace
        cws->res_m = qp_ws->qp_step->m->pa;

        // update approximation of qcqp as qp for absolute step
        d_ocp_qcqp_update_qp_abs_step(qcqp, qcqp_sol, qp, qcqp_ws);

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
            d_ocp_qp_ipm_abs_step(kk, qp, qp_sol, qp_arg, qp_ws);
            // blasfeo_print_exp_tran_dvec(cws->nc, qp_sol->lam, 0);
            d_ocp_qp_sol_conv_qcqp_sol(qp_sol, qcqp_sol);

            // update approximation of qcqp as qp for absolute step
            d_ocp_qcqp_update_qp_abs_step(qcqp, qcqp_sol, qp, qcqp_ws);

            // compute mu
            mu = dvecmuldot(cws->nc, qcqp_sol->lam, 0, qcqp_sol->t, 0, qp_ws->tmp_m, 0);
            mu /= cws->nc;
            cws->mu = mu;
            if (kk + 1 < stat_max)
                stat[stat_m * (kk + 1) + 5] = mu;
        }

        if (qp_arg->comp_res_exit & qp_arg->comp_dual_sol_eq) {
            // compute residuals
            d_ocp_qcqp_res_compute(qcqp, qcqp_sol, qcqp_res, qcqp_res_ws);
            if (qp_ws->mask_constr) {
                // mask out disregarded constraints
                for (ii = 0; ii <= N; ii++)
                    dvecmul(2 * ns[ii], qp->d_mask + ii, 2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii], qcqp_res->res_g + ii, nu[ii] + nx[ii], qcqp_res->res_g + ii, nu[ii] + nx[ii]);
                dvecmul(cws->nc, qp->d_mask, 0, qcqp_res->res_d, 0, qcqp_res->res_d, 0);
                dvecmul(cws->nc, qp->d_mask, 0, qcqp_res->res_m, 0, qcqp_res->res_m, 0);
            }
            d_ocp_qcqp_res_compute_inf_norm(qcqp_res);
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

    } else {
        // compute residuals
        d_ocp_qcqp_res_compute(qcqp, qcqp_sol, qcqp_res, qcqp_res_ws);
        if (qp_ws->mask_constr) {
            // mask out disregarded constraints
            for (ii = 0; ii <= N; ii++)
                dvecmul(2 * ns[ii], qp->d_mask + ii, 2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii], qcqp_res->res_g + ii, nu[ii] + nx[ii], qcqp_res->res_g + ii, nu[ii] + nx[ii]);
            dvecmul(cws->nc, qp->d_mask, 0, qcqp_res->res_d, 0, qcqp_res->res_d, 0);
            dvecmul(cws->nc, qp->d_mask, 0, qcqp_res->res_m, 0, qcqp_res->res_m, 0);
        }
        d_ocp_qcqp_res_compute_inf_norm(qcqp_res);
        d_ocp_qcqp_res_conv_qp_res(qcqp_res, qp_ws->res);
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
             kk<qcqp_arg->iter_max & cws->alpha> qcqp_arg->alpha_min &
             (qcqp_res_max[0] > qcqp_arg->res_g_max |
              qcqp_res_max[1] > qcqp_arg->res_b_max |
              qcqp_res_max[2] > qcqp_arg->res_d_max |
              qcqp_res_max[3] > qcqp_arg->res_m_max);
             kk++) {

            // hessian is updated with quad constr: can not reuse hessian factorization !!!
            // XXX is it in ocp ws ?????
            for (ii = 0; ii <= N; ii++)
                qp_ws->use_hess_fact[ii] = 0;

            // compute delta step
            d_ocp_qp_ipm_delta_step(kk, qp, qp_sol, qp_arg, qp_ws);
            // blasfeo_print_exp_tran_dvec(cws->nc, qp_sol->lam, 0);
            d_ocp_qp_sol_conv_qcqp_sol(qp_sol, qcqp_sol);
            // XXX maybe not needed
            if (qp_ws->mask_constr) {
                // mask out disregarded constraints
                dvecmul(cws->nc, qp->d_mask, 0, qcqp_sol->lam, 0, qcqp_sol->lam, 0);
            }

            // update approximation of qcqp as qp
            d_ocp_qcqp_update_qp(qcqp, qcqp_sol, qp, qcqp_ws);

            // compute residuals
            d_ocp_qcqp_res_compute(qcqp, qcqp_sol, qcqp_res, qcqp_res_ws);
            if (qp_ws->mask_constr) {
                // mask out disregarded constraints
                for (ii = 0; ii <= N; ii++)
                    dvecmul(2 * ns[ii], qp->d_mask + ii, 2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii], qcqp_res->res_g + ii, nu[ii] + nx[ii], qcqp_res->res_g + ii, nu[ii] + nx[ii]);
                dvecmul(cws->nc, qp->d_mask, 0, qcqp_res->res_d, 0, qcqp_res->res_d, 0);
                dvecmul(cws->nc, qp->d_mask, 0, qcqp_res->res_m, 0, qcqp_res->res_m, 0);
            }
            d_ocp_qcqp_res_compute_inf_norm(qcqp_res);
            d_ocp_qcqp_res_conv_qp_res(qcqp_res, qp_ws->res);
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
    }

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
