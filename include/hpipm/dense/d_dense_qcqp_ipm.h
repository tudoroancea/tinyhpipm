#ifndef HPIPM_D_DENSE_QCQP_IPM_H_
#define HPIPM_D_DENSE_QCQP_IPM_H_

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/common.h"
#include "hpipm/dense/d_dense_qcqp.h"
#include "hpipm/dense/d_dense_qcqp_dim.h"
#include "hpipm/dense/d_dense_qcqp_res.h"
#include "hpipm/dense/d_dense_qcqp_sol.h"


#ifdef __cplusplus
extern "C" {
#endif


struct d_dense_qcqp_ipm_arg {
    struct d_dense_qp_ipm_arg* qp_arg;
    double mu0;  // initial value for duality measure
    double alpha_min;  // exit cond on step length
    double res_g_max;  // exit cond on inf norm of residuals
    double res_b_max;  // exit cond on inf norm of residuals
    double res_d_max;  // exit cond on inf norm of residuals
    double res_m_max;  // exit cond on inf norm of residuals
    double reg_prim;  // reg of primal hessian
    double reg_dual;  // reg of dual hessian
    double lam_min;  // min value in lam vector
    double t_min;  // min value in t vector
    int iter_max;  // exit cond in iter number
    int stat_max;  // iterations saved in stat
    int pred_corr;  // Mehrotra's predictor-corrector IPM algirthm
    int cond_pred_corr;  // conditional Mehrotra's predictor-corrector
    int scale;  // scale hessian
    int itref_pred_max;  // max number of iterative refinement steps for predictor step
    int itref_corr_max;  // max number of iterative refinement steps for corrector step
    int warm_start;  // 0 no warm start, 1 warm start primal sol, 2 warm start primal and dual sol
    int lq_fact;  // 0 syrk+potrf, 1 mix, 2 lq
    int abs_form;  // absolute IPM formulation
    int comp_res_exit;  // compute residuals on exit (only for abs_form==1)
    int comp_res_pred;  // compute residuals of prediction
    int split_step;  // use different step for primal and dual variables
    int t_lam_min;  // clip t and lam: 0 no, 1 in Gamma computation, 2 in solution
    int mode;
    hpipm_size_t memsize;
};


struct d_dense_qcqp_ipm_ws {
    struct d_dense_qp_ipm_ws* qp_ws;
    struct d_dense_qp* qp;
    struct d_dense_qp_sol* qp_sol;
    struct d_dense_qcqp_res_ws* qcqp_res_ws;
    struct d_dense_qcqp_res* qcqp_res;
    struct blasfeo_dvec* tmp_nv;
    int iter;  // iteration number
    int status;
    hpipm_size_t memsize;  // memory size (in bytes) of workspace
};


//
hpipm_size_t d_dense_qcqp_ipm_arg_strsize();
//
hpipm_size_t d_dense_qcqp_ipm_arg_memsize(struct d_dense_qcqp_dim* dim);
//
void d_dense_qcqp_ipm_arg_create(struct d_dense_qcqp_dim* dim, struct d_dense_qcqp_ipm_arg* arg, void* mem);
//
void d_dense_qcqp_ipm_arg_set_default(enum hpipm_mode mode, struct d_dense_qcqp_ipm_arg* arg);
//
void d_dense_qcqp_ipm_arg_set(char* field, void* value, struct d_dense_qcqp_ipm_arg* arg);
//
void d_dense_qcqp_ipm_arg_set_iter_max(int* iter_max, struct d_dense_qcqp_ipm_arg* arg);
//
void d_dense_qcqp_ipm_arg_set_alpha_min(double* alpha_min, struct d_dense_qcqp_ipm_arg* arg);
//
void d_dense_qcqp_ipm_arg_set_mu0(double* mu0, struct d_dense_qcqp_ipm_arg* arg);
//
void d_dense_qcqp_ipm_arg_set_tol_stat(double* tol_stat, struct d_dense_qcqp_ipm_arg* arg);
//
void d_dense_qcqp_ipm_arg_set_tol_eq(double* tol_eq, struct d_dense_qcqp_ipm_arg* arg);
//
void d_dense_qcqp_ipm_arg_set_tol_ineq(double* tol_ineq, struct d_dense_qcqp_ipm_arg* arg);
//
void d_dense_qcqp_ipm_arg_set_tol_comp(double* tol_comp, struct d_dense_qcqp_ipm_arg* arg);
//
void d_dense_qcqp_ipm_arg_set_reg_prim(double* reg, struct d_dense_qcqp_ipm_arg* arg);
//
void d_dense_qcqp_ipm_arg_set_reg_dual(double* reg, struct d_dense_qcqp_ipm_arg* arg);
//
void d_dense_qcqp_ipm_arg_set_warm_start(int* warm_start, struct d_dense_qcqp_ipm_arg* arg);
//
void d_dense_qcqp_ipm_arg_set_pred_corr(int* pred_corr, struct d_dense_qcqp_ipm_arg* arg);
//
void d_dense_qcqp_ipm_arg_set_cond_pred_corr(int* cond_pred_corr, struct d_dense_qcqp_ipm_arg* arg);
//
void d_dense_qcqp_ipm_arg_set_comp_res_pred(int* comp_res_pred, struct d_dense_qcqp_ipm_arg* arg);
//
void d_dense_qcqp_ipm_arg_set_comp_res_exit(int* comp_res_exit, struct d_dense_qcqp_ipm_arg* arg);
//
void d_dense_qcqp_ipm_arg_set_lam_min(double* value, struct d_dense_qcqp_ipm_arg* arg);
//
void d_dense_qcqp_ipm_arg_set_t_min(double* value, struct d_dense_qcqp_ipm_arg* arg);
//
void d_dense_qcqp_ipm_arg_set_split_step(int* value, struct d_dense_qcqp_ipm_arg* arg);
//
void d_dense_qcqp_ipm_arg_set_t_lam_min(int* value, struct d_dense_qcqp_ipm_arg* arg);


//
hpipm_size_t d_dense_qcqp_ipm_ws_strsize();
//
hpipm_size_t d_dense_qcqp_ipm_ws_memsize(struct d_dense_qcqp_dim* qp_dim, struct d_dense_qcqp_ipm_arg* arg);
//
void d_dense_qcqp_ipm_ws_create(struct d_dense_qcqp_dim* qp_dim, struct d_dense_qcqp_ipm_arg* arg, struct d_dense_qcqp_ipm_ws* ws, void* mem);
//
void d_dense_qcqp_ipm_get(char* field, struct d_dense_qcqp_ipm_ws* ws, void* value);
//
void d_dense_qcqp_ipm_get_status(struct d_dense_qcqp_ipm_ws* ws, int* status);
//
void d_dense_qcqp_ipm_get_iter(struct d_dense_qcqp_ipm_ws* ws, int* iter);
//
void d_dense_qcqp_ipm_get_max_res_stat(struct d_dense_qcqp_ipm_ws* ws, double* res_stat);
//
void d_dense_qcqp_ipm_get_max_res_eq(struct d_dense_qcqp_ipm_ws* ws, double* res_eq);
//
void d_dense_qcqp_ipm_get_max_res_ineq(struct d_dense_qcqp_ipm_ws* ws, double* res_ineq);
//
void d_dense_qcqp_ipm_get_max_res_comp(struct d_dense_qcqp_ipm_ws* ws, double* res_comp);
//
void d_dense_qcqp_ipm_get_obj(struct d_dense_qcqp_ipm_ws* ws, double* obj);
//
void d_dense_qcqp_ipm_get_stat(struct d_dense_qcqp_ipm_ws* ws, double** stat);
//
void d_dense_qcqp_ipm_get_stat_m(struct d_dense_qcqp_ipm_ws* ws, int* stat_m);
//
void d_dense_qcqp_init_var(struct d_dense_qcqp* qp, struct d_dense_qcqp_sol* qp_sol, struct d_dense_qcqp_ipm_arg* arg, struct d_dense_qcqp_ipm_ws* ws);
//
void d_dense_qcqp_ipm_solve(struct d_dense_qcqp* qp, struct d_dense_qcqp_sol* qp_sol, struct d_dense_qcqp_ipm_arg* arg, struct d_dense_qcqp_ipm_ws* ws);
#if 0
//
void d_dense_qcqp_ipm_predict(struct d_dense_qcqp *qp, struct d_dense_qcqp_sol *qp_sol, struct d_dense_qcqp_ipm_arg *arg, struct d_dense_qcqp_ipm_ws *ws);
//
void d_dense_qcqp_ipm_sens(struct d_dense_qcqp *qp, struct d_dense_qcqp_sol *qp_sol, struct d_dense_qcqp_ipm_arg *arg, struct d_dense_qcqp_ipm_ws *ws);
#endif


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif  // HPIPM_D_DENSE_QCQP_IPM_H_
