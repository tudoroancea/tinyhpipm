#ifndef HPIPM_D_OCP_qp_ipm_H_
#define HPIPM_D_OCP_qp_ipm_H_

#include "tinyhpipm/blas.h"

#include "tinyhpipm/common.h"
#include "tinyhpipm/ocp/d_ocp_qp.h"
#include "tinyhpipm/ocp/d_ocp_qp_dim.h"
#include "tinyhpipm/ocp/d_ocp_qp_res.h"
#include "tinyhpipm/ocp/d_ocp_qp_sol.h"


#ifdef __cplusplus
extern "C" {
#endif


struct d_ocp_qp_ipm_arg {
    double mu0;  // initial value for complementarity slackness
    double alpha_min;  // exit cond on step length
    double res_g_max;  // exit cond on inf norm of residuals
    double res_b_max;  // exit cond on inf norm of residuals
    double res_d_max;  // exit cond on inf norm of residuals
    double res_m_max;  // exit cond on inf norm of residuals
    double reg_prim;  // reg of primal hessian
    double lam_min;  // min value in lam vector
    double t_min;  // min value in t vector
    double tau_min;  // min value of barrier parameter
    int iter_max;  // exit cond in iter number
    int stat_max;  // iterations saved in stat
    int pred_corr;  // use Mehrotra's predictor-corrector IPM algirthm
    int cond_pred_corr;  // conditional Mehrotra's predictor-corrector
    int itref_pred_max;  // max number of iterative refinement steps for predictor step
    int itref_corr_max;  // max number of iterative refinement steps for corrector step
    int warm_start;  // 0 no warm start, 1 warm start primal sol, 2 warm start primal and dual sol
    int square_root_alg;  // 0 classical Riccati, 1 square-root Riccati
    int lq_fact;  // 0 syrk+potrf, 1 mix, 2 lq (for square_root_alg==1)
    int abs_form;  // absolute IPM formulation
    int comp_dual_sol_eq;  // dual solution of equality constraints (only for abs_form==1)
    int comp_res_exit;  // compute residuals on exit (only for abs_form==1 and comp_dual_sol_eq==1)
    int comp_res_pred;  // compute residuals of prediction
    int split_step;  // use different steps for primal and dual variables
    int var_init_scheme;  // variables initialization scheme
    int t_lam_min;  // clip t and lam: 0 no, 1 in Gamma computation, 2 in solution
    int mode;
    hpipm_size_t memsize;
};


struct d_ocp_qp_ipm_ws {
    double qp_res[4];  // infinity norm of residuals
    struct d_core_qp_ipm_workspace* core_workspace;
    struct d_ocp_qp_dim* dim;  // cache dim
    struct d_ocp_qp_res_ws* res_workspace;
    struct d_ocp_qp_sol* sol_step;
    struct d_ocp_qp_sol* sol_itref;
    struct d_ocp_qp* qp_step;
    struct d_ocp_qp* qp_itref;
    struct d_ocp_qp_res* res_itref;
    struct d_ocp_qp_res* res;
    struct vec* Gamma;  // hessian update
    struct vec* gamma;  // hessian update
    struct vec* tmp_nuxM;  // work space of size nxM
    struct vec* tmp_nbgM;  // work space of size nbM+ngM
    struct vec* tmp_nsM;  // work space of size nsM
    struct vec* Pb;  // Pb
    struct vec* Zs_inv;
    struct vec* tmp_m;
    struct vec* l;  // cache linear part for _get_ric_xxx
    struct mat* L;
    struct mat* Ls;
    struct mat* P;
    struct mat* Lh;
    struct mat* AL;
    struct mat* lq0;
    struct mat* tmp_nxM_nxM;
    double* stat;  // convergence statistics
    int* use_hess_fact;
    void* lq_work0;
    int iter;  // iteration number
    int stat_max;  // iterations saved in stat
    int stat_m;  // number of recorded stat per IPM iter
    int use_Pb;
    int status;  // solver status
    int square_root_alg;  // cache from arg
    int lq_fact;  // cache from arg
    int mask_constr;  // use constr mask
    int valid_ric_vec;  // meaningful riccati vectors
    int valid_ric_p;  // form of riccati p: 0 p*inv(L), 1 p
    hpipm_size_t memsize;
};


//
hpipm_size_t d_ocp_qp_ipm_arg_strsize();
//
hpipm_size_t d_ocp_qp_ipm_arg_memsize(struct d_ocp_qp_dim* ocp_dim);
//
void d_ocp_qp_ipm_arg_create(struct d_ocp_qp_dim* ocp_dim, struct d_ocp_qp_ipm_arg* arg, void* mem);
//
void d_ocp_qp_ipm_arg_set_default(enum tinyhpipm_mode mode, struct d_ocp_qp_ipm_arg* arg);
//
void d_ocp_qp_ipm_arg_set(char* field, void* value, struct d_ocp_qp_ipm_arg* arg);
// set maximum number of iterations
void d_ocp_qp_ipm_arg_set_iter_max(int* iter_max, struct d_ocp_qp_ipm_arg* arg);
// set minimum step lenght
void d_ocp_qp_ipm_arg_set_alpha_min(double* alpha_min, struct d_ocp_qp_ipm_arg* arg);
// set initial value of barrier parameter
void d_ocp_qp_ipm_arg_set_mu0(double* mu0, struct d_ocp_qp_ipm_arg* arg);
// set exit tolerance on stationarity condition
void d_ocp_qp_ipm_arg_set_tol_stat(double* tol_stat, struct d_ocp_qp_ipm_arg* arg);
// set exit tolerance on equality constr
void d_ocp_qp_ipm_arg_set_tol_eq(double* tol_eq, struct d_ocp_qp_ipm_arg* arg);
// set exit tolerance on inequality constr
void d_ocp_qp_ipm_arg_set_tol_ineq(double* tol_ineq, struct d_ocp_qp_ipm_arg* arg);
// set exit tolerance on complementarity condition
void d_ocp_qp_ipm_arg_set_tol_comp(double* tol_comp, struct d_ocp_qp_ipm_arg* arg);
// set regularization of primal variables
void d_ocp_qp_ipm_arg_set_reg_prim(double* tol_comp, struct d_ocp_qp_ipm_arg* arg);
// set warm start: 0 no warm start, 1 primal var
void d_ocp_qp_ipm_arg_set_warm_start(int* warm_start, struct d_ocp_qp_ipm_arg* arg);
// Mehrotra's predictor-corrector IPM algorithm: 0 no predictor-corrector, 1 use predictor-corrector
void d_ocp_qp_ipm_arg_set_pred_corr(int* pred_corr, struct d_ocp_qp_ipm_arg* arg);
// conditional predictor-corrector: 0 no conditinal predictor-corrector, 1 conditional predictor-corrector
void d_ocp_qp_ipm_arg_set_cond_pred_corr(int* value, struct d_ocp_qp_ipm_arg* arg);
// set riccati algorithm: 0 classic, 1 square-root
void d_ocp_qp_ipm_arg_set_ric_alg(int* value, struct d_ocp_qp_ipm_arg* arg);
// dual solution of equality constraints (only for abs_form==1)
void d_ocp_qp_ipm_arg_set_comp_dual_sol_eq(int* value, struct d_ocp_qp_ipm_arg* arg);
// compute residuals after solution
void d_ocp_qp_ipm_arg_set_comp_res_exit(int* value, struct d_ocp_qp_ipm_arg* arg);
// compute residuals of prediction
void d_ocp_qp_ipm_arg_set_comp_res_pred(int* value, struct d_ocp_qp_ipm_arg* arg);
// min value of lam in the solution
void d_ocp_qp_ipm_arg_set_lam_min(double* value, struct d_ocp_qp_ipm_arg* arg);
// min value of t in the solution
void d_ocp_qp_ipm_arg_set_t_min(double* value, struct d_ocp_qp_ipm_arg* arg);
// min value of tau in the solution
void d_ocp_qp_ipm_arg_set_tau_min(double* value, struct d_ocp_qp_ipm_arg* arg);
// set split step: 0 same step, 1 different step for primal and dual variables
void d_ocp_qp_ipm_arg_set_split_step(int* value, struct d_ocp_qp_ipm_arg* arg);
// variables initialization scheme
void d_ocp_qp_ipm_arg_set_var_init_scheme(int* value, struct d_ocp_qp_ipm_arg* arg);
// clip t and lam: 0 no, 1 in Gamma computation, 2 in solution
void d_ocp_qp_ipm_arg_set_t_lam_min(int* value, struct d_ocp_qp_ipm_arg* arg);

//
hpipm_size_t d_ocp_qp_ipm_ws_strsize();
//
hpipm_size_t d_ocp_qp_ipm_ws_memsize(struct d_ocp_qp_dim* ocp_dim, struct d_ocp_qp_ipm_arg* arg);
//
void d_ocp_qp_ipm_ws_create(struct d_ocp_qp_dim* ocp_dim, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws, void* mem);
//
void d_ocp_qp_ipm_get(char* field, struct d_ocp_qp_ipm_ws* ws, void* value);
//
void d_ocp_qp_ipm_get_status(struct d_ocp_qp_ipm_ws* ws, int* status);
//
void d_ocp_qp_ipm_get_iter(struct d_ocp_qp_ipm_ws* ws, int* iter);
//
void d_ocp_qp_ipm_get_max_res_stat(struct d_ocp_qp_ipm_ws* ws, double* res_stat);
//
void d_ocp_qp_ipm_get_max_res_eq(struct d_ocp_qp_ipm_ws* ws, double* res_eq);
//
void d_ocp_qp_ipm_get_max_res_ineq(struct d_ocp_qp_ipm_ws* ws, double* res_ineq);
//
void d_ocp_qp_ipm_get_max_res_comp(struct d_ocp_qp_ipm_ws* ws, double* res_comp);
//
void d_ocp_qp_ipm_get_obj(struct d_ocp_qp_ipm_ws* ws, double* obj);
//
void d_ocp_qp_ipm_get_stat(struct d_ocp_qp_ipm_ws* ws, double** stat);
//
void d_ocp_qp_ipm_get_stat_m(struct d_ocp_qp_ipm_ws* ws, int* stat_m);
//
void d_ocp_qp_ipm_get_ric_Lr(struct d_ocp_qp* qp, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws, int stage, double* Lr);
//
void d_ocp_qp_ipm_get_ric_Ls(struct d_ocp_qp* qp, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws, int stage, double* Ls);
//
void d_ocp_qp_ipm_get_ric_P(struct d_ocp_qp* qp, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws, int stage, double* P);
//
void d_ocp_qp_ipm_get_ric_lr(struct d_ocp_qp* qp, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws, int stage, double* lr);
//
void d_ocp_qp_ipm_get_ric_p(struct d_ocp_qp* qp, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws, int stage, double* p);
// feedback control gain in the form u = K x + k
// void d_ocp_qp_ipm_get_ric_K(struct d_ocp_qp* qp, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws, int stage, double* K);
// feedback control gain in the form u = K x + k
// void d_ocp_qp_ipm_get_ric_k(struct d_ocp_qp* qp, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws, int stage, double* k);
//
void d_ocp_qp_init_var(struct d_ocp_qp* qp, struct d_ocp_qp_sol* qp_sol, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws);
//
void d_ocp_qp_ipm_abs_step(int kk, struct d_ocp_qp* qp, struct d_ocp_qp_sol* qp_sol, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws);
//
void d_ocp_qp_ipm_delta_step(int kk, struct d_ocp_qp* qp, struct d_ocp_qp_sol* qp_sol, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws);
//
void d_ocp_qp_ipm_solve(struct d_ocp_qp* qp, struct d_ocp_qp_sol* qp_sol, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws);
//
void d_ocp_qp_ipm_predict(struct d_ocp_qp* qp, struct d_ocp_qp_sol* qp_sol, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws);
//
void d_ocp_qp_ipm_sens(struct d_ocp_qp* qp, struct d_ocp_qp_sol* qp_sol, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws);


#ifdef __cplusplus
}  // #extern "C"
#endif


#endif  // HPIPM_D_OCP_qp_ipm_H_
