#ifndef HPIPM_S_DENSE_QP_KKT_H_
#define HPIPM_S_DENSE_QP_KKT_H_

#include "hpipm/common.h"
#include "hpipm/dense/s_dense_qp.h"
#include "hpipm/dense/s_dense_qp_ipm.h"
#include "hpipm/dense/s_dense_qp_sol.h"


#ifdef __cplusplus
extern "C" {
#endif


//
void s_fact_solve_kkt_unconstr_dense_qp(struct s_dense_qp* qp, struct s_dense_qp_sol* qp_sol, struct s_dense_qp_ipm_arg* arg, struct s_dense_qp_ipm_ws* ws);
//
void s_fact_solve_kkt_step_dense_qp(struct s_dense_qp* qp, struct s_dense_qp_sol* qp_sol, struct s_dense_qp_ipm_arg* arg, struct s_dense_qp_ipm_ws* ws);
//
void s_fact_lq_solve_kkt_step_dense_qp(struct s_dense_qp* qp, struct s_dense_qp_sol* qp_sol, struct s_dense_qp_ipm_arg* arg, struct s_dense_qp_ipm_ws* ws);
//
void s_solve_kkt_step_dense_qp(struct s_dense_qp* qp, struct s_dense_qp_sol* qp_sol, struct s_dense_qp_ipm_arg* arg, struct s_dense_qp_ipm_ws* ws);
//
void s_dense_qp_remove_lin_dep_eq(struct s_dense_qp* qp, struct s_dense_qp_ipm_arg* arg, struct s_dense_qp_ipm_ws* ws);
//
void s_dense_qp_restore_lin_dep_eq(struct s_dense_qp* qp, struct s_dense_qp_ipm_arg* arg, struct s_dense_qp_ipm_ws* ws);
//
void s_dense_qp_compute_obj(struct s_dense_qp* qp, struct s_dense_qp_sol* qp_sol, struct s_dense_qp_ipm_arg* arg, struct s_dense_qp_ipm_ws* ws);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif  // HPIPM_S_DENSE_QP_KKT_H_
