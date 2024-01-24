#ifndef HPIPM_D_d_dense_qp_KKT_H_
#define HPIPM_D_d_dense_qp_KKT_H_

#include "hpipm/common.h"
#include "hpipm/dense/d_dense_qp.h"
#include "hpipm/dense/d_dense_qp_ipm.h"
#include "hpipm/dense/d_dense_qp_sol.h"


#ifdef __cplusplus
extern "C" {
#endif


//
void d_fact_solve_kkt_unconstr_dense_qp(struct d_dense_qp* qp, struct d_dense_qp_sol* qp_sol, struct d_dense_qp_ipm_arg* arg, struct d_dense_qp_ipm_ws* ws);
//
void d_fact_solve_kkt_step_dense_qp(struct d_dense_qp* qp, struct d_dense_qp_sol* qp_sol, struct d_dense_qp_ipm_arg* arg, struct d_dense_qp_ipm_ws* ws);
//
void d_fact_lq_solve_kkt_step_dense_qp(struct d_dense_qp* qp, struct d_dense_qp_sol* qp_sol, struct d_dense_qp_ipm_arg* arg, struct d_dense_qp_ipm_ws* ws);
//
void d_solve_kkt_step_dense_qp(struct d_dense_qp* qp, struct d_dense_qp_sol* qp_sol, struct d_dense_qp_ipm_arg* arg, struct d_dense_qp_ipm_ws* ws);
//
void d_dense_qp_remove_lin_dep_eq(struct d_dense_qp* qp, struct d_dense_qp_ipm_arg* arg, struct d_dense_qp_ipm_ws* ws);
//
void d_dense_qp_restore_lin_dep_eq(struct d_dense_qp* qp, struct d_dense_qp_ipm_arg* arg, struct d_dense_qp_ipm_ws* ws);
//
void d_dense_qp_compute_obj(struct d_dense_qp* qp, struct d_dense_qp_sol* qp_sol, struct d_dense_qp_ipm_arg* arg, struct d_dense_qp_ipm_ws* ws);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif  // HPIPM_D_d_dense_qp_KKT_H_
