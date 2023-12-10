#ifndef HPIPM_S_OCP_QP_KKT_H_
#define HPIPM_S_OCP_QP_KKT_H_

#include "hpipm/common.h"
#include "hpipm/ocp/d_ocp_qp_ipm.h"
#include "hpipm/ocp/s_ocp_qp.h"
#include "hpipm/ocp/s_ocp_qp_sol.h"


#ifdef __cplusplus
extern "C" {
#endif

//
void s_ocp_qp_fact_solve_kkt_unconstr(struct s_ocp_qp* qp, struct s_ocp_qp_sol* qp_sol, struct s_ocp_qp_ipm_arg* arg, struct s_ocp_qp_ipm_ws* ws);
//
void s_ocp_qp_fact_solve_kkt_step(struct s_ocp_qp* qp, struct s_ocp_qp_sol* qp_sol, struct s_ocp_qp_ipm_arg* arg, struct s_ocp_qp_ipm_ws* ws);
//
void s_ocp_qp_fact_lq_solve_kkt_step(struct s_ocp_qp* qp, struct s_ocp_qp_sol* qp_sol, struct s_ocp_qp_ipm_arg* arg, struct s_ocp_qp_ipm_ws* ws);
//
void s_ocp_qp_solve_kkt_step(struct s_ocp_qp* qp, struct s_ocp_qp_sol* qp_sol, struct s_ocp_qp_ipm_arg* arg, struct s_ocp_qp_ipm_ws* ws);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif  // HPIPM_S_OCP_QP_KKT_H_
