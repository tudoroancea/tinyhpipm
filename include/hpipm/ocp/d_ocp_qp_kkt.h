#ifndef HPIPM_D_OCP_qp_KKT_H_
#define HPIPM_D_OCP_qp_KKT_H_

#include "hpipm/common.h"
#include "hpipm/ocp/d_ocp_qp.h"
#include "hpipm/ocp/d_ocp_qp_ipm.h"
#include "hpipm/ocp/d_ocp_qp_sol.h"


#ifdef __cplusplus
extern "C" {
#endif


//
void d_ocp_qp_fact_solve_kkt_unconstr(struct d_ocp_qp* qp, struct d_ocp_qp_sol* qp_sol, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws);
//
void d_ocp_qp_fact_solve_kkt_step(struct d_ocp_qp* qp, struct d_ocp_qp_sol* qp_sol, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws);
//
void d_ocp_qp_fact_lq_solve_kkt_step(struct d_ocp_qp* qp, struct d_ocp_qp_sol* qp_sol, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws);
//
void d_ocp_qp_solve_kkt_step(struct d_ocp_qp* qp, struct d_ocp_qp_sol* qp_sol, struct d_ocp_qp_ipm_arg* arg, struct d_ocp_qp_ipm_ws* ws);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif  // HPIPM_D_OCP_qp_KKT_H_
