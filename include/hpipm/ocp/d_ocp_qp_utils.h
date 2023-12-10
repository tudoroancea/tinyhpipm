#ifndef HPIPM_D_OCP_QP_UTILS_H_
#define HPIPM_D_OCP_QP_UTILS_H_

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/ocp/d_ocp_qp.h"
#include "hpipm/ocp/d_ocp_qp_dim.h"
#include "hpipm/ocp/d_ocp_qp_ipm.h"
#include "hpipm/ocp/d_ocp_qp_sol.h"


#ifdef __cplusplus
extern "C" {
#endif


//
void d_ocp_qp_dim_print(struct d_ocp_qp_dim* qp_dim);
//
void d_ocp_qp_dim_codegen(char* file_name, char* mode, struct d_ocp_qp_dim* qp_dim);
//
void d_ocp_qp_print(struct d_ocp_qp_dim* qp_dim, struct d_ocp_qp* qp);
//
void d_ocp_qp_codegen(char* file_name, char* mode, struct d_ocp_qp_dim* qp_dim, struct d_ocp_qp* qp);
//
void d_ocp_qp_sol_print(struct d_ocp_qp_dim* qp_dim, struct d_ocp_qp_sol* ocp_qp_sol);
//
void d_ocp_qp_ipm_arg_print(struct d_ocp_qp_dim* qp_dim, struct d_ocp_qp_ipm_arg* arg);
//
void d_ocp_qp_ipm_arg_codegen(char* file_name, char* mode, struct d_ocp_qp_dim* qp_dim, struct d_ocp_qp_ipm_arg* arg);
//
void d_ocp_qp_res_print(struct d_ocp_qp_dim* qp_dim, struct d_ocp_qp_res* ocp_qp_res);


#ifdef __cplusplus
}  // #extern "C"
#endif


#endif  // HPIPM_D_OCP_QP_UTILS_H_
