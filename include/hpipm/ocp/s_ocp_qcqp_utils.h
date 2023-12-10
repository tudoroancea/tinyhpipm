#ifndef HPIPM_S_OCP_QCQP_UTILS_H_
#define HPIPM_S_OCP_QCQP_UTILS_H_

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_target.h"

#include "hpipm/ocp/s_ocp_qcqp_dim.h"
#include "hpipm/ocp/s_ocp_qcqp_ipm.h"
#include "hpipm/ocp/s_ocp_qcqp_sol.h"
#include "hpipm/ocp/s_ocp_qp.h"


#ifdef __cplusplus
extern "C" {
#endif


//
void s_ocp_qcqp_dim_print(struct s_ocp_qcqp_dim* qcqp_dim);
//
void s_ocp_qcqp_dim_codegen(char* file_name, char* mode, struct s_ocp_qcqp_dim* qcqp_dim);
//
void s_ocp_qcqp_print(struct s_ocp_qcqp_dim* qcqp_dim, struct s_ocp_qcqp* qp);
//
void s_ocp_qcqp_codegen(char* file_name, char* mode, struct s_ocp_qcqp_dim* qcqp_dim, struct s_ocp_qcqp* qp);
//
void s_ocp_qcqp_sol_print(struct s_ocp_qcqp_dim* qcqp_dim, struct s_ocp_qcqp_sol* ocp_qcqp_sol);
//
void s_ocp_qcqp_ipm_arg_codegen(char* file_name, char* mode, struct s_ocp_qcqp_dim* qcqp_dim, struct s_ocp_qcqp_ipm_arg* arg);
//
void s_ocp_qcqp_res_print(struct s_ocp_qcqp_dim* qcqp_dim, struct s_ocp_qcqp_res* ocp_qcqp_res);


#ifdef __cplusplus
}  // #extern "C"
#endif


#endif  // HPIPM_S_OCP_QCQP_UTILS_H_
