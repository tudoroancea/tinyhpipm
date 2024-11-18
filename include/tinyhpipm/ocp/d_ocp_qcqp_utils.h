#ifndef HPIPM_D_OCP_qcqp_UTILS_H_
#define HPIPM_D_OCP_qcqp_UTILS_H_

#include "tinyhpipm/blas.h"
#include "tinyhpipm/ocp/d_ocp_qcqp_dim.h"
#include "tinyhpipm/ocp/d_ocp_qcqp_ipm.h"
#include "tinyhpipm/ocp/d_ocp_qcqp_sol.h"
#include "tinyhpipm/ocp/d_ocp_qp.h"


#ifdef __cplusplus
extern "C" {
#endif


//
void d_ocp_qcqp_dim_print(struct d_ocp_qcqp_dim* qcqp_dim);
//
void d_ocp_qcqp_dim_codegen(char* file_name, char* mode, struct d_ocp_qcqp_dim* qcqp_dim);
//
void d_ocp_qcqp_print(struct d_ocp_qcqp_dim* qcqp_dim, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_codegen(char* file_name, char* mode, struct d_ocp_qcqp_dim* qcqp_dim, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_sol_print(struct d_ocp_qcqp_dim* qcqp_dim, struct d_ocp_qcqp_sol* ocp_qcqp_sol);
//
void d_ocp_qcqp_ipm_arg_codegen(char* file_name, char* mode, struct d_ocp_qcqp_dim* qcqp_dim, struct d_ocp_qcqp_ipm_arg* arg);
//
void d_ocp_qcqp_res_print(struct d_ocp_qcqp_dim* qcqp_dim, struct d_ocp_qcqp_res* ocp_qcqp_res);


#ifdef __cplusplus
}  // #extern "C"
#endif


#endif  // HPIPM_D_OCP_qcqp_UTILS_H_
