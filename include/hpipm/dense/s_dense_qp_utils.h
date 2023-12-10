#ifndef HPIPM_S_DENSE_QP_UTILS_H_
#define HPIPM_S_DENSE_QP_UTILS_H_

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_target.h"

#include "hpipm/dense/s_dense_qp.h"
#include "hpipm/dense/s_dense_qp_dim.h"
#include "hpipm/dense/s_dense_qp_ipm.h"
#include "hpipm/dense/s_dense_qp_sol.h"


#ifdef __cplusplus
extern "C" {
#endif


//
void s_dense_qp_dim_print(struct s_dense_qp_dim* qp_dim);
//
void s_dense_qp_dim_codegen(char* file_name, char* mode, struct s_dense_qp_dim* qp_dim);
//
void s_dense_qp_print(struct s_dense_qp_dim* qp_dim, struct s_dense_qp* qp);
//
void s_dense_qp_codegen(char* file_name, char* mode, struct s_dense_qp_dim* qp_dim, struct s_dense_qp* qp);
//
void s_dense_qp_sol_print(struct s_dense_qp_dim* qp_dim, struct s_dense_qp_sol* dense_qp_sol);
//
void s_dense_qp_ipm_arg_codegen(char* file_name, char* mode, struct s_dense_qp_dim* qp_dim, struct s_dense_qp_ipm_arg* arg);
//
void s_dense_qp_res_print(struct s_dense_qp_dim* qp_dim, struct s_dense_qp_res* dense_qp_res);
//
void s_dense_qp_arg_print(struct s_dense_qp_dim* qp_dim, struct s_dense_qp_ipm_arg* qp_ipm_arg);


#ifdef __cplusplus
}  // #extern "C"
#endif


#endif  // HPIPM_S_DENSE_QP_UTILS_H_
