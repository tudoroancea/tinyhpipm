#ifndef HPIPM_D_d_dense_qp_UTILS_H_
#define HPIPM_D_d_dense_qp_UTILS_H_

#include "hpipm/blas.h"

#include "hpipm/dense/d_dense_qp.h"
#include "hpipm/dense/d_dense_qp_dim.h"
#include "hpipm/dense/d_dense_qp_ipm.h"
#include "hpipm/dense/d_dense_qp_sol.h"


#ifdef __cplusplus
extern "C" {
#endif


//
void d_dense_qp_dim_print(struct d_dense_qp_dim* qp_dim);
//
void d_dense_qp_dim_codegen(char* file_name, char* mode, struct d_dense_qp_dim* qp_dim);
//
void d_dense_qp_print(struct d_dense_qp_dim* qp_dim, struct d_dense_qp* qp);
//
void d_dense_qp_codegen(char* file_name, char* mode, struct d_dense_qp_dim* qp_dim, struct d_dense_qp* qp);
//
void d_dense_qp_sol_print(struct d_dense_qp_dim* qp_dim, struct d_dense_qp_sol* dense_qp_sol);
//
void d_dense_qp_ipm_arg_codegen(char* file_name, char* mode, struct d_dense_qp_dim* qp_dim, struct d_dense_qp_ipm_arg* arg);
//
void d_dense_qp_res_print(struct d_dense_qp_dim* qp_dim, struct d_dense_qp_res* dense_qp_res);
//
void d_dense_qp_arg_print(struct d_dense_qp_dim* qp_dim, struct d_dense_qp_ipm_arg* qp_ipm_arg);


#ifdef __cplusplus
}  // #extern "C"
#endif


#endif  // HPIPM_D_d_dense_qp_UTILS_H_
