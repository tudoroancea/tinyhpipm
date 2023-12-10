#ifndef HPIPM_D_DENSE_QCQP_UTILS_H_
#define HPIPM_D_DENSE_QCQP_UTILS_H_

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_target.h"

#include "hpipm/dense/d_dense_qcqp.h"
#include "hpipm/dense/d_dense_qcqp_dim.h"
#include "hpipm/dense/d_dense_qcqp_ipm.h"
#include "hpipm/dense/d_dense_qcqp_sol.h"


#ifdef __cplusplus
extern "C" {
#endif


//
void d_dense_qcqp_dim_print(struct d_dense_qcqp_dim* qp_dim);
//
// void d_dense_qcqp_dim_codegen(char *file_name, char *mode, struct d_dense_qcqp_dim *qp_dim);
//
void d_dense_qcqp_print(struct d_dense_qcqp_dim* qp_dim, struct d_dense_qcqp* qp);
//
// void d_dense_qcqp_codegen(char *file_name, char *mode, struct d_dense_qcqp_dim *qp_dim, struct d_dense_qcqp *qp);
//
void d_dense_qcqp_sol_print(struct d_dense_qcqp_dim* qp_dim, struct d_dense_qcqp_sol* dense_qcqp_sol);
//
// void d_dense_qcqp_ipm_arg_codegen(char *file_name, char *mode, struct d_dense_qcqp_dim *qp_dim, struct d_dense_qcqp_ipm_arg *arg);
//
void d_dense_qcqp_res_print(struct d_dense_qcqp_dim* qp_dim, struct d_dense_qcqp_res* dense_qcqp_res);


#ifdef __cplusplus
}  // #extern "C"
#endif


#endif  // HPIPM_D_DENSE_QCQP_UTILS_H_
