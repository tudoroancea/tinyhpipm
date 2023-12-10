#ifndef HPIPM_S_DENSE_QCQP_UTILS_H_
#define HPIPM_S_DENSE_QCQP_UTILS_H_

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_target.h"

#include "hpipm/dense/s_dense_qcqp.h"
#include "hpipm/dense/s_dense_qcqp_dim.h"
#include "hpipm/dense/s_dense_qcqp_sol.h"
// #include "hpipm/dense/s_dense_qcqp_ipm.h"


#ifdef __cplusplus
extern "C" {
#endif


//
void s_dense_qcqp_dim_print(struct s_dense_qcqp_dim* qp_dim);
//
// void s_dense_qcqp_dim_codegen(char *file_name, char *mode, struct s_dense_qcqp_dim *qp_dim);
//
void s_dense_qcqp_print(struct s_dense_qcqp_dim* qp_dim, struct s_dense_qcqp* qp);
//
// void s_dense_qcqp_codegen(char *file_name, char *mode, struct s_dense_qcqp_dim *qp_dim, struct s_dense_qcqp *qp);
//
void s_dense_qcqp_sol_print(struct s_dense_qcqp_dim* qp_dim, struct s_dense_qcqp_sol* dense_qcqp_sol);
//
// void s_dense_qcqp_ipm_arg_codegen(char *file_name, char *mode, struct s_dense_qcqp_dim *qp_dim, struct s_dense_qcqp_ipm_arg *arg);
//
void s_dense_qcqp_res_print(struct s_dense_qcqp_dim* qp_dim, struct s_dense_qcqp_res* dense_qcqp_res);


#ifdef __cplusplus
}  // #extern "C"
#endif


#endif  // HPIPM_S_DENSE_QCQP_UTILS_H_
