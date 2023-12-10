#ifndef HPIPM_M_DENSE_QP_DIM_H_
#define HPIPM_M_DENSE_QP_DIM_H_

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_target.h"

#include "hpipm/dense/d_dense_qp_dim.h"
#include "hpipm/dense/s_dense_qp_dim.h"


#ifdef __cplusplus
extern "C" {
#endif


void cvt_d2s_dense_qp_dim(struct d_dense_qp_dim* qpd, struct s_dense_qp_dim* qps);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif  // HPIPM_M_DENSE_QP_DIM_H_
