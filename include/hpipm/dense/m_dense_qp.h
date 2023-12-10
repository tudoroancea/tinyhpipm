#ifndef HPIPM_M_DENSE_QP_H_
#define HPIPM_M_DENSE_QP_H_

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_target.h"

#include "hpipm/dense/d_dense_qp.h"
#include "hpipm/dense/s_dense_qp.h"

#ifdef __cplusplus
extern "C" {
#endif

void cvt_d2s_dense_qp(struct d_dense_qp* qpd, struct s_dense_qp* qps);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // HPIPM_M_DENSE_QP_H_
