#ifndef HPIPM_D_CAST_QCQP_H_
#define HPIPM_D_CAST_QCQP_H_


#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/dense/d_dense_qcqp.h"
#include "hpipm/dense/d_dense_qcqp_sol.h"
#include "hpipm/ocp/d_ocp_qcqp.h"
#include "hpipm/ocp/d_ocp_qcqp_dim.h"
#include "hpipm/ocp/d_ocp_qcqp_sol.h"

#ifdef __cplusplus
extern "C" {
#endif


//
void d_cast_qcqp_compute_dim(struct d_ocp_qcqp_dim* ocp_dim, struct d_dense_qcqp_dim* dense_dim);
//
void d_cast_qcqp_cond(struct d_ocp_qcqp* ocp_qp, struct d_dense_qcqp* dense_qp);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif  // HPIPM_D_CAST_QCQP_H_
