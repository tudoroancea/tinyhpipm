#include <stdio.h>
#include <stdlib.h>

#include <blasfeo/blasfeo_common.h>
#include <blasfeo/blasfeo_d_aux.h>
#include <blasfeo/blasfeo_d_blas.h>
#include <blasfeo/blasfeo_target.h>

#include "hpipm/dense/d_dense_qcqp.h"
#include "hpipm/dense/d_dense_qcqp_dim.h"
#include "hpipm/dense/d_dense_qcqp_sol.h"
#include "hpipm/ocp/d_ocp_qcqp.h"
#include "hpipm/ocp/d_ocp_qcqp_dim.h"
#include "hpipm/ocp/d_ocp_qcqp_sol.h"


#define DENSE_QCQP d_dense_qcqp
#define DENSE_QCQP_DIM d_dense_qcqp_dim
#define DENSE_QCQP_DIM_SET_NV d_dense_qcqp_dim_set_nv
#define DENSE_QCQP_DIM_SET_NE d_dense_qcqp_dim_set_ne
#define DENSE_QCQP_DIM_SET_NB d_dense_qcqp_dim_set_nb
#define DENSE_QCQP_DIM_SET_NG d_dense_qcqp_dim_set_ng
#define DENSE_QCQP_DIM_SET_NQ d_dense_qcqp_dim_set_nq
#define DENSE_QCQP_DIM_SET_NS d_dense_qcqp_dim_set_ns
#define DENSE_QCQP_DIM_SET_NSB d_dense_qcqp_dim_set_nsb
#define DENSE_QCQP_DIM_SET_NSG d_dense_qcqp_dim_set_nsg
#define DENSE_QCQP_DIM_SET_NSQ d_dense_qcqp_dim_set_nsq
#define DENSE_QCQP_SOL d_dense_qcqp_sol
#define DENSE_QP d_dense_qp
#define DENSE_QP_DIM d_dense_qp_dim
#define DENSE_QP_SOL d_dense_qp_sol
#define DIARE blasfeo_ddiare
#define GECP blasfeo_dgecp
#define GETR blasfeo_dgetr
#define OCP_QCQP d_ocp_qcqp
#define OCP_QCQP_DIM d_ocp_qcqp_dim
#define OCP_QCQP_SOL d_ocp_qcqp_sol
#define OCP_QP d_ocp_qp
#define OCP_QP_DIM d_ocp_qp_dim
#define OCP_QP_SOL d_ocp_qp_sol
#define SIZE_STRMAT blasfeo_memsize_dmat
#define SIZE_STRVEC blasfeo_memsize_dvec
#define STRMAT blasfeo_dmat
#define STRVEC blasfeo_dvec
#define VECCP blasfeo_dveccp

#define CAST_QCQP_COMPUTE_DIM d_cast_qcqp_compute_dim
#define CAST_QCQP_COND d_cast_qcqp_cond


#include "x_cast_qcqp.c"
