#include <stdlib.h>
#include <stdio.h>

#include <blasfeo/blasfeo_common.h>
#include <blasfeo/blasfeo_s_aux.h>
#include <blasfeo/blasfeo_s_blas.h>
#include <blasfeo/blasfeo_target.h>

#include "hpipm/dense/s_dense_qcqp.h"
#include "hpipm/dense/s_dense_qcqp_dim.h"
#include "hpipm/dense/s_dense_qcqp_sol.h"
#include "hpipm/ocp/s_ocp_qcqp.h"
#include "hpipm/ocp/s_ocp_qcqp_dim.h"
#include "hpipm/ocp/s_ocp_qcqp_sol.h"



#define DENSE_QCQP s_dense_qcqp
#define DENSE_QCQP_DIM s_dense_qcqp_dim
#define DENSE_QCQP_DIM_SET_NV s_dense_qcqp_dim_set_nv
#define DENSE_QCQP_DIM_SET_NE s_dense_qcqp_dim_set_ne
#define DENSE_QCQP_DIM_SET_NB s_dense_qcqp_dim_set_nb
#define DENSE_QCQP_DIM_SET_NG s_dense_qcqp_dim_set_ng
#define DENSE_QCQP_DIM_SET_NQ s_dense_qcqp_dim_set_nq
#define DENSE_QCQP_DIM_SET_NS s_dense_qcqp_dim_set_ns
#define DENSE_QCQP_DIM_SET_NSB s_dense_qcqp_dim_set_nsb
#define DENSE_QCQP_DIM_SET_NSG s_dense_qcqp_dim_set_nsg
#define DENSE_QCQP_DIM_SET_NSQ s_dense_qcqp_dim_set_nsq
#define DENSE_QCQP_SOL s_dense_qcqp_sol
#define DENSE_QP s_dense_qp
#define DENSE_QP_DIM s_dense_qp_dim
#define DENSE_QP_SOL s_dense_qp_sol
#define DIARE blasfeo_sdiare
#define GECP blasfeo_sgecp
#define GETR blasfeo_sgetr
#define OCP_QCQP s_ocp_qcqp
#define OCP_QCQP_DIM s_ocp_qcqp_dim
#define OCP_QCQP_SOL s_ocp_qcqp_sol
#define OCP_QP s_ocp_qp
#define OCP_QP_DIM s_ocp_qp_dim
#define OCP_QP_SOL s_ocp_qp_sol
#define SIZE_STRMAT blasfeo_memsize_smat
#define SIZE_STRVEC blasfeo_memsize_svec
#define STRMAT blasfeo_smat
#define STRVEC blasfeo_svec
#define VECCP blasfeo_sveccp

#define CAST_QCQP_COMPUTE_DIM s_cast_qcqp_compute_dim
#define CAST_QCQP_COND s_cast_qcqp_cond



#include "x_cast_qcqp.c"



