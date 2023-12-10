#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_d_aux.h"
#include "blasfeo/blasfeo_d_blas.h"
#include "blasfeo/blasfeo_target.h"

#include "hpipm/cond/d_cond.h"
#include "hpipm/cond/d_cond_aux.h"
#include "hpipm/cond/d_part_cond.h"
#include "hpipm/dense/d_dense_qp.h"
#include "hpipm/dense/d_dense_qp_sol.h"
#include "hpipm/ocp/d_ocp_qp.h"
#include "hpipm/ocp/d_ocp_qp_dim.h"
#include "hpipm/ocp/d_ocp_qp_sol.h"


#define SINGLE_PRECISION


#define COLAD blasfeo_scolad
#define COND_DCT s_cond_DCt
#define COND_DCTD s_cond_DCtd
#define COND_D s_cond_d
#define COND_B s_cond_b
#define COND_BABT s_cond_BAbt
#define COND_BAT s_cond_BAt
#define COND_RQ s_cond_rq
#define COND_RSQ s_cond_RSQ
#define COND_RSQRQ s_cond_RSQrq
#define COND_QCQP_ARG s_cond_qcqp_arg
#define COND_QCQP_WS s_cond_qcqp_ws
#define COND_QP_ARG s_cond_qp_arg
#define COND_QP_ARG_MEMSIZE s_cond_qp_arg_memsize
#define COND_QP_ARG_CREATE s_cond_qp_arg_create
#define COND_QP_WS s_cond_qp_ws
#define COND_QP_WS_MEMSIZE s_cond_qp_ws_memsize
#define COND_QP_WS_CREATE s_cond_qp_ws_create
#define CREATE_STRMAT blasfeo_create_smat
#define CREATE_STRVEC blasfeo_create_svec
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
#define DOT blasfeo_sdot
#define EXPAND_SOL s_expand_sol
#define EXPAND_PRIMAL_SOL s_expand_primal_sol
#define GEAD blasfeo_sgead
#define GECP blasfeo_sgecp
#define GEMM_NN blasfeo_sgemm_nn
#define GEMV_N blasfeo_sgemv_n
#define GESE blasfeo_sgese
#define OCP_QCQP s_ocp_qcqp
#define OCP_QCQP_DIM s_ocp_qcqp_dim
#define OCP_QCQP_SOL s_ocp_qcqp_sol
#define OCP_QP s_ocp_qp
#define OCP_QP_DIM s_ocp_qp_dim
#define OCP_QP_SOL s_ocp_qp_sol
#define REAL float
#define ROWEX blasfeo_srowex
#define SIZE_STRMAT blasfeo_memsize_smat
#define SIZE_STRVEC blasfeo_memsize_svec
#define STRMAT blasfeo_smat
#define STRVEC blasfeo_svec
#define SYMV_L blasfeo_ssymv_l
#define SYRK_LN blasfeo_ssyrk_ln
#define SYRK_LN_MN blasfeo_ssyrk_ln_mn
#define TRCP_L blasfeo_strcp_l
#define TRTR_L blasfeo_strtr_l

#define COND_QCQP_COMPUTE_DIM s_cond_qcqp_compute_dim
#define COND_QCQP_ARG_MEMSIZE s_cond_qcqp_arg_memsize
#define COND_QCQP_ARG_CREATE s_cond_qcqp_arg_create
#define COND_QCQP_ARG_SET_DEFAULT s_cond_qcqp_arg_set_default
#define COND_QCQP_ARG_SET_RIC_ALG s_cond_qcqp_arg_set_ric_alg
#define COND_QCQP_ARG_SET_COND_LAST_STAGE s_cond_qcqp_arg_set_cond_last_stage
#define COND_QCQP_WS_MEMSIZE s_cond_qcqp_ws_memsize
#define COND_QCQP_WS_CREATE s_cond_qcqp_ws_create
#define COND_QCQP_QC s_cond_qcqp_qc
#define COND_QCQP_QC_LHS s_cond_qcqp_qc_lhs
#define COND_QCQP_QC_RHS s_cond_qcqp_qc_rhs
#define COND_QCQP_COND s_cond_qcqp_cond
#define COND_QCQP_COND_LHS s_cond_qcqp_cond_lhs
#define COND_QCQP_COND_RHS s_cond_qcqp_cond_rhs
#define COND_QCQP_EXPAND_SOL s_cond_qcqp_expand_sol


#include "x_cond_qcqp.c"
