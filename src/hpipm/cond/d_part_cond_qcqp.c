#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_d_aux.h"
#include "blasfeo/blasfeo_d_blas.h"
#include "blasfeo/blasfeo_target.h"

#include "hpipm/cond/d_cond.h"
#include "hpipm/cond/d_cond_aux.h"
#include "hpipm/cond/d_cond_qcqp.h"
#include "hpipm/cond/d_part_cond_qcqp.h"
#include "hpipm/dense/d_dense_qcqp.h"
#include "hpipm/dense/d_dense_qcqp_sol.h"
#include "hpipm/dense/d_dense_qp.h"
#include "hpipm/dense/d_dense_qp_sol.h"
#include "hpipm/ocp/d_ocp_qcqp.h"
#include "hpipm/ocp/d_ocp_qcqp_dim.h"
#include "hpipm/ocp/d_ocp_qcqp_sol.h"
#include "hpipm/ocp/d_ocp_qp.h"
#include "hpipm/ocp/d_ocp_qp_dim.h"
#include "hpipm/ocp/d_ocp_qp_sol.h"


#define COLAD blasfeo_dcolad
#define COND_B d_cond_b
#define COND_BABT d_cond_BAbt
#define COND_BAT d_cond_BAt
#define COND_D d_cond_d
#define COND_DCT d_cond_DCt
#define COND_DCTD d_cond_DCtd
#define COND_RQ d_cond_rq
#define COND_RSQ d_cond_RSQ
#define COND_RSQRQ d_cond_RSQrq
#define COND_QCQP_ARG d_cond_qcqp_arg
#define COND_QCQP_ARG_CREATE d_cond_qcqp_arg_create
#define COND_QCQP_ARG_MEMSIZE d_cond_qcqp_arg_memsize
#define COND_QCQP_ARG_SET_COND_LAST_STAGE d_cond_qcqp_arg_set_cond_last_stage
#define COND_QCQP_ARG_SET_DEFAULT d_cond_qcqp_arg_set_default
#define COND_QCQP_ARG_SET_RIC_ALG d_cond_qcqp_arg_set_ric_alg
#define COND_QCQP_ARG_WS d_cond_qcqp_ws
#define COND_QCQP_QC d_cond_qcqp_qc
#define COND_QCQP_QC_LHS d_cond_qcqp_qc_lhs
#define COND_QCQP_QC_RHS d_cond_qcqp_qc_rhs
#define PART_COND_QCQP_ARG d_part_cond_qcqp_arg
#define PART_COND_QCQP_WS d_part_cond_qcqp_ws
#define COND_QCQP_WS_CREATE d_cond_qcqp_ws_create
#define COND_QCQP_WS_MEMSIZE d_cond_qcqp_ws_memsize
#define CREATE_STRVEC blasfeo_create_dvec
#define DENSE_QCQP d_dense_qcqp
#define DENSE_QCQP_SOL d_dense_qcqp_sol
#define DENSE_QP d_dense_qp
#define DENSE_QP_SOL d_dense_qp_sol
#define EXPAND_SOL d_expand_sol
#define GECP blasfeo_dgecp
#define OCP_QCQP d_ocp_qcqp
#define OCP_QCQP_DIM d_ocp_qcqp_dim
#define OCP_QCQP_DIM_SET_NX d_ocp_qcqp_dim_set_nx
#define OCP_QCQP_DIM_SET_NU d_ocp_qcqp_dim_set_nu
#define OCP_QCQP_DIM_SET_NBX d_ocp_qcqp_dim_set_nbx
#define OCP_QCQP_DIM_SET_NBU d_ocp_qcqp_dim_set_nbu
#define OCP_QCQP_DIM_SET_NG d_ocp_qcqp_dim_set_ng
#define OCP_QCQP_DIM_SET_NQ d_ocp_qcqp_dim_set_nq
#define OCP_QCQP_DIM_SET_NSBX d_ocp_qcqp_dim_set_nsbx
#define OCP_QCQP_DIM_SET_NSBU d_ocp_qcqp_dim_set_nsbu
#define OCP_QCQP_DIM_SET_NSG d_ocp_qcqp_dim_set_nsg
#define OCP_QCQP_DIM_SET_NSQ d_ocp_qcqp_dim_set_nsq
#define OCP_QCQP_DIM_SET_NS d_ocp_qcqp_dim_set_ns
#define OCP_QCQP_SOL d_ocp_qcqp_sol
#define OCP_QP d_ocp_qp
#define OCP_QP_DIM d_ocp_qp_dim
#define OCP_QP_SOL d_ocp_qp_sol
#define STRVEC blasfeo_dvec
#define SYMV_L blasfeo_dsymv_l
#define UPDATE_COND_BABT d_update_cond_BAbt
#define UPDATE_COND_DCTD d_update_cond_DCtd
#define UPDATE_COND_RSQRQ_N2NX3 d_update_cond_RSQrq_N2nx3
#define VECCP blasfeo_dveccp

#define PART_COND_QCQP_COMPUTE_BLOCK_SIZE d_part_cond_qcqp_compute_block_size
#define PART_COND_QCQP_COMPUTE_DIM d_part_cond_qcqp_compute_dim
#define PART_COND_QCQP_ARG_MEMSIZE d_part_cond_qcqp_arg_memsize
#define PART_COND_QCQP_ARG_CREATE d_part_cond_qcqp_arg_create
#define PART_COND_QCQP_ARG_SET_DEFAULT d_part_cond_qcqp_arg_set_default
#define PART_COND_QCQP_ARG_SET_RIC_ALG d_part_cond_qcqp_arg_set_ric_alg
#define PART_COND_QCQP_WS_MEMSIZE d_part_cond_qcqp_ws_memsize
#define PART_COND_QCQP_WS_CREATE d_part_cond_qcqp_ws_create
#define PART_COND_QCQP_COND d_part_cond_qcqp_cond
#define PART_COND_QCQP_COND_LHS d_part_cond_qcqp_cond_lhs
#define PART_COND_QCQP_COND_RHS d_part_cond_qcqp_cond_rhs
#define PART_COND_QCQP_EXPAND_SOL d_part_cond_qcqp_expand_sol


#include "x_part_cond_qcqp.c"