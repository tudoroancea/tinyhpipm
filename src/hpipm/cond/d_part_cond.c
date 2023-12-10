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


#define COND_B d_cond_b
#define COND_BABT d_cond_BAbt
#define COND_BAT d_cond_BAt
#define COND_D d_cond_d
#define COND_DCT d_cond_DCt
#define COND_DCTD d_cond_DCtd
#define COND_RQ d_cond_rq
#define COND_RSQ d_cond_RSQ
#define COND_RSQRQ d_cond_RSQrq
#define COND_QP_ARG d_cond_qp_arg
#define COND_QP_ARG_CREATE d_cond_qp_arg_create
#define COND_QP_ARG_MEMSIZE d_cond_qp_arg_memsize
#define COND_QP_ARG_SET_COMP_DUAL_SOL_EQ d_cond_qp_arg_set_comp_dual_sol_eq
#define COND_QP_ARG_SET_COMP_DUAL_SOL_INEQ d_cond_qp_arg_set_comp_dual_sol_ineq
#define COND_QP_ARG_SET_COMP_PRIM_SOL d_cond_qp_arg_set_comp_prim_sol
#define COND_QP_ARG_SET_COND_LAST_STAGE d_cond_qp_arg_set_cond_last_stage
#define COND_QP_ARG_SET_DEFAULT d_cond_qp_arg_set_default
#define COND_QP_ARG_SET_RIC_ALG d_cond_qp_arg_set_ric_alg
#define COND_QP_ARG_WS d_cond_qp_ws
#define COND_QP_WS_CREATE d_cond_qp_ws_create
#define COND_QP_WS_MEMSIZE d_cond_qp_ws_memsize
#define CREATE_STRVEC blasfeo_create_dvec
#define DENSE_QP d_dense_qp
#define DENSE_QP_SOL d_dense_qp_sol
#define EXPAND_SOL d_expand_sol
#define GECP_LIBSTR blasfeo_dgecp
#define OCP_QP d_ocp_qp
#define OCP_QP_DIM d_ocp_qp_dim
#define OCP_QP_SOL d_ocp_qp_sol
#define PART_COND_QP_ARG d_part_cond_qp_arg
#define PART_COND_QP_WS d_part_cond_qp_ws
#define STRVEC blasfeo_dvec
#define UPDATE_COND_BABT d_update_cond_BAbt
#define UPDATE_COND_DCTD d_update_cond_DCtd
#define UPDATE_COND_RSQRQ_N2NX3 d_update_cond_RSQrq_N2nx3
#define VECCP_LIBSTR blasfeo_dveccp

#define PART_COND_QP_COMPUTE_BLOCK_SIZE d_part_cond_qp_compute_block_size
#define PART_COND_QP_COMPUTE_DIM d_part_cond_qp_compute_dim
#define PART_COND_QP_ARG_MEMSIZE d_part_cond_qp_arg_memsize
#define PART_COND_QP_ARG_CREATE d_part_cond_qp_arg_create
#define PART_COND_QP_ARG_SET_DEFAULT d_part_cond_qp_arg_set_default
#define PART_COND_QP_ARG_SET_RIC_ALG d_part_cond_qp_arg_set_ric_alg
#define PART_COND_QP_ARG_SET_COMP_PRIM_SOL d_part_cond_qp_arg_set_comp_prim_sol
#define PART_COND_QP_ARG_SET_COMP_DUAL_SOL_EQ d_part_cond_qp_arg_set_comp_dual_sol_eq
#define PART_COND_QP_ARG_SET_COMP_DUAL_SOL_INEQ d_part_cond_qp_arg_set_comp_dual_sol_ineq
#define PART_COND_QP_WS_MEMSIZE d_part_cond_qp_ws_memsize
#define PART_COND_QP_WS_CREATE d_part_cond_qp_ws_create
#define PART_COND_QP_COND d_part_cond_qp_cond
#define PART_COND_QP_COND_LHS d_part_cond_qp_cond_lhs
#define PART_COND_QP_COND_RHS d_part_cond_qp_cond_rhs
#define PART_COND_QP_EXPAND_SOL d_part_cond_qp_expand_sol
#define PART_COND_QP_UPDATE d_part_cond_qp_update


#include "x_part_cond.c"
