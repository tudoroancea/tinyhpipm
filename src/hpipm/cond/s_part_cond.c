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



#define COND_B s_cond_b
#define COND_BABT s_cond_BAbt
#define COND_BAT s_cond_BAt
#define COND_D s_cond_d
#define COND_DCTD s_cond_DCtd
#define COND_DCT s_cond_DCt
#define COND_RQ s_cond_rq
#define COND_RSQ s_cond_RSQ
#define COND_RSQRQ s_cond_RSQrq
#define COND_QP_ARG s_cond_qp_arg
#define COND_QP_ARG_CREATE s_cond_qp_arg_create
#define COND_QP_ARG_MEMSIZE s_cond_qp_arg_memsize
#define COND_QP_ARG_SET_COMP_DUAL_SOL_EQ s_cond_qp_arg_set_comp_dual_sol_eq
#define COND_QP_ARG_SET_COMP_DUAL_SOL_INEQ s_cond_qp_arg_set_comp_dual_sol_ineq
#define COND_QP_ARG_SET_COMP_PRIM_SOL s_cond_qp_arg_set_comp_prim_sol
#define COND_QP_ARG_SET_COND_LAST_STAGE s_cond_qp_arg_set_cond_last_stage
#define COND_QP_ARG_SET_RIC_ALG s_cond_qp_arg_set_ric_alg
#define COND_QP_ARG_SET_DEFAULT s_cond_qp_arg_set_default
#define COND_QP_ARG_WS s_cond_qp_ws
#define COND_QP_WS_CREATE s_cond_qp_ws_create
#define COND_QP_WS_MEMSIZE s_cond_qp_ws_memsize
#define CREATE_STRVEC blasfeo_create_svec
#define DENSE_QP s_dense_qp
#define DENSE_QP_SOL s_dense_qp_sol
#define EXPAND_SOL s_expand_sol
#define GECP_LIBSTR blasfeo_sgecp
#define OCP_QP s_ocp_qp
#define OCP_QP_DIM s_ocp_qp_dim
#define OCP_QP_SOL s_ocp_qp_sol
#define PART_COND_QP_ARG s_part_cond_qp_arg
#define PART_COND_QP_WS s_part_cond_qp_ws
#define STRVEC blasfeo_svec
#define UPDATE_COND_BABT s_update_cond_BAbt
#define UPDATE_COND_DCTD s_update_cond_DCtd
#define UPDATE_COND_RSQRQ_N2NX3 s_update_cond_RSQrq_N2nx3
#define VECCP_LIBSTR blasfeo_sveccp

#define PART_COND_QP_COMPUTE_BLOCK_SIZE s_part_cond_qp_compute_block_size
#define PART_COND_QP_COMPUTE_DIM s_part_cond_qp_compute_dim
#define PART_COND_QP_ARG_MEMSIZE s_part_cond_qp_arg_memsize
#define PART_COND_QP_ARG_CREATE s_part_cond_qp_arg_create
#define PART_COND_QP_ARG_SET_DEFAULT s_part_cond_qp_arg_set_default
#define PART_COND_QP_ARG_SET_RIC_ALG s_part_cond_qp_arg_set_ric_alg
#define PART_COND_QP_ARG_SET_COMP_PRIM_SOL s_part_cond_qp_arg_set_comp_prim_sol
#define PART_COND_QP_ARG_SET_COMP_DUAL_SOL_EQ s_part_cond_qp_arg_set_comp_dual_sol_eq
#define PART_COND_QP_ARG_SET_COMP_DUAL_SOL_INEQ s_part_cond_qp_arg_set_comp_dual_sol_ineq
#define PART_COND_QP_WS_MEMSIZE s_part_cond_qp_ws_memsize
#define PART_COND_QP_WS_CREATE s_part_cond_qp_ws_create
#define PART_COND_QP_COND s_part_cond_qp_cond
#define PART_COND_QP_COND_LHS s_part_cond_qp_cond_lhs
#define PART_COND_QP_COND_RHS s_part_cond_qp_cond_rhs
#define PART_COND_QP_EXPAND_SOL s_part_cond_qp_expand_sol
#define PART_COND_QP_UPDATE s_part_cond_qp_update



#include "x_part_cond.c"

