#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_s_aux.h"
#include "blasfeo/blasfeo_s_blas.h"
#include "blasfeo/blasfeo_target.h"

#include "hpipm/auxiliary/aux_mem.h"
#include "hpipm/cond/s_cond.h"
#include "hpipm/cond/s_cond_aux.h"
#include "hpipm/dense/s_dense_qp.h"
#include "hpipm/dense/s_dense_qp_sol.h"
#include "hpipm/ocp/s_ocp_qp.h"
#include "hpipm/ocp/s_ocp_qp_dim.h"
#include "hpipm/ocp/s_ocp_qp_sol.h"



#define COND_DCT s_cond_DCt
#define COND_DCTD s_cond_DCtd
#define COND_D s_cond_d
#define COND_B s_cond_b
#define COND_BABT s_cond_BAbt
#define COND_BAT s_cond_BAt
#define COND_RQ s_cond_rq
#define COND_RSQ s_cond_RSQ
#define COND_RSQRQ s_cond_RSQrq
#define COND_QP_ARG s_cond_qp_arg
#define COND_QP_WS s_cond_qp_ws
#define CREATE_STRMAT blasfeo_create_smat
#define CREATE_STRVEC blasfeo_create_svec
#define DENSE_QP s_dense_qp
#define DENSE_QP_DIM s_dense_qp_dim
#define DENSE_QP_SOL s_dense_qp_sol
#define EXPAND_SOL s_expand_sol
#define EXPAND_PRIMAL_SOL s_expand_primal_sol
#define OCP_QP s_ocp_qp
#define OCP_QP_DIM s_ocp_qp_dim
#define OCP_QP_SOL s_ocp_qp_sol
#define REAL float
#define SIZE_STRMAT blasfeo_memsize_smat
#define SIZE_STRVEC blasfeo_memsize_svec
#define STRMAT blasfeo_smat
#define STRVEC blasfeo_svec
#define UPDATE_COND_DCTD s_update_cond_DCtd
#define UPDATE_COND_BABT s_update_cond_BAbt
#define UPDATE_COND_RSQRQ_N2NX3 s_update_cond_RSQrq_N2nx3

#define COND_QP_COMPUTE_DIM s_cond_qp_compute_dim
#define COND_QP_ARG_MEMSIZE s_cond_qp_arg_memsize
#define COND_QP_ARG_CREATE s_cond_qp_arg_create
#define COND_QP_ARG_SET_DEFAULT s_cond_qp_arg_set_default
#define COND_QP_ARG_SET_COND_ALG s_cond_qp_arg_set_cond_alg
#define COND_QP_ARG_SET_RIC_ALG s_cond_qp_arg_set_ric_alg
#define COND_QP_ARG_SET_COND_LAST_STAGE s_cond_qp_arg_set_cond_last_stage
#define COND_QP_ARG_SET_COMP_PRIM_SOL s_cond_qp_arg_set_comp_prim_sol
#define COND_QP_ARG_SET_COMP_DUAL_SOL_EQ s_cond_qp_arg_set_comp_dual_sol_eq
#define COND_QP_ARG_SET_COMP_DUAL_SOL_INEQ s_cond_qp_arg_set_comp_dual_sol_ineq
#define COND_QP_WS_MEMSIZE s_cond_qp_ws_memsize
#define COND_QP_WS_CREATE s_cond_qp_ws_create
#define COND_QP_COND s_cond_qp_cond
#define COND_QP_COND_LHS s_cond_qp_cond_lhs
#define COND_QP_COND_RHS s_cond_qp_cond_rhs
#define COND_QP_EXPAND_SOL s_cond_qp_expans_sol
#define COND_QP_EXPAND_PRIMAL_SOL s_cond_qp_expans_primal_sol
#define COND_QP_UPDATE s_cond_qp_update



#include "x_cond.c"

