#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_s_aux.h"
#include "blasfeo/blasfeo_s_blas.h"
#include "blasfeo/blasfeo_target.h"

#include "hpipm/auxiliary/aux_mem.h"
#include "hpipm/cond/s_cond.h"
#include "hpipm/cond/s_cons_aux.h"
#include "hpipm/dense/s_dense_qp.h"
#include "hpipm/dense/s_dense_qp_sol.h"
#include "hpipm/ocp/s_ocp_qp.h"
#include "hpipm/ocp/s_ocp_qp_dim.h"
#include "hpipm/ocp/s_ocp_qp_sol.h"



#define CONs_DCT s_cons_DCt
#define CONs_DCTD s_cons_DCtd
#define CONs_D s_cons_d
#define CONs_B s_cons_b
#define CONs_BABT s_cons_BAbt
#define CONs_BAT s_cons_BAt
#define CONs_RQ s_cons_rq
#define CONs_RSQ s_cons_RSQ
#define CONs_RSQRQ s_cons_RSQrq
#define CONs_QP_ARG s_cons_qp_arg
#define CONs_QP_WS s_cons_qp_ws
#define CREATE_STRMAT blasfeo_create_smat
#define CREATE_STRVEC blasfeo_create_svec
#define DENSE_QP s_dense_qp
#define DENSE_QP_DIM s_dense_qp_dim
#define DENSE_QP_SOL s_dense_qp_sol
#define EXPANs_SOL s_expans_sol
#define EXPANs_PRIMAL_SOL s_expans_primal_sol
#define OCP_QP s_ocp_qp
#define OCP_QP_DIM s_ocp_qp_dim
#define OCP_QP_SOL s_ocp_qp_sol
#define REAL float
#define SIZE_STRMAT blasfeo_memsize_smat
#define SIZE_STRVEC blasfeo_memsize_svec
#define STRMAT blasfeo_smat
#define STRVEC blasfeo_svec
#define UPDATE_CONs_DCTD s_update_cons_DCtd
#define UPDATE_CONs_BABT s_update_cons_BAbt
#define UPDATE_CONs_RSQRQ_N2NX3 s_update_cons_RSQrq_N2nx3

#define CONs_QP_COMPUTE_DIM s_cons_qp_compute_dim
#define CONs_QP_ARG_MEMSIZE s_cons_qp_arg_memsize
#define CONs_QP_ARG_CREATE s_cons_qp_arg_create
#define CONs_QP_ARG_SET_DEFAULT s_cons_qp_arg_set_default
#define CONs_QP_ARG_SET_CONs_ALG s_cons_qp_arg_set_cons_alg
#define CONs_QP_ARG_SET_RIC_ALG s_cons_qp_arg_set_ric_alg
#define CONs_QP_ARG_SET_CONs_LAST_STAGE s_cons_qp_arg_set_cons_last_stage
#define CONs_QP_ARG_SET_COMP_PRIM_SOL s_cons_qp_arg_set_comp_prim_sol
#define CONs_QP_ARG_SET_COMP_DUAL_SOL_EQ s_cons_qp_arg_set_comp_dual_sol_eq
#define CONs_QP_ARG_SET_COMP_DUAL_SOL_INEQ s_cons_qp_arg_set_comp_dual_sol_ineq
#define CONs_QP_WS_MEMSIZE s_cons_qp_ws_memsize
#define CONs_QP_WS_CREATE s_cons_qp_ws_create
#define CONs_QP_COND s_cons_qp_cond
#define CONs_QP_CONs_LHS s_cons_qp_cons_lhs
#define CONs_QP_CONs_RHS s_cons_qp_cons_rhs
#define CONs_QP_EXPANs_SOL s_cons_qp_expans_sol
#define CONs_QP_EXPANs_PRIMAL_SOL s_cons_qp_expans_primal_sol
#define CONs_QP_UPDATE s_cons_qp_update



#include "x_cond.c"

