#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_s_aux.h"
#include "blasfeo/blasfeo_s_blas.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/auxiliary/aux_mem.h"
#include "hpipm/auxiliary/aux_string.h"
#include "hpipm/ocp/s_ocp_qp.h"
#include "hpipm/ocp/s_ocp_qp_dim.h"
#include "hpipm/ocp/s_ocp_qp_red.h"
#include "hpipm/ocp/s_ocp_qp_sol.h"


#define SINGLE_PRECISION


#define AXPY blasfeo_saxpy
#define CREATE_STRVEC blasfeo_create_svec
#define GECP blasfeo_sgecp
#define GEMV_N blasfeo_sgemv_n
#define GEMV_T blasfeo_sgemv_t
#define OCP_QP s_ocp_qp
#define OCP_QP_SOL s_ocp_qp_sol
#define OCP_QP_DIM s_ocp_qp_dim
#define OCP_QP_REDUCE_EQ_DOF_ARG s_ocp_qp_reduce_eq_dof_arg
#define OCP_QP_REDUCE_EQ_DOF_WS s_ocp_qp_reduce_eq_dof_ws
#define MATEL BLASFEO_SMATEL
#define REAL float
#define SIZE_STRVEC blasfeo_memsize_svec
#define STRVEC blasfeo_svec
#define SYMV_L blasfeo_ssymv_l
#define VECAD_SP blasfeo_svecad_sp
#define VECCP blasfeo_sveccp
#define VECEL BLASFEO_SVECEL
#define VECSE blasfeo_svecse

#define OCP_QP_DIM_REDUCE_EQ_DOF s_ocp_qp_dim_reduce_eq_dof
#define OCP_QP_REDUCE_EQ_DOF_ARG_MEMSIZE s_ocp_qp_reduce_eq_dof_arg_memsize
#define OCP_QP_REDUCE_EQ_DOF_ARG_CREATE s_ocp_qp_reduce_eq_dof_arg_create
#define OCP_QP_REDUCE_EQ_DOF_ARG_SET_DEFAULT s_ocp_qp_reduce_eq_dof_arg_set_default
#define OCP_QP_REDUCE_EQ_DOF_ARG_SET_ALIAS_UNCHANGED s_ocp_qp_reduce_eq_dof_arg_set_alias_unchanged
#define OCP_QP_REDUCE_EQ_DOF_ARG_SET_COMP_PRIM_SOL s_ocp_qp_reduce_eq_dof_arg_set_comp_prim_sol
#define OCP_QP_REDUCE_EQ_DOF_ARG_SET_COMP_DUAL_SOL_EQ s_ocp_qp_reduce_eq_dof_arg_set_comp_dual_sol_eq
#define OCP_QP_REDUCE_EQ_DOF_ARG_SET_COMP_DUAL_SOL_INEQ s_ocp_qp_reduce_eq_dof_arg_set_comp_dual_sol_ineq
#define OCP_QP_REDUCE_EQ_DOF_WS_MEMSIZE s_ocp_qp_reduce_eq_dof_ws_memsize
#define OCP_QP_REDUCE_EQ_DOF_WS_CREATE s_ocp_qp_reduce_eq_dof_ws_create
#define OCP_QP_REDUCE_EQ_DOF s_ocp_qp_reduce_eq_dof
#define OCP_QP_REDUCE_EQ_DOF_LHS s_ocp_qp_reduce_eq_dof_lhs
#define OCP_QP_REDUCE_EQ_DOF_RHS s_ocp_qp_reduce_eq_dof_rhs
#define OCP_QP_RESTORE_EQ_DOF s_ocp_qp_restore_eq_dof


#include "x_ocp_qp_red.c"
