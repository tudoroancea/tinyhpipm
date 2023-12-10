#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_d_aux.h"
#include "blasfeo/blasfeo_d_blas.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/auxiliary/aux_mem.h"
#include "hpipm/auxiliary/aux_string.h"
#include "hpipm/ocp/d_ocp_qp.h"
#include "hpipm/ocp/d_ocp_qp_dim.h"
#include "hpipm/ocp/d_ocp_qp_red.h"
#include "hpipm/ocp/d_ocp_qp_sol.h"


#define DOUBLE_PRECISION


#define AXPY blasfeo_daxpy
#define CREATE_STRVEC blasfeo_create_dvec
#define GECP blasfeo_dgecp
#define GEMV_N blasfeo_dgemv_n
#define GEMV_T blasfeo_dgemv_t
#define OCP_QP d_ocp_qp
#define OCP_QP_DIM d_ocp_qp_dim
#define OCP_QP_REDUCE_EQ_DOF_ARG d_ocp_qp_reduce_eq_dof_arg
#define OCP_QP_REDUCE_EQ_DOF_WS d_ocp_qp_reduce_eq_dof_ws
#define OCP_QP_SOL d_ocp_qp_sol
#define MATEL BLASFEO_DMATEL
#define REAL double
#define SIZE_STRVEC blasfeo_memsize_dvec
#define STRVEC blasfeo_dvec
#define SYMV_L blasfeo_dsymv_l
#define VECAD_SP blasfeo_dvecad_sp
#define VECCP blasfeo_dveccp
#define VECEL BLASFEO_DVECEL
#define VECSE blasfeo_dvecse

#define OCP_QP_DIM_REDUCE_EQ_DOF d_ocp_qp_dim_reduce_eq_dof
#define OCP_QP_REDUCE_EQ_DOF_ARG_MEMSIZE d_ocp_qp_reduce_eq_dof_arg_memsize
#define OCP_QP_REDUCE_EQ_DOF_ARG_CREATE d_ocp_qp_reduce_eq_dof_arg_create
#define OCP_QP_REDUCE_EQ_DOF_ARG_SET_DEFAULT d_ocp_qp_reduce_eq_dof_arg_set_default
#define OCP_QP_REDUCE_EQ_DOF_ARG_SET_ALIAS_UNCHANGED d_ocp_qp_reduce_eq_dof_arg_set_alias_unchanged
#define OCP_QP_REDUCE_EQ_DOF_ARG_SET_COMP_PRIM_SOL d_ocp_qp_reduce_eq_dof_arg_set_comp_prim_sol
#define OCP_QP_REDUCE_EQ_DOF_ARG_SET_COMP_DUAL_SOL_EQ d_ocp_qp_reduce_eq_dof_arg_set_comp_dual_sol_eq
#define OCP_QP_REDUCE_EQ_DOF_ARG_SET_COMP_DUAL_SOL_INEQ d_ocp_qp_reduce_eq_dof_arg_set_comp_dual_sol_ineq
#define OCP_QP_REDUCE_EQ_DOF_WS_MEMSIZE d_ocp_qp_reduce_eq_dof_ws_memsize
#define OCP_QP_REDUCE_EQ_DOF_WS_CREATE d_ocp_qp_reduce_eq_dof_ws_create
#define OCP_QP_REDUCE_EQ_DOF d_ocp_qp_reduce_eq_dof
#define OCP_QP_REDUCE_EQ_DOF_LHS d_ocp_qp_reduce_eq_dof_lhs
#define OCP_QP_REDUCE_EQ_DOF_RHS d_ocp_qp_reduce_eq_dof_rhs
#define OCP_QP_RESTORE_EQ_DOF d_ocp_qp_restore_eq_dof


#include "x_ocp_qp_red.c"
