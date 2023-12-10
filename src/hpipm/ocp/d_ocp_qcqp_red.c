#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_d_aux.h"
#include "blasfeo/blasfeo_d_blas.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/auxiliary/aux_mem.h"
#include "hpipm/auxiliary/aux_string.h"
#include "hpipm/ocp/d_ocp_qcqp.h"
#include "hpipm/ocp/d_ocp_qcqp_dim.h"
#include "hpipm/ocp/d_ocp_qcqp_red.h"
#include "hpipm/ocp/d_ocp_qcqp_sol.h"


#define DOUBLE_PRECISION


#define AXPY blasfeo_daxpy
#define COLEX blasfeo_dcolex
#define CREATE_STRVEC blasfeo_create_dvec
#define DOT blasfeo_ddot
#define GECP blasfeo_dgecp
#define GEMV_N blasfeo_dgemv_n
#define GEMV_T blasfeo_dgemv_t
#define OCP_QCQP d_ocp_qcqp
#define OCP_QCQP_DIM d_ocp_qcqp_dim
#define OCP_QCQP_DIM_SET_NU d_ocp_qcqp_dim_set_nu
#define OCP_QCQP_DIM_SET_NX d_ocp_qcqp_dim_set_nx
#define OCP_QCQP_DIM_SET_NBU d_ocp_qcqp_dim_set_nbu
#define OCP_QCQP_DIM_SET_NBX d_ocp_qcqp_dim_set_nbx
#define OCP_QCQP_DIM_SET_NG d_ocp_qcqp_dim_set_ng
#define OCP_QCQP_DIM_SET_NQ d_ocp_qcqp_dim_set_nq
#define OCP_QCQP_DIM_SET_NS d_ocp_qcqp_dim_set_ns
#define OCP_QCQP_DIM_SET_NSBU d_ocp_qcqp_dim_set_nsbu
#define OCP_QCQP_DIM_SET_NSBX d_ocp_qcqp_dim_set_nsbx
#define OCP_QCQP_DIM_SET_NSG d_ocp_qcqp_dim_set_nsg
#define OCP_QCQP_DIM_SET_NSQ d_ocp_qcqp_dim_set_nsq
#define OCP_QCQP_DIM_SET_NBUE d_ocp_qcqp_dim_set_nbue
#define OCP_QCQP_DIM_SET_NBXE d_ocp_qcqp_dim_set_nbxe
#define OCP_QCQP_DIM_SET_NGE d_ocp_qcqp_dim_set_nge
#define OCP_QCQP_DIM_SET_NQE d_ocp_qcqp_dim_set_nqe
#define OCP_QCQP_REDUCE_EQ_DOF_ARG d_ocp_qcqp_reduce_eq_dof_arg
#define OCP_QCQP_REDUCE_EQ_DOF_WS d_ocp_qcqp_reduce_eq_dof_ws
#define OCP_QCQP_SOL d_ocp_qcqp_sol
#define MATEL BLASFEO_DMATEL
#define REAL double
#define SIZE_STRVEC blasfeo_memsize_dvec
#define STRVEC blasfeo_dvec
#define SYMV_L blasfeo_dsymv_l
#define VECAD_SP blasfeo_dvecad_sp
#define VECCP blasfeo_dveccp
#define VECEL BLASFEO_DVECEL
#define VECSE blasfeo_dvecse

#define OCP_QCQP_DIM_REDUCE_EQ_DOF d_ocp_qcqp_dim_reduce_eq_dof
#define OCP_QCQP_REDUCE_EQ_DOF_ARG_MEMSIZE d_ocp_qcqp_reduce_eq_dof_arg_memsize
#define OCP_QCQP_REDUCE_EQ_DOF_ARG_CREATE d_ocp_qcqp_reduce_eq_dof_arg_create
#define OCP_QCQP_REDUCE_EQ_DOF_ARG_SET_DEFAULT d_ocp_qcqp_reduce_eq_dof_arg_set_default
#define OCP_QCQP_REDUCE_EQ_DOF_ARG_SET_ALIAS_UNCHANGED d_ocp_qcqp_reduce_eq_dof_arg_set_alias_unchanged
#define OCP_QCQP_REDUCE_EQ_DOF_ARG_SET_COMP_PRIM_SOL d_ocp_qcqp_reduce_eq_dof_arg_set_comp_prim_sol
#define OCP_QCQP_REDUCE_EQ_DOF_ARG_SET_COMP_DUAL_SOL_EQ d_ocp_qcqp_reduce_eq_dof_arg_set_comp_dual_sol_eq
#define OCP_QCQP_REDUCE_EQ_DOF_ARG_SET_COMP_DUAL_SOL_INEQ d_ocp_qcqp_reduce_eq_dof_arg_set_comp_dual_sol_ineq
#define OCP_QCQP_REDUCE_EQ_DOF_WS_MEMSIZE d_ocp_qcqp_reduce_eq_dof_ws_memsize
#define OCP_QCQP_REDUCE_EQ_DOF_WS_CREATE d_ocp_qcqp_reduce_eq_dof_ws_create
#define OCP_QCQP_REDUCE_EQ_DOF d_ocp_qcqp_reduce_eq_dof
#define OCP_QCQP_REDUCE_EQ_DOF_LHS d_ocp_qcqp_reduce_eq_dof_lhs
#define OCP_QCQP_REDUCE_EQ_DOF_RHS d_ocp_qcqp_reduce_eq_dof_rhs
#define OCP_QCQP_RESTORE_EQ_DOF d_ocp_qcqp_restore_eq_dof


#include "x_ocp_qcqp_red.c"