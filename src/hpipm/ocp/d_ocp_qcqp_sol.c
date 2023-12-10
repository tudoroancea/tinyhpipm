#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_d_aux.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/auxiliary/aux_mem.h"
#include "hpipm/auxiliary/aux_string.h"
#include "hpipm/ocp/d_ocp_qcqp_dim.h"
#include "hpipm/ocp/d_ocp_qcqp_sol.h"
#include "hpipm/ocp/d_ocp_qp.h"


#define CREATE_STRVEC blasfeo_create_dvec
#define UNPACK_VEC blasfeo_unpack_dvec
#define PACK_VEC blasfeo_pack_dvec
#define OCP_QP d_ocp_qp
#define OCP_QCQP_DIM d_ocp_qcqp_dim
#define OCP_QCQP_SOL d_ocp_qcqp_sol
#define REAL double
#define STRVEC blasfeo_dvec
#define SIZE_STRVEC blasfeo_memsize_dvec
#define VECCP blasfeo_dveccp
#define VECSE blasfeo_dvecse

#define OCP_QCQP_SOL_STRSIZE d_ocp_qcqp_sol_strsize
#define OCP_QCQP_SOL_MEMSIZE d_ocp_qcqp_sol_memsize
#define OCP_QCQP_SOL_CREATE d_ocp_qcqp_sol_create
#define OCP_QCQP_SOL_COPY_ALL d_ocp_qcqp_sol_copy_all
#define OCP_QCQP_SOL_GET_ALL d_ocp_qcqp_sol_get_all
#define OCP_QCQP_SOL_SET_ALL d_ocp_qcqp_sol_set_all
#define OCP_QCQP_SOL_GET d_ocp_qcqp_sol_get
#define OCP_QCQP_SOL_GET_U d_ocp_qcqp_sol_get_u
#define OCP_QCQP_SOL_GET_X d_ocp_qcqp_sol_get_x
#define OCP_QCQP_SOL_GET_SL d_ocp_qcqp_sol_get_sl
#define OCP_QCQP_SOL_GET_SU d_ocp_qcqp_sol_get_su
#define OCP_QCQP_SOL_GET_PI d_ocp_qcqp_sol_get_pi
#define OCP_QCQP_SOL_GET_LAM_LB d_ocp_qcqp_sol_get_lam_lb
#define OCP_QCQP_SOL_GET_LAM_UB d_ocp_qcqp_sol_get_lam_ub
#define OCP_QCQP_SOL_GET_LAM_LG d_ocp_qcqp_sol_get_lam_lg
#define OCP_QCQP_SOL_GET_LAM_UG d_ocp_qcqp_sol_get_lam_ug
#define OCP_QCQP_SOL_SET d_ocp_qcqp_sol_set
#define OCP_QCQP_SOL_SET_U d_ocp_qcqp_sol_set_u
#define OCP_QCQP_SOL_SET_X d_ocp_qcqp_sol_set_x
#define OCP_QCQP_SOL_SET_SL d_ocp_qcqp_sol_set_sl
#define OCP_QCQP_SOL_SET_SU d_ocp_qcqp_sol_set_su


#include "x_ocp_qcqp_sol.c"
