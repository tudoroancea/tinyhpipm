#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_s_aux.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/auxiliary/aux_mem.h"
#include "hpipm/auxiliary/aux_string.h"
#include "hpipm/ocp/s_ocp_qcqp_dim.h"
#include "hpipm/ocp/s_ocp_qcqp_sol.h"
#include "hpipm/ocp/s_ocp_qp.h"


#define CREATE_STRVEC blasfeo_create_svec
#define UNPACK_VEC blasfeo_unpack_svec
#define PACK_VEC blasfeo_pack_svec
#define OCP_QP s_ocp_qp
#define OCP_QCQP_DIM s_ocp_qcqp_dim
#define OCP_QCQP_SOL s_ocp_qcqp_sol
#define REAL float
#define STRVEC blasfeo_svec
#define SIZE_STRVEC blasfeo_memsize_svec
#define VECCP blasfeo_sveccp
#define VECSE blasfeo_svecse

#define OCP_QCQP_SOL_STRSIZE s_ocp_qcqp_sol_strsize
#define OCP_QCQP_SOL_MEMSIZE s_ocp_qcqp_sol_memsize
#define OCP_QCQP_SOL_CREATE s_ocp_qcqp_sol_create
#define OCP_QCQP_SOL_COPY_ALL s_ocp_qcqp_sol_copy_all
#define OCP_QCQP_SOL_GET_ALL s_ocp_qcqp_sol_get_all
#define OCP_QCQP_SOL_SET_ALL s_ocp_qcqp_sol_set_all
#define OCP_QCQP_SOL_GET s_ocp_qcqp_sol_get
#define OCP_QCQP_SOL_GET_U s_ocp_qcqp_sol_get_u
#define OCP_QCQP_SOL_GET_X s_ocp_qcqp_sol_get_x
#define OCP_QCQP_SOL_GET_SL s_ocp_qcqp_sol_get_sl
#define OCP_QCQP_SOL_GET_SU s_ocp_qcqp_sol_get_su
#define OCP_QCQP_SOL_GET_PI s_ocp_qcqp_sol_get_pi
#define OCP_QCQP_SOL_GET_LAM_LB s_ocp_qcqp_sol_get_lam_lb
#define OCP_QCQP_SOL_GET_LAM_UB s_ocp_qcqp_sol_get_lam_ub
#define OCP_QCQP_SOL_GET_LAM_LG s_ocp_qcqp_sol_get_lam_lg
#define OCP_QCQP_SOL_GET_LAM_UG s_ocp_qcqp_sol_get_lam_ug
#define OCP_QCQP_SOL_SET s_ocp_qcqp_sol_set
#define OCP_QCQP_SOL_SET_U s_ocp_qcqp_sol_set_u
#define OCP_QCQP_SOL_SET_X s_ocp_qcqp_sol_set_x
#define OCP_QCQP_SOL_SET_SL s_ocp_qcqp_sol_set_sl
#define OCP_QCQP_SOL_SET_SU s_ocp_qcqp_sol_set_su


#include "x_ocp_qcqp_sol.c"
