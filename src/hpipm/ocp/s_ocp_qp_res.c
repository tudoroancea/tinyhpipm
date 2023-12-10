#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_s_aux.h"
#include "blasfeo/blasfeo_s_blas.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/auxiliary/aux_mem.h"
#include "hpipm/ocp/s_ocp_qp_dim.h"
#include "hpipm/ocp/s_ocp_qp_res.h"


#define BLASFEO_VECEL BLASFEO_SVECEL


#define AXPY blasfeo_saxpy
#define CREATE_STRVEC blasfeo_create_svec
#define DOT blasfeo_sdot
#define UNPACK_VEC blasfeo_unpack_svec
#define GEMV_DIAG blasfeo_sgemv_d
#define GEMV_NT blasfeo_sgemv_nt
#define OCP_QP s_ocp_qp
#define OCP_QP_DIM s_ocp_qp_dim
#define OCP_QP_RES s_ocp_qp_res
#define OCP_QP_RES_WS s_ocp_qp_res_ws
#define OCP_QP_SOL s_ocp_qp_sol
#define REAL float
#define SIZE_STRVEC blasfeo_memsize_svec
#define STRMAT blasfeo_smat
#define STRVEC blasfeo_svec
#define SYMV_L blasfeo_ssymv_l
#define VECAD_SP blasfeo_svecad_sp
#define VECCP blasfeo_sveccp
#define VECEX_SP blasfeo_svecex_sp
#define VECMULACC blasfeo_svecmulacc
#define VECMULDOT blasfeo_svecmuldot
#define VECNRM_INF blasfeo_svecnrm_inf


#define OCP_QP_RES_MEMSIZE s_ocp_qp_res_memsize
#define OCP_QP_RES_CREATE s_ocp_qp_res_create
#define OCP_QP_RES_WS_MEMSIZE s_ocp_qp_res_ws_memsize
#define OCP_QP_RES_WS_CREATE s_ocp_qp_res_ws_create
#define OCP_QP_RES_COMPUTE s_ocp_qp_res_compute
#define OCP_QP_RES_COMPUTE_LIN s_ocp_qp_res_compute_lin
#define OCP_QP_RES_COMPUTE_INF_NORM s_ocp_qp_res_compute_inf_norm
#define OCP_QP_RES_GET_ALL s_ocp_qp_res_get_all
#define OCP_QP_RES_GET_MAX_RES_STAT s_ocp_qp_res_get_max_res_stat
#define OCP_QP_RES_GET_MAX_RES_EQ s_ocp_qp_res_get_max_res_eq
#define OCP_QP_RES_GET_MAX_RES_INEQ s_ocp_qp_res_get_max_res_ineq
#define OCP_QP_RES_GET_MAX_RES_COMP s_ocp_qp_res_get_max_res_comp


#include "x_ocp_qp_res.c"
