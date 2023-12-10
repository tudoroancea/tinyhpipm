#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_s_aux.h"
#include "blasfeo/blasfeo_s_blas.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/auxiliary/aux_mem.h"
#include "hpipm/ocp/s_ocp_qcqp_dim.h"
#include "hpipm/ocp/s_ocp_qcqp_res.h"


#define SINGLE_PRECISION
#define BLASFEO_VECEL BLASFEO_SVECEL


#define AXPY blasfeo_saxpy
#define COLEX blasfeo_scolex
#define CREATE_STRVEC blasfeo_create_svec
#define CVT_STRVEC2VEC blasfeo_unpack_svec
#define DOT blasfeo_sdot
#define GEMV_DIAG blasfeo_sgemv_d
#define GEMV_NT blasfeo_sgemv_nt
#define OCP_QCQP s_ocp_qcqp
#define OCP_QCQP_DIM s_ocp_qcqp_dim
#define OCP_QCQP_RES s_ocp_qcqp_res
#define OCP_QCQP_RES_WS s_ocp_qcqp_res_ws
#define OCP_QCQP_SOL s_ocp_qcqp_sol
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


#define OCP_QCQP_RES_MEMSIZE s_ocp_qcqp_res_memsize
#define OCP_QCQP_RES_CREATE s_ocp_qcqp_res_create
#define OCP_QCQP_RES_WS_MEMSIZE s_ocp_qcqp_res_ws_memsize
#define OCP_QCQP_RES_WS_CREATE s_ocp_qcqp_res_ws_create
#define OCP_QCQP_RES_COMPUTE s_ocp_qcqp_res_compute
#define OCP_QCQP_RES_COMPUTE_LIN s_ocp_qcqp_res_compute_lin
#define OCP_QCQP_RES_COMPUTE_INF_NORM s_ocp_qcqp_res_compute_inf_norm
#define OCP_QCQP_RES_GET_ALL s_ocp_qcqp_res_get_all
#define OCP_QCQP_RES_GET_MAX_RES_STAT s_ocp_qcqp_res_get_max_res_stat
#define OCP_QCQP_RES_GET_MAX_RES_EQ s_ocp_qcqp_res_get_max_res_eq
#define OCP_QCQP_RES_GET_MAX_RES_INEQ s_ocp_qcqp_res_get_max_res_ineq
#define OCP_QCQP_RES_GET_MAX_RES_COMP s_ocp_qcqp_res_get_max_res_comp


#include "x_ocp_qcqp_res.c"
