#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_d_aux.h"
#include "blasfeo/blasfeo_d_blas.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/auxiliary/aux_mem.h"
#include "hpipm/ocp/d_ocp_qcqp_dim.h"
#include "hpipm/ocp/d_ocp_qcqp_res.h"


#define DOUBLE_PRECISION
#define BLASFEO_VECEL BLASFEO_DVECEL


#define AXPY blasfeo_daxpy
#define COLEX blasfeo_dcolex
#define CREATE_STRVEC blasfeo_create_dvec
#define CVT_STRVEC2VEC blasfeo_unpack_dvec
#define DOT blasfeo_ddot
#define GEMV_DIAG blasfeo_dgemv_d
#define GEMV_NT blasfeo_dgemv_nt
#define OCP_QCQP d_ocp_qcqp
#define OCP_QCQP_DIM d_ocp_qcqp_dim
#define OCP_QCQP_RES d_ocp_qcqp_res
#define OCP_QCQP_RES_WS d_ocp_qcqp_res_ws
#define OCP_QCQP_SOL d_ocp_qcqp_sol
#define REAL double
#define SIZE_STRVEC blasfeo_memsize_dvec
#define STRMAT blasfeo_dmat
#define STRVEC blasfeo_dvec
#define SYMV_L blasfeo_dsymv_l
#define VECAD_SP blasfeo_dvecad_sp
#define VECCP blasfeo_dveccp
#define VECEX_SP blasfeo_dvecex_sp
#define VECMULACC blasfeo_dvecmulacc
#define VECMULDOT blasfeo_dvecmuldot
#define VECNRM_INF blasfeo_dvecnrm_inf


#define OCP_QCQP_RES_MEMSIZE d_ocp_qcqp_res_memsize
#define OCP_QCQP_RES_CREATE d_ocp_qcqp_res_create
#define OCP_QCQP_RES_WS_MEMSIZE d_ocp_qcqp_res_ws_memsize
#define OCP_QCQP_RES_WS_CREATE d_ocp_qcqp_res_ws_create
#define OCP_QCQP_RES_COMPUTE d_ocp_qcqp_res_compute
#define OCP_QCQP_RES_COMPUTE_LIN d_ocp_qcqp_res_compute_lin
#define OCP_QCQP_RES_COMPUTE_INF_NORM d_ocp_qcqp_res_compute_inf_norm
#define OCP_QCQP_RES_GET_ALL d_ocp_qcqp_res_get_all
#define OCP_QCQP_RES_GET_MAX_RES_STAT d_ocp_qcqp_res_get_max_res_stat
#define OCP_QCQP_RES_GET_MAX_RES_EQ d_ocp_qcqp_res_get_max_res_eq
#define OCP_QCQP_RES_GET_MAX_RES_INEQ d_ocp_qcqp_res_get_max_res_ineq
#define OCP_QCQP_RES_GET_MAX_RES_COMP d_ocp_qcqp_res_get_max_res_comp


#include "x_ocp_qcqp_res.c"
