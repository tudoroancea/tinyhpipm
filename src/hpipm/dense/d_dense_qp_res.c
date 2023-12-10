#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_d_aux.h"
#include "blasfeo/blasfeo_d_blas.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/auxiliary/aux_mem.h"
#include "hpipm/dense/d_dense_qp.h"
#include "hpipm/dense/d_dense_qp_dim.h"
#include "hpipm/dense/d_dense_qp_res.h"
#include "hpipm/dense/d_dense_qp_sol.h"


#define DOUBLE_PRECISION


#define AXPY blasfeo_daxpy
#define CREATE_STRVEC blasfeo_create_dvec
#define UNPACK_VEC blasfeo_unpack_dvec
#define DENSE_QP d_dense_qp
#define DENSE_QP_DIM d_dense_qp_dim
#define DENSE_QP_RES d_dense_qp_res
#define DENSE_QP_RES_WS d_dense_qp_res_ws
#define DENSE_QP_SOL d_dense_qp_sol
#define DOT blasfeo_ddot
#define GEMV_DIAG blasfeo_dgemv_d
#define GEMV_NT blasfeo_dgemv_nt
#define REAL double
#define SIZE_STRVEC blasfeo_memsize_dvec
#define STRMAT blasfeo_dmat
#define STRVEC blasfeo_dvec
#define SYMV_L blasfeo_dsymv_l
#define VECAD_SP blasfeo_dvecad_sp
#define VECCP blasfeo_dveccp
#define VECCPSC blasfeo_dveccpsc
#define VECEX_SP blasfeo_dvecex_sp
#define VECMUL blasfeo_dvecmul
#define VECMULACC blasfeo_dvecmulacc
#define VECMULDOT blasfeo_dvecmuldot
#define VECNRM_INF blasfeo_dvecnrm_inf


#define DENSE_QP_RES_MEMSIZE d_dense_qp_res_memsize
#define DENSE_QP_RES_CREATE d_dense_qp_res_create
#define DENSE_QP_RES_WS_MEMSIZE d_dense_qp_res_ws_memsize
#define DENSE_QP_RES_WS_CREATE d_dense_qp_res_ws_create
#define DENSE_QP_RES_COMPUTE d_dense_qp_res_compute
#define DENSE_QP_RES_COMPUTE_LIN d_dense_qp_res_compute_lin
#define DENSE_QP_RES_COMPUTE_INF_NORM d_dense_qp_res_compute_inf_norm
#define DENSE_QP_RES_GET_ALL d_dense_qp_res_get_all


#include "x_dense_qp_res.c"
