#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_d_aux.h"
#include "blasfeo/blasfeo_d_blas.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/auxiliary/aux_mem.h"
#include "hpipm/dense/d_dense_qcqp.h"
#include "hpipm/dense/d_dense_qcqp_dim.h"
#include "hpipm/dense/d_dense_qcqp_res.h"
#include "hpipm/dense/d_dense_qcqp_sol.h"


#define DOUBLE_PRECISION


#define AXPY blasfeo_daxpy
#define CREATE_STRVEC blasfeo_create_dvec
#define COLEX blasfeo_dcolex
#define CVT_STRVEC2VEC blasfeo_unpack_dvec
#define DENSE_QCQP d_dense_qcqp
#define DENSE_QCQP_DIM d_dense_qcqp_dim
#define DENSE_QCQP_RES d_dense_qcqp_res
#define DENSE_QCQP_RES_WS d_dense_qcqp_res_ws
#define DENSE_QCQP_SOL d_dense_qcqp_sol
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
#define VECEX_SP blasfeo_dvecex_sp
#define VECMULACC blasfeo_dvecmulacc
#define VECMULDOT blasfeo_dvecmuldot
#define VECNRM_INF blasfeo_dvecnrm_inf


#define DENSE_QCQP_RES_MEMSIZE d_dense_qcqp_res_memsize
#define DENSE_QCQP_RES_CREATE d_dense_qcqp_res_create
#define DENSE_QCQP_RES_WS_MEMSIZE d_dense_qcqp_res_ws_memsize
#define DENSE_QCQP_RES_WS_CREATE d_dense_qcqp_res_ws_create
#define DENSE_QCQP_RES_COMPUTE d_dense_qcqp_res_compute
#define DENSE_QCQP_RES_COMPUTE_INF_NORM d_dense_qcqp_res_compute_inf_norm


#include "x_dense_qcqp_res.c"
