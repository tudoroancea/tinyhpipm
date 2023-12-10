#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_s_aux.h"
#include "blasfeo/blasfeo_s_blas.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/auxiliary/aux_mem.h"
#include "hpipm/dense/s_dense_qp.h"
#include "hpipm/dense/s_dense_qp_dim.h"
#include "hpipm/dense/s_dense_qp_res.h"
#include "hpipm/dense/s_dense_qp_sol.h"


#define SINGLE_PRECISION



#define AXPY blasfeo_saxpy
#define CREATE_STRVEC blasfeo_create_svec
#define UNPACK_VEC blasfeo_unpack_svec
#define DENSE_QP s_dense_qp
#define DENSE_QP_DIM s_dense_qp_dim
#define DENSE_QP_RES s_dense_qp_res
#define DENSE_QP_RES_WS s_dense_qp_res_ws
#define DENSE_QP_SOL s_dense_qp_sol
#define DOT blasfeo_sdot
#define GEMV_DIAG blasfeo_sgemv_d
#define GEMV_NT blasfeo_sgemv_nt
#define REAL float
#define SIZE_STRVEC blasfeo_memsize_svec
#define STRMAT blasfeo_smat
#define STRVEC blasfeo_svec
#define SYMV_L blasfeo_ssymv_l
#define VECAD_SP blasfeo_svecad_sp
#define VECCP blasfeo_sveccp
#define VECCPSC blasfeo_sveccpsc
#define VECEX_SP blasfeo_svecex_sp
#define VECMUL blasfeo_svecmul
#define VECMULACC blasfeo_svecmulacc
#define VECMULDOT blasfeo_svecmuldot
#define VECNRM_INF blasfeo_svecnrm_inf



#define DENSE_QP_RES_MEMSIZE s_dense_qp_res_memsize
#define DENSE_QP_RES_CREATE s_dense_qp_res_create
#define DENSE_QP_RES_WS_MEMSIZE s_dense_qp_res_ws_memsize
#define DENSE_QP_RES_WS_CREATE s_dense_qp_res_ws_create
#define DENSE_QP_RES_COMPUTE s_dense_qp_res_compute
#define DENSE_QP_RES_COMPUTE_LIN s_dense_qp_res_compute_lin
#define DENSE_QP_RES_COMPUTE_INF_NORM s_dense_qp_res_compute_inf_norm
#define DENSE_QP_RES_GET_ALL s_dense_qp_res_get_all



#include "x_dense_qp_res.c"


