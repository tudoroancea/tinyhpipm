#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_s_aux.h"
#include "blasfeo/blasfeo_s_blas.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/auxiliary/aux_mem.h"
#include "hpipm/dense/s_dense_qcqp.h"
#include "hpipm/dense/s_dense_qcqp_dim.h"
#include "hpipm/dense/s_dense_qcqp_res.h"
#include "hpipm/dense/s_dense_qcqp_sol.h"


#define SINGLE_PRECISION



#define AXPY blasfeo_saxpy
#define CREATE_STRVEC blasfeo_create_svec
#define COLEX blasfeo_scolex
#define CVT_STRVEC2VEC blasfeo_unpack_svec
#define DENSE_QCQP s_dense_qcqp
#define DENSE_QCQP_DIM s_dense_qcqp_dim
#define DENSE_QCQP_RES s_dense_qcqp_res
#define DENSE_QCQP_RES_WS s_dense_qcqp_res_ws
#define DENSE_QCQP_SOL s_dense_qcqp_sol
#define DOT blasfeo_sdot
#define GEMV_DIAG blasfeo_sgemv_d
#define GEMV_NT blasfeo_sgemv_nt
#define REAL double
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



#define DENSE_QCQP_RES_MEMSIZE s_dense_qcqp_res_memsize
#define DENSE_QCQP_RES_CREATE s_dense_qcqp_res_create
#define DENSE_QCQP_RES_WS_MEMSIZE s_dense_qcqp_res_ws_memsize
#define DENSE_QCQP_RES_WS_CREATE s_dense_qcqp_res_ws_create
#define DENSE_QCQP_RES_COMPUTE s_dense_qcqp_res_compute
#define DENSE_QCQP_RES_COMPUTE_INF_NORM s_dense_qcqp_res_compute_inf_norm



#include "x_dense_qcqp_res.c"



