#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_s_aux.h"
#include "blasfeo/blasfeo_s_blas.h"
#include "blasfeo/blasfeo_target.h"

#include "hpipm/cond/s_cond.h"
#include "hpipm/dense/s_dense_qp.h"
#include "hpipm/dense/s_dense_qp_sol.h"
#include "hpipm/ocp/s_ocp_qp.h"
#include "hpipm/ocp/s_ocp_qp_sol.h"

#define SINGLE_PRECISION


#define AXPY blasfeo_saxpy
#define COND_QP_ARG s_cond_qp_arg
#define COND_QP_WS s_cond_qp_ws
#define DENSE_QP_SOL s_dense_qp_sol
#define DIAEX blasfeo_sdiaex
#define GEAD blasfeo_sgead
#define GECP blasfeo_sgecp
#define GEEX1 blasfeo_sgeex1
#define GESE blasfeo_sgese
#define GEMM_ND blasfeo_sgemm_nd
#define GEMM_NN blasfeo_sgemm_nn
#define GEMM_NT blasfeo_sgemm_nt
#define GEMV_D blasfeo_sgemv_d
#define GEMV_N blasfeo_sgemv_n
#define GEMV_T blasfeo_sgemv_t
#define OCP_QP s_ocp_qp
#define OCP_QP_SOL s_ocp_qp_sol
#define POTRF_L blasfeo_spotrf_l
#define POTRF_L_MN blasfeo_spotrf_l_mn
#define PRINT_MAT blasfeo_print_smat
#define REAL float
#define ROWAD blasfeo_srowad
#define ROWEX blasfeo_srowex
#define ROWIN blasfeo_srowin
#define STRMAT blasfeo_smat
#define STRVEC blasfeo_svec
#define SYMV_L blasfeo_ssymv_l
#define SYRK_LN_MN blasfeo_ssyrk_ln_mn
#define TRCP_L blasfeo_strcp_l
#define TRTR_L blasfeo_strtr_l
#define TRMM_RLNN blasfeo_strmm_rlnn
#define VECAD_SP blasfeo_svecad_sp
#define VECCP blasfeo_sveccp
#define VECSE blasfeo_svecse

#define COND_BABT s_cond_BAbt
#define COND_BAT s_cond_BAt
#define COND_B s_cond_b
#define COND_RSQRQ s_cond_RSQrq
#define COND_RSQ s_cond_RSQ
#define COND_RQ s_cond_rq
#define COND_DCTD s_cond_DCtd
#define COND_DCT s_cond_DCt
#define COND_D s_cond_d
#define EXPAND_SOL s_expand_sol
#define EXPAND_PRIMAL_SOL s_expand_primal_sol
#define UPDATE_COND_BABT s_update_cond_BAbt
#define UPDATE_COND_RSQRQ_N2NX3 s_update_cond_RSQrq_N2nx3
#define UPDATE_COND_DCTD s_update_cond_DCtd


#include "x_cond_aux.c"
