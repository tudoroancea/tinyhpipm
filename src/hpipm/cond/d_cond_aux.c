#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_d_aux.h"
#include "blasfeo/blasfeo_d_blas.h"
#include "blasfeo/blasfeo_target.h"

#include "hpipm/cond/d_cond.h"
#include "hpipm/dense/d_dense_qp.h"
#include "hpipm/dense/d_dense_qp_sol.h"
#include "hpipm/ocp/d_ocp_qp.h"
#include "hpipm/ocp/d_ocp_qp_sol.h"


#define DOUBLE_PRECISION


#define AXPY blasfeo_daxpy
#define COND_QP_ARG d_cond_qp_arg
#define COND_QP_WS d_cond_qp_ws
#define DENSE_QP_SOL d_dense_qp_sol
#define DIAEX blasfeo_ddiaex
#define GEAD blasfeo_dgead
#define GECP blasfeo_dgecp
#define GEEX1 blasfeo_dgeex1
#define GESE blasfeo_dgese
#define GEMM_ND blasfeo_dgemm_nd
#define GEMM_NN blasfeo_dgemm_nn
#define GEMM_NT blasfeo_dgemm_nt
#define GEMV_D blasfeo_dgemv_d
#define GEMV_N blasfeo_dgemv_n
#define GEMV_T blasfeo_dgemv_t
#define OCP_QP d_ocp_qp
#define OCP_QP_SOL d_ocp_qp_sol
#define POTRF_L blasfeo_dpotrf_l
#define POTRF_L_MN blasfeo_dpotrf_l_mn
#define PRINT_MAT blasfeo_print_dmat
#define REAL double
#define ROWAD blasfeo_drowad
#define ROWEX blasfeo_drowex
#define ROWIN blasfeo_drowin
#define STRMAT blasfeo_dmat
#define STRVEC blasfeo_dvec
#define SYMV_L blasfeo_dsymv_l
#define SYRK_LN_MN blasfeo_dsyrk_ln_mn
#define TRCP_L blasfeo_dtrcp_l
#define TRTR_L blasfeo_dtrtr_l
#define TRMM_RLNN blasfeo_dtrmm_rlnn
#define VECAD_SP blasfeo_dvecad_sp
#define VECCP blasfeo_dveccp
#define VECSE blasfeo_dvecse

#define COND_BABT d_cond_BAbt
#define COND_BAT d_cond_BAt
#define COND_B d_cond_b
#define COND_RSQRQ d_cond_RSQrq
#define COND_RSQ d_cond_RSQ
#define COND_RQ d_cond_rq
#define COND_DCTD d_cond_DCtd
#define COND_DCT d_cond_DCt
#define COND_D d_cond_d
#define EXPAND_SOL d_expand_sol
#define EXPAND_PRIMAL_SOL d_expand_primal_sol
#define UPDATE_COND_BABT d_update_cond_BAbt
#define UPDATE_COND_RSQRQ_N2NX3 d_update_cond_RSQrq_N2nx3
#define UPDATE_COND_DCTD d_update_cond_DCtd


#include "x_cond_aux.c"
