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
#define CONs_QP_ARG s_cons_qp_arg
#define CONs_QP_WS s_cons_qp_ws
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
#define VECAs_SP blasfeo_svecas_sp
#define VECCP blasfeo_sveccp
#define VECSE blasfeo_svecse

#define CONs_BABT s_cons_BAbt
#define CONs_BAT s_cons_BAt
#define CONs_B s_cons_b
#define CONs_RSQRQ s_cons_RSQrq
#define CONs_RSQ s_cons_RSQ
#define CONs_RQ s_cons_rq
#define CONs_DCTD s_cons_DCtd
#define CONs_DCT s_cons_DCt
#define CONs_D s_cons_d
#define EXPANs_SOL s_expans_sol
#define EXPANs_PRIMAL_SOL s_expans_primal_sol
#define UPDATE_CONs_BABT s_update_cons_BAbt
#define UPDATE_CONs_RSQRQ_N2NX3 s_update_cons_RSQrq_N2nx3
#define UPDATE_CONs_DCTD s_update_cons_DCtd



#include "x_cons_aux.c"

