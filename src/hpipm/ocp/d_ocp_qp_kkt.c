#include <math.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_d_aux.h"
#include "blasfeo/blasfeo_d_blas.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/ipm_core/d_core_qp_ipm.h"
#include "hpipm/ipm_core/d_core_qp_ipm_aux.h"
#include "hpipm/ocp/d_ocp_qp.h"
#include "hpipm/ocp/d_ocp_qp_dim.h"
#include "hpipm/ocp/d_ocp_qp_ipm.h"
#include "hpipm/ocp/d_ocp_qp_sol.h"


#define DOUBLE_PRECISION
#define BLASFEO_MATEL BLASFEO_DMATEL
#define BLASFEO_VECEL BLASFEO_DVECEL


#define AXPY blasfeo_daxpy
#define COLSC blasfeo_dcolsc
#define COMPUTE_LAM_T_QP d_compute_lam_t_qp
#define COMPUTE_GAMMA_GAMMA_QP d_compute_Gamma_gamma_qp
#define COMPUTE_GAMMA_QP d_compute_gamma_qp
#define CORE_QP_IPM_WORKSPACE d_core_qp_ipm_workspace
#define DIAAD_SP blasfeo_ddiaad_sp
#define DIARE blasfeo_ddiare
#define GEAD blasfeo_dgead
#define GECP blasfeo_dgecp
#define GELQF blasfeo_dgelqf
#define GELQF_PD blasfeo_dgelqf_pd
#define GELQF_PD_LA blasfeo_dgelqf_pd_la
#define GELQF_PD_LLA blasfeo_dgelqf_pd_lla
#define GEMM_NT blasfeo_dgemm_nt
#define GEMM_R_DIAG blasfeo_dgemm_nd
#define GEMV_N blasfeo_dgemv_n
#define GEMV_T blasfeo_dgemv_t
#define GESE blasfeo_dgese
#define OCP_QP d_ocp_qp
#define OCP_QP_IPM_ARG d_ocp_qp_ipm_arg
#define OCP_QP_IPM_WS d_ocp_qp_ipm_ws
#define OCP_QP_DIM d_ocp_qp_dim
#define OCP_QP_SOL d_ocp_qp_sol
#define POTRF_L blasfeo_dpotrf_l
#define POTRF_L_MN blasfeo_dpotrf_l_mn
#define PRINT_E_MAT d_print_exp_mat
#define PRINT_E_STRVEC blasfeo_print_exp_dvec
#define PRINT_E_TRAN_STRVEC blasfeo_print_exp_tran_dvec
#define PRINT_STRMAT d_print_strmat
#define PRINT_STRVEC blasfeo_print_dvec
#define PRINT_TRAN_STRVEC blasfeo_print_tran_dvec
#define REAL double
#define ROWAD_SP blasfeo_drowad_sp
#define ROWEX blasfeo_drowex
#define ROWIN blasfeo_drowin
#define STRMAT blasfeo_dmat
#define STRVEC blasfeo_dvec
#define SYRK_LN blasfeo_dsyrk_ln
#define SYRK_LN_MN blasfeo_dsyrk_ln_mn
#define SYRK_POTRF_LN_MN blasfeo_dsyrk_dpotrf_ln_mn
#define TRCP_L blasfeo_dtrcp_l
#define TRMM_RLNN blasfeo_dtrmm_rlnn
#define TRMV_LNN blasfeo_dtrmv_lnn
#define TRMV_LTN blasfeo_dtrmv_ltn
#define TRSV_LNN blasfeo_dtrsv_lnn
#define TRSV_LNN_MN blasfeo_dtrsv_lnn_mn
#define TRSV_LTN blasfeo_dtrsv_ltn
#define TRSV_LTN_MN blasfeo_dtrsv_ltn_mn
#define TRTR_L blasfeo_dtrtr_l
#define VECAD_SP blasfeo_dvecad_sp
#define VECCP blasfeo_dveccp
#define VECCPSC blasfeo_dveccpsc
#define VECEX_SP blasfeo_dvecex_sp
#define VECSC blasfeo_dvecsc


#define OCP_QP_FACT_SOLVE_KKT_UNCONSTR d_ocp_qp_fact_solve_kkt_unconstr
#define OCP_QP_FACT_SOLVE_KKT_STEP d_ocp_qp_fact_solve_kkt_step
#define OCP_QP_FACT_LQ_SOLVE_KKT_STEP d_ocp_qp_fact_lq_solve_kkt_step
#define OCP_QP_SOLVE_KKT_STEP d_ocp_qp_solve_kkt_step


#include "x_ocp_qp_kkt.c"
