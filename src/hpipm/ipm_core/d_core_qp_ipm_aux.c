#include "hpipm/ipm_core/d_core_qp_ipm.h"


#define CORE_QP_IPM_WORKSPACE d_core_qp_ipm_workspace
#define REAL double

#define COMPUTE_GAMMA_GAMMA_QP d_compute_Gamma_gamma_qp
#define COMPUTE_GAMMA_QP d_compute_gamma_qp
#define COMPUTE_LAM_T_QP d_compute_lam_t_qp
#define COMPUTE_ALPHA_QP d_compute_alpha_qp
#define UPDATE_VAR_QP d_update_var_qp
#define COMPUTE_MU_AFF_QP d_compute_mu_aff_qp
#define BACKUP_RES_M d_backup_res_m
#define COMPUTE_CENTERING_CORRECTION_QP d_compute_centering_correction_qp
#define COMPUTE_CENTERING_QP d_compute_centering_qp
#define COMPUTE_TAU_MIN_QP d_compute_tau_min_qp


#include "x_core_qp_ipm_aux.c"
