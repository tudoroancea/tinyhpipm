#ifndef HPIPM_S_CORE_QP_IPM_AUX_
#define HPIPM_S_CORE_QP_IPM_AUX_

#ifdef __cplusplus
extern "C" {
#endif

//
void s_compute_Gamma_gamma_qp(float* res_d, float* res_m, struct s_core_qp_ipm_workspace* rws);
//
void s_compute_gamma_qp(float* res_d, float* res_m, struct s_core_qp_ipm_workspace* rws);
//
void s_compute_lam_t_qp(float* res_d, float* res_m, float* dlam, float* dt, struct s_core_qp_ipm_workspace* rws);
//
void s_compute_alpha_qp(struct s_core_qp_ipm_workspace* rws);
//
void s_update_var_qp(struct s_core_qp_ipm_workspace* rws);
//
void s_compute_mu_aff_qp(struct s_core_qp_ipm_workspace* rws);
//
void s_backup_res_m(struct s_core_qp_ipm_workspace* rws);
//
void s_compute_centering_correction_qp(struct s_core_qp_ipm_workspace* rws);
//
void s_compute_centering_qp(struct s_core_qp_ipm_workspace* rws);
//
void s_compute_tau_min_qp(struct s_core_qp_ipm_workspace* rws);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // HPIPM_S_CORE_QP_IPM_AUX_
