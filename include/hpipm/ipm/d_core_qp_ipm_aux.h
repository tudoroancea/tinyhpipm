#ifndef HPIPM_D_CORE_QP_IPM_AUX_
#define HPIPM_D_CORE_QP_IPM_AUX_

#ifdef __cplusplus
extern "C" {
#endif

//
void d_compute_Gamma_gamma_qp(double* res_d, double* res_m, struct d_core_qp_ipm_workspace* rws);
//
void d_compute_gamma_qp(double* res_d, double* res_m, struct d_core_qp_ipm_workspace* rws);
//
void d_compute_lam_t_qp(double* res_d, double* res_m, double* dlam, double* dt, struct d_core_qp_ipm_workspace* rws);
//
void d_compute_alpha_qp(struct d_core_qp_ipm_workspace* rws);
//
void d_update_var_qp(struct d_core_qp_ipm_workspace* rws);
//
void d_compute_mu_aff_qp(struct d_core_qp_ipm_workspace* rws);
//
void d_backup_res_m(struct d_core_qp_ipm_workspace* rws);
//
void d_compute_centering_correction_qp(struct d_core_qp_ipm_workspace* rws);
//
void d_compute_centering_qp(struct d_core_qp_ipm_workspace* rws);
//
void d_compute_tau_min_qp(struct d_core_qp_ipm_workspace* rws);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // HPIPM_D_CORE_QP_IPM_AUX_
