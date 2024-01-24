#include "hpipm/ipm_core/d_core_qp_ipm.h"


#define d_core_qp_ipm_workspace d_core_qp_ipm_workspace


#define d_compute_Gamma_gamma_qp d_compute_Gamma_gamma_qp
#define d_compute_gamma_qp d_compute_gamma_qp
#define d_compute_lam_t_qp d_compute_lam_t_qp
#define d_compute_alpha_qp d_compute_alpha_qp
#define d_update_var_qp d_update_var_qp
#define d_compute_mu_aff_qp d_compute_mu_aff_qp
#define d_backup_res_m d_backup_res_m
#define d_compute_centering_correction_qp d_compute_centering_correction_qp
#define d_compute_centering_qp d_compute_centering_qp
#define d_compute_tau_min_qp d_compute_tau_min_qp


void d_compute_Gamma_gamma_qp(double* res_d, double* res_m, struct d_core_qp_ipm_workspace* cws) {

    int nc = cws->nc;

    double* lam = cws->lam;
    double* t = cws->t;
    //	double *res_d = cws->res_d; // TODO rename d ???
    //	double *res_m = cws->res_m; // TODO rename m ???
    double* t_inv = cws->t_inv;
    double* Gamma = cws->Gamma;
    double* gamma = cws->gamma;
    double lam_min = cws->lam_min;
    double t_min = cws->t_min;
    double t_min_inv = cws->t_min_inv;
    double lam0, t0, t_inv_tmp, lam_tmp;

    // local variables
    int ii;

    // printf("\ngamma gamma\n");
    if (cws->t_lam_min == 1) {
        for (ii = 0; ii < nc; ii++) {
            lam0 = lam[ii];
            t0 = t[ii];
            t_inv[ii] = 1.0 / t0;
            t_inv_tmp = t0 < t_min ? t_min_inv : t_inv[ii];
            lam_tmp = lam0 < lam_min ? lam_min : lam0;
            Gamma[ii] = t_inv_tmp * lam_tmp;
            gamma[ii] = t_inv[ii] * (res_m[ii] - lam0 * res_d[ii]);
        }
    } else {
        for (ii = 0; ii < nc; ii++) {
            lam0 = lam[ii];
            t0 = t[ii];
            t_inv[ii] = 1.0 / t0;
            Gamma[ii] = t_inv[ii] * lam0;
            gamma[ii] = t_inv[ii] * (res_m[ii] - lam0 * res_d[ii]);
        }
    }

    
}


void d_compute_gamma_qp(double* res_d, double* res_m, struct d_core_qp_ipm_workspace* cws) {

    int nc = cws->nc;

    double* lam = cws->lam;
    //	double *res_m = cws->res_m;
    //	double *res_d = cws->res_d;
    double* t_inv = cws->t_inv;
    double* gamma = cws->gamma;
    double lam_min = cws->lam_min;
    double lam0;

    // local variables
    int ii;

    for (ii = 0; ii < nc; ii++) {
        lam0 = lam[ii];
        //		lam0 = lam0<lam_min ? lam_min : lam0;
        gamma[ii] = t_inv[ii] * (res_m[ii] - lam0 * res_d[ii]);
    }

    
}

void d_compute_lam_t_qp(double* res_d, double* res_m, double* dlam, double* dt, struct d_core_qp_ipm_workspace* cws) {

    int nc = cws->nc;

    double* lam = cws->lam;
    double* t_inv = cws->t_inv;
    double lam_min = cws->lam_min;
    double lam0;

    // local variables
    int ii;

    for (ii = 0; ii < nc; ii++) {
        lam0 = lam[ii];
        //		lam0 = lam0<lam_min ? lam_min : lam0;
        dlam[ii] = -t_inv[ii] * (res_m[ii] + (lam0 * dt[ii]) - (lam0 * res_d[ii]));
        dt[ii] -= res_d[ii];
        // TODO compute lamda alone ???
        //		dlam[ii] = - t_inv[ii] * (lam0*dt[ii] + res_m[ii]);
    }

    
}


void d_compute_alpha_qp(struct d_core_qp_ipm_workspace* cws) {

    // extract workspace members
    int nc = cws->nc;

    double* lam = cws->lam;
    double* t = cws->t;
    double* dlam = cws->dlam;
    double* dt = cws->dt;

    double alpha_prim = -1.0;
    double alpha_dual = -1.0;
    double alpha = -1.0;

#if 1


    // local variables
    int ii;

    for (ii = 0; ii < nc; ii++) {

        if (alpha_dual * dlam[ii + 0] > lam[ii + 0]) {
            alpha_dual = lam[ii + 0] / dlam[ii + 0];
        }
        if (alpha_prim * dt[ii + 0] > t[ii + 0]) {
            alpha_prim = t[ii + 0] / dt[ii + 0];
        }
    }

#else  // fraction to the boundary

    double mu = cws->mu;
    double tau = 1.0;
    //	double tau = 0.995;

    tau = tau > (1 - mu) ? tau : 1 - mu;

    // local variables
    int ii;

    for (ii = 0; ii < nc; ii++) {

        if (alpha_dual * dlam[ii + 0] > tau * lam[ii + 0]) {
            alpha_dual = tau * lam[ii + 0] / dlam[ii + 0];
        }
        if (alpha_prim * dt[ii + 0] > tau * t[ii + 0]) {
            alpha_prim = tau * t[ii + 0] / dt[ii + 0];
        }
    }

#endif
    alpha = alpha_prim > alpha_dual ? alpha_prim : alpha_dual;

    // store alpha
    cws->alpha_prim = -alpha_prim;
    cws->alpha_dual = -alpha_dual;
    cws->alpha = -alpha;

    
}


void d_update_var_qp(struct d_core_qp_ipm_workspace* cws) {

    // extract workspace members
    int nv = cws->nv;
    int ne = cws->ne;
    int nc = cws->nc;

    double* v = cws->v;
    double* pi = cws->pi;
    double* lam = cws->lam;
    double* t = cws->t;
    double* v_bkp = cws->v_bkp;
    double* pi_bkp = cws->pi_bkp;
    double* lam_bkp = cws->lam_bkp;
    double* t_bkp = cws->t_bkp;
    double* dv = cws->dv;
    double* dpi = cws->dpi;
    double* dlam = cws->dlam;
    double* dt = cws->dt;
    double alpha = cws->alpha;
    double alpha_prim = cws->alpha_prim;
    double alpha_dual = cws->alpha_dual;
    double lam_min = cws->lam_min;
    double t_min = cws->t_min;

    double tmp_alpha_prim, tmp_alpha_dual;

#if 0
	if(alpha<1.0)
		alpha *= 0.995;
#else
    if (alpha < 1.0) {
        alpha_prim = alpha_prim * ((1.0 - alpha_prim) * 0.99 + alpha_prim * 0.9999999);
        alpha_dual = alpha_dual * ((1.0 - alpha_dual) * 0.99 + alpha_dual * 0.9999999);
        alpha = alpha * ((1.0 - alpha) * 0.99 + alpha * 0.9999999);
    }
#endif

    // local variables
    int ii;

    if (cws->split_step == 0) {
        tmp_alpha_prim = alpha;
        tmp_alpha_dual = alpha;
    } else {
        tmp_alpha_prim = alpha_prim;
        tmp_alpha_dual = alpha_dual;
    }

    // update v
    for (ii = 0; ii < nv; ii++) {
        v_bkp[ii] = v[ii];
        v[ii] += tmp_alpha_prim * dv[ii];
    }

    // update pi
    for (ii = 0; ii < ne; ii++) {
        pi_bkp[ii] = pi[ii];
        pi[ii] += tmp_alpha_dual * dpi[ii];
    }

    if (cws->t_lam_min == 2) {
        // update lam
        for (ii = 0; ii < nc; ii++) {
            lam_bkp[ii] = lam[ii];
            lam[ii] += tmp_alpha_dual * dlam[ii];
            lam[ii] = lam[ii] <= lam_min ? lam_min : lam[ii];
        }

        // update t
        for (ii = 0; ii < nc; ii++) {
            t_bkp[ii] = t[ii];
            t[ii] += tmp_alpha_prim * dt[ii];
            t[ii] = t[ii] <= t_min ? t_min : t[ii];
        }
    } else {
        // update lam
        for (ii = 0; ii < nc; ii++) {
            lam_bkp[ii] = lam[ii];
            lam[ii] += tmp_alpha_dual * dlam[ii];
        }

        // update t
        for (ii = 0; ii < nc; ii++) {
            t_bkp[ii] = t[ii];
            t[ii] += tmp_alpha_prim * dt[ii];
        }
    }

    
}


void d_compute_mu_aff_qp(struct d_core_qp_ipm_workspace* cws) {

    int ii;

    // extract workspace members
    int nc = cws->nc;

    double* ptr_lam = cws->lam;
    double* ptr_t = cws->t;
    double* ptr_dlam = cws->dlam;
    double* ptr_dt = cws->dt;
    double alpha = cws->alpha;
    // this affects the minimum value of signa !!!
    //		alpha *= 0.99;

    double mu = 0;

    for (ii = 0; ii < nc; ii++) {
        mu += (ptr_lam[ii + 0] + alpha * ptr_dlam[ii + 0]) * (ptr_t[ii + 0] + alpha * ptr_dt[ii + 0]);
    }

    cws->mu_aff = mu * cws->nc_mask_inv;

    
}


void d_backup_res_m(struct d_core_qp_ipm_workspace* cws) {

    int ii;

    // extract workspace members
    int nc = cws->nc;

    double* ptr_res_m = cws->res_m;
    double* ptr_res_m_bkp = cws->res_m_bkp;

    for (ii = 0; ii < nc; ii++) {
        ptr_res_m_bkp[ii + 0] = ptr_res_m[ii + 0];
    }

    
}


void d_compute_centering_correction_qp(struct d_core_qp_ipm_workspace* cws) {

    int ii;

    // extract workspace members
    int nc = cws->nc;

    double* ptr_dlam = cws->dlam;
    double* ptr_dt = cws->dt;
    double* ptr_res_m = cws->res_m;
    double* ptr_res_m_bkp = cws->res_m_bkp;

    double sigma_mu = cws->sigma * cws->mu;
    sigma_mu = sigma_mu > cws->tau_min ? sigma_mu : cws->tau_min;

    for (ii = 0; ii < nc; ii++) {
        ptr_res_m[ii + 0] = ptr_res_m_bkp[ii + 0] + ptr_dt[ii + 0] * ptr_dlam[ii + 0] - sigma_mu;
    }

    
}


void d_compute_centering_qp(struct d_core_qp_ipm_workspace* cws) {

    int ii;

    // extract workspace members
    int nc = cws->nc;

    double* ptr_res_m = cws->res_m;
    double* ptr_res_m_bkp = cws->res_m_bkp;

    double sigma_mu = cws->sigma * cws->mu;
    sigma_mu = sigma_mu > cws->tau_min ? sigma_mu : cws->tau_min;

    for (ii = 0; ii < nc; ii++) {
        ptr_res_m[ii + 0] = ptr_res_m_bkp[ii + 0] - sigma_mu;
    }

    
}


void d_compute_tau_min_qp(struct d_core_qp_ipm_workspace* cws) {

    int ii;

    // extract workspace members
    int nc = cws->nc;

    double* ptr_res_m = cws->res_m;
    double* ptr_res_m_bkp = cws->res_m_bkp;

    double tau_min = cws->tau_min;

    for (ii = 0; ii < nc; ii++) {
        ptr_res_m[ii + 0] = ptr_res_m_bkp[ii + 0] - tau_min;
    }

    
}

