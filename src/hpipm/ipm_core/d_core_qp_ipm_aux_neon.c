#include <arm_neon.h>
#include <math.h>

#include "hpipm/ipm_core/d_core_qp_ipm.h"
#include "hpipm/ipm_core/d_core_qp_ipm_aux.h"

void d_compute_Gamma_gamma_qp(double* res_d, double* res_m, struct d_core_qp_ipm_workspace* cws) {
    int nc = cws->nc;

    double* lam = cws->lam;
    double* t = cws->t;
    double* t_inv = cws->t_inv;
    double* Gamma = cws->Gamma;
    double* gamma = cws->gamma;
    double lam_min = cws->lam_min;
    double t_min = cws->t_min;
    double t_min_inv = cws->t_min_inv;

    float64x2_t lam_vec, t_vec, res_d_vec, res_m_vec, t_inv_vec, lam_min_vec, t_min_vec, t_min_inv_vec, t_inv_tmp_vec, lam_tmp_vec, Gamma_vec,
            gamma_vec;
    t_min_vec = vdupq_n_f64(t_min);
    t_min_inv_vec = vdupq_n_f64(t_min_inv);
    lam_min_vec = vdupq_n_f64(lam_min);
    int ii = 0;

#if 0
    if (cws->t_lam_min == 1) {
        for (; ii < nc - 1; ii += 2) {
            // load lam and t in neon vectors
            lam_vec = vld1q_f64(&lam[ii]);
            t_vec = vld1q_f64(&t[ii]);

            // compute t_inv
            t_inv_vec = vrecpeq_f64(t_vec);
            // t_inv_vec = vmulq_f64(vrecpsq_f64(t_vec, t_inv_vec), t_inv_vec);  // Refine estimate
            vst1q_f64(&t_inv[ii], t_inv_vec);

            // compare lam and t with lam_min and t_min, select only the
            t_inv_tmp_vec = vbslq_f64(vcltq_f64(t_vec, t_min_vec), t_min_inv_vec, t_inv_vec);
            lam_tmp_vec = vbslq_f64(vcltq_f64(lam_vec, lam_min_vec), lam_min_vec, lam_vec);

            // compute Gamma and store it back in memory
            Gamma_vec = vmulq_f64(t_inv_tmp_vec, lam_tmp_vec);
            vst1q_f64(&Gamma[ii], Gamma_vec);

            // compute gamma and store it back in memory
            res_d_vec = vld1q_f64(&res_d[ii]);
            res_m_vec = vld1q_f64(&res_m[ii]);
            gamma_vec = vmulq_f64(t_inv_vec, vsubq_f64(res_m_vec, vmulq_f64(lam_vec, res_d_vec)));
            vst1q_f64(&gamma[ii], gamma_vec);
        }
    } else {
        for (; ii < nc - 1; ii += 2) {
            // load lam and t in neon vectors
            lam_vec = vld1q_f64(&lam[ii]);
            t_vec = vld1q_f64(&t[ii]);

            // compute t_inv
            t_inv_vec = vrecpeq_f64(t_vec);
            t_inv_vec = vmulq_f64(vrecpsq_f64(t_vec, t_inv_vec), t_inv_vec);  // Refine estimate
            vst1q_f64(&t_inv[ii], t_inv_vec);
            // TODO: do we need to store t_inv_vec back to memory (t_inv)?

            // compute Gamma and store it back in memory
            Gamma_vec = vmulq_f64(t_inv_vec, lam_vec);
            vst1q_f64(&Gamma[ii], Gamma_vec);

            // compute gamma and store it back in memory
            res_d_vec = vld1q_f64(&res_d[ii]);
            res_m_vec = vld1q_f64(&res_m[ii]);
            gamma_vec = vmulq_f64(t_inv_vec, vsubq_f64(res_m_vec, vmulq_f64(lam_vec, res_d_vec)));
            vst1q_f64(&gamma[ii], gamma_vec);
        }
    }
#endif
    for (; ii < nc; ii++) {
        t_inv[ii] = 1.0 / t[ii];
        Gamma[ii] = t_inv[ii] * lam[ii];
        gamma[ii] = t_inv[ii] * (res_m[ii] - lam[ii] * res_d[ii]);
    }
}


void d_compute_gamma_qp(double* res_d, double* res_m, struct d_core_qp_ipm_workspace* cws) {
    int nc = cws->nc;
    double* lam = cws->lam;
    double* t_inv = cws->t_inv;
    double* gamma = cws->gamma;

    float64x2_t lam_vec, t_inv_vec, gamma_vec;
    int ii = 0;
    // #if 0
    for (; ii < nc - 1; ii += 2) {
        // load lam and t_inv in neon vectors
        lam_vec = vld1q_f64(&lam[ii]);
        t_inv_vec = vld1q_f64(&t_inv[ii]);

        // compute gamma and store it back in memory
        gamma_vec = vmulq_f64(t_inv_vec, vsubq_f64(vld1q_f64(&res_m[ii]), vmulq_f64(lam_vec, vld1q_f64(&res_d[ii]))));
        vst1q_f64(&gamma[ii], gamma_vec);
    }
    // #endif
    for (; ii < nc; ii++) {
        gamma[ii] = t_inv[ii] * (res_m[ii] - lam[ii] * res_d[ii]);
    }
}


void d_compute_lam_t_qp(double* res_d, double* res_m, double* dlam, double* dt, struct d_core_qp_ipm_workspace* cws) {
    int nc = cws->nc;
    double* lam = cws->lam;
    double* t_inv = cws->t_inv;


    float64x2_t lam_vec, t_inv_vec, dlam_vec, dt_vec, res_d_vec, res_m_vec;
    float lam0;
    int ii = 0;

    // #if 0
    for (; ii < nc - 1; ii += 2) {
        lam_vec = vld1q_f64(&lam[ii]);
        dt_vec = vld1q_f64(&dt[ii]);
        res_d_vec = vld1q_f64(&res_d[ii]);
        res_m_vec = vld1q_f64(&res_m[ii]);
        t_inv_vec = vld1q_f64(&t_inv[ii]);
        dlam_vec = vaddq_f64(res_m_vec, vmulq_f64(lam_vec, dt_vec));
        dlam_vec = vsubq_f64(dlam_vec, vmulq_f64(lam_vec, res_d_vec));
        dlam_vec = vmulq_f64(dlam_vec, t_inv_vec);
        dlam_vec = vnegq_f64(dlam_vec);
        vst1q_f64(&dlam[ii], dlam_vec);

        dt_vec = vsubq_f64(dt_vec, res_d_vec);
        vst1q_f64(&dt[ii], dt_vec);
    }
    // #endif
    for (; ii < nc; ii++) {
        lam0 = lam[ii];
        dlam[ii] = -t_inv[ii] * (res_m[ii] + (lam0 * dt[ii]) - (lam0 * res_d[ii]));
        dt[ii] -= res_d[ii];
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

    float64x2_t dlam_vec, lam_vec, t_vec, dt_vec, tpr_vec;
    uint32x2_t mask;
    double tpr;
    int ii = 0;
#if 0
    for (; ii < nc - 1; ii += 2) {
        dlam_vec = vld1q_f64(&dlam[ii]);
        lam_vec = vld1q_f64(&lam[ii]);
        t_vec = vld1q_f64(&t[ii]);
        dt_vec = vld1q_f64(&dt[ii]);
        // alpha_dual = fmin(alpha_dual, vminvq_f64(vdivq_f64(lam_vec, dlam_vec)));
        // alpha_prim = fmin(alpha_prim, vminvq_f64(vdivq_f64(t_vec, dt_vec)));
        tpr_vec = vdupq_n_f64(alpha_dual);
        mask = vcltq_f64(lam_vec, vmulq_f64(dlam_vec, tpr_vec));

        tpr = vminvq_f64(vdivq_f64(lam_vec, dlam_vec));
        alpha_dual = alpha_dual < tpr ? alpha_dual : tpr;
        tpr = vminvq_f64(vdivq_f64(t_vec, dt_vec));
        alpha_prim = alpha_prim < tpr ? alpha_prim : tpr;
    }
#endif
    for (; ii < nc; ii++) {

        if (alpha_dual * dlam[ii] > lam[ii]) {
            alpha_dual = lam[ii] / dlam[ii];
        }
        if (alpha_prim * dt[ii] > t[ii]) {
            alpha_prim = t[ii] / dt[ii];
        }
    }

    alpha = alpha_prim > alpha_dual ? alpha_prim : alpha_dual;

    // store alpha
    cws->alpha_prim = -alpha_prim;
    cws->alpha_dual = -alpha_dual;
    cws->alpha = -alpha;

    return;
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

    double tmp_alpha_prim, tmp_alpha_dual;

    if (alpha < 1.0) {
        alpha_prim = alpha_prim * ((1.0 - alpha_prim) * 0.99 + alpha_prim * 0.9999999);
        alpha_dual = alpha_dual * ((1.0 - alpha_dual) * 0.99 + alpha_dual * 0.9999999);
        alpha = alpha * ((1.0 - alpha) * 0.99 + alpha * 0.9999999);
    }

    if (cws->split_step == 0) {
        tmp_alpha_prim = alpha;
        tmp_alpha_dual = alpha;
    } else {
        tmp_alpha_prim = alpha_prim;
        tmp_alpha_dual = alpha_dual;
    }


    float64x2_t tmp0, tmp1, tmp_alpha_prim_vec = vdupq_n_f64(tmp_alpha_prim), tmp_alpha_dual_vec = vdupq_n_f64(tmp_alpha_dual);
    int ii;

    // update v
    ii = 0;
    // #if 0
    for (; ii < nv - 1; ii += 2) {
        tmp0 = vld1q_f64(&v[ii]);
        vst1q_f64(&v_bkp[ii], tmp0);
        tmp0 = vaddq_f64(tmp0, vmulq_f64(vld1q_f64(&dv[ii]), tmp_alpha_prim_vec));
        vst1q_f64(&v[ii], tmp0);
    }
    // #endif
    for (; ii < nv; ii++) {
        v_bkp[ii] = v[ii];
        v[ii] += tmp_alpha_prim * dv[ii];
    }

    // update pi
    ii = 0;
    // #if 0
    for (; ii < ne - 1; ii += 2) {
        tmp0 = vld1q_f64(&pi[ii]);
        vst1q_f64(&pi_bkp[ii], tmp0);
        tmp0 = vaddq_f64(tmp0, vmulq_f64(vld1q_f64(&dpi[ii]), tmp_alpha_dual_vec));
        vst1q_f64(&pi[ii], tmp0);
    }
    // #endif
    for (; ii < ne; ii++) {
        pi_bkp[ii] = pi[ii];
        pi[ii] += tmp_alpha_dual * dpi[ii];
    }

    // update lam and t
    ii = 0;
    if (cws->t_lam_min == 2) {
        // #if 0
        float64x2_t lam_min_vec = vdupq_n_f64(cws->lam_min), t_min_vec = vdupq_n_f64(cws->t_min);
        for (; ii < nc - 1; ii += 2) {
            tmp0 = vld1q_f64(&lam[ii]);
            vst1q_f64(&lam_bkp[ii], tmp0);
            tmp1 = vld1q_f64(dlam + ii);  // TODO: is this syntax better than &dlam[ii]?
            tmp0 = vaddq_f64(tmp0, vmulq_f64(tmp1, tmp_alpha_dual_vec));
            tmp0 = vmaxq_f64(tmp0, lam_min_vec);  // TODO: does not preserve NaN ?
            vst1q_f64(&lam[ii], tmp0);

            tmp0 = vld1q_f64(&t[ii]);
            vst1q_f64(&t_bkp[ii], tmp0);
            tmp1 = vld1q_f64(dt + ii);
            tmp0 = vaddq_f64(tmp0, vmulq_f64(tmp1, tmp_alpha_prim_vec));
            tmp0 = vmaxq_f64(tmp0, t_min_vec);
            vst1q_f64(&t[ii], tmp0);
        }
        // #endif
        for (; ii < nc; ii++) {
            lam_bkp[ii] = lam[ii];
            lam[ii] += tmp_alpha_dual * dlam[ii];
            lam[ii] = lam[ii] <= cws->lam_min ? cws->lam_min : lam[ii];

            t_bkp[ii] = t[ii];
            t[ii] += tmp_alpha_prim * dt[ii];
            t[ii] = t[ii] <= cws->t_min ? cws->t_min : t[ii];
        }
    } else {  // split step
        // #if 0
        for (; ii < nc - 1; ii += 2) {
            tmp0 = vld1q_f64(&lam[ii]);
            vst1q_f64(&lam_bkp[ii], tmp0);
            tmp1 = vld1q_f64(dlam + ii);  // TODO: is this syntax better than &dlam[ii]?
            tmp0 = vaddq_f64(tmp0, vmulq_f64(tmp1, tmp_alpha_dual_vec));
            vst1q_f64(&lam[ii], tmp0);

            tmp0 = vld1q_f64(&t[ii]);
            vst1q_f64(&t_bkp[ii], tmp0);
            tmp1 = vld1q_f64(dt + ii);
            tmp0 = vaddq_f64(tmp0, vmulq_f64(tmp1, tmp_alpha_prim_vec));
            vst1q_f64(&t[ii], tmp0);
        }
        // #endif 0
        for (; ii < nc; ii++) {
            lam_bkp[ii] = lam[ii];
            lam[ii] += tmp_alpha_dual * dlam[ii];

            t_bkp[ii] = t[ii];
            t[ii] += tmp_alpha_prim * dt[ii];
        }
    }
}


void d_compute_mu_aff_qp(struct d_core_qp_ipm_workspace* cws) {
    // extract workspace members
    int nc = cws->nc;

    double* lam = cws->lam;
    double* t = cws->t;
    double* dlam = cws->dlam;
    double* dt = cws->dt;
    double alpha = cws->alpha;
    // this affects the minimum value of signa !!!
    //		alpha *= 0.99;

    double mu = 0;
    float64x2_t lam_vec, dlam_vec, t_vec, dt_vec, alpha_vec = vdupq_n_f64(alpha), tpr;

    int ii = 0;
    // #if 0
    for (; ii < nc - 1; ii += 2) {
        lam_vec = vld1q_f64(&lam[ii]);
        dlam_vec = vld1q_f64(&dlam[ii]);
        t_vec = vld1q_f64(&t[ii]);
        dt_vec = vld1q_f64(&dt[ii]);
        tpr = vmulq_f64(vaddq_f64(lam_vec, vmulq_f64(alpha_vec, dlam_vec)), vaddq_f64(t_vec, vmulq_f64(alpha_vec, dt_vec)));
        mu += vaddvq_f64(tpr);
    }
    // #endif
    for (; ii < nc; ii++) {
        mu += (lam[ii] + alpha * dlam[ii]) * (t[ii] + alpha * dt[ii]);
    }

    cws->mu_aff = mu * cws->nc_mask_inv;
}


void d_backup_res_m(struct d_core_qp_ipm_workspace* cws) {
    // extract workspace members
    int nc = cws->nc;

    double* res_m = cws->res_m;
    double* res_m_bkp = cws->res_m_bkp;

    float64x2_t res_m_vec;

    int ii = 0;
    // #if 0
    for (; ii < nc - 1; ii += 2) {
        res_m_vec = vld1q_f64(&res_m[ii]);
        vst1q_f64(&res_m_bkp[ii], res_m_vec);
    }
    // #endif
    for (; ii < nc; ii++) {
        res_m_bkp[ii] = res_m[ii];
    }
}


void d_compute_centering_correction_qp(struct d_core_qp_ipm_workspace* cws) {
    // extract workspace members
    int nc = cws->nc;

    double* dlam = cws->dlam;
    double* dt = cws->dt;
    double* res_m = cws->res_m;
    double* res_m_bkp = cws->res_m_bkp;

    float64x2_t dt_vec, dlam_vec, res_m_vec, res_m_bkp_vec;

    double sigma_mu = cws->sigma * cws->mu;
    sigma_mu = sigma_mu > cws->tau_min ? sigma_mu : cws->tau_min;
    float64x2_t sigma_mu_vec = vdupq_n_f64(sigma_mu);

    int ii = 0;
    // #if 0
    for (; ii < nc - 1; ii += 2) {
        res_m_vec = vld1q_f64(&res_m[ii]);
        res_m_bkp_vec = vld1q_f64(&res_m_bkp[ii]);
        dt_vec = vld1q_f64(&dt[ii]);
        dlam_vec = vld1q_f64(&dlam[ii]);
        res_m_vec = vsubq_f64(vaddq_f64(res_m_bkp_vec, vmulq_f64(dt_vec, dlam_vec)), sigma_mu_vec);
        vst1q_f64(&res_m[ii], res_m_vec);
    }
    // #endif
    for (; ii < nc; ii++) {
        res_m[ii] = res_m_bkp[ii] + dt[ii] * dlam[ii] - sigma_mu;
    }
}


void d_compute_centering_qp(struct d_core_qp_ipm_workspace* cws) {
    // extract workspace members
    int nc = cws->nc;

    double* res_m = cws->res_m;
    double* res_m_bkp = cws->res_m_bkp;

    double sigma_mu = cws->sigma * cws->mu;
    sigma_mu = sigma_mu > cws->tau_min ? sigma_mu : cws->tau_min;

    float64x2_t res_m_vec, res_m_bkp_vec, sigma_mu_vec = vdupq_n_f64(sigma_mu);

    int ii = 0;
    // #if 0
    for (; ii < nc - 1; ii += 2) {
        res_m_vec = vld1q_f64(&res_m[ii]);
        res_m_bkp_vec = vld1q_f64(&res_m_bkp[ii]);
        res_m_vec = vsubq_f64(res_m_bkp_vec, sigma_mu_vec);
        vst1q_f64(&res_m[ii], res_m_vec);
    }
    // #endif
    for (; ii < nc; ii++) {
        res_m[ii] = res_m_bkp[ii] - sigma_mu;
    }
}


void d_compute_tau_min_qp(struct d_core_qp_ipm_workspace* cws) {
    // extract workspace members
    int nc = cws->nc;

    double* res_m = cws->res_m;
    double* res_m_bkp = cws->res_m_bkp;
    double tau_min = cws->tau_min;

    float64x2_t res_m_vec, res_m_bkp_vec, tau_min_vec = vdupq_n_f64(tau_min);

    int ii = 0;
    // #if 0
    for (; ii < nc - 1; ii += 2) {
        res_m_vec = vld1q_f64(&res_m[ii]);
        res_m_bkp_vec = vld1q_f64(&res_m_bkp[ii]);
        res_m_vec = vsubq_f64(res_m_bkp_vec, tau_min_vec);
        vst1q_f64(&res_m[ii], res_m_vec);
    }
    // #endif
    for (; ii < nc; ii++) {
        res_m[ii] = res_m_bkp[ii] - tau_min;
    }
}
