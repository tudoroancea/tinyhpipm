#ifndef HPIPM_D_COND_AUX_H_
#define HPIPM_D_COND_AUX_H_

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_target.h"

#ifdef __cplusplus
extern "C" {
#endif


//
void d_cond_BAbt(struct d_ocp_qp* ocp_qp, struct blasfeo_dmat* BAbt2, struct blasfeo_dvec* b, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws);
//
void d_cond_BAt(struct d_ocp_qp* ocp_qp, struct blasfeo_dmat* BAbt2, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws);
//
void d_cond_b(struct d_ocp_qp* ocp_qp, struct blasfeo_dvec* b, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws);
//
void d_cond_RSQrq(struct d_ocp_qp* ocp_qp, struct blasfeo_dmat* RSQrq2, struct blasfeo_dvec* rq, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws);
//
void d_cond_RSQ(struct d_ocp_qp* ocp_qp, struct blasfeo_dmat* RSQrq2, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws);
//
void d_cond_rq(struct d_ocp_qp* ocp_qp, struct blasfeo_dvec* rq, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws);
//
void d_cond_DCtd(struct d_ocp_qp* ocp_qp, int* idxb2, struct blasfeo_dmat* DCt2, struct blasfeo_dvec* d2, struct blasfeo_dvec* d_mask2, int* idxs_rev2, struct blasfeo_dvec* Z2, struct blasfeo_dvec* z, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws);
//
void d_cond_DCt(struct d_ocp_qp* ocp_qp, int* idxb2, struct blasfeo_dmat* DCt2, int* idxs_rev2, struct blasfeo_dvec* Z2, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws);
//
void d_cond_d(struct d_ocp_qp* ocp_qp, struct blasfeo_dvec* d2, struct blasfeo_dvec* d_mask2, struct blasfeo_dvec* z, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws);
//
void d_expand_sol(struct d_ocp_qp* ocp_qp, struct d_dense_qp_sol* dense_qp_sol, struct d_ocp_qp_sol* ocp_qp_so, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws);
//
void d_expand_primal_sol(struct d_ocp_qp* ocp_qp, struct d_dense_qp_sol* dense_qp_sol, struct d_ocp_qp_sol* ocp_qp_so, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws);

//
void d_update_cond_BAbt(int* idxc, struct d_ocp_qp* ocp_qp, struct blasfeo_dmat* BAbt2, struct blasfeo_dvec* b, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws);
//
void d_update_cond_RSQrq_N2nx3(int* idxc, struct d_ocp_qp* ocp_qp, struct blasfeo_dmat* RSQrq2, struct blasfeo_dvec* rq, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws);
//
void d_update_cond_DCtd(int* idxc, struct d_ocp_qp* ocp_qp, int* idxb2, struct blasfeo_dmat* DCt2, struct blasfeo_dvec* d2, int* idxs2, struct blasfeo_dvec* Z2, struct blasfeo_dvec* z, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif  // HPIPM_D_COND_AUX_H_
