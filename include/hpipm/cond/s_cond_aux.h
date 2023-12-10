#ifndef HPIPM_S_COND_AUX_H_
#define HPIPM_S_COND_AUX_H_

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_target.h"

#ifdef __cplusplus
extern "C" {
#endif


//
void s_cond_BAbt(struct s_ocp_qp* ocp_qp, struct blasfeo_smat* BAbt2, struct blasfeo_svec* b2, struct s_cond_qp_arg* cond_arg, struct s_cond_qp_ws* cond_ws);
//
void s_cond_BAt(struct s_ocp_qp* ocp_qp, struct blasfeo_smat* BAbt2, struct s_cond_qp_arg* cond_arg, struct s_cond_qp_ws* cond_ws);
//
void s_cond_b(struct s_ocp_qp* ocp_qp, struct blasfeo_svec* b2, struct s_cond_qp_arg* cond_arg, struct s_cond_qp_ws* cond_ws);
//
void s_cond_RSQrq(struct s_ocp_qp* ocp_qp, struct blasfeo_smat* RSQrq2, struct blasfeo_svec* rq2, struct s_cond_qp_arg* cond_arg, struct s_cond_qp_ws* cond_ws);
//
void s_cond_RSQ(struct s_ocp_qp* ocp_qp, struct blasfeo_smat* RSQrq2, struct s_cond_qp_arg* cond_arg, struct s_cond_qp_ws* cond_ws);
//
void s_cond_rq(struct s_ocp_qp* ocp_qp, struct blasfeo_svec* rq2, struct s_cond_qp_arg* cond_arg, struct s_cond_qp_ws* cond_ws);
//
void s_cond_DCtd(struct s_ocp_qp* ocp_qp, int* idxb2, struct blasfeo_smat* DCt2, struct blasfeo_svec* d2, struct blasfeo_svec* d_mask2, int* idxs_rev2, struct blasfeo_svec* Z2, struct blasfeo_svec* z2, struct s_cond_qp_arg* cond_arg, struct s_cond_qp_ws* cond_ws);
//
void s_cond_DCt(struct s_ocp_qp* ocp_qp, int* idxb2, struct blasfeo_smat* DCt2, int* idxs_rev2, struct blasfeo_svec* Z2, struct s_cond_qp_arg* cond_arg, struct s_cond_qp_ws* cond_ws);
//
void s_cond_d(struct s_ocp_qp* ocp_qp, struct blasfeo_svec* d2, struct blasfeo_svec* d_mask2, struct blasfeo_svec* z2, struct s_cond_qp_arg* cond_arg, struct s_cond_qp_ws* cond_ws);
//
void s_expand_sol(struct s_ocp_qp* ocp_qp, struct s_dense_qp_sol* dense_qp_sol, struct s_ocp_qp_sol* ocp_qp_sol, struct s_cond_qp_arg* cond_arg, struct s_cond_qp_ws* cond_ws);
//
void s_expand_primal_sol(struct s_ocp_qp* ocp_qp, struct s_dense_qp_sol* dense_qp_sol, struct s_ocp_qp_sol* ocp_qp_sol, struct s_cond_qp_arg* cond_arg, struct s_cond_qp_ws* cond_ws);

//
void s_update_cond_BAbt(int* idxc, struct s_ocp_qp* ocp_qp, struct blasfeo_smat* BAbt2, struct blasfeo_svec* b2, struct s_cond_qp_arg* cond_arg, struct s_cond_qp_ws* cond_ws);
//
void s_update_cond_RSQrq_N2nx3(int* idxc, struct s_ocp_qp* ocp_qp, struct blasfeo_smat* RSQrq2, struct blasfeo_svec* rq2, struct s_cond_qp_arg* cond_arg, struct s_cond_qp_ws* cond_ws);
//
void s_update_cond_DCtd(int* idxc, struct s_ocp_qp* ocp_qp, int* idxb2, struct blasfeo_smat* DCt2, struct blasfeo_svec* d2, int* idxs2, struct blasfeo_svec* Z2, struct blasfeo_svec* z2, struct s_cond_qp_arg* cond_arg, struct s_cond_qp_ws* cond_ws);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif  // HPIPM_S_COND_AUX_H_
