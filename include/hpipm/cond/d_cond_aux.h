#ifndef HPIPM_D_COND_AUX_H_
#define HPIPM_D_COND_AUX_H_

#include "hpipm/blas.h"
#include "hpipm/common.h"

#ifdef __cplusplus
extern "C" {
#endif


//
void d_cond_BAbt(struct d_ocp_qp* ocp_qp, struct mat* BAbt2, struct vec* b, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws);
//
void d_cond_BAt(struct d_ocp_qp* ocp_qp, struct mat* BAbt2, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws);
//
void d_cond_b(struct d_ocp_qp* ocp_qp, struct vec* b, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws);
//
void d_cond_RSQrq(struct d_ocp_qp* ocp_qp, struct mat* RSQrq2, struct vec* rq, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws);
//
void d_cond_RSQ(struct d_ocp_qp* ocp_qp, struct mat* RSQrq2, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws);
//
void d_cond_rq(struct d_ocp_qp* ocp_qp, struct vec* rq, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws);
//
void d_cond_DCtd(struct d_ocp_qp* ocp_qp, int* idxb2, struct mat* DCt2, struct vec* d2, struct vec* d_mask2, int* idxs_rev2, struct vec* Z2, struct vec* z, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws);
//
void d_cond_DCt(struct d_ocp_qp* ocp_qp, int* idxb2, struct mat* DCt2, int* idxs_rev2, struct vec* Z2, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws);
//
void d_cond_d(struct d_ocp_qp* ocp_qp, struct vec* d2, struct vec* d_mask2, struct vec* z, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws);
//
void d_expand_sol(struct d_ocp_qp* ocp_qp, struct d_dense_qp_sol* dense_qp_sol, struct d_ocp_qp_sol* ocp_qp_so, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws);
//
void d_expand_primal_sol(struct d_ocp_qp* ocp_qp, struct d_dense_qp_sol* dense_qp_sol, struct d_ocp_qp_sol* ocp_qp_so, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws);

//
void d_update_cond_BAbt(int* idxc, struct d_ocp_qp* ocp_qp, struct mat* BAbt2, struct vec* b, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws);
//
void d_update_cond_RSQrq_N2nx3(int* idxc, struct d_ocp_qp* ocp_qp, struct mat* RSQrq2, struct vec* rq, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws);
//
void d_update_cond_DCtd(int* idxc, struct d_ocp_qp* ocp_qp, int* idxb2, struct mat* DCt2, struct vec* d2, int* idxs2, struct vec* Z2, struct vec* z, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif  // HPIPM_D_COND_AUX_H_
