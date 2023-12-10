#ifndef HPIPM_D_COND_QCQP_H_
#define HPIPM_D_COND_QCQP_H_

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_target.h"

#include "hpipm/dense/d_dense_qcqp.h"
#include "hpipm/dense/d_dense_qcqp_sol.h"
#include "hpipm/ocp/d_ocp_qcqp.h"
#include "hpipm/ocp/d_ocp_qcqp_dim.h"
#include "hpipm/ocp/d_ocp_qcqp_sol.h"

#ifdef __cplusplus
extern "C" {
#endif


struct d_cond_qcqp_arg {
    struct d_cond_qp_arg* qp_arg;
    int cond_last_stage;  // condense last stage
    //	int cond_variant; // TODO
    int comp_prim_sol;  // primal solution (v)
    int comp_dual_sol_eq;  // dual solution equality constr (pi)
    int comp_dual_sol_ineq;  // dual solution equality constr (lam t)
    int square_root_alg;  // square root algorithm (faster but requires RSQ>0)
    hpipm_size_t memsize;
};


struct d_cond_qcqp_ws {
    struct d_cond_qp_ws* qp_ws;
    struct blasfeo_dmat* hess_array;  // TODO remove
    struct blasfeo_dmat* zero_hess;  // TODO remove
    struct blasfeo_dvec* grad_array;  // TODO remove
    struct blasfeo_dvec* zero_grad;  // TODO remove
    struct blasfeo_dvec* tmp_nvc;
    struct blasfeo_dvec* tmp_nuxM;
    struct blasfeo_dmat* GammaQ;
    struct blasfeo_dmat* tmp_DCt;
    struct blasfeo_dmat* tmp_nuM_nxM;
    //	struct blasfeo_dvec *d_qp;
    //	struct blasfeo_dvec *d_mask_qp;
    hpipm_size_t memsize;
};


//
hpipm_size_t d_cond_qcqp_arg_memsize();
//
void d_cond_qcqp_arg_create(struct d_cond_qcqp_arg* cond_arg, void* mem);
//
void d_cond_qcqp_arg_set_default(struct d_cond_qcqp_arg* cond_arg);
// set riccati-like algorithm: 0 classical, 1 square-root
void d_cond_qcqp_arg_set_ric_alg(int ric_alg, struct d_cond_qcqp_arg* cond_arg);
// condense last stage: 0 last stage disregarded, 1 last stage condensed too
void d_cond_qcqp_arg_set_cond_last_stage(int cond_last_stage, struct d_cond_qcqp_arg* cond_arg);

//
void d_cond_qcqp_compute_dim(struct d_ocp_qcqp_dim* ocp_dim, struct d_dense_qcqp_dim* dense_dim);
//
hpipm_size_t d_cond_qcqp_ws_memsize(struct d_ocp_qcqp_dim* ocp_dim, struct d_cond_qcqp_arg* cond_arg);
//
void d_cond_qcqp_ws_create(struct d_ocp_qcqp_dim* ocp_dim, struct d_cond_qcqp_arg* cond_arg, struct d_cond_qcqp_ws* cond_ws, void* mem);
//
void d_cond_qcqp_qc(struct d_ocp_qcqp* ocp_qp, struct blasfeo_dmat* Hq2, int* Hq_nzero2, struct blasfeo_dmat* Ct2, struct blasfeo_dvec* d2, struct d_cond_qcqp_arg* cond_arg, struct d_cond_qcqp_ws* cond_ws);
//
void d_cond_qcqp_qc_lhs(struct d_ocp_qcqp* ocp_qp, struct blasfeo_dmat* Hq2, int* Hq_nzero2, struct blasfeo_dmat* Ct2, struct d_cond_qcqp_arg* cond_arg, struct d_cond_qcqp_ws* cond_ws);
//
void d_cond_qcqp_qc_rhs(struct d_ocp_qcqp* ocp_qp, struct blasfeo_dvec* d2, struct d_cond_qcqp_arg* cond_arg, struct d_cond_qcqp_ws* cond_ws);
//
void d_cond_qcqp_cond(struct d_ocp_qcqp* ocp_qp, struct d_dense_qcqp* dense_qp, struct d_cond_qcqp_arg* cond_arg, struct d_cond_qcqp_ws* cond_ws);
//
void d_cond_qcqp_cond_rhs(struct d_ocp_qcqp* ocp_qp, struct d_dense_qcqp* dense_qp, struct d_cond_qcqp_arg* cond_arg, struct d_cond_qcqp_ws* cond_ws);
//
void d_cond_qcqp_cond_lhs(struct d_ocp_qcqp* ocp_qp, struct d_dense_qcqp* dense_qp, struct d_cond_qcqp_arg* cond_arg, struct d_cond_qcqp_ws* cond_ws);
//
void d_cond_qcqp_expand_sol(struct d_ocp_qcqp* ocp_qp, struct d_dense_qcqp_sol* dense_qp_sol, struct d_ocp_qcqp_sol* ocp_qp_sol, struct d_cond_qcqp_arg* cond_arg, struct d_cond_qcqp_ws* cond_ws);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif  // HPIPM_D_COND_QCQP_H_
