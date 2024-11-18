#ifndef HPIPM_D_COND_qcqp_H_
#define HPIPM_D_COND_qcqp_H_

#include "tinyhpipm/blas.h"
#include "tinyhpipm/common.h"

#include "tinyhpipm/dense/d_dense_qcqp.h"
#include "tinyhpipm/dense/d_dense_qcqp_sol.h"
#include "tinyhpipm/ocp/d_ocp_qcqp.h"
#include "tinyhpipm/ocp/d_ocp_qcqp_dim.h"
#include "tinyhpipm/ocp/d_ocp_qcqp_sol.h"

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
    struct mat* hess_array;  // TODO remove
    struct mat* zero_hess;  // TODO remove
    struct vec* grad_array;  // TODO remove
    struct vec* zero_grad;  // TODO remove
    struct vec* tmp_nvc;
    struct vec* tmp_nuxM;
    struct mat* GammaQ;
    struct mat* tmp_DCt;
    struct mat* tmp_nuM_nxM;
    //	struct vec *d_qp;
    //	struct vec *d_mask_qp;
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
void d_cond_qcqp_qc(struct d_ocp_qcqp* ocp_qp, struct mat* Hq2, int* Hq_nzero2, struct mat* Ct2, struct vec* d2, struct d_cond_qcqp_arg* cond_arg, struct d_cond_qcqp_ws* cond_ws);
//
void d_cond_qcqp_qc_lhs(struct d_ocp_qcqp* ocp_qp, struct mat* Hq2, int* Hq_nzero2, struct mat* Ct2, struct d_cond_qcqp_arg* cond_arg, struct d_cond_qcqp_ws* cond_ws);
//
void d_cond_qcqp_qc_rhs(struct d_ocp_qcqp* ocp_qp, struct vec* d2, struct d_cond_qcqp_arg* cond_arg, struct d_cond_qcqp_ws* cond_ws);
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


#endif  // HPIPM_D_COND_qcqp_H_
