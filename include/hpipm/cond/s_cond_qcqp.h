#ifndef HPIPM_S_COND_QCQP_H_
#define HPIPM_S_COND_QCQP_H_

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_target.h"

#include "hpipm/dense/s_dense_qcqp.h"
#include "hpipm/dense/s_dense_qcqp_sol.h"
#include "hpipm/ocp/s_ocp_qcqp.h"
#include "hpipm/ocp/s_ocp_qcqp_dim.h"
#include "hpipm/ocp/s_ocp_qcqp_sol.h"

#ifdef __cplusplus
extern "C" {
#endif


struct s_cond_qcqp_arg {
    struct s_cond_qp_arg* qp_arg;
    int cond_last_stage;  // condense last stage
    //	int cond_variant; // TODO
    int comp_prim_sol;  // primal solution (v)
    int comp_dual_sol_eq;  // dual solution equality constr (pi)
    int comp_dual_sol_ineq;  // dual solution equality constr (lam t)
    int square_root_alg;  // square root algorithm (faster but requires RSQ>0)
    hpipm_size_t memsize;
};


struct s_cond_qcqp_ws {
    struct s_cond_qp_ws* qp_ws;
    struct blasfeo_smat* hess_array;  // TODO remove
    struct blasfeo_smat* zero_hess;  // TODO remove
    struct blasfeo_svec* zero_grad;  // TODO remove
    struct blasfeo_svec* grad_array;  // TODO remove
    struct blasfeo_svec* tmp_nvc;
    struct blasfeo_svec* tmp_nuxM;
    struct blasfeo_smat* GammaQ;
    struct blasfeo_smat* tmp_DCt;
    struct blasfeo_smat* tmp_nuM_nxM;
    //	struct blasfeo_svec *d_qp;
    //	struct blasfeo_svec *d_mask_qp;
    hpipm_size_t memsize;
};


//
hpipm_size_t s_cond_qcqp_arg_memsize();
//
void s_cond_qcqp_arg_create(struct s_cond_qcqp_arg* cond_arg, void* mem);
//
void s_cond_qcqp_arg_set_default(struct s_cond_qcqp_arg* cond_arg);
// set riccati-like algorithm: 0 classical, 1 square-root
void s_cond_qcqp_arg_set_ric_alg(int ric_alg, struct s_cond_qcqp_arg* cond_arg);
// condense last stage: 0 last stage disregarded, 1 last stage condensed too
void s_cond_qcqp_arg_set_cond_last_stage(int cond_last_stage, struct s_cond_qcqp_arg* cond_arg);

//
void s_cond_qcqp_compute_dim(struct s_ocp_qcqp_dim* ocp_dim, struct s_dense_qcqp_dim* dense_dim);
//
hpipm_size_t s_cond_qcqp_ws_memsize(struct s_ocp_qcqp_dim* ocp_dim, struct s_cond_qcqp_arg* cond_arg);
//
void s_cond_qcqp_ws_create(struct s_ocp_qcqp_dim* ocp_dim, struct s_cond_qcqp_arg* cond_arg, struct s_cond_qcqp_ws* cond_ws, void* mem);
//
void s_cond_qcqp_qc(struct s_ocp_qcqp* ocp_qp, struct blasfeo_smat* Hq2, int* Hq_nzero2, struct blasfeo_smat* Ct2, struct blasfeo_svec* d2, struct s_cond_qcqp_arg* cond_arg, struct s_cond_qcqp_ws* cond_ws);
//
void s_cond_qcqp_qc_lhs(struct s_ocp_qcqp* ocp_qp, struct blasfeo_smat* Hq2, int* Hq_nzero2, struct blasfeo_smat* Ct2, struct s_cond_qcqp_arg* cond_arg, struct s_cond_qcqp_ws* cond_ws);
//
void s_cond_qcqp_qc_rhs(struct s_ocp_qcqp* ocp_qp, struct blasfeo_svec* d2, struct s_cond_qcqp_arg* cond_arg, struct s_cond_qcqp_ws* cond_ws);
//
void s_cond_qcqp_cond(struct s_ocp_qcqp* ocp_qp, struct s_dense_qcqp* dense_qp, struct s_cond_qcqp_arg* cond_arg, struct s_cond_qcqp_ws* cond_ws);
//
void s_cond_qcqp_cond_rhs(struct s_ocp_qcqp* ocp_qp, struct s_dense_qcqp* dense_qp, struct s_cond_qcqp_arg* cond_arg, struct s_cond_qcqp_ws* cond_ws);
//
void s_cond_qcqp_cond_lhs(struct s_ocp_qcqp* ocp_qp, struct s_dense_qcqp* dense_qp, struct s_cond_qcqp_arg* cond_arg, struct s_cond_qcqp_ws* cond_ws);
//
void s_cond_qcqp_expand_sol(struct s_ocp_qcqp* ocp_qp, struct s_dense_qcqp_sol* dense_qp_sol, struct s_ocp_qcqp_sol* ocp_qp_sol, struct s_cond_qcqp_arg* cond_arg, struct s_cond_qcqp_ws* cond_ws);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif  // HPIPM_S_COND_QCQP_H_
