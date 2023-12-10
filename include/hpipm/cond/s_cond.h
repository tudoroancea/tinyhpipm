#ifndef HPIPM_S_COND_H_
#define HPIPM_S_COND_H_

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_target.h"

#include "hpipm/dense/s_dense_qp.h"
#include "hpipm/dense/s_dense_qp_sol.h"
#include "hpipm/ocp/s_ocp_qp.h"
#include "hpipm/ocp/s_ocp_qp_dim.h"
#include "hpipm/ocp/s_ocp_qp_sol.h"


#ifdef __cplusplus
extern "C" {
#endif


struct s_cond_qp_arg {
    int cond_last_stage;  // condense last stage
    int cond_alg;  // condensing algorithm: 0 N2-nx3, 1 N3-nx2
    int comp_prim_sol;  // primal solution (v)
    int comp_dual_sol_eq;  // dual solution equality constr (pi)
    int comp_dual_sol_ineq;  // dual solution inequality constr (lam t)
    int square_root_alg;  // square root algorithm (faster but requires RSQ>0)
    hpipm_size_t memsize;
};


struct s_cond_qp_ws {
    struct blasfeo_smat* Gamma;
    struct blasfeo_smat* GammaQ;
    struct blasfeo_smat* L;
    struct blasfeo_smat* Lx;
    struct blasfeo_smat* AL;
    struct blasfeo_svec* Gammab;
    struct blasfeo_svec* l;
    struct blasfeo_svec* tmp_nbgM;
    struct blasfeo_svec* tmp_nuxM;
    int bs;  // block size
    hpipm_size_t memsize;
};


//
hpipm_size_t s_cond_qp_arg_memsize();
//
void s_cond_qp_arg_create(struct s_cond_qp_arg* cond_arg, void* mem);
//
void s_cond_qp_arg_set_default(struct s_cond_qp_arg* cond_arg);
// condensing algorithm: 0 N2-nx3, 1 N3-nx2
void s_cond_qp_arg_set_cond_alg(int cond_alg, struct s_cond_qp_arg* cond_arg);
// set riccati-like algorithm: 0 classical, 1 square-root
void s_cond_qp_arg_set_ric_alg(int ric_alg, struct s_cond_qp_arg* cond_arg);
// condense last stage: 0 last stage disregarded, 1 last stage condensed too
void s_cond_qp_arg_set_cond_last_stage(int cond_last_stage, struct s_cond_qp_arg* cond_arg);
//
void s_cond_qp_arg_set_comp_prim_sol(int value, struct s_cond_qp_arg* cond_arg);
//
void s_cond_qp_arg_set_comp_dual_sol_eq(int value, struct s_cond_qp_arg* cond_arg);
//
void s_cond_qp_arg_set_comp_dual_sol_ineq(int value, struct s_cond_qp_arg* cond_arg);

//
void s_cond_qp_compute_dim(struct s_ocp_qp_dim* ocp_dim, struct s_dense_qp_dim* dense_dim);
//
hpipm_size_t s_cond_qp_ws_memsize(struct s_ocp_qp_dim* ocp_dim, struct s_cond_qp_arg* cond_arg);
//
void s_cond_qp_ws_create(struct s_ocp_qp_dim* ocp_dim, struct s_cond_qp_arg* cond_arg, struct s_cond_qp_ws* cond_ws, void* mem);
//
void s_cond_qp_cond(struct s_ocp_qp* ocp_qp, struct s_dense_qp* dense_qp, struct s_cond_qp_arg* cond_arg, struct s_cond_qp_ws* cond_ws);
//
void s_cond_qp_cond_lhs(struct s_ocp_qp* ocp_qp, struct s_dense_qp* dense_qp, struct s_cond_qp_arg* cond_arg, struct s_cond_qp_ws* cond_ws);
//
void s_cond_qp_cond_rhs(struct s_ocp_qp* ocp_qp, struct s_dense_qp* dense_qp, struct s_cond_qp_arg* cond_arg, struct s_cond_qp_ws* cond_ws);
//
void s_cond_qp_expand_sol(struct s_ocp_qp* ocp_qp, struct s_dense_qp_sol* dense_qp_sol, struct s_ocp_qp_sol* ocp_qp_sol, struct s_cond_qp_arg* cond_arg, struct s_cond_qp_ws* cond_ws);
// TODO remove
void s_cond_qp_expand_primal_sol(struct s_ocp_qp* ocp_qp, struct s_dense_qp_sol* dense_qp_sol, struct s_ocp_qp_sol* ocp_qp_sol, struct s_cond_qp_arg* cond_arg, struct s_cond_qp_ws* cond_ws);

//
void s_cond_qp_update(int* idxc, struct s_ocp_qp* ocp_qp, struct s_dense_qp* dense_qp, struct s_cond_qp_arg* cond_arg, struct s_cond_qp_ws* cond_ws);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif  // HPIPM_S_COND_H_
