#ifndef HPIPM_D_PART_COND_H_
#define HPIPM_D_PART_COND_H_

#include "tinyhpipm/blas.h"
#include "tinyhpipm/common.h"

#include "tinyhpipm/common.h"
#include "tinyhpipm/cond/d_cond.h"


#ifdef __cplusplus
extern "C" {
#endif


struct d_part_cond_qp_arg {
    struct d_cond_qp_arg* cond_arg;
    int N2;
    hpipm_size_t memsize;
};


struct d_part_cond_qp_ws {
    struct d_cond_qp_ws* cond_workspace;
    hpipm_size_t memsize;
};


//
hpipm_size_t d_part_cond_qp_arg_memsize(int N2);
//
void d_part_cond_qp_arg_create(int N2, struct d_part_cond_qp_arg* cond_arg, void* mem);
//
void d_part_cond_qp_arg_set_default(struct d_part_cond_qp_arg* cond_arg);
// set riccati-like algorithm: 0 classical, 1 squre-root
void d_part_cond_qp_arg_set_ric_alg(int ric_alg, struct d_part_cond_qp_arg* cond_arg);
//
void d_part_cond_qp_arg_set_comp_prim_sol(int value, struct d_part_cond_qp_arg* cond_arg);
//
void d_part_cond_qp_arg_set_comp_dual_sol_eq(int value, struct d_part_cond_qp_arg* cond_arg);
//
void d_part_cond_qp_arg_set_comp_dual_sol_ineq(int value, struct d_part_cond_qp_arg* cond_arg);

//
void d_part_cond_qp_compute_block_size(int N, int N2, int* block_size);
//
void d_part_cond_qp_compute_dim(struct d_ocp_qp_dim* ocp_dim, int* block_size, struct d_ocp_qp_dim* part_dense_dim);
//
hpipm_size_t d_part_cond_qp_ws_memsize(struct d_ocp_qp_dim* ocp_dim, int* block_size, struct d_ocp_qp_dim* part_dense_dim, struct d_part_cond_qp_arg* cond_arg);
//
void d_part_cond_qp_ws_create(struct d_ocp_qp_dim* ocp_dim, int* block_size, struct d_ocp_qp_dim* part_dense_dim, struct d_part_cond_qp_arg* cond_arg, struct d_part_cond_qp_ws* cond_ws, void* mem);
//
void d_part_cond_qp_cond(struct d_ocp_qp* ocp_qp, struct d_ocp_qp* part_dense_qp, struct d_part_cond_qp_arg* cond_arg, struct d_part_cond_qp_ws* cond_ws);
//
void d_part_cond_qp_cond_lhs(struct d_ocp_qp* ocp_qp, struct d_ocp_qp* part_dense_qp, struct d_part_cond_qp_arg* cond_arg, struct d_part_cond_qp_ws* cond_ws);
//
void d_part_cond_qp_cond_rhs(struct d_ocp_qp* ocp_qp, struct d_ocp_qp* part_dense_qp, struct d_part_cond_qp_arg* cond_arg, struct d_part_cond_qp_ws* cond_ws);
//
void d_part_cond_qp_expand_sol(struct d_ocp_qp* ocp_qp, struct d_ocp_qp* part_dense_qp, struct d_ocp_qp_sol* part_dense_qp_sol, struct d_ocp_qp_sol* ocp_qp_sol, struct d_part_cond_qp_arg* cond_arg, struct d_part_cond_qp_ws* cond_ws);

//
void d_part_cond_qp_update(int* idxc, struct d_ocp_qp* ocp_qp, struct d_ocp_qp* part_dense_qp, struct d_part_cond_qp_arg* cond_arg, struct d_part_cond_qp_ws* cond_ws);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif  // HPIPM_D_PART_COND_H_
