#ifndef HPIPM_D_OCP_qcqp_RED_H_
#define HPIPM_D_OCP_qcqp_RED_H_

#include "hpipm/blas.h"
#include "hpipm/ocp/d_ocp_qcqp_dim.h"
#include "hpipm/ocp/d_ocp_qcqp_sol.h"
#include "hpipm/ocp/d_ocp_qp.h"


#ifdef __cplusplus
extern "C" {
#endif


struct d_ocp_qcqp_reduce_eq_dof_arg {
    double lam_min;
    double t_min;
    int alias_unchanged;  // do not keep copy unchanged stage
    int comp_prim_sol;  // primal solution (v)
    int comp_dual_sol_eq;  // dual solution equality constr (pi)
    int comp_dual_sol_ineq;  // dual solution inequality constr (lam t)
    hpipm_size_t memsize;  // memory size in bytes
};


struct d_ocp_qcqp_reduce_eq_dof_ws {
    struct vec* tmp_nuxM;
    struct vec* tmp_nbgqM;
    int* e_imask_ux;
    int* e_imask_d;
    hpipm_size_t memsize;  // memory size in bytes
};


//
void d_ocp_qcqp_dim_reduce_eq_dof(struct d_ocp_qcqp_dim* dim, struct d_ocp_qcqp_dim* dim_red);
//
hpipm_size_t d_ocp_qcqp_reduce_eq_dof_arg_memsize();
//
void d_ocp_qcqp_reduce_eq_dof_arg_create(struct d_ocp_qcqp_reduce_eq_dof_arg* arg, void* mem);
//
void d_ocp_qcqp_reduce_eq_dof_arg_set_default(struct d_ocp_qcqp_reduce_eq_dof_arg* arg);
//
void d_ocp_qcqp_reduce_eq_dof_arg_set_alias_unchanged(struct d_ocp_qcqp_reduce_eq_dof_arg* arg, int value);
//
void d_ocp_qcqp_reduce_eq_dof_arg_set_comp_prim_sol(struct d_ocp_qcqp_reduce_eq_dof_arg* arg, int value);
//
void d_ocp_qcqp_reduce_eq_dof_arg_set_comp_dual_sol_eq(struct d_ocp_qcqp_reduce_eq_dof_arg* arg, int value);
//
void d_ocp_qcqp_reduce_eq_dof_arg_set_comp_dual_sol_ineq(struct d_ocp_qcqp_reduce_eq_dof_arg* arg, int value);
//
hpipm_size_t d_ocp_qcqp_reduce_eq_dof_ws_memsize(struct d_ocp_qcqp_dim* dim);
//
void d_ocp_qcqp_reduce_eq_dof_ws_create(struct d_ocp_qcqp_dim* dim, struct d_ocp_qcqp_reduce_eq_dof_ws* work, void* mem);
//
void d_ocp_qcqp_reduce_eq_dof(struct d_ocp_qcqp* qp, struct d_ocp_qcqp* qp_red, struct d_ocp_qcqp_reduce_eq_dof_arg* arg, struct d_ocp_qcqp_reduce_eq_dof_ws* work);
//
void d_ocp_qcqp_reduce_eq_dof_lhs(struct d_ocp_qcqp* qp, struct d_ocp_qcqp* qp_red, struct d_ocp_qcqp_reduce_eq_dof_arg* arg, struct d_ocp_qcqp_reduce_eq_dof_ws* work);
//
void d_ocp_qcqp_reduce_eq_dof_rhs(struct d_ocp_qcqp* qp, struct d_ocp_qcqp* qp_red, struct d_ocp_qcqp_reduce_eq_dof_arg* arg, struct d_ocp_qcqp_reduce_eq_dof_ws* work);
//
void d_ocp_qcqp_restore_eq_dof(struct d_ocp_qcqp* qp, struct d_ocp_qcqp_sol* qp_sol_red, struct d_ocp_qcqp_sol* qp_sol, struct d_ocp_qcqp_reduce_eq_dof_arg* arg, struct d_ocp_qcqp_reduce_eq_dof_ws* work);


#ifdef __cplusplus
}  // #extern "C"
#endif


#endif  // HPIPM_D_OCP_qcqp_RED_H_
