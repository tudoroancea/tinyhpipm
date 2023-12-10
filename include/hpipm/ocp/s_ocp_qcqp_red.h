#ifndef HPIPM_S_OCP_QCQP_RED_H_
#define HPIPM_S_OCP_QCQP_RED_H_

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_target.h"

#include "hpipm/ocp/s_ocp_qcqp_dim.h"
#include "hpipm/ocp/s_ocp_qcqp_sol.h"
#include "hpipm/ocp/s_ocp_qp.h"


#ifdef __cplusplus
extern "C" {
#endif


struct s_ocp_qcqp_reduce_eq_dof_arg {
    float lam_min;
    float t_min;
    int alias_unchanged;  // do not keep copy unchanged stage
    int comp_prim_sol;  // primal solution (v)
    int comp_dual_sol_eq;  // dual solution equality constr (pi)
    int comp_dual_sol_ineq;  // dual solution inequality constr (lam t)
    hpipm_size_t memsize;  // memory size in bytes
};


struct s_ocp_qcqp_reduce_eq_dof_ws {
    struct blasfeo_svec* tmp_nuxM;
    struct blasfeo_svec* tmp_nbgqM;
    int* e_imask_ux;
    int* e_imask_d;
    hpipm_size_t memsize;  // memory size in bytes
};


//
void s_ocp_qcqp_dim_reduce_eq_dof(struct s_ocp_qcqp_dim* dim, struct s_ocp_qcqp_dim* dim_red);
//
hpipm_size_t s_ocp_qcqp_reduce_eq_dof_arg_memsize();
//
void s_ocp_qcqp_reduce_eq_dof_arg_create(struct s_ocp_qcqp_reduce_eq_dof_arg* arg, void* mem);
//
void s_ocp_qcqp_reduce_eq_dof_arg_set_default(struct s_ocp_qcqp_reduce_eq_dof_arg* arg);
//
void s_ocp_qcqp_reduce_eq_dof_arg_set_alias_unchanged(struct s_ocp_qcqp_reduce_eq_dof_arg* arg, int value);
//
void s_ocp_qcqp_reduce_eq_dof_arg_set_comp_prim_sol(struct s_ocp_qcqp_reduce_eq_dof_arg* arg, int value);
//
void s_ocp_qcqp_reduce_eq_dof_arg_set_comp_dual_sol_eq(struct s_ocp_qcqp_reduce_eq_dof_arg* arg, int value);
//
void s_ocp_qcqp_reduce_eq_dof_arg_set_comp_dual_sol_ineq(struct s_ocp_qcqp_reduce_eq_dof_arg* arg, int value);
//
hpipm_size_t s_ocp_qcqp_reduce_eq_dof_ws_memsize(struct s_ocp_qcqp_dim* dim);
//
void s_ocp_qcqp_reduce_eq_dof_ws_create(struct s_ocp_qcqp_dim* dim, struct s_ocp_qcqp_reduce_eq_dof_ws* work, void* mem);
//
void s_ocp_qcqp_reduce_eq_dof(struct s_ocp_qcqp* qp, struct s_ocp_qcqp* qp_red, struct s_ocp_qcqp_reduce_eq_dof_arg* arg, struct s_ocp_qcqp_reduce_eq_dof_ws* work);
//
void s_ocp_qcqp_reduce_eq_dof_lhs(struct s_ocp_qcqp* qp, struct s_ocp_qcqp* qp_red, struct s_ocp_qcqp_reduce_eq_dof_arg* arg, struct s_ocp_qcqp_reduce_eq_dof_ws* work);
//
void s_ocp_qcqp_reduce_eq_dof_rhs(struct s_ocp_qcqp* qp, struct s_ocp_qcqp* qp_red, struct s_ocp_qcqp_reduce_eq_dof_arg* arg, struct s_ocp_qcqp_reduce_eq_dof_ws* work);
//
void s_ocp_qcqp_restore_eq_dof(struct s_ocp_qcqp* qp, struct s_ocp_qcqp_sol* qp_sol_red, struct s_ocp_qcqp_sol* qp_sol, struct s_ocp_qcqp_reduce_eq_dof_arg* arg, struct s_ocp_qcqp_reduce_eq_dof_ws* work);


#ifdef __cplusplus
}  // #extern "C"
#endif


#endif  // HPIPM_S_OCP_QCQP_RED_H_
