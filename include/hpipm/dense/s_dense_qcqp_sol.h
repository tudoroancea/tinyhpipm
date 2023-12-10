#ifndef HPIPM_S_DENSE_QCQP_SOL_H_
#define HPIPM_S_DENSE_QCQP_SOL_H_

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_target.h"

#include "hpipm/dense/s_dense_qcqp_dim.h"


#ifdef __cplusplus
extern "C" {
#endif


struct s_dense_qcqp_sol {
    struct s_dense_qcqp_dim* dim;
    struct blasfeo_svec* v;
    struct blasfeo_svec* pi;
    struct blasfeo_svec* lam;
    struct blasfeo_svec* t;
    void* misc;
    hpipm_size_t memsize;
};


//
hpipm_size_t s_dense_qcqp_sol_memsize(struct s_dense_qcqp_dim* dim);
//
void s_dense_qcqp_sol_create(struct s_dense_qcqp_dim* dim, struct s_dense_qcqp_sol* qp_sol, void* memory);
//
void s_dense_qcqp_sol_get(char* field, struct s_dense_qcqp_sol* qp_sol, void* value);
//
void s_dense_qcqp_sol_get_v(struct s_dense_qcqp_sol* qp_sol, float* v);
//
void s_dense_qcqp_sol_get_pi(struct s_dense_qcqp_sol* qp_sol, float* pi);
//
void s_dense_qcqp_sol_get_lam_lb(struct s_dense_qcqp_sol* qp_sol, float* lam_lb);
//
void s_dense_qcqp_sol_get_lam_ub(struct s_dense_qcqp_sol* qp_sol, float* lam_ub);
//
void s_dense_qcqp_sol_get_lam_lg(struct s_dense_qcqp_sol* qp_sol, float* lam_lg);
//
void s_dense_qcqp_sol_get_lam_ug(struct s_dense_qcqp_sol* qp_sol, float* lam_ug);
//
void s_dense_qcqp_sol_get_lam_uq(struct s_dense_qcqp_sol* qp_sol, float* lam_uq);
//
void s_dense_qcqp_sol_set_v(float* v, struct s_dense_qcqp_sol* qp_sol);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif  // HPIPM_S_DENSE_QCQP_SOL_H_
