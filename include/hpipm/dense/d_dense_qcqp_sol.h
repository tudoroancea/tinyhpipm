#ifndef HPIPM_D_DENSE_QCQP_SOL_H_
#define HPIPM_D_DENSE_QCQP_SOL_H_

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/dense/d_dense_qcqp_dim.h"


#ifdef __cplusplus
extern "C" {
#endif


struct d_dense_qcqp_sol {
    struct d_dense_qcqp_dim* dim;
    struct blasfeo_dvec* v;
    struct blasfeo_dvec* pi;
    struct blasfeo_dvec* lam;
    struct blasfeo_dvec* t;
    void* misc;
    hpipm_size_t memsize;
};


//
hpipm_size_t d_dense_qcqp_sol_memsize(struct d_dense_qcqp_dim* dim);
//
void d_dense_qcqp_sol_create(struct d_dense_qcqp_dim* dim, struct d_dense_qcqp_sol* qp_sol, void* memory);
//
void d_dense_qcqp_sol_get(char* field, struct d_dense_qcqp_sol* qp_sol, void* value);
//
void d_dense_qcqp_sol_get_v(struct d_dense_qcqp_sol* qp_sol, double* v);
//
void d_dense_qcqp_sol_get_pi(struct d_dense_qcqp_sol* qp_sol, double* pi);
//
void d_dense_qcqp_sol_get_lam_lb(struct d_dense_qcqp_sol* qp_sol, double* lam_lb);
//
void d_dense_qcqp_sol_get_lam_ub(struct d_dense_qcqp_sol* qp_sol, double* lam_ub);
//
void d_dense_qcqp_sol_get_lam_lg(struct d_dense_qcqp_sol* qp_sol, double* lam_lg);
//
void d_dense_qcqp_sol_get_lam_ug(struct d_dense_qcqp_sol* qp_sol, double* lam_ug);
//
void d_dense_qcqp_sol_get_lam_uq(struct d_dense_qcqp_sol* qp_sol, double* lam_uq);
//
void d_dense_qcqp_sol_set_v(double* v, struct d_dense_qcqp_sol* qp_sol);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif  // HPIPM_D_DENSE_QCQP_SOL_H_
