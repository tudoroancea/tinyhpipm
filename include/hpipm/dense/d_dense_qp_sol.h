#ifndef HPIPM_D_DENSE_QP_SOL_H_
#define HPIPM_D_DENSE_QP_SOL_H_

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_target.h"

#include "hpipm/dense/d_dense_qp_dim.h"


#ifdef __cplusplus
extern "C" {
#endif


struct d_dense_qp_sol {
    struct d_dense_qp_dim* dim;
    struct blasfeo_dvec* v;
    struct blasfeo_dvec* pi;
    struct blasfeo_dvec* lam;
    struct blasfeo_dvec* t;
    void* misc;
    double obj;
    int valid_obj;
    hpipm_size_t memsize;
};


//
hpipm_size_t d_dense_qp_sol_strsize();
//
hpipm_size_t d_dense_qp_sol_memsize(struct d_dense_qp_dim* dim);
//
void d_dense_qp_sol_create(struct d_dense_qp_dim* dim, struct d_dense_qp_sol* qp_sol, void* memory);
//
void d_dense_qp_sol_get_all(struct d_dense_qp_sol* qp_sol, double* v, double* ls, double* us, double* pi, double* lam_lb, double* lam_ub, double* lam_lg, double* lam_ug, double* lam_ls, double* lam_us);
//
void d_dense_qp_sol_get(char* field, struct d_dense_qp_sol* sol, void* value);
//
void d_dense_qp_sol_get_v(struct d_dense_qp_sol* sol, double* v);
//
void d_dense_qp_sol_get_pi(struct d_dense_qp_sol* sol, double* pi);
//
void d_dense_qp_sol_get_lam_lb(struct d_dense_qp_sol* sol, double* lam_lb);
//
void d_dense_qp_sol_get_lam_ub(struct d_dense_qp_sol* sol, double* lam_ub);
//
void d_dense_qp_sol_get_lam_lg(struct d_dense_qp_sol* sol, double* lam_lg);
//
void d_dense_qp_sol_get_lam_ug(struct d_dense_qp_sol* sol, double* lam_ug);
//
void d_dense_qp_sol_get_valid_obj(struct d_dense_qp_sol* sol, int* valid_obj);
//
void d_dense_qp_sol_get_obj(struct d_dense_qp_sol* sol, double* obj);
//
void d_dense_qp_sol_set_v(double* v, struct d_dense_qp_sol* sol);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif  // HPIPM_D_DENSE_QP_SOL_H_
