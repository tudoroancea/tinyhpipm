#ifndef HPIPM_S_DENSE_QP_SOL_H_
#define HPIPM_S_DENSE_QP_SOL_H_

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/dense/s_dense_qp_dim.h"


#ifdef __cplusplus
extern "C" {
#endif


struct s_dense_qp_sol {
    struct s_dense_qp_dim* dim;
    struct blasfeo_svec* v;
    struct blasfeo_svec* pi;
    struct blasfeo_svec* lam;
    struct blasfeo_svec* t;
    void* misc;
    float obj;
    int valid_obj;
    hpipm_size_t memsize;
};


//
hpipm_size_t s_dense_qp_sol_strsize();
//
hpipm_size_t s_dense_qp_sol_memsize(struct s_dense_qp_dim* dim);
//
void s_dense_qp_sol_create(struct s_dense_qp_dim* dim, struct s_dense_qp_sol* qp_sol, void* memory);
//
void s_dense_qp_sol_get_all(struct s_dense_qp_sol* qp_sol, float* v, float* ls, float* us, float* pi, float* lam_lb, float* lam_ub, float* lam_lg, float* lam_ug, float* lam_ls, float* lam_us);
//
void s_dense_qp_sol_get(char* field, struct s_dense_qp_sol* sol, void* value);
//
void s_dense_qp_sol_get_v(struct s_dense_qp_sol* sol, float* v);
//
void s_dense_qp_sol_get_pi(struct s_dense_qp_sol* sol, float* pi);
//
void s_dense_qp_sol_get_lam_lb(struct s_dense_qp_sol* sol, float* lam_lb);
//
void s_dense_qp_sol_get_lam_ub(struct s_dense_qp_sol* sol, float* lam_ub);
//
void s_dense_qp_sol_get_lam_lg(struct s_dense_qp_sol* sol, float* lam_lg);
//
void s_dense_qp_sol_get_lam_ug(struct s_dense_qp_sol* sol, float* lam_ug);
//
void s_dense_qp_sol_get_valid_obj(struct s_dense_qp_sol* sol, int* valid_obj);
//
void s_dense_qp_sol_get_obj(struct s_dense_qp_sol* sol, float* obj);

//
void s_dense_qp_sol_set_v(float* v, struct s_dense_qp_sol* sol);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif  // HPIPM_S_DENSE_QP_SOL_H_
