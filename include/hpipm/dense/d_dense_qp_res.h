#ifndef HPIPM_D_d_dense_qp_res_H_
#define HPIPM_D_d_dense_qp_res_H_


#include "hpipm/blas.h"
#include "hpipm/dense/d_dense_qp.h"
#include "hpipm/dense/d_dense_qp_dim.h"
#include "hpipm/dense/d_dense_qp_sol.h"


#ifdef __cplusplus
extern "C" {
#endif


struct d_dense_qp_res {
    struct d_dense_qp_dim* dim;
    struct vec* res_g;  // q-residuals
    struct vec* res_b;  // b-residuals
    struct vec* res_d;  // d-residuals
    struct vec* res_m;  // m-residuals
    double res_max[4];  // max of residuals
    double res_mu;  // mu-residual
    double obj;  // (primal) objective
    hpipm_size_t memsize;
};


struct d_dense_qp_res_ws {
    struct vec* tmp_nbg;  // work space of size nbM+ngM
    struct vec* tmp_ns;  // work space of size nsM
    hpipm_size_t memsize;
};


//
hpipm_size_t d_dense_qp_res_memsize(struct d_dense_qp_dim* dim);
//
void d_dense_qp_res_create(struct d_dense_qp_dim* dim, struct d_dense_qp_res* res, void* mem);
//
hpipm_size_t d_dense_qp_res_ws_memsize(struct d_dense_qp_dim* dim);
//
void d_dense_qp_res_ws_create(struct d_dense_qp_dim* dim, struct d_dense_qp_res_ws* workspace, void* mem);
//
void d_dense_qp_res_compute(struct d_dense_qp* qp, struct d_dense_qp_sol* qp_sol, struct d_dense_qp_res* res, struct d_dense_qp_res_ws* ws);
//
void d_dense_qp_res_compute_lin(struct d_dense_qp* qp, struct d_dense_qp_sol* qp_sol, struct d_dense_qp_sol* qp_step, struct d_dense_qp_res* res, struct d_dense_qp_res_ws* ws);
//
void d_dense_qp_res_compute_inf_norm(struct d_dense_qp_res* res);
//
void d_dense_qp_res_get_all(struct d_dense_qp_res* res, double* res_g, double* res_ls, double* res_us, double* res_b, double* res_d_lb, double* res_d_ub, double* res_d_lg, double* res_d_ug, double* res_d_ls, double* res_d_us, double* res_m_lb, double* res_m_ub, double* res_m_lg, double* res_m_ug, double* res_m_ls, double* res_m_us);


#ifdef __cplusplus
}  // #extern "C"
#endif


#endif  // HPIPM_D_d_dense_qp_res_H_
