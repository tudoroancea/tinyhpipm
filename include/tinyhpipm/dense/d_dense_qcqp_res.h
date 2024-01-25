#ifndef HPIPM_D_d_dense_qcqp_res_H_
#define HPIPM_D_d_dense_qcqp_res_H_

#include "tinyhpipm/blas.h"
#include "tinyhpipm/common.h"
#include "tinyhpipm/dense/d_dense_qcqp.h"
#include "tinyhpipm/dense/d_dense_qcqp_dim.h"
#include "tinyhpipm/dense/d_dense_qcqp_sol.h"


#ifdef __cplusplus
extern "C" {
#endif


struct d_dense_qcqp_res {
    struct d_dense_qcqp_dim* dim;
    struct vec* res_g;  // q-residuals
    struct vec* res_b;  // b-residuals
    struct vec* res_d;  // d-residuals
    struct vec* res_m;  // m-residuals
    double res_max[4];  // infinity norm of residuals
    double res_mu;  // mu-residual
    double obj;  // (primal) objective
    hpipm_size_t memsize;
};


struct d_dense_qcqp_res_ws {
    struct vec* tmp_nv;  // work space of size nv
    struct vec* tmp_nbgq;  // work space of size nbM+ngM+nqM
    struct vec* tmp_ns;  // work space of size nsM
    struct vec* q_fun;  // value for evaluation of quadr constr
    struct vec* q_adj;  // value for adjoint of quadr constr
    int use_q_fun;  // reuse cached value for evaluation of quadr constr
    int use_q_adj;  // reuse cached value for adjoint of quadr constr
    hpipm_size_t memsize;
};


//
hpipm_size_t d_dense_qcqp_res_memsize(struct d_dense_qcqp_dim* dim);
//
void d_dense_qcqp_res_create(struct d_dense_qcqp_dim* dim, struct d_dense_qcqp_res* res, void* mem);
//
hpipm_size_t d_dense_qcqp_res_ws_memsize(struct d_dense_qcqp_dim* dim);
//
void d_dense_qcqp_res_ws_create(struct d_dense_qcqp_dim* dim, struct d_dense_qcqp_res_ws* workspace, void* mem);
//
void d_dense_qcqp_res_compute(struct d_dense_qcqp* qp, struct d_dense_qcqp_sol* qp_sol, struct d_dense_qcqp_res* res, struct d_dense_qcqp_res_ws* ws);
//
void d_dense_qcqp_res_compute_inf_norm(struct d_dense_qcqp_res* res);


#ifdef __cplusplus
}  // #extern "C"
#endif


#endif  // HPIPM_D_d_dense_qcqp_res_H_
