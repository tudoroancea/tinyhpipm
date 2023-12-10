#ifndef HPIPM_D_DENSE_QCQP_RES_H_
#define HPIPM_D_DENSE_QCQP_RES_H_

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/dense/d_dense_qcqp.h"
#include "hpipm/dense/d_dense_qcqp_dim.h"
#include "hpipm/dense/d_dense_qcqp_sol.h"


#ifdef __cplusplus
extern "C" {
#endif


struct d_dense_qcqp_res {
    struct d_dense_qcqp_dim* dim;
    struct blasfeo_dvec* res_g;  // q-residuals
    struct blasfeo_dvec* res_b;  // b-residuals
    struct blasfeo_dvec* res_d;  // d-residuals
    struct blasfeo_dvec* res_m;  // m-residuals
    double res_max[4];  // infinity norm of residuals
    double res_mu;  // mu-residual
    double obj;  // (primal) objective
    hpipm_size_t memsize;
};


struct d_dense_qcqp_res_ws {
    struct blasfeo_dvec* tmp_nv;  // work space of size nv
    struct blasfeo_dvec* tmp_nbgq;  // work space of size nbM+ngM+nqM
    struct blasfeo_dvec* tmp_ns;  // work space of size nsM
    struct blasfeo_dvec* q_fun;  // value for evaluation of quadr constr
    struct blasfeo_dvec* q_adj;  // value for adjoint of quadr constr
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


#endif  // HPIPM_D_DENSE_QCQP_RES_H_
