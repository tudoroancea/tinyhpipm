#ifndef HPIPM_S_DENSE_QCQP_RES_H_
#define HPIPM_S_DENSE_QCQP_RES_H_

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/dense/s_dense_qcqp.h"
#include "hpipm/dense/s_dense_qcqp_dim.h"
#include "hpipm/dense/s_dense_qcqp_sol.h"


#ifdef __cplusplus
extern "C" {
#endif


struct s_dense_qcqp_res {
    struct s_dense_qcqp_dim* dim;
    struct blasfeo_svec* res_g;  // q-residuals
    struct blasfeo_svec* res_b;  // b-residuals
    struct blasfeo_svec* res_d;  // d-residuals
    struct blasfeo_svec* res_m;  // m-residuals
    float res_max[4];  // infinity norm of residuals
    float res_mu;  // mu-residual
    float obj;  // (primal) objective
    hpipm_size_t memsize;
};


struct s_dense_qcqp_res_ws {
    struct blasfeo_svec* tmp_nv;  // work space of size nv
    struct blasfeo_svec* tmp_nbgq;  // work space of size nbM+ngM+nqM
    struct blasfeo_svec* tmp_ns;  // work space of size nsM
    struct blasfeo_svec* q_fun;  // value for evaluation of quadr constr
    struct blasfeo_svec* q_adj;  // value for adjoint of quadr constr
    int use_q_fun;  // reuse cached value for evaluation of quadr constr
    int use_q_adj;  // reuse cached value for adjoint of quadr constr
    hpipm_size_t memsize;
};


//
hpipm_size_t s_dense_qcqp_res_memsize(struct s_dense_qcqp_dim* dim);
//
void s_dense_qcqp_res_create(struct s_dense_qcqp_dim* dim, struct s_dense_qcqp_res* res, void* mem);
//
hpipm_size_t s_dense_qcqp_res_ws_memsize(struct s_dense_qcqp_dim* dim);
//
void s_dense_qcqp_res_ws_create(struct s_dense_qcqp_dim* dim, struct s_dense_qcqp_res_ws* workspace, void* mem);
//
void s_dense_qcqp_res_compute(struct s_dense_qcqp* qp, struct s_dense_qcqp_sol* qp_sol, struct s_dense_qcqp_res* res, struct s_dense_qcqp_res_ws* ws);
//
void s_dense_qcqp_res_compute_inf_norm(struct s_dense_qcqp_res* res);


#ifdef __cplusplus
}  // #extern "C"
#endif


#endif  // HPIPM_S_DENSE_QCQP_RES_H_
