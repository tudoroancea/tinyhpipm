#ifndef HPIPM_D_OCP_qcqp_res_H_
#define HPIPM_D_OCP_qcqp_res_H_

#include "hpipm/blas.h"
#include "hpipm/common.h"
#include "hpipm/ocp/d_ocp_qcqp.h"
#include "hpipm/ocp/d_ocp_qcqp_dim.h"
#include "hpipm/ocp/d_ocp_qcqp_sol.h"


#ifdef __cplusplus
extern "C" {
#endif


struct d_ocp_qcqp_res {
    struct d_ocp_qcqp_dim* dim;
    struct vec* res_g;  // q-residuals
    struct vec* res_b;  // b-residuals
    struct vec* res_d;  // d-residuals
    struct vec* res_m;  // m-residuals
    double res_max[4];  // max of residuals
    double res_mu;  // mu-residual
    double obj;  // (primal) objective
    hpipm_size_t memsize;
};


struct d_ocp_qcqp_res_ws {
    struct vec* tmp_nuxM;  // work space of size nuM+nxM
    struct vec* tmp_nbgqM;  // work space of size nbM+ngM+nqM
    struct vec* tmp_nsM;  // work space of size nsM
    struct vec* q_fun;  // value for evaluation of quadr constr
    struct vec* q_adj;  // value for adjoint of quadr constr
    int use_q_fun;  // reuse cached value for evaluation of quadr constr
    int use_q_adj;  // reuse cached value for adjoint of quadr constr
    hpipm_size_t memsize;
};


//
hpipm_size_t d_ocp_qcqp_res_memsize(struct d_ocp_qcqp_dim* ocp_dim);
//
void d_ocp_qcqp_res_create(struct d_ocp_qcqp_dim* ocp_dim, struct d_ocp_qcqp_res* res, void* mem);
//
hpipm_size_t d_ocp_qcqp_res_ws_memsize(struct d_ocp_qcqp_dim* ocp_dim);
//
void d_ocp_qcqp_res_ws_create(struct d_ocp_qcqp_dim* ocp_dim, struct d_ocp_qcqp_res_ws* workspace, void* mem);
//
void d_ocp_qcqp_res_compute(struct d_ocp_qcqp* qp, struct d_ocp_qcqp_sol* qp_sol, struct d_ocp_qcqp_res* res, struct d_ocp_qcqp_res_ws* ws);
//
void d_ocp_qcqp_res_compute_inf_norm(struct d_ocp_qcqp_res* res);
//
void d_ocp_qcqp_res_get_max_res_stat(struct d_ocp_qcqp_res* res, double* value);
//
void d_ocp_qcqp_res_get_max_res_eq(struct d_ocp_qcqp_res* res, double* value);
//
void d_ocp_qcqp_res_get_max_res_ineq(struct d_ocp_qcqp_res* res, double* value);
//
void d_ocp_qcqp_res_get_max_res_comp(struct d_ocp_qcqp_res* res, double* value);


#ifdef __cplusplus
}  // #extern "C"
#endif


#endif  // HPIPM_D_OCP_qcqp_res_H_
