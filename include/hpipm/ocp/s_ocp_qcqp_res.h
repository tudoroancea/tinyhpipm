#ifndef HPIPM_S_OCP_QCQP_RES_H_
#define HPIPM_S_OCP_QCQP_RES_H_

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_target.h"

#include "hpipm/common.h"
#include "hpipm/ocp/s_ocp_qcqp.h"
#include "hpipm/ocp/s_ocp_qcqp_dim.h"
#include "hpipm/ocp/s_ocp_qcqp_sol.h"


#ifdef __cplusplus
extern "C" {
#endif


struct s_ocp_qcqp_res {
    struct s_ocp_qcqp_dim* dim;
    struct blasfeo_svec* res_g;  // q-residuals
    struct blasfeo_svec* res_b;  // b-residuals
    struct blasfeo_svec* res_d;  // d-residuals
    struct blasfeo_svec* res_m;  // m-residuals
    float res_max[4];  // max of residuals
    float res_mu;  // mu-residual
    float obj;  // (primal) objective
    hpipm_size_t memsize;
};


struct s_ocp_qcqp_res_ws {
    struct blasfeo_svec* tmp_nuxM;  // work space of size nuM+nxM
    struct blasfeo_svec* tmp_nbgqM;  // work space of size nbM+ngM+nqM
    struct blasfeo_svec* tmp_nsM;  // work space of size nsM
    struct blasfeo_svec* q_fun;  // value for evaluation of quadr constr
    struct blasfeo_svec* q_adj;  // value for adjoint of quadr constr
    int use_q_fun;  // reuse cached value for evaluation of quadr constr
    int use_q_adj;  // reuse cached value for adjoint of quadr constr
    hpipm_size_t memsize;
};


//
hpipm_size_t s_ocp_qcqp_res_memsize(struct s_ocp_qcqp_dim* ocp_dim);
//
void s_ocp_qcqp_res_create(struct s_ocp_qcqp_dim* ocp_dim, struct s_ocp_qcqp_res* res, void* mem);
//
hpipm_size_t s_ocp_qcqp_res_ws_memsize(struct s_ocp_qcqp_dim* ocp_dim);
//
void s_ocp_qcqp_res_ws_create(struct s_ocp_qcqp_dim* ocp_dim, struct s_ocp_qcqp_res_ws* workspace, void* mem);
//
void s_ocp_qcqp_res_compute(struct s_ocp_qcqp* qp, struct s_ocp_qcqp_sol* qp_sol, struct s_ocp_qcqp_res* res, struct s_ocp_qcqp_res_ws* ws);
//
void s_ocp_qcqp_res_compute_inf_norm(struct s_ocp_qcqp_res* res);
//
void s_ocp_qcqp_res_get_max_res_stat(struct s_ocp_qcqp_res* res, float* value);
//
void s_ocp_qcqp_res_get_max_res_eq(struct s_ocp_qcqp_res* res, float* value);
//
void s_ocp_qcqp_res_get_max_res_ineq(struct s_ocp_qcqp_res* res, float* value);
//
void s_ocp_qcqp_res_get_max_res_comp(struct s_ocp_qcqp_res* res, float* value);


#ifdef __cplusplus
}  // #extern "C"
#endif


#endif  // HPIPM_S_OCP_QCQP_RES_H_
