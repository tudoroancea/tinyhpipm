#ifndef HPIPM_S_OCP_QP_RES_H_
#define HPIPM_S_OCP_QP_RES_H_

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/common.h"
#include "hpipm/ocp/s_ocp_qp.h"
#include "hpipm/ocp/s_ocp_qp_dim.h"
#include "hpipm/ocp/s_ocp_qp_sol.h"


#ifdef __cplusplus
extern "C" {
#endif


struct s_ocp_qp_res {
    struct s_ocp_qp_dim* dim;
    struct blasfeo_svec* res_g;  // q-residuals
    struct blasfeo_svec* res_b;  // b-residuals
    struct blasfeo_svec* res_d;  // d-residuals
    struct blasfeo_svec* res_m;  // m-residuals
    float res_max[4];  // max of residuals
    float res_mu;  // mu-residual
    float obj;  // (primal) objective
    hpipm_size_t memsize;
};


struct s_ocp_qp_res_ws {
    struct blasfeo_svec* tmp_nbgM;  // work space of size nbM+ngM
    struct blasfeo_svec* tmp_nsM;  // work space of size nsM
    hpipm_size_t memsize;
};


//
hpipm_size_t s_ocp_qp_res_memsize(struct s_ocp_qp_dim* ocp_dim);
//
void s_ocp_qp_res_create(struct s_ocp_qp_dim* ocp_dim, struct s_ocp_qp_res* res, void* mem);
//
hpipm_size_t s_ocp_qp_res_ws_memsize(struct s_ocp_qp_dim* ocp_dim);
//
void s_ocp_qp_res_ws_create(struct s_ocp_qp_dim* ocp_dim, struct s_ocp_qp_res_ws* workspace, void* mem);
//
void s_ocp_qp_res_compute(struct s_ocp_qp* qp, struct s_ocp_qp_sol* qp_sol, struct s_ocp_qp_res* res, struct s_ocp_qp_res_ws* ws);
//
void s_ocp_qp_res_compute_lin(struct s_ocp_qp* qp, struct s_ocp_qp_sol* qp_sol, struct s_ocp_qp_sol* qp_step, struct s_ocp_qp_res* res, struct s_ocp_qp_res_ws* ws);
//
void s_ocp_qp_res_compute_inf_norm(struct s_ocp_qp_res* res);
//
void s_ocp_qp_res_get_all(struct s_ocp_qp_res* res, float** res_r, float** res_q, float** res_ls, float** res_us, float** res_b, float** res_d_lb, float** res_d_ub, float** res_d_lg, float** res_d_ug, float** res_d_ls, float** res_d_us, float** res_m_lb, float** res_m_ub, float** res_m_lg, float** res_m_ug, float** res_m_ls, float** res_m_us);
//
void s_ocp_qp_res_get_max_res_stat(struct s_ocp_qp_res* res, float* value);
//
void s_ocp_qp_res_get_max_res_eq(struct s_ocp_qp_res* res, float* value);
//
void s_ocp_qp_res_get_max_res_ineq(struct s_ocp_qp_res* res, float* value);
//
void s_ocp_qp_res_get_max_res_comp(struct s_ocp_qp_res* res, float* value);


#ifdef __cplusplus
}  // #extern "C"
#endif


#endif  // HPIPM_S_OCP_QP_RES_H_
