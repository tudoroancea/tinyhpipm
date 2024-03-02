#ifndef HPIPM_D_OCP_qp_res_H_
#define HPIPM_D_OCP_qp_res_H_

#include "tinyhpipm/blas.h"
#include "tinyhpipm/common.h"
#include "tinyhpipm/ocp/d_ocp_qp.h"
#include "tinyhpipm/ocp/d_ocp_qp_dim.h"
#include "tinyhpipm/ocp/d_ocp_qp_sol.h"


#ifdef __cplusplus
extern "C" {
#endif


struct d_ocp_qp_res {
    struct d_ocp_qp_dim* dim;
    struct vec* res_g;  // q-residuals
    struct vec* res_b;  // b-residuals
    struct vec* res_d;  // d-residuals
    struct vec* res_m;  // m-residuals
    double res_max[4];  // max of residuals
    double res_mu;  // mu-residual
    double obj;  // (primal) objective
    hpipm_size_t memsize;
};


struct d_ocp_qp_res_ws {
    struct vec* tmp_nbgM;  // work space of size nbM+ngM
    struct vec* tmp_nsM;  // work space of size nsM
    hpipm_size_t memsize;
};


//
hpipm_size_t d_ocp_qp_res_memsize(struct d_ocp_qp_dim* ocp_dim);
//
void d_ocp_qp_res_create(struct d_ocp_qp_dim* ocp_dim, struct d_ocp_qp_res* res, void* mem);
//
hpipm_size_t d_ocp_qp_res_ws_memsize(struct d_ocp_qp_dim* ocp_dim);
//
void d_ocp_qp_res_ws_create(struct d_ocp_qp_dim* ocp_dim, struct d_ocp_qp_res_ws* workspace, void* mem);
//
void d_ocp_qp_res_compute(struct d_ocp_qp* qp, struct d_ocp_qp_sol* qp_sol, struct d_ocp_qp_res* res, struct d_ocp_qp_res_ws* ws);
//
void d_ocp_qp_res_compute_lin(struct d_ocp_qp* qp, struct d_ocp_qp_sol* qp_sol, struct d_ocp_qp_sol* qp_step, struct d_ocp_qp_res* res, struct d_ocp_qp_res_ws* ws);
//
void d_ocp_qp_res_compute_inf_norm(struct d_ocp_qp_res* res);
//
void d_ocp_qp_res_get_all(struct d_ocp_qp_res* res, double** res_r, double** res_q, double** res_ls, double** res_us, double** res_b, double** res_d_lb, double** res_d_ub, double** res_d_lg, double** res_d_ug, double** res_d_ls, double** res_d_us, double** res_m_lb, double** res_m_ub, double** res_m_lg, double** res_m_ug, double** res_m_ls, double** res_m_us);
//
void d_ocp_qp_res_get_max_res_stat(struct d_ocp_qp_res* res, double* value);
//
void d_ocp_qp_res_get_max_res_eq(struct d_ocp_qp_res* res, double* value);
//
void d_ocp_qp_res_get_max_res_ineq(struct d_ocp_qp_res* res, double* value);
//
void d_ocp_qp_res_get_max_res_comp(struct d_ocp_qp_res* res, double* value);


#ifdef __cplusplus
}  // #extern "C"
#endif


#endif  // HPIPM_D_OCP_qp_res_H_
