#ifndef HPIPM_D_OCP_QCQP_SOL_H_
#define HPIPM_D_OCP_QCQP_SOL_H_

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/ocp/d_ocp_qcqp_dim.h"


#ifdef __cplusplus
extern "C" {
#endif


struct d_ocp_qcqp_sol {
    struct d_ocp_qcqp_dim* dim;
    struct blasfeo_dvec* ux;
    struct blasfeo_dvec* pi;
    struct blasfeo_dvec* lam;
    struct blasfeo_dvec* t;
    hpipm_size_t memsize;  // memory size in bytes
};


//
hpipm_size_t d_ocp_qcqp_sol_strsize();
//
hpipm_size_t d_ocp_qcqp_sol_memsize(struct d_ocp_qcqp_dim* dim);
//
void d_ocp_qcqp_sol_create(struct d_ocp_qcqp_dim* dim, struct d_ocp_qcqp_sol* qp_sol, void* memory);
//
void d_ocp_qcqp_sol_copy_all(struct d_ocp_qcqp_sol* qp_sol_orig, struct d_ocp_qcqp_sol* qp_sol_dest);
//
void d_ocp_qcqp_sol_get(char* field, int stage, struct d_ocp_qcqp_sol* qp_sol, double* vec);
//
void d_ocp_qcqp_sol_get_u(int stage, struct d_ocp_qcqp_sol* qp_sol, double* vec);
//
void d_ocp_qcqp_sol_get_x(int stage, struct d_ocp_qcqp_sol* qp_sol, double* vec);
//
void d_ocp_qcqp_sol_get_sl(int stage, struct d_ocp_qcqp_sol* qp_sol, double* vec);
//
void d_ocp_qcqp_sol_get_su(int stage, struct d_ocp_qcqp_sol* qp_sol, double* vec);
//
void d_ocp_qcqp_sol_get_pi(int stage, struct d_ocp_qcqp_sol* qp_sol, double* vec);
//
void d_ocp_qcqp_sol_get_lam_lb(int stage, struct d_ocp_qcqp_sol* qp_sol, double* vec);
//
void d_ocp_qcqp_sol_get_lam_ub(int stage, struct d_ocp_qcqp_sol* qp_sol, double* vec);
//
void d_ocp_qcqp_sol_get_lam_lg(int stage, struct d_ocp_qcqp_sol* qp_sol, double* vec);
//
void d_ocp_qcqp_sol_get_lam_ug(int stage, struct d_ocp_qcqp_sol* qp_sol, double* vec);
//
void d_ocp_qcqp_sol_set(char* field, int stage, double* vec, struct d_ocp_qcqp_sol* qp_sol);
//
void d_ocp_qcqp_sol_set_u(int stage, double* vec, struct d_ocp_qcqp_sol* qp_sol);
//
void d_ocp_qcqp_sol_set_x(int stage, double* vec, struct d_ocp_qcqp_sol* qp_sol);
//
void d_ocp_qcqp_sol_set_sl(int stage, double* vec, struct d_ocp_qcqp_sol* qp_sol);
//
void d_ocp_qcqp_sol_set_su(int stage, double* vec, struct d_ocp_qcqp_sol* qp_sol);


#ifdef __cplusplus
}  // #extern "C"
#endif


#endif  // HPIPM_D_OCP_QCQP_SOL_H_
