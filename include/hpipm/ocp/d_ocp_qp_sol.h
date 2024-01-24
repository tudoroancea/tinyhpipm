#ifndef HPIPM_D_OCP_qp_sol_H_
#define HPIPM_D_OCP_qp_sol_H_

#include "hpipm/blas.h"
#include "hpipm/common.h"
#include "hpipm/ocp/d_ocp_qp_dim.h"


#ifdef __cplusplus
extern "C" {
#endif


struct d_ocp_qp_sol {
    struct d_ocp_qp_dim* dim;
    struct vec* ux;
    struct vec* pi;
    struct vec* lam;
    struct vec* t;
    void* misc;
    hpipm_size_t memsize;  // memory size in bytes
};


//
hpipm_size_t d_ocp_qp_sol_strsize();
//
hpipm_size_t d_ocp_qp_sol_memsize(struct d_ocp_qp_dim* dim);
//
void d_ocp_qp_sol_create(struct d_ocp_qp_dim* dim, struct d_ocp_qp_sol* qp_sol, void* memory);
//
void d_ocp_qp_sol_copy_all(struct d_ocp_qp_sol* qp_sol_orig, struct d_ocp_qp_sol* qp_sol_dest);
//
void d_ocp_qp_sol_get_all(struct d_ocp_qp_sol* qp_sol, double** u, double** x, double** ls, double** us, double** pi, double** lam_lb, double** lam_ub, double** lam_lg, double** lam_ug, double** lam_ls, double** lam_us);
//
void d_ocp_qp_sol_get_all_rowmaj(struct d_ocp_qp_sol* qp_sol, double** u, double** x, double** ls, double** us, double** pi, double** lam_lb, double** lam_ub, double** lam_lg, double** lam_ug, double** lam_ls, double** lam_us);
//
void d_ocp_qp_sol_set_all(double** u, double** x, double** ls, double** us, double** pi, double** lam_lb, double** lam_ub, double** lam_lg, double** lam_ug, double** lam_ls, double** lam_us, struct d_ocp_qp_sol* qp_sol);
//
void d_ocp_qp_sol_get(char* field, int stage, struct d_ocp_qp_sol* qp_sol, double* vec);
//
void d_ocp_qp_sol_get_u(int stage, struct d_ocp_qp_sol* qp_sol, double* vec);
//
void d_ocp_qp_sol_get_x(int stage, struct d_ocp_qp_sol* qp_sol, double* vec);
//
void d_ocp_qp_sol_get_sl(int stage, struct d_ocp_qp_sol* qp_sol, double* vec);
//
void d_ocp_qp_sol_get_su(int stage, struct d_ocp_qp_sol* qp_sol, double* vec);
//
void d_ocp_qp_sol_get_pi(int stage, struct d_ocp_qp_sol* qp_sol, double* vec);
//
void d_ocp_qp_sol_get_lam_lb(int stage, struct d_ocp_qp_sol* qp_sol, double* vec);
//
void d_ocp_qp_sol_get_lam_lbu(int stage, struct d_ocp_qp_sol* qp_sol, double* vec);
//
void d_ocp_qp_sol_get_lam_lbx(int stage, struct d_ocp_qp_sol* qp_sol, double* vec);
//
void d_ocp_qp_sol_get_lam_ub(int stage, struct d_ocp_qp_sol* qp_sol, double* vec);
//
void d_ocp_qp_sol_get_lam_ubu(int stage, struct d_ocp_qp_sol* qp_sol, double* vec);
//
void d_ocp_qp_sol_get_lam_ubx(int stage, struct d_ocp_qp_sol* qp_sol, double* vec);
//
void d_ocp_qp_sol_get_lam_lg(int stage, struct d_ocp_qp_sol* qp_sol, double* vec);
//
void d_ocp_qp_sol_get_lam_ug(int stage, struct d_ocp_qp_sol* qp_sol, double* vec);
//
void d_ocp_qp_sol_set(char* field, int stage, double* vec, struct d_ocp_qp_sol* qp_sol);
//
void d_ocp_qp_sol_set_u(int stage, double* vec, struct d_ocp_qp_sol* qp_sol);
//
void d_ocp_qp_sol_set_x(int stage, double* vec, struct d_ocp_qp_sol* qp_sol);
//
void d_ocp_qp_sol_set_sl(int stage, double* vec, struct d_ocp_qp_sol* qp_sol);
//
void d_ocp_qp_sol_set_su(int stage, double* vec, struct d_ocp_qp_sol* qp_sol);


#ifdef __cplusplus
}  // #extern "C"
#endif


#endif  // HPIPM_D_OCP_qp_sol_H_
