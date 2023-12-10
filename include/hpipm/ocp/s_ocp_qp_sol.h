#ifndef HPIPM_S_OCP_QP_SOL_H_
#define HPIPM_S_OCP_QP_SOL_H_

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/ocp/s_ocp_qp_dim.h"


#ifdef __cplusplus
extern "C" {
#endif


struct s_ocp_qp_sol {
    struct s_ocp_qp_dim* dim;
    struct blasfeo_svec* ux;
    struct blasfeo_svec* pi;
    struct blasfeo_svec* lam;
    struct blasfeo_svec* t;
    void* misc;
    hpipm_size_t memsize;  // memory size in bytes
};


//
hpipm_size_t s_ocp_qp_sol_strsize();
//
hpipm_size_t s_ocp_qp_sol_memsize(struct s_ocp_qp_dim* dim);
//
void s_ocp_qp_sol_create(struct s_ocp_qp_dim* dim, struct s_ocp_qp_sol* qp_sol, void* memory);
//
void s_ocp_qp_sol_copy_all(struct s_ocp_qp_sol* qp_sol_orig, struct s_ocp_qp_sol* qp_sol_dest);
//
void s_qp_sol_get_all(struct s_ocp_qp_sol* qp_sol, float** u, float** x, float** ls, float** us, float** pi, float** lam_lb, float** lam_ub, float** lam_lg, float** lam_ug, float** lam_ls, float** lam_us);
//
void s_qp_sol_get_all_rowmaj(struct s_ocp_qp_sol* qp_sol, float** u, float** x, float** ls, float** us, float** pi, float** lam_lb, float** lam_ub, float** lam_lg, float** lam_ug, float** lam_ls, float** lam_us);
//
void s_ocp_qp_sol_set_all(float** u, float** x, float** ls, float** us, float** pi, float** lam_lb, float** lam_ub, float** lam_lg, float** lam_ug, float** lam_ls, float** lam_us, struct s_ocp_qp_sol* qp_sol);
//
void s_ocp_qp_sol_get(char* field, int stage, struct s_ocp_qp_sol* qp_sol, float* vec);
//
void s_ocp_qp_sol_get_u(int stage, struct s_ocp_qp_sol* qp_sol, float* vec);
//
void s_ocp_qp_sol_get_x(int stage, struct s_ocp_qp_sol* qp_sol, float* vec);
//
void s_ocp_qp_sol_get_sl(int stage, struct s_ocp_qp_sol* qp_sol, float* vec);
//
void s_ocp_qp_sol_get_su(int stage, struct s_ocp_qp_sol* qp_sol, float* vec);
//
void s_ocp_qp_sol_get_pi(int stage, struct s_ocp_qp_sol* qp_sol, float* vec);
//
void s_ocp_qp_sol_get_lam_lb(int stage, struct s_ocp_qp_sol* qp_sol, float* vec);
//
void s_ocp_qp_sol_get_lam_lbu(int stage, struct s_ocp_qp_sol* qp_sol, float* vec);
//
void s_ocp_qp_sol_get_lam_lbx(int stage, struct s_ocp_qp_sol* qp_sol, float* vec);
//
void s_ocp_qp_sol_get_lam_ub(int stage, struct s_ocp_qp_sol* qp_sol, float* vec);
//
void s_ocp_qp_sol_get_lam_ubu(int stage, struct s_ocp_qp_sol* qp_sol, float* vec);
//
void s_ocp_qp_sol_get_lam_ubx(int stage, struct s_ocp_qp_sol* qp_sol, float* vec);
//
void s_ocp_qp_sol_get_lam_lg(int stage, struct s_ocp_qp_sol* qp_sol, float* vec);
//
void s_ocp_qp_sol_get_lam_ug(int stage, struct s_ocp_qp_sol* qp_sol, float* vec);
//
void s_ocp_qp_sol_set(char* field, int stage, float* vec, struct s_ocp_qp_sol* qp_sol);
//
void s_ocp_qp_sol_set_u(int stage, float* vec, struct s_ocp_qp_sol* qp_sol);
//
void s_ocp_qp_sol_set_x(int stage, float* vec, struct s_ocp_qp_sol* qp_sol);
//
void s_ocp_qp_sol_set_sl(int stage, float* vec, struct s_ocp_qp_sol* qp_sol);
//
void s_ocp_qp_sol_set_su(int stage, float* vec, struct s_ocp_qp_sol* qp_sol);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif  // HPIPM_S_OCP_QP_SOL_H_
