#ifndef HPIPM_S_DENSE_QP_H_
#define HPIPM_S_DENSE_QP_H_

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_target.h"

#include "hpipm/dense/s_dense_qp_dim.h"


#ifdef __cplusplus
extern "C" {
#endif


struct s_dense_qp {
    struct s_dense_qp_dim* dim;
    struct blasfeo_smat* Hv;  // hessian & gradient
    struct blasfeo_smat* A;  // dynamics matrix
    struct blasfeo_smat* Ct;  // constraints matrix
    struct blasfeo_svec* gz;  // gradient & gradient of slacks
    struct blasfeo_svec* b;  // dynamics vector
    struct blasfeo_svec* d;  // constraints vector
    struct blasfeo_svec* d_mask;  // inequality constraints mask vector
    struct blasfeo_svec* m;  // rhs of complementarity condition
    struct blasfeo_svec* Z;  // (diagonal) hessian of slacks
    int* idxb;  // indices of box constrained variables within [u; x]
    int* idxs_rev;  // index of soft constraints (reverse storage)
    hpipm_size_t memsize;  // memory size in bytes
};


//
hpipm_size_t s_dense_qp_memsize(struct s_dense_qp_dim* dim);
//
void s_dense_qp_create(struct s_dense_qp_dim* dim, struct s_dense_qp* qp, void* memory);

// setters - colmaj
//
void s_dense_qp_set_all(float* H, float* g, float* A, float* b, int* idxb, float* d_lb, float* d_ub, float* C, float* d_lg, float* d_ug, float* Zl, float* Zu, float* zl, float* zu, int* idxs, float* d_ls, float* d_us, struct s_dense_qp* qp);
//
void s_dense_qp_get_all(struct s_dense_qp* qp, float* H, float* g, float* A, float* b, int* idxb, float* d_lb, float* d_ub, float* C, float* d_lg, float* d_ug, float* Zl, float* Zu, float* zl, float* zu, int* idxs, float* d_ls, float* d_us);
//
void s_dense_qp_set(char* field, void* value, struct s_dense_qp* qp);
//
void s_dense_qp_set_H(float* H, struct s_dense_qp* qp);
//
void s_dense_qp_set_g(float* g, struct s_dense_qp* qp);
//
void s_dense_qp_set_A(float* A, struct s_dense_qp* qp);
//
void s_dense_qp_set_b(float* b, struct s_dense_qp* qp);
//
void s_dense_qp_set_idxb(int* idxb, struct s_dense_qp* qp);
//
void s_dense_qp_set_Jb(float* Jb, struct s_dense_qp* qp);
//
void s_dense_qp_set_lb(float* lb, struct s_dense_qp* qp);
//
void s_dense_qp_set_lb_mask(float* lb, struct s_dense_qp* qp);
//
void s_dense_qp_set_ub(float* ub, struct s_dense_qp* qp);
//
void s_dense_qp_set_ub_mask(float* ub, struct s_dense_qp* qp);
//
void s_dense_qp_set_C(float* C, struct s_dense_qp* qp);
//
void s_dense_qp_set_lg(float* lg, struct s_dense_qp* qp);
//
void s_dense_qp_set_lg_mask(float* lg, struct s_dense_qp* qp);
//
void s_dense_qp_set_ug(float* ug, struct s_dense_qp* qp);
//
void s_dense_qp_set_ug_mask(float* ug, struct s_dense_qp* qp);
//
void s_dense_qp_set_idxs(int* idxs, struct s_dense_qp* qp);
//
void s_dense_qp_set_idxs_rev(int* idxs_rev, struct s_dense_qp* qp);
//
void s_dense_qp_set_Jsb(float* Jsb, struct s_dense_qp* qp);
//
void s_dense_qp_set_Jsg(float* Jsg, struct s_dense_qp* qp);
//
void s_dense_qp_set_Zl(float* Zl, struct s_dense_qp* qp);
//
void s_dense_qp_set_Zu(float* Zu, struct s_dense_qp* qp);
//
void s_dense_qp_set_zl(float* zl, struct s_dense_qp* qp);
//
void s_dense_qp_set_zu(float* zu, struct s_dense_qp* qp);
//
void s_dense_qp_set_lls(float* ls, struct s_dense_qp* qp);
//
void s_dense_qp_set_lls_mask(float* ls, struct s_dense_qp* qp);
//
void s_dense_qp_set_lus(float* us, struct s_dense_qp* qp);
//
void s_dense_qp_set_lus_mask(float* us, struct s_dense_qp* qp);

// getters - colmaj
//
void s_dense_qp_get_H(struct s_dense_qp* qp, float* H);
//
void s_dense_qp_get_g(struct s_dense_qp* qp, float* g);
//
void s_dense_qp_get_A(struct s_dense_qp* qp, float* A);
//
void s_dense_qp_get_b(struct s_dense_qp* qp, float* b);
//
void s_dense_qp_get_idxb(struct s_dense_qp* qp, int* idxb);
//
void s_dense_qp_get_lb(struct s_dense_qp* qp, float* lb);
//
void s_dense_qp_get_lb_mask(struct s_dense_qp* qp, float* lb);
//
void s_dense_qp_get_ub(struct s_dense_qp* qp, float* ub);
//
void s_dense_qp_get_ub_mask(struct s_dense_qp* qp, float* ub);
//
void s_dense_qp_get_C(struct s_dense_qp* qp, float* C);
//
void s_dense_qp_get_lg(struct s_dense_qp* qp, float* lg);
//
void s_dense_qp_get_lg_mask(struct s_dense_qp* qp, float* lg);
//
void s_dense_qp_get_ug(struct s_dense_qp* qp, float* ug);
//
void s_dense_qp_get_ug_mask(struct s_dense_qp* qp, float* ug);
//
void s_dense_qp_get_idxs(struct s_dense_qp* qp, int* idxs);
//
void s_dense_qp_get_idxs_rev(struct s_dense_qp* qp, int* idxs_rev);
//
void s_dense_qp_get_Zl(struct s_dense_qp* qp, float* Zl);
//
void s_dense_qp_get_Zu(struct s_dense_qp* qp, float* Zu);
//
void s_dense_qp_get_zl(struct s_dense_qp* qp, float* zl);
//
void s_dense_qp_get_zu(struct s_dense_qp* qp, float* zu);
//
void s_dense_qp_get_ls(struct s_dense_qp* qp, float* ls);
//
void s_dense_qp_get_ls_mask(struct s_dense_qp* qp, float* ls);
//
void s_dense_qp_get_us(struct s_dense_qp* qp, float* us);
//
void s_dense_qp_get_us_mask(struct s_dense_qp* qp, float* us);

// setters - rowmaj
//
void s_dense_qp_set_all_rowmaj(float* H, float* g, float* A, float* b, int* idxb, float* d_lb, float* d_ub, float* C, float* d_lg, float* d_ug, float* Zl, float* Zu, float* zl, float* zu, int* idxs, float* d_ls, float* d_us, struct s_dense_qp* qp);

// getters - rowmaj
//
void s_dense_qp_get_all_rowmaj(struct s_dense_qp* qp, float* H, float* g, float* A, float* b, int* idxb, float* d_lb, float* d_ub, float* C, float* d_lg, float* d_ug, float* Zl, float* Zu, float* zl, float* zu, int* idxs, float* d_ls, float* d_us);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif  // HPIPM_S_DENSE_QP_H_
