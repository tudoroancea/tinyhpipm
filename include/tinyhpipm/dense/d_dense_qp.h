#ifndef HPIPM_D_d_dense_qp_H_
#define HPIPM_D_d_dense_qp_H_

#include "tinyhpipm/blas.h"
#include "tinyhpipm/common.h"
#include "tinyhpipm/dense/d_dense_qp_dim.h"


#ifdef __cplusplus
extern "C" {
#endif


struct d_dense_qp {
    struct d_dense_qp_dim* dim;
    struct mat* Hv;  // hessian of cost & vector work space
    struct mat* A;  // equality constraint matrix
    struct mat* Ct;  // inequality constraints matrix
    struct vec* gz;  // gradient of cost & gradient of slacks
    struct vec* b;  // equality constraint vector
    struct vec* d;  // inequality constraints vector
    struct vec* d_mask;  // inequality constraints mask vector
    struct vec* m;  // rhs of complementarity condition
    struct vec* Z;  // (diagonal) hessian of slacks
    int* idxb;  // indices of box constrained variables within [u; x]
    int* idxs_rev;  // index of soft constraints (reverse storage)
    hpipm_size_t memsize;  // memory size in bytes
};


//
hpipm_size_t d_dense_qp_strsize();
//
hpipm_size_t d_dense_qp_memsize(struct d_dense_qp_dim* dim);
//
void d_dense_qp_create(struct d_dense_qp_dim* dim, struct d_dense_qp* qp, void* memory);

// setters - colmaj
//
void d_dense_qp_set_all(double* H, double* g, double* A, double* b, int* idxb, double* d_lb, double* d_ub, double* C, double* d_lg, double* d_ug, double* Zl, double* Zu, double* zl, double* zu, int* idxs, double* d_ls, double* d_us, struct d_dense_qp* qp);
//
void d_dense_qp_get_all(struct d_dense_qp* qp, double* H, double* g, double* A, double* b, int* idxb, double* d_lb, double* d_ub, double* C, double* d_lg, double* d_ug, double* Zl, double* Zu, double* zl, double* zu, int* idxs, double* d_ls, double* d_us);
//
void d_dense_qp_set(char* field, void* value, struct d_dense_qp* qp);
//
void d_dense_qp_set_H(double* H, struct d_dense_qp* qp);
//
void d_dense_qp_set_g(double* g, struct d_dense_qp* qp);
//
void d_dense_qp_set_A(double* A, struct d_dense_qp* qp);
//
void d_dense_qp_set_b(double* b, struct d_dense_qp* qp);
//
void d_dense_qp_set_idxb(int* idxb, struct d_dense_qp* qp);
//
void d_dense_qp_set_Jb(double* Jb, struct d_dense_qp* qp);
//
void d_dense_qp_set_lb(double* lb, struct d_dense_qp* qp);
//
void d_dense_qp_set_lb_mask(double* lb, struct d_dense_qp* qp);
//
void d_dense_qp_set_ub(double* ub, struct d_dense_qp* qp);
//
void d_dense_qp_set_ub_mask(double* ub, struct d_dense_qp* qp);
//
void d_dense_qp_set_C(double* C, struct d_dense_qp* qp);
//
void d_dense_qp_set_lg(double* lg, struct d_dense_qp* qp);
//
void d_dense_qp_set_lg_mask(double* lg, struct d_dense_qp* qp);
//
void d_dense_qp_set_ug(double* ug, struct d_dense_qp* qp);
//
void d_dense_qp_set_ug_mask(double* ug, struct d_dense_qp* qp);
//
void d_dense_qp_set_idxs(int* idxs, struct d_dense_qp* qp);
//
void d_dense_qp_set_idxs_rev(int* idxs_rev, struct d_dense_qp* qp);
//
void d_dense_qp_set_Jsb(double* Jsb, struct d_dense_qp* qp);
//
void d_dense_qp_set_Jsg(double* Jsg, struct d_dense_qp* qp);
//
void d_dense_qp_set_Zl(double* Zl, struct d_dense_qp* qp);
//
void d_dense_qp_set_Zu(double* Zu, struct d_dense_qp* qp);
//
void d_dense_qp_set_zl(double* zl, struct d_dense_qp* qp);
//
void d_dense_qp_set_zu(double* zu, struct d_dense_qp* qp);
//
void d_dense_qp_set_lls(double* ls, struct d_dense_qp* qp);
//
void d_dense_qp_set_lls_mask(double* ls, struct d_dense_qp* qp);
//
void d_dense_qp_set_lus(double* us, struct d_dense_qp* qp);
//
void d_dense_qp_set_lus_mask(double* us, struct d_dense_qp* qp);

// getters - colmaj
//
void d_dense_qp_get_H(struct d_dense_qp* qp, double* H);
//
void d_dense_qp_get_g(struct d_dense_qp* qp, double* g);
//
void d_dense_qp_get_A(struct d_dense_qp* qp, double* A);
//
void d_dense_qp_get_b(struct d_dense_qp* qp, double* b);
//
void d_dense_qp_get_idxb(struct d_dense_qp* qp, int* idxb);
//
void d_dense_qp_get_lb(struct d_dense_qp* qp, double* lb);
//
void d_dense_qp_get_lb_mask(struct d_dense_qp* qp, double* lb);
//
void d_dense_qp_get_ub(struct d_dense_qp* qp, double* ub);
//
void d_dense_qp_get_ub_mask(struct d_dense_qp* qp, double* ub);
//
void d_dense_qp_get_C(struct d_dense_qp* qp, double* C);
//
void d_dense_qp_get_lg(struct d_dense_qp* qp, double* lg);
//
void d_dense_qp_get_lg_mask(struct d_dense_qp* qp, double* lg);
//
void d_dense_qp_get_ug(struct d_dense_qp* qp, double* ug);
//
void d_dense_qp_get_ug_mask(struct d_dense_qp* qp, double* ug);
//
void d_dense_qp_get_idxs(struct d_dense_qp* qp, int* idxs);
//
void d_dense_qp_get_idxs_rev(struct d_dense_qp* qp, int* idxs_rev);
//
void d_dense_qp_get_Zl(struct d_dense_qp* qp, double* Zl);
//
void d_dense_qp_get_Zu(struct d_dense_qp* qp, double* Zu);
//
void d_dense_qp_get_zl(struct d_dense_qp* qp, double* zl);
//
void d_dense_qp_get_zu(struct d_dense_qp* qp, double* zu);
//
void d_dense_qp_get_ls(struct d_dense_qp* qp, double* ls);
//
void d_dense_qp_get_ls_mask(struct d_dense_qp* qp, double* ls);
//
void d_dense_qp_get_us(struct d_dense_qp* qp, double* us);
//
void d_dense_qp_get_us_mask(struct d_dense_qp* qp, double* us);

// setters - rowmaj
//
void d_dense_qp_set_all_rowmaj(double* H, double* g, double* A, double* b, int* idxb, double* d_lb, double* d_ub, double* C, double* d_lg, double* d_ug, double* Zl, double* Zu, double* zl, double* zu, int* idxs, double* d_ls, double* d_us, struct d_dense_qp* qp);

// getters - rowmaj
//
void d_dense_qp_get_all_rowmaj(struct d_dense_qp* qp, double* H, double* g, double* A, double* b, int* idxb, double* d_lb, double* d_ub, double* C, double* d_lg, double* d_ug, double* Zl, double* Zu, double* zl, double* zu, int* idxs, double* d_ls, double* d_us);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif  // HPIPM_D_d_dense_qp_H_
