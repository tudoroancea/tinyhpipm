#ifndef HPIPM_D_d_dense_qcqp_H_
#define HPIPM_D_d_dense_qcqp_H_

#include "hpipm/blas.h"
#include "hpipm/common.h"

#include "hpipm/dense/d_dense_qcqp_dim.h"

#ifdef __cplusplus
extern "C" {
#endif


struct d_dense_qcqp {
    struct d_dense_qcqp_dim* dim;
    struct mat* Hv;  // hessian of cost & vector work space
    struct mat* A;  // equality constraint matrix
    struct mat* Ct;  // inequality constraints matrix
    struct mat* Hq;  // hessians of quadratic constraints
    struct vec* gz;  // gradient of cost & gradient of slacks
    struct vec* b;  // equality constraint vector
    struct vec* d;  // inequality constraints vector
    struct vec* d_mask;  // inequality constraints mask vector
    struct vec* m;  // rhs of complementarity condition
    struct vec* Z;  // (diagonal) hessian of slacks
    int* idxb;  // indices of box constrained variables within [u; x]
    int* idxs_rev;  // index of soft constraints (reverse storage)
    int* Hq_nzero;  // for each int, the last 3 bits ...abc, {a,b,c}=0 => {R,S,Q}=0
    hpipm_size_t memsize;  // memory size in bytes
};


//
hpipm_size_t d_dense_qcqp_strsize();
//
hpipm_size_t d_dense_qcqp_memsize(struct d_dense_qcqp_dim* dim);
//
void d_dense_qcqp_create(struct d_dense_qcqp_dim* dim, struct d_dense_qcqp* qp, void* memory);

//
void d_dense_qcqp_set(char* field, void* value, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_H(double* H, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_g(double* g, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_A(double* A, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_b(double* b, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_idxb(int* idxb, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_Jb(double* Jb, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_lb(double* lb, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_lb_mask(double* lb, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_ub(double* ub, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_ub_mask(double* ub, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_C(double* C, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_lg(double* lg, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_lg_mask(double* lg, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_ug(double* ug, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_ug_mask(double* ug, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_Hq(double* Hq, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_gq(double* gq, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_uq(double* uq, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_uq_mask(double* uq, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_idxs(int* idxs, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_idxs_rev(int* idxs_rev, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_Jsb(double* Jsb, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_Jsg(double* Jsg, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_Jsq(double* Jsq, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_Zl(double* Zl, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_Zu(double* Zu, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_zl(double* zl, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_zu(double* zu, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_ls(double* ls, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_ls_mask(double* ls, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_us(double* us, struct d_dense_qcqp* qp);
//
void d_dense_qcqp_set_us_mask(double* us, struct d_dense_qcqp* qp);

// getters (COLMAJ)

void d_dense_qcqp_get_H(struct d_dense_qcqp* qp, double* H);
//
void d_dense_qcqp_get_g(struct d_dense_qcqp* qp, double* g);
//
void d_dense_qcqp_get_A(struct d_dense_qcqp* qp, double* A);
//
void d_dense_qcqp_get_b(struct d_dense_qcqp* qp, double* b);
//
void d_dense_qcqp_get_idxb(struct d_dense_qcqp* qp, int* idxb);
//
void d_dense_qcqp_get_lb(struct d_dense_qcqp* qp, double* lb);
//
void d_dense_qcqp_get_lb_mask(struct d_dense_qcqp* qp, double* lb);
//
void d_dense_qcqp_get_ub(struct d_dense_qcqp* qp, double* ub);
//
void d_dense_qcqp_get_ub_mask(struct d_dense_qcqp* qp, double* ub);
//
void d_dense_qcqp_get_C(struct d_dense_qcqp* qp, double* C);
//
void d_dense_qcqp_get_lg(struct d_dense_qcqp* qp, double* lg);
//
void d_dense_qcqp_get_lg_mask(struct d_dense_qcqp* qp, double* lg);
//
void d_dense_qcqp_get_ug(struct d_dense_qcqp* qp, double* ug);
//
void d_dense_qcqp_get_ug_mask(struct d_dense_qcqp* qp, double* ug);
//
void d_dense_qcqp_get_idxs(struct d_dense_qcqp* qp, int* idxs);
//
void d_dense_qcqp_get_idxs_rev(struct d_dense_qcqp* qp, int* idxs_rev);
//
void d_dense_qcqp_get_Zl(struct d_dense_qcqp* qp, double* Zl);
//
void d_dense_qcqp_get_Zu(struct d_dense_qcqp* qp, double* Zu);
//
void d_dense_qcqp_get_zl(struct d_dense_qcqp* qp, double* zl);
//
void d_dense_qcqp_get_zu(struct d_dense_qcqp* qp, double* zu);
//
void d_dense_qcqp_get_ls(struct d_dense_qcqp* qp, double* ls);
//
void d_dense_qcqp_get_ls_mask(struct d_dense_qcqp* qp, double* ls);
//
void d_dense_qcqp_get_us(struct d_dense_qcqp* qp, double* us);
//
void d_dense_qcqp_get_us_mask(struct d_dense_qcqp* qp, double* us);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif  // HPIPM_D_d_dense_qcqp_H_
