#ifndef HPIPM_S_DENSE_QCQP_H_
#define HPIPM_S_DENSE_QCQP_H_


#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_target.h"

#include "hpipm/dense/s_dense_qcqp_dim.h"


#ifdef __cplusplus
extern "C" {
#endif


struct s_dense_qcqp {
    struct s_dense_qcqp_dim* dim;
    struct blasfeo_smat* Hv;  // hessian of cost & vector work space
    struct blasfeo_smat* A;  // equality constraint matrix
    struct blasfeo_smat* Ct;  // inequality constraints matrix
    struct blasfeo_smat* Hq;  // hessians of quadratic constraints
    struct blasfeo_svec* gz;  // gradient of cost & gradient of slacks
    struct blasfeo_svec* b;  // equality constraint vector
    struct blasfeo_svec* d;  // inequality constraints vector
    struct blasfeo_svec* d_mask;  // inequality constraints mask vector
    struct blasfeo_svec* m;  // rhs of complementarity condition
    struct blasfeo_svec* Z;  // (diagonal) hessian of slacks
    int* idxb;  // indices of box constrained variables within [u; x]
    int* idxs_rev;  // index of soft constraints (reverse storage)
    int* Hq_nzero;  // for each int, the last 3 bits ...abc, {a,b,c}=0 => {R,S,Q}=0
    hpipm_size_t memsize;  // memory size in bytes
};


//
hpipm_size_t s_dense_qcqp_strsize();
//
hpipm_size_t s_dense_qcqp_memsize(struct s_dense_qcqp_dim* dim);
//
void s_dense_qcqp_create(struct s_dense_qcqp_dim* dim, struct s_dense_qcqp* qp, void* memory);

//
void s_dense_qcqp_set(char* field, void* value, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_H(float* H, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_g(float* g, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_A(float* A, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_b(float* b, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_idxb(int* idxb, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_Jb(float* Jb, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_lb(float* lb, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_lb_mask(float* lb, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_ub(float* ub, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_ub_mask(float* ub, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_C(float* C, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_lg(float* lg, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_lg_mask(float* lg, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_ug(float* ug, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_ug_mask(float* ug, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_Hq(float* Hq, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_gq(float* gq, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_uq(float* uq, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_uq_mask(float* uq, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_idxs(int* idxs, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_idxs_rev(int* idxs_rev, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_Jsb(float* Jsb, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_Jsg(float* Jsg, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_Jsq(float* Jsq, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_Zl(float* Zl, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_Zu(float* Zu, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_zl(float* zl, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_zu(float* zu, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_ls(float* ls, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_ls_mask(float* ls, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_us(float* us, struct s_dense_qcqp* qp);
//
void s_dense_qcqp_set_us_mask(float* us, struct s_dense_qcqp* qp);

// getters (COLMAJ)

void s_dense_qcqp_get_H(struct s_dense_qcqp* qp, float* H);
//
void s_dense_qcqp_get_g(struct s_dense_qcqp* qp, float* g);
//
void s_dense_qcqp_get_A(struct s_dense_qcqp* qp, float* A);
//
void s_dense_qcqp_get_b(struct s_dense_qcqp* qp, float* b);
//
void s_dense_qcqp_get_idxb(struct s_dense_qcqp* qp, int* idxb);
//
void s_dense_qcqp_get_lb(struct s_dense_qcqp* qp, float* lb);
//
void s_dense_qcqp_get_lb_mask(struct s_dense_qcqp* qp, float* lb);
//
void s_dense_qcqp_get_ub(struct s_dense_qcqp* qp, float* ub);
//
void s_dense_qcqp_get_ub_mask(struct s_dense_qcqp* qp, float* ub);
//
void s_dense_qcqp_get_C(struct s_dense_qcqp* qp, float* C);
//
void s_dense_qcqp_get_lg(struct s_dense_qcqp* qp, float* lg);
//
void s_dense_qcqp_get_lg_mask(struct s_dense_qcqp* qp, float* lg);
//
void s_dense_qcqp_get_ug(struct s_dense_qcqp* qp, float* ug);
//
void s_dense_qcqp_get_ug_mask(struct s_dense_qcqp* qp, float* ug);
//
void s_dense_qcqp_get_idxs(struct s_dense_qcqp* qp, int* idxs);
//
void s_dense_qcqp_get_idxs_rev(struct s_dense_qcqp* qp, int* idxs_rev);
//
void s_dense_qcqp_get_Zl(struct s_dense_qcqp* qp, float* Zl);
//
void s_dense_qcqp_get_Zu(struct s_dense_qcqp* qp, float* Zu);
//
void s_dense_qcqp_get_zl(struct s_dense_qcqp* qp, float* zl);
//
void s_dense_qcqp_get_zu(struct s_dense_qcqp* qp, float* zu);
//
void s_dense_qcqp_get_ls(struct s_dense_qcqp* qp, float* ls);
//
void s_dense_qcqp_get_ls_mask(struct s_dense_qcqp* qp, float* ls);
//
void s_dense_qcqp_get_us(struct s_dense_qcqp* qp, float* us);
//
void s_dense_qcqp_get_us_mask(struct s_dense_qcqp* qp, float* us);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif  // HPIPM_S_DENSE_QCQP_H_
