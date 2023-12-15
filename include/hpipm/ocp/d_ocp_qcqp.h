#ifndef HPIPM_D_OCP_QCQP_H_
#define HPIPM_D_OCP_QCQP_H_

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/ocp/d_ocp_qcqp_dim.h"


#ifdef __cplusplus
extern "C" {
#endif


struct d_ocp_qcqp {
    struct d_ocp_qcqp_dim* dim;
    struct blasfeo_dmat* BAbt;  // dynamics matrix & vector work space
    struct blasfeo_dmat* RSQrq;  // hessian of cost & vector work space
    struct blasfeo_dmat* DCt;  // inequality constraints matrix
    struct blasfeo_dmat** Hq;  // hessians of quadratic constraints
    struct blasfeo_dvec* b;  // dynamics vector
    struct blasfeo_dvec* rqz;  // gradient of cost & gradient of slacks
    struct blasfeo_dvec* d;  // inequality constraints vector
    struct blasfeo_dvec* d_mask;  // inequality constraints mask vector
    struct blasfeo_dvec* m;  // rhs of complementarity condition
    struct blasfeo_dvec* Z;  // (diagonal) hessian of slacks
    int** idxb;  // indices of box constrained variables within [u; x]
    int** idxs_rev;  // index of soft constraints (reverse storage)
    int** idxe;  // indices of constraints within [bu, bx, g, q] that are equalities, subset of [0, ..., nbu+nbx+ng+nq-1]
    int** Hq_nzero;  // for each int, the last 3 bits ...abc, {a,b,c}=0 => {R,S,Q}=0
    hpipm_size_t memsize;  // memory size in bytes
};


//
hpipm_size_t d_ocp_qcqp_strsize();
//
hpipm_size_t d_ocp_qcqp_memsize(struct d_ocp_qcqp_dim* dim);
//
void d_ocp_qcqp_create(struct d_ocp_qcqp_dim* dim, struct d_ocp_qcqp* qp, void* memory);
//
void d_ocp_qcqp_copy_all(struct d_ocp_qcqp* qp_orig, struct d_ocp_qcqp* qp_dest);

// setters
//
void d_ocp_qcqp_set_all_zero(struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_rhs_zero(struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set(char* field_name, int stage, void* value, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_el(char* field_name, int stage, int index, void* value, struct d_ocp_qcqp* qp);
// inter-stage equality constraints (dynamics)
//
void d_ocp_qcqp_set_A(int stage, double* mat, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_B(int stage, double* mat, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_b(int stage, double* vec, struct d_ocp_qcqp* qp);
// costs
//
void d_ocp_qcqp_set_Q(int stage, double* mat, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_S(int stage, double* mat, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_R(int stage, double* mat, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_q(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_r(int stage, double* vec, struct d_ocp_qcqp* qp);
// box constraints
//
void d_ocp_qcqp_set_lb(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_lb_mask(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_ub(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_ub_mask(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_lbx(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_lbx_mask(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_el_lbx(int stage, int index, double* elem, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_ubx(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_ubx_mask(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_el_ubx(int stage, int index, double* elem, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_lbu(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_lbu_mask(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_ubu(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_ubu_mask(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_idxb(int stage, int* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_idxbu(int stage, int* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_Jbu(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_idxbx(int stage, int* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_Jbx(int stage, double* vec, struct d_ocp_qcqp* qp);
// general linear inequality constraints
//
void d_ocp_qcqp_set_C(int stage, double* mat, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_D(int stage, double* mat, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_lg(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_lg_mask(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_ug(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_ug_mask(int stage, double* vec, struct d_ocp_qcqp* qp);
// quadratic inequality constraints
//
void d_ocp_qcqp_set_Qq(int stage, double* mat, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_Sq(int stage, double* mat, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_Rq(int stage, double* mat, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_qq(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_rq(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_uq(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_uq_mask(int stage, double* vec, struct d_ocp_qcqp* qp);
// softed constraints
//
void d_ocp_qcqp_set_Zl(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_Zu(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_zl(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_zu(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_lls(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_lls_mask(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_lus(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_lus_mask(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_idxs(int stage, int* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_idxs_rev(int stage, int* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_idxsbu(int stage, int* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_Jsbu(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_idxsbx(int stage, int* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_Jsbx(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_Jsg(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_Jsq(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_idxe(int stage, int* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_idxbxe(int stage, int* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_idxbue(int stage, int* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_idxge(int stage, int* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_idxqe(int stage, int* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_Jbxe(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_Jbue(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_Jge(int stage, double* vec, struct d_ocp_qcqp* qp);
//
void d_ocp_qcqp_set_Jqe(int stage, double* vec, struct d_ocp_qcqp* qp);

// getters
//
void d_ocp_qcqp_get(char* field, int stage, struct d_ocp_qcqp* qp, void* value);
//
void d_ocp_qcqp_get_A(int stage, struct d_ocp_qcqp* qp, double* mat);
//
void d_ocp_qcqp_get_B(int stage, struct d_ocp_qcqp* qp, double* mat);
//
void d_ocp_qcqp_get_b(int stage, struct d_ocp_qcqp* qp, double* vec);
//
void d_ocp_qcqp_get_Q(int stage, struct d_ocp_qcqp* qp, double* mat);
//
void d_ocp_qcqp_get_S(int stage, struct d_ocp_qcqp* qp, double* mat);
//
void d_ocp_qcqp_get_R(int stage, struct d_ocp_qcqp* qp, double* mat);
//
void d_ocp_qcqp_get_q(int stage, struct d_ocp_qcqp* qp, double* vec);
//
void d_ocp_qcqp_get_r(int stage, struct d_ocp_qcqp* qp, double* vec);
//
void d_ocp_qcqp_get_ub(int stage, struct d_ocp_qcqp* qp, double* vec);
//
void d_ocp_qcqp_get_ub_mask(int stage, struct d_ocp_qcqp* qp, double* vec);
//
void d_ocp_qcqp_get_lb(int stage, struct d_ocp_qcqp* qp, double* vec);
//
void d_ocp_qcqp_get_lb_mask(int stage, struct d_ocp_qcqp* qp, double* vec);
//
void d_ocp_qcqp_get_lbx(int stage, struct d_ocp_qcqp* qp, double* vec);
//
void d_ocp_qcqp_get_lbx_mask(int stage, struct d_ocp_qcqp* qp, double* vec);
//
void d_ocp_qcqp_get_ubx(int stage, struct d_ocp_qcqp* qp, double* vec);
//
void d_ocp_qcqp_get_ubx_mask(int stage, struct d_ocp_qcqp* qp, double* vec);
//
void d_ocp_qcqp_get_lbu(int stage, struct d_ocp_qcqp* qp, double* vec);
//
void d_ocp_qcqp_get_lbu_mask(int stage, struct d_ocp_qcqp* qp, double* vec);
//
void d_ocp_qcqp_get_ubu(int stage, struct d_ocp_qcqp* qp, double* vec);
//
void d_ocp_qcqp_get_ubu_mask(int stage, struct d_ocp_qcqp* qp, double* vec);
//
void d_ocp_qcqp_get_idxb(int stage, struct d_ocp_qcqp* qp, int* vec);
//
// void d_ocp_qcqp_get_idxbx(int stage, struct d_ocp_qcqp *qp, int *vec);
//
// void d_ocp_qcqp_get_Jbx(int stage, struct d_ocp_qcqp *qp, double *vec);
//
// void d_ocp_qcqp_get_idxbu(int stage, struct d_ocp_qcqp *qp, int *vec);
//
// void d_ocp_qcqp_get_Jbu(int stage, struct d_ocp_qcqp *qp, double *vec);
//
void d_ocp_qcqp_get_C(int stage, struct d_ocp_qcqp* qp, double* mat);
//
void d_ocp_qcqp_get_D(int stage, struct d_ocp_qcqp* qp, double* mat);
//
void d_ocp_qcqp_get_lg(int stage, struct d_ocp_qcqp* qp, double* vec);
//
void d_ocp_qcqp_get_lg_mask(int stage, struct d_ocp_qcqp* qp, double* vec);
//
void d_ocp_qcqp_get_ug(int stage, struct d_ocp_qcqp* qp, double* vec);
//
void d_ocp_qcqp_get_ug_mask(int stage, struct d_ocp_qcqp* qp, double* vec);
//
void d_ocp_qcqp_get_Zl(int stage, struct d_ocp_qcqp* qp, double* vec);
//
void d_ocp_qcqp_get_Zu(int stage, struct d_ocp_qcqp* qp, double* vec);
//
void d_ocp_qcqp_get_zl(int stage, struct d_ocp_qcqp* qp, double* vec);
//
void d_ocp_qcqp_get_zu(int stage, struct d_ocp_qcqp* qp, double* vec);
//
void d_ocp_qcqp_get_lls(int stage, struct d_ocp_qcqp* qp, double* vec);
//
void d_ocp_qcqp_get_lls_mask(int stage, struct d_ocp_qcqp* qp, double* vec);
//
void d_ocp_qcqp_get_lus(int stage, struct d_ocp_qcqp* qp, double* vec);
//
void d_ocp_qcqp_get_lus_mask(int stage, struct d_ocp_qcqp* qp, double* vec);
// XXX only valid if there is one slack per softed constraint !!!
void d_ocp_qcqp_get_idxs(int stage, struct d_ocp_qcqp* qp, int* vec);
//
void d_ocp_qcqp_get_idxs_rev(int stage, struct d_ocp_qcqp* qp, int* vec);
//
// void d_ocp_qcqp_get_Jsbu(int stage, struct d_ocp_qcqp *qp, float *vec);
//
// void d_ocp_qcqp_get_Jsbx(int stage, struct d_ocp_qcqp *qp, float *vec);
//
// void d_ocp_qcqp_get_Jsg(int stage, struct d_ocp_qcqp *qp, float *vec);


#ifdef __cplusplus
}  // #extern "C"
#endif


#endif  // HPIPM_D_OCP_QCQP_H_
