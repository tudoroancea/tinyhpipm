#ifndef HPIPM_S_OCP_QCQP_DIM_H_
#define HPIPM_S_OCP_QCQP_DIM_H_

#include "hpipm/common.h"

#ifdef __cplusplus
extern "C" {
#endif


struct s_ocp_qcqp_dim {
    struct s_ocp_qp_dim* qp_dim;  // dim of qp approximation
    int* nx;  // number of states
    int* nu;  // number of inputs
    int* nb;  // number of box constraints
    int* nbx;  // number of (two-sided) state box constraints
    int* nbu;  // number of (two-sided) input box constraints
    int* ng;  // number of (two-sided) general constraints
    int* nq;  // number of (upper) quadratic constraints
    int* ns;  // number of soft constraints
    int* nsbx;  // number of soft state box constraints
    int* nsbu;  // number of soft input box constraints
    int* nsg;  // number of soft general constraints
    int* nsq;  // number of soft quadratic constraints
    int* nbxe;  // number of state box constraints which are equality
    int* nbue;  // number of input box constraints which are equality
    int* nge;  // number of general constraints which are equality
    int* nqe;  // number of quadratic constraints which are equality
    int N;  // horizon length
    hpipm_size_t memsize;
};


//
hpipm_size_t s_ocp_qcqp_dim_strsize();
//
hpipm_size_t s_ocp_qcqp_dim_memsize(int N);
//
void s_ocp_qcqp_dim_create(int N, struct s_ocp_qcqp_dim* qp_dim, void* memory);
//
void s_ocp_qcqp_dim_copy_all(struct s_ocp_qcqp_dim* dim_orig, struct s_ocp_qcqp_dim* dim_dest);
//
void s_ocp_qcqp_dim_set(char* field, int stage, int value, struct s_ocp_qcqp_dim* dim);
//
void s_ocp_qcqp_dim_set_nx(int stage, int value, struct s_ocp_qcqp_dim* dim);
//
void s_ocp_qcqp_dim_set_nu(int stage, int value, struct s_ocp_qcqp_dim* dim);
//
void s_ocp_qcqp_dim_set_nbx(int stage, int value, struct s_ocp_qcqp_dim* dim);
//
void s_ocp_qcqp_dim_set_nbu(int stage, int value, struct s_ocp_qcqp_dim* dim);
//
void s_ocp_qcqp_dim_set_ng(int stage, int value, struct s_ocp_qcqp_dim* dim);
//
void s_ocp_qcqp_dim_set_nq(int stage, int value, struct s_ocp_qcqp_dim* dim);
//
void s_ocp_qcqp_dim_set_ns(int stage, int value, struct s_ocp_qcqp_dim* dim);
//
void s_ocp_qcqp_dim_set_nsbx(int stage, int value, struct s_ocp_qcqp_dim* dim);
//
void s_ocp_qcqp_dim_set_nsbu(int stage, int value, struct s_ocp_qcqp_dim* dim);
//
void s_ocp_qcqp_dim_set_nsg(int stage, int value, struct s_ocp_qcqp_dim* dim);
//
void s_ocp_qcqp_dim_set_nsq(int stage, int value, struct s_ocp_qcqp_dim* dim);
//
void s_ocp_qcqp_dim_get(struct s_ocp_qcqp_dim* dim, char* field, int stage, int* value);
//
void s_ocp_qcqp_dim_set_nbxe(int stage, int value, struct s_ocp_qcqp_dim* dim);
//
void s_ocp_qcqp_dim_set_nbue(int stage, int value, struct s_ocp_qcqp_dim* dim);
//
void s_ocp_qcqp_dim_set_nge(int stage, int value, struct s_ocp_qcqp_dim* dim);
//
void s_ocp_qcqp_dim_set_nqe(int stage, int value, struct s_ocp_qcqp_dim* dim);
//
void s_ocp_qcqp_dim_get_N(struct s_ocp_qcqp_dim* dim, int* value);
//
void s_ocp_qcqp_dim_get_nx(struct s_ocp_qcqp_dim* dim, int stage, int* value);
//
void s_ocp_qcqp_dim_get_nu(struct s_ocp_qcqp_dim* dim, int stage, int* value);


#ifdef __cplusplus
}  // #extern "C"
#endif


#endif  // HPIPM_S_OCP_QCQP_DIM_H_
