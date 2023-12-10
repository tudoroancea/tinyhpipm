#ifndef HPIPM_S_DENSE_QCQP_DIM_H_
#define HPIPM_S_DENSE_QCQP_DIM_H_

#include "hpipm/common.h"

#ifdef __cplusplus
extern "C" {
#endif


struct s_dense_qcqp_dim {
    struct s_dense_qp_dim* qp_dim;  // dim of qp approximation
    int nv;  // number of variables
    int ne;  // number of equality constraints
    int nb;  // number of box constraints
    int ng;  // number of general constraints
    int nq;  // number of quadratic constraints
    int nsb;  // number of softened box constraints
    int nsg;  // number of softened general constraints
    int nsq;  // number of softened quadratic constraints
    int ns;  // number of softened constraints (nsb+nsg+nsq) TODO number of slacks
    hpipm_size_t memsize;
};


//
hpipm_size_t s_dense_qcqp_dim_strsize();
//
hpipm_size_t s_dense_qcqp_dim_memsize();
//
void s_dense_qcqp_dim_create(struct s_dense_qcqp_dim* dim, void* memory);
//
void s_dense_qcqp_dim_set(char* fiels_name, int value, struct s_dense_qcqp_dim* dim);
//
void s_dense_qcqp_dim_set_nv(int value, struct s_dense_qcqp_dim* dim);
//
void s_dense_qcqp_dim_set_ne(int value, struct s_dense_qcqp_dim* dim);
//
void s_dense_qcqp_dim_set_nb(int value, struct s_dense_qcqp_dim* dim);
//
void s_dense_qcqp_dim_set_ng(int value, struct s_dense_qcqp_dim* dim);
//
void s_dense_qcqp_dim_set_nq(int value, struct s_dense_qcqp_dim* dim);
//
void s_dense_qcqp_dim_set_nsb(int value, struct s_dense_qcqp_dim* dim);
//
void s_dense_qcqp_dim_set_nsg(int value, struct s_dense_qcqp_dim* dim);
//
void s_dense_qcqp_dim_set_nsq(int value, struct s_dense_qcqp_dim* dim);
//
void s_dense_qcqp_dim_set_ns(int value, struct s_dense_qcqp_dim* dim);

// getters
//
void s_dense_qcqp_dim_get_nv(struct s_dense_qcqp_dim* dim, int* value);
//
void s_dense_qcqp_dim_get_ne(struct s_dense_qcqp_dim* dim, int* value);
//
void s_dense_qcqp_dim_get_nb(struct s_dense_qcqp_dim* dim, int* value);
//
void s_dense_qcqp_dim_get_ng(struct s_dense_qcqp_dim* dim, int* value);
//
void s_dense_qcqp_dim_get_nq(struct s_dense_qcqp_dim* dim, int* value);


#ifdef __cplusplus
}  // #extern "C"
#endif


#endif  // HPIPM_S_DENSE_QCQP_DIM_H_
