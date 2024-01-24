#ifndef HPIPM_D_d_dense_qcqp_dim_H_
#define HPIPM_D_d_dense_qcqp_dim_H_

#include "hpipm/common.h"

#ifdef __cplusplus
extern "C" {
#endif


struct d_dense_qcqp_dim {
    struct d_dense_qp_dim* qp_dim;  // dim of qp approximation
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
hpipm_size_t d_dense_qcqp_dim_strsize();
//
hpipm_size_t d_dense_qcqp_dim_memsize();
//
void d_dense_qcqp_dim_create(struct d_dense_qcqp_dim* dim, void* memory);
//
void d_dense_qcqp_dim_set(char* field_name, int value, struct d_dense_qcqp_dim* dim);
//
void d_dense_qcqp_dim_set_nv(int value, struct d_dense_qcqp_dim* dim);
//
void d_dense_qcqp_dim_set_ne(int value, struct d_dense_qcqp_dim* dim);
//
void d_dense_qcqp_dim_set_nb(int value, struct d_dense_qcqp_dim* dim);
//
void d_dense_qcqp_dim_set_ng(int value, struct d_dense_qcqp_dim* dim);
//
void d_dense_qcqp_dim_set_nq(int value, struct d_dense_qcqp_dim* dim);
//
void d_dense_qcqp_dim_set_nsb(int value, struct d_dense_qcqp_dim* dim);
//
void d_dense_qcqp_dim_set_nsg(int value, struct d_dense_qcqp_dim* dim);
//
void d_dense_qcqp_dim_set_nsq(int value, struct d_dense_qcqp_dim* dim);
//
void d_dense_qcqp_dim_set_ns(int value, struct d_dense_qcqp_dim* dim);

// getters
//
void d_dense_qcqp_dim_get_nv(struct d_dense_qcqp_dim* dim, int* value);
//
void d_dense_qcqp_dim_get_ne(struct d_dense_qcqp_dim* dim, int* value);
//
void d_dense_qcqp_dim_get_nb(struct d_dense_qcqp_dim* dim, int* value);
//
void d_dense_qcqp_dim_get_ng(struct d_dense_qcqp_dim* dim, int* value);
//
void d_dense_qcqp_dim_get_nq(struct d_dense_qcqp_dim* dim, int* value);


#ifdef __cplusplus
}  // #extern "C"
#endif


#endif  // HPIPM_D_d_dense_qcqp_dim_H_
