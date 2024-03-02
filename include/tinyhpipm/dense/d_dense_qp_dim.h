#ifndef HPIPM_D_d_dense_qp_dim_H_
#define HPIPM_D_d_dense_qp_dim_H_

#include "tinyhpipm/common.h"

#ifdef __cplusplus
extern "C" {
#endif


struct d_dense_qp_dim {
    int nv;  // number of variables
    int ne;  // number of equality constraints
    int nb;  // number of box constraints
    int ng;  // number of general constraints
    int nsb;  // number of softened box constraints
    int nsg;  // number of softened general constraints
    int ns;  // number of softened constraints (nsb+nsg)
    hpipm_size_t memsize;
};


//
hpipm_size_t d_dense_qp_dim_strsize();
//
hpipm_size_t d_dense_qp_dim_memsize();
//
void d_dense_qp_dim_create(struct d_dense_qp_dim* qp_dim, void* memory);
//
void d_dense_qp_dim_set_all(int nv, int ne, int nb, int ng, int nsb, int nsg, struct d_dense_qp_dim* dim);
//
void d_dense_qp_dim_set(char* field_name, int value, struct d_dense_qp_dim* dim);
//
void d_dense_qp_dim_set_nv(int value, struct d_dense_qp_dim* dim);
//
void d_dense_qp_dim_set_ne(int value, struct d_dense_qp_dim* dim);
//
void d_dense_qp_dim_set_nb(int value, struct d_dense_qp_dim* dim);
//
void d_dense_qp_dim_set_ng(int value, struct d_dense_qp_dim* dim);
//
void d_dense_qp_dim_set_nsb(int value, struct d_dense_qp_dim* dim);
//
void d_dense_qp_dim_set_nsg(int value, struct d_dense_qp_dim* dim);
//
void d_dense_qp_dim_set_ns(int value, struct d_dense_qp_dim* dim);

// getters
//
void d_dense_qp_dim_get_nv(struct d_dense_qp_dim* dim, int* value);
//
void d_dense_qp_dim_get_ne(struct d_dense_qp_dim* dim, int* value);
//
void d_dense_qp_dim_get_nb(struct d_dense_qp_dim* dim, int* value);
//
void d_dense_qp_dim_get_ng(struct d_dense_qp_dim* dim, int* value);


#ifdef __cplusplus
}  // #extern "C"
#endif


#endif  // HPIPM_D_d_dense_qp_dim_H_
