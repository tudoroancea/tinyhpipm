#ifndef HPIPM_S_DENSE_QP_DIM_H_
#define HPIPM_S_DENSE_QP_DIM_H_

#include "hpipm/common.h"

#ifdef __cplusplus
extern "C" {
#endif


struct s_dense_qp_dim {
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
hpipm_size_t s_dense_qp_dim_strsize();
//
hpipm_size_t s_dense_qp_dim_memsize();
//
void s_dense_qp_dim_create(struct s_dense_qp_dim* qp_dim, void* memory);
//
void s_dense_qp_dim_set_all(int nv, int ne, int nb, int ng, int nsb, int nsg, struct s_dense_qp_dim* dim);
//
void s_dense_qp_dim_set(char* fiels_name, int value, struct s_dense_qp_dim* dim);
//
void s_dense_qp_dim_set_nv(int value, struct s_dense_qp_dim* dim);
//
void s_dense_qp_dim_set_ne(int value, struct s_dense_qp_dim* dim);
//
void s_dense_qp_dim_set_nb(int value, struct s_dense_qp_dim* dim);
//
void s_dense_qp_dim_set_ng(int value, struct s_dense_qp_dim* dim);
//
void s_dense_qp_dim_set_nsb(int value, struct s_dense_qp_dim* dim);
//
void s_dense_qp_dim_set_nsg(int value, struct s_dense_qp_dim* dim);
//
void s_dense_qp_dim_set_ns(int value, struct s_dense_qp_dim* dim);

// getters
//
void s_dense_qp_dim_get_nv(struct s_dense_qp_dim* dim, int* value);
//
void s_dense_qp_dim_get_ne(struct s_dense_qp_dim* dim, int* value);
//
void s_dense_qp_dim_get_nb(struct s_dense_qp_dim* dim, int* value);
//
void s_dense_qp_dim_get_ng(struct s_dense_qp_dim* dim, int* value);


#ifdef __cplusplus
}  // #extern "C"
#endif


#endif  // HPIPM_S_DENSE_QP_DIM_H_
