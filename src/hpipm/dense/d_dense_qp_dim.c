#include <stdio.h>
#include <stdlib.h>

#include "hpipm/common.h"
#include "hpipm/dense/d_dense_qp_dim.h"


hpipm_size_t d_dense_qp_dim_strsize() {

    hpipm_size_t size = 0;

    size += sizeof(struct d_dense_qp_dim);

    return size;
}


hpipm_size_t d_dense_qp_dim_memsize() {

    hpipm_size_t size = 0;

    size = (size + 8 - 1) / 8 * 8;

    return size;
}


void d_dense_qp_dim_create(struct d_dense_qp_dim* dim, void* mem) {

    // zero memory (to avoid corrupted memory like e.g. NaN)
    //	hpipm_size_t memsize = d_dense_qp_dim_memsize();
    //	hpipm_zero_memset(memsize, mem);

    dim->memsize = d_dense_qp_dim_memsize();

    // initialize dims to zero by default

    dim->nv = 0;
    dim->ne = 0;
    dim->nb = 0;
    dim->ng = 0;
    dim->ns = 0;
    dim->nsb = 0;
    dim->nsg = 0;
}


void d_dense_qp_dim_set_all(int nv, int ne, int nb, int ng, int nsb, int nsg, struct d_dense_qp_dim* dim) {

    dim->nv = nv;
    dim->ne = ne;
    dim->nb = nb;
    dim->ng = ng;
    dim->ns = nsb + nsg;
    dim->nsb = nsb;
    dim->nsg = nsg;
}


void d_dense_qp_dim_set(char* field_name, int value, struct d_dense_qp_dim* dim) {
    if (hpipm_strcmp(field_name, "nv")) {
        d_dense_qp_dim_set_nv(value, dim);
    } else if (hpipm_strcmp(field_name, "ne")) {
        d_dense_qp_dim_set_ne(value, dim);
    } else if (hpipm_strcmp(field_name, "nb")) {
        d_dense_qp_dim_set_nb(value, dim);
    } else if (hpipm_strcmp(field_name, "ng")) {
        d_dense_qp_dim_set_ng(value, dim);
    } else if (hpipm_strcmp(field_name, "nsb")) {
        d_dense_qp_dim_set_nsb(value, dim);
    } else if (hpipm_strcmp(field_name, "nsg")) {
        d_dense_qp_dim_set_nsg(value, dim);
    } else if (hpipm_strcmp(field_name, "ns")) {
        d_dense_qp_dim_set_ns(value, dim);
    } else {
        printf("error: set_d_ocp_qp_dim: wrong field %s\n", field_name);
        exit(1);
    }
}


void d_dense_qp_dim_set_nv(int value, struct d_dense_qp_dim* dim) {
    dim->nv = value;
}


void d_dense_qp_dim_set_ne(int value, struct d_dense_qp_dim* dim) {
    dim->ne = value;
}


void d_dense_qp_dim_set_nb(int value, struct d_dense_qp_dim* dim) {
    dim->nb = value;
}


void d_dense_qp_dim_set_ng(int value, struct d_dense_qp_dim* dim) {
    dim->ng = value;
}


void d_dense_qp_dim_set_nsb(int value, struct d_dense_qp_dim* dim) {
    dim->nsb = value;
    dim->ns = dim->nsb + dim->nsg;
}


void d_dense_qp_dim_set_nsg(int value, struct d_dense_qp_dim* dim) {
    dim->nsg = value;
    dim->ns = dim->nsb + dim->nsg;
}


void d_dense_qp_dim_set_ns(int value, struct d_dense_qp_dim* dim) {
    dim->ns = value;
}


void d_dense_qp_dim_get_nv(struct d_dense_qp_dim* dim, int* value) {
    *value = dim->nv;
}


void d_dense_qp_dim_get_ne(struct d_dense_qp_dim* dim, int* value) {
    *value = dim->ne;
}


void d_dense_qp_dim_get_nb(struct d_dense_qp_dim* dim, int* value) {
    *value = dim->nb;
}


void d_dense_qp_dim_get_ng(struct d_dense_qp_dim* dim, int* value) {
    *value = dim->ng;
}
