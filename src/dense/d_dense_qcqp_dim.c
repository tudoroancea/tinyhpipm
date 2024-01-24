#include <stdio.h>
#include <stdlib.h>

#include "hpipm/common.h"
#include "hpipm/dense/d_dense_qcqp_dim.h"
#include "hpipm/dense/d_dense_qp_dim.h"


hpipm_size_t d_dense_qcqp_dim_strsize() {
    return sizeof(struct d_dense_qcqp_dim);
}


hpipm_size_t d_dense_qcqp_dim_memsize() {
    hpipm_size_t size = 0;
    size += 1 * sizeof(struct d_dense_qp_dim);
    size += 1 * d_dense_qp_dim_memsize();
    size = (size + 63) / 64 * 64;  // make multiple of typical cache line size
    size += 1 * 64;  // align once to typical cache line size
    return size;
}


void d_dense_qcqp_dim_create(struct d_dense_qcqp_dim* dim, void* mem) {
    // zero memory (to avoid corrupted memory like e.g. NaN)
    hpipm_size_t memsize = d_dense_qcqp_dim_memsize();
    hpipm_zero_memset(memsize, mem);
    // qp_dim struct
    struct d_dense_qp_dim* dim_ptr = mem;
    dim->qp_dim = dim_ptr;
    dim_ptr += 1;
    // align to typical cache line size
    hpipm_size_t s_ptr = (hpipm_size_t) dim_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;
    // void
    char* c_ptr = (char*) s_ptr;
    d_dense_qp_dim_create(dim->qp_dim, c_ptr);
    c_ptr += dim->qp_dim->memsize;
    dim->memsize = d_dense_qcqp_dim_memsize();

#if defined(RUNTIME_CHECKS)
    if (c_ptr > ((char*) mem) + dim->memsize) {
        printf("\nerror: d_dense_qcqp_dim_create: outside memory bounds!\n\n");
        exit(1);
    }
#endif

    // initialize dims to zero by default
    dim->nv = 0;
    dim->ne = 0;
    dim->nb = 0;
    dim->ng = 0;
    dim->nq = 0;
    dim->ns = 0;
    dim->nsb = 0;
    dim->nsg = 0;
    dim->nsq = 0;
}


void d_dense_qcqp_dim_set(char* field_name, int value, struct d_dense_qcqp_dim* dim) {
    if (hpipm_strcmp(field_name, "nv")) {
        d_dense_qcqp_dim_set_nv(value, dim);
    } else if (hpipm_strcmp(field_name, "ne")) {
        d_dense_qcqp_dim_set_ne(value, dim);
    } else if (hpipm_strcmp(field_name, "nb")) {
        d_dense_qcqp_dim_set_nb(value, dim);
    } else if (hpipm_strcmp(field_name, "ng")) {
        d_dense_qcqp_dim_set_ng(value, dim);
    } else if (hpipm_strcmp(field_name, "nq")) {
        d_dense_qcqp_dim_set_nq(value, dim);
    } else if (hpipm_strcmp(field_name, "nsb")) {
        d_dense_qcqp_dim_set_nsb(value, dim);
    } else if (hpipm_strcmp(field_name, "nsg")) {
        d_dense_qcqp_dim_set_nsg(value, dim);
    } else if (hpipm_strcmp(field_name, "nsq")) {
        d_dense_qcqp_dim_set_nsq(value, dim);
    } else if (hpipm_strcmp(field_name, "ns")) {
        d_dense_qcqp_dim_set_ns(value, dim);
    } else {
        printf("error: set_d_ocp_qcqp_dim: wrong field %s\n", field_name);
        exit(1);
    }
}


void d_dense_qcqp_dim_set_nv(int value, struct d_dense_qcqp_dim* dim) {
    dim->nv = value;
    d_dense_qp_dim_set_nv(dim->nv, dim->qp_dim);
}


void d_dense_qcqp_dim_set_ne(int value, struct d_dense_qcqp_dim* dim) {
    dim->ne = value;
    d_dense_qp_dim_set_ne(dim->ne, dim->qp_dim);
}


void d_dense_qcqp_dim_set_nb(int value, struct d_dense_qcqp_dim* dim) {
    dim->nb = value;
    d_dense_qp_dim_set_nb(dim->nb, dim->qp_dim);
}


void d_dense_qcqp_dim_set_ng(int value, struct d_dense_qcqp_dim* dim) {
    dim->ng = value;
    d_dense_qp_dim_set_ng(dim->ng + dim->nq, dim->qp_dim);
}


void d_dense_qcqp_dim_set_nq(int value, struct d_dense_qcqp_dim* dim) {
    dim->nq = value;
    d_dense_qp_dim_set_ng(dim->ng + dim->nq, dim->qp_dim);
}


void d_dense_qcqp_dim_set_nsb(int value, struct d_dense_qcqp_dim* dim) {
    dim->nsb = value;
    dim->ns = dim->nsb + dim->nsg + dim->nsq;
    d_dense_qp_dim_set_nsb(dim->nsb, dim->qp_dim);
    d_dense_qp_dim_set_ns(dim->ns, dim->qp_dim);
}


void d_dense_qcqp_dim_set_nsg(int value, struct d_dense_qcqp_dim* dim) {
    dim->nsg = value;
    dim->ns = dim->nsb + dim->nsg + dim->nsq;
    d_dense_qp_dim_set_nsg(dim->nsg + dim->nsq, dim->qp_dim);
    d_dense_qp_dim_set_ns(dim->ns, dim->qp_dim);
}


void d_dense_qcqp_dim_set_nsq(int value, struct d_dense_qcqp_dim* dim) {
    dim->nsq = value;
    dim->ns = dim->nsb + dim->nsg + dim->nsq;
    d_dense_qp_dim_set_nsg(dim->nsg + dim->nsq, dim->qp_dim);
    d_dense_qp_dim_set_ns(dim->ns, dim->qp_dim);
}


void d_dense_qcqp_dim_set_ns(int value, struct d_dense_qcqp_dim* dim) {
    dim->ns = value;
    d_dense_qp_dim_set_ns(dim->ns, dim->qp_dim);
}

void d_dense_qcqp_dim_get_nv(struct d_dense_qcqp_dim* dim, int* value) {
    *value = dim->nv;
}


void d_dense_qcqp_dim_get_ne(struct d_dense_qcqp_dim* dim, int* value) {
    *value = dim->ne;
}


void d_dense_qcqp_dim_get_nb(struct d_dense_qcqp_dim* dim, int* value) {
    *value = dim->nb;
}


void d_dense_qcqp_dim_get_ng(struct d_dense_qcqp_dim* dim, int* value) {
    *value = dim->ng;
}


void d_dense_qcqp_dim_get_nq(struct d_dense_qcqp_dim* dim, int* value) {
    *value = dim->nq;
}
