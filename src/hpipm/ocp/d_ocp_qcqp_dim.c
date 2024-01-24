#include <stdio.h>
#include <stdlib.h>

#include "hpipm/blas.h"
#include "hpipm/common.h"
#include "hpipm/ocp/d_ocp_qcqp_dim.h"
#include "hpipm/ocp/d_ocp_qp_dim.h"


hpipm_size_t d_ocp_qcqp_dim_strsize() {

    hpipm_size_t size = 0;

    size += sizeof(struct d_ocp_qcqp_dim);

    return size;
}


hpipm_size_t d_ocp_qcqp_dim_memsize(int N) {

    hpipm_size_t size = 0;

    size += 16 * (N + 1) * sizeof(int);

    size += 1 * sizeof(struct d_ocp_qp_dim);
    size += 1 * d_ocp_qp_dim_memsize(N);

    size = (size + 63) / 64 * 64;  // make multiple of typical cache line size
    size += 1 * 64;  // align once to typical cache line size

    return size;
}


void d_ocp_qcqp_dim_create(int N, struct d_ocp_qcqp_dim* dim, void* mem) {

    // loop index
    int ii;

    // zero memory (to avoid corrupted memory like e.g. NaN)
    hpipm_size_t memsize = d_ocp_qcqp_dim_memsize(N);
    hpipm_zero_memset(memsize, mem);

    // qp_dim struct
    struct d_ocp_qp_dim* dim_ptr = mem;

    dim->qp_dim = dim_ptr;
    dim_ptr += 1;

    // align to typical cache line size
    hpipm_size_t s_ptr = (hpipm_size_t) dim_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;

    // void
    char* c_ptr = (char*) s_ptr;

    d_ocp_qp_dim_create(N, dim->qp_dim, c_ptr);
    c_ptr += dim->qp_dim->memsize;

    // nx
    dim->nx = (int*) c_ptr;
    c_ptr += (N + 1) * sizeof(int);
    // nu
    dim->nu = (int*) c_ptr;
    c_ptr += (N + 1) * sizeof(int);
    // nb
    dim->nb = (int*) c_ptr;
    c_ptr += (N + 1) * sizeof(int);
    // nbx
    dim->nbx = (int*) c_ptr;
    c_ptr += (N + 1) * sizeof(int);
    // nbu
    dim->nbu = (int*) c_ptr;
    c_ptr += (N + 1) * sizeof(int);
    // ng
    dim->ng = (int*) c_ptr;
    c_ptr += (N + 1) * sizeof(int);
    // nq
    dim->nq = (int*) c_ptr;
    c_ptr += (N + 1) * sizeof(int);
    // ns
    dim->ns = (int*) c_ptr;
    c_ptr += (N + 1) * sizeof(int);
    // nsbx
    dim->nsbx = (int*) c_ptr;
    c_ptr += (N + 1) * sizeof(int);
    // nsbu
    dim->nsbu = (int*) c_ptr;
    c_ptr += (N + 1) * sizeof(int);
    // nsg
    dim->nsg = (int*) c_ptr;
    c_ptr += (N + 1) * sizeof(int);
    // nsq
    dim->nsq = (int*) c_ptr;
    c_ptr += (N + 1) * sizeof(int);
    // nbxe
    dim->nbxe = (int*) c_ptr;
    c_ptr += (N + 1) * sizeof(int);
    // nbue
    dim->nbue = (int*) c_ptr;
    c_ptr += (N + 1) * sizeof(int);
    // nge
    dim->nge = (int*) c_ptr;
    c_ptr += (N + 1) * sizeof(int);
    // nqe
    dim->nqe = (int*) c_ptr;
    c_ptr += (N + 1) * sizeof(int);

    // N
    dim->N = N;

    // initialize dims to zero by default
    // XXX already zero-initialized while zeroing all memory

    dim->memsize = memsize;  // d_ocp_qcqp_dim_memsize(N);
}


void d_ocp_qcqp_dim_copy_all(struct d_ocp_qcqp_dim* dim_orig, struct d_ocp_qcqp_dim* dim_dest) {

#if defined(RUNTIME_CHECKS)
    if (dim_orig->N != dim_dest->N) {
        printf("\nerror: d_ocp_qcqp_dim_copy_all: dim_orig->N != dim_dest->N\n");
        exit(1);
    }
#endif

    // loop index
    int ii;

    // N
    int N = dim_orig->N;

    // copy qp dim
    for (ii = 0; ii <= N; ii++)
        d_ocp_qcqp_dim_set_nx(ii, dim_orig->nx[ii], dim_dest);
    for (ii = 0; ii <= N; ii++)
        d_ocp_qcqp_dim_set_nu(ii, dim_orig->nu[ii], dim_dest);
    for (ii = 0; ii <= N; ii++)
        d_ocp_qcqp_dim_set_nbx(ii, dim_orig->nbx[ii], dim_dest);
    for (ii = 0; ii <= N; ii++)
        d_ocp_qcqp_dim_set_nbu(ii, dim_orig->nbu[ii], dim_dest);
    for (ii = 0; ii <= N; ii++)
        d_ocp_qcqp_dim_set_nq(ii, dim_orig->nq[ii], dim_dest);
    for (ii = 0; ii <= N; ii++)
        d_ocp_qcqp_dim_set_ng(ii, dim_orig->ng[ii], dim_dest);
    for (ii = 0; ii <= N; ii++)
        d_ocp_qcqp_dim_set_nsbx(ii, dim_orig->nsbx[ii], dim_dest);
    for (ii = 0; ii <= N; ii++)
        d_ocp_qcqp_dim_set_nsbu(ii, dim_orig->nsbu[ii], dim_dest);
    for (ii = 0; ii <= N; ii++)
        d_ocp_qcqp_dim_set_nsg(ii, dim_orig->nsg[ii], dim_dest);
    for (ii = 0; ii <= N; ii++)
        d_ocp_qcqp_dim_set_nsq(ii, dim_orig->nsq[ii], dim_dest);
    for (ii = 0; ii <= N; ii++)
        d_ocp_qcqp_dim_set_nbxe(ii, dim_orig->nbxe[ii], dim_dest);
    for (ii = 0; ii <= N; ii++)
        d_ocp_qcqp_dim_set_nbue(ii, dim_orig->nbue[ii], dim_dest);
    for (ii = 0; ii <= N; ii++)
        d_ocp_qcqp_dim_set_nge(ii, dim_orig->nge[ii], dim_dest);
    for (ii = 0; ii <= N; ii++)
        d_ocp_qcqp_dim_set_nq(ii, dim_orig->nqe[ii], dim_dest);
}


void d_ocp_qcqp_dim_set(char* field_name, int stage, int value, struct d_ocp_qcqp_dim* dim) {
    if (hpipm_strcmp(field_name, "nx")) {
        d_ocp_qcqp_dim_set_nx(stage, value, dim);
    } else if (hpipm_strcmp(field_name, "nu")) {
        d_ocp_qcqp_dim_set_nu(stage, value, dim);
    } else if (hpipm_strcmp(field_name, "nbx")) {
        d_ocp_qcqp_dim_set_nbx(stage, value, dim);
    } else if (hpipm_strcmp(field_name, "nbu")) {
        d_ocp_qcqp_dim_set_nbu(stage, value, dim);
    } else if (hpipm_strcmp(field_name, "ng")) {
        d_ocp_qcqp_dim_set_ng(stage, value, dim);
    } else if (hpipm_strcmp(field_name, "nq")) {
        d_ocp_qcqp_dim_set_nq(stage, value, dim);
    } else if (hpipm_strcmp(field_name, "ns")) {
        d_ocp_qcqp_dim_set_ns(stage, value, dim);
    } else if (hpipm_strcmp(field_name, "nsbx")) {
        d_ocp_qcqp_dim_set_nsbx(stage, value, dim);
    } else if (hpipm_strcmp(field_name, "nsbu")) {
        d_ocp_qcqp_dim_set_nsbu(stage, value, dim);
    } else if (hpipm_strcmp(field_name, "nsg")) {
        d_ocp_qcqp_dim_set_nsg(stage, value, dim);
    } else if (hpipm_strcmp(field_name, "nsq")) {
        d_ocp_qcqp_dim_set_nsq(stage, value, dim);
    } else if (hpipm_strcmp(field_name, "nbxe")) {
        d_ocp_qcqp_dim_set_nbxe(stage, value, dim);
    } else if (hpipm_strcmp(field_name, "nbue")) {
        d_ocp_qcqp_dim_set_nbue(stage, value, dim);
    } else if (hpipm_strcmp(field_name, "nge")) {
        d_ocp_qcqp_dim_set_nge(stage, value, dim);
    } else if (hpipm_strcmp(field_name, "nqe")) {
        d_ocp_qcqp_dim_set_nq(stage, value, dim);
    } else {
        printf("error: d_ocp_qcqp_dim_set: wrong field %s\n", field_name);
        exit(1);
    }
}


void d_ocp_qcqp_dim_set_nx(int stage, int value, struct d_ocp_qcqp_dim* dim) {
    dim->nx[stage] = value;
    d_ocp_qp_dim_set_nx(stage, dim->nx[stage], dim->qp_dim);
}


void d_ocp_qcqp_dim_set_nu(int stage, int value, struct d_ocp_qcqp_dim* dim) {
    dim->nu[stage] = value;
    d_ocp_qp_dim_set_nu(stage, dim->nu[stage], dim->qp_dim);
}


void d_ocp_qcqp_dim_set_nbx(int stage, int value, struct d_ocp_qcqp_dim* dim) {
    dim->nbx[stage] = value;
    dim->nb[stage] = dim->nbx[stage] + dim->nbu[stage];
    d_ocp_qp_dim_set_nbx(stage, dim->nbx[stage], dim->qp_dim);
}


void d_ocp_qcqp_dim_set_nbu(int stage, int value, struct d_ocp_qcqp_dim* dim) {
    dim->nbu[stage] = value;
    dim->nb[stage] = dim->nbx[stage] + dim->nbu[stage];
    d_ocp_qp_dim_set_nbu(stage, dim->nbu[stage], dim->qp_dim);
}


void d_ocp_qcqp_dim_set_ng(int stage, int value, struct d_ocp_qcqp_dim* dim) {
    dim->ng[stage] = value;
    d_ocp_qp_dim_set_ng(stage, dim->ng[stage] + dim->nq[stage], dim->qp_dim);
}


void d_ocp_qcqp_dim_set_nq(int stage, int value, struct d_ocp_qcqp_dim* dim) {
    dim->nq[stage] = value;
    d_ocp_qp_dim_set_ng(stage, dim->ng[stage] + dim->nq[stage], dim->qp_dim);
}


void d_ocp_qcqp_dim_set_ns(int stage, int value, struct d_ocp_qcqp_dim* dim) {
    dim->ns[stage] = value;
    d_ocp_qp_dim_set_ns(stage, dim->ns[stage], dim->qp_dim);
}


void d_ocp_qcqp_dim_set_nsbx(int stage, int value, struct d_ocp_qcqp_dim* dim) {
    dim->nsbx[stage] = value;
    dim->ns[stage] = dim->nsbx[stage] + dim->nsbu[stage] + dim->nsg[stage] + dim->nsq[stage];
    d_ocp_qp_dim_set_nsbx(stage, dim->nsbx[stage], dim->qp_dim);
}


void d_ocp_qcqp_dim_set_nsbu(int stage, int value, struct d_ocp_qcqp_dim* dim) {
    dim->nsbu[stage] = value;
    dim->ns[stage] = dim->nsbx[stage] + dim->nsbu[stage] + dim->nsg[stage] + dim->nsq[stage];
    d_ocp_qp_dim_set_nsbu(stage, dim->nsbu[stage], dim->qp_dim);
}


void d_ocp_qcqp_dim_set_nsg(int stage, int value, struct d_ocp_qcqp_dim* dim) {
    dim->nsg[stage] = value;
    dim->ns[stage] = dim->nsbx[stage] + dim->nsbu[stage] + dim->nsg[stage] + dim->nsq[stage];
    d_ocp_qp_dim_set_nsg(stage, dim->nsg[stage] + dim->nsq[stage], dim->qp_dim);
}


void d_ocp_qcqp_dim_set_nsq(int stage, int value, struct d_ocp_qcqp_dim* dim) {
    dim->nsq[stage] = value;
    dim->ns[stage] = dim->nsbx[stage] + dim->nsbu[stage] + dim->nsg[stage] + dim->nsq[stage];
    d_ocp_qp_dim_set_nsg(stage, dim->nsg[stage] + dim->nsq[stage], dim->qp_dim);
}


void d_ocp_qcqp_dim_set_nbxe(int stage, int value, struct d_ocp_qcqp_dim* dim) {
    dim->nbxe[stage] = value;
    d_ocp_qp_dim_set_nbxe(stage, dim->nbxe[stage], dim->qp_dim);
}


void d_ocp_qcqp_dim_set_nbue(int stage, int value, struct d_ocp_qcqp_dim* dim) {
    dim->nbue[stage] = value;
    d_ocp_qp_dim_set_nbue(stage, dim->nbue[stage], dim->qp_dim);
}


void d_ocp_qcqp_dim_set_nge(int stage, int value, struct d_ocp_qcqp_dim* dim) {
    dim->nge[stage] = value;
    d_ocp_qp_dim_set_nge(stage, dim->nge[stage] + dim->nqe[stage], dim->qp_dim);
}


void d_ocp_qcqp_dim_set_nqe(int stage, int value, struct d_ocp_qcqp_dim* dim) {
    dim->nqe[stage] = value;
    d_ocp_qp_dim_set_nge(stage, dim->nge[stage] + dim->nqe[stage], dim->qp_dim);
}


void d_ocp_qcqp_dim_get(struct d_ocp_qcqp_dim* dim, char* field_name, int stage, int* value) {
    if (hpipm_strcmp(field_name, "nx")) {
        d_ocp_qcqp_dim_get_nx(dim, stage, value);
    } else if (hpipm_strcmp(field_name, "nu")) {
        d_ocp_qcqp_dim_get_nu(dim, stage, value);
    } else {
        printf("error: d_ocp_qcqp_dim_get: wrong field %s\n", field_name);
        exit(1);
    }
}


void d_ocp_qcqp_dim_get_N(struct d_ocp_qcqp_dim* dim, int* value) {
    *value = dim->N;
}


void d_ocp_qcqp_dim_get_nx(struct d_ocp_qcqp_dim* dim, int stage, int* value) {
    *value = dim->nx[stage];
}


void d_ocp_qcqp_dim_get_nu(struct d_ocp_qcqp_dim* dim, int stage, int* value) {
    *value = dim->nu[stage];
}
