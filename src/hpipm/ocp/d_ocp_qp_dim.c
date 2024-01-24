#include <stdio.h>
#include <stdlib.h>

#include "hpipm/ocp/d_ocp_qp_dim.h"


hpipm_size_t d_ocp_qp_dim_strsize() {
    return sizeof(struct d_ocp_qp_dim);
}


hpipm_size_t d_ocp_qp_dim_memsize(int N) {
    hpipm_size_t size = 0;
    size += 13 * (N + 1) * sizeof(int);
    size = (size + 8 - 1) / 8 * 8;  // make multiple of 8 bytes
    return size;
}


void d_ocp_qp_dim_create(int N, struct d_ocp_qp_dim* dim, void* mem) {
    // loop index
    int ii;

    // zero memory (to avoid corrupted memory like e.g. NaN)
    hpipm_size_t memsize = d_ocp_qp_dim_memsize(N);
    hpipm_zero_memset(memsize, mem);
    char* c_ptr = mem;

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
    // nbxe
    dim->nbxe = (int*) c_ptr;
    c_ptr += (N + 1) * sizeof(int);
    // nbue
    dim->nbue = (int*) c_ptr;
    c_ptr += (N + 1) * sizeof(int);
    // nge
    dim->nge = (int*) c_ptr;
    c_ptr += (N + 1) * sizeof(int);

    // N
    dim->N = N;

    // initialize dims to zero by default
    // XXX already zero-initialized while zeroing all memory

    dim->memsize = d_ocp_qp_dim_memsize(N);
}


void d_ocp_qp_dim_copy_all(struct d_ocp_qp_dim* dim_orig, struct d_ocp_qp_dim* dim_dest) {

#if defined(RUNTIME_CHECKS)
    if (dim_orig->N != dim_dest->N) {
        printf("\nerror: d_ocp_qp_dim_copy_all: dim_orig->N != dim_dest->N\n");
        exit(1);
    }
#endif

    // loop index
    int ii;

    // N
    int N = dim_orig->N;

    // copy qp dim
    for (ii = 0; ii <= N; ii++)
        dim_dest->nx[ii] = dim_orig->nx[ii];
    for (ii = 0; ii <= N; ii++)
        dim_dest->nu[ii] = dim_orig->nu[ii];
    for (ii = 0; ii <= N; ii++)
        dim_dest->nb[ii] = dim_orig->nb[ii];
    for (ii = 0; ii <= N; ii++)
        dim_dest->nbx[ii] = dim_orig->nbx[ii];
    for (ii = 0; ii <= N; ii++)
        dim_dest->nbu[ii] = dim_orig->nbu[ii];
    for (ii = 0; ii <= N; ii++)
        dim_dest->ng[ii] = dim_orig->ng[ii];
    for (ii = 0; ii <= N; ii++)
        dim_dest->ns[ii] = dim_orig->ns[ii];
    for (ii = 0; ii <= N; ii++)
        dim_dest->nsbx[ii] = dim_orig->nsbx[ii];
    for (ii = 0; ii <= N; ii++)
        dim_dest->nsbu[ii] = dim_orig->nsbu[ii];
    for (ii = 0; ii <= N; ii++)
        dim_dest->nsg[ii] = dim_orig->nsg[ii];
    for (ii = 0; ii <= N; ii++)
        dim_dest->nbxe[ii] = dim_orig->nbxe[ii];
    for (ii = 0; ii <= N; ii++)
        dim_dest->nbue[ii] = dim_orig->nbue[ii];
    for (ii = 0; ii <= N; ii++)
        dim_dest->nge[ii] = dim_orig->nge[ii];
}


// TODO deprecated and remove ???
void d_ocp_qp_dim_set_all(int* nx, int* nu, int* nbx, int* nbu, int* ng, int* nsbx, int* nsbu, int* nsg, struct d_ocp_qp_dim* dim) {

    // loop index
    int ii;

    // N
    int N = dim->N;

    // copy qp dim
    for (ii = 0; ii <= N; ii++)
        dim->nx[ii] = nx[ii];
    for (ii = 0; ii <= N; ii++)
        dim->nu[ii] = nu[ii];
    for (ii = 0; ii <= N; ii++)
        dim->nb[ii] = nbx[ii] + nbu[ii];
    for (ii = 0; ii <= N; ii++)
        dim->nbx[ii] = nbx[ii];
    for (ii = 0; ii <= N; ii++)
        dim->nbu[ii] = nbu[ii];
    for (ii = 0; ii <= N; ii++)
        dim->ng[ii] = ng[ii];
    for (ii = 0; ii <= N; ii++)
        dim->ns[ii] = nsbx[ii] + nsbu[ii] + nsg[ii];
    for (ii = 0; ii <= N; ii++)
        dim->nsbx[ii] = nsbx[ii];
    for (ii = 0; ii <= N; ii++)
        dim->nsbu[ii] = nsbu[ii];
    for (ii = 0; ii <= N; ii++)
        dim->nsg[ii] = nsg[ii];
}


void d_ocp_qp_dim_set(char* field_name, int stage, int value, struct d_ocp_qp_dim* dim) {
    if (hpipm_strcmp(field_name, "nx")) {
        d_ocp_qp_dim_set_nx(stage, value, dim);
    } else if (hpipm_strcmp(field_name, "nu")) {
        d_ocp_qp_dim_set_nu(stage, value, dim);
    } else if (hpipm_strcmp(field_name, "nbx")) {
        d_ocp_qp_dim_set_nbx(stage, value, dim);
    } else if (hpipm_strcmp(field_name, "nbu")) {
        d_ocp_qp_dim_set_nbu(stage, value, dim);
    } else if (hpipm_strcmp(field_name, "ng")) {
        d_ocp_qp_dim_set_ng(stage, value, dim);
    } else if (hpipm_strcmp(field_name, "ns")) {
        d_ocp_qp_dim_set_ns(stage, value, dim);
    } else if (hpipm_strcmp(field_name, "nsbx")) {
        d_ocp_qp_dim_set_nsbx(stage, value, dim);
    } else if (hpipm_strcmp(field_name, "nsbu")) {
        d_ocp_qp_dim_set_nsbu(stage, value, dim);
    } else if (hpipm_strcmp(field_name, "nsg")) {
        d_ocp_qp_dim_set_nsg(stage, value, dim);
    } else if (hpipm_strcmp(field_name, "nbxe")) {
        d_ocp_qp_dim_set_nbxe(stage, value, dim);
    } else if (hpipm_strcmp(field_name, "nbue")) {
        d_ocp_qp_dim_set_nbue(stage, value, dim);
    } else if (hpipm_strcmp(field_name, "nge")) {
        d_ocp_qp_dim_set_nge(stage, value, dim);
    } else {
        printf("error: d_ocp_qp_dim_set: wrong field %s\n", field_name);
        exit(1);
    }
}


void d_ocp_qp_dim_set_nx(int stage, int value, struct d_ocp_qp_dim* dim) {
    dim->nx[stage] = value;
}


void d_ocp_qp_dim_set_nu(int stage, int value, struct d_ocp_qp_dim* dim) {
    dim->nu[stage] = value;
}


void d_ocp_qp_dim_set_nbx(int stage, int value, struct d_ocp_qp_dim* dim) {
    dim->nbx[stage] = value;
    dim->nb[stage] = dim->nbx[stage] + dim->nbu[stage];
}


void d_ocp_qp_dim_set_nbu(int stage, int value, struct d_ocp_qp_dim* dim) {
    dim->nbu[stage] = value;
    dim->nb[stage] = dim->nbx[stage] + dim->nbu[stage];
}


void d_ocp_qp_dim_set_ng(int stage, int value, struct d_ocp_qp_dim* dim) {
    dim->ng[stage] = value;
}


void d_ocp_qp_dim_set_ns(int stage, int value, struct d_ocp_qp_dim* dim) {
    dim->ns[stage] = value;
}


void d_ocp_qp_dim_set_nsbx(int stage, int value, struct d_ocp_qp_dim* dim) {
    dim->nsbx[stage] = value;
    dim->ns[stage] = dim->nsbx[stage] + dim->nsbu[stage] + dim->nsg[stage];
}


void d_ocp_qp_dim_set_nsbu(int stage, int value, struct d_ocp_qp_dim* dim) {
    dim->nsbu[stage] = value;
    dim->ns[stage] = dim->nsbx[stage] + dim->nsbu[stage] + dim->nsg[stage];
}


void d_ocp_qp_dim_set_nsg(int stage, int value, struct d_ocp_qp_dim* dim) {
    dim->nsg[stage] = value;
    dim->ns[stage] = dim->nsbx[stage] + dim->nsbu[stage] + dim->nsg[stage];
}


void d_ocp_qp_dim_set_nbxe(int stage, int value, struct d_ocp_qp_dim* dim) {
    dim->nbxe[stage] = value;
}


void d_ocp_qp_dim_set_nbue(int stage, int value, struct d_ocp_qp_dim* dim) {
    dim->nbue[stage] = value;
}


void d_ocp_qp_dim_set_nge(int stage, int value, struct d_ocp_qp_dim* dim) {
    dim->nge[stage] = value;
}


void d_ocp_qp_dim_get(struct d_ocp_qp_dim* dim, char* field_name, int stage, int* value) {
    if (hpipm_strcmp(field_name, "nx")) {
        d_ocp_qp_dim_get_nx(dim, stage, value);
    } else if (hpipm_strcmp(field_name, "nu")) {
        d_ocp_qp_dim_get_nu(dim, stage, value);
    } else if (hpipm_strcmp(field_name, "nbx")) {
        d_ocp_qp_dim_get_nbx(dim, stage, value);
    } else if (hpipm_strcmp(field_name, "nbu")) {
        d_ocp_qp_dim_get_nbu(dim, stage, value);
    } else if (hpipm_strcmp(field_name, "ng")) {
        d_ocp_qp_dim_get_ng(dim, stage, value);
    } else if (hpipm_strcmp(field_name, "ns")) {
        d_ocp_qp_dim_get_ns(dim, stage, value);
    } else if (hpipm_strcmp(field_name, "nsbx")) {
        d_ocp_qp_dim_get_nsbx(dim, stage, value);
    } else if (hpipm_strcmp(field_name, "nsbu")) {
        d_ocp_qp_dim_get_nsbu(dim, stage, value);
    } else if (hpipm_strcmp(field_name, "nsg")) {
        d_ocp_qp_dim_get_nsg(dim, stage, value);
    } else if (hpipm_strcmp(field_name, "nbxe")) {
        d_ocp_qp_dim_get_nbxe(dim, stage, value);
    } else if (hpipm_strcmp(field_name, "nbue")) {
        d_ocp_qp_dim_get_nbue(dim, stage, value);
    } else if (hpipm_strcmp(field_name, "nge")) {
        d_ocp_qp_dim_get_nge(dim, stage, value);
    } else {
        printf("error: d_ocp_qp_dim_get: wrong field %s\n", field_name);
        exit(1);
    }
}


void d_ocp_qp_dim_get_N(struct d_ocp_qp_dim* dim, int* value) {
    *value = dim->N;
}


void d_ocp_qp_dim_get_nx(struct d_ocp_qp_dim* dim, int stage, int* value) {
    *value = dim->nx[stage];
}


void d_ocp_qp_dim_get_nu(struct d_ocp_qp_dim* dim, int stage, int* value) {
    *value = dim->nu[stage];
}


void d_ocp_qp_dim_get_nbx(struct d_ocp_qp_dim* dim, int stage, int* value) {
    *value = dim->nbx[stage];
}


void d_ocp_qp_dim_get_nbu(struct d_ocp_qp_dim* dim, int stage, int* value) {
    *value = dim->nbu[stage];
}


void d_ocp_qp_dim_get_ng(struct d_ocp_qp_dim* dim, int stage, int* value) {
    *value = dim->ng[stage];
}


void d_ocp_qp_dim_get_ns(struct d_ocp_qp_dim* dim, int stage, int* value) {
    *value = dim->ns[stage];
}


void d_ocp_qp_dim_get_nsbx(struct d_ocp_qp_dim* dim, int stage, int* value) {
    *value = dim->nsbx[stage];
}


void d_ocp_qp_dim_get_nsbu(struct d_ocp_qp_dim* dim, int stage, int* value) {
    *value = dim->nsbu[stage];
}


void d_ocp_qp_dim_get_nsg(struct d_ocp_qp_dim* dim, int stage, int* value) {
    *value = dim->nsg[stage];
}


void d_ocp_qp_dim_get_nbxe(struct d_ocp_qp_dim* dim, int stage, int* value) {
    *value = dim->nbxe[stage];
}


void d_ocp_qp_dim_get_nbue(struct d_ocp_qp_dim* dim, int stage, int* value) {
    *value = dim->nbue[stage];
}


void d_ocp_qp_dim_get_nge(struct d_ocp_qp_dim* dim, int stage, int* value) {
    *value = dim->nge[stage];
}
