#include <stdio.h>
#include <stdlib.h>

#include "hpipm/blas.h"
#include "hpipm/common.h"
#include "hpipm/ocp/d_ocp_qcqp.h"
#include "hpipm/ocp/d_ocp_qcqp_dim.h"

hpipm_size_t d_ocp_qcqp_strsize() {
    return sizeof(struct d_ocp_qcqp);
}


hpipm_size_t d_ocp_qcqp_memsize(struct d_ocp_qcqp_dim* dim) {
    // extract dim
    int N = dim->N;
    int* nx = dim->nx;
    int* nu = dim->nu;
    int* nb = dim->nb;
    int* ng = dim->ng;
    int* nq = dim->nq;
    int* ns = dim->ns;
    int* nbue = dim->nbue;
    int* nbxe = dim->nbxe;
    int* nge = dim->nge;
    int* nqe = dim->nqe;

    // loop index
    int ii;

    // compute core qp size
    int nvt = 0;
    int net = 0;
    int nct = 0;
    int nqt = 0;
    for (ii = 0; ii < N; ii++) {
        nvt += nx[ii] + nu[ii] + 2 * ns[ii];
        net += nx[ii + 1];
        nct += 2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii];
        nqt += nq[ii];
    }
    nvt += nx[ii] + nu[ii] + 2 * ns[ii];
    nct += 2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii];
    nqt += nq[ii];

    hpipm_size_t size = 0;

    size += (N + 1) * sizeof(struct mat*);  // Hq

    //	size += 5*(N+1)*sizeof(int); // nx nu nb ng ns
    size += 4 * (N + 1) * sizeof(int*);  // idxb idxs_rev idxe Hq_nzero
    size += (2 * (N + 1) + nqt) * sizeof(struct mat);  // RSqrq DCt Hq
    size += 1 * N * sizeof(struct mat);  // BAbt
    size += 5 * (N + 1) * sizeof(struct vec);  // rqz d m Z d_mask
    size += 1 * N * sizeof(struct vec);  // b

    for (ii = 0; ii < N; ii++) {
        size += nb[ii] * sizeof(int);  // idxb
        size += (nb[ii] + ng[ii] + nq[ii]) * sizeof(int);  // idxs_rev
        size += (nbue[ii] + nbxe[ii] + nge[ii] + nqe[ii]) * sizeof(int);  // idxe
        size += nq[ii] * sizeof(int);  // Hq_nzero
        size += memsize_mat(nu[ii] + nx[ii] + 1, nx[ii + 1]);  // BAbt
        size += memsize_mat(nu[ii] + nx[ii] + 1, nu[ii] + nx[ii]);  // RSQrq
        size += memsize_mat(nu[ii] + nx[ii], ng[ii] + nq[ii]);  // DCt
        size += memsize_vec(2 * ns[ii]);  // Z
        size += nq[ii] * memsize_mat(nu[ii] + nx[ii] + 1, nu[ii] + nx[ii]);  // Hq
    }
    ii = N;
    size += nb[ii] * sizeof(int);  // idxb
    size += nq[ii] * sizeof(int);  // Hq_nzero
    size += (nb[ii] + ng[ii] + nq[ii]) * sizeof(int);  // idxs_rev
    size += (nbue[ii] + nbxe[ii] + nge[ii] + nqe[ii]) * sizeof(int);  // idxe
    size += memsize_mat(nu[ii] + nx[ii] + 1, nu[ii] + nx[ii]);  // RSQrq
    size += memsize_mat(nu[ii] + nx[ii], ng[ii] + nq[ii]);  // DCt
    size += memsize_vec(2 * ns[ii]);  // Z
    size += nq[ii] * memsize_mat(nu[ii] + nx[ii] + 1, nu[ii] + nx[ii]);  // Hq

    size += 1 * memsize_vec(nvt);  // rqz
    size += 1 * memsize_vec(net);  // b
    size += 3 * memsize_vec(nct);  // d m d_mask

    size = (size + 63) / 64 * 64;  // make multiple of typical cache line size
    size += 64;  // align to typical cache line size
    // TODO: replace with actual CACHE_LINE_SIZE variable

    return size;
}


void d_ocp_qcqp_create(struct d_ocp_qcqp_dim* dim, struct d_ocp_qcqp* qp, void* mem) {
    // loop index
    int ii, jj;

    // zero memory (to avoid corrupted memory like e.g. NaN)
    hpipm_size_t memsize = d_ocp_qcqp_memsize(dim);
    hpipm_zero_memset(memsize, mem);

    // extract dim
    int N = dim->N;
    int* nx = dim->nx;
    int* nu = dim->nu;
    int* nb = dim->nb;
    int* ng = dim->ng;
    int* nq = dim->nq;
    int* ns = dim->ns;
    int* nbxe = dim->nbxe;
    int* nbue = dim->nbue;
    int* nge = dim->nge;
    int* nqe = dim->nqe;

    // compute core qp size
    int nvt = 0;
    int net = 0;
    int nct = 0;
    for (ii = 0; ii < N; ii++) {
        nvt += nx[ii] + nu[ii] + 2 * ns[ii];
        net += nx[ii + 1];
        nct += 2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii];
    }
    nvt += nx[ii] + nu[ii] + 2 * ns[ii];
    nct += 2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii];


    // int pointer stuff
    int** ip_ptr;
    ip_ptr = (int**) mem;

    // idxb
    qp->idxb = ip_ptr;
    ip_ptr += N + 1;
    // idxs_rev
    qp->idxs_rev = ip_ptr;
    ip_ptr += N + 1;
    // idxe
    qp->idxe = ip_ptr;
    ip_ptr += N + 1;
    // Hq_nzero
    qp->Hq_nzero = ip_ptr;
    ip_ptr += N + 1;


    // matrix struct pointer stuff
    struct mat** smp_ptr = (struct mat**) ip_ptr;

    // Hq
    qp->Hq = smp_ptr;
    smp_ptr += N + 1;


    // matrix struct stuff
    struct mat* sm_ptr = (struct mat*) smp_ptr;

    // BAbt
    qp->BAbt = sm_ptr;
    sm_ptr += N;
    // RSQrq
    qp->RSQrq = sm_ptr;
    sm_ptr += N + 1;
    // DCt
    qp->DCt = sm_ptr;
    sm_ptr += N + 1;
    // Hq
    for (ii = 0; ii <= N; ii++) {
        qp->Hq[ii] = sm_ptr;
        sm_ptr += nq[ii];
    }


    // vector struct stuff
    struct vec* sv_ptr = (struct vec*) sm_ptr;

    // b
    qp->b = sv_ptr;
    sv_ptr += N;
    // rqz
    qp->rqz = sv_ptr;
    sv_ptr += N + 1;
    // d
    qp->d = sv_ptr;
    sv_ptr += N + 1;
    // d_mask
    qp->d_mask = sv_ptr;
    sv_ptr += N + 1;
    // m
    qp->m = sv_ptr;
    sv_ptr += N + 1;
    // Z
    qp->Z = sv_ptr;
    sv_ptr += N + 1;


    // integer stuff
    int* i_ptr;
    i_ptr = (int*) sv_ptr;

    // idxb
    for (ii = 0; ii <= N; ii++) {
        (qp->idxb)[ii] = i_ptr;
        i_ptr += nb[ii];
        for (jj = 0; jj < nb[ii]; jj++)
            qp->idxb[ii][jj] = 0;
    }
    // idxs_rev
    for (ii = 0; ii <= N; ii++) {
        (qp->idxs_rev)[ii] = i_ptr;
        i_ptr += nb[ii] + ng[ii] + nq[ii];
        for (jj = 0; jj < nb[ii] + ng[ii] + nq[ii]; jj++)
            qp->idxs_rev[ii][jj] = -1;
    }
    // idxe
    for (ii = 0; ii <= N; ii++) {
        (qp->idxe)[ii] = i_ptr;
        i_ptr += nbue[ii] + nbxe[ii] + nge[ii] + nqe[ii];
    }
    // Hq_nzero
    for (ii = 0; ii <= N; ii++) {
        (qp->Hq_nzero)[ii] = i_ptr;
        i_ptr += nq[ii];
        for (jj = 0; jj < nq[ii]; jj++)
            qp->Hq_nzero[ii][jj] = 0;
    }


    // align to typical cache line size
    hpipm_size_t l_ptr = (hpipm_size_t) i_ptr;
    l_ptr = (l_ptr + 63) / 64 * 64;


    // floating point stuff
    char* c_ptr;
    c_ptr = (char*) l_ptr;

    char* tmp_ptr;

    // BAbt
    for (ii = 0; ii < N; ii++) {
        create_mat(nu[ii] + nx[ii] + 1, nx[ii + 1], qp->BAbt + ii, c_ptr);
        c_ptr += (qp->BAbt + ii)->memsize;
        //		dgese(nu[ii]+nx[ii]+1, nx[ii+1], 0.0, qp->BAbt+ii, 0, 0);
    }
    // RSQrq
    for (ii = 0; ii <= N; ii++) {
        create_mat(nu[ii] + nx[ii] + 1, nu[ii] + nx[ii], qp->RSQrq + ii, c_ptr);
        c_ptr += (qp->RSQrq + ii)->memsize;
        //		dgese(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], 0.0, qp->RSQrq+ii, 0, 0);
    }
    // DCt
    for (ii = 0; ii <= N; ii++) {
        create_mat(nu[ii] + nx[ii], ng[ii] + nq[ii], qp->DCt + ii, c_ptr);
        c_ptr += (qp->DCt + ii)->memsize;
        //		dgese(nu[ii]+nx[ii], ng[ii], 0.0, qp->DCt+ii, 0, 0);
    }
    // Hq
    for (ii = 0; ii <= N; ii++) {
        for (jj = 0; jj < nq[ii]; jj++) {
            create_mat(nu[ii] + nx[ii] + 1, nu[ii] + nx[ii], &qp->Hq[ii][jj], c_ptr);
            c_ptr += qp->Hq[ii][jj].memsize;
            //			dgese(nu[ii]+nx[ii], ng[ii], 0.0, qp->DCt+ii, 0, 0);
        }
    }
    // Z
    for (ii = 0; ii <= N; ii++) {
        create_vec(2 * ns[ii], qp->Z + ii, c_ptr);
        c_ptr += (qp->Z + ii)->memsize;
        //		dvecse(2*ns[ii], 0.0, qp->Z+ii, 0);
    }
    // g
    tmp_ptr = c_ptr;
    c_ptr += memsize_vec(nvt);
    for (ii = 0; ii <= N; ii++) {
        create_vec(nu[ii] + nx[ii] + 2 * ns[ii], qp->rqz + ii, tmp_ptr);
        tmp_ptr += nu[ii] * sizeof(double);
        tmp_ptr += nx[ii] * sizeof(double);
        tmp_ptr += ns[ii] * sizeof(double);
        tmp_ptr += ns[ii] * sizeof(double);
        //		dvecse(nu[ii]+nx[ii]+2*ns[ii], 0.0, qp->rqz+ii, 0);
    }
    // b
    tmp_ptr = c_ptr;
    c_ptr += memsize_vec(net);
    for (ii = 0; ii < N; ii++) {
        create_vec(nx[ii + 1], qp->b + ii, tmp_ptr);
        tmp_ptr += nx[ii + 1] * sizeof(double);
        //		dvecse(nx[ii+1], 0.0, qp->b+ii, 0);
    }
    // d
    tmp_ptr = c_ptr;
    c_ptr += memsize_vec(nct);
    for (ii = 0; ii <= N; ii++) {
        create_vec(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], qp->d + ii, tmp_ptr);
        tmp_ptr += nb[ii] * sizeof(double);  // lb
        tmp_ptr += ng[ii] * sizeof(double);  // lg
        tmp_ptr += nq[ii] * sizeof(double);  // lq
        tmp_ptr += nb[ii] * sizeof(double);  // ub
        tmp_ptr += ng[ii] * sizeof(double);  // ug
        tmp_ptr += nq[ii] * sizeof(double);  // uq
        tmp_ptr += ns[ii] * sizeof(double);  // ls
        tmp_ptr += ns[ii] * sizeof(double);  // us
        //		dvecse(2*nb[ii]+2*ng[ii]+2*nq[ii]+2*ns[ii], 0.0, qp->d+ii, 0);
    }
    // d_mask
    tmp_ptr = c_ptr;
    c_ptr += memsize_vec(nct);
    for (ii = 0; ii <= N; ii++) {
        create_vec(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], qp->d_mask + ii, tmp_ptr);
        tmp_ptr += nb[ii] * sizeof(double);  // lb
        tmp_ptr += ng[ii] * sizeof(double);  // lg
        tmp_ptr += nq[ii] * sizeof(double);  // lq
        tmp_ptr += nb[ii] * sizeof(double);  // ub
        tmp_ptr += ng[ii] * sizeof(double);  // ug
        tmp_ptr += nq[ii] * sizeof(double);  // uq
        tmp_ptr += ns[ii] * sizeof(double);  // ls
        tmp_ptr += ns[ii] * sizeof(double);  // us
        dvecse(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], 1.0, qp->d_mask + ii, 0);
    }
    // m
    tmp_ptr = c_ptr;
    c_ptr += memsize_vec(nct);
    for (ii = 0; ii <= N; ii++) {
        create_vec(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], qp->m + ii, tmp_ptr);
        tmp_ptr += nb[ii] * sizeof(double);  // lb
        tmp_ptr += ng[ii] * sizeof(double);  // lg
        tmp_ptr += nq[ii] * sizeof(double);  // lq
        tmp_ptr += nb[ii] * sizeof(double);  // ub
        tmp_ptr += ng[ii] * sizeof(double);  // ug
        tmp_ptr += nq[ii] * sizeof(double);  // uq
        tmp_ptr += ns[ii] * sizeof(double);  // ls
        tmp_ptr += ns[ii] * sizeof(double);  // us
        //		dvecse(2*nb[ii]+2*ng[ii]+2*nq[ii]+2*ns[ii], 0.0, qp->m+ii, 0);
    }

    qp->dim = dim;

    qp->memsize = memsize;  // d_ocp_qcqp_memsize(dim);


#if defined(RUNTIME_CHECKS)
    if (c_ptr > ((char*) mem) + qp->memsize) {
        printf("\nerror: d_ocp_qcqp_create: outside memory bounds!\n\n");
        exit(1);
    }
#endif
}


void d_ocp_qcqp_copy_all(struct d_ocp_qcqp* qp_orig, struct d_ocp_qcqp* qp_dest) {

    // extract dim
    int N = qp_orig->dim->N;
    int* nx = qp_orig->dim->nx;
    int* nu = qp_orig->dim->nu;
    int* nb = qp_orig->dim->nb;
    int* nbx = qp_orig->dim->nbx;
    int* nbu = qp_orig->dim->nbu;
    int* ng = qp_orig->dim->ng;
    int* nq = qp_orig->dim->nq;
    int* ns = qp_orig->dim->ns;
    int* nbxe = qp_orig->dim->nbxe;
    int* nbue = qp_orig->dim->nbue;
    int* nge = qp_orig->dim->nge;
    int* nqe = qp_orig->dim->nqe;

    int ii, jj;

    // copy dim pointer
    //	qp_dest->dim = qp_orig->dim;

    for (ii = 0; ii < N; ii++) {
        dgecp(nu[ii] + nx[ii] + 1, nx[ii + 1], qp_orig->BAbt + ii, 0, 0, qp_dest->BAbt + ii, 0, 0);
        dveccp(nx[ii + 1], qp_orig->b + ii, 0, qp_dest->b + ii, 0);
    }

    for (ii = 0; ii <= N; ii++) {
        dgecp(nu[ii] + nx[ii] + 1, nu[ii] + nx[ii], qp_orig->RSQrq + ii, 0, 0, qp_dest->RSQrq + ii, 0, 0);
        dveccp(2 * ns[ii], qp_orig->Z + ii, 0, qp_dest->Z + ii, 0);
        dveccp(nu[ii] + nx[ii] + 2 * ns[ii], qp_orig->rqz + ii, 0, qp_dest->rqz + ii, 0);
        for (jj = 0; jj < nb[ii]; jj++)
            qp_dest->idxb[ii][jj] = qp_orig->idxb[ii][jj];
        dgecp(nu[ii] + nx[ii], ng[ii] + nq[ii], qp_orig->DCt + ii, 0, 0, qp_dest->DCt + ii, 0, 0);
        for (jj = 0; jj < nq[ii]; jj++)
            dgecp(nu[ii] + nx[ii], nu[ii] + nx[ii], &qp_orig->Hq[ii][jj], 0, 0, &qp_dest->Hq[ii][jj], 0, 0);
        dveccp(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], qp_orig->d + ii, 0, qp_dest->d + ii, 0);
        dveccp(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], qp_orig->d_mask + ii, 0, qp_dest->d_mask + ii, 0);
        dveccp(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], qp_orig->m + ii, 0, qp_dest->m + ii, 0);
        for (jj = 0; jj < nb[ii] + ng[ii] + nq[ii]; jj++)
            qp_dest->idxs_rev[ii][jj] = qp_orig->idxs_rev[ii][jj];
        for (jj = 0; jj < nq[ii]; jj++)
            qp_dest->Hq_nzero[ii][jj] = qp_orig->Hq_nzero[ii][jj];
        for (jj = 0; jj < nbue[ii] + nbxe[ii] + nge[ii] + nqe[ii]; jj++)
            qp_dest->idxe[ii][jj] = qp_orig->idxe[ii][jj];
    }
}


void d_ocp_qcqp_set_all_zero(struct d_ocp_qcqp* qp) {

    // extract dim
    int N = qp->dim->N;
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;
    int* nb = qp->dim->nb;
    int* nbx = qp->dim->nbx;
    int* nbu = qp->dim->nbu;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;
    int* ns = qp->dim->ns;
    int* nbxe = qp->dim->nbxe;
    int* nbue = qp->dim->nbue;
    int* nge = qp->dim->nge;
    int* nqe = qp->dim->nqe;

    int ii, jj;

    for (ii = 0; ii < N; ii++) {
        dgese(nu[ii] + nx[ii] + 1, nx[ii + 1], 0.0, qp->BAbt + ii, 0, 0);
        dvecse(nx[ii + 1], 0.0, qp->b + ii, 0);
    }

    for (ii = 0; ii <= N; ii++) {
        dgese(nu[ii] + nx[ii] + 1, nu[ii] + nx[ii], 0.0, qp->RSQrq + ii, 0, 0);
        dvecse(2 * ns[ii], 0.0, qp->Z + ii, 0);
        dvecse(nu[ii] + nx[ii] + 2 * ns[ii], 0.0, qp->rqz + ii, 0);
        for (jj = 0; jj < nb[ii]; jj++)
            qp->idxb[ii][jj] = 0;
        dgese(nu[ii] + nx[ii], ng[ii] + nq[ii], 0.0, qp->DCt + ii, 0, 0);
        for (jj = 0; jj < nq[ii]; jj++)
            dgese(nu[ii] + nx[ii], nu[ii] + nx[ii], 0.0, &qp->Hq[ii][jj], 0, 0);
        dvecse(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], 0.0, qp->d + ii, 0);
        dvecse(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], 1.0, qp->d_mask + ii, 0);
        dvecse(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], 0.0, qp->m + ii, 0);
        for (jj = 0; jj < nb[ii] + ng[ii] + nq[ii]; jj++)
            qp->idxs_rev[ii][jj] = -1;
        for (jj = 0; jj < nq[ii]; jj++)
            qp->Hq_nzero[ii][jj] = 0;
        for (jj = 0; jj < nbue[ii] + nbxe[ii] + nge[ii] + nqe[ii]; jj++)
            qp->idxe[ii][jj] = 0;
    }
}


void d_ocp_qcqp_set_rhs_zero(struct d_ocp_qcqp* qp) {

    // extract dim
    int N = qp->dim->N;
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;
    int* nb = qp->dim->nb;
    int* nbx = qp->dim->nbx;
    int* nbu = qp->dim->nbu;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;
    int* ns = qp->dim->ns;

    int ii, jj;

    for (ii = 0; ii < N; ii++) {
        dvecse(nx[ii + 1], 0.0, qp->b + ii, 0);
    }

    for (ii = 0; ii <= N; ii++) {
        dvecse(2 * ns[ii], 0.0, qp->Z + ii, 0);
        dvecse(nu[ii] + nx[ii] + 2 * ns[ii], 0.0, qp->rqz + ii, 0);
        dvecse(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], 0.0, qp->d + ii, 0);
        dvecse(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], 1.0, qp->d_mask + ii, 0);
        dvecse(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], 0.0, qp->m + ii, 0);
    }
}


void d_ocp_qcqp_set(char* field, int stage, void* value, struct d_ocp_qcqp* qp) {
    double* r_ptr;
    int* i_ptr;

    if (hpipm_strcmp(field, "A")) {
        d_ocp_qcqp_set_A(stage, value, qp);
    } else if (hpipm_strcmp(field, "B")) {
        d_ocp_qcqp_set_B(stage, value, qp);
    } else if (hpipm_strcmp(field, "Q")) {
        d_ocp_qcqp_set_Q(stage, value, qp);
    } else if (hpipm_strcmp(field, "S")) {
        d_ocp_qcqp_set_S(stage, value, qp);
    } else if (hpipm_strcmp(field, "R")) {
        d_ocp_qcqp_set_R(stage, value, qp);
    } else if (hpipm_strcmp(field, "C")) {
        d_ocp_qcqp_set_C(stage, value, qp);
    } else if (hpipm_strcmp(field, "D")) {
        d_ocp_qcqp_set_D(stage, value, qp);
    } else if (hpipm_strcmp(field, "b")) {
        d_ocp_qcqp_set_b(stage, value, qp);
    } else if (hpipm_strcmp(field, "q")) {
        d_ocp_qcqp_set_q(stage, value, qp);
    } else if (hpipm_strcmp(field, "r")) {
        d_ocp_qcqp_set_r(stage, value, qp);
    } else if (hpipm_strcmp(field, "lb")) {
        d_ocp_qcqp_set_lb(stage, value, qp);
    } else if (hpipm_strcmp(field, "lb_mask")) {
        d_ocp_qcqp_set_lb_mask(stage, value, qp);
    } else if (hpipm_strcmp(field, "lbu") | hpipm_strcmp(field, "lu")) {
        d_ocp_qcqp_set_lbu(stage, value, qp);
    } else if (hpipm_strcmp(field, "lbu_mask")) {
        d_ocp_qcqp_set_lbu_mask(stage, value, qp);
    } else if (hpipm_strcmp(field, "lbx") | hpipm_strcmp(field, "lx")) {
        d_ocp_qcqp_set_lbx(stage, value, qp);
    } else if (hpipm_strcmp(field, "lbx_mask")) {
        d_ocp_qcqp_set_lbx_mask(stage, value, qp);
    } else if (hpipm_strcmp(field, "ub")) {
        d_ocp_qcqp_set_ub(stage, value, qp);
    } else if (hpipm_strcmp(field, "ub_mask")) {
        d_ocp_qcqp_set_ub_mask(stage, value, qp);
    } else if (hpipm_strcmp(field, "ubu") | hpipm_strcmp(field, "uu")) {
        d_ocp_qcqp_set_ubu(stage, value, qp);
    } else if (hpipm_strcmp(field, "ubu_mask")) {
        d_ocp_qcqp_set_ubu_mask(stage, value, qp);
    } else if (hpipm_strcmp(field, "ubx") | hpipm_strcmp(field, "ux")) {
        d_ocp_qcqp_set_ubx(stage, value, qp);
    } else if (hpipm_strcmp(field, "ubx_mask")) {
        d_ocp_qcqp_set_ubx_mask(stage, value, qp);
    } else if (hpipm_strcmp(field, "lg")) {
        d_ocp_qcqp_set_lg(stage, value, qp);
    } else if (hpipm_strcmp(field, "lg_mask")) {
        d_ocp_qcqp_set_lg_mask(stage, value, qp);
    } else if (hpipm_strcmp(field, "ug")) {
        d_ocp_qcqp_set_ug(stage, value, qp);
    } else if (hpipm_strcmp(field, "ug_mask")) {
        d_ocp_qcqp_set_ug_mask(stage, value, qp);
    } else if (hpipm_strcmp(field, "Qq")) {
        d_ocp_qcqp_set_Qq(stage, value, qp);
    } else if (hpipm_strcmp(field, "Sq")) {
        d_ocp_qcqp_set_Sq(stage, value, qp);
    } else if (hpipm_strcmp(field, "Rq")) {
        d_ocp_qcqp_set_Rq(stage, value, qp);
    } else if (hpipm_strcmp(field, "qq")) {
        d_ocp_qcqp_set_qq(stage, value, qp);
    } else if (hpipm_strcmp(field, "rq")) {
        d_ocp_qcqp_set_rq(stage, value, qp);
    } else if (hpipm_strcmp(field, "uq")) {
        d_ocp_qcqp_set_uq(stage, value, qp);
    } else if (hpipm_strcmp(field, "uq_mask")) {
        d_ocp_qcqp_set_uq_mask(stage, value, qp);
    } else if (hpipm_strcmp(field, "Zl")) {
        d_ocp_qcqp_set_Zl(stage, value, qp);
    } else if (hpipm_strcmp(field, "Zu")) {
        d_ocp_qcqp_set_Zu(stage, value, qp);
    } else if (hpipm_strcmp(field, "zl")) {
        d_ocp_qcqp_set_zl(stage, value, qp);
    } else if (hpipm_strcmp(field, "zu")) {
        d_ocp_qcqp_set_zu(stage, value, qp);
    } else if (hpipm_strcmp(field, "lls")) {
        d_ocp_qcqp_set_lls(stage, value, qp);
    } else if (hpipm_strcmp(field, "lls_mask")) {
        d_ocp_qcqp_set_lls_mask(stage, value, qp);
    } else if (hpipm_strcmp(field, "lus")) {
        d_ocp_qcqp_set_lus(stage, value, qp);
    } else if (hpipm_strcmp(field, "lus_mask")) {
        d_ocp_qcqp_set_lus_mask(stage, value, qp);
    }
    // int
    else if (hpipm_strcmp(field, "idxb")) {
        d_ocp_qcqp_set_idxb(stage, value, qp);
    } else if (hpipm_strcmp(field, "idxbx")) {
        d_ocp_qcqp_set_idxbx(stage, value, qp);
    } else if (hpipm_strcmp(field, "Jbx") | hpipm_strcmp(field, "Jx")) {
        d_ocp_qcqp_set_Jbx(stage, value, qp);
    } else if (hpipm_strcmp(field, "idxbu")) {
        d_ocp_qcqp_set_idxbu(stage, value, qp);
    } else if (hpipm_strcmp(field, "Jbu") | hpipm_strcmp(field, "Ju")) {
        d_ocp_qcqp_set_Jbu(stage, value, qp);
    } else if (hpipm_strcmp(field, "idxs")) {
        d_ocp_qcqp_set_idxs(stage, value, qp);
    } else if (hpipm_strcmp(field, "idxs_rev")) {
        d_ocp_qcqp_set_idxs_rev(stage, value, qp);
    } else if (hpipm_strcmp(field, "Jsbu") | hpipm_strcmp(field, "Jsu")) {
        d_ocp_qcqp_set_Jsbu(stage, value, qp);
    } else if (hpipm_strcmp(field, "Jsbx") | hpipm_strcmp(field, "Jsx")) {
        d_ocp_qcqp_set_Jsbx(stage, value, qp);
    } else if (hpipm_strcmp(field, "Jsg")) {
        d_ocp_qcqp_set_Jsg(stage, value, qp);
    } else if (hpipm_strcmp(field, "Jsq")) {
        d_ocp_qcqp_set_Jsq(stage, value, qp);
    } else if (hpipm_strcmp(field, "idxe")) {
        d_ocp_qcqp_set_idxe(stage, value, qp);
    } else if (hpipm_strcmp(field, "idxbue")) {
        d_ocp_qcqp_set_idxbue(stage, value, qp);
    } else if (hpipm_strcmp(field, "idxbxe")) {
        d_ocp_qcqp_set_idxbxe(stage, value, qp);
    } else if (hpipm_strcmp(field, "idxge")) {
        d_ocp_qcqp_set_idxge(stage, value, qp);
    } else if (hpipm_strcmp(field, "idxqe")) {
        d_ocp_qcqp_set_idxqe(stage, value, qp);
    } else if (hpipm_strcmp(field, "Jbue")) {
        d_ocp_qcqp_set_Jbue(stage, value, qp);
    } else if (hpipm_strcmp(field, "Jbxe")) {
        d_ocp_qcqp_set_Jbxe(stage, value, qp);
    } else if (hpipm_strcmp(field, "Jge")) {
        d_ocp_qcqp_set_Jge(stage, value, qp);
    } else if (hpipm_strcmp(field, "Jqe")) {
        d_ocp_qcqp_set_Jqe(stage, value, qp);
    } else {
        printf("error: d_ocp_qcqp_set: wrong field name '%s'. Exiting.\n", field);
        exit(1);
    }
}


void d_ocp_qcqp_set_el(char* field, int stage, int index, void* elem, struct d_ocp_qcqp* qp) {
    double* r_ptr;
    int* i_ptr;

    // matrices
    if (hpipm_strcmp(field, "lbx") | hpipm_strcmp(field, "lx")) {
        d_ocp_qcqp_set_el_lbx(stage, index, elem, qp);
    } else if (hpipm_strcmp(field, "ubx") | hpipm_strcmp(field, "ux")) {
        d_ocp_qcqp_set_el_ubx(stage, index, elem, qp);
    } else {
        printf("error: d_ocp_qcqp_set_el: wrong field%s\n", field);
        exit(1);
    }
}


void d_ocp_qcqp_set_A(int stage, double* A, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;

    pack_tran_mat(nx[stage + 1], nx[stage], A, nx[stage + 1], qp->BAbt + stage, nu[stage], 0);
}


void d_ocp_qcqp_set_B(int stage, double* B, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;

    pack_tran_mat(nx[stage + 1], nu[stage], B, nx[stage + 1], qp->BAbt + stage, 0, 0);
}


void d_ocp_qcqp_set_b(int stage, double* b, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;

    pack_tran_mat(nx[stage + 1], 1, b, nx[stage + 1], &(qp->BAbt[stage]), nu[stage] + nx[stage], 0);  // TODO remove ???
    pack_vec(nx[stage + 1], b, 1, qp->b + stage, 0);
}


void d_ocp_qcqp_set_Q(int stage, double* Q, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;

    pack_mat(nx[stage], nx[stage], Q, nx[stage], qp->RSQrq + stage, nu[stage], nu[stage]);
}


void d_ocp_qcqp_set_S(int stage, double* S, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;

    pack_tran_mat(nu[stage], nx[stage], S, nu[stage], qp->RSQrq + stage, nu[stage], 0);
}


void d_ocp_qcqp_set_R(int stage, double* R, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;

    pack_mat(nu[stage], nu[stage], R, nu[stage], qp->RSQrq + stage, 0, 0);
}


void d_ocp_qcqp_set_q(int stage, double* q, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;

    pack_tran_mat(nx[stage], 1, q, nx[stage], &(qp->RSQrq[stage]), nu[stage] + nx[stage], nu[stage]);  // TODO remove ???
    pack_vec(nx[stage], q, 1, qp->rqz + stage, nu[stage]);
}


void d_ocp_qcqp_set_r(int stage, double* r, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;

    pack_tran_mat(nu[stage], 1, r, nu[stage], &(qp->RSQrq[stage]), nu[stage] + nx[stage], 0);  // TODO remove ???
    pack_vec(nu[stage], r, 1, qp->rqz + stage, 0);
}


void d_ocp_qcqp_set_lb(int stage, double* lb, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nb = qp->dim->nb;

    pack_vec(nb[stage], lb, 1, qp->d + stage, 0);
}


void d_ocp_qcqp_set_lb_mask(int stage, double* lb_mask, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nb = qp->dim->nb;

    pack_vec(nb[stage], lb_mask, 1, qp->d_mask + stage, 0);
}


void d_ocp_qcqp_set_lbx(int stage, double* lbx, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nbu = qp->dim->nbu;
    int* nbx = qp->dim->nbx;

    pack_vec(nbx[stage], lbx, 1, qp->d + stage, nbu[stage]);
}


void d_ocp_qcqp_set_lbx_mask(int stage, double* lbx_mask, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nbu = qp->dim->nbu;
    int* nbx = qp->dim->nbx;

    pack_vec(nbx[stage], lbx_mask, 1, qp->d_mask + stage, nbu[stage]);
}


void d_ocp_qcqp_set_el_lbx(int stage, int index, double* elem, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nbu = qp->dim->nbu;

    VECEL(qp->d + stage, nbu[stage] + index) = *elem;
}


void d_ocp_qcqp_set_lbu(int stage, double* lbu, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nbu = qp->dim->nbu;

    pack_vec(nbu[stage], lbu, 1, qp->d + stage, 0);
}


void d_ocp_qcqp_set_lbu_mask(int stage, double* lbu_mask, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nbu = qp->dim->nbu;

    pack_vec(nbu[stage], lbu_mask, 1, qp->d_mask + stage, 0);
}


void d_ocp_qcqp_set_ub(int stage, double* ub, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;

    pack_vec(nb[stage], ub, 1, qp->d + stage, nb[stage] + ng[stage] + nq[stage]);
    dvecsc(nb[stage], -1.0, qp->d + stage, nb[stage] + ng[stage] + nq[stage]);
}


void d_ocp_qcqp_set_ub_mask(int stage, double* ub_mask, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;

    pack_vec(nb[stage], ub_mask, 1, qp->d_mask + stage, nb[stage] + ng[stage] + nq[stage]);
}


void d_ocp_qcqp_set_ubx(int stage, double* ubx, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nb = qp->dim->nb;
    int* nbx = qp->dim->nbx;
    int* nbu = qp->dim->nbu;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;

    pack_vec(nbx[stage], ubx, 1, qp->d + stage, nb[stage] + ng[stage] + nq[stage] + nbu[stage]);
    dvecsc(nbx[stage], -1.0, qp->d + stage, nb[stage] + ng[stage] + nq[stage] + nbu[stage]);
}


void d_ocp_qcqp_set_ubx_mask(int stage, double* ubx_mask, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nb = qp->dim->nb;
    int* nbx = qp->dim->nbx;
    int* nbu = qp->dim->nbu;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;

    pack_vec(nbx[stage], ubx_mask, 1, qp->d_mask + stage, nb[stage] + ng[stage] + nq[stage] + nbu[stage]);
}


void d_ocp_qcqp_set_el_ubx(int stage, int index, double* elem, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nb = qp->dim->nb;
    int* nbu = qp->dim->nbu;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;

    VECEL(qp->d + stage, nb[stage] + ng[stage] + nq[stage] + nbu[stage] + index) = -*elem;
}


void d_ocp_qcqp_set_ubu(int stage, double* ubu, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nb = qp->dim->nb;
    int* nbu = qp->dim->nbu;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;

    pack_vec(nbu[stage], ubu, 1, qp->d + stage, nb[stage] + ng[stage] + nq[stage]);
    dvecsc(nbu[stage], -1.0, qp->d + stage, nb[stage] + ng[stage] + nq[stage]);
}


void d_ocp_qcqp_set_ubu_mask(int stage, double* ubu_mask, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nb = qp->dim->nb;
    int* nbu = qp->dim->nbu;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;

    pack_vec(nbu[stage], ubu_mask, 1, qp->d_mask + stage, nb[stage] + ng[stage] + nq[stage]);
}


void d_ocp_qcqp_set_idxb(int stage, int* idxb, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nb = qp->dim->nb;

    int ii;
    for (ii = 0; ii < nb[stage]; ii++)
        qp->idxb[stage][ii] = idxb[ii];
}


void d_ocp_qcqp_set_idxbx(int stage, int* idxbx, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nu = qp->dim->nu;
    int* nbx = qp->dim->nbx;
    int* nbu = qp->dim->nbu;

    int ii;
    for (ii = 0; ii < nbx[stage]; ii++) {
        qp->idxb[stage][nbu[stage] + ii] = nu[stage] + idxbx[ii];
    }
}


void d_ocp_qcqp_set_Jbx(int stage, double* Jbx, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;
    int* nbx = qp->dim->nbx;
    int* nbu = qp->dim->nbu;

    int ii, jj, jj0;
    for (ii = 0; ii < nbx[stage]; ii++) {
        jj0 = -1;
        for (jj = 0; jj < nx[stage] & jj0 == -1; jj++) {
            if (Jbx[ii + jj * nbx[stage]] != 0.0) {
                jj0 = jj;
                qp->idxb[stage][nbu[stage] + ii] = nu[stage] + jj;
            }
        }
    }
}


void d_ocp_qcqp_set_idxbu(int stage, int* idxbx, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nbu = qp->dim->nbu;

    int ii;
    for (ii = 0; ii < nbu[stage]; ii++) {
        qp->idxb[stage][ii] = idxbx[ii];
    }
}


void d_ocp_qcqp_set_Jbu(int stage, double* Jbu, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nu = qp->dim->nu;
    int* nbu = qp->dim->nbu;

    int ii, jj, jj0;
    for (ii = 0; ii < nbu[stage]; ii++) {
        jj0 = -1;
        for (jj = 0; jj < nu[stage] & jj0 == -1; jj++) {
            if (Jbu[ii + jj * nbu[stage]] != 0.0) {
                jj0 = jj;
                qp->idxb[stage][ii] = jj;
            }
        }
    }
}


void d_ocp_qcqp_set_C(int stage, double* C, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;
    int* ng = qp->dim->ng;

    pack_tran_mat(ng[stage], nx[stage], C, ng[stage], qp->DCt + stage, nu[stage], 0);
}


void d_ocp_qcqp_set_D(int stage, double* D, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nu = qp->dim->nu;
    int* ng = qp->dim->ng;

    pack_tran_mat(ng[stage], nu[stage], D, ng[stage], qp->DCt + stage, 0, 0);
}


void d_ocp_qcqp_set_lg(int stage, double* lg, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;

    pack_vec(ng[stage], lg, 1, qp->d + stage, nb[stage]);
}


void d_ocp_qcqp_set_lg_mask(int stage, double* lg_mask, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;

    pack_vec(ng[stage], lg_mask, 1, qp->d_mask + stage, nb[stage]);
}


void d_ocp_qcqp_set_ug(int stage, double* ug, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;

    pack_vec(ng[stage], ug, 1, qp->d + stage, 2 * nb[stage] + ng[stage] + nq[stage]);
    dvecsc(ng[stage], -1.0, qp->d + stage, 2 * nb[stage] + ng[stage] + nq[stage]);
}


void d_ocp_qcqp_set_ug_mask(int stage, double* ug_mask, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;

    pack_vec(ng[stage], ug_mask, 1, qp->d_mask + stage, 2 * nb[stage] + ng[stage] + nq[stage]);
}


void d_ocp_qcqp_set_Qq(int stage, double* Qq, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;
    int* nq = qp->dim->nq;

    int ii, jj;

    int nzero = 0;

    for (ii = 0; ii < nq[stage]; ii++) {
        for (jj = 0; jj < nx[stage] * nx[stage]; jj++) {
            if (Qq[ii * nx[stage] * nx[stage] + jj] != 0.0) {
                nzero = 1;
                break;
            }
        }
        pack_mat(nx[stage], nx[stage], Qq + ii * nx[stage] * nx[stage], nx[stage], &qp->Hq[stage][ii], nu[stage], nu[stage]);
        // print to check
        // blasfeo_print_dmat(nx[stage] + nu[stage], nx[stage] + nu[stage], &qp->Hq[stage][ii], 0, 0);
        if (nzero) {
            qp->Hq_nzero[stage][ii] |= 1;
        }
    }
}


void d_ocp_qcqp_set_Sq(int stage, double* Sq, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;
    int* nq = qp->dim->nq;

    int ii, jj;

    int nzero = 0;

    for (ii = 0; ii < nq[stage]; ii++) {
        for (jj = 0; jj < nu[stage] * nx[stage]; jj++) {
            if (Sq[ii * nu[stage] * nx[stage] + jj] != 0.0) {
                nzero = 1;
                break;
            }
        }
        pack_tran_mat(nu[stage], nx[stage], Sq + ii * nu[stage] * nx[stage], nu[stage], &qp->Hq[stage][ii], nu[stage], 0);
        if (nzero) {
            qp->Hq_nzero[stage][ii] |= 2;
        }
    }
}


void d_ocp_qcqp_set_Rq(int stage, double* Rq, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;
    int* nq = qp->dim->nq;

    int ii, jj;

    int nzero = 0;

    for (ii = 0; ii < nq[stage]; ii++) {
        for (jj = 0; jj < nu[stage] * nu[stage]; jj++) {
            if (Rq[ii * nu[stage] * nu[stage] + jj] != 0.0) {
                nzero = 1;
                break;
            }
        }
        pack_mat(nu[stage], nu[stage], Rq + ii * nu[stage] * nu[stage], nu[stage], &qp->Hq[stage][ii], 0, 0);
        if (nzero) {
            qp->Hq_nzero[stage][ii] |= 4;
        }
    }
}


void d_ocp_qcqp_set_qq(int stage, double* qq, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;

    int ii;

    pack_mat(nx[stage], nq[stage], qq, nx[stage], qp->DCt + stage, nu[stage], ng[stage]);
}


void d_ocp_qcqp_set_rq(int stage, double* rq, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;

    int ii;

    pack_mat(nu[stage], nq[stage], rq, nu[stage], qp->DCt + stage, 0, ng[stage]);
}


void d_ocp_qcqp_set_uq(int stage, double* value, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;

    pack_vec(nq[stage], value, 1, qp->d + stage, 2 * nb[stage] + 2 * ng[stage] + nq[stage]);
    dvecsc(nq[stage], -1.0, qp->d + stage, 2 * nb[stage] + 2 * ng[stage] + nq[stage]);
}


void d_ocp_qcqp_set_uq_mask(int stage, double* value, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;
    pack_vec(nq[stage], value, 1, qp->d_mask + stage, 2 * nb[stage] + 2 * ng[stage] + nq[stage]);
}


void d_ocp_qcqp_set_Zl(int stage, double* Zl, struct d_ocp_qcqp* qp) {
    // extract dim
    int* ns = qp->dim->ns;

    pack_vec(ns[stage], Zl, 1, qp->Z + stage, 0);
}


void d_ocp_qcqp_set_Zu(int stage, double* Zu, struct d_ocp_qcqp* qp) {
    // extract dim
    int* ns = qp->dim->ns;

    pack_vec(ns[stage], Zu, 1, qp->Z + stage, ns[stage]);
}


void d_ocp_qcqp_set_zl(int stage, double* zl, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nu = qp->dim->nu;
    int* nx = qp->dim->nx;
    int* ns = qp->dim->ns;

    pack_vec(ns[stage], zl, 1, qp->rqz + stage, nu[stage] + nx[stage]);
}


void d_ocp_qcqp_set_zu(int stage, double* zu, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nu = qp->dim->nu;
    int* nx = qp->dim->nx;
    int* ns = qp->dim->ns;

    pack_vec(ns[stage], zu, 1, qp->rqz + stage, nu[stage] + nx[stage] + ns[stage]);
}


void d_ocp_qcqp_set_idxs(int stage, int* idxs, struct d_ocp_qcqp* qp) {
    // extract dim
    int* ns = qp->dim->ns;

    int ii;
    for (ii = 0; ii < ns[stage]; ii++) {
        qp->idxs_rev[stage][idxs[ii]] = ii;
    }
}


void d_ocp_qcqp_set_idxs_rev(int stage, int* idxs_rev, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;

    int ii;
    for (ii = 0; ii < nb[stage] + ng[stage] + nq[stage]; ii++) {
        qp->idxs_rev[stage][ii] = idxs_rev[ii];
    }
}


void d_ocp_qcqp_set_Jsbu(int stage, double* Jsbu, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;
    int* nb = qp->dim->nb;
    int* nbx = qp->dim->nbx;
    int* nbu = qp->dim->nbu;
    int* ng = qp->dim->ng;
    int* ns = qp->dim->ns;

    int ii, jj, jj0, idx_tmp;
    // compute nbu part of idxs_rev
    for (ii = 0; ii < nbu[stage]; ii++) {
        jj0 = -1;
        for (jj = 0; jj < ns[stage] & jj0 == -1; jj++) {
            if (Jsbu[ii + jj * nbu[stage]] != 0.0) {
                jj0 = jj;
                qp->idxs_rev[stage][0 + ii] = jj;
            }
        }
    }
}


void d_ocp_qcqp_set_Jsbx(int stage, double* Jsbx, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;
    int* nb = qp->dim->nb;
    int* nbx = qp->dim->nbx;
    int* nbu = qp->dim->nbu;
    int* ng = qp->dim->ng;
    int* ns = qp->dim->ns;

    int ii, jj, jj0, idx_tmp;
    // compute nbx part of idxs_rev
    for (ii = 0; ii < nbx[stage]; ii++) {
        jj0 = -1;
        for (jj = 0; jj < ns[stage] & jj0 == -1; jj++) {
            if (Jsbx[ii + jj * nbx[stage]] != 0.0) {
                jj0 = jj;
                qp->idxs_rev[stage][nbu[stage] + ii] = jj;
            }
        }
    }
}


void d_ocp_qcqp_set_Jsg(int stage, double* Jsg, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;
    int* nb = qp->dim->nb;
    int* nbx = qp->dim->nbx;
    int* nbu = qp->dim->nbu;
    int* ng = qp->dim->ng;
    int* ns = qp->dim->ns;

    int ii, jj, jj0, idx_tmp;
    // compute ng part of idxs_rev
    for (ii = 0; ii < ng[stage]; ii++) {
        jj0 = -1;
        for (jj = 0; jj < ns[stage] & jj0 == -1; jj++) {
            if (Jsg[ii + jj * ng[stage]] != 0.0) {
                jj0 = jj;
                qp->idxs_rev[stage][nb[stage] + ii] = jj;
            }
        }
    }
}


void d_ocp_qcqp_set_Jsq(int stage, double* Jsq, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;
    int* nb = qp->dim->nb;
    int* nbx = qp->dim->nbx;
    int* nbu = qp->dim->nbu;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;
    int* ns = qp->dim->ns;

    int ii, jj, jj0, idx_tmp;
    // compute nq part of idxs_rev
    for (ii = 0; ii < nq[stage]; ii++) {
        jj0 = -1;
        for (jj = 0; jj < ns[stage] & jj0 == -1; jj++) {
            if (Jsq[ii + jj * nq[stage]] != 0.0) {
                jj0 = jj;
                qp->idxs_rev[stage][nb[stage] + ng[stage] + ii] = jj;
            }
        }
    }
}


void d_ocp_qcqp_set_lls(int stage, double* ls, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;
    int* ns = qp->dim->ns;

    pack_vec(ns[stage], ls, 1, qp->d + stage, 2 * nb[stage] + 2 * ng[stage] + 2 * nq[stage]);
}


void d_ocp_qcqp_set_lls_mask(int stage, double* ls_mask, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;
    int* ns = qp->dim->ns;

    pack_vec(ns[stage], ls_mask, 1, qp->d_mask + stage, 2 * nb[stage] + 2 * ng[stage] + 2 * nq[stage]);
}


void d_ocp_qcqp_set_lus(int stage, double* us, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;
    int* ns = qp->dim->ns;

    pack_vec(ns[stage], us, 1, qp->d + stage, 2 * nb[stage] + 2 * ng[stage] + 2 * nq[stage] + ns[stage]);
}


void d_ocp_qcqp_set_lus_mask(int stage, double* lus_mask, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;
    int* ns = qp->dim->ns;

    pack_vec(ns[stage], lus_mask, 1, qp->d_mask + stage, 2 * nb[stage] + 2 * ng[stage] + 2 * nq[stage] + ns[stage]);
}


void d_ocp_qcqp_set_idxe(int stage, int* idxe, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nbxe = qp->dim->nbxe;
    int* nbue = qp->dim->nbue;
    int* nge = qp->dim->nge;
    int* nqe = qp->dim->nqe;

    int ii;
    for (ii = 0; ii < nbxe[stage] + nbue[stage] + nge[stage] + nqe[stage]; ii++)
        qp->idxe[stage][ii] = idxe[ii];
}


void d_ocp_qcqp_set_idxbue(int stage, int* idxbue, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nbxe = qp->dim->nbxe;
    int* nbue = qp->dim->nbue;
    int* nge = qp->dim->nge;
    int* nqe = qp->dim->nqe;

    int ii;
    for (ii = 0; ii < nbue[stage]; ii++)
        qp->idxe[stage][ii] = idxbue[ii];
}


void d_ocp_qcqp_set_idxbxe(int stage, int* idxbxe, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nbu = qp->dim->nbu;
    int* nbxe = qp->dim->nbxe;
    int* nbue = qp->dim->nbue;
    int* nge = qp->dim->nge;
    int* nqe = qp->dim->nqe;

    int ii;
    for (ii = 0; ii < nbxe[stage]; ii++)
        qp->idxe[stage][nbue[stage] + ii] = nbu[stage] + idxbxe[ii];
}


void d_ocp_qcqp_set_idxge(int stage, int* idxge, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nbu = qp->dim->nbu;
    int* nbx = qp->dim->nbx;
    int* nbxe = qp->dim->nbxe;
    int* nbue = qp->dim->nbue;
    int* nge = qp->dim->nge;
    int* nqe = qp->dim->nqe;

    int ii;
    for (ii = 0; ii < nge[stage]; ii++)
        qp->idxe[stage][nbue[stage] + nbxe[stage] + ii] = nbu[stage] + nbx[stage] + idxge[ii];
}


void d_ocp_qcqp_set_idxqe(int stage, int* idxqe, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nbu = qp->dim->nbu;
    int* nbx = qp->dim->nbx;
    int* ng = qp->dim->ng;
    int* nbxe = qp->dim->nbxe;
    int* nbue = qp->dim->nbue;
    int* nge = qp->dim->nge;
    int* nqe = qp->dim->nqe;

    int ii;
    for (ii = 0; ii < nqe[stage]; ii++)
        qp->idxe[stage][nbue[stage] + nbxe[stage] + nge[stage] + ii] = nbu[stage] + nbx[stage] + ng[stage] + idxqe[ii];
}


void d_ocp_qcqp_set_Jbue(int stage, double* Jbue, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nbu = qp->dim->nbu;
    int* nbue = qp->dim->nbue;

    int ii, idx;
    idx = 0;
    for (ii = 0; ii < nbu[stage]; ii++) {
        if (idx < nbue[stage] | Jbue[ii] != 0.0) {
            qp->idxe[stage][idx] = ii;
            idx++;
        }
    }
}


void d_ocp_qcqp_set_Jbxe(int stage, double* Jbxe, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nbu = qp->dim->nbu;
    int* nbx = qp->dim->nbx;
    int* nbue = qp->dim->nbue;
    int* nbxe = qp->dim->nbxe;

    int ii, idx;
    idx = 0;
    for (ii = 0; ii < nbx[stage]; ii++) {
        if (idx < nbxe[stage] | Jbxe[ii] != 0.0) {
            qp->idxe[stage][nbue[stage] + idx] = nbu[stage] + ii;
            idx++;
        }
    }
}


void d_ocp_qcqp_set_Jge(int stage, double* Jge, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nbu = qp->dim->nbu;
    int* nbx = qp->dim->nbx;
    int* ng = qp->dim->ng;
    int* nbue = qp->dim->nbue;
    int* nbxe = qp->dim->nbxe;
    int* nge = qp->dim->nge;

    int ii, idx;
    idx = 0;
    for (ii = 0; ii < ng[stage]; ii++) {
        if (idx < nge[stage] | Jge[ii] != 0.0) {
            qp->idxe[stage][nbue[stage] + nbxe[stage] + idx] = nbu[stage] + nbx[stage] + ii;
            idx++;
        }
    }
}


void d_ocp_qcqp_set_Jqe(int stage, double* Jqe, struct d_ocp_qcqp* qp) {
    // extract dim
    int* nbu = qp->dim->nbu;
    int* nbx = qp->dim->nbx;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;
    int* nbue = qp->dim->nbue;
    int* nbxe = qp->dim->nbxe;
    int* nge = qp->dim->nge;
    int* nqe = qp->dim->nqe;

    int ii, idx;
    idx = 0;
    for (ii = 0; ii < nq[stage]; ii++) {
        if (idx < nqe[stage] | Jqe[ii] != 0.0) {
            qp->idxe[stage][nbue[stage] + nbxe[stage] + nge[stage] + idx] = nbu[stage] + nbx[stage] + ng[stage] + ii;
            idx++;
        }
    }
}


void d_ocp_qcqp_get(char* field, int stage, struct d_ocp_qcqp* qp, void* value) {
    // matrices
    if (hpipm_strcmp(field, "A")) {
        d_ocp_qcqp_get_A(stage, qp, value);
    }
    // vectors
    else if (hpipm_strcmp(field, "lbx") | hpipm_strcmp(field, "lx")) {
        d_ocp_qcqp_get_lbx(stage, qp, value);
    } else if (hpipm_strcmp(field, "ubx") | hpipm_strcmp(field, "ux")) {
        d_ocp_qcqp_get_ubx(stage, qp, value);
    }
    // int
    else {
        printf("error: d_ocp_qcqp_get: wrong field %s\n", field);
        exit(1);
    }
}


void d_ocp_qcqp_get_A(int stage, struct d_ocp_qcqp* qp, double* A) {
    // extract dim
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;

    unpack_tran_mat(nx[stage], nx[stage + 1], qp->BAbt + stage, nu[stage], 0, A, nx[stage + 1]);
}


void d_ocp_qcqp_get_B(int stage, struct d_ocp_qcqp* qp, double* B) {
    // extract dim
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;

    unpack_tran_mat(nu[stage], nx[stage + 1], qp->BAbt + stage, 0, 0, B, nx[stage + 1]);
}


void d_ocp_qcqp_get_b(int stage, struct d_ocp_qcqp* qp, double* b) {
    // extract dim
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;

    unpack_vec(nx[stage + 1], qp->b + stage, 0, b, 1);
}


void d_ocp_qcqp_get_Q(int stage, struct d_ocp_qcqp* qp, double* Q) {
    // extract dim
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;

    unpack_mat(nx[stage], nx[stage], qp->RSQrq + stage, nu[stage], nu[stage], Q, nx[stage]);
}


void d_ocp_qcqp_get_S(int stage, struct d_ocp_qcqp* qp, double* S) {
    // extract dim
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;

    unpack_tran_mat(nx[stage], nu[stage], qp->RSQrq + stage, nu[stage], 0, S, nu[stage]);
}


void d_ocp_qcqp_get_R(int stage, struct d_ocp_qcqp* qp, double* R) {
    // extract dim
    int* nu = qp->dim->nu;

    unpack_mat(nu[stage], nu[stage], qp->RSQrq + stage, 0, 0, R, nu[stage]);
}


void d_ocp_qcqp_get_q(int stage, struct d_ocp_qcqp* qp, double* q) {
    // extract dim
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;

    unpack_vec(nx[stage], qp->rqz + stage, nu[stage], q, 1);
}


void d_ocp_qcqp_get_r(int stage, struct d_ocp_qcqp* qp, double* r) {
    // extract dim
    int* nu = qp->dim->nu;

    unpack_vec(nu[stage], qp->rqz + stage, 0, r, 1);
}


void d_ocp_qcqp_get_lb(int stage, struct d_ocp_qcqp* qp, double* lb) {
    // extract dim
    int* nb = qp->dim->nb;

    int i;

    unpack_vec(nb[stage], qp->d + stage, 0, lb, 1);
}


void d_ocp_qcqp_get_lb_mask(int stage, struct d_ocp_qcqp* qp, double* lb_mask) {
    // extract dim
    int* nb = qp->dim->nb;

    int i;

    unpack_vec(nb[stage], qp->d_mask + stage, 0, lb_mask, 1);
}


void d_ocp_qcqp_get_ub(int stage, struct d_ocp_qcqp* qp, double* ub) {
    // extract dim
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;

    int i;

    unpack_vec(nb[stage], qp->d + stage, nb[stage] + ng[stage] + nq[stage], ub, 1);
    for (i = 0; i < nb[stage]; i++) {
        ub[i] = -ub[i];
    }
}


void d_ocp_qcqp_get_ub_mask(int stage, struct d_ocp_qcqp* qp, double* ub_mask) {
    // extract dim
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;

    int i;

    unpack_vec(nb[stage], qp->d_mask + stage, nb[stage] + ng[stage] + nq[stage], ub_mask, 1);
}


void d_ocp_qcqp_get_lbx(int stage, struct d_ocp_qcqp* qp, double* lbx) {
    // extract dim
    int* nbu = qp->dim->nbu;
    int* nbx = qp->dim->nbx;

    unpack_vec(nbx[stage], qp->d + stage, nbu[stage], lbx, 1);
}


void d_ocp_qcqp_get_lbx_mask(int stage, struct d_ocp_qcqp* qp, double* lbx_mask) {
    // extract dim
    int* nbu = qp->dim->nbu;
    int* nbx = qp->dim->nbx;

    unpack_vec(nbx[stage], qp->d_mask + stage, nbu[stage], lbx_mask, 1);
}


void d_ocp_qcqp_get_ubx(int stage, struct d_ocp_qcqp* qp, double* ubx) {
    // extract dim
    int* nb = qp->dim->nb;
    int* nbx = qp->dim->nbx;
    int* nbu = qp->dim->nbu;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;

    int i;

    unpack_vec(nbx[stage], qp->d + stage, nb[stage] + ng[stage] + nq[stage] + nbu[stage], ubx, 1);
    for (i = 0; i < nbx[stage]; i++) {
        ubx[i] = -ubx[i];
    }
}


void d_ocp_qcqp_get_ubx_mask(int stage, struct d_ocp_qcqp* qp, double* ubx_mask) {
    // extract dim
    int* nb = qp->dim->nb;
    int* nbx = qp->dim->nbx;
    int* nbu = qp->dim->nbu;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;

    int i;

    unpack_vec(nbx[stage], qp->d_mask + stage, nb[stage] + ng[stage] + nq[stage] + nbu[stage], ubx_mask, 1);
}


void d_ocp_qcqp_get_lbu(int stage, struct d_ocp_qcqp* qp, double* lbu) {
    // extract dim
    int* nbu = qp->dim->nbu;

    unpack_vec(nbu[stage], qp->d + stage, 0, lbu, 1);
}


void d_ocp_qcqp_get_lbu_mask(int stage, struct d_ocp_qcqp* qp, double* lbu_mask) {
    // extract dim
    int* nbu = qp->dim->nbu;

    unpack_vec(nbu[stage], qp->d_mask + stage, 0, lbu_mask, 1);
}


void d_ocp_qcqp_get_ubu(int stage, struct d_ocp_qcqp* qp, double* ubu) {
    // extract dim
    int* nb = qp->dim->nb;
    int* nbu = qp->dim->nbu;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;

    int i;

    unpack_vec(nbu[stage], qp->d + stage, nb[stage] + ng[stage] + nq[stage], ubu, 1);
    for (i = 0; i < nbu[stage]; i++) {
        ubu[i] = -ubu[i];
    }
}


void d_ocp_qcqp_get_ubu_mask(int stage, struct d_ocp_qcqp* qp, double* ubu_mask) {
    // extract dim
    int* nb = qp->dim->nb;
    int* nbu = qp->dim->nbu;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;

    int i;

    unpack_vec(nbu[stage], qp->d_mask + stage, nb[stage] + ng[stage] + nq[stage], ubu_mask, 1);
}


void d_ocp_qcqp_get_idxb(int stage, struct d_ocp_qcqp* qp, int* idxb) {
    // extract dim
    int* nb = qp->dim->nb;

    int ii;
    for (ii = 0; ii < nb[stage]; ii++)
        idxb[ii] = qp->idxb[stage][ii];
}


// void d_ocp_qcqp_get_idxbx(int stage, struct d_ocp_qcqp *qp, int *idxb)
//	{
//	TODO
//
//	}


// void d_ocp_qcqp_get_Jbx(int stage, struct d_ocp_qcqp *qp, int *Jbx)
//	{
//	TODO
//
//	}


// void d_ocp_qcqp_get_idxbu(int stage, struct d_ocp_qcqp *qp, int *idxbu)
//	{
//	TODO
//
//	}


// void d_ocp_qcqp_get_Jbu(int stage, struct d_ocp_qcqp *qp, int *Jbu)
//	{
//	TODO
//
//	}


void d_ocp_qcqp_get_C(int stage, struct d_ocp_qcqp* qp, double* C) {
    // extract dim
    int* nx = qp->dim->nx;
    int* nu = qp->dim->nu;
    int* ng = qp->dim->ng;

    unpack_tran_mat(nx[stage], ng[stage], qp->DCt + stage, nu[stage], 0, C, ng[stage]);
}


void d_ocp_qcqp_get_D(int stage, struct d_ocp_qcqp* qp, double* D) {
    // extract dim
    int* nu = qp->dim->nu;
    int* ng = qp->dim->ng;

    unpack_tran_mat(nu[stage], ng[stage], qp->DCt + stage, 0, 0, D, ng[stage]);
}


void d_ocp_qcqp_get_lg(int stage, struct d_ocp_qcqp* qp, double* lg) {
    // extract dim
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;

    unpack_vec(ng[stage], qp->d + stage, nb[stage], lg, 1);
}


void d_ocp_qcqp_get_lg_mask(int stage, struct d_ocp_qcqp* qp, double* lg_mask) {
    // extract dim
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;

    unpack_vec(ng[stage], qp->d_mask + stage, nb[stage], lg_mask, 1);
}


void d_ocp_qcqp_get_ug(int stage, struct d_ocp_qcqp* qp, double* ug) {
    // extract dim
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;

    int i;

    unpack_vec(ng[stage], qp->d + stage, 2 * nb[stage] + ng[stage] + nq[stage], ug, 1);
    for (i = 0; i < ng[stage]; i++) {
        ug[i] = -ug[i];
    }
}


void d_ocp_qcqp_get_ug_mask(int stage, struct d_ocp_qcqp* qp, double* ug_mask) {
    // extract dim
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;

    int i;

    unpack_vec(ng[stage], qp->d_mask + stage, 2 * nb[stage] + ng[stage] + nq[stage], ug_mask, 1);
}


void d_ocp_qcqp_get_Zl(int stage, struct d_ocp_qcqp* qp, double* Zl) {
    // extract dim
    int* ns = qp->dim->ns;

    unpack_vec(ns[stage], qp->Z + stage, 0, Zl, 1);
}


void d_ocp_qcqp_get_Zu(int stage, struct d_ocp_qcqp* qp, double* Zu) {
    // extract dim
    int* ns = qp->dim->ns;

    unpack_vec(ns[stage], qp->Z + stage, ns[stage], Zu, 1);
}


void d_ocp_qcqp_get_zl(int stage, struct d_ocp_qcqp* qp, double* zl) {
    // extract dim
    int* nu = qp->dim->nu;
    int* nx = qp->dim->nx;
    int* ns = qp->dim->ns;

    unpack_vec(ns[stage], qp->rqz + stage, nu[stage] + nx[stage], zl, 1);
}


void d_ocp_qcqp_get_zu(int stage, struct d_ocp_qcqp* qp, double* zu) {
    // extract dim
    int* nu = qp->dim->nu;
    int* nx = qp->dim->nx;
    int* ns = qp->dim->ns;

    unpack_vec(ns[stage], qp->rqz + stage, nu[stage] + nx[stage] + ns[stage], zu, 1);
}


// XXX only valid if there is one slack per softed constraint !!!
void d_ocp_qcqp_get_idxs(int stage, struct d_ocp_qcqp* qp, int* idxs) {
    // extract dim
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;
    int* ns = qp->dim->ns;

    int ii, idx_tmp;

    for (ii = 0; ii < nb[stage] + ng[stage] + nq[stage]; ii++) {
        idx_tmp = qp->idxs_rev[stage][ii];
        if (idx_tmp != -1) {
            idxs[idx_tmp] = ii;
        }
    }
}


void d_ocp_qcqp_get_idxs_rev(int stage, struct d_ocp_qcqp* qp, int* idxs_rev) {
    // extract dim
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;

    int ii;
    for (ii = 0; ii < nb[stage] + ng[stage] + nq[stage]; ii++)
        idxs_rev[ii] = qp->idxs_rev[stage][ii];
}


// void d_ocp_qcqp_get_JSBX(int stage, struct d_ocp_qcqp *qp, int *Jsbx)
//	{
//	TODO
//
//	}


// void d_ocp_qcqp_get_JSBX(int stage, struct d_ocp_qcqp *qp, int *Jsbx)
//	{
//	TODO
//
//	}


// void d_ocp_qcqp_get_JSG(int stage, struct d_ocp_qcqp *qp, int *Jsg)
//	{
//	TODO
//
//	}


void d_ocp_qcqp_get_lls(int stage, struct d_ocp_qcqp* qp, double* ls) {
    // extract dim
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;
    int* ns = qp->dim->ns;

    unpack_vec(ns[stage], qp->d + stage, 2 * nb[stage] + 2 * ng[stage] + 2 * nq[stage], ls, 1);
}


void d_ocp_qcqp_get_lls_mask(int stage, struct d_ocp_qcqp* qp, double* ls_mask) {
    // extract dim
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;
    int* ns = qp->dim->ns;

    unpack_vec(ns[stage], qp->d_mask + stage, 2 * nb[stage] + 2 * ng[stage] + 2 * nq[stage], ls_mask, 1);
}


void d_ocp_qcqp_get_lus(int stage, struct d_ocp_qcqp* qp, double* us) {
    // extract dim
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;
    int* ns = qp->dim->ns;

    int i;

    unpack_vec(ns[stage], qp->d + stage, 2 * nb[stage] + 2 * ng[stage] + 2 * nq[stage] + ns[stage], us, 1);
}


void d_ocp_qcqp_get_lus_mask(int stage, struct d_ocp_qcqp* qp, double* us_mask) {
    // extract dim
    int* nb = qp->dim->nb;
    int* ng = qp->dim->ng;
    int* nq = qp->dim->nq;
    int* ns = qp->dim->ns;

    int i;

    unpack_vec(ns[stage], qp->d_mask + stage, 2 * nb[stage] + 2 * ng[stage] + 2 * nq[stage] + ns[stage], us_mask, 1);
}
