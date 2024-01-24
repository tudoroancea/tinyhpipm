#include <stdio.h>
#include <stdlib.h>

#include "hpipm/blas.h"
#include "hpipm/common.h"

#include "hpipm/cond/d_cond.h"
#include "hpipm/cond/d_cond_aux.h"
#include "hpipm/cond/d_cond_qcqp.h"
#include "hpipm/cond/d_part_cond_qcqp.h"
#include "hpipm/dense/d_dense_qcqp.h"
#include "hpipm/dense/d_dense_qcqp_sol.h"
#include "hpipm/dense/d_dense_qp.h"
#include "hpipm/dense/d_dense_qp_sol.h"
#include "hpipm/ocp/d_ocp_qcqp.h"
#include "hpipm/ocp/d_ocp_qcqp_dim.h"
#include "hpipm/ocp/d_ocp_qcqp_sol.h"
#include "hpipm/ocp/d_ocp_qp.h"
#include "hpipm/ocp/d_ocp_qp_dim.h"
#include "hpipm/ocp/d_ocp_qp_sol.h"


void d_part_cond_qcqp_compute_block_size(int N, int N2, int* block_size) {

    int ii;

    int bs0 = N / N2;  // (floor) size of small blocks

    // the first blocks have size bs0+1
    for (ii = 0; ii < N - N2 * bs0; ii++)
        block_size[ii] = bs0 + 1;
    // the following blocks have size bs0
    for (; ii < N2; ii++)
        block_size[ii] = bs0;
    // the last block has size 0
    block_size[N2] = 0;
}


void d_part_cond_qcqp_compute_dim(struct d_ocp_qcqp_dim* ocp_dim, int* block_size, struct d_ocp_qcqp_dim* part_dense_dim) {

    // TODO run time check on sum(block_size) = N

    int N = ocp_dim->N;
    int* nx = ocp_dim->nx;
    int* nu = ocp_dim->nu;
    int* nb = ocp_dim->nb;
    int* nbx = ocp_dim->nbx;
    int* nbu = ocp_dim->nbu;
    int* ng = ocp_dim->ng;
    int* nq = ocp_dim->nq;
    int* ns = ocp_dim->ns;
    int* nsbx = ocp_dim->nsbx;
    int* nsbu = ocp_dim->nsbu;
    int* nsg = ocp_dim->nsg;
    int* nsq = ocp_dim->nsq;

    int N2 = part_dense_dim->N;
    //	int *nx2 = part_dense_dim->nx;
    //	int *nu2 = part_dense_dim->nu;
    //	int *nb2 = part_dense_dim->nb;
    //	int *nbx2 = part_dense_dim->nbx;
    //	int *nbu2 = part_dense_dim->nbu;
    //	int *ng2 = part_dense_dim->ng;
    //	int *nq2 = part_dense_dim->nq;
    //	int *ns2 = part_dense_dim->ns;
    //	int *nsbx2 = part_dense_dim->nsbx;
    //	int *nsbu2 = part_dense_dim->nsbu;
    //	int *nsg2 = part_dense_dim->nsg;
    //	int *nsq2 = part_dense_dim->nsq;

    int nx2, nu2, nb2, nbx2, nbu2, ng2, nq2, ns2, nsbu2, nsbx2, nsg2, nsq2;

    int ii, jj;

    // TODO equality constraints !!!!!!!!!

    int nbb;  // box constr that remain box constr
    int nbg;  // box constr that becomes general constr
    int N_tmp = 0;  // temporary sum of block size
    // first stages
    for (ii = 0; ii < N2; ii++) {
        nx2 = nx[N_tmp + 0];
        nu2 = nu[N_tmp + 0];
        nbx2 = nbx[N_tmp + 0];
        nbu2 = nbu[N_tmp + 0];
        nb2 = nb[N_tmp + 0];
        ng2 = ng[N_tmp + 0];
        nq2 = nq[N_tmp + 0];
        ns2 = ns[N_tmp + 0];
        nsbx2 = nsbx[N_tmp + 0];
        nsbu2 = nsbu[N_tmp + 0];
        nsg2 = nsg[N_tmp + 0];
        nsq2 = nsq[N_tmp + 0];
        for (jj = 1; jj < block_size[ii]; jj++) {
            nx2 += 0;
            nu2 += nu[N_tmp + jj];
            nbx2 += 0;
            nbu2 += nbu[N_tmp + jj];
            nb2 += nbu[N_tmp + jj];
            ng2 += ng[N_tmp + jj] + nbx[N_tmp + jj];
            nq2 += nq[N_tmp + jj];
            ns2 += ns[N_tmp + jj];
            nsbx2 += 0;
            nsbu2 += nsbu[N_tmp + jj];
            nsg2 += nsg[N_tmp + jj] + nsbx[N_tmp + jj];
            nsq2 += nsq[N_tmp + jj];
        }
        N_tmp += block_size[ii];
        // XXX must use setters to correctly set qp ones too !
        d_ocp_qcqp_dim_set_nx(ii, nx2, part_dense_dim);
        d_ocp_qcqp_dim_set_nu(ii, nu2, part_dense_dim);
        d_ocp_qcqp_dim_set_nbx(ii, nbx2, part_dense_dim);
        d_ocp_qcqp_dim_set_nbu(ii, nbu2, part_dense_dim);
        d_ocp_qcqp_dim_set_ng(ii, ng2, part_dense_dim);
        d_ocp_qcqp_dim_set_nq(ii, nq2, part_dense_dim);
        d_ocp_qcqp_dim_set_ns(ii, ns2, part_dense_dim);
        d_ocp_qcqp_dim_set_nsbx(ii, nsbx2, part_dense_dim);
        d_ocp_qcqp_dim_set_nsbu(ii, nsbu2, part_dense_dim);
        d_ocp_qcqp_dim_set_nsg(ii, nsg2, part_dense_dim);
        d_ocp_qcqp_dim_set_nsq(ii, nsq2, part_dense_dim);
    }
    // last stage: condense also following stage
    ii = N2;
    nx2 = nx[N_tmp + 0];
    nu2 = nu[N_tmp + 0];
    nbx2 = nbx[N_tmp + 0];
    nbu2 = nbu[N_tmp + 0];
    nb2 = nb[N_tmp + 0];
    ng2 = ng[N_tmp + 0];
    nq2 = nq[N_tmp + 0];
    ns2 = ns[N_tmp + 0];
    nsbx2 = nsbx[N_tmp + 0];
    nsbu2 = nsbu[N_tmp + 0];
    nsg2 = nsg[N_tmp + 0];
    nsq2 = nsq[N_tmp + 0];
    for (jj = 1; jj < block_size[ii] + 1; jj++) {
        nx2 += 0;
        nu2 += nu[N_tmp + jj];
        nbx2 += 0;
        nbu2 += nbu[N_tmp + jj];
        nb2 += nbu[N_tmp + jj];
        ng2 += ng[N_tmp + jj] + nbx[N_tmp + jj];
        nq2 += nq[N_tmp + jj];
        ns2 += ns[N_tmp + jj];
        nsbx2 += 0;
        nsbu2 += nsbu[N_tmp + jj];
        //		nsbx2 = nsbx[N_tmp+0];
        //		nsbu2 = nsbu[N_tmp+0];
        nsg2 += nsg[N_tmp + jj] + nsbx[N_tmp + jj];
        nsq2 += nsq[N_tmp + jj];
    }
    // XXX must use setters to correctly set qp ones too !
    d_ocp_qcqp_dim_set_nx(ii, nx2, part_dense_dim);
    d_ocp_qcqp_dim_set_nu(ii, nu2, part_dense_dim);
    d_ocp_qcqp_dim_set_nbx(ii, nbx2, part_dense_dim);
    d_ocp_qcqp_dim_set_nbu(ii, nbu2, part_dense_dim);
    d_ocp_qcqp_dim_set_ng(ii, ng2, part_dense_dim);
    d_ocp_qcqp_dim_set_nq(ii, nq2, part_dense_dim);
    d_ocp_qcqp_dim_set_ns(ii, ns2, part_dense_dim);
    d_ocp_qcqp_dim_set_nsbx(ii, nsbx2, part_dense_dim);
    d_ocp_qcqp_dim_set_nsbu(ii, nsbu2, part_dense_dim);
    d_ocp_qcqp_dim_set_nsg(ii, nsg2, part_dense_dim);
    d_ocp_qcqp_dim_set_nsq(ii, nsq2, part_dense_dim);
}


hpipm_size_t d_part_cond_qcqp_arg_memsize(int N2) {

    int ii;

    hpipm_size_t size = 0;

    size += (N2 + 1) * sizeof(struct d_cond_qcqp_arg);

    for (ii = 0; ii <= N2; ii++) {

        size += d_cond_qcqp_arg_memsize();
    }

    size = (size + 63) / 64 * 64;  // make multiple of typical cache line size
    size += 1 * 64;  // align once to typical cache line size

    return size;
}


void d_part_cond_qcqp_arg_create(int N2, struct d_part_cond_qcqp_arg* part_cond_arg, void* mem) {

    int ii;

    // cond workspace struct
    struct d_cond_qcqp_arg* cws_ptr = mem;
    part_cond_arg->cond_arg = cws_ptr;
    cws_ptr += N2 + 1;

    // align to typical cache line size
    hpipm_size_t s_ptr = (hpipm_size_t) cws_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;

    char* c_ptr = (char*) s_ptr;

    for (ii = 0; ii <= N2; ii++) {

        d_cond_qcqp_arg_create(part_cond_arg->cond_arg + ii, c_ptr);
        c_ptr += (part_cond_arg->cond_arg + ii)->memsize;
    }

    part_cond_arg->N2 = N2;

    part_cond_arg->memsize = d_part_cond_qcqp_arg_memsize(N2);

#if defined(RUNTIME_CHECKS)
    if (c_ptr > ((char*) mem) + part_cond_arg->memsize) {
        printf("\nCreate_cond_qcqp_ocp2ocp_arg: outside memory bounds!\n\n");
        exit(1);
    }
#endif
}


void d_part_cond_qcqp_arg_set_default(struct d_part_cond_qcqp_arg* part_cond_arg) {

    int ii;

    int N2 = part_cond_arg->N2;

    for (ii = 0; ii <= N2; ii++) {

        d_cond_qcqp_arg_set_default(part_cond_arg->cond_arg + ii);
        d_cond_qcqp_arg_set_cond_last_stage(0, part_cond_arg->cond_arg + ii);
    }
    // cond_last_stage at last stage
    d_cond_qcqp_arg_set_cond_last_stage(1, part_cond_arg->cond_arg + N2);
}


void d_part_cond_qcqp_arg_set_ric_alg(int ric_alg, struct d_part_cond_qcqp_arg* part_cond_arg) {

    int ii;

    int N2 = part_cond_arg->N2;

    for (ii = 0; ii <= N2; ii++) {
        d_cond_qcqp_arg_set_ric_alg(ric_alg, part_cond_arg->cond_arg + ii);
    }
}


hpipm_size_t d_part_cond_qcqp_ws_memsize(struct d_ocp_qcqp_dim* ocp_dim, int* block_size, struct d_ocp_qcqp_dim* part_dense_dim, struct d_part_cond_qcqp_arg* part_cond_arg) {

    struct d_ocp_qcqp_dim tmp_ocp_qcqp_dim;
    struct d_ocp_qp_dim tmp_ocp_qp_dim;

    int ii;

    int N = ocp_dim->N;
    int N2 = part_dense_dim->N;

    hpipm_size_t size = 0;

    size += (N2 + 1) * sizeof(struct d_cond_qcqp_ws);

    int N_tmp = 0;  // temporary sum of horizons
    for (ii = 0; ii <= N2; ii++) {

        // alias ocp_qcqp_dim
        tmp_ocp_qcqp_dim.N = block_size[ii];
        tmp_ocp_qcqp_dim.nx = ocp_dim->nx + N_tmp;
        tmp_ocp_qcqp_dim.nu = ocp_dim->nu + N_tmp;
        tmp_ocp_qcqp_dim.nbx = ocp_dim->nbx + N_tmp;
        tmp_ocp_qcqp_dim.nbu = ocp_dim->nbu + N_tmp;
        tmp_ocp_qcqp_dim.nb = ocp_dim->nb + N_tmp;
        tmp_ocp_qcqp_dim.ng = ocp_dim->ng + N_tmp;
        tmp_ocp_qcqp_dim.nq = ocp_dim->nq + N_tmp;
        tmp_ocp_qcqp_dim.nsbx = ocp_dim->nsbx + N_tmp;
        tmp_ocp_qcqp_dim.nsbu = ocp_dim->nsbu + N_tmp;
        tmp_ocp_qcqp_dim.nsg = ocp_dim->nsg + N_tmp;
        tmp_ocp_qcqp_dim.nsq = ocp_dim->nsq + N_tmp;
        tmp_ocp_qcqp_dim.ns = ocp_dim->ns + N_tmp;

        // alias ocp_qcqp_dim
        tmp_ocp_qp_dim.N = block_size[ii];
        tmp_ocp_qp_dim.nx = ocp_dim->qp_dim->nx + N_tmp;
        tmp_ocp_qp_dim.nu = ocp_dim->qp_dim->nu + N_tmp;
        tmp_ocp_qp_dim.nbx = ocp_dim->qp_dim->nbx + N_tmp;
        tmp_ocp_qp_dim.nbu = ocp_dim->qp_dim->nbu + N_tmp;
        tmp_ocp_qp_dim.nb = ocp_dim->qp_dim->nb + N_tmp;
        tmp_ocp_qp_dim.ng = ocp_dim->qp_dim->ng + N_tmp;
        tmp_ocp_qp_dim.nsbx = ocp_dim->qp_dim->nsbx + N_tmp;
        tmp_ocp_qp_dim.nsbu = ocp_dim->qp_dim->nsbu + N_tmp;
        tmp_ocp_qp_dim.nsg = ocp_dim->qp_dim->nsg + N_tmp;
        tmp_ocp_qp_dim.ns = ocp_dim->qp_dim->ns + N_tmp;

        tmp_ocp_qcqp_dim.qp_dim = &tmp_ocp_qp_dim;

        size += d_cond_qcqp_ws_memsize(&tmp_ocp_qcqp_dim, part_cond_arg->cond_arg + ii);

        N_tmp += block_size[ii];
    }

    size = (size + 63) / 64 * 64;  // make multiple of typical cache line size
    size += 1 * 64;  // align once to typical cache line size

    return size;
}


void d_part_cond_qcqp_ws_create(struct d_ocp_qcqp_dim* ocp_dim, int* block_size, struct d_ocp_qcqp_dim* part_dense_dim, struct d_part_cond_qcqp_arg* part_cond_arg, struct d_part_cond_qcqp_ws* part_cond_ws, void* mem) {

    struct d_ocp_qcqp_dim tmp_ocp_qcqp_dim;
    struct d_ocp_qp_dim tmp_ocp_qp_dim;

    int ii;

    int N = ocp_dim->N;
    int N2 = part_dense_dim->N;

    // cond workspace struct
    struct d_cond_qcqp_ws* cws_ptr = mem;
    part_cond_ws->cond_ws = cws_ptr;
    cws_ptr += N2 + 1;

    // align to typical cache line size
    hpipm_size_t s_ptr = (hpipm_size_t) cws_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;

    char* c_ptr = (char*) s_ptr;

    int N_tmp = 0;  // temporary sum of horizons
    for (ii = 0; ii <= N2; ii++) {

        // alias ocp_qcqp_dim
        tmp_ocp_qcqp_dim.N = block_size[ii];
        tmp_ocp_qcqp_dim.nx = ocp_dim->nx + N_tmp;
        tmp_ocp_qcqp_dim.nu = ocp_dim->nu + N_tmp;
        tmp_ocp_qcqp_dim.nbx = ocp_dim->nbx + N_tmp;
        tmp_ocp_qcqp_dim.nbu = ocp_dim->nbu + N_tmp;
        tmp_ocp_qcqp_dim.nb = ocp_dim->nb + N_tmp;
        tmp_ocp_qcqp_dim.ng = ocp_dim->ng + N_tmp;
        tmp_ocp_qcqp_dim.nq = ocp_dim->nq + N_tmp;
        tmp_ocp_qcqp_dim.nsbx = ocp_dim->nsbx + N_tmp;
        tmp_ocp_qcqp_dim.nsbu = ocp_dim->nsbu + N_tmp;
        tmp_ocp_qcqp_dim.nsg = ocp_dim->nsg + N_tmp;
        tmp_ocp_qcqp_dim.nsq = ocp_dim->nsq + N_tmp;
        tmp_ocp_qcqp_dim.ns = ocp_dim->ns + N_tmp;

        // alias ocp_qcqp_dim
        tmp_ocp_qp_dim.N = block_size[ii];
        tmp_ocp_qp_dim.nx = ocp_dim->qp_dim->nx + N_tmp;
        tmp_ocp_qp_dim.nu = ocp_dim->qp_dim->nu + N_tmp;
        tmp_ocp_qp_dim.nbx = ocp_dim->qp_dim->nbx + N_tmp;
        tmp_ocp_qp_dim.nbu = ocp_dim->qp_dim->nbu + N_tmp;
        tmp_ocp_qp_dim.nb = ocp_dim->qp_dim->nb + N_tmp;
        tmp_ocp_qp_dim.ng = ocp_dim->qp_dim->ng + N_tmp;
        tmp_ocp_qp_dim.nsbx = ocp_dim->qp_dim->nsbx + N_tmp;
        tmp_ocp_qp_dim.nsbu = ocp_dim->qp_dim->nsbu + N_tmp;
        tmp_ocp_qp_dim.nsg = ocp_dim->qp_dim->nsg + N_tmp;
        tmp_ocp_qp_dim.ns = ocp_dim->qp_dim->ns + N_tmp;

        tmp_ocp_qcqp_dim.qp_dim = &tmp_ocp_qp_dim;

        d_cond_qcqp_ws_create(&tmp_ocp_qcqp_dim, part_cond_arg->cond_arg + ii, part_cond_ws->cond_ws + ii, c_ptr);
        c_ptr += (part_cond_ws->cond_ws + ii)->memsize;

        N_tmp += block_size[ii];
    }

    part_cond_ws->memsize = d_part_cond_qcqp_ws_memsize(ocp_dim, block_size, part_dense_dim, part_cond_arg);

#if defined(RUNTIME_CHECKS)
    if (c_ptr > ((char*) mem) + part_cond_ws->memsize) {
        printf("\nCreate_cond_qp_ocp2ocp: outside memory bounds!\n\n");
        exit(1);
    }
#endif
}


void d_part_cond_qcqp_cond(struct d_ocp_qcqp* ocp_qp, struct d_ocp_qcqp* part_dense_qp, struct d_part_cond_qcqp_arg* part_cond_arg, struct d_part_cond_qcqp_ws* part_cond_ws) {

    struct d_ocp_qp_dim tmp_ocp_dim;
    struct d_ocp_qp tmp_ocp_qp;

    struct d_ocp_qcqp_dim tmp_ocp_qcqp_dim;
    struct d_ocp_qcqp tmp_ocp_qcqp;

    int ii;

    int N = ocp_qp->dim->N;
    int N2 = part_dense_qp->dim->N;
    int bs;  // horizon of current block

    int N_tmp = 0;  // temporary sum of horizons
    for (ii = 0; ii <= N2; ii++) {

        bs = part_cond_ws->cond_ws[ii].qp_ws->bs;

        // alias ocp_dim
        tmp_ocp_dim.N = bs;
        tmp_ocp_dim.nx = ocp_qp->dim->qp_dim->nx + N_tmp;
        tmp_ocp_dim.nu = ocp_qp->dim->qp_dim->nu + N_tmp;
        tmp_ocp_dim.nbx = ocp_qp->dim->qp_dim->nbx + N_tmp;
        tmp_ocp_dim.nbu = ocp_qp->dim->qp_dim->nbu + N_tmp;
        tmp_ocp_dim.nb = ocp_qp->dim->qp_dim->nb + N_tmp;
        tmp_ocp_dim.ng = ocp_qp->dim->qp_dim->ng + N_tmp;
        tmp_ocp_dim.nsbx = ocp_qp->dim->qp_dim->nsbx + N_tmp;
        tmp_ocp_dim.nsbu = ocp_qp->dim->qp_dim->nsbu + N_tmp;
        tmp_ocp_dim.nsg = ocp_qp->dim->qp_dim->nsg + N_tmp;
        tmp_ocp_dim.ns = ocp_qp->dim->qp_dim->ns + N_tmp;
        // TODO equality constraints !!!!!!!!!!!!!!!!!!!!!!!!!!!

        // alias ocp_qp
        tmp_ocp_qp.dim = &tmp_ocp_dim;
        tmp_ocp_qp.idxb = ocp_qp->idxb + N_tmp;
        tmp_ocp_qp.BAbt = ocp_qp->BAbt + N_tmp;
        tmp_ocp_qp.b = ocp_qp->b + N_tmp;
        tmp_ocp_qp.RSQrq = ocp_qp->RSQrq + N_tmp;
        tmp_ocp_qp.rqz = ocp_qp->rqz + N_tmp;
        tmp_ocp_qp.DCt = ocp_qp->DCt + N_tmp;
        tmp_ocp_qp.d = ocp_qp->d + N_tmp;
        tmp_ocp_qp.d_mask = ocp_qp->d_mask + N_tmp;
        tmp_ocp_qp.Z = ocp_qp->Z + N_tmp;
        tmp_ocp_qp.idxs_rev = ocp_qp->idxs_rev + N_tmp;

        d_cond_BAbt(&tmp_ocp_qp, part_dense_qp->BAbt + ii, part_dense_qp->b + ii, part_cond_arg->cond_arg[ii].qp_arg, part_cond_ws->cond_ws[ii].qp_ws);

        d_cond_RSQrq(&tmp_ocp_qp, part_dense_qp->RSQrq + ii, part_dense_qp->rqz + ii, part_cond_arg->cond_arg[ii].qp_arg, part_cond_ws->cond_ws[ii].qp_ws);

        d_cond_DCtd(&tmp_ocp_qp, part_dense_qp->idxb[ii], part_dense_qp->DCt + ii, part_dense_qp->d + ii, part_dense_qp->d_mask + ii, part_dense_qp->idxs_rev[ii], part_dense_qp->Z + ii, part_dense_qp->rqz + ii, part_cond_arg->cond_arg[ii].qp_arg, part_cond_ws->cond_ws[ii].qp_ws);

        // alias ocp_dim
        tmp_ocp_qcqp_dim.N = bs;
        tmp_ocp_qcqp_dim.nx = ocp_qp->dim->nx + N_tmp;
        tmp_ocp_qcqp_dim.nu = ocp_qp->dim->nu + N_tmp;
        tmp_ocp_qcqp_dim.nbx = ocp_qp->dim->nbx + N_tmp;
        tmp_ocp_qcqp_dim.nbu = ocp_qp->dim->nbu + N_tmp;
        tmp_ocp_qcqp_dim.nb = ocp_qp->dim->nb + N_tmp;
        tmp_ocp_qcqp_dim.ng = ocp_qp->dim->ng + N_tmp;
        tmp_ocp_qcqp_dim.nq = ocp_qp->dim->nq + N_tmp;
        tmp_ocp_qcqp_dim.nsbx = ocp_qp->dim->nsbx + N_tmp;
        tmp_ocp_qcqp_dim.nsbu = ocp_qp->dim->nsbu + N_tmp;
        tmp_ocp_qcqp_dim.nsg = ocp_qp->dim->nsg + N_tmp;
        tmp_ocp_qcqp_dim.nsq = ocp_qp->dim->nsq + N_tmp;
        tmp_ocp_qcqp_dim.ns = ocp_qp->dim->ns + N_tmp;

        // alias ocp_qp
        tmp_ocp_qcqp.dim = &tmp_ocp_qcqp_dim;
        tmp_ocp_qcqp.idxb = ocp_qp->idxb + N_tmp;
        tmp_ocp_qcqp.BAbt = ocp_qp->BAbt + N_tmp;
        tmp_ocp_qcqp.b = ocp_qp->b + N_tmp;
        tmp_ocp_qcqp.RSQrq = ocp_qp->RSQrq + N_tmp;
        tmp_ocp_qcqp.rqz = ocp_qp->rqz + N_tmp;
        tmp_ocp_qcqp.DCt = ocp_qp->DCt + N_tmp;
        tmp_ocp_qcqp.d = ocp_qp->d + N_tmp;
        tmp_ocp_qcqp.d_mask = ocp_qp->d_mask + N_tmp;
        tmp_ocp_qcqp.Z = ocp_qp->Z + N_tmp;
        tmp_ocp_qcqp.idxs_rev = ocp_qp->idxs_rev + N_tmp;
        tmp_ocp_qcqp.Hq = ocp_qp->Hq + N_tmp;
        tmp_ocp_qcqp.Hq_nzero = ocp_qp->Hq_nzero + N_tmp;

        d_cond_qcqp_qc(&tmp_ocp_qcqp, part_dense_qp->Hq[ii], part_dense_qp->Hq_nzero[ii], part_dense_qp->DCt + ii, part_dense_qp->d + ii, part_cond_arg->cond_arg + ii, part_cond_ws->cond_ws + ii);

        N_tmp += bs;
    }

#if 0
	// copy last stage
	int *nx = ocp_qp->dim->nx;
	int *nu = ocp_qp->dim->nu;
	int *nb = ocp_qp->dim->nb;
	int *ng = ocp_qp->dim->ng;
	int *ns = ocp_qp->dim->ns;

	dgecp(nu[N]+nx[N]+1, nu[N]+nx[N], ocp_qp->RSQrq+N, 0, 0, part_dense_qp->RSQrq+N2, 0, 0);
	dveccp(nu[N]+nx[N], ocp_qp->rq+N, 0, part_dense_qp->rq+N2, 0);
	dgecp(nu[N]+nx[N], ng[N], ocp_qp->DCt+N, 0, 0, part_dense_qp->DCt+N2, 0, 0);
	dveccp(2*nb[N]+2*ng[N], ocp_qp->d+N, 0, part_dense_qp->d+N2, 0);
	for(ii=0; ii<nb[N]; ii++) part_dense_qp->idxb[N2][ii] = ocp_qp->idxb[N][ii];
	dveccp(2*ns[N], ocp_qp->Z+N, 0, part_dense_qp->Z+N2, 0);
	dveccp(2*ns[N], ocp_qp->z+N, 0, part_dense_qp->z+N2, 0);
	for(ii=0; ii<ns[N]; ii++) part_dense_qp->idxs_rev[N2][ii] = ocp_qp->idxs_rev[N][ii];
#endif
}


void d_part_cond_qcqp_cond_lhs(struct d_ocp_qcqp* ocp_qp, struct d_ocp_qcqp* part_dense_qp, struct d_part_cond_qcqp_arg* part_cond_arg, struct d_part_cond_qcqp_ws* part_cond_ws) {

    struct d_ocp_qp_dim tmp_ocp_dim;
    struct d_ocp_qp tmp_ocp_qp;

    struct d_ocp_qcqp_dim tmp_ocp_qcqp_dim;
    struct d_ocp_qcqp tmp_ocp_qcqp;

    int ii;

    int N = ocp_qp->dim->N;
    int N2 = part_dense_qp->dim->N;
    int bs;  // horizon of current block

    int N_tmp = 0;  // temporary sum of horizons
    for (ii = 0; ii <= N2; ii++) {

        bs = part_cond_ws->cond_ws[ii].qp_ws->bs;

        // alias ocp_dim
        tmp_ocp_dim.N = bs;
        tmp_ocp_dim.nx = ocp_qp->dim->qp_dim->nx + N_tmp;
        tmp_ocp_dim.nu = ocp_qp->dim->qp_dim->nu + N_tmp;
        tmp_ocp_dim.nbx = ocp_qp->dim->qp_dim->nbx + N_tmp;
        tmp_ocp_dim.nbu = ocp_qp->dim->qp_dim->nbu + N_tmp;
        tmp_ocp_dim.nb = ocp_qp->dim->qp_dim->nb + N_tmp;
        tmp_ocp_dim.ng = ocp_qp->dim->qp_dim->ng + N_tmp;
        tmp_ocp_dim.nsbx = ocp_qp->dim->qp_dim->nsbx + N_tmp;
        tmp_ocp_dim.nsbu = ocp_qp->dim->qp_dim->nsbu + N_tmp;
        tmp_ocp_dim.nsg = ocp_qp->dim->qp_dim->nsg + N_tmp;
        tmp_ocp_dim.ns = ocp_qp->dim->qp_dim->ns + N_tmp;

        // alias ocp_qp
        tmp_ocp_qp.dim = &tmp_ocp_dim;
        tmp_ocp_qp.idxb = ocp_qp->idxb + N_tmp;
        tmp_ocp_qp.BAbt = ocp_qp->BAbt + N_tmp;
        tmp_ocp_qp.b = ocp_qp->b + N_tmp;
        tmp_ocp_qp.RSQrq = ocp_qp->RSQrq + N_tmp;
        tmp_ocp_qp.rqz = ocp_qp->rqz + N_tmp;
        tmp_ocp_qp.DCt = ocp_qp->DCt + N_tmp;
        tmp_ocp_qp.d = ocp_qp->d + N_tmp;
        tmp_ocp_qp.d_mask = ocp_qp->d_mask + N_tmp;
        tmp_ocp_qp.Z = ocp_qp->Z + N_tmp;
        tmp_ocp_qp.idxs_rev = ocp_qp->idxs_rev + N_tmp;

        d_cond_BAt(&tmp_ocp_qp, part_dense_qp->BAbt + ii, part_cond_arg->cond_arg[ii].qp_arg, part_cond_ws->cond_ws[ii].qp_ws);

        d_cond_RSQ(&tmp_ocp_qp, part_dense_qp->RSQrq + ii, part_cond_arg->cond_arg[ii].qp_arg, part_cond_ws->cond_ws[ii].qp_ws);

        d_cond_DCt(&tmp_ocp_qp, part_dense_qp->idxb[ii], part_dense_qp->DCt + ii, part_dense_qp->idxs_rev[ii], part_dense_qp->Z + ii, part_cond_arg->cond_arg[ii].qp_arg, part_cond_ws->cond_ws[ii].qp_ws);

        // alias ocp_dim
        tmp_ocp_qcqp_dim.N = bs;
        tmp_ocp_qcqp_dim.nx = ocp_qp->dim->nx + N_tmp;
        tmp_ocp_qcqp_dim.nu = ocp_qp->dim->nu + N_tmp;
        tmp_ocp_qcqp_dim.nbx = ocp_qp->dim->nbx + N_tmp;
        tmp_ocp_qcqp_dim.nbu = ocp_qp->dim->nbu + N_tmp;
        tmp_ocp_qcqp_dim.nb = ocp_qp->dim->nb + N_tmp;
        tmp_ocp_qcqp_dim.ng = ocp_qp->dim->ng + N_tmp;
        tmp_ocp_qcqp_dim.nq = ocp_qp->dim->nq + N_tmp;
        tmp_ocp_qcqp_dim.nsbx = ocp_qp->dim->nsbx + N_tmp;
        tmp_ocp_qcqp_dim.nsbu = ocp_qp->dim->nsbu + N_tmp;
        tmp_ocp_qcqp_dim.nsg = ocp_qp->dim->nsg + N_tmp;
        tmp_ocp_qcqp_dim.nsq = ocp_qp->dim->nsq + N_tmp;
        tmp_ocp_qcqp_dim.ns = ocp_qp->dim->ns + N_tmp;

        // alias ocp_qp
        tmp_ocp_qcqp.dim = &tmp_ocp_qcqp_dim;
        tmp_ocp_qcqp.idxb = ocp_qp->idxb + N_tmp;
        tmp_ocp_qcqp.BAbt = ocp_qp->BAbt + N_tmp;
        tmp_ocp_qcqp.b = ocp_qp->b + N_tmp;
        tmp_ocp_qcqp.RSQrq = ocp_qp->RSQrq + N_tmp;
        tmp_ocp_qcqp.rqz = ocp_qp->rqz + N_tmp;
        tmp_ocp_qcqp.DCt = ocp_qp->DCt + N_tmp;
        tmp_ocp_qcqp.d = ocp_qp->d + N_tmp;
        tmp_ocp_qcqp.d_mask = ocp_qp->d_mask + N_tmp;
        tmp_ocp_qcqp.Z = ocp_qp->Z + N_tmp;
        tmp_ocp_qcqp.idxs_rev = ocp_qp->idxs_rev + N_tmp;
        tmp_ocp_qcqp.Hq = ocp_qp->Hq + N_tmp;
        tmp_ocp_qcqp.Hq_nzero = ocp_qp->Hq_nzero + N_tmp;

        d_cond_qcqp_qc_lhs(&tmp_ocp_qcqp, part_dense_qp->Hq[ii], part_dense_qp->Hq_nzero[ii], part_dense_qp->DCt + ii, part_cond_arg->cond_arg + ii, part_cond_ws->cond_ws + ii);

        N_tmp += bs;
    }

#if 0
	// copy last stage
	int *nx = ocp_qp->dim->nx;
	int *nu = ocp_qp->dim->nu;
	int *nb = ocp_qp->dim->nb;
	int *ng = ocp_qp->dim->ng;
	int *ns = ocp_qp->dim->ns;

	dgecp(nu[N]+nx[N]+1, nu[N]+nx[N], ocp_qp->RSQrq+N, 0, 0, part_dense_qp->RSQrq+N2, 0, 0);
	dveccp(nu[N]+nx[N], ocp_qp->rq+N, 0, part_dense_qp->rq+N2, 0);
	dgecp(nu[N]+nx[N], ng[N], ocp_qp->DCt+N, 0, 0, part_dense_qp->DCt+N2, 0, 0);
	dveccp(2*nb[N]+2*ng[N], ocp_qp->d+N, 0, part_dense_qp->d+N2, 0);
	for(ii=0; ii<nb[N]; ii++) part_dense_qp->idxb[N2][ii] = ocp_qp->idxb[N][ii];
	dveccp(2*ns[N], ocp_qp->Z+N, 0, part_dense_qp->Z+N2, 0);
	dveccp(2*ns[N], ocp_qp->z+N, 0, part_dense_qp->z+N2, 0);
	for(ii=0; ii<ns[N]; ii++) part_dense_qp->idxs_rev[N2][ii] = ocp_qp->idxs_rev[N][ii];
#endif
}


void d_part_cond_qcqp_cond_rhs(struct d_ocp_qcqp* ocp_qp, struct d_ocp_qcqp* part_dense_qp, struct d_part_cond_qcqp_arg* part_cond_arg, struct d_part_cond_qcqp_ws* part_cond_ws) {

    struct d_ocp_qp_dim tmp_ocp_dim;
    struct d_ocp_qp tmp_ocp_qp;

    struct d_ocp_qcqp_dim tmp_ocp_qcqp_dim;
    struct d_ocp_qcqp tmp_ocp_qcqp;

    int ii;

    int N = ocp_qp->dim->N;
    int N2 = part_dense_qp->dim->N;
    int bs;  // horizon of current block

    int N_tmp = 0;  // temporary sum of horizons
    for (ii = 0; ii <= N2; ii++) {

        bs = part_cond_ws->cond_ws[ii].qp_ws->bs;

        // alias ocp_dim
        tmp_ocp_dim.N = bs;
        tmp_ocp_dim.nx = ocp_qp->dim->nx + N_tmp;
        tmp_ocp_dim.nu = ocp_qp->dim->nu + N_tmp;
        tmp_ocp_dim.nbx = ocp_qp->dim->nbx + N_tmp;
        tmp_ocp_dim.nbu = ocp_qp->dim->nbu + N_tmp;
        tmp_ocp_dim.nb = ocp_qp->dim->nb + N_tmp;
        tmp_ocp_dim.ng = ocp_qp->dim->ng + N_tmp;
        tmp_ocp_dim.nsbx = ocp_qp->dim->nsbx + N_tmp;
        tmp_ocp_dim.nsbu = ocp_qp->dim->nsbu + N_tmp;
        tmp_ocp_dim.nsg = ocp_qp->dim->nsg + N_tmp;
        tmp_ocp_dim.ns = ocp_qp->dim->ns + N_tmp;

        // alias ocp_qp
        tmp_ocp_qp.dim = &tmp_ocp_dim;
        tmp_ocp_qp.idxb = ocp_qp->idxb + N_tmp;
        tmp_ocp_qp.BAbt = ocp_qp->BAbt + N_tmp;
        tmp_ocp_qp.b = ocp_qp->b + N_tmp;
        tmp_ocp_qp.RSQrq = ocp_qp->RSQrq + N_tmp;
        tmp_ocp_qp.rqz = ocp_qp->rqz + N_tmp;
        tmp_ocp_qp.DCt = ocp_qp->DCt + N_tmp;
        tmp_ocp_qp.d = ocp_qp->d + N_tmp;
        tmp_ocp_qp.d_mask = ocp_qp->d_mask + N_tmp;
        tmp_ocp_qp.Z = ocp_qp->Z + N_tmp;
        tmp_ocp_qp.idxs_rev = ocp_qp->idxs_rev + N_tmp;

        d_cond_b(&tmp_ocp_qp, part_dense_qp->b + ii, part_cond_arg->cond_arg[ii].qp_arg, part_cond_ws->cond_ws[ii].qp_ws);

        d_cond_rq(&tmp_ocp_qp, part_dense_qp->rqz + ii, part_cond_arg->cond_arg[ii].qp_arg, part_cond_ws->cond_ws[ii].qp_ws);

        d_cond_d(&tmp_ocp_qp, part_dense_qp->d + ii, part_dense_qp->d_mask + ii, part_dense_qp->rqz + ii, part_cond_arg->cond_arg[ii].qp_arg, part_cond_ws->cond_ws[ii].qp_ws);

        // alias ocp_dim
        tmp_ocp_qcqp_dim.N = bs;
        tmp_ocp_qcqp_dim.nx = ocp_qp->dim->nx + N_tmp;
        tmp_ocp_qcqp_dim.nu = ocp_qp->dim->nu + N_tmp;
        tmp_ocp_qcqp_dim.nbx = ocp_qp->dim->nbx + N_tmp;
        tmp_ocp_qcqp_dim.nbu = ocp_qp->dim->nbu + N_tmp;
        tmp_ocp_qcqp_dim.nb = ocp_qp->dim->nb + N_tmp;
        tmp_ocp_qcqp_dim.ng = ocp_qp->dim->ng + N_tmp;
        tmp_ocp_qcqp_dim.nq = ocp_qp->dim->nq + N_tmp;
        tmp_ocp_qcqp_dim.nsbx = ocp_qp->dim->nsbx + N_tmp;
        tmp_ocp_qcqp_dim.nsbu = ocp_qp->dim->nsbu + N_tmp;
        tmp_ocp_qcqp_dim.nsg = ocp_qp->dim->nsg + N_tmp;
        tmp_ocp_qcqp_dim.nsq = ocp_qp->dim->nsq + N_tmp;
        tmp_ocp_qcqp_dim.ns = ocp_qp->dim->ns + N_tmp;

        // alias ocp_qp
        tmp_ocp_qcqp.dim = &tmp_ocp_qcqp_dim;
        tmp_ocp_qcqp.idxb = ocp_qp->idxb + N_tmp;
        tmp_ocp_qcqp.BAbt = ocp_qp->BAbt + N_tmp;
        tmp_ocp_qcqp.b = ocp_qp->b + N_tmp;
        tmp_ocp_qcqp.RSQrq = ocp_qp->RSQrq + N_tmp;
        tmp_ocp_qcqp.rqz = ocp_qp->rqz + N_tmp;
        tmp_ocp_qcqp.DCt = ocp_qp->DCt + N_tmp;
        tmp_ocp_qcqp.d = ocp_qp->d + N_tmp;
        tmp_ocp_qcqp.d_mask = ocp_qp->d_mask + N_tmp;
        tmp_ocp_qcqp.Z = ocp_qp->Z + N_tmp;
        tmp_ocp_qcqp.idxs_rev = ocp_qp->idxs_rev + N_tmp;
        tmp_ocp_qcqp.Hq = ocp_qp->Hq + N_tmp;
        tmp_ocp_qcqp.Hq_nzero = ocp_qp->Hq_nzero + N_tmp;

        d_cond_qcqp_qc_rhs(&tmp_ocp_qcqp, part_dense_qp->d + ii, part_cond_arg->cond_arg + ii, part_cond_ws->cond_ws + ii);

        N_tmp += bs;
    }

#if 0
	// copy last stage
	int *nx = ocp_qp->dim->nx;
	int *nu = ocp_qp->dim->nu;
	int *nb = ocp_qp->dim->nb;
	int *ng = ocp_qp->dim->ng;
	int *ns = ocp_qp->dim->ns;

	dveccp(nu[N]+nx[N], ocp_qp->rq+N, 0, part_dense_qp->rq+N2, 0);
	dveccp(2*nb[N]+2*ng[N], ocp_qp->d+N, 0, part_dense_qp->d+N2, 0);
	dveccp(2*ns[N], ocp_qp->z+N, 0, part_dense_qp->z+N2, 0);
#endif
}


void d_part_cond_qcqp_expand_sol(struct d_ocp_qcqp* ocp_qp, struct d_ocp_qcqp* part_dense_qp, struct d_ocp_qcqp_sol* part_dense_qp_sol, struct d_ocp_qcqp_sol* ocp_qp_sol, struct d_part_cond_qcqp_arg* part_cond_arg, struct d_part_cond_qcqp_ws* part_cond_ws) {

    struct d_ocp_qp_dim tmp_ocp_dim;
    struct d_ocp_qp tmp_ocp_qp;
    struct d_ocp_qp_sol tmp_ocp_qp_sol;
    struct d_dense_qp_sol dense_qp_sol;

    int bkp_comp_prim_sol;
    int bkp_comp_dual_sol_eq;
    int bkp_comp_dual_sol_ineq;

    int* nx = ocp_qp->dim->nx;
    int* nu = ocp_qp->dim->nu;
    int* nb = ocp_qp->dim->nb;
    int* ng = ocp_qp->dim->ng;
    int* nq = ocp_qp->dim->nq;
    int* ns = ocp_qp->dim->ns;

    int ii, jj, kk;

    int N = ocp_qp->dim->N;
    int N2 = part_dense_qp->dim->N;
    int bs;  // horizon of current block

    int N_tmp = 0;  // temporary sum of horizons
    for (ii = 0; ii <= N2; ii++) {

        bs = part_cond_ws->cond_ws[ii].qp_ws->bs;

        // alias ocp_dim
        tmp_ocp_dim.N = bs;
        tmp_ocp_dim.nx = ocp_qp->dim->qp_dim->nx + N_tmp;
        tmp_ocp_dim.nu = ocp_qp->dim->qp_dim->nu + N_tmp;
        tmp_ocp_dim.nbx = ocp_qp->dim->qp_dim->nbx + N_tmp;
        tmp_ocp_dim.nbu = ocp_qp->dim->qp_dim->nbu + N_tmp;
        tmp_ocp_dim.nb = ocp_qp->dim->qp_dim->nb + N_tmp;
        tmp_ocp_dim.ng = ocp_qp->dim->qp_dim->ng + N_tmp;
        tmp_ocp_dim.nsbx = ocp_qp->dim->qp_dim->nsbx + N_tmp;
        tmp_ocp_dim.nsbu = ocp_qp->dim->qp_dim->nsbu + N_tmp;
        tmp_ocp_dim.nsg = ocp_qp->dim->qp_dim->nsg + N_tmp;
        tmp_ocp_dim.ns = ocp_qp->dim->qp_dim->ns + N_tmp;

        // alias ocp_qp
        tmp_ocp_qp.dim = &tmp_ocp_dim;
        tmp_ocp_qp.idxb = ocp_qp->idxb + N_tmp;
        tmp_ocp_qp.BAbt = ocp_qp->BAbt + N_tmp;
        tmp_ocp_qp.b = ocp_qp->b + N_tmp;
        tmp_ocp_qp.RSQrq = ocp_qp->RSQrq + N_tmp;
        tmp_ocp_qp.rqz = ocp_qp->rqz + N_tmp;
        tmp_ocp_qp.DCt = ocp_qp->DCt + N_tmp;
        tmp_ocp_qp.d = ocp_qp->d + N_tmp;
        tmp_ocp_qp.d_mask = ocp_qp->d_mask + N_tmp;
        tmp_ocp_qp.Z = ocp_qp->Z + N_tmp;
        tmp_ocp_qp.idxs_rev = ocp_qp->idxs_rev + N_tmp;

        // alias ocp qp sol
        tmp_ocp_qp_sol.dim = &tmp_ocp_dim;
        tmp_ocp_qp_sol.ux = ocp_qp_sol->ux + N_tmp;
        tmp_ocp_qp_sol.pi = ocp_qp_sol->pi + N_tmp;
        tmp_ocp_qp_sol.lam = ocp_qp_sol->lam + N_tmp;
        tmp_ocp_qp_sol.t = ocp_qp_sol->t + N_tmp;

        // alias ocp qp sol
        dense_qp_sol.v = part_dense_qp_sol->ux + ii;
        dense_qp_sol.pi = part_dense_qp_sol->pi + ii;
        dense_qp_sol.lam = part_dense_qp_sol->lam + ii;
        dense_qp_sol.t = part_dense_qp_sol->t + ii;

        bkp_comp_prim_sol = part_cond_arg->cond_arg[ii].qp_arg->comp_prim_sol;
        bkp_comp_dual_sol_eq = part_cond_arg->cond_arg[ii].qp_arg->comp_dual_sol_eq;
        bkp_comp_dual_sol_ineq = part_cond_arg->cond_arg[ii].qp_arg->comp_dual_sol_ineq;

        part_cond_arg->cond_arg[ii].qp_arg->comp_prim_sol = 1 & bkp_comp_prim_sol;
        part_cond_arg->cond_arg[ii].qp_arg->comp_dual_sol_eq = 0 & bkp_comp_dual_sol_eq;
        part_cond_arg->cond_arg[ii].qp_arg->comp_dual_sol_ineq = 1 & bkp_comp_dual_sol_ineq;

        d_expand_sol(&tmp_ocp_qp, &dense_qp_sol, &tmp_ocp_qp_sol, part_cond_arg->cond_arg[ii].qp_arg, part_cond_ws->cond_ws[ii].qp_ws);

        // linearize quadr constr
        for (jj = N_tmp; jj <= N_tmp + bs; jj++) {
            dgecp(nu[jj] + nx[jj], ng[jj] + nq[jj], ocp_qp->DCt + jj, 0, 0, part_cond_ws->cond_ws[ii].tmp_DCt + (jj - N_tmp), 0, 0);
            for (kk = 0; kk < nq[jj]; kk++) {
                dsymv_l(nu[jj] + nx[jj], 1.0, ocp_qp->Hq[jj] + kk, 0, 0, ocp_qp_sol->ux + jj, 0, 0.0, part_cond_ws->cond_ws[ii].tmp_nuxM, 0, part_cond_ws->cond_ws[ii].tmp_nuxM, 0);
                dcolad(nu[jj] + nx[jj], 1.0, part_cond_ws->cond_ws[ii].tmp_nuxM, 0, part_cond_ws->cond_ws[ii].tmp_DCt + (jj - N_tmp), 0, ng[jj] + kk);
            }
        }

        tmp_ocp_qp.DCt = part_cond_ws->cond_ws[ii].tmp_DCt + 0;

        part_cond_arg->cond_arg[ii].qp_arg->comp_prim_sol = 0 & bkp_comp_prim_sol;
        part_cond_arg->cond_arg[ii].qp_arg->comp_dual_sol_eq = 1 & bkp_comp_dual_sol_eq;
        part_cond_arg->cond_arg[ii].qp_arg->comp_dual_sol_ineq = 0 & bkp_comp_dual_sol_ineq;

        d_expand_sol(&tmp_ocp_qp, &dense_qp_sol, &tmp_ocp_qp_sol, part_cond_arg->cond_arg[ii].qp_arg, part_cond_ws->cond_ws[ii].qp_ws);

        part_cond_arg->cond_arg[ii].qp_arg->comp_prim_sol = bkp_comp_prim_sol;
        part_cond_arg->cond_arg[ii].qp_arg->comp_dual_sol_eq = bkp_comp_dual_sol_eq;
        part_cond_arg->cond_arg[ii].qp_arg->comp_dual_sol_ineq = bkp_comp_dual_sol_ineq;

        N_tmp += bs;
    }

#if 0
	// copy last stage
	dveccp(nu[N]+nx[N]+2*ns[N], part_dense_qp_sol->ux+N2, 0, ocp_qp_sol->ux+N, 0);
	dveccp(2*nb[N]+2*ng[N]+2*ns[N], part_dense_qp_sol->lam+N2, 0, ocp_qp_sol->lam+N, 0);
	dveccp(2*nb[N]+2*ng[N]+2*ns[N], part_dense_qp_sol->t+N2, 0, ocp_qp_sol->t+N, 0);
#endif
}
