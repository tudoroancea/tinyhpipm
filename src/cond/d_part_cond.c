#include <stdio.h>
#include <stdlib.h>

#include "tinyhpipm/blas.h"
#include "tinyhpipm/common.h"

#include "tinyhpipm/cond/d_cond.h"
#include "tinyhpipm/cond/d_cond_aux.h"
#include "tinyhpipm/cond/d_part_cond.h"
#include "tinyhpipm/dense/d_dense_qp.h"
#include "tinyhpipm/dense/d_dense_qp_sol.h"
#include "tinyhpipm/ocp/d_ocp_qp.h"
#include "tinyhpipm/ocp/d_ocp_qp_dim.h"
#include "tinyhpipm/ocp/d_ocp_qp_sol.h"


void d_part_cond_qp_compute_block_size(int N, int N2, int* block_size) {

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


void d_part_cond_qp_compute_dim(struct d_ocp_qp_dim* ocp_dim, int* block_size, struct d_ocp_qp_dim* part_dense_dim) {

    // TODO run time check on sum(block_size) = N

    int N = ocp_dim->N;
    int* nx = ocp_dim->nx;
    int* nu = ocp_dim->nu;
    int* nb = ocp_dim->nb;
    int* nbx = ocp_dim->nbx;
    int* nbu = ocp_dim->nbu;
    int* ng = ocp_dim->ng;
    int* ns = ocp_dim->ns;
    int* nsbx = ocp_dim->nsbx;
    int* nsbu = ocp_dim->nsbu;
    int* nsg = ocp_dim->nsg;

    int N2 = part_dense_dim->N;
    int* nx2 = part_dense_dim->nx;
    int* nu2 = part_dense_dim->nu;
    int* nb2 = part_dense_dim->nb;
    int* nbx2 = part_dense_dim->nbx;
    int* nbu2 = part_dense_dim->nbu;
    int* ng2 = part_dense_dim->ng;
    int* ns2 = part_dense_dim->ns;
    int* nsbx2 = part_dense_dim->nsbx;
    int* nsbu2 = part_dense_dim->nsbu;
    int* nsg2 = part_dense_dim->nsg;

    int ii, jj;

    // TODO equality constraints !!!!!!!!!

    //	int nbb; // box constr that remain box constr
    //	int nbg; // box constr that becomes general constr
    int N_tmp = 0;  // temporary sum of block size
    // first stages
    for (ii = 0; ii < N2; ii++) {
        nx2[ii] = nx[N_tmp + 0];
        nu2[ii] = nu[N_tmp + 0];
        nbx2[ii] = nbx[N_tmp + 0];
        nbu2[ii] = nbu[N_tmp + 0];
        nb2[ii] = nb[N_tmp + 0];
        ng2[ii] = ng[N_tmp + 0];
        ns2[ii] = ns[N_tmp + 0];
        nsbx2[ii] = nsbx[N_tmp + 0];
        nsbu2[ii] = nsbu[N_tmp + 0];
        nsg2[ii] = nsg[N_tmp + 0];
        for (jj = 1; jj < block_size[ii]; jj++) {
            nx2[ii] += 0;
            nu2[ii] += nu[N_tmp + jj];
            nbx2[ii] += 0;
            nbu2[ii] += nbu[N_tmp + jj];
            nb2[ii] += nbu[N_tmp + jj];
            ng2[ii] += ng[N_tmp + jj] + nbx[N_tmp + jj];
            ns2[ii] += ns[N_tmp + jj];
            nsbx2[ii] += 0;
            nsbu2[ii] += nsbu[N_tmp + jj];
            nsg2[ii] += nsg[N_tmp + jj] + nsbx[N_tmp + jj];
        }
        N_tmp += block_size[ii];
    }
    // last stage: condense also following stage
    ii = N2;
    nx2[ii] = nx[N_tmp + 0];
    nu2[ii] = nu[N_tmp + 0];
    nbx2[ii] = nbx[N_tmp + 0];
    nbu2[ii] = nbu[N_tmp + 0];
    nb2[ii] = nb[N_tmp + 0];
    ng2[ii] = ng[N_tmp + 0];
    ns2[ii] = ns[N_tmp + 0];
    nsbx2[ii] = nsbx[N_tmp + 0];
    nsbu2[ii] = nsbu[N_tmp + 0];
    nsg2[ii] = nsg[N_tmp + 0];
    for (jj = 1; jj < block_size[ii] + 1; jj++) {
        nx2[ii] += 0;
        nu2[ii] += nu[N_tmp + jj];
        nbx2[ii] += 0;
        nbu2[ii] += nbu[N_tmp + jj];
        nb2[ii] += nbu[N_tmp + jj];
        ng2[ii] += ng[N_tmp + jj] + nbx[N_tmp + jj];
        ns2[ii] += ns[N_tmp + jj];
        nsbx2[ii] += 0;
        nsbu2[ii] += nsbu[N_tmp + jj];
        //		nsbx2[ii] = nsbx[N_tmp+0];
        //		nsbu2[ii] = nsbu[N_tmp+0];
        nsg2[ii] += nsg[N_tmp + jj] + nsbx[N_tmp + jj];
    }
}


hpipm_size_t d_part_cond_qp_arg_memsize(int N2) {

    int ii;

    hpipm_size_t size = 0;

    size += (N2 + 1) * sizeof(struct d_cond_qp_arg);

    for (ii = 0; ii <= N2; ii++) {

        size += d_cond_qp_arg_memsize();
    }

    size = (size + 63) / 64 * 64;  // make multiple of typical cache line size
    size += 1 * 64;  // align once to typical cache line size

    return size;
}


void d_part_cond_qp_arg_create(int N2, struct d_part_cond_qp_arg* part_cond_arg, void* mem) {

    int ii;

    // cond workspace struct
    struct d_cond_qp_arg* cws_ptr = mem;
    part_cond_arg->cond_arg = cws_ptr;
    cws_ptr += N2 + 1;

    // align to typical cache line size
    hpipm_size_t s_ptr = (hpipm_size_t) cws_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;

    char* c_ptr = (char*) s_ptr;

    for (ii = 0; ii <= N2; ii++) {

        d_cond_qp_arg_create(part_cond_arg->cond_arg + ii, c_ptr);
        c_ptr += (part_cond_arg->cond_arg + ii)->memsize;
    }

    part_cond_arg->N2 = N2;

    part_cond_arg->memsize = d_part_cond_qp_arg_memsize(N2);

#if defined(RUNTIME_CHECKS)
    if (c_ptr > ((char*) mem) + part_cond_arg->memsize) {
        printf("\nCreate_cond_qp_ocp2ocp_arg: outside memory bounds!\n\n");
        exit(1);
    }
#endif
}


void d_part_cond_qp_arg_set_default(struct d_part_cond_qp_arg* part_cond_arg) {

    int ii;

    int N2 = part_cond_arg->N2;

    for (ii = 0; ii <= N2; ii++) {

        d_cond_qp_arg_set_default(part_cond_arg->cond_arg + ii);
        d_cond_qp_arg_set_cond_last_stage(0, part_cond_arg->cond_arg + ii);
    }
    // cond_last_stage at last stage
    d_cond_qp_arg_set_cond_last_stage(1, part_cond_arg->cond_arg + N2);
}


void d_part_cond_qp_arg_set_ric_alg(int ric_alg, struct d_part_cond_qp_arg* part_cond_arg) {

    int ii;

    int N2 = part_cond_arg->N2;

    for (ii = 0; ii <= N2; ii++) {
        d_cond_qp_arg_set_ric_alg(ric_alg, part_cond_arg->cond_arg + ii);
    }
}


void d_part_cond_qp_arg_set_comp_prim_sol(int value, struct d_part_cond_qp_arg* part_cond_arg) {

    int ii;

    int N2 = part_cond_arg->N2;

    for (ii = 0; ii <= N2; ii++) {
        d_cond_qp_arg_set_comp_prim_sol(value, part_cond_arg->cond_arg + ii);
    }
}


void d_part_cond_qp_arg_set_comp_dual_sol_eq(int value, struct d_part_cond_qp_arg* part_cond_arg) {

    int ii;

    int N2 = part_cond_arg->N2;

    for (ii = 0; ii <= N2; ii++) {
        d_cond_qp_arg_set_comp_dual_sol_eq(value, part_cond_arg->cond_arg + ii);
    }
}


void d_part_cond_qp_arg_set_comp_dual_sol_ineq(int value, struct d_part_cond_qp_arg* part_cond_arg) {

    int ii;

    int N2 = part_cond_arg->N2;

    for (ii = 0; ii <= N2; ii++) {
        d_cond_qp_arg_set_comp_dual_sol_ineq(value, part_cond_arg->cond_arg + ii);
    }
}


hpipm_size_t d_part_cond_qp_ws_memsize(struct d_ocp_qp_dim* ocp_dim, int* block_size, struct d_ocp_qp_dim* part_dense_dim, struct d_part_cond_qp_arg* part_cond_arg) {

    struct d_ocp_qp_dim tmp_ocp_dim;

    int ii;

    int N = ocp_dim->N;
    int N2 = part_dense_dim->N;

    hpipm_size_t size = 0;

    size += (N2 + 1) * sizeof(struct d_cond_qp_ws);

    int N_tmp = 0;  // temporary sum of horizons
    for (ii = 0; ii <= N2; ii++) {

        // alias ocp_dim
        tmp_ocp_dim.N = block_size[ii];
        tmp_ocp_dim.nx = ocp_dim->nx + N_tmp;
        tmp_ocp_dim.nu = ocp_dim->nu + N_tmp;
        tmp_ocp_dim.nbx = ocp_dim->nbx + N_tmp;
        tmp_ocp_dim.nbu = ocp_dim->nbu + N_tmp;
        tmp_ocp_dim.nb = ocp_dim->nb + N_tmp;
        tmp_ocp_dim.ng = ocp_dim->ng + N_tmp;
        tmp_ocp_dim.nsbx = ocp_dim->nsbx + N_tmp;
        tmp_ocp_dim.nsbu = ocp_dim->nsbu + N_tmp;
        tmp_ocp_dim.nsg = ocp_dim->nsg + N_tmp;
        tmp_ocp_dim.ns = ocp_dim->ns + N_tmp;

        size += d_cond_qp_ws_memsize(&tmp_ocp_dim, part_cond_arg->cond_arg + ii);

        N_tmp += block_size[ii];
    }

    size = (size + 63) / 64 * 64;  // make multiple of typical cache line size
    size += 1 * 64;  // align once to typical cache line size

    return size;
}


void d_part_cond_qp_ws_create(struct d_ocp_qp_dim* ocp_dim, int* block_size, struct d_ocp_qp_dim* part_dense_dim, struct d_part_cond_qp_arg* part_cond_arg, struct d_part_cond_qp_ws* part_cond_ws, void* mem) {

    struct d_ocp_qp_dim tmp_ocp_dim;

    int ii;

    int N = ocp_dim->N;
    int N2 = part_dense_dim->N;

    // cond workspace struct
    struct d_cond_qp_ws* cws_ptr = mem;
    part_cond_ws->cond_workspace = cws_ptr;
    cws_ptr += N2 + 1;

    // align to typical cache line size
    hpipm_size_t s_ptr = (hpipm_size_t) cws_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;

    char* c_ptr = (char*) s_ptr;

    int N_tmp = 0;  // temporary sum of horizons
    for (ii = 0; ii <= N2; ii++) {

        // alias ocp_dim
        tmp_ocp_dim.N = block_size[ii];
        tmp_ocp_dim.nx = ocp_dim->nx + N_tmp;
        tmp_ocp_dim.nu = ocp_dim->nu + N_tmp;
        tmp_ocp_dim.nbx = ocp_dim->nbx + N_tmp;
        tmp_ocp_dim.nbu = ocp_dim->nbu + N_tmp;
        tmp_ocp_dim.nb = ocp_dim->nb + N_tmp;
        tmp_ocp_dim.ng = ocp_dim->ng + N_tmp;
        tmp_ocp_dim.nsbx = ocp_dim->nsbx + N_tmp;
        tmp_ocp_dim.nsbu = ocp_dim->nsbu + N_tmp;
        tmp_ocp_dim.nsg = ocp_dim->nsg + N_tmp;
        tmp_ocp_dim.ns = ocp_dim->ns + N_tmp;
        // TODO equality constraints !!!!!!!!!!!!!!!!!!!!!!!!!!!

        d_cond_qp_ws_create(&tmp_ocp_dim, part_cond_arg->cond_arg + ii, part_cond_ws->cond_workspace + ii, c_ptr);
        c_ptr += (part_cond_ws->cond_workspace + ii)->memsize;

        N_tmp += block_size[ii];
    }

    part_cond_ws->memsize = d_part_cond_qp_ws_memsize(ocp_dim, block_size, part_dense_dim, part_cond_arg);

#if defined(RUNTIME_CHECKS)
    if (c_ptr > ((char*) mem) + part_cond_ws->memsize) {
        printf("\nCreate_cond_qp_ocp2ocp: outside memory bounds!\n\n");
        exit(1);
    }
#endif
}


void d_part_cond_qp_cond_lhs(struct d_ocp_qp* ocp_qp, struct d_ocp_qp* part_dense_qp, struct d_part_cond_qp_arg* part_cond_arg, struct d_part_cond_qp_ws* part_cond_ws) {

    struct d_ocp_qp_dim tmp_ocp_dim;
    struct d_ocp_qp tmp_ocp_qp;

    int ii;

    int N = ocp_qp->dim->N;
    int N2 = part_dense_qp->dim->N;
    int bs;  // horizon of current block

    int N_tmp = 0;  // temporary sum of horizons
    for (ii = 0; ii <= N2; ii++) {

        bs = part_cond_ws->cond_workspace[ii].bs;

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
        tmp_ocp_qp.diag_H_flag = ocp_qp->diag_H_flag + N_tmp;

        d_cond_BAt(&tmp_ocp_qp, part_dense_qp->BAbt + ii, part_cond_arg->cond_arg + ii, part_cond_ws->cond_workspace + ii);

        d_cond_RSQ(&tmp_ocp_qp, part_dense_qp->RSQrq + ii, part_cond_arg->cond_arg + ii, part_cond_ws->cond_workspace + ii);

        d_cond_DCt(&tmp_ocp_qp, part_dense_qp->idxb[ii], part_dense_qp->DCt + ii, part_dense_qp->idxs_rev[ii], part_dense_qp->Z + ii, part_cond_arg->cond_arg + ii, part_cond_ws->cond_workspace + ii);

        N_tmp += bs;
    }
}


void d_part_cond_qp_cond(struct d_ocp_qp* ocp_qp, struct d_ocp_qp* part_dense_qp, struct d_part_cond_qp_arg* part_cond_arg, struct d_part_cond_qp_ws* part_cond_ws) {

    struct d_ocp_qp_dim tmp_ocp_dim;
    struct d_ocp_qp tmp_ocp_qp;

    int ii;

    int N = ocp_qp->dim->N;
    int N2 = part_dense_qp->dim->N;
    int bs;  // horizon of current block

    int N_tmp = 0;  // temporary sum of horizons
    for (ii = 0; ii <= N2; ii++) {

        bs = part_cond_ws->cond_workspace[ii].bs;

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
        tmp_ocp_qp.diag_H_flag = ocp_qp->diag_H_flag + N_tmp;

        d_cond_BAbt(&tmp_ocp_qp, part_dense_qp->BAbt + ii, part_dense_qp->b + ii, part_cond_arg->cond_arg + ii, part_cond_ws->cond_workspace + ii);

        d_cond_RSQrq(&tmp_ocp_qp, part_dense_qp->RSQrq + ii, part_dense_qp->rqz + ii, part_cond_arg->cond_arg + ii, part_cond_ws->cond_workspace + ii);

        d_cond_DCtd(&tmp_ocp_qp, part_dense_qp->idxb[ii], part_dense_qp->DCt + ii, part_dense_qp->d + ii, part_dense_qp->d_mask + ii, part_dense_qp->idxs_rev[ii], part_dense_qp->Z + ii, part_dense_qp->rqz + ii, part_cond_arg->cond_arg + ii, part_cond_ws->cond_workspace + ii);

        N_tmp += bs;
    }
}


void d_part_cond_qp_cond_rhs(struct d_ocp_qp* ocp_qp, struct d_ocp_qp* part_dense_qp, struct d_part_cond_qp_arg* part_cond_arg, struct d_part_cond_qp_ws* part_cond_ws) {

    struct d_ocp_qp_dim tmp_ocp_dim;
    struct d_ocp_qp tmp_ocp_qp;

    int ii;

    int N = ocp_qp->dim->N;
    int N2 = part_dense_qp->dim->N;
    int bs;  // horizon of current block

    int N_tmp = 0;  // temporary sum of horizons
    for (ii = 0; ii <= N2; ii++) {

        bs = part_cond_ws->cond_workspace[ii].bs;

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
        tmp_ocp_qp.diag_H_flag = ocp_qp->diag_H_flag + N_tmp;

        d_cond_b(&tmp_ocp_qp, part_dense_qp->b + ii, part_cond_arg->cond_arg + ii, part_cond_ws->cond_workspace + ii);

        d_cond_rq(&tmp_ocp_qp, part_dense_qp->rqz + ii, part_cond_arg->cond_arg + ii, part_cond_ws->cond_workspace + ii);

        d_cond_d(&tmp_ocp_qp, part_dense_qp->d + ii, part_dense_qp->d_mask + ii, part_dense_qp->rqz + ii, part_cond_arg->cond_arg + ii, part_cond_ws->cond_workspace + ii);

        N_tmp += bs;
    }
}


void d_part_cond_qp_expand_sol(struct d_ocp_qp* ocp_qp, struct d_ocp_qp* part_dense_qp, struct d_ocp_qp_sol* part_dense_qp_sol, struct d_ocp_qp_sol* ocp_qp_sol, struct d_part_cond_qp_arg* part_cond_arg, struct d_part_cond_qp_ws* part_cond_ws) {

    struct d_ocp_qp_dim tmp_ocp_dim;
    struct d_ocp_qp tmp_ocp_qp;
    struct d_ocp_qp_sol tmp_ocp_qp_sol;
    struct d_dense_qp_sol dense_qp_sol;

    int* nx = ocp_qp->dim->nx;
    int* nu = ocp_qp->dim->nu;
    int* nb = ocp_qp->dim->nb;
    int* ng = ocp_qp->dim->ng;
    int* ns = ocp_qp->dim->ns;

    int ii;

    int N = ocp_qp->dim->N;
    int N2 = part_dense_qp->dim->N;
    int bs;  // horizon of current block

    int N_tmp = 0;  // temporary sum of horizons
    for (ii = 0; ii <= N2; ii++) {

        bs = part_cond_ws->cond_workspace[ii].bs;

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
        tmp_ocp_qp.diag_H_flag = ocp_qp->diag_H_flag + N_tmp;

        // alias ocp qp sol
        tmp_ocp_qp_sol.ux = ocp_qp_sol->ux + N_tmp;
        tmp_ocp_qp_sol.pi = ocp_qp_sol->pi + N_tmp;
        tmp_ocp_qp_sol.lam = ocp_qp_sol->lam + N_tmp;
        tmp_ocp_qp_sol.t = ocp_qp_sol->t + N_tmp;

        // alias ocp qp sol
        dense_qp_sol.v = part_dense_qp_sol->ux + ii;
        dense_qp_sol.pi = part_dense_qp_sol->pi + ii;
        dense_qp_sol.lam = part_dense_qp_sol->lam + ii;
        dense_qp_sol.t = part_dense_qp_sol->t + ii;

        d_expand_sol(&tmp_ocp_qp, &dense_qp_sol, &tmp_ocp_qp_sol, part_cond_arg->cond_arg + ii, part_cond_ws->cond_workspace + ii);

        N_tmp += bs;
    }
}

/************************************************
 * update cond
 ************************************************/

void d_part_cond_qp_update(int* idxc, struct d_ocp_qp* ocp_qp, struct d_ocp_qp* part_dense_qp, struct d_part_cond_qp_arg* part_cond_arg, struct d_part_cond_qp_ws* part_cond_ws) {

    struct d_ocp_qp_dim tmp_ocp_dim;
    struct d_ocp_qp tmp_ocp_qp;

    int ii;

    int N = ocp_qp->dim->N;
    int N2 = part_dense_qp->dim->N;
    int bs;  // horizon of current block

    int N_tmp = 0;  // temporary sum of horizons
    for (ii = 0; ii <= N2; ii++) {

        bs = part_cond_ws->cond_workspace[ii].bs;

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
        tmp_ocp_qp.diag_H_flag = ocp_qp->diag_H_flag + N_tmp;

        d_update_cond_BAbt(idxc + N_tmp, &tmp_ocp_qp, part_dense_qp->BAbt + ii, part_dense_qp->b + ii, part_cond_arg->cond_arg + ii, part_cond_ws->cond_workspace + ii);

        d_update_cond_RSQrq_N2nx3(idxc + N_tmp, &tmp_ocp_qp, part_dense_qp->RSQrq + ii, part_dense_qp->rqz + ii, part_cond_arg->cond_arg + ii, part_cond_ws->cond_workspace + ii);

        d_update_cond_DCtd(idxc + N_tmp, &tmp_ocp_qp, part_dense_qp->idxb[ii], part_dense_qp->DCt + ii, part_dense_qp->d + ii, part_dense_qp->idxs_rev[ii], part_dense_qp->Z + ii, part_dense_qp->rqz + ii, part_cond_arg->cond_arg + ii, part_cond_ws->cond_workspace + ii);

        N_tmp += bs;
    }
}
