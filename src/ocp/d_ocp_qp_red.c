#include <stdio.h>
#include <stdlib.h>

#include "hpipm/blas.h"
#include "hpipm/common.h"
#include "hpipm/ocp/d_ocp_qp.h"
#include "hpipm/ocp/d_ocp_qp_dim.h"
#include "hpipm/ocp/d_ocp_qp_red.h"
#include "hpipm/ocp/d_ocp_qp_sol.h"

void d_ocp_qp_dim_reduce_eq_dof(struct d_ocp_qp_dim* dim, struct d_ocp_qp_dim* dim_red) {
    int ii;

    // first stage: DOF: inputs and states
    dim_red->nu[0] = dim->nu[0] - dim->nbue[0];
    dim_red->nx[0] = dim->nx[0] - dim->nbxe[0];
    dim_red->nbu[0] = dim->nbu[0] - dim->nbue[0];
    dim_red->nbx[0] = dim->nbx[0] - dim->nbxe[0];
    dim_red->nb[0] = dim_red->nbu[0] + dim_red->nbx[0];
    dim_red->ng[0] = dim->ng[0];
    dim_red->ns[0] = dim->ns[0];
    dim_red->nsbu[0] = dim->nsbu[0];
    dim_red->nsbx[0] = dim->nsbx[0];
    dim_red->nsg[0] = dim->nsg[0];
    dim_red->nbue[0] = 0;
    dim_red->nbxe[0] = 0;
    dim_red->nge[0] = dim->nge[0];
    // other stages: DOF: inputs
    for (ii = 1; ii <= dim->N; ii++) {
        dim_red->nu[ii] = dim->nu[ii] - dim->nbue[ii];
        dim_red->nx[ii] = dim->nx[ii];
        dim_red->nbu[ii] = dim->nbu[ii] - dim->nbue[ii];
        dim_red->nbx[ii] = dim->nbx[ii];
        dim_red->nb[ii] = dim_red->nbu[ii] + dim_red->nbx[ii];
        dim_red->ng[ii] = dim->ng[ii];
        dim_red->ns[ii] = dim->ns[ii];
        dim_red->nsbu[ii] = dim->nsbu[ii];
        dim_red->nsbx[ii] = dim->nsbx[ii];
        dim_red->nsg[ii] = dim->nsg[ii];
        dim_red->nbue[ii] = 0;
        dim_red->nbxe[ii] = dim->nbxe[ii];
        dim_red->nge[ii] = dim->nge[ii];
    }
}


hpipm_size_t d_ocp_qp_reduce_eq_dof_arg_memsize() {

    return 0;
}


void d_ocp_qp_reduce_eq_dof_arg_create(struct d_ocp_qp_reduce_eq_dof_arg* arg, void* mem) {

    arg->memsize = d_ocp_qp_reduce_eq_dof_arg_memsize();
}


void d_ocp_qp_reduce_eq_dof_arg_set_default(struct d_ocp_qp_reduce_eq_dof_arg* arg) {

    arg->lam_min = 1e-16;
    arg->t_min = 1e-16;
    arg->alias_unchanged = 0;
    arg->comp_prim_sol = 1;
    arg->comp_dual_sol_eq = 1;
    arg->comp_dual_sol_ineq = 1;
}


void d_ocp_qp_reduce_eq_dof_arg_set_alias_unchanged(struct d_ocp_qp_reduce_eq_dof_arg* arg, int value) {

    arg->alias_unchanged = value;
}


void d_ocp_qp_reduce_eq_dof_arg_set_comp_prim_sol(struct d_ocp_qp_reduce_eq_dof_arg* arg, int value) {

    arg->comp_prim_sol = value;
}


void d_ocp_qp_reduce_eq_dof_arg_set_comp_dual_sol_eq(struct d_ocp_qp_reduce_eq_dof_arg* arg, int value) {

    arg->comp_dual_sol_eq = value;
}


void d_ocp_qp_reduce_eq_dof_arg_set_comp_dual_sol_ineq(struct d_ocp_qp_reduce_eq_dof_arg* arg, int value) {

    arg->comp_dual_sol_ineq = value;
}


hpipm_size_t d_ocp_qp_reduce_eq_dof_ws_memsize(struct d_ocp_qp_dim* dim) {

    int ii;

    // extract dim
    int N = dim->N;
    int* nu = dim->nu;
    int* nx = dim->nx;
    int* nb = dim->nb;
    int* ng = dim->ng;

    int nuxM = nu[0] + nx[0];
    int nbgM = nb[0] + ng[0];
    for (ii = 1; ii <= N; ii++) {
        nuxM = nu[ii] + nx[ii] > nuxM ? nu[ii] + nx[ii] : nuxM;
        nbgM = nb[ii] + ng[ii] > nbgM ? nb[ii] + ng[ii] : nbgM;
    }

    hpipm_size_t size = 0;

    size += nuxM * sizeof(int);  // e_imask_ux
    size += nbgM * sizeof(int);  // e_imask_d

    size += 3 * sizeof(struct vec);  // tmp_nuxM(0,1) tmp_nbgM

    size += 2 * memsize_vec(nuxM);  // tmp_nuxM(0,1)
    size += 1 * memsize_vec(nbgM);  // tmp_nbgM

    size = (size + 63) / 64 * 64;  // make multiple of typical cache line size
    size += 64;  // align to typical cache line size

    return size;
}


void d_ocp_qp_reduce_eq_dof_ws_create(struct d_ocp_qp_dim* dim, struct d_ocp_qp_reduce_eq_dof_ws* work, void* mem) {

    int ii;

    // zero memory (to avoid corrupted memory like e.g. NaN)
    hpipm_size_t memsize = d_ocp_qp_reduce_eq_dof_ws_memsize(dim);
    hpipm_zero_memset(memsize, mem);

    // extract dim
    int N = dim->N;
    int* nu = dim->nu;
    int* nx = dim->nx;
    int* nb = dim->nb;
    int* ng = dim->ng;

    int nuxM = nu[0] + nx[0];
    int nbgM = nb[0] + ng[0];
    for (ii = 1; ii <= N; ii++) {
        nuxM = nu[ii] + nx[ii] > nuxM ? nu[ii] + nx[ii] : nuxM;
        nbgM = nb[ii] + ng[ii] > nbgM ? nb[ii] + ng[ii] : nbgM;
    }

    // vector struct stuff
    struct vec* sv_ptr = (struct vec*) mem;

    // tmp_nuxM
    work->tmp_nuxM = sv_ptr;
    sv_ptr += 2;
    // tmp_nbgM
    work->tmp_nbgM = sv_ptr;
    sv_ptr += 1;

    // integer stuff
    int* i_ptr;
    i_ptr = (int*) sv_ptr;

    // e_imask_ux
    work->e_imask_ux = i_ptr;
    i_ptr += nuxM;
    // e_imask_d
    work->e_imask_d = i_ptr;
    i_ptr += nbgM;

    // align to typical cache line size
    hpipm_size_t l_ptr = (hpipm_size_t) i_ptr;
    l_ptr = (l_ptr + 63) / 64 * 64;

    // floating point stuff
    char* c_ptr;
    c_ptr = (char*) l_ptr;

    // tmp_nuxM
    create_vec(nuxM, work->tmp_nuxM + 0, c_ptr);
    c_ptr += (work->tmp_nuxM + 0)->memsize;
    create_vec(nuxM, work->tmp_nuxM + 1, c_ptr);
    c_ptr += (work->tmp_nuxM + 1)->memsize;
    // tmp_nbgM
    create_vec(nbgM, work->tmp_nbgM, c_ptr);
    c_ptr += (work->tmp_nbgM)->memsize;

    work->memsize = memsize;
}


void d_ocp_qp_reduce_eq_dof(struct d_ocp_qp* qp, struct d_ocp_qp* qp_red, struct d_ocp_qp_reduce_eq_dof_arg* arg, struct d_ocp_qp_reduce_eq_dof_ws* work) {

    int ii, jj, kk, idx0, idx1;

    struct d_ocp_qp_dim* dim = qp->dim;
    int N = dim->N;
    int* nx = dim->nx;
    int* nu = dim->nu;
    int* nb = dim->nb;
    int* nbx = dim->nbx;
    int* nbu = dim->nbu;
    int* ng = dim->ng;
    int* ns = dim->ns;
    int* nbue = dim->nbue;
    int* nbxe = dim->nbxe;
    int* nge = dim->nge;

    struct d_ocp_qp_dim* dim_red = qp_red->dim;
    int* nx_red = dim_red->nx;
    int* nu_red = dim_red->nu;
    int* nb_red = dim_red->nb;
    int* ng_red = dim_red->ng;
    int* ns_red = dim_red->ns;

    // TODO handle case of softed equalities !!!!!!!!!!!!!!!!

    int ne_thr;

    for (ii = 0; ii <= N; ii++) {
        if (ii == 0)
            ne_thr = nbue[ii] + nbxe[ii];
        else
            ne_thr = nbue[ii];
        if (ne_thr > 0)  // reduce inputs and/or states
        {
            dvecse(nu[ii] + nx[ii], 0.0, work->tmp_nuxM + 0, 0);
            for (jj = 0; jj < nu[ii] + nx[ii]; jj++)
                work->e_imask_ux[jj] = 0;
            for (jj = 0; jj < nbu[ii] + nbx[ii]; jj++)
                work->e_imask_d[jj] = 0;
            for (jj = 0; jj < ne_thr; jj++)  // set 1s for both inputs and states
            {
                VECEL(work->tmp_nuxM + 0, qp->idxb[ii][qp->idxe[ii][jj]]) = VECEL(qp->d + ii, qp->idxe[ii][jj]);
                work->e_imask_ux[qp->idxb[ii][qp->idxe[ii][jj]]] = 1;
                work->e_imask_d[qp->idxe[ii][jj]] = 1;
            }
            // TODO check first and last non-zero in e_mask and only multiply between them
            if (ii < N) {
                // BAt
                idx0 = 0;
                for (jj = 0; jj < nu[ii] + nx[ii]; jj++) {
                    if (work->e_imask_ux[jj] == 0) {
                        dgecp(1, nx[ii + 1], qp->BAbt + ii, jj, 0, qp_red->BAbt + ii, idx0, 0);
                        idx0++;
                    }
                }
                // b
                dgemv_t(nu[ii] + nx[ii], nx[ii + 1], 1.0, qp->BAbt + ii, 0, 0, work->tmp_nuxM + 0, 0, 1.0, qp->b + ii, 0, qp_red->b + ii, 0);
            }
            // RSQ
            idx0 = 0;
            for (jj = 0; jj < nu[ii] + nx[ii]; jj++) {
                if (work->e_imask_ux[jj] == 0) {
                    //					idx1 = 0;
                    idx1 = idx0;
                    //					for(kk=0; kk<nu[ii]+nx[ii]; kk++)
                    for (kk = jj; kk < nu[ii] + nx[ii]; kk++) {
                        if (work->e_imask_ux[kk] == 0) {
                            MATEL(qp_red->RSQrq + ii, idx1, idx0) = MATEL(qp->RSQrq + ii, kk, jj);

                            idx1++;
                        }
                    }
                    idx0++;
                }
            }
            // rq
            dsymv_l(nu[ii] + nx[ii], 1.0, qp->RSQrq + ii, 0, 0, work->tmp_nuxM + 0, 0, 1.0, qp->rqz + ii, 0, work->tmp_nuxM + 1, 0);
            idx0 = 0;
            for (jj = 0; jj < nu[ii] + nx[ii]; jj++) {
                if (work->e_imask_ux[jj] == 0) {
                    VECEL(qp_red->rqz + ii, idx0) = VECEL(work->tmp_nuxM + 1, jj);
                    idx0++;
                }
            }
            // d d_mask m idxb idxs_rev
            idx0 = 0;
            for (jj = 0; jj < nb[ii]; jj++) {
                if (work->e_imask_d[jj] == 0) {
                    VECEL(qp_red->d + ii, idx0) = VECEL(qp->d + ii, jj);
                    VECEL(qp_red->d + ii, nb_red[ii] + ng_red[ii] + idx0) = VECEL(qp->d + ii, nb[ii] + ng[ii] + jj);
                    VECEL(qp_red->d_mask + ii, idx0) = VECEL(qp->d_mask + ii, jj);
                    VECEL(qp_red->d_mask + ii, nb_red[ii] + ng_red[ii] + idx0) = VECEL(qp->d_mask + ii, nb[ii] + ng[ii] + jj);
                    VECEL(qp_red->m + ii, idx0) = VECEL(qp->m + ii, jj);
                    VECEL(qp_red->m + ii, nb_red[ii] + ng_red[ii] + idx0) = VECEL(qp->m + ii, nb[ii] + ng[ii] + jj);
                    qp_red->idxb[ii][idx0] = qp->idxb[ii][jj];
                    qp_red->idxs_rev[ii][idx0] = qp->idxs_rev[ii][jj];  // keep softed inequality constr with same slack
                    idx0++;
                }
            }
            dveccp(ng[ii], qp->d + ii, nb[ii], qp_red->d + ii, nb_red[ii]);
            dveccp(ng[ii] + 2 * ns[ii], qp->d + ii, 2 * nb[ii] + ng[ii], qp_red->d + ii, 2 * nb_red[ii] + ng_red[ii]);
            dveccp(ng[ii], qp->d_mask + ii, nb[ii], qp_red->d_mask + ii, nb_red[ii]);
            dveccp(ng[ii] + 2 * ns[ii], qp->d_mask + ii, 2 * nb[ii] + ng[ii], qp_red->d_mask + ii, 2 * nb_red[ii] + ng_red[ii]);
            dveccp(ng[ii], qp->m + ii, nb[ii], qp_red->m + ii, nb_red[ii]);
            dveccp(ng[ii] + 2 * ns[ii], qp->m + ii, 2 * nb[ii] + ng[ii], qp_red->m + ii, 2 * nb_red[ii] + ng_red[ii]);
            dgemv_t(nu[ii] + nx[ii], ng[ii], 1.0, qp->DCt + ii, 0, 0, work->tmp_nuxM + 0, 0, 0.0, work->tmp_nbgM, 0, work->tmp_nbgM, 0);
            daxpy(ng[ii], -1.0, work->tmp_nbgM, 0, qp->d + ii, nb[ii], qp_red->d + ii, nb_red[ii]);
            daxpy(ng[ii], 1.0, work->tmp_nbgM, 0, qp->d + ii, 2 * nb[ii] + ng[ii], qp_red->d + ii, 2 * nb_red[ii] + ng_red[ii]);
            // DCt
            idx0 = 0;
            for (jj = 0; jj < nu[ii] + nx[ii]; jj++) {
                if (work->e_imask_ux[jj] == 0) {
                    dgecp(1, ng[ii], qp->DCt + ii, jj, 0, qp_red->DCt + ii, idx0, 0);
                    idx0++;
                }
            }
            // soft constraints
            for (jj = 0; jj < ng[ii]; jj++)
                qp_red->idxs_rev[ii][nb_red[ii] + jj] = qp->idxs_rev[ii][nb[ii] + jj];  // keep softed inequality constr with same slack
            dveccp(2 * ns[ii], qp->Z + ii, 0, qp_red->Z + ii, 0);
            dveccp(2 * ns[ii], qp->rqz + ii, nu[ii] + nx[ii], qp_red->rqz + ii, nu_red[ii] + nx_red[ii]);
            qp_red->diag_H_flag[ii] = qp->diag_H_flag[ii];
            // TODO idxe !!!!!!!!!!!!!!!
        } else  // copy everything
        {
            // copy vectors which are contiguous in the QP (e.g. to alias to res)
            if (ii < N) {
                dveccp(nx[ii + 1], qp->b + ii, 0, qp_red->b + ii, 0);
            }
            dveccp(nu[ii] + nx[ii] + 2 * ns[ii], qp->rqz + ii, 0, qp_red->rqz + ii, 0);
            dveccp(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], qp->d + ii, 0, qp_red->d + ii, 0);
            dveccp(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], qp->d_mask + ii, 0, qp_red->d_mask + ii, 0);
            dveccp(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], qp->m + ii, 0, qp_red->m + ii, 0);
            qp_red->diag_H_flag[ii] = qp->diag_H_flag[ii];
            if (arg->alias_unchanged) {
                if (ii < N) {
                    qp_red->BAbt[ii] = qp->BAbt[ii];
                }
                qp_red->RSQrq[ii] = qp->RSQrq[ii];
                qp_red->Z[ii] = qp->Z[ii];
                qp_red->idxb[ii] = qp->idxb[ii];
                qp_red->DCt[ii] = qp->DCt[ii];
                qp_red->idxs_rev[ii] = qp->idxs_rev[ii];
                qp_red->idxe[ii] = qp->idxe[ii];
            } else {
                if (ii < N) {
                    dgecp(nu[ii] + nx[ii] + 1, nx[ii + 1], qp->BAbt + ii, 0, 0, qp_red->BAbt + ii, 0, 0);
                }
                dgecp(nu[ii] + nx[ii] + 1, nu[ii] + nx[ii], qp->RSQrq + ii, 0, 0, qp_red->RSQrq + ii, 0, 0);
                dveccp(2 * ns[ii], qp->Z + ii, 0, qp_red->Z + ii, 0);
                for (jj = 0; jj < nb[ii]; jj++)
                    qp_red->idxb[ii][jj] = qp->idxb[ii][jj];
                dgecp(nu[ii] + nx[ii], ng[ii], qp->DCt + ii, 0, 0, qp_red->DCt + ii, 0, 0);
                for (jj = 0; jj < nb[ii] + ng[ii]; jj++)
                    qp_red->idxs_rev[ii][jj] = qp->idxs_rev[ii][jj];
                for (jj = 0; jj < nbue[ii] + nbxe[ii] + nge[ii]; jj++)
                    qp_red->idxe[ii][jj] = qp->idxe[ii][jj];
            }
        }
    }
}


void d_ocp_qp_reduce_eq_dof_lhs(struct d_ocp_qp* qp, struct d_ocp_qp* qp_red, struct d_ocp_qp_reduce_eq_dof_arg* arg, struct d_ocp_qp_reduce_eq_dof_ws* work) {

    int ii, jj, kk, idx0, idx1;

    struct d_ocp_qp_dim* dim = qp->dim;
    int N = dim->N;
    int* nx = dim->nx;
    int* nu = dim->nu;
    int* nb = dim->nb;
    int* nbx = dim->nbx;
    int* nbu = dim->nbu;
    int* ng = dim->ng;
    int* ns = dim->ns;
    int* nbue = dim->nbue;
    int* nbxe = dim->nbxe;
    int* nge = dim->nge;

    struct d_ocp_qp_dim* dim_red = qp_red->dim;
    int* nx_red = dim_red->nx;
    int* nu_red = dim_red->nu;
    int* nb_red = dim_red->nb;
    int* ng_red = dim_red->ng;
    int* ns_red = dim_red->ns;

    // TODO handle case of softed equalities !!!!!!!!!!!!!!!!

    int ne_thr;

    for (ii = 0; ii <= N; ii++) {
        if (ii == 0)
            ne_thr = nbue[ii] + nbxe[ii];
        else
            ne_thr = nbue[ii];
        if (ne_thr > 0)  // reduce inputs and/or states
        {
            dvecse(nu[ii] + nx[ii], 0.0, work->tmp_nuxM + 0, 0);
            for (jj = 0; jj < nu[ii] + nx[ii]; jj++)
                work->e_imask_ux[jj] = 0;
            for (jj = 0; jj < nbu[ii] + nbx[ii]; jj++)
                work->e_imask_d[jj] = 0;
            for (jj = 0; jj < ne_thr; jj++)  // set 1s for both inputs and states
            {
                VECEL(work->tmp_nuxM + 0, qp->idxb[ii][qp->idxe[ii][jj]]) = VECEL(qp->d + ii, qp->idxe[ii][jj]);
                work->e_imask_ux[qp->idxb[ii][qp->idxe[ii][jj]]] = 1;
                work->e_imask_d[qp->idxe[ii][jj]] = 1;
            }
            // TODO check first and last non-zero in e_mask and only multiply between them
            if (ii < N) {
                // BAt
                idx0 = 0;
                for (jj = 0; jj < nu[ii] + nx[ii]; jj++) {
                    if (work->e_imask_ux[jj] == 0) {
                        dgecp(1, nx[ii + 1], qp->BAbt + ii, jj, 0, qp_red->BAbt + ii, idx0, 0);
                        idx0++;
                    }
                }
            }
            // RSQ
            idx0 = 0;
            for (jj = 0; jj < nu[ii] + nx[ii]; jj++) {
                if (work->e_imask_ux[jj] == 0) {
                    //					idx1 = 0;
                    idx1 = idx0;
                    //					for(kk=0; kk<nu[ii]+nx[ii]; kk++)
                    for (kk = jj; kk < nu[ii] + nx[ii]; kk++) {
                        if (work->e_imask_ux[kk] == 0) {
                            MATEL(qp_red->RSQrq + ii, idx1, idx0) = MATEL(qp->RSQrq + ii, kk, jj);

                            idx1++;
                        }
                    }
                    idx0++;
                }
            }
            // idxb idxs_rev
            idx0 = 0;
            for (jj = 0; jj < nb[ii]; jj++) {
                if (work->e_imask_d[jj] == 0) {
                    qp_red->idxb[ii][idx0] = qp->idxb[ii][jj];
                    qp_red->idxs_rev[ii][idx0] = qp->idxs_rev[ii][jj];  // keep softed inequality constr with same slack
                    idx0++;
                }
            }
            // DCt
            idx0 = 0;
            for (jj = 0; jj < nu[ii] + nx[ii]; jj++) {
                if (work->e_imask_ux[jj] == 0) {
                    dgecp(1, ng[ii], qp->DCt + ii, jj, 0, qp_red->DCt + ii, idx0, 0);
                    idx0++;
                }
            }
            // soft constraints
            for (jj = 0; jj < ng[ii]; jj++)
                qp_red->idxs_rev[ii][nb_red[ii] + jj] = qp->idxs_rev[ii][nb[ii] + jj];  // keep softed inequality constr with same slack
            dveccp(2 * ns[ii], qp->Z + ii, 0, qp_red->Z + ii, 0);
            qp_red->diag_H_flag[ii] = qp->diag_H_flag[ii];
            // TODO idxe !!!!!!!!!!!!!!!
        } else  // copy everything
        {
            qp_red->diag_H_flag[ii] = qp->diag_H_flag[ii];
            if (arg->alias_unchanged) {
                if (ii < N) {
                    qp_red->BAbt[ii] = qp->BAbt[ii];
                }
                qp_red->RSQrq[ii] = qp->RSQrq[ii];
                qp_red->Z[ii] = qp->Z[ii];
                qp_red->idxb[ii] = qp->idxb[ii];
                qp_red->DCt[ii] = qp->DCt[ii];
                qp_red->idxs_rev[ii] = qp->idxs_rev[ii];
                qp_red->idxe[ii] = qp->idxe[ii];
            } else {
                if (ii < N) {
                    dgecp(nu[ii] + nx[ii] + 1, nx[ii + 1], qp->BAbt + ii, 0, 0, qp_red->BAbt + ii, 0, 0);
                }
                dgecp(nu[ii] + nx[ii] + 1, nu[ii] + nx[ii], qp->RSQrq + ii, 0, 0, qp_red->RSQrq + ii, 0, 0);
                dveccp(2 * ns[ii], qp->Z + ii, 0, qp_red->Z + ii, 0);
                for (jj = 0; jj < nb[ii]; jj++)
                    qp_red->idxb[ii][jj] = qp->idxb[ii][jj];
                dgecp(nu[ii] + nx[ii], ng[ii], qp->DCt + ii, 0, 0, qp_red->DCt + ii, 0, 0);
                for (jj = 0; jj < nb[ii] + ng[ii]; jj++)
                    qp_red->idxs_rev[ii][jj] = qp->idxs_rev[ii][jj];
                for (jj = 0; jj < nbue[ii] + nbxe[ii] + nge[ii]; jj++)
                    qp_red->idxe[ii][jj] = qp->idxe[ii][jj];
            }
        }
    }
}


void d_ocp_qp_reduce_eq_dof_rhs(struct d_ocp_qp* qp, struct d_ocp_qp* qp_red, struct d_ocp_qp_reduce_eq_dof_arg* arg, struct d_ocp_qp_reduce_eq_dof_ws* work) {

    int ii, jj, kk, idx0, idx1;

    struct d_ocp_qp_dim* dim = qp->dim;
    int N = dim->N;
    int* nx = dim->nx;
    int* nu = dim->nu;
    int* nb = dim->nb;
    int* nbx = dim->nbx;
    int* nbu = dim->nbu;
    int* ng = dim->ng;
    int* ns = dim->ns;
    int* nbue = dim->nbue;
    int* nbxe = dim->nbxe;
    int* nge = dim->nge;

    struct d_ocp_qp_dim* dim_red = qp_red->dim;
    int* nx_red = dim_red->nx;
    int* nu_red = dim_red->nu;
    int* nb_red = dim_red->nb;
    int* ng_red = dim_red->ng;
    int* ns_red = dim_red->ns;

    // TODO handle case of softed equalities !!!!!!!!!!!!!!!!

    int ne_thr;

    for (ii = 0; ii <= N; ii++) {
        if (ii == 0)
            ne_thr = nbue[ii] + nbxe[ii];
        else
            ne_thr = nbue[ii];
        if (ne_thr > 0)  // reduce inputs and/or states
        {
            dvecse(nu[ii] + nx[ii], 0.0, work->tmp_nuxM + 0, 0);
            for (jj = 0; jj < nu[ii] + nx[ii]; jj++)
                work->e_imask_ux[jj] = 0;
            for (jj = 0; jj < nbu[ii] + nbx[ii]; jj++)
                work->e_imask_d[jj] = 0;
            for (jj = 0; jj < ne_thr; jj++)  // set 1s for both inputs and states
            {
                VECEL(work->tmp_nuxM + 0, qp->idxb[ii][qp->idxe[ii][jj]]) = VECEL(qp->d + ii, qp->idxe[ii][jj]);
                work->e_imask_ux[qp->idxb[ii][qp->idxe[ii][jj]]] = 1;
                work->e_imask_d[qp->idxe[ii][jj]] = 1;
            }
            // TODO check first and last non-zero in e_mask and only multiply between them
            if (ii < N) {
                // b
                dgemv_t(nu[ii] + nx[ii], nx[ii + 1], 1.0, qp->BAbt + ii, 0, 0, work->tmp_nuxM + 0, 0, 1.0, qp->b + ii, 0, qp_red->b + ii, 0);
            }
            // rq
            dsymv_l(nu[ii] + nx[ii], 1.0, qp->RSQrq + ii, 0, 0, work->tmp_nuxM + 0, 0, 1.0, qp->rqz + ii, 0, work->tmp_nuxM + 1, 0);
            idx0 = 0;
            for (jj = 0; jj < nu[ii] + nx[ii]; jj++) {
                if (work->e_imask_ux[jj] == 0) {
                    VECEL(qp_red->rqz + ii, idx0) = VECEL(work->tmp_nuxM + 1, jj);
                    idx0++;
                }
            }
            // d d_mask m idxb idxs_rev
            idx0 = 0;
            for (jj = 0; jj < nb[ii]; jj++) {
                if (work->e_imask_d[jj] == 0) {
                    VECEL(qp_red->d + ii, idx0) = VECEL(qp->d + ii, jj);
                    VECEL(qp_red->d + ii, nb_red[ii] + ng_red[ii] + idx0) = VECEL(qp->d + ii, nb[ii] + ng[ii] + jj);
                    VECEL(qp_red->d_mask + ii, idx0) = VECEL(qp->d_mask + ii, jj);
                    VECEL(qp_red->d_mask + ii, nb_red[ii] + ng_red[ii] + idx0) = VECEL(qp->d_mask + ii, nb[ii] + ng[ii] + jj);
                    VECEL(qp_red->m + ii, idx0) = VECEL(qp->m + ii, jj);
                    VECEL(qp_red->m + ii, nb_red[ii] + ng_red[ii] + idx0) = VECEL(qp->m + ii, nb[ii] + ng[ii] + jj);
                    idx0++;
                }
            }
            dveccp(ng[ii], qp->d + ii, nb[ii], qp_red->d + ii, nb_red[ii]);
            dveccp(ng[ii] + 2 * ns[ii], qp->d + ii, 2 * nb[ii] + ng[ii], qp_red->d + ii, 2 * nb_red[ii] + ng_red[ii]);
            dveccp(ng[ii], qp->d_mask + ii, nb[ii], qp_red->d_mask + ii, nb_red[ii]);
            dveccp(ng[ii] + 2 * ns[ii], qp->d_mask + ii, 2 * nb[ii] + ng[ii], qp_red->d_mask + ii, 2 * nb_red[ii] + ng_red[ii]);
            dveccp(ng[ii], qp->m + ii, nb[ii], qp_red->m + ii, nb_red[ii]);
            dveccp(ng[ii] + 2 * ns[ii], qp->m + ii, 2 * nb[ii] + ng[ii], qp_red->m + ii, 2 * nb_red[ii] + ng_red[ii]);
            dgemv_t(nu[ii] + nx[ii], ng[ii], 1.0, qp->DCt + ii, 0, 0, work->tmp_nuxM + 0, 0, 0.0, work->tmp_nbgM, 0, work->tmp_nbgM, 0);
            daxpy(ng[ii], -1.0, work->tmp_nbgM, 0, qp->d + ii, nb[ii], qp_red->d + ii, nb_red[ii]);
            daxpy(ng[ii], 1.0, work->tmp_nbgM, 0, qp->d + ii, 2 * nb[ii] + ng[ii], qp_red->d + ii, 2 * nb_red[ii] + ng_red[ii]);
            // soft constraints
            dveccp(2 * ns[ii], qp->rqz + ii, nu[ii] + nx[ii], qp_red->rqz + ii, nu_red[ii] + nx_red[ii]);
            // TODO idxe !!!!!!!!!!!!!!!
        } else  // copy everything
        {
            // copy vectors which are contiguous in the QP (e.g. to alias to res)
            if (ii < N) {
                dveccp(nx[ii + 1], qp->b + ii, 0, qp_red->b + ii, 0);
            }
            dveccp(nu[ii] + nx[ii] + 2 * ns[ii], qp->rqz + ii, 0, qp_red->rqz + ii, 0);
            dveccp(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], qp->d + ii, 0, qp_red->d + ii, 0);
            dveccp(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], qp->d_mask + ii, 0, qp_red->d_mask + ii, 0);
            dveccp(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], qp->m + ii, 0, qp_red->m + ii, 0);
        }
    }
}


void d_ocp_qp_restore_eq_dof(struct d_ocp_qp* qp, struct d_ocp_qp_sol* qp_sol_red, struct d_ocp_qp_sol* qp_sol, struct d_ocp_qp_reduce_eq_dof_arg* arg, struct d_ocp_qp_reduce_eq_dof_ws* work) {

    int ii, jj, idx0;

    struct d_ocp_qp_dim* dim = qp->dim;
    int N = dim->N;
    int* nx = dim->nx;
    int* nu = dim->nu;
    int* nb = dim->nb;
    int* nbx = dim->nbx;
    int* nbu = dim->nbu;
    int* ng = dim->ng;
    int* ns = dim->ns;
    int* nbue = dim->nbue;
    int* nbxe = dim->nbxe;
    int* nge = dim->nge;

    struct d_ocp_qp_dim* dim_red = qp_sol_red->dim;
    int* nx_red = dim_red->nx;
    int* nu_red = dim_red->nu;
    int* nb_red = dim_red->nb;
    int* ng_red = dim_red->ng;

    int ne_thr;

    double tmp;

    for (ii = 0; ii <= N; ii++) {
        if (ii == 0)
            ne_thr = nbue[ii] + nbxe[ii];
        else
            ne_thr = nbue[ii];
        if (ne_thr > 0)  // restore inputs and/or states
        {
            dvecse(nu[ii] + nx[ii], 0.0, work->tmp_nuxM + 0, 0);
            for (jj = 0; jj < nu[ii] + nx[ii]; jj++)
                work->e_imask_ux[jj] = 0;
            for (jj = 0; jj < nb[ii]; jj++)
                work->e_imask_d[jj] = 0;
            for (jj = 0; jj < ne_thr; jj++)  // set 1s for both inputs and states
            {
                work->e_imask_ux[qp->idxb[ii][qp->idxe[ii][jj]]] = 1;
                work->e_imask_d[qp->idxe[ii][jj]] = 1;
            }
            // ux
            if (arg->comp_prim_sol) {
                idx0 = 0;
                for (jj = 0; jj < nu[ii] + nx[ii]; jj++) {
                    if (work->e_imask_ux[jj] == 0) {
                        VECEL(qp_sol->ux + ii, jj) = VECEL(qp_sol_red->ux + ii, idx0);
                        idx0++;
                    }
                }
                for (jj = 0; jj < ne_thr; jj++) {
                    VECEL(qp_sol->ux + ii, qp->idxb[ii][qp->idxe[ii][jj]]) = VECEL(qp->d + ii, qp->idxe[ii][jj]);
                }
                // TODO update based on choices on reduce !!!!!!!!!!!!!
                dveccp(2 * ns[ii], qp_sol_red->ux + ii, nu_red[ii] + nx_red[ii], qp_sol->ux + ii, nu[ii] + nx[ii]);
            }
            if (arg->comp_dual_sol_eq) {
                // pi
                if (ii < N)
                    dveccp(nx[ii + 1], qp_sol_red->pi + ii, 0, qp_sol->pi + ii, 0);
            }
            if (arg->comp_dual_sol_ineq) {
                // lam t
                idx0 = 0;
                for (jj = 0; jj < nb[ii]; jj++) {
                    if (work->e_imask_d[jj] == 0) {
                        VECEL(qp_sol->lam + ii, jj) = VECEL(qp_sol_red->lam + ii, idx0);
                        VECEL(qp_sol->lam + ii, nb[ii] + ng[ii] + jj) = VECEL(qp_sol_red->lam + ii, nb_red[ii] + ng_red[ii] + idx0);
                        VECEL(qp_sol->t + ii, jj) = VECEL(qp_sol_red->t + ii, idx0);
                        VECEL(qp_sol->t + ii, nb[ii] + ng[ii] + jj) = VECEL(qp_sol_red->t + ii, nb_red[ii] + ng_red[ii] + idx0);
                        idx0++;
                    } else {
                        // lam
                        // t
                        VECEL(qp_sol->lam + ii, jj) = arg->lam_min;
                        VECEL(qp_sol->lam + ii, nb[ii] + ng[ii] + jj) = arg->lam_min;
                        VECEL(qp_sol->t + ii, jj) = arg->t_min;
                        VECEL(qp_sol->t + ii, nb[ii] + ng[ii] + jj) = arg->t_min;
                    }
                }
                // TODO update based on choices on reduce !!!!!!!!!!!!!
                dveccp(ng[ii], qp_sol_red->lam + ii, nb_red[ii], qp_sol->lam + ii, nb[ii]);
                dveccp(ng[ii] + 2 * ns[ii], qp_sol_red->lam + ii, 2 * nb_red[ii] + ng_red[ii], qp_sol->lam + ii, 2 * nb[ii] + ng[ii]);
                dveccp(ng[ii], qp_sol_red->t + ii, nb_red[ii], qp_sol->t + ii, nb[ii]);
                dveccp(ng[ii] + 2 * ns[ii], qp_sol_red->t + ii, 2 * nb_red[ii] + ng_red[ii], qp_sol->t + ii, 2 * nb[ii] + ng[ii]);
                // update lam_lb for removed eq, keep lam_ub to zero
                dveccp(nu[ii] + nx[ii], qp->rqz + ii, 0, work->tmp_nuxM, 0);
                daxpy(nb[ii] + ng[ii], -1.0, qp_sol->lam + ii, 0, qp_sol->lam + ii, nb[ii] + ng[ii], work->tmp_nbgM, 0);
                dvecad_sp(nb[ii], 1.0, work->tmp_nbgM, 0, qp->idxb[ii], work->tmp_nuxM, 0);
                dsymv_l(nu[ii] + nx[ii], 1.0, qp->RSQrq + ii, 0, 0, qp_sol->ux + ii, 0, 1.0, work->tmp_nuxM, 0, work->tmp_nuxM, 0);
                if (ii < N)
                    dgemv_n(nu[ii] + nx[ii], nx[ii + 1], 1.0, qp->BAbt + ii, 0, 0, qp_sol_red->pi + ii, 0, 1.0, work->tmp_nuxM, 0, work->tmp_nuxM, 0);
                dgemv_n(nu[ii] + nx[ii], ng[ii], 1.0, qp->DCt + ii, 0, 0, work->tmp_nbgM, nb[ii], 1.0, work->tmp_nuxM, 0, work->tmp_nuxM, 0);
                for (jj = 0; jj < nb[ii]; jj++) {
                    if (work->e_imask_d[jj] != 0) {
                        tmp = VECEL(work->tmp_nuxM, qp->idxb[ii][jj]);
                        if (tmp >= 0)
                            VECEL(qp_sol->lam + ii, jj) = tmp;
                        else
                            VECEL(qp_sol->lam + ii, nb[ii] + ng[ii] + jj) = -tmp;
                    }
                }
            }
        } else  // copy
        {
            if (arg->comp_prim_sol) {
                // ux
                dveccp(nu[ii] + nx[ii] + 2 * ns[ii], qp_sol_red->ux + ii, 0, qp_sol->ux + ii, 0);
            }
            if (arg->comp_dual_sol_eq) {
                // pi
                if (ii < N)
                    dveccp(nx[ii + 1], qp_sol_red->pi + ii, 0, qp_sol->pi + ii, 0);
            }
            if (arg->comp_dual_sol_ineq) {
                // lam
                dveccp(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], qp_sol_red->lam + ii, 0, qp_sol->lam + ii, 0);
                // t
                dveccp(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], qp_sol_red->t + ii, 0, qp_sol->t + ii, 0);
            }
        }
    }
}