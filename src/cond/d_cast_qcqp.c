#include <stdio.h>
#include <stdlib.h>

#include "tinyhpipm/blas.h"
#include "tinyhpipm/dense/d_dense_qcqp.h"
#include "tinyhpipm/dense/d_dense_qcqp_dim.h"
#include "tinyhpipm/dense/d_dense_qcqp_sol.h"
#include "tinyhpipm/ocp/d_ocp_qcqp.h"
#include "tinyhpipm/ocp/d_ocp_qcqp_dim.h"
#include "tinyhpipm/ocp/d_ocp_qcqp_sol.h"

void d_cast_qcqp_compute_dim(struct d_ocp_qcqp_dim* ocp_dim, struct d_dense_qcqp_dim* dense_dim) {

    int N = ocp_dim->N;
    int* nx = ocp_dim->nx;
    int* nu = ocp_dim->nu;
    int* nbx = ocp_dim->nbx;
    int* nbu = ocp_dim->nbu;
    int* ng = ocp_dim->ng;
    int* nq = ocp_dim->nq;
    int* ns = ocp_dim->ns;
    int* nsbx = ocp_dim->nsbx;
    int* nsbu = ocp_dim->nsbu;
    int* nsg = ocp_dim->nsg;
    int* nsq = ocp_dim->nsq;


    int nvc = 0;
    int nec = 0;
    int nbc = 0;
    int ngc = 0;
    int nqc = 0;
    int nsc = 0;
    int nsbc = 0;
    int nsgc = 0;
    int nsqc = 0;

    // first stage
    nvc += nx[0] + nu[0];
    nec = 0;
    nbc += nbx[0] + nbu[0];
    ngc += ng[0];
    nqc += nq[0];
    nsc += ns[0];
    nsbc += nsbx[0] + nsbu[0];
    nsgc += nsg[0];
    nsqc += nsq[0];
    // remaining stages
    for (int ii = 1; ii <= N; ii++) {
        nvc += nx[ii] + nu[ii];
        nec += nx[ii];
        nbc += nbx[ii] + nbu[ii];
        ngc += ng[ii];
        nqc += nq[ii];
        nsc += ns[ii];
        nsbc += nsbx[ii] + nsbu[ii];
        nsgc += nsg[ii];
        nsqc += nsq[ii];
    }

    // XXX must use setters to correctly set qp ones too !
    d_dense_qcqp_dim_set_nv(nvc, dense_dim);
    d_dense_qcqp_dim_set_ne(nec, dense_dim);
    d_dense_qcqp_dim_set_nb(nbc, dense_dim);
    d_dense_qcqp_dim_set_ng(ngc, dense_dim);
    d_dense_qcqp_dim_set_nq(nqc, dense_dim);
    d_dense_qcqp_dim_set_ns(nsc, dense_dim);
    d_dense_qcqp_dim_set_nsb(nsbc, dense_dim);
    d_dense_qcqp_dim_set_nsg(nsgc, dense_dim);
    d_dense_qcqp_dim_set_nsq(nsqc, dense_dim);
}


void d_cast_qcqp_cond(struct d_ocp_qcqp* ocp_qp, struct d_dense_qcqp* dense_qp) {

    int ii, jj;

    int N = ocp_qp->dim->N;
    int* nu = ocp_qp->dim->nu;
    int* nx = ocp_qp->dim->nx;
    int* nb = ocp_qp->dim->nb;
    int* ng = ocp_qp->dim->ng;
    int* nq = ocp_qp->dim->nq;

    //	int nvc = dense_qp->dim->nv;
    //	int nec = dense_qp->dim->ne;
    int nbc = dense_qp->dim->nb;
    int ngc = dense_qp->dim->ng;
    int nqc = dense_qp->dim->nq;

    int idxc, idxr, idxq;

    // cost
    idxr = 0;
    idxc = 0;
    for (ii = 0; ii <= N; ii++) {
        dgecp(nu[ii] + nx[ii], nu[ii] + nx[ii], ocp_qp->RSQrq + ii, 0, 0, dense_qp->Hv, idxr, idxc);
        dveccp(nu[ii] + nx[ii], ocp_qp->rqz + ii, 0, dense_qp->gz, idxr);
        idxr += nu[ii] + nx[ii];
        idxc += nu[ii] + nx[ii];
    }

    // dynamics
    idxr = 0;
    idxc = 0;
    for (ii = 0; ii < N; ii++) {
        dgetr(nu[ii] + nx[ii], nx[ii + 1], ocp_qp->BAbt + ii, 0, 0, dense_qp->A, idxr, idxc);
        ddiare(nx[ii + 1], -1.0, dense_qp->A, idxr, idxc + nu[ii] + nx[ii] + nu[ii + 1]);
        dveccp(nx[ii + 1], ocp_qp->b + ii, 0, dense_qp->b, idxr);
        idxr += nx[ii + 1];
        idxc += nu[ii] + nx[ii];
    }

    // box constraints
    idxr = 0;
    idxc = 0;
    for (ii = 0; ii <= N; ii++) {
        dveccp(nb[ii], ocp_qp->d + ii, 0, dense_qp->d, idxr);
        dveccp(nb[ii], ocp_qp->d + ii, nb[ii] + ng[ii] + nq[ii], dense_qp->d, idxr + nbc + ngc + nqc);
        for (jj = 0; jj < nb[ii]; jj++)
            dense_qp->idxb[idxr + jj] = ocp_qp->idxb[ii][jj] + idxc;
        idxr += nb[ii];
        idxc += nu[ii] + nx[ii];
    }

    // general constraints
    idxr = 0;
    idxc = 0;
    for (ii = 0; ii <= N; ii++) {
        dveccp(ng[ii], ocp_qp->d + ii, nb[ii], dense_qp->d, idxr + nbc);
        dveccp(ng[ii], ocp_qp->d + ii, 2 * nb[ii] + ng[ii] + nq[ii], dense_qp->d, idxr + 2 * nbc + ngc + nqc);
        dgecp(nu[ii] + nx[ii], ng[ii], ocp_qp->DCt + ii, 0, 0, dense_qp->Ct, idxc, idxr);
        idxr += ng[ii];
        idxc += nu[ii] + nx[ii];
    }

    // quadratic constraints
    idxr = 0;
    idxc = 0;
    idxq = 0;
    for (ii = 0; ii <= N; ii++) {
        for (jj = 0; jj < nq[ii]; jj++) {
            dgecp(nu[ii] + nx[ii], nu[ii] + nx[ii], ocp_qp->Hq[ii] + jj, 0, 0, dense_qp->Hq + idxq, idxr, idxc);
            idxq++;
        }
        idxr += nu[ii] + nx[ii];
        idxc += nu[ii] + nx[ii];
    }
    idxr = 0;
    idxc = 0;
    for (ii = 0; ii <= N; ii++) {
        dgecp(nu[ii] + nx[ii], nq[ii], ocp_qp->DCt + ii, 0, ng[ii], dense_qp->Ct, idxr, ngc + idxc);
        idxr += nu[ii] + nx[ii];
        idxc += nq[ii];
    }
    idxr = 0;
    idxc = 0;
    for (ii = 0; ii <= N; ii++) {
        //		dveccp(nb[ii], ocp_qp->d+ii, 0, dense_qp->d, idxr);
        dveccp(nq[ii], ocp_qp->d + ii, 2 * nb[ii] + 2 * ng[ii] + nq[ii], dense_qp->d, idxr + 2 * nbc + 2 * ngc + nqc);
        idxr += nq[ii];
    }

    // TODO soft constraints
}
