#include <stdio.h>
#include <stdlib.h>

#include "hpipm/blas.h"
#include "hpipm/common.h"
#include "hpipm/dense/d_dense_qcqp.h"
#include "hpipm/dense/d_dense_qcqp_dim.h"
#include "hpipm/dense/d_dense_qcqp_res.h"
#include "hpipm/dense/d_dense_qcqp_sol.h"


void d_dense_qcqp_dim_print(struct d_dense_qcqp_dim* qp_dim) {
    int ii;

    int nv = qp_dim->nv;
    int ne = qp_dim->ne;
    int nb = qp_dim->nb;
    int ng = qp_dim->ng;
    int nq = qp_dim->nq;
    int nsb = qp_dim->nsb;
    int nsg = qp_dim->nsg;
    int nsq = qp_dim->nsq;
    int ns = qp_dim->ns;

    printf("nv = %d\n\n", nv);
    printf("ne = %d\n\n", ne);
    printf("nb = %d\n\n", nb);
    printf("ng = %d\n\n", ng);
    printf("nq = %d\n\n", nq);
    printf("nsb = %d\n\n", nsb);
    printf("nsg = %d\n\n", nsg);
    printf("nsq = %d\n\n", nsq);
    printf("ns = %d\n\n", ns);
}


void d_dense_qcqp_print(struct d_dense_qcqp_dim* qp_dim, struct d_dense_qcqp* qp) {
    int ii;

    int nv = qp_dim->nv;
    int ne = qp_dim->ne;
    int nb = qp_dim->nb;
    int ng = qp_dim->ng;
    int nq = qp_dim->nq;
    int ns = qp_dim->ns;

    printf("H = \n");
    print_mat(nv, nv, qp->Hv, 0, 0);

    printf("A = \n");
    print_mat(ne, nv, qp->A, 0, 0);

    printf("Ct = \n");
    print_mat(nv, ng, qp->Ct, 0, 0);

    printf("Hq = \n");
    for (ii = 0; ii < nq; ii++)
        print_mat(nv, nv, qp->Hq + ii, 0, 0);

    printf("Hq_nzero = \n");
    int_print_mat(1, nq, qp->Hq_nzero, 1);

    printf("idxb = \n");
    int_print_mat(1, nb, qp->idxb, 1);

    printf("gz = \n");
    print_tran_vec(nv + 2 * ns, qp->gz, 0);

    printf("b = \n");
    print_tran_vec(ne, qp->b, 0);

    printf("d = \n");
    print_tran_vec(2 * nb + 2 * ng + 2 * nq + 2 * ns, qp->d, 0);

    printf("d_mask = \n");
    print_tran_vec(2 * nb + 2 * ng + 2 * nq + 2 * ns, qp->d_mask, 0);

    printf("m = \n");
    print_tran_vec(2 * nb + 2 * ng + 2 * nq + 2 * ns, qp->m, 0);

    printf("Z = \n");
    print_tran_vec(2 * ns, qp->Z, 0);

    printf("gq = \n");
    for (ii = 0; ii < nq; ii++)
        print_tran_mat(nv, 1, qp->Ct, 0, ng + ii);

    printf("idxs_rev = \n");
    int_print_mat(1, nb + ng + nq, qp->idxs_rev, 1);
}


void d_dense_qcqp_sol_print(struct d_dense_qcqp_dim* qp_dim, struct d_dense_qcqp_sol* qp_sol) {
    int ii;

    int nv = qp_dim->nv;
    int ne = qp_dim->ne;
    int nb = qp_dim->nb;
    int ng = qp_dim->ng;
    int nq = qp_dim->nq;
    int ns = qp_dim->ns;

    printf("v = \n");
    print_tran_vec(nv + 2 * ns, qp_sol->v, 0);

    printf("pi = \n");
    print_tran_vec(ne, qp_sol->pi, 0);

    printf("lam = \n");
    print_tran_vec(2 * nb + 2 * ng + 2 * nq + 2 * ns, qp_sol->lam, 0);

    printf("t = \n");
    print_tran_vec(2 * nb + 2 * ng + 2 * nq + 2 * ns, qp_sol->t, 0);
}


void d_dense_qcqp_res_print(struct d_dense_qcqp_dim* qp_dim, struct d_dense_qcqp_res* qp_res) {
    int ii;

    int nv = qp_dim->nv;
    int ne = qp_dim->ne;
    int nb = qp_dim->nb;
    int ng = qp_dim->ng;
    int nq = qp_dim->nq;
    int ns = qp_dim->ns;

    printf("res_g = \n");
    print_tran_vec(nv + 2 * ns, qp_res->res_g, 0);

    printf("res_b = \n");
    print_tran_vec(ne, qp_res->res_b, 0);

    printf("res_d = \n");
    print_tran_vec(2 * nb + 2 * ng + 2 * nq + 2 * ns, qp_res->res_d, 0);

    printf("res_m = \n");
    print_tran_vec(2 * nb + 2 * ng + 2 * nq + 2 * ns, qp_res->res_m, 0);
}
