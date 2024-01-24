#include <stdio.h>
#include <stdlib.h>

#include "hpipm/blas.h"
#include "hpipm/common.h"
#include "hpipm/ocp/d_ocp_qcqp_dim.h"
#include "hpipm/ocp/d_ocp_qcqp_ipm.h"
#include "hpipm/ocp/d_ocp_qcqp_sol.h"
#include "hpipm/ocp/d_ocp_qp.h"

void d_ocp_qcqp_dim_print(struct d_ocp_qcqp_dim* qp_dim) {
    int ii;

    int N = qp_dim->N;
    int* nx = qp_dim->nx;
    int* nu = qp_dim->nu;
    int* nbx = qp_dim->nbx;
    int* nbu = qp_dim->nbu;
    int* ng = qp_dim->ng;
    int* nq = qp_dim->nq;
    int* ns = qp_dim->ns;
    int* nsbx = qp_dim->nsbx;
    int* nsbu = qp_dim->nsbu;
    int* nsg = qp_dim->nsg;
    int* nsq = qp_dim->nsq;
    int* nbxe = qp_dim->nbxe;
    int* nbue = qp_dim->nbue;
    int* nge = qp_dim->nge;
    int* nqe = qp_dim->nqe;

    printf("N = %d\n\n", N);

    printf("nx =\n");
    for (ii = 0; ii <= N; ii++)
        printf("\t%d", nx[ii]);
    printf("\n\n");

    printf("nu =\n");
    for (ii = 0; ii <= N; ii++)
        printf("\t%d", nu[ii]);
    printf("\n\n");

    printf("nbx =\n");
    for (ii = 0; ii <= N; ii++)
        printf("\t%d", nbx[ii]);
    printf("\n\n");

    printf("nbu =\n");
    for (ii = 0; ii <= N; ii++)
        printf("\t%d", nbu[ii]);
    printf("\n\n");

    printf("ng =\n");
    for (ii = 0; ii <= N; ii++)
        printf("\t%d", ng[ii]);
    printf("\n\n");

    printf("nq =\n");
    for (ii = 0; ii <= N; ii++)
        printf("\t%d", nq[ii]);
    printf("\n\n");

    printf("ns =\n");
    for (ii = 0; ii <= N; ii++)
        printf("\t%d", ns[ii]);
    printf("\n\n");

    printf("nsbx =\n");
    for (ii = 0; ii <= N; ii++)
        printf("\t%d", nsbx[ii]);
    printf("\n\n");

    printf("nsbu =\n");
    for (ii = 0; ii <= N; ii++)
        printf("\t%d", nsbu[ii]);
    printf("\n\n");

    printf("nsg =\n");
    for (ii = 0; ii <= N; ii++)
        printf("\t%d", nsg[ii]);
    printf("\n\n");

    printf("nsq =\n");
    for (ii = 0; ii <= N; ii++)
        printf("\t%d", nsq[ii]);
    printf("\n\n");

    printf("nbxe =\n");
    for (ii = 0; ii <= N; ii++)
        printf("\t%d", nbxe[ii]);
    printf("\n\n");

    printf("nbue =\n");
    for (ii = 0; ii <= N; ii++)
        printf("\t%d", nbue[ii]);
    printf("\n\n");

    printf("nge =\n");
    for (ii = 0; ii <= N; ii++)
        printf("\t%d", nge[ii]);
    printf("\n\n");

    printf("nqe =\n");
    for (ii = 0; ii <= N; ii++)
        printf("\t%d", nqe[ii]);
    printf("\n\n");
}


void d_ocp_qcqp_dim_codegen(char* file_name, char* mode, struct d_ocp_qcqp_dim* qp_dim) {
    int ii;

    FILE* file = fopen(file_name, mode);

    int N = qp_dim->N;
    int* nx = qp_dim->nx;
    int* nu = qp_dim->nu;
    int* nbx = qp_dim->nbx;
    int* nbu = qp_dim->nbu;
    int* ng = qp_dim->ng;
    int* nq = qp_dim->nq;
    int* ns = qp_dim->ns;
    int* nsbx = qp_dim->nsbx;
    int* nsbu = qp_dim->nsbu;
    int* nsg = qp_dim->nsg;
    int* nsq = qp_dim->nsq;
    int* nbxe = qp_dim->nbxe;
    int* nbue = qp_dim->nbue;
    int* nge = qp_dim->nge;
    int* nqe = qp_dim->nqe;

    fprintf(file, "/***************\n* dim\n***************/\n");

    // N
    fprintf(file, "/* N */\n");
    fprintf(file, "int N = %d;\n", N);
    // nx
    fprintf(file, "/* nx */\n");
    fprintf(file, "static int nnx[] = {");
    for (ii = 0; ii <= N; ii++)
        fprintf(file, "%d, ", nx[ii]);
    fprintf(file, "};\n");
    fprintf(file, "int *nx = nnx;\n");
    // nu
    fprintf(file, "/* nu */\n");
    fprintf(file, "static int nnu[] = {");
    for (ii = 0; ii <= N; ii++)
        fprintf(file, "%d, ", nu[ii]);
    fprintf(file, "};\n");
    fprintf(file, "int *nu = nnu;\n");
    // nbx
    fprintf(file, "/* nbx */\n");
    fprintf(file, "static int nnbx[] = {");
    for (ii = 0; ii <= N; ii++)
        fprintf(file, "%d, ", nbx[ii]);
    fprintf(file, "};\n");
    fprintf(file, "int *nbx = nnbx;\n");
    // nbu
    fprintf(file, "/* nbu */\n");
    fprintf(file, "static int nnbu[] = {");
    for (ii = 0; ii <= N; ii++)
        fprintf(file, "%d, ", nbu[ii]);
    fprintf(file, "};\n");
    fprintf(file, "int *nbu = nnbu;\n");
    // ng
    fprintf(file, "/* ng */\n");
    fprintf(file, "static int nng[] = {");
    for (ii = 0; ii <= N; ii++)
        fprintf(file, "%d, ", ng[ii]);
    fprintf(file, "};\n");
    fprintf(file, "int *ng = nng;\n");
    // nq
    fprintf(file, "/* nq */\n");
    fprintf(file, "static int nnq[] = {");
    for (ii = 0; ii <= N; ii++)
        fprintf(file, "%d, ", nq[ii]);
    fprintf(file, "};\n");
    fprintf(file, "int *nq = nnq;\n");
    // ns
    fprintf(file, "/* ns */\n");
    fprintf(file, "static int nns[] = {");
    for (ii = 0; ii <= N; ii++)
        fprintf(file, "%d, ", ns[ii]);
    fprintf(file, "};\n");
    fprintf(file, "int *ns = nns;\n");
    // nsbx
    fprintf(file, "/* nsbx */\n");
    fprintf(file, "static int nnsbx[] = {");
    for (ii = 0; ii <= N; ii++)
        fprintf(file, "%d, ", nsbx[ii]);
    fprintf(file, "};\n");
    fprintf(file, "int *nsbx = nnsbx;\n");
    // nsbu
    fprintf(file, "/* nsbu */\n");
    fprintf(file, "static int nnsbu[] = {");
    for (ii = 0; ii <= N; ii++)
        fprintf(file, "%d, ", nsbu[ii]);
    fprintf(file, "};\n");
    fprintf(file, "int *nsbu = nnsbu;\n");
    // nsg
    fprintf(file, "/* nsg */\n");
    fprintf(file, "static int nnsg[] = {");
    for (ii = 0; ii <= N; ii++)
        fprintf(file, "%d, ", nsg[ii]);
    fprintf(file, "};\n");
    fprintf(file, "int *nsg = nnsg;\n");
    // nsq
    fprintf(file, "/* nsq */\n");
    fprintf(file, "static int nnsq[] = {");
    for (ii = 0; ii <= N; ii++)
        fprintf(file, "%d, ", nsq[ii]);
    fprintf(file, "};\n");
    fprintf(file, "int *nsq = nnsq;\n");
    // nbxe
    fprintf(file, "/* nbxe */\n");
    fprintf(file, "static int nnbxe[] = {");
    for (ii = 0; ii <= N; ii++)
        fprintf(file, "%d, ", nbxe[ii]);
    fprintf(file, "};\n");
    fprintf(file, "int *nbxe = nnbxe;\n");
    // nbue
    fprintf(file, "/* nbue */\n");
    fprintf(file, "static int nnbue[] = {");
    for (ii = 0; ii <= N; ii++)
        fprintf(file, "%d, ", nbue[ii]);
    fprintf(file, "};\n");
    fprintf(file, "int *nbue = nnbue;\n");
    // nge
    fprintf(file, "/* nge */\n");
    fprintf(file, "static int nnge[] = {");
    for (ii = 0; ii <= N; ii++)
        fprintf(file, "%d, ", nge[ii]);
    fprintf(file, "};\n");
    fprintf(file, "int *nge = nnge;\n");
    // nqe
    fprintf(file, "/* nqe */\n");
    fprintf(file, "static int nnqe[] = {");
    for (ii = 0; ii <= N; ii++)
        fprintf(file, "%d, ", nqe[ii]);
    fprintf(file, "};\n");
    fprintf(file, "int *nqe = nnqe;\n");

    fclose(file);
}


void d_ocp_qcqp_print(struct d_ocp_qcqp_dim* dim, struct d_ocp_qcqp* qp) {
    int ii, jj;

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

    printf("BAt =\n");
    for (ii = 0; ii < N; ii++)
        print_mat(nu[ii] + nx[ii], nx[ii + 1], qp->BAbt + ii, 0, 0);

    printf("b =\n");
    for (ii = 0; ii < N; ii++)
        print_tran_vec(nx[ii + 1], qp->b + ii, 0);

    printf("RSQ =\n");
    for (ii = 0; ii <= N; ii++)
        print_mat(nu[ii] + nx[ii], nu[ii] + nx[ii], qp->RSQrq + ii, 0, 0);

    printf("Z =\n");
    for (ii = 0; ii <= N; ii++)
        print_tran_vec(2 * ns[ii], qp->Z + ii, 0);

    printf("rqz =\n");
    for (ii = 0; ii <= N; ii++)
        print_tran_vec(nu[ii] + nx[ii] + 2 * ns[ii], qp->rqz + ii, 0);

    printf("idxb = \n");
    for (ii = 0; ii <= N; ii++)
        int_print_mat(1, nb[ii], qp->idxb[ii], 1);

    printf("d =\n");
    for (ii = 0; ii <= N; ii++)
        print_tran_vec(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], qp->d + ii, 0);

    printf("d_mask =\n");
    for (ii = 0; ii <= N; ii++)
        print_tran_vec(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], qp->d_mask + ii, 0);

    printf("DCt =\n");
    for (ii = 0; ii <= N; ii++)
        print_mat(nu[ii] + nx[ii], ng[ii], qp->DCt + ii, 0, 0);

    printf("Hq =\n");
    for (ii = 0; ii <= N; ii++)
        if (nq[ii] == 0)
            printf("\n\n");
        else
            for (jj = 0; jj < nq[ii]; jj++)
                print_mat(nu[ii] + nx[ii], nu[ii] + nx[ii], &qp->Hq[ii][jj], 0, 0);

    printf("Hq_nzero = \n");
    for (ii = 0; ii <= N; ii++)
        int_print_mat(1, nq[ii], qp->Hq_nzero[ii], 1);

    printf("gq =\n");
    for (ii = 0; ii <= N; ii++)
        if (nq[ii] == 0)
            printf("\n\n");
        else
            for (jj = 0; jj < nq[ii]; jj++)
                print_tran_mat(nu[ii] + nx[ii], 1, qp->DCt + ii, 0, ng[ii] + jj);

    printf("idxs_rev = \n");
    for (ii = 0; ii <= N; ii++)
        int_print_mat(1, nb[ii] + ng[ii] + nq[ii], qp->idxs_rev[ii], 1);

    printf("m =\n");
    for (ii = 0; ii <= N; ii++)
        print_tran_vec(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], qp->m + ii, 0);

    printf("idxe = \n");
    for (ii = 0; ii <= N; ii++)
        int_print_mat(1, nbue[ii] + nbxe[ii] + nge[ii] + nqe[ii], qp->idxe[ii], 1);
}


void d_ocp_qcqp_codegen(char* file_name, char* mode, struct d_ocp_qcqp_dim* dim, struct d_ocp_qcqp* qp) {
    int nn, ii, jj, kk, idx_tmp;

    FILE* file = fopen(file_name, mode);

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

    fprintf(file, "/***************\n* qp\n***************/\n");

    // A
    fprintf(file, "/* A */\n");
    for (nn = 0; nn < N; nn++) {

        fprintf(file, "static double A%d[] = {", nn);


        for (ii = 0; ii < nx[nn]; ii++) {
            for (jj = 0; jj < nx[nn + 1]; jj++) {

                fprintf(file, "%18.15e, ", MATEL(qp->BAbt + nn, nu[nn] + ii, jj));
            }
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *AA[] = {");


    for (nn = 0; nn < N; nn++)
        fprintf(file, "A%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hA = AA;\n");


    // B
    fprintf(file, "/* B */\n");
    for (nn = 0; nn < N; nn++) {

        fprintf(file, "static double B%d[] = {", nn);


        for (ii = 0; ii < nu[nn]; ii++) {
            for (jj = 0; jj < nx[nn + 1]; jj++) {

                fprintf(file, "%18.15e, ", MATEL(qp->BAbt + nn, ii, jj));
            }
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *BB[] = {");


    for (nn = 0; nn < N; nn++)
        fprintf(file, "B%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hB = BB;\n");


    // b
    fprintf(file, "/* b */\n");
    for (nn = 0; nn < N; nn++) {

        fprintf(file, "static double b%d[] = {", nn);


        for (jj = 0; jj < nx[nn + 1]; jj++) {

            fprintf(file, "%18.15e, ", VECEL(qp->b + nn, jj));
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *bb[] = {");


    for (nn = 0; nn < N; nn++)
        fprintf(file, "b%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hb = bb;\n");


    // Q
    fprintf(file, "/* Q */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double Q%d[] = {", nn);


        for (jj = 0; jj < nx[nn]; jj++) {
            for (ii = 0; ii < nx[nn]; ii++) {

                fprintf(file, "%18.15e, ", MATEL(qp->RSQrq + nn, nu[nn] + ii, nu[nn] + jj));
            }
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *QQ[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "Q%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hQ = QQ;\n");


    // S
    fprintf(file, "/* S */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double S%d[] = {", nn);


        for (ii = 0; ii < nx[nn]; ii++) {
            for (jj = 0; jj < nu[nn]; jj++) {

                fprintf(file, "%18.15e, ", MATEL(qp->RSQrq + nn, nu[nn] + ii, jj));
            }
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *SS[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "S%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hS = SS;\n");


    // R
    fprintf(file, "/* R */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double R%d[] = {", nn);


        for (jj = 0; jj < nu[nn]; jj++) {
            for (ii = 0; ii < nu[nn]; ii++) {

                fprintf(file, "%18.15e, ", MATEL(qp->RSQrq + nn, ii, jj));
            }
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *RR[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "R%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hR = RR;\n");


    // r
    fprintf(file, "/* r */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double r%d[] = {", nn);


        for (jj = 0; jj < nu[nn]; jj++) {

            fprintf(file, "%18.15e, ", VECEL(qp->rqz + nn, jj));
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *rr[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "r%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hr = rr;\n");


    // q
    fprintf(file, "/* q */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double q%d[] = {", nn);


        for (jj = 0; jj < nx[nn]; jj++) {

            fprintf(file, "%18.15e, ", VECEL(qp->rqz + nn, nu[nn] + jj));
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *qq[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "q%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hq = qq;\n");


    // idxbu
    fprintf(file, "/* idxbu */\n");
    for (nn = 0; nn <= N; nn++) {
        fprintf(file, "static int idxbu%d[] = {", nn);
        for (jj = 0; jj < nb[nn]; jj++) {
            if (qp->idxb[nn][jj] < nu[nn]) {
                fprintf(file, "%d, ", qp->idxb[nn][jj]);
            }
        }
        fprintf(file, "};\n");
    }
    fprintf(file, "static int *iidxbu[] = {");
    for (nn = 0; nn <= N; nn++)
        fprintf(file, "idxbu%d, ", nn);
    fprintf(file, "};\n");
    fprintf(file, "int **hidxbu = iidxbu;\n");

    // lbu
    fprintf(file, "/* lbu */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double lbu%d[] = {", nn);


        for (jj = 0; jj < nb[nn]; jj++) {
            if (qp->idxb[nn][jj] < nu[nn]) {
                fprintf(file, "%18.15e, ", VECEL(qp->d + nn, jj));
            }
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *llbu[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "lbu%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hlbu = llbu;\n");


    // lbu_mask
    fprintf(file, "/* lbu_mask */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double lbu_mask%d[] = {", nn);


        for (jj = 0; jj < nb[nn]; jj++) {
            if (qp->idxb[nn][jj] < nu[nn]) {
                fprintf(file, "%18.15e, ", VECEL(qp->d_mask + nn, jj));
            }
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *llbu_mask[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "lbu_mask%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hlbu_mask = llbu_mask;\n");


    // ubu
    fprintf(file, "/* ubu */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double ubu%d[] = {", nn);


        for (jj = 0; jj < nb[nn]; jj++) {
            if (qp->idxb[nn][jj] < nu[nn]) {
                fprintf(file, "%18.15e, ", -VECEL(qp->d + nn, nb[nn] + ng[nn] + nq[nn] + jj));
            }
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *uubu[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "ubu%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hubu = uubu;\n");


    // ubu_mask
    fprintf(file, "/* ubu_mask */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double ubu_mask%d[] = {", nn);


        for (jj = 0; jj < nb[nn]; jj++) {
            if (qp->idxb[nn][jj] < nu[nn]) {
                fprintf(file, "%18.15e, ", VECEL(qp->d_mask + nn, nb[nn] + ng[nn] + nq[nn] + jj));
            }
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *uubu_mask[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "ubu_mask%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hubu_mask = uubu_mask;\n");


    // idxbx
    fprintf(file, "/* idxbx */\n");
    for (nn = 0; nn <= N; nn++) {
        fprintf(file, "static int idxbx%d[] = {", nn);
        for (jj = 0; jj < nb[nn]; jj++) {
            if (qp->idxb[nn][jj] >= nu[nn]) {
                fprintf(file, "%d, ", qp->idxb[nn][jj] - nu[nn]);
            }
        }
        fprintf(file, "};\n");
    }
    fprintf(file, "static int *iidxbx[] = {");
    for (nn = 0; nn <= N; nn++)
        fprintf(file, "idxbx%d, ", nn);
    fprintf(file, "};\n");
    fprintf(file, "int **hidxbx = iidxbx;\n");

    // lbx
    fprintf(file, "/* lbx */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double lbx%d[] = {", nn);


        for (jj = 0; jj < nb[nn]; jj++) {
            if (qp->idxb[nn][jj] >= nu[nn]) {
                fprintf(file, "%18.15e, ", VECEL(qp->d + nn, jj));
            }
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *llbx[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "lbx%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hlbx = llbx;\n");


    // lbx_mask
    fprintf(file, "/* lbx_mask */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double lbx_mask%d[] = {", nn);


        for (jj = 0; jj < nb[nn]; jj++) {
            if (qp->idxb[nn][jj] >= nu[nn]) {
                fprintf(file, "%18.15e, ", VECEL(qp->d_mask + nn, jj));
            }
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *llbx_mask[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "lbx_mask%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hlbx_mask = llbx_mask;\n");


    // ubx
    fprintf(file, "/* ubx */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double ubx%d[] = {", nn);


        for (jj = 0; jj < nb[nn]; jj++) {
            if (qp->idxb[nn][jj] >= nu[nn]) {
                fprintf(file, "%18.15e, ", -VECEL(qp->d + nn, nb[nn] + ng[nn] + nq[nn] + jj));
            }
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *uubx[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "ubx%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hubx = uubx;\n");


    // ubx_mask
    fprintf(file, "/* ubx_mask */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double ubx_mask%d[] = {", nn);


        for (jj = 0; jj < nb[nn]; jj++) {
            if (qp->idxb[nn][jj] >= nu[nn]) {
                fprintf(file, "%18.15e, ", VECEL(qp->d_mask + nn, nb[nn] + ng[nn] + nq[nn] + jj));
            }
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *uubx_mask[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "ubx_mask%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hubx_mask = uubx_mask;\n");


    // C
    fprintf(file, "/* C */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double C%d[] = {", nn);


        for (ii = 0; ii < nx[nn]; ii++) {
            for (jj = 0; jj < ng[nn]; jj++) {

                fprintf(file, "%18.15e, ", MATEL(qp->DCt + nn, nu[nn] + ii, jj));
            }
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *CC[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "C%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hC = CC;\n");


    // D
    fprintf(file, "/* D */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double D%d[] = {", nn);


        for (ii = 0; ii < nu[nn]; ii++) {
            for (jj = 0; jj < ng[nn]; jj++) {

                fprintf(file, "%18.15e, ", MATEL(qp->DCt + nn, ii, jj));
            }
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *DD[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "D%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hD = DD;\n");


    // lg
    fprintf(file, "/* lg */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double lg%d[] = {", nn);


        for (jj = 0; jj < ng[nn]; jj++) {
            fprintf(file, "%18.15e, ", VECEL(qp->d + nn, nb[nn] + jj));
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *llg[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "lg%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hlg = llg;\n");


    // lg_mask
    fprintf(file, "/* lg_mask */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double lg_mask%d[] = {", nn);


        for (jj = 0; jj < ng[nn]; jj++) {
            fprintf(file, "%18.15e, ", VECEL(qp->d_mask + nn, nb[nn] + jj));
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *llg_mask[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "lg_mask%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hlg_mask = llg_mask;\n");


    // ug
    fprintf(file, "/* ug */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double ug%d[] = {", nn);


        for (jj = 0; jj < ng[nn]; jj++) {
            fprintf(file, "%18.15e, ", -VECEL(qp->d + nn, 2 * nb[nn] + ng[nn] + nq[nn] + jj));
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *uug[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "ug%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hug = uug;\n");


    // ug_mask
    fprintf(file, "/* ug_mask */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double ug_mask%d[] = {", nn);


        // TODO: why no consideration of floating point precision here?
        for (jj = 0; jj < ng[nn]; jj++) {
            fprintf(file, "%18.15e, ", VECEL(qp->d_mask + nn, 2 * nb[nn] + ng[nn] + nq[nn] + jj));
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *uug_mask[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "ug_mask%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hug_mask = uug_mask;\n");


    // Qq
    fprintf(file, "/* Qq */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double Qq%d[] = {", nn);


        for (kk = 0; kk < nq[nn]; kk++) {
            for (jj = 0; jj < nx[nn]; jj++) {
                for (ii = 0; ii < nx[nn]; ii++) {

                    fprintf(file, "%18.15e, ", MATEL(&qp->Hq[nn][kk], nu[nn] + ii, nu[nn] + jj));
                }
            }
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *QQq[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "Qq%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hQq = QQq;\n");


    // Sq
    fprintf(file, "/* Sq */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double Sq%d[] = {", nn);


        for (kk = 0; kk < nq[nn]; kk++) {
            for (ii = 0; ii < nx[nn]; ii++) {
                for (jj = 0; jj < nu[nn]; jj++) {

                    fprintf(file, "%18.15e, ", MATEL(&qp->Hq[nn][kk], nu[nn] + ii, jj));
                }
            }
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *SSq[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "Sq%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hSq = SSq;\n");


    // Rq
    fprintf(file, "/* Rq */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double Rq%d[] = {", nn);


        for (kk = 0; kk < nq[nn]; kk++) {
            for (jj = 0; jj < nu[nn]; jj++) {
                for (ii = 0; ii < nu[nn]; ii++) {

                    fprintf(file, "%18.15e, ", MATEL(&qp->Hq[nn][kk], ii, jj));
                }
            }
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *RRq[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "Rq%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hRq = RRq;\n");


    // qq
    fprintf(file, "/* qq */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double qq%d[] = {", nn);


        for (kk = 0; kk < nq[nn]; kk++) {
            for (jj = 0; jj < nx[nn]; jj++) {

                //				fprintf(file, "%18.15e, ", VECEL(&qp->gq[nn][kk], nu[nn]+jj));
                fprintf(file, "%18.15e, ", MATEL(qp->DCt + nn, nu[nn] + jj, ng[nn] + kk));
            }
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *qqq[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "qq%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hqq = qqq;\n");


    // rq
    fprintf(file, "/* rq */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double rq%d[] = {", nn);


        for (kk = 0; kk < nq[nn]; kk++) {
            for (jj = 0; jj < nu[nn]; jj++) {
                fprintf(file, "%18.15e, ", MATEL(qp->DCt + nn, jj, ng[nn] + kk));
            }
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *rrq[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "rq%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hrq = rrq;\n");


    // uq
    fprintf(file, "/* uq */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double uq%d[] = {", nn);


        for (jj = 0; jj < nq[nn]; jj++) {
            fprintf(file, "%18.15e, ", -VECEL(qp->d + nn, 2 * nb[nn] + 2 * ng[nn] + nq[nn] + jj));
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *uuq[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "uq%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **huq = uuq;\n");


    // uq_mask
    fprintf(file, "/* uq_mask */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double uq_mask%d[] = {", nn);


        for (jj = 0; jj < nq[nn]; jj++) {
            fprintf(file, "%18.15e, ", VECEL(qp->d_mask + nn, 2 * nb[nn] + 2 * ng[nn] + nq[nn] + jj));
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *uuq_mask[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "uq_mask%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **huq_mask = uuq_mask;\n");


    // Zl
    fprintf(file, "/* Zl */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double Zl%d[] = {", nn);


        for (jj = 0; jj < ns[nn]; jj++) {
            fprintf(file, "%18.15e, ", VECEL(qp->Z + nn, jj));
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *ZZl[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "Zl%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hZl = ZZl;\n");


    // Zu
    fprintf(file, "/* Zu */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double Zu%d[] = {", nn);


        for (jj = 0; jj < ns[nn]; jj++) {
            fprintf(file, "%18.15e, ", VECEL(qp->Z + nn, ns[nn] + jj));
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *ZZu[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "Zu%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hZu = ZZu;\n");


    // zl
    fprintf(file, "/* zl */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double zl%d[] = {", nn);


        for (jj = 0; jj < ns[nn]; jj++) {
            fprintf(file, "%18.15e, ", VECEL(qp->rqz + nn, nu[nn] + nx[nn] + jj));
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *zzl[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "zl%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hzl = zzl;\n");


    // zu
    fprintf(file, "/* zu */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double zu%d[] = {", nn);


        for (jj = 0; jj < ns[nn]; jj++) {
            fprintf(file, "%18.15e, ", VECEL(qp->rqz + nn, nu[nn] + nx[nn] + ns[nn] + jj));
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *zzu[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "zu%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hzu = zzu;\n");


    // idxs_rev
    fprintf(file, "/* idxs_rev */\n");
    for (nn = 0; nn <= N; nn++) {
        fprintf(file, "static int idxs_rev%d[] = {", nn);
        for (jj = 0; jj < nb[nn] + ng[nn] + nq[nn]; jj++) {
            fprintf(file, "%d, ", qp->idxs_rev[nn][jj]);
        }
        fprintf(file, "};\n");
    }
    fprintf(file, "static int *iidxs_rev[] = {");
    for (nn = 0; nn <= N; nn++)
        fprintf(file, "idxs_rev%d, ", nn);
    fprintf(file, "};\n");
    fprintf(file, "int **hidxs_rev = iidxs_rev;\n");

    // idxs
    fprintf(file, "/* idxs */\n");
    for (nn = 0; nn <= N; nn++) {
        fprintf(file, "static int idxs%d[] = {", nn);
        for (jj = 0; jj < nb[nn] + ng[nn] + nq[nn]; jj++) {
            idx_tmp = qp->idxs_rev[nn][jj];
            if (idx_tmp != -1) {
                fprintf(file, "%d, ", jj);
            }
        }
        fprintf(file, "};\n");
    }
    fprintf(file, "static int *iidxs[] = {");
    for (nn = 0; nn <= N; nn++)
        fprintf(file, "idxs%d, ", nn);
    fprintf(file, "};\n");
    fprintf(file, "int **hidxs = iidxs;\n");

    // lls
    fprintf(file, "/* lls */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double lls%d[] = {", nn);


        for (jj = 0; jj < ns[nn]; jj++) {
            fprintf(file, "%18.15e, ", VECEL(qp->d + nn, 2 * nb[nn] + 2 * ng[nn] + 2 * nq[nn] + jj));
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *llls[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "lls%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hlls = llls;\n");


    // lus
    fprintf(file, "/* lus */\n");
    for (nn = 0; nn <= N; nn++) {

        fprintf(file, "static double lus%d[] = {", nn);


        for (jj = 0; jj < ns[nn]; jj++) {
            fprintf(file, "%18.15e, ", VECEL(qp->d + nn, 2 * nb[nn] + 2 * ng[nn] + 2 * nq[nn] + ns[nn] + jj));
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *llus[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "lus%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hlus = llus;\n");


    // idxe
    fprintf(file, "/* idxe */\n");
    for (nn = 0; nn <= N; nn++) {
        fprintf(file, "static int idxe%d[] = {", nn);
        for (jj = 0; jj < nbue[nn] + nbxe[nn] + nge[nn]; jj++) {
            fprintf(file, "%d, ", qp->idxe[nn][jj]);
        }
        fprintf(file, "};\n");
    }
    fprintf(file, "static int *iidxe[] = {");
    for (nn = 0; nn <= N; nn++)
        fprintf(file, "idxe%d, ", nn);
    fprintf(file, "};\n");
    fprintf(file, "int **hidxe = iidxe;\n");

    // XXX what follows is not part of the QP !!!

    // u_guess
    //	fprintf(file, "/* u_guess */\n");
    //	fprintf(file, "double **hu_guess;\n");

    // x_guess
    //	fprintf(file, "/* x_guess */\n");
    //	fprintf(file, "double **hx_guess;\n");

    // sl_guess
    //	fprintf(file, "/* sl_guess */\n");
    //	fprintf(file, "double **hsl_guess;\n");

    // su_guess
    //	fprintf(file, "/* su_guess */\n");
    //	fprintf(file, "double **hsu_guess;\n");

    fclose(file);
}


void d_ocp_qcqp_sol_print(struct d_ocp_qcqp_dim* qp_dim, struct d_ocp_qcqp_sol* qp_sol) {
    int ii;

    int N = qp_dim->N;
    int* nx = qp_dim->nx;
    int* nu = qp_dim->nu;
    int* nb = qp_dim->nb;
    int* ng = qp_dim->ng;
    int* nq = qp_dim->nq;
    int* ns = qp_dim->ns;

    printf("uxs =\n");
    for (ii = 0; ii <= N; ii++)
        print_tran_vec(nu[ii] + nx[ii] + 2 * ns[ii], &qp_sol->ux[ii], 0);

    printf("pi =\n");
    for (ii = 0; ii < N; ii++)
        print_tran_vec(nx[ii + 1], &qp_sol->pi[ii], 0);

    printf("lam =\n");
    for (ii = 0; ii <= N; ii++)
        print_tran_vec(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], &qp_sol->lam[ii], 0);

    printf("t =\n");
    for (ii = 0; ii <= N; ii++)
        print_tran_vec(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], &qp_sol->t[ii], 0);
}


void d_ocp_qcqp_ipm_arg_codegen(char* file_name, char* mode, struct d_ocp_qcqp_dim* qp_dim, struct d_ocp_qcqp_ipm_arg* arg) {
    int ii;

    FILE* file = fopen(file_name, mode);

    fprintf(file, "/***************\n* arg\n***************/\n");

    // iter_max
    fprintf(file, "/* mode */\n");
    fprintf(file, "int mode = %d;\n", arg->mode);
    // iter_max
    fprintf(file, "/* iter_max */\n");
    fprintf(file, "int iter_max = %d;\n", arg->iter_max);
    // alpha_min
    fprintf(file, "/* alpha_min */\n");
    fprintf(file, "double alpha_min = %18.15e;\n", arg->alpha_min);
    // mu0
    fprintf(file, "/* mu0 */\n");
    fprintf(file, "double mu0 = %18.15e;\n", arg->mu0);
    // tol_stat
    fprintf(file, "/* tol_stat */\n");
    fprintf(file, "double tol_stat = %18.15e;\n", arg->res_g_max);
    // tol_eq
    fprintf(file, "/* tol_eq */\n");
    fprintf(file, "double tol_eq = %18.15e;\n", arg->res_b_max);
    // tol_ineq
    fprintf(file, "/* tol_ineq */\n");
    fprintf(file, "double tol_ineq = %18.15e;\n", arg->res_d_max);
    // tol_comp
    fprintf(file, "/* tol_comp */\n");
    fprintf(file, "double tol_comp = %18.15e;\n", arg->res_m_max);
    // reg_prim
    fprintf(file, "/* reg_prim */\n");
    fprintf(file, "double reg_prim = %18.15e;\n", arg->reg_prim);
    // warm_start
    fprintf(file, "/* warm_start */\n");
    fprintf(file, "int warm_start = %d;\n", arg->warm_start);
    // pred_corr
    fprintf(file, "/* pred_corr */\n");
    fprintf(file, "int pred_corr = %d;\n", arg->pred_corr);
    // ric_alg
    fprintf(file, "/* ric_alg */\n");
    fprintf(file, "int ric_alg = %d;\n", arg->square_root_alg);
    // split_step
    fprintf(file, "/* split_step */\n");
    fprintf(file, "int split_step = %d;\n", arg->split_step);

    fclose(file);
}


void d_ocp_qcqp_res_print(struct d_ocp_qcqp_dim* qp_dim, struct d_ocp_qcqp_res* qp_res) {
    int ii;

    int N = qp_dim->N;
    int* nx = qp_dim->nx;
    int* nu = qp_dim->nu;
    int* nb = qp_dim->nb;
    int* ng = qp_dim->ng;
    int* nq = qp_dim->nq;
    int* ns = qp_dim->ns;

    printf("res_g =\n");
    for (ii = 0; ii <= N; ii++)
        print_exp_tran_vec(nu[ii] + nx[ii] + 2 * ns[ii], &qp_res->res_g[ii], 0);

    printf("res_b =\n");
    for (ii = 0; ii < N; ii++)
        print_exp_tran_vec(nx[ii + 1], &qp_res->res_b[ii], 0);

    printf("res_d =\n");
    for (ii = 0; ii <= N; ii++)
        print_exp_tran_vec(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], &qp_res->res_d[ii], 0);

    printf("res_m =\n");
    for (ii = 0; ii <= N; ii++)
        print_exp_tran_vec(2 * nb[ii] + 2 * ng[ii] + 2 * nq[ii] + 2 * ns[ii], &qp_res->res_m[ii], 0);
}
