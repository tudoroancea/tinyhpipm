#include <stdio.h>
#include <stdlib.h>

#include "tinyhpipm/blas.h"
#include "tinyhpipm/ocp/d_ocp_qp.h"
#include "tinyhpipm/ocp/d_ocp_qp_dim.h"
#include "tinyhpipm/ocp/d_ocp_qp_ipm.h"
#include "tinyhpipm/ocp/d_ocp_qp_sol.h"


void d_ocp_qp_dim_print(struct d_ocp_qp_dim* qp_dim) {
    int ii;

    int N = qp_dim->N;
    int* nx = qp_dim->nx;
    int* nu = qp_dim->nu;
    int* nbx = qp_dim->nbx;
    int* nbu = qp_dim->nbu;
    int* ng = qp_dim->ng;
    int* nsbx = qp_dim->nsbx;
    int* nsbu = qp_dim->nsbu;
    int* nsg = qp_dim->nsg;
    int* nbxe = qp_dim->nbxe;
    int* nbue = qp_dim->nbue;
    int* nge = qp_dim->nge;

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
}


void d_ocp_qp_dim_codegen(char* file_name, char* mode, struct d_ocp_qp_dim* qp_dim) {
    int ii;

    FILE* file = fopen(file_name, mode);

    int N = qp_dim->N;
    int* nx = qp_dim->nx;
    int* nu = qp_dim->nu;
    int* nbx = qp_dim->nbx;
    int* nbu = qp_dim->nbu;
    int* ng = qp_dim->ng;
    int* nsbx = qp_dim->nsbx;
    int* nsbu = qp_dim->nsbu;
    int* nsg = qp_dim->nsg;
    int* nbxe = qp_dim->nbxe;
    int* nbue = qp_dim->nbue;
    int* nge = qp_dim->nge;

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

    fclose(file);
}


void d_ocp_qp_print(struct d_ocp_qp_dim* dim, struct d_ocp_qp* qp) {
    int ii;

    int N = dim->N;
    int* nx = dim->nx;
    int* nu = dim->nu;
    int* nb = dim->nb;
    int* ng = dim->ng;
    int* ns = dim->ns;
    int* nbxe = dim->nbxe;
    int* nbue = dim->nbue;
    int* nge = dim->nge;

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
        print_tran_vec(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], qp->d + ii, 0);

    printf("d_mask =\n");
    for (ii = 0; ii <= N; ii++)
        print_tran_vec(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], qp->d_mask + ii, 0);

    printf("DCt =\n");
    for (ii = 0; ii <= N; ii++)
        print_mat(nu[ii] + nx[ii], ng[ii], qp->DCt + ii, 0, 0);

    printf("idxs_rev = \n");
    for (ii = 0; ii <= N; ii++)
        int_print_mat(1, nb[ii] + ng[ii], qp->idxs_rev[ii], 1);

    printf("m =\n");
    for (ii = 0; ii <= N; ii++)
        print_tran_vec(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], qp->m + ii, 0);

    printf("idxe = \n");
    for (ii = 0; ii <= N; ii++)
        int_print_mat(1, nbue[ii] + nbxe[ii] + nge[ii], qp->idxe[ii], 1);
}


void d_ocp_qp_codegen(char* file_name, char* mode, struct d_ocp_qp_dim* dim, struct d_ocp_qp* qp) {
    int nn, ii, jj, idx_tmp;

    FILE* file = fopen(file_name, mode);

    int N = dim->N;
    int* nx = dim->nx;
    int* nu = dim->nu;
    int* nb = dim->nb;
    int* ng = dim->ng;
    int* ns = dim->ns;
    int* nbxe = dim->nbxe;
    int* nbue = dim->nbue;
    int* nge = dim->nge;

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
                fprintf(file, "%18.15e, ", -VECEL(qp->d + nn, nb[nn] + ng[nn] + jj));
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
                fprintf(file, "%18.15e, ", VECEL(qp->d_mask + nn, nb[nn] + ng[nn] + jj));
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
                fprintf(file, "%18.15e, ", -VECEL(qp->d + nn, nb[nn] + ng[nn] + jj));
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
                fprintf(file, "%18.15e, ", VECEL(qp->d_mask + nn, nb[nn] + ng[nn] + jj));
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
            fprintf(file, "%18.15e, ", -VECEL(qp->d + nn, 2 * nb[nn] + ng[nn] + jj));
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


        for (jj = 0; jj < ng[nn]; jj++) {
            fprintf(file, "%18.15e, ", VECEL(qp->d_mask + nn, 2 * nb[nn] + ng[nn] + jj));
        }
        fprintf(file, "};\n");
    }

    fprintf(file, "static double *uug_mask[] = {");


    for (nn = 0; nn <= N; nn++)
        fprintf(file, "ug_mask%d, ", nn);
    fprintf(file, "};\n");

    fprintf(file, "double **hug_mask = uug_mask;\n");


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
        for (jj = 0; jj < nb[nn] + ng[nn]; jj++) {
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
        for (jj = 0; jj < nb[nn] + ng[nn]; jj++) {
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
            fprintf(file, "%18.15e, ", VECEL(qp->d + nn, 2 * nb[nn] + 2 * ng[nn] + jj));
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
            fprintf(file, "%18.15e, ", VECEL(qp->d + nn, 2 * nb[nn] + 2 * ng[nn] + ns[nn] + jj));
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


void d_ocp_qp_sol_print(struct d_ocp_qp_dim* qp_dim, struct d_ocp_qp_sol* qp_sol) {
    int ii;

    int N = qp_dim->N;
    int* nx = qp_dim->nx;
    int* nu = qp_dim->nu;
    int* nb = qp_dim->nb;
    int* ng = qp_dim->ng;
    int* ns = qp_dim->ns;

    printf("uxs =\n");
    for (ii = 0; ii <= N; ii++)
        print_tran_vec(nu[ii] + nx[ii] + 2 * ns[ii], &qp_sol->ux[ii], 0);

    printf("pi =\n");
    for (ii = 0; ii < N; ii++)
        print_tran_vec(nx[ii + 1], &qp_sol->pi[ii], 0);

    printf("lam =\n");
    for (ii = 0; ii <= N; ii++)
        print_tran_vec(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], &qp_sol->lam[ii], 0);

    printf("t =\n");
    for (ii = 0; ii <= N; ii++)
        print_tran_vec(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], &qp_sol->t[ii], 0);
}


void d_ocp_qp_ipm_arg_print(struct d_ocp_qp_dim* qp_dim, struct d_ocp_qp_ipm_arg* arg) {
    int ii;

    // mode
    printf("/* mode */\n");
    printf("int mode = %d;\n", arg->mode);
    // iter_max
    printf("/* iter_max */\n");
    printf("int iter_max = %d;\n", arg->iter_max);
    // alpha_min
    printf("/* alpha_min */\n");
    printf("double alpha_min = %18.15e;\n", arg->alpha_min);
    // mu0
    printf("/* mu0 */\n");
    printf("double mu0 = %18.15e;\n", arg->mu0);
    // tol_stat
    printf("/* tol_stat */\n");
    printf("double tol_stat = %18.15e;\n", arg->res_g_max);
    // tol_eq
    printf("/* tol_eq */\n");
    printf("double tol_eq = %18.15e;\n", arg->res_b_max);
    // tol_ineq
    printf("/* tol_ineq */\n");
    printf("double tol_ineq = %18.15e;\n", arg->res_d_max);
    // tol_comp
    printf("/* tol_comp */\n");
    printf("double tol_comp = %18.15e;\n", arg->res_m_max);
    // reg_prim
    printf("/* reg_prim */\n");
    printf("double reg_prim = %18.15e;\n", arg->reg_prim);
    // warm_start
    printf("/* warm_start */\n");
    printf("int warm_start = %d;\n", arg->warm_start);
    // pred_corr
    printf("/* pred_corr */\n");
    printf("int pred_corr = %d;\n", arg->pred_corr);
    // ric_alg
    printf("/* ric_alg */\n");
    printf("int ric_alg = %d;\n", arg->square_root_alg);
    // split_step
    printf("/* split_step */\n");
    printf("int split_step = %d;\n", arg->split_step);
}


void d_ocp_qp_ipm_arg_codegen(char* file_name, char* mode, struct d_ocp_qp_dim* qp_dim, struct d_ocp_qp_ipm_arg* arg) {
    int ii;

    FILE* file = fopen(file_name, mode);

    fprintf(file, "/***************\n* arg\n***************/\n");

    // mode
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


void d_ocp_qp_res_print(struct d_ocp_qp_dim* qp_dim, struct d_ocp_qp_res* qp_res) {
    int ii;

    int N = qp_dim->N;
    int* nx = qp_dim->nx;
    int* nu = qp_dim->nu;
    int* nb = qp_dim->nb;
    int* ng = qp_dim->ng;
    int* ns = qp_dim->ns;

    printf("res_g =\n");
    for (ii = 0; ii <= N; ii++)
        print_tran_vec(nu[ii] + nx[ii] + 2 * ns[ii], &qp_res->res_g[ii], 0);

    printf("res_b =\n");
    for (ii = 0; ii < N; ii++)
        print_tran_vec(nx[ii + 1], &qp_res->res_b[ii], 0);

    printf("res_d =\n");
    for (ii = 0; ii <= N; ii++)
        print_tran_vec(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], &qp_res->res_d[ii], 0);

    printf("res_m =\n");
    for (ii = 0; ii <= N; ii++)
        print_tran_vec(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], &qp_res->res_m[ii], 0);
}
