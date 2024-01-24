#include <stdio.h>
#include <stdlib.h>

#include "hpipm/blas.h"
#include "hpipm/common.h"
#include "hpipm/dense/d_dense_qp.h"
#include "hpipm/dense/d_dense_qp_dim.h"
#include "hpipm/dense/d_dense_qp_ipm.h"
#include "hpipm/dense/d_dense_qp_res.h"
#include "hpipm/dense/d_dense_qp_sol.h"


void d_dense_qp_dim_print(struct d_dense_qp_dim* qp_dim) {
    int ii;

    int nv = qp_dim->nv;
    int ne = qp_dim->ne;
    int nb = qp_dim->nb;
    int ng = qp_dim->ng;
    int nsb = qp_dim->nsb;
    int nsg = qp_dim->nsg;
    int ns = qp_dim->ns;

    printf("nv = %d\n\n", nv);
    printf("ne = %d\n\n", ne);
    printf("nb = %d\n\n", nb);
    printf("ng = %d\n\n", ng);
    printf("nsb = %d\n\n", nsb);
    printf("nsg = %d\n\n", nsg);
    printf("ns = %d\n\n", ns);
}


void d_dense_qp_dim_codegen(char* file_name, char* mode, struct d_dense_qp_dim* qp_dim) {
    int ii;

    FILE* file = fopen(file_name, mode);

    int nv = qp_dim->nv;
    int ne = qp_dim->ne;
    int nb = qp_dim->nb;
    int ng = qp_dim->ng;
    int nsb = qp_dim->nsb;
    int nsg = qp_dim->nsg;
    int ns = qp_dim->ns;

    fprintf(file, "/***************\n* dim\n***************/\n");

    // nv
    fprintf(file, "/* nv */\n");
    fprintf(file, "int nv = %d;\n", nv);
    // ne
    fprintf(file, "/* ne */\n");
    fprintf(file, "int ne = %d;\n", ne);
    // nb
    fprintf(file, "/* nb */\n");
    fprintf(file, "int nb = %d;\n", nb);
    // ng
    fprintf(file, "/* ng */\n");
    fprintf(file, "int ng = %d;\n", ng);
    // nsb
    fprintf(file, "/* nsb */\n");
    fprintf(file, "int nsb = %d;\n", nsb);
    // nsg
    fprintf(file, "/* nsg */\n");
    fprintf(file, "int nsg = %d;\n", nsg);
    // ns
    fprintf(file, "/* ns */\n");
    fprintf(file, "int ns = %d;\n", ns);

    fclose(file);
}


void d_dense_qp_print(struct d_dense_qp_dim* qp_dim, struct d_dense_qp* qp) {
    int ii;

    int nv = qp_dim->nv;
    int ne = qp_dim->ne;
    int nb = qp_dim->nb;
    int ng = qp_dim->ng;
    int ns = qp_dim->ns;

    printf("H = \n");
    print_mat(nv, nv, qp->Hv, 0, 0);

    printf("A = \n");
    print_mat(ne, nv, qp->A, 0, 0);

    printf("Ct = \n");
    print_mat(nv, ng, qp->Ct, 0, 0);

    printf("idxb = \n");
    int_print_mat(1, nb, qp->idxb, 1);

    printf("gz = \n");
    print_tran_vec(nv + 2 * ns, qp->gz, 0);

    printf("b = \n");
    print_tran_vec(ne, qp->b, 0);

    printf("d = \n");
    print_tran_vec(2 * nb + 2 * ng + 2 * ns, qp->d, 0);

    printf("d_mask = \n");
    print_tran_vec(2 * nb + 2 * ng + 2 * ns, qp->d_mask, 0);

    printf("m = \n");
    print_tran_vec(2 * nb + 2 * ng + 2 * ns, qp->m, 0);

    printf("Z = \n");
    print_tran_vec(2 * ns, qp->Z, 0);

    printf("idxs_rev = \n");
    int_print_mat(1, nb + ng, qp->idxs_rev, 1);
}


void d_dense_qp_codegen(char* file_name, char* mode, struct d_dense_qp_dim* qp_dim, struct d_dense_qp* qp) {
    int ii, jj;

    FILE* file = fopen(file_name, mode);

    int nv = qp_dim->nv;
    int ne = qp_dim->ne;
    int nb = qp_dim->nb;
    int ng = qp_dim->ng;
    int ns = qp_dim->ns;

    fprintf(file, "/***************\n* qp\n***************/\n");

    // H
    fprintf(file, "/* H */\n");
    fprintf(file, "static double HH[] = {");
    for (jj = 0; jj < nv; jj++) {
        for (ii = 0; ii < nv; ii++) {
            fprintf(file, "%18.15e, ", MATEL(qp->Hv, ii, jj));
        }
    }
    fprintf(file, "};\n");
    fprintf(file, "double *H = HH;\n");

    // A
    fprintf(file, "/* A */\n");
    fprintf(file, "static double AA[] = {");
    for (jj = 0; jj < nv; jj++) {
        for (ii = 0; ii < ne; ii++) {
            fprintf(file, "%18.15e, ", MATEL(qp->A, ii, jj));
        }
    }
    fprintf(file, "};\n");
    fprintf(file, "double *A = AA;\n");

    // C
    fprintf(file, "/* C */\n");
    fprintf(file, "static double CC[] = {");
    for (ii = 0; ii < nv; ii++) {
        for (jj = 0; jj < ng; jj++) {
            fprintf(file, "%18.15e, ", MATEL(qp->Ct, ii, jj));
        }
    }
    fprintf(file, "};\n");
    fprintf(file, "double *C = CC;\n");

    // idxb
    fprintf(file, "/* idxb */\n");
    fprintf(file, "static int iidxb[] = {");
    for (ii = 0; ii < nb; ii++)
        fprintf(file, "%d, ", qp->idxb[ii]);
    fprintf(file, "};\n");
    fprintf(file, "int *idxb = iidxb;\n");

    // g
    fprintf(file, "/* g */\n");
    fprintf(file, "static double gg[] = {");
    for (ii = 0; ii < nv; ii++) {
        fprintf(file, "%18.15e, ", VECEL(qp->gz, ii));
    }
    fprintf(file, "};\n");
    fprintf(file, "double *g = gg;\n");

    // zl
    fprintf(file, "/* zl */\n");
    fprintf(file, "static double zzl[] = {");
    for (ii = 0; ii < ns; ii++) {
        fprintf(file, "%18.15e, ", VECEL(qp->gz, nv + ii));
    }
    fprintf(file, "};\n");
    fprintf(file, "double *zl = zzl;\n");

    // zu
    fprintf(file, "/* zu */\n");
    fprintf(file, "static double zzu[] = {");
    for (ii = 0; ii < ns; ii++) {
        fprintf(file, "%18.15e, ", VECEL(qp->gz, nv + ns + ii));
    }
    fprintf(file, "};\n");
    fprintf(file, "double *zu = zzu;\n");

    // b
    fprintf(file, "/* b */\n");
    fprintf(file, "static double bb[] = {");
    for (ii = 0; ii < ne; ii++) {
        fprintf(file, "%18.15e, ", VECEL(qp->b, ii));
    }
    fprintf(file, "};\n");
    fprintf(file, "double *b = bb;\n");

    // lb
    fprintf(file, "/* lb */\n");
    fprintf(file, "static double llb[] = {");
    for (ii = 0; ii < nb; ii++) {
        fprintf(file, "%18.15e, ", VECEL(qp->d, ii));
    }
    fprintf(file, "};\n");
    fprintf(file, "double *lb = llb;\n");

    // lb_mask
    fprintf(file, "/* lb_mask */\n");
    fprintf(file, "static double llb_mask[] = {");
    for (ii = 0; ii < nb; ii++) {
        fprintf(file, "%18.15e, ", VECEL(qp->d_mask, ii));
    }
    fprintf(file, "};\n");
    fprintf(file, "double *lb_mask = llb_mask;\n");

    // ub
    fprintf(file, "/* ub */\n");
    fprintf(file, "static double uub[] = {");
    for (ii = 0; ii < nb; ii++) {
        fprintf(file, "%18.15e, ", -VECEL(qp->d, nb + ng + ii));
    }
    fprintf(file, "};\n");
    fprintf(file, "double *ub = uub;\n");

    // ub_mask
    fprintf(file, "/* ub_mask */\n");
    fprintf(file, "static double uub_mask[] = {");
    for (ii = 0; ii < nb; ii++) {
        fprintf(file, "%18.15e, ", VECEL(qp->d_mask, nb + ng + ii));
    }
    fprintf(file, "};\n");
    fprintf(file, "double *ub_mask = uub_mask;\n");

    // lg
    fprintf(file, "/* lg */\n");
    fprintf(file, "static double llg[] = {");
    for (ii = 0; ii < ng; ii++) {
        fprintf(file, "%18.15e, ", VECEL(qp->d, nb + ii));
    }
    fprintf(file, "};\n");
    fprintf(file, "double *lg = llg;\n");

    // lg_mask
    fprintf(file, "/* lg_mask */\n");
    fprintf(file, "static double llg_mask[] = {");
    for (ii = 0; ii < ng; ii++) {
        fprintf(file, "%18.15e, ", VECEL(qp->d_mask, nb + ii));
    }
    fprintf(file, "};\n");
    fprintf(file, "double *lg_mask = llg_mask;\n");

    // ug
    fprintf(file, "/* ug */\n");
    fprintf(file, "static double uug[] = {");
    for (ii = 0; ii < ng; ii++) {
        fprintf(file, "%18.15e, ", -VECEL(qp->d, 2 * nb + ng + ii));
    }
    fprintf(file, "};\n");
    fprintf(file, "double *ug = uug;\n");

    // ug_mask
    fprintf(file, "/* ug_mask */\n");
    fprintf(file, "static double uug_mask[] = {");
    for (ii = 0; ii < ng; ii++) {
        fprintf(file, "%18.15e, ", VECEL(qp->d_mask, 2 * nb + ng + ii));
    }
    fprintf(file, "};\n");
    fprintf(file, "double *ug_mask = uug_mask;\n");

    // lls
    fprintf(file, "/* lls */\n");
    fprintf(file, "static double llls[] = {");
    for (ii = 0; ii < ns; ii++) {
        fprintf(file, "%18.15e, ", VECEL(qp->d, 2 * nb + 2 * ng + ii));
    }
    fprintf(file, "};\n");
    fprintf(file, "double *lls = llls;\n");

    // lls_mask
    fprintf(file, "/* lls_mask */\n");
    fprintf(file, "static double llls_mask[] = {");
    for (ii = 0; ii < ns; ii++) {
        fprintf(file, "%18.15e, ", VECEL(qp->d_mask, 2 * nb + 2 * ng + ii));
    }
    fprintf(file, "};\n");
    fprintf(file, "double *lls_mask = llls_mask;\n");

    // lus
    fprintf(file, "/* lus */\n");
    fprintf(file, "static double llus[] = {");
    for (ii = 0; ii < ns; ii++) {
        fprintf(file, "%18.15e, ", VECEL(qp->d, 2 * nb + 2 * ng + ns + ii));
    }
    fprintf(file, "};\n");
    fprintf(file, "double *lus = llus;\n");

    // lus_mask
    fprintf(file, "/* lus_mask */\n");
    fprintf(file, "static double llus_mask[] = {");
    for (ii = 0; ii < ns; ii++) {
        fprintf(file, "%18.15e, ", VECEL(qp->d_mask, 2 * nb + 2 * ng + ns + ii));
    }
    fprintf(file, "};\n");
    fprintf(file, "double *lus_mask = llus_mask;\n");

    printf("Z = \n");
    print_tran_vec(2 * ns, qp->Z, 0);
    // zl
    fprintf(file, "/* Zl */\n");
    fprintf(file, "static double ZZl[] = {");
    for (ii = 0; ii < ns; ii++) {
        fprintf(file, "%18.15e, ", VECEL(qp->Z, ii));
    }
    fprintf(file, "};\n");
    fprintf(file, "double *Zl = ZZl;\n");

    // Zu
    fprintf(file, "/* Zu */\n");
    fprintf(file, "static double ZZu[] = {");
    for (ii = 0; ii < ns; ii++) {
        fprintf(file, "%18.15e, ", VECEL(qp->Z, ns + ii));
    }
    fprintf(file, "};\n");
    fprintf(file, "double *Zu = ZZu;\n");

    // idxs_rev
    fprintf(file, "/* idxs_rev */\n");
    fprintf(file, "static int iidxs_rev[] = {");
    for (ii = 0; ii < nb + ng; ii++)
        fprintf(file, "%d, ", qp->idxs_rev[ii]);
    fprintf(file, "};\n");
    fprintf(file, "int *idxs_rev = iidxs_rev;\n");

    fclose(file);
}


void d_dense_qp_sol_print(struct d_dense_qp_dim* qp_dim, struct d_dense_qp_sol* qp_sol) {
    int ii;

    int nv = qp_dim->nv;
    int ne = qp_dim->ne;
    int nb = qp_dim->nb;
    int ng = qp_dim->ng;
    int ns = qp_dim->ns;

    printf("v = \n");
    print_tran_vec(nv + 2 * ns, qp_sol->v, 0);

    printf("pi = \n");
    print_tran_vec(ne, qp_sol->pi, 0);

    printf("lam = \n");
    print_tran_vec(2 * nb + 2 * ng + 2 * ns, qp_sol->lam, 0);

    printf("t = \n");
    print_tran_vec(2 * nb + 2 * ng + 2 * ns, qp_sol->t, 0);
}


void d_dense_qp_res_print(struct d_dense_qp_dim* qp_dim, struct d_dense_qp_res* qp_res) {
    int ii;

    int nv = qp_dim->nv;
    int ne = qp_dim->ne;
    int nb = qp_dim->nb;
    int ng = qp_dim->ng;
    int ns = qp_dim->ns;

    printf("res_g = \n");
    print_tran_vec(nv + 2 * ns, qp_res->res_g, 0);

    printf("res_b = \n");
    print_tran_vec(ne, qp_res->res_b, 0);

    printf("res_d = \n");
    print_tran_vec(2 * nb + 2 * ng + 2 * ns, qp_res->res_d, 0);

    printf("res_m = \n");
    print_tran_vec(2 * nb + 2 * ng + 2 * ns, qp_res->res_m, 0);
}


void d_dense_qp_ipm_arg_print(struct d_dense_qp_dim* qp_dim, struct d_dense_qp_ipm_arg* arg) {
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
    // reg_dual
    printf("/* reg_dual */\n");
    printf("double reg_dual = %18.15e;\n", arg->reg_dual);
    // warm_start
    printf("/* warm_start */\n");
    printf("int warm_start = %d;\n", arg->warm_start);
    // pred_corr
    printf("/* pred_corr */\n");
    printf("int pred_corr = %d;\n", arg->pred_corr);
    // split_step
    printf("/* split_step */\n");
    printf("int split_step = %d;\n", arg->split_step);
}


void d_dense_qp_ipm_arg_codegen(char* file_name, char* mode, struct d_dense_qp_dim* qp_dim, struct d_dense_qp_ipm_arg* arg) {
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
    // reg_dual
    fprintf(file, "/* reg_dual */\n");
    fprintf(file, "double reg_dual = %18.15e;\n", arg->reg_dual);
    // warm_start
    fprintf(file, "/* warm_start */\n");
    fprintf(file, "int warm_start = %d;\n", arg->warm_start);
    // pred_corr
    fprintf(file, "/* pred_corr */\n");
    fprintf(file, "int pred_corr = %d;\n", arg->pred_corr);
    // split_step
    fprintf(file, "/* split_step */\n");
    fprintf(file, "int split_step = %d;\n", arg->split_step);

    fclose(file);
}
