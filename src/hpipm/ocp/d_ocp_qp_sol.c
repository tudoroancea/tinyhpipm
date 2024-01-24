#include <stdio.h>
#include <stdlib.h>

#include "hpipm/blas.h"
#include "hpipm/common.h"
#include "hpipm/ocp/d_ocp_qp.h"
#include "hpipm/ocp/d_ocp_qp_dim.h"
#include "hpipm/ocp/d_ocp_qp_sol.h"


hpipm_size_t d_ocp_qp_sol_strsize() {
    return sizeof(struct d_ocp_qp_sol);
}


hpipm_size_t d_ocp_qp_sol_memsize(struct d_ocp_qp_dim* dim) {

    // extract dim
    int N = dim->N;
    int* nx = dim->nx;
    int* nu = dim->nu;
    int* nb = dim->nb;
    int* ng = dim->ng;
    int* ns = dim->ns;

    // loop index
    int ii;

    int nvt = 0;
    int net = 0;
    int nct = 0;
    for (ii = 0; ii < N; ii++) {
        nvt += nu[ii] + nx[ii] + 2 * ns[ii];
        net += nx[ii + 1];
        nct += 2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii];
    }
    nvt += nu[ii] + nx[ii] + 2 * ns[ii];
    nct += 2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii];

    hpipm_size_t size = 0;

    size += 3 * (N + 1) * sizeof(struct vec);  // ux lam t
    size += 1 * N * sizeof(struct vec);  // pi

    size += 1 * memsize_vec(nvt);  // ux
    size += 1 * memsize_vec(net);  // pi
    size += 2 * memsize_vec(nct);  // lam t

    size = (size + 63) / 64 * 64;  // make multiple of typical cache line size
    size += 64;  // align to typical cache line size

    return size;
}


void d_ocp_qp_sol_create(struct d_ocp_qp_dim* dim, struct d_ocp_qp_sol* qp_sol, void* mem) {

    // loop index
    int ii;

    // zero memory (to avoid corrupted memory like e.g. NaN)
    hpipm_size_t memsize = d_ocp_qp_sol_memsize(dim);
    hpipm_zero_memset(memsize, mem);

    // extract dim
    int N = dim->N;
    int* nx = dim->nx;
    int* nu = dim->nu;
    int* nb = dim->nb;
    int* ng = dim->ng;
    int* ns = dim->ns;

    int nvt = 0;
    int net = 0;
    int nct = 0;
    for (ii = 0; ii < N; ii++) {
        nvt += nu[ii] + nx[ii] + 2 * ns[ii];
        net += nx[ii + 1];
        nct += 2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii];
    }
    nvt += nu[ii] + nx[ii] + 2 * ns[ii];
    nct += 2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii];


    // vector struct stuff
    struct vec* sv_ptr = (struct vec*) mem;

    qp_sol->ux = sv_ptr;
    sv_ptr += N + 1;
    qp_sol->pi = sv_ptr;
    sv_ptr += N;
    qp_sol->lam = sv_ptr;
    sv_ptr += N + 1;
    qp_sol->t = sv_ptr;
    sv_ptr += N + 1;


    // align to typical cache line size
    hpipm_size_t l_ptr = (hpipm_size_t) sv_ptr;
    l_ptr = (l_ptr + 63) / 64 * 64;


    // double stuff
    char* c_ptr;
    c_ptr = (char*) l_ptr;

    char* tmp_ptr;

    // ux
    tmp_ptr = c_ptr;
    c_ptr += memsize_vec(nvt);
    for (ii = 0; ii <= N; ii++) {
        create_vec(nu[ii] + nx[ii] + 2 * ns[ii], qp_sol->ux + ii, tmp_ptr);
        tmp_ptr += nu[ii] * sizeof(double);  // u
        tmp_ptr += nx[ii] * sizeof(double);  // x
        tmp_ptr += ns[ii] * sizeof(double);  // s_ls
        tmp_ptr += ns[ii] * sizeof(double);  // s_us
        dvecse(nu[ii] + nx[ii] + 2 * ns[ii], 0.0, qp_sol->ux + ii, 0);
    }
    // pi
    tmp_ptr = c_ptr;
    c_ptr += memsize_vec(net);
    for (ii = 0; ii < N; ii++) {
        create_vec(nx[ii + 1], qp_sol->pi + ii, tmp_ptr);
        tmp_ptr += nx[ii + 1] * sizeof(double);  // pi
        dvecse(nx[ii + 1], 0.0, qp_sol->pi + ii, 0);
    }
    // lam
    tmp_ptr = c_ptr;
    c_ptr += memsize_vec(nct);
    for (ii = 0; ii <= N; ii++) {
        create_vec(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], qp_sol->lam + ii, tmp_ptr);
        tmp_ptr += nb[ii] * sizeof(double);  // lb
        tmp_ptr += ng[ii] * sizeof(double);  // lg
        tmp_ptr += nb[ii] * sizeof(double);  // ub
        tmp_ptr += ng[ii] * sizeof(double);  // ug
        tmp_ptr += ns[ii] * sizeof(double);  // ls
        tmp_ptr += ns[ii] * sizeof(double);  // us
        dvecse(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], 0.0, qp_sol->lam + ii, 0);
    }
    // t
    tmp_ptr = c_ptr;
    c_ptr += memsize_vec(nct);
    for (ii = 0; ii <= N; ii++) {
        create_vec(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], qp_sol->t + ii, tmp_ptr);
        tmp_ptr += nb[ii] * sizeof(double);  // lb
        tmp_ptr += ng[ii] * sizeof(double);  // lg
        tmp_ptr += nb[ii] * sizeof(double);  // ub
        tmp_ptr += ng[ii] * sizeof(double);  // ug
        tmp_ptr += ns[ii] * sizeof(double);  // ls
        tmp_ptr += ns[ii] * sizeof(double);  // us
        dvecse(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], 0.0, qp_sol->t + ii, 0);
    }

    qp_sol->dim = dim;

    qp_sol->memsize = d_ocp_qp_sol_memsize(dim);


#if defined(RUNTIME_CHECKS)
    if (c_ptr > ((char*) mem) + qp_sol->memsize) {
        printf("\nCreate_ocp_qp_sol: outsize memory bounds!\n\n");
        exit(1);
    }
#endif
}


void d_ocp_qp_sol_copy_all(struct d_ocp_qp_sol* qp_sol_orig, struct d_ocp_qp_sol* qp_sol_dest) {

    int N = qp_sol_orig->dim->N;
    int* nx = qp_sol_orig->dim->nx;
    int* nu = qp_sol_orig->dim->nu;
    int* nb = qp_sol_orig->dim->nb;
    int* ng = qp_sol_orig->dim->ng;
    int* ns = qp_sol_orig->dim->ns;

    int ii;

    // copy dim pointer
    //	qp_sol_dest->dim = qp_sol_orig->dim;

    for (ii = 0; ii < N; ii++) {
        dveccp(nu[ii] + nx[ii] + 2 * ns[ii], qp_sol_orig->ux + ii, 0, qp_sol_dest->ux + ii, 0);
        dveccp(nx[ii + 1], qp_sol_orig->pi + ii, 0, qp_sol_dest->pi + ii, 0);
        dveccp(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], qp_sol_orig->lam + ii, 0, qp_sol_dest->lam + ii, 0);
        dveccp(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], qp_sol_orig->t + ii, 0, qp_sol_dest->t + ii, 0);
    }
    ii = N;
    dveccp(nu[ii] + nx[ii] + 2 * ns[ii], qp_sol_orig->ux + ii, 0, qp_sol_dest->ux + ii, 0);
    dveccp(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], qp_sol_orig->lam + ii, 0, qp_sol_dest->lam + ii, 0);
    dveccp(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], qp_sol_orig->t + ii, 0, qp_sol_dest->t + ii, 0);
}


void d_ocp_qp_sol_get_all(struct d_ocp_qp_sol* qp_sol, double** u, double** x, double** ls, double** us, double** pi, double** lam_lb, double** lam_ub, double** lam_lg, double** lam_ug, double** lam_ls, double** lam_us) {

    int N = qp_sol->dim->N;
    int* nx = qp_sol->dim->nx;
    int* nu = qp_sol->dim->nu;
    int* nb = qp_sol->dim->nb;
    int* ng = qp_sol->dim->ng;
    int* ns = qp_sol->dim->ns;

    int ii;

    for (ii = 0; ii < N; ii++) {
        unpack_vec(nx[ii + 1], qp_sol->pi + ii, 0, pi[ii], 1);
    }

    for (ii = 0; ii <= N; ii++) {
        unpack_vec(nu[ii], qp_sol->ux + ii, 0, u[ii], 1);
        unpack_vec(nx[ii], qp_sol->ux + ii, nu[ii], x[ii], 1);
        if (nb[ii] > 0) {
            unpack_vec(nb[ii], qp_sol->lam + ii, 0, lam_lb[ii], 1);
            unpack_vec(nb[ii], qp_sol->lam + ii, nb[ii] + ng[ii], lam_ub[ii], 1);
        }
        if (ng[ii] > 0) {
            unpack_vec(ng[ii], qp_sol->lam + ii, nb[ii], lam_lg[ii], 1);
            unpack_vec(ng[ii], qp_sol->lam + ii, 2 * nb[ii] + ng[ii], lam_ug[ii], 1);
        }
        if (ns[ii] > 0) {
            unpack_vec(ns[ii], qp_sol->ux + ii, nu[ii] + nx[ii], ls[ii], 1);
            unpack_vec(ns[ii], qp_sol->ux + ii, nu[ii] + nx[ii] + ns[ii], us[ii], 1);
            unpack_vec(ns[ii], qp_sol->lam + ii, 2 * nb[ii] + 2 * ng[ii], lam_ls[ii], 1);
            unpack_vec(ns[ii], qp_sol->lam + ii, 2 * nb[ii] + 2 * ng[ii] + ns[ii], lam_us[ii], 1);
        }
    }
}


void d_ocp_qp_sol_get_all_rowmaj(struct d_ocp_qp_sol* qp_sol, double** u, double** x, double** ls, double** us, double** pi, double** lam_lb, double** lam_ub, double** lam_lg, double** lam_ug, double** lam_ls, double** lam_us) {
    d_ocp_qp_sol_get_all(qp_sol, u, x, ls, us, pi, lam_lb, lam_ub, lam_lg, lam_ug, lam_ls, lam_us);
}


void d_ocp_qp_sol_set_all(double** u, double** x, double** ls, double** us, double** pi, double** lam_lb, double** lam_ub, double** lam_lg, double** lam_ug, double** lam_ls, double** lam_us, struct d_ocp_qp_sol* qp_sol) {

    int N = qp_sol->dim->N;
    int* nx = qp_sol->dim->nx;
    int* nu = qp_sol->dim->nu;
    int* nb = qp_sol->dim->nb;
    int* ng = qp_sol->dim->ng;
    int* ns = qp_sol->dim->ns;

    int ii;

    for (ii = 0; ii < N; ii++) {
        pack_vec(nx[ii + 1], pi[ii], 1, qp_sol->pi + ii, 0);
    }

    for (ii = 0; ii <= N; ii++) {
        pack_vec(nu[ii], u[ii], 1, qp_sol->ux + ii, 0);
        pack_vec(nx[ii], x[ii], 1, qp_sol->ux + ii, nu[ii]);
        if (nb[ii] > 0) {
            pack_vec(nb[ii], lam_lb[ii], 1, qp_sol->lam + ii, 0);
            pack_vec(nb[ii], lam_ub[ii], 1, qp_sol->lam + ii, nb[ii] + ng[ii]);
        }
        if (ng[ii] > 0) {
            pack_vec(ng[ii], lam_lg[ii], 1, qp_sol->lam + ii, nb[ii]);
            pack_vec(ng[ii], lam_ug[ii], 1, qp_sol->lam + ii, 2 * nb[ii] + ng[ii]);
        }
        if (ns[ii] > 0) {
            pack_vec(ns[ii], ls[ii], 1, qp_sol->ux + ii, nu[ii] + nx[ii]);
            pack_vec(ns[ii], us[ii], 1, qp_sol->ux + ii, nu[ii] + nx[ii] + ns[ii]);
            pack_vec(ns[ii], lam_ls[ii], 1, qp_sol->lam + ii, 2 * nb[ii] + 2 * ng[ii]);
            pack_vec(ns[ii], lam_us[ii], 1, qp_sol->lam + ii, 2 * nb[ii] + 2 * ng[ii] + ns[ii]);
        }
    }
}


void d_ocp_qp_sol_get(char* field, int stage, struct d_ocp_qp_sol* qp_sol, double* vec) {
    if (hpipm_strcmp(field, "u")) {
        d_ocp_qp_sol_get_u(stage, qp_sol, vec);
    } else if (hpipm_strcmp(field, "x")) {
        d_ocp_qp_sol_get_x(stage, qp_sol, vec);
    } else if (hpipm_strcmp(field, "sl")) {
        d_ocp_qp_sol_get_sl(stage, qp_sol, vec);
    } else if (hpipm_strcmp(field, "su")) {
        d_ocp_qp_sol_get_su(stage, qp_sol, vec);
    } else if (hpipm_strcmp(field, "pi")) {
        d_ocp_qp_sol_get_pi(stage, qp_sol, vec);
    } else if (hpipm_strcmp(field, "lam_lb")) {
        d_ocp_qp_sol_get_lam_lb(stage, qp_sol, vec);
    } else if (hpipm_strcmp(field, "lam_lbu")) {
        d_ocp_qp_sol_get_lam_lbu(stage, qp_sol, vec);
    } else if (hpipm_strcmp(field, "lam_lbx")) {
        d_ocp_qp_sol_get_lam_lbx(stage, qp_sol, vec);
    } else if (hpipm_strcmp(field, "lam_ub")) {
        d_ocp_qp_sol_get_lam_ub(stage, qp_sol, vec);
    } else if (hpipm_strcmp(field, "lam_ubu")) {
        d_ocp_qp_sol_get_lam_ubu(stage, qp_sol, vec);
    } else if (hpipm_strcmp(field, "lam_ubx")) {
        d_ocp_qp_sol_get_lam_ubx(stage, qp_sol, vec);
    } else if (hpipm_strcmp(field, "lam_lg")) {
        d_ocp_qp_sol_get_lam_lg(stage, qp_sol, vec);
    } else if (hpipm_strcmp(field, "lam_ug")) {
        d_ocp_qp_sol_get_lam_ug(stage, qp_sol, vec);
    } else {
        printf("error [d_ocp_qp_dim_get]: unknown field name '%s'. Exiting.\n", field);
        exit(1);
    }
}


void d_ocp_qp_sol_get_u(int stage, struct d_ocp_qp_sol* qp_sol, double* vec) {
    int* nu = qp_sol->dim->nu;
    unpack_vec(nu[stage], qp_sol->ux + stage, 0, vec, 1);
}


void d_ocp_qp_sol_get_x(int stage, struct d_ocp_qp_sol* qp_sol, double* vec) {
    int* nx = qp_sol->dim->nx;
    int* nu = qp_sol->dim->nu;
    unpack_vec(nx[stage], qp_sol->ux + stage, nu[stage], vec, 1);
}


void d_ocp_qp_sol_get_sl(int stage, struct d_ocp_qp_sol* qp_sol, double* vec) {
    int* nu = qp_sol->dim->nu;
    int* nx = qp_sol->dim->nx;
    int* ns = qp_sol->dim->ns;
    unpack_vec(ns[stage], qp_sol->ux + stage, nu[stage] + nx[stage], vec, 1);
}


void d_ocp_qp_sol_get_su(int stage, struct d_ocp_qp_sol* qp_sol, double* vec) {
    int* nu = qp_sol->dim->nu;
    int* nx = qp_sol->dim->nx;
    int* ns = qp_sol->dim->ns;
    unpack_vec(ns[stage], qp_sol->ux + stage, nu[stage] + nx[stage] + ns[stage], vec, 1);
}


void d_ocp_qp_sol_get_pi(int stage, struct d_ocp_qp_sol* qp_sol, double* vec) {
    int* nx = qp_sol->dim->nx;
    unpack_vec(nx[stage + 1], qp_sol->pi + stage, 0, vec, 1);
}


void d_ocp_qp_sol_get_lam_lb(int stage, struct d_ocp_qp_sol* qp_sol, double* vec) {
    int* nb = qp_sol->dim->nb;
    unpack_vec(nb[stage], qp_sol->lam + stage, 0, vec, 1);
}


void d_ocp_qp_sol_get_lam_lbu(int stage, struct d_ocp_qp_sol* qp_sol, double* vec) {
    int* nbu = qp_sol->dim->nbu;
    unpack_vec(nbu[stage], qp_sol->lam + stage, 0, vec, 1);
}


void d_ocp_qp_sol_get_lam_lbx(int stage, struct d_ocp_qp_sol* qp_sol, double* vec) {
    int* nbu = qp_sol->dim->nbu;
    int* nbx = qp_sol->dim->nbx;
    unpack_vec(nbx[stage], qp_sol->lam + stage, nbu[stage], vec, 1);
}


void d_ocp_qp_sol_get_lam_ub(int stage, struct d_ocp_qp_sol* qp_sol, double* vec) {
    int* nb = qp_sol->dim->nb;
    int* ng = qp_sol->dim->ng;
    unpack_vec(nb[stage], qp_sol->lam + stage, nb[stage] + ng[stage], vec, 1);
}


void d_ocp_qp_sol_get_lam_ubu(int stage, struct d_ocp_qp_sol* qp_sol, double* vec) {
    int* nb = qp_sol->dim->nb;
    int* nbu = qp_sol->dim->nbu;
    int* ng = qp_sol->dim->ng;
    unpack_vec(nbu[stage], qp_sol->lam + stage, nb[stage] + ng[stage], vec, 1);
}


void d_ocp_qp_sol_get_lam_ubx(int stage, struct d_ocp_qp_sol* qp_sol, double* vec) {
    int* nb = qp_sol->dim->nb;
    int* nbu = qp_sol->dim->nbu;
    int* nbx = qp_sol->dim->nbx;
    int* ng = qp_sol->dim->ng;
    unpack_vec(nbx[stage], qp_sol->lam + stage, nb[stage] + ng[stage] + nbu[stage], vec, 1);
}


void d_ocp_qp_sol_get_lam_lg(int stage, struct d_ocp_qp_sol* qp_sol, double* vec) {
    int* nb = qp_sol->dim->nb;
    int* ng = qp_sol->dim->ng;
    unpack_vec(ng[stage], qp_sol->lam + stage, nb[stage], vec, 1);
}


void d_ocp_qp_sol_get_lam_ug(int stage, struct d_ocp_qp_sol* qp_sol, double* vec) {
    int* nb = qp_sol->dim->nb;
    int* ng = qp_sol->dim->ng;
    unpack_vec(ng[stage], qp_sol->lam + stage, 2 * nb[stage] + ng[stage], vec, 1);
}


void d_ocp_qp_sol_set(char* field, int stage, double* vec, struct d_ocp_qp_sol* qp_sol) {
    if (hpipm_strcmp(field, "u")) {
        d_ocp_qp_sol_set_u(stage, vec, qp_sol);
    } else if (hpipm_strcmp(field, "x")) {
        d_ocp_qp_sol_set_x(stage, vec, qp_sol);
    } else if (hpipm_strcmp(field, "sl")) {
        d_ocp_qp_sol_set_sl(stage, vec, qp_sol);
    } else if (hpipm_strcmp(field, "su")) {
        d_ocp_qp_sol_set_su(stage, vec, qp_sol);
    } else {
        printf("error [d_ocp_qp_dim_set]: unknown field name '%s'. Exiting.\n", field);
        exit(1);
    }
}


void d_ocp_qp_sol_set_u(int stage, double* vec, struct d_ocp_qp_sol* qp_sol) {
    int* nu = qp_sol->dim->nu;
    pack_vec(nu[stage], vec, 1, qp_sol->ux + stage, 0);
}


void d_ocp_qp_sol_set_x(int stage, double* vec, struct d_ocp_qp_sol* qp_sol) {
    int* nu = qp_sol->dim->nu;
    int* nx = qp_sol->dim->nx;
    pack_vec(nx[stage], vec, 1, qp_sol->ux + stage, nu[stage]);
}


void d_ocp_qp_sol_set_sl(int stage, double* vec, struct d_ocp_qp_sol* qp_sol) {
    int* nu = qp_sol->dim->nu;
    int* nx = qp_sol->dim->nx;
    int* ns = qp_sol->dim->ns;
    pack_vec(ns[stage], vec, 1, qp_sol->ux + stage, nu[stage] + nx[stage]);
}


void d_ocp_qp_sol_set_su(int stage, double* vec, struct d_ocp_qp_sol* qp_sol) {
    int* nu = qp_sol->dim->nu;
    int* nx = qp_sol->dim->nx;
    int* ns = qp_sol->dim->ns;
    pack_vec(ns[stage], vec, 1, qp_sol->ux + stage, nu[stage] + nx[stage] + ns[stage]);
}
