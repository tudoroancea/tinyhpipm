#include <stdio.h>
#include <stdlib.h>

#include "tinyhpipm/blas.h"
#include "tinyhpipm/common.h"
#include "tinyhpipm/dense/d_dense_qp.h"
#include "tinyhpipm/dense/d_dense_qp_dim.h"


hpipm_size_t d_dense_qp_strsize() {
    return sizeof(struct d_dense_qp);
}


hpipm_size_t d_dense_qp_memsize(struct d_dense_qp_dim* dim) {

    int nv = dim->nv;
    int ne = dim->ne;
    int nb = dim->nb;
    int ng = dim->ng;
    int ns = dim->ns;

    hpipm_size_t size = 0;

    size += 6 * sizeof(struct vec);  // gz b d m Z d_mask
    size += 3 * sizeof(struct mat);  // Hv A Ct

    size += 1 * memsize_vec(nv + 2 * ns);  // g
    size += 1 * memsize_vec(ne);  // b
    size += 3 * memsize_vec(2 * nb + 2 * ng + 2 * ns);  // d m d_mask
    size += 1 * memsize_vec(2 * ns);  // Z
    size += 1 * nb * sizeof(int);  // idxb
    size += 1 * (nb + ng) * sizeof(int);  // idxs_rev

    size += 1 * memsize_mat(nv + 1, nv);  // Hv
    size += 1 * memsize_mat(ne, nv);  // A
    size += 1 * memsize_mat(nv, ng);  // Ct

    size = (size + 63) / 64 * 64;  // make multiple of typical cache line size
    size += 1 * 64;  // align once to typical cache line size

    return size;
}


void d_dense_qp_create(struct d_dense_qp_dim* dim, struct d_dense_qp* qp, void* mem) {

    int ii;

    // zero memory (to avoid corrupted memory like e.g. NaN)
    hpipm_size_t memsize = d_dense_qp_memsize(dim);
    hpipm_zero_memset(memsize, mem);

    // extract dim
    int nv = dim->nv;
    int ne = dim->ne;
    int nb = dim->nb;
    int ng = dim->ng;
    int ns = dim->ns;


    // matrix struct stuff
    struct mat* sm_ptr = (struct mat*) mem;

    qp->Hv = sm_ptr;
    sm_ptr += 1;

    qp->A = sm_ptr;
    sm_ptr += 1;

    qp->Ct = sm_ptr;
    sm_ptr += 1;


    // vector struct stuff
    struct vec* sv_ptr = (struct vec*) sm_ptr;

    qp->gz = sv_ptr;
    sv_ptr += 1;

    qp->b = sv_ptr;
    sv_ptr += 1;

    qp->d = sv_ptr;
    sv_ptr += 1;

    qp->d_mask = sv_ptr;
    sv_ptr += 1;

    qp->m = sv_ptr;
    sv_ptr += 1;

    qp->Z = sv_ptr;
    sv_ptr += 1;


    // int stuff
    int* i_ptr;
    i_ptr = (int*) sv_ptr;

    // idxb
    qp->idxb = i_ptr;
    i_ptr += nb;

    // idxs_rev
    qp->idxs_rev = i_ptr;
    i_ptr += nb + ng;
    for (ii = 0; ii < nb + ng; ii++)
        qp->idxs_rev[ii] = -1;


    // align to typical cache line size
    hpipm_size_t s_ptr = (hpipm_size_t) i_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;


    //  stuff
    char* c_ptr;
    c_ptr = (char*) s_ptr;

    create_mat(nv + 1, nv, qp->Hv, c_ptr);
    c_ptr += qp->Hv->memsize;

    create_mat(ne, nv, qp->A, c_ptr);
    c_ptr += qp->A->memsize;

    create_mat(nv, ng, qp->Ct, c_ptr);
    c_ptr += qp->Ct->memsize;

    create_vec(nv + 2 * ns, qp->gz, c_ptr);
    c_ptr += qp->gz->memsize;

    create_vec(ne, qp->b, c_ptr);
    c_ptr += qp->b->memsize;

    create_vec(2 * nb + 2 * ng + 2 * ns, qp->d, c_ptr);
    c_ptr += qp->d->memsize;

    create_vec(2 * nb + 2 * ng + 2 * ns, qp->d_mask, c_ptr);
    c_ptr += qp->d_mask->memsize;

    create_vec(2 * nb + 2 * ng + 2 * ns, qp->m, c_ptr);
    c_ptr += qp->m->memsize;

    create_vec(2 * ns, qp->Z, c_ptr);
    c_ptr += qp->Z->memsize;


    // default initialization
    dvecse(2 * nb + 2 * ng + 2 * ns, 1.0, qp->d_mask, 0);


    qp->dim = dim;

    qp->memsize = d_dense_qp_memsize(dim);


#if defined(RUNTIME_CHECKS)
    if (c_ptr > ((char*) mem) + qp->memsize) {
        printf("\nCreate_ocp_qp: outside memory bounds!\n\n");
        exit(1);
    }
#endif
}


void d_dense_qp_set_all(double* H, double* g, double* A, double* b, int* idxb, double* d_lb, double* d_ub, double* C, double* d_lg, double* d_ug, double* Zl, double* Zu, double* zl, double* zu, int* idxs, double* d_ls, double* d_us, struct d_dense_qp* qp) {

    int ii, jj;

    int nv = qp->dim->nv;
    int ne = qp->dim->ne;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;
    int ns = qp->dim->ns;

    pack_mat(nv, nv, H, nv, qp->Hv, 0, 0);
    pack_vec(nv, g, 1, qp->gz, 0);
    if (ne > 0) {
        pack_mat(ne, nv, A, ne, qp->A, 0, 0);
        pack_vec(ne, b, 1, qp->b, 0);
    }
    if (nb > 0) {
        for (ii = 0; ii < nb; ii++) qp->idxb[ii] = idxb[ii];
        pack_vec(nb, d_lb, 1, qp->d, 0);
        pack_vec(nb, d_ub, 1, qp->d, nb + ng);
        dvecsc(nb, -1.0, qp->d, nb + ng);
        dvecse(nb, 0.0, qp->m, 0);
        dvecse(nb, 0.0, qp->m, nb + ng);
    }
    if (ng > 0) {
        pack_tran_mat(ng, nv, C, ng, qp->Ct, 0, 0);
        pack_vec(ng, d_lg, 1, qp->d, nb);
        pack_vec(ng, d_ug, 1, qp->d, 2 * nb + ng);
        dvecsc(ng, -1.0, qp->d, 2 * nb + ng);
        dvecse(ng, 0.0, qp->m, nb);
        dvecse(ng, 0.0, qp->m, 2 * nb + ng);
    }
    if (ns > 0) {
        for (ii = 0; ii < ns; ii++) {
            qp->idxs_rev[idxs[ii]] = ii;
        }
        pack_vec(ns, Zl, 1, qp->Z, 0);
        pack_vec(ns, Zu, 1, qp->Z, ns);
        pack_vec(ns, zl, 1, qp->gz, nv);
        pack_vec(ns, zu, 1, qp->gz, nv + ns);
        pack_vec(ns, d_ls, 1, qp->d, 2 * nb + 2 * ng);
        pack_vec(ns, d_us, 1, qp->d, 2 * nb + 2 * ng + ns);
        dvecse(ns, 0.0, qp->m, 2 * nb + 2 * ng);
        dvecse(ns, 0.0, qp->m, 2 * nb + 2 * ng + ns);
    }
}


void d_dense_qp_get_all(struct d_dense_qp* qp, double* H, double* g, double* A, double* b, int* idxb, double* d_lb, double* d_ub, double* C, double* d_lg, double* d_ug, double* Zl, double* Zu, double* zl, double* zu, int* idxs, double* d_ls, double* d_us) {

    int ii, idx_tmp;

    int nv = qp->dim->nv;
    int ne = qp->dim->ne;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;
    int ns = qp->dim->ns;

    unpack_mat(nv, nv, qp->Hv, 0, 0, H, nv);
    unpack_vec(nv, qp->gz, 0, g, 1);
    if (ne > 0) {
        unpack_mat(ne, nv, qp->A, 0, 0, A, ne);
        unpack_vec(ne, qp->b, 0, b, 1);
    }
    if (nb > 0) {
        for (ii = 0; ii < nb; ii++) idxb[ii] = qp->idxb[ii];
        unpack_vec(nb, qp->d, 0, d_lb, 1);
        unpack_vec(nb, qp->d, nb + ng, d_ub, 1);
        for (ii = 0; ii < nb; ii++) d_ub[ii] = -d_ub[ii];
    }
    if (ng > 0) {
        unpack_tran_mat(nv, ng, qp->Ct, 0, 0, C, ng);
        unpack_vec(ng, qp->d, nb, d_lg, 1);
        unpack_vec(ng, qp->d, 2 * nb + ng, d_ug, 1);
        for (ii = 0; ii < ng; ii++) d_ug[ii] = -d_ug[ii];
    }
    if (ns > 0) {
        // TODO only valid if there is one slack variable per soft constraint !!!
        for (ii = 0; ii < nb + ng; ii++) {
            idx_tmp = qp->idxs_rev[ii];
            if (idx_tmp != -1) {
                idxs[idx_tmp] = ii;
            }
        }
        unpack_vec(ns, qp->Z, 0, Zl, 1);
        unpack_vec(ns, qp->Z, ns, Zu, 1);
        unpack_vec(ns, qp->gz, nv, zl, 1);
        unpack_vec(ns, qp->gz, nv + ns, zu, 1);
        unpack_vec(ns, qp->d, 2 * nb + 2 * ng, d_ls, 1);
        unpack_vec(ns, qp->d, 2 * nb + 2 * ng + ns, d_us, 1);
    }
}


void d_dense_qp_set(char* field, void* value, struct d_dense_qp* qp) {
    double* r_ptr;
    int* i_ptr;

    // matrices
    if (hpipm_strcmp(field, "H")) {
        d_dense_qp_set_H(value, qp);
    } else if (hpipm_strcmp(field, "A")) {
        d_dense_qp_set_A(value, qp);
    } else if (hpipm_strcmp(field, "C")) {
        d_dense_qp_set_C(value, qp);
    }
    // vectors
    else if (hpipm_strcmp(field, "g")) {
        d_dense_qp_set_g(value, qp);
    } else if (hpipm_strcmp(field, "b")) {
        d_dense_qp_set_b(value, qp);
    } else if (hpipm_strcmp(field, "lb")) {
        d_dense_qp_set_lb(value, qp);
    } else if (hpipm_strcmp(field, "lb_mask")) {
        d_dense_qp_set_lb_mask(value, qp);
    } else if (hpipm_strcmp(field, "ub")) {
        d_dense_qp_set_ub(value, qp);
    } else if (hpipm_strcmp(field, "ub_mask")) {
        d_dense_qp_set_ub_mask(value, qp);
    } else if (hpipm_strcmp(field, "lg")) {
        d_dense_qp_set_lg(value, qp);
    } else if (hpipm_strcmp(field, "lg_mask")) {
        d_dense_qp_set_lg_mask(value, qp);
    } else if (hpipm_strcmp(field, "ug")) {
        d_dense_qp_set_ug(value, qp);
    } else if (hpipm_strcmp(field, "ug_mask")) {
        d_dense_qp_set_ug_mask(value, qp);
    } else if (hpipm_strcmp(field, "Zl")) {
        d_dense_qp_set_Zl(value, qp);
    } else if (hpipm_strcmp(field, "Zu")) {
        d_dense_qp_set_Zu(value, qp);
    } else if (hpipm_strcmp(field, "zl")) {
        d_dense_qp_set_zl(value, qp);
    } else if (hpipm_strcmp(field, "zu")) {
        d_dense_qp_set_zu(value, qp);
    } else if (hpipm_strcmp(field, "lls")) {
        d_dense_qp_set_lls(value, qp);
    } else if (hpipm_strcmp(field, "lls_mask")) {
        d_dense_qp_set_lls_mask(value, qp);
    } else if (hpipm_strcmp(field, "lus")) {
        d_dense_qp_set_lus(value, qp);
    } else if (hpipm_strcmp(field, "lus_mask")) {
        d_dense_qp_set_lus_mask(value, qp);
    }
    // int
    else if (hpipm_strcmp(field, "idxb")) {
        d_dense_qp_set_idxb(value, qp);
    } else if (hpipm_strcmp(field, "Jb")) {
        d_dense_qp_set_Jb(value, qp);
    } else if (hpipm_strcmp(field, "idxs")) {
        d_dense_qp_set_idxs(value, qp);
    } else if (hpipm_strcmp(field, "idxs_rev")) {
        d_dense_qp_set_idxs_rev(value, qp);
    } else if (hpipm_strcmp(field, "Jsb")) {
        d_dense_qp_set_Jsb(value, qp);
    } else if (hpipm_strcmp(field, "Jsg")) {
        d_dense_qp_set_Jsg(value, qp);
    } else {
        printf("error: d_dense_qp_set: wrong field name '%s'. Exiting.\n", field);
        exit(1);
    }
}


void d_dense_qp_set_H(double* H, struct d_dense_qp* qp) {

    int nv = qp->dim->nv;

    pack_mat(nv, nv, H, nv, qp->Hv, 0, 0);
}


void d_dense_qp_set_g(double* g, struct d_dense_qp* qp) {

    int nv = qp->dim->nv;

    pack_vec(nv, g, 1, qp->gz, 0);
}


void d_dense_qp_set_A(double* A, struct d_dense_qp* qp) {

    int nv = qp->dim->nv;
    int ne = qp->dim->ne;

    pack_mat(ne, nv, A, ne, qp->A, 0, 0);
}


void d_dense_qp_set_b(double* b, struct d_dense_qp* qp) {

    int ne = qp->dim->ne;

    pack_vec(ne, b, 1, qp->b, 0);
}


void d_dense_qp_set_idxb(int* idxb, struct d_dense_qp* qp) {

    int ii;
    int nb = qp->dim->nb;

    for (ii = 0; ii < nb; ii++) qp->idxb[ii] = idxb[ii];
}


void d_dense_qp_set_Jb(double* Jb, struct d_dense_qp* qp) {
    // extract dim
    int nv = qp->dim->nv;
    int nb = qp->dim->nb;

    int ii, jj, jj0;
    for (ii = 0; ii < nb; ii++) {
        jj0 = -1;
        for (jj = 0; jj < nv & jj0 == -1; jj++) {
            if (Jb[ii + jj * nb] != 0.0) {
                jj0 = jj;
                qp->idxb[ii] = jj;
            }
        }
    }
}


void d_dense_qp_set_lb(double* lb, struct d_dense_qp* qp) {

    int nb = qp->dim->nb;

    pack_vec(nb, lb, 1, qp->d, 0);
    dvecse(nb, 0.0, qp->m, 0);
}


void d_dense_qp_set_lb_mask(double* lb_mask, struct d_dense_qp* qp) {

    int nb = qp->dim->nb;

    int ii;

    for (ii = 0; ii < nb; ii++)
        if (lb_mask[ii] == 0.0)
            VECEL(qp->d_mask, ii) = 0.0;
        else
            VECEL(qp->d_mask, ii) = 1.0;
}


void d_dense_qp_set_ub(double* ub, struct d_dense_qp* qp) {

    int nb = qp->dim->nb;
    int ng = qp->dim->ng;

    pack_vec(nb, ub, 1, qp->d, nb + ng);
    dvecsc(nb, -1.0, qp->d, nb + ng);
    dvecse(nb, 0.0, qp->m, nb + ng);
}


void d_dense_qp_set_ub_mask(double* ub_mask, struct d_dense_qp* qp) {

    int nb = qp->dim->nb;
    int ng = qp->dim->ng;

    int ii;

    for (ii = 0; ii < nb; ii++)
        if (ub_mask[ii] == 0.0)
            VECEL(qp->d_mask, nb + ng + ii) = 0.0;
        else
            VECEL(qp->d_mask, nb + ng + ii) = 1.0;
}


void d_dense_qp_set_C(double* C, struct d_dense_qp* qp) {

    int nv = qp->dim->nv;
    int ng = qp->dim->ng;

    pack_tran_mat(ng, nv, C, ng, qp->Ct, 0, 0);
}


void d_dense_qp_set_lg(double* lg, struct d_dense_qp* qp) {

    int nb = qp->dim->nb;
    int ng = qp->dim->ng;

    pack_vec(ng, lg, 1, qp->d, nb);
    dvecse(ng, 0.0, qp->m, nb);
}


void d_dense_qp_set_lg_mask(double* lg_mask, struct d_dense_qp* qp) {

    int nb = qp->dim->nb;
    int ng = qp->dim->ng;

    int ii;

    for (ii = 0; ii < ng; ii++)
        if (lg_mask[ii] == 0.0)
            VECEL(qp->d_mask, nb + ii) = 0.0;
        else
            VECEL(qp->d_mask, nb + ii) = 1.0;
}


void d_dense_qp_set_ug(double* ug, struct d_dense_qp* qp) {

    int nb = qp->dim->nb;
    int ng = qp->dim->ng;

    pack_vec(ng, ug, 1, qp->d, 2 * nb + ng);
    dvecsc(ng, -1.0, qp->d, 2 * nb + ng);
    dvecse(ng, 0.0, qp->m, 2 * nb + ng);
}


void d_dense_qp_set_ug_mask(double* ug_mask, struct d_dense_qp* qp) {

    int nb = qp->dim->nb;
    int ng = qp->dim->ng;

    int ii;

    for (ii = 0; ii < ng; ii++)
        if (ug_mask[ii] == 0.0)
            VECEL(qp->d_mask, 2 * nb + ng + ii) = 0.0;
        else
            VECEL(qp->d_mask, 2 * nb + ng + ii) = 1.0;
}


void d_dense_qp_set_idxs(int* idxs, struct d_dense_qp* qp) {

    int ii;
    int ns = qp->dim->ns;

    for (ii = 0; ii < ns; ii++) {
        qp->idxs_rev[idxs[ii]] = ii;
    }
}


void d_dense_qp_set_idxs_rev(int* idxs_rev, struct d_dense_qp* qp) {

    int ii;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;

    for (ii = 0; ii < nb + ng; ii++) {
        qp->idxs_rev[ii] = idxs_rev[ii];
    }
}


void d_dense_qp_set_Jsb(double* Jsb, struct d_dense_qp* qp) {
    // extract dim
    int nv = qp->dim->nv;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;
    int ns = qp->dim->ns;

    int ii, jj, jj0, idx_tmp;
    // compute nb part of idxs_rev
    for (ii = 0; ii < nb; ii++) {
        jj0 = -1;
        for (jj = 0; jj < ns & jj0 == -1; jj++) {
            if (Jsb[ii + jj * nb] != 0.0) {
                jj0 = jj;
                qp->idxs_rev[0 + ii] = jj;
            }
        }
    }
}


void d_dense_qp_set_Jsg(double* Jsg, struct d_dense_qp* qp) {
    // extract dim
    int nv = qp->dim->nv;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;
    int ns = qp->dim->ns;

    int ii, jj, jj0, idx_tmp;
    // compute ng part of idxs_rev
    for (ii = 0; ii < ng; ii++) {
        jj0 = -1;
        for (jj = 0; jj < ns & jj0 == -1; jj++) {
            if (Jsg[ii + jj * ng] != 0.0) {
                jj0 = jj;
                qp->idxs_rev[nb + ii] = jj;
            }
        }
    }
}


void d_dense_qp_set_Zl(double* Zl, struct d_dense_qp* qp) {

    int ns = qp->dim->ns;

    pack_vec(ns, Zl, 1, qp->Z, 0);
}


void d_dense_qp_set_Zu(double* Zu, struct d_dense_qp* qp) {

    int ns = qp->dim->ns;

    pack_vec(ns, Zu, 1, qp->Z, ns);
}


void d_dense_qp_set_zl(double* zl, struct d_dense_qp* qp) {

    int ns = qp->dim->ns;
    int nv = qp->dim->nv;

    pack_vec(ns, zl, 1, qp->gz, nv);
}


void d_dense_qp_set_zu(double* zu, struct d_dense_qp* qp) {

    int ns = qp->dim->ns;
    int nv = qp->dim->nv;

    pack_vec(ns, zu, 1, qp->gz, nv + ns);
}


void d_dense_qp_set_lls(double* ls, struct d_dense_qp* qp) {

    int ns = qp->dim->ns;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;

    pack_vec(ns, ls, 1, qp->d, 2 * nb + 2 * ng);
    dvecse(ns, 0.0, qp->m, 2 * nb + 2 * ng);
}


void d_dense_qp_set_lls_mask(double* ls_mask, struct d_dense_qp* qp) {

    int nb = qp->dim->nb;
    int ng = qp->dim->ng;
    int ns = qp->dim->ns;

    int ii;

    for (ii = 0; ii < ns; ii++)
        if (ls_mask[ii] == 0.0)
            VECEL(qp->d_mask, 2 * nb + 2 * ng + ii) = 0.0;
        else
            VECEL(qp->d_mask, 2 * nb + 2 * ng + ii) = 1.0;
}


void d_dense_qp_set_lus(double* us, struct d_dense_qp* qp) {

    int ns = qp->dim->ns;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;

    pack_vec(ns, us, 1, qp->d, 2 * nb + 2 * ng + ns);
    dvecse(ns, 0.0, qp->m, 2 * nb + 2 * ng + ns);
}


void d_dense_qp_set_lus_mask(double* us_mask, struct d_dense_qp* qp) {

    int nb = qp->dim->nb;
    int ng = qp->dim->ng;
    int ns = qp->dim->ns;

    int ii;

    for (ii = 0; ii < ns; ii++)
        if (us_mask[ii] == 0.0)
            VECEL(qp->d_mask, 2 * nb + 2 * ng + ns + ii) = 0.0;
        else
            VECEL(qp->d_mask, 2 * nb + 2 * ng + ns + ii) = 1.0;
}


void d_dense_qp_get_H(struct d_dense_qp* qp, double* H) {

    int nv = qp->dim->nv;

    unpack_mat(nv, nv, qp->Hv, 0, 0, H, nv);
}


void d_dense_qp_get_g(struct d_dense_qp* qp, double* g) {

    int nv = qp->dim->nv;

    unpack_vec(nv, qp->gz, 0, g, 1);
}


void d_dense_qp_get_A(struct d_dense_qp* qp, double* A) {

    int nv = qp->dim->nv;
    int ne = qp->dim->ne;

    unpack_mat(ne, nv, qp->A, 0, 0, A, ne);
}


void d_dense_qp_get_b(struct d_dense_qp* qp, double* b) {

    int ne = qp->dim->ne;

    unpack_vec(ne, qp->b, 0, b, 1);
}


void d_dense_qp_get_idxb(struct d_dense_qp* qp, int* idxb) {

    int ii;
    int nb = qp->dim->nb;

    for (ii = 0; ii < nb; ii++) idxb[ii] = qp->idxb[ii];
}


void d_dense_qp_get_lb(struct d_dense_qp* qp, double* lb) {

    int nb = qp->dim->nb;

    unpack_vec(nb, qp->d, 0, lb, 1);
}


void d_dense_qp_get_lb_mask(struct d_dense_qp* qp, double* lb_mask) {

    int nb = qp->dim->nb;

    unpack_vec(nb, qp->d_mask, 0, lb_mask, 1);
}


void d_dense_qp_get_ub(struct d_dense_qp* qp, double* ub) {

    int ii;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;

    unpack_vec(nb, qp->d, nb + ng, ub, 1);
    for (ii = 0; ii < nb; ii++) ub[ii] = -ub[ii];
}


void d_dense_qp_get_ub_mask(struct d_dense_qp* qp, double* ub_mask) {

    int ii;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;

    unpack_vec(nb, qp->d_mask, nb + ng, ub_mask, 1);
}


void d_dense_qp_get_C(struct d_dense_qp* qp, double* C) {

    int nv = qp->dim->nv;
    int ng = qp->dim->ng;

    unpack_tran_mat(nv, ng, qp->Ct, 0, 0, C, ng);
}


void d_dense_qp_get_lg(struct d_dense_qp* qp, double* lg) {

    int nb = qp->dim->nb;
    int ng = qp->dim->ng;

    unpack_vec(ng, qp->d, nb, lg, 1);
}


void d_dense_qp_get_lg_mask(struct d_dense_qp* qp, double* lg_mask) {

    int nb = qp->dim->nb;
    int ng = qp->dim->ng;

    unpack_vec(ng, qp->d_mask, nb, lg_mask, 1);
}


void d_dense_qp_get_ug(struct d_dense_qp* qp, double* ug) {

    int ii;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;

    unpack_vec(ng, qp->d, 2 * nb + ng, ug, 1);
    for (ii = 0; ii < ng; ii++) ug[ii] = -ug[ii];
}


void d_dense_qp_get_ug_mask(struct d_dense_qp* qp, double* ug_mask) {

    int ii;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;

    unpack_vec(ng, qp->d_mask, 2 * nb + ng, ug_mask, 1);
}


void d_dense_qp_get_idxs(struct d_dense_qp* qp, int* idxs) {

    int ii, idx_tmp;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;
    int ns = qp->dim->ns;

    for (ii = 0; ii < nb + ng; ii++) {
        idx_tmp = qp->idxs_rev[ii];
        if (idx_tmp != -1) {
            idxs[idx_tmp] = ii;
        }
    }
}


void d_dense_qp_get_idxs_rev(struct d_dense_qp* qp, int* idxs_rev) {

    int ii;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;

    for (ii = 0; ii < nb + ng; ii++) {
        idxs_rev[ii] = qp->idxs_rev[ii];
    }
}


void d_dense_qp_get_Zl(struct d_dense_qp* qp, double* Zl) {

    int ns = qp->dim->ns;

    unpack_vec(ns, qp->Z, 0, Zl, 1);
}


void d_dense_qp_get_Zu(struct d_dense_qp* qp, double* Zu) {

    int ns = qp->dim->ns;

    unpack_vec(ns, qp->Z, ns, Zu, 1);
}


void d_dense_qp_get_zl(struct d_dense_qp* qp, double* zl) {

    int ns = qp->dim->ns;
    int nv = qp->dim->nv;

    unpack_vec(ns, qp->gz, nv, zl, 1);
}


void d_dense_qp_get_zu(struct d_dense_qp* qp, double* zu) {

    int ns = qp->dim->ns;
    int nv = qp->dim->nv;

    unpack_vec(ns, qp->gz, nv + ns, zu, 1);
}


void d_dense_qp_get_ls(struct d_dense_qp* qp, double* ls) {

    int ns = qp->dim->ns;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;

    unpack_vec(ns, qp->d, 2 * nb + 2 * ng, ls, 1);
}


void d_dense_qp_get_ls_mask(struct d_dense_qp* qp, double* ls_mask) {

    int ns = qp->dim->ns;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;

    unpack_vec(ns, qp->d_mask, 2 * nb + 2 * ng, ls_mask, 1);
}


void d_dense_qp_get_us(struct d_dense_qp* qp, double* us) {

    int ns = qp->dim->ns;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;

    unpack_vec(ns, qp->d, 2 * nb + 2 * ng + ns, us, 1);
}


void d_dense_qp_get_us_mask(struct d_dense_qp* qp, double* us_mask) {

    int ns = qp->dim->ns;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;

    unpack_vec(ns, qp->d_mask, 2 * nb + 2 * ng + ns, us_mask, 1);
}


void d_dense_qp_set_all_rowmaj(double* H, double* g, double* A, double* b, int* idxb, double* d_lb, double* d_ub, double* C, double* d_lg, double* d_ug, double* Zl, double* Zu, double* zl, double* zu, int* idxs, double* d_ls, double* d_us, struct d_dense_qp* qp) {

    int ii;

    int nv = qp->dim->nv;
    int ne = qp->dim->ne;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;
    int ns = qp->dim->ns;

    pack_tran_mat(nv, nv, H, nv, qp->Hv, 0, 0);
    pack_vec(nv, g, 1, qp->gz, 0);
    if (ne > 0) {
        pack_tran_mat(nv, ne, A, nv, qp->A, 0, 0);
        pack_vec(ne, b, 1, qp->b, 0);
    }
    if (nb > 0) {
        for (ii = 0; ii < nb; ii++) qp->idxb[ii] = idxb[ii];
        pack_vec(nb, d_lb, 1, qp->d, 0);
        pack_vec(nb, d_ub, 1, qp->d, nb + ng);
        dvecsc(nb, -1.0, qp->d, nb + ng);
        dvecse(nb, 0.0, qp->m, 0);
        dvecse(nb, 0.0, qp->m, nb + ng);
    }
    if (ng > 0) {
        pack_mat(nv, ng, C, nv, qp->Ct, 0, 0);
        pack_vec(ng, d_lg, 1, qp->d, nb);
        pack_vec(ng, d_ug, 1, qp->d, 2 * nb + ng);
        dvecsc(ng, -1.0, qp->d, 2 * nb + ng);
        dvecse(ng, 0.0, qp->m, nb);
        dvecse(ng, 0.0, qp->m, 2 * nb + ng);
    }
    if (ns > 0) {
        for (ii = 0; ii < ns; ii++) {
            qp->idxs_rev[idxs[ii]] = ii;
        }
        pack_vec(ns, Zl, 1, qp->Z, 0);
        pack_vec(ns, Zu, 1, qp->Z, ns);
        pack_vec(ns, zl, 1, qp->gz, nv);
        pack_vec(ns, zu, 1, qp->gz, nv + ns);
        pack_vec(ns, d_ls, 1, qp->d, 2 * nb + 2 * ng);
        pack_vec(ns, d_us, 1, qp->d, 2 * nb + 2 * ng + ns);
        dvecse(ns, 0.0, qp->m, 2 * nb + 2 * ng);
        dvecse(ns, 0.0, qp->m, 2 * nb + 2 * ng + ns);
    }
}


void d_dense_qp_get_all_rowmaj(struct d_dense_qp* qp, double* H, double* g, double* A, double* b, int* idxb, double* d_lb, double* d_ub, double* C, double* d_lg, double* d_ug, double* Zl, double* Zu, double* zl, double* zu, int* idxs, double* d_ls, double* d_us) {

    int ii, idx_tmp;

    int nv = qp->dim->nv;
    int ne = qp->dim->ne;
    int nb = qp->dim->nb;
    int ng = qp->dim->ng;
    int ns = qp->dim->ns;

    unpack_tran_mat(nv, nv, qp->Hv, 0, 0, H, nv);
    unpack_vec(nv, qp->gz, 0, g, 1);
    if (ne > 0) {
        unpack_tran_mat(ne, nv, qp->A, 0, 0, A, nv);
        unpack_vec(ne, qp->b, 0, b, 1);
    }
    if (nb > 0) {
        for (ii = 0; ii < nb; ii++) idxb[ii] = qp->idxb[ii];
        unpack_vec(nb, qp->d, 0, d_lb, 1);
        unpack_vec(nb, qp->d, nb + ng, d_ub, 1);
        for (ii = 0; ii < nb; ii++) d_ub[ii] = -d_ub[ii];
    }
    if (ng > 0) {
        unpack_mat(nv, ng, qp->Ct, 0, 0, C, nv);
        unpack_vec(ng, qp->d, nb, d_lg, 1);
        unpack_vec(ng, qp->d, 2 * nb + ng, d_ug, 1);
        for (ii = 0; ii < ng; ii++) d_ug[ii] = -d_ug[ii];
    }
    if (ns > 0) {
        // TODO only valid if there is one slack variable per soft constraint !!!
        for (ii = 0; ii < nb + ng; ii++) {
            idx_tmp = qp->idxs_rev[ii];
            if (idx_tmp != -1) {
                idxs[idx_tmp] = ii;
            }
        }
        unpack_vec(ns, qp->Z, 0, Zl, 1);
        unpack_vec(ns, qp->Z, ns, Zu, 1);
        unpack_vec(ns, qp->gz, nv, zl, 1);
        unpack_vec(ns, qp->gz, nv + ns, zu, 1);
        unpack_vec(ns, qp->d, 2 * nb + 2 * ng, d_ls, 1);
        unpack_vec(ns, qp->d, 2 * nb + 2 * ng + ns, d_us, 1);
    }
}
