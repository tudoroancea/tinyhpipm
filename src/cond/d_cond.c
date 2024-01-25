#include <stdio.h>
#include <stdlib.h>

#include "tinyhpipm/blas.h"
#include "tinyhpipm/common.h"
#include "tinyhpipm/cond/d_cond.h"
#include "tinyhpipm/cond/d_cond_aux.h"
#include "tinyhpipm/dense/d_dense_qp.h"
#include "tinyhpipm/dense/d_dense_qp_sol.h"
#include "tinyhpipm/ocp/d_ocp_qp.h"
#include "tinyhpipm/ocp/d_ocp_qp_dim.h"
#include "tinyhpipm/ocp/d_ocp_qp_sol.h"


void d_cond_qp_compute_dim(struct d_ocp_qp_dim* ocp_dim, struct d_dense_qp_dim* dense_dim) {

    int N = ocp_dim->N;
    int* nx = ocp_dim->nx;
    int* nu = ocp_dim->nu;
    int* nbx = ocp_dim->nbx;
    int* nbu = ocp_dim->nbu;
    int* ng = ocp_dim->ng;
    int* ns = ocp_dim->ns;
    int* nsbx = ocp_dim->nsbx;
    int* nsbu = ocp_dim->nsbu;
    int* nsg = ocp_dim->nsg;

    int ii;

    int nvc = 0;
    int nec = 0;
    int nbc = 0;
    int ngc = 0;
    int nsc = 0;
    int nsbc = 0;
    int nsgc = 0;

    // first stage
    nvc += nx[0] + nu[0];
    nbc += nbx[0] + nbu[0];
    ngc += ng[0];
    nsc += ns[0];
    nsbc += nsbx[0] + nsbu[0];
    nsgc += nsg[0];
    // remaining stages
    for (ii = 1; ii <= N; ii++) {
        nvc += nu[ii];
        nbc += nbu[ii];
        ngc += nbx[ii] + ng[ii];
        nsc += ns[ii];
        nsbc += nsbu[ii];
        nsgc += nsbx[ii] + nsg[ii];
    }

    dense_dim->nv = nvc;
    dense_dim->ne = nec;
    dense_dim->nb = nbc;
    dense_dim->ng = ngc;
    dense_dim->ns = nsc;
    dense_dim->nsb = nsbc;
    dense_dim->nsg = nsgc;
}


hpipm_size_t d_cond_qp_arg_memsize() {

    hpipm_size_t size = 0;

    return size;
}


void d_cond_qp_arg_create(struct d_cond_qp_arg* cond_arg, void* mem) {

    cond_arg->memsize = d_cond_qp_arg_memsize();
}


void d_cond_qp_arg_set_default(struct d_cond_qp_arg* cond_arg) {

    cond_arg->cond_alg = 0;  // condensing algorithm
    cond_arg->cond_last_stage = 1;  // condense last stage
    cond_arg->comp_prim_sol = 1;  // compute primal solution (v)
    cond_arg->comp_dual_sol_eq = 1;  // compute dual solution equality constr (pi)
    cond_arg->comp_dual_sol_ineq = 1;  // compute dual solution inequality constr (lam t)
    cond_arg->square_root_alg = 1;  // square root algorithm (faster but requires RSQ>0)
}


void d_cond_qp_arg_set_cond_alg(int cond_alg, struct d_cond_qp_arg* cond_arg) {

    cond_arg->cond_alg = cond_alg;
}


void d_cond_qp_arg_set_ric_alg(int ric_alg, struct d_cond_qp_arg* cond_arg) {

    cond_arg->square_root_alg = ric_alg;
}


void d_cond_qp_arg_set_cond_last_stage(int cond_last_stage, struct d_cond_qp_arg* cond_arg) {

    cond_arg->cond_last_stage = cond_last_stage;
}


void d_cond_qp_arg_set_comp_prim_sol(int value, struct d_cond_qp_arg* cond_arg) {

    cond_arg->comp_prim_sol = value;
}


void d_cond_qp_arg_set_comp_dual_sol_eq(int value, struct d_cond_qp_arg* cond_arg) {

    cond_arg->comp_dual_sol_eq = value;
}


void d_cond_qp_arg_set_comp_dual_sol_ineq(int value, struct d_cond_qp_arg* cond_arg) {

    cond_arg->comp_dual_sol_ineq = value;
}


hpipm_size_t d_cond_qp_ws_memsize(struct d_ocp_qp_dim* ocp_dim, struct d_cond_qp_arg* cond_arg) {

    int ii;
    int nu_tmp;

    int N = ocp_dim->N;
    int* nx = ocp_dim->nx;
    int* nu = ocp_dim->nu;
    int* nb = ocp_dim->nb;
    int* ng = ocp_dim->ng;

    // compute core qp size and max size
    int nvt = 0;
    int net = 0;
    int nbt = 0;
    int ngt = 0;
    int nxM = 0;
    int nuM = 0;
    int nbM = 0;
    int ngM = 0;

    for (ii = 0; ii < N; ii++) {
        nvt += nx[ii] + nu[ii];
        net += nx[ii + 1];
        nbt += nb[ii];
        ngt += ng[ii];
        nxM = nx[ii] > nxM ? nx[ii] : nxM;
        nuM = nu[ii] > nuM ? nu[ii] : nuM;
        nbM = nb[ii] > nbM ? nb[ii] : nbM;
        ngM = ng[ii] > ngM ? ng[ii] : ngM;
    }
    ii = N;
    nvt += nx[ii] + nu[ii];
    nbt += nb[ii];
    ngt += ng[ii];
    nxM = nx[ii] > nxM ? nx[ii] : nxM;
    nuM = nu[ii] > nuM ? nu[ii] : nuM;
    nbM = nb[ii] > nbM ? nb[ii] : nbM;
    ngM = ng[ii] > ngM ? ng[ii] : ngM;

    hpipm_size_t size = 0;

    size += 1 * N * sizeof(struct mat);  // Gamma
    size += 1 * N * sizeof(struct vec);  // Gammab
    size += 2 * sizeof(struct vec);  // tmp_nbgM tmp_nuxM
    if (cond_arg->cond_alg == 0) {
        size += 1 * (N + 1) * sizeof(struct mat);  // L
        size += 2 * sizeof(struct mat);  // Lx AL
        size += 1 * (N + 1) * sizeof(struct vec);  // l
    } else {
        size += 1 * sizeof(struct mat);  // GammaQ
    }

    nu_tmp = 0;
    for (ii = 0; ii < N; ii++) {
        nu_tmp += nu[ii];
        size += memsize_mat(nu_tmp + nx[0] + 1, nx[ii + 1]);  // Gamma
    }
    if (cond_arg->cond_alg == 0) {
        for (ii = 0; ii <= N; ii++)
            size += memsize_mat(nu[ii] + nx[ii] + 1, nu[ii] + nx[ii]);  // L
        size += memsize_mat(nxM + 1, nxM);  // Lx
        size += memsize_mat(nuM + nxM + 1, nxM);  // AL
        for (ii = 0; ii <= N; ii++)
            size += memsize_vec(nu[ii] + nx[ii]);  // l
    } else {
        nu_tmp = 0;
        for (ii = 0; ii < N; ii++)
            nu_tmp += nu[ii];
        size += memsize_mat(nu_tmp + nxM + 1, nxM);  // GammaQ
    }

    for (ii = 0; ii < N; ii++)
        size += 1 * memsize_vec(nx[ii + 1]);  // Gammab
    size += 1 * memsize_vec(nbM + ngM);  // tmp_nbgM
    size += 1 * memsize_vec(nuM + nxM);  // tmp_nuxM

    size += 1 * 64;  // align to typical cache line size
    size += 1 * 8;  // align to 64 bits

    size = (size + 63) / 64 * 64;  // make multiple of typical cache line size

    return size;
}


void d_cond_qp_ws_create(struct d_ocp_qp_dim* ocp_dim, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws, void* mem) {

    int ii;
    int nu_tmp;

    size_t s_ptr;

    // zero memory (to avoid corrupted memory like e.g. NaN)
    hpipm_size_t memsize = d_cond_qp_ws_memsize(ocp_dim, cond_arg);
    hpipm_zero_memset(memsize, mem);

    int N = ocp_dim->N;
    int* nx = ocp_dim->nx;
    int* nu = ocp_dim->nu;
    int* nb = ocp_dim->nb;
    int* ng = ocp_dim->ng;

    // compute core qp dim and max dim
    int nvt = 0;
    int net = 0;
    int nbt = 0;
    int ngt = 0;
    int nxM = 0;
    int nuM = 0;
    int nbM = 0;
    int ngM = 0;

    for (ii = 0; ii < N; ii++) {
        nvt += nx[ii] + nu[ii];
        net += nx[ii + 1];
        nbt += nb[ii];
        ngt += ng[ii];
        nxM = nx[ii] > nxM ? nx[ii] : nxM;
        nuM = nu[ii] > nuM ? nu[ii] : nuM;
        nbM = nb[ii] > nbM ? nb[ii] : nbM;
        ngM = ng[ii] > ngM ? ng[ii] : ngM;
    }
    ii = N;
    nvt += nx[ii] + nu[ii];
    nbt += nb[ii];
    ngt += ng[ii];
    nxM = nx[ii] > nxM ? nx[ii] : nxM;
    nuM = nu[ii] > nuM ? nu[ii] : nuM;
    nbM = nb[ii] > nbM ? nb[ii] : nbM;
    ngM = ng[ii] > ngM ? ng[ii] : ngM;


    // align to typical cache line size
    s_ptr = (size_t) mem;
    s_ptr = (s_ptr + 7) / 8 * 8;


    // matrix struct
    struct mat* sm_ptr = (struct mat*) s_ptr;

    cond_ws->Gamma = sm_ptr;
    sm_ptr += N;
    if (cond_arg->cond_alg == 0) {
        cond_ws->L = sm_ptr;
        sm_ptr += N + 1;
        cond_ws->Lx = sm_ptr;
        sm_ptr += 1;
        cond_ws->AL = sm_ptr;
        sm_ptr += 1;
    } else {
        cond_ws->GammaQ = sm_ptr;
        sm_ptr += 1;
    }


    // vector struct
    struct vec* sv_ptr = (struct vec*) sm_ptr;

    cond_ws->Gammab = sv_ptr;
    sv_ptr += N;
    cond_ws->tmp_nbgM = sv_ptr;
    sv_ptr += 1;
    cond_ws->tmp_nuxM = sv_ptr;
    sv_ptr += 1;
    if (cond_arg->cond_alg == 0) {
        cond_ws->l = sv_ptr;
        sv_ptr += N + 1;
    }


    // align to typical cache line size
    s_ptr = (size_t) sv_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;


    // void stuf
    char* c_ptr = (char*) s_ptr;
    char* c_tmp;

    nu_tmp = 0;
    for (ii = 0; ii < N; ii++) {
        nu_tmp += nu[ii];
        create_mat(nu_tmp + nx[0] + 1, nx[ii + 1], cond_ws->Gamma + ii, c_ptr);
        c_ptr += (cond_ws->Gamma + ii)->memsize;
    }
    if (cond_arg->cond_alg == 0) {
        for (ii = 0; ii <= N; ii++) {
            create_mat(nu[ii] + nx[ii] + 1, nu[ii] + nx[ii], cond_ws->L + ii, c_ptr);
            c_ptr += (cond_ws->L + ii)->memsize;
        }
        create_mat(nxM + 1, nxM, cond_ws->Lx, c_ptr);
        c_ptr += cond_ws->Lx->memsize;
        create_mat(nuM + nxM + 1, nxM, cond_ws->AL, c_ptr);
        c_ptr += cond_ws->AL->memsize;
        for (ii = 0; ii <= N; ii++) {
            create_vec(nu[ii] + nx[ii], cond_ws->l + ii, c_ptr);
            c_ptr += (cond_ws->l + ii)->memsize;
        }
    } else {
        nu_tmp = 0;
        for (ii = 0; ii < N; ii++)
            nu_tmp += nu[ii];
        create_mat(nu_tmp + nxM + 1, nxM, cond_ws->GammaQ, c_ptr);
        c_ptr += cond_ws->GammaQ->memsize;
    }
    for (ii = 0; ii < N; ii++) {
        create_vec(nx[ii + 1], cond_ws->Gammab + ii, c_ptr);
        c_ptr += (cond_ws->Gammab + ii)->memsize;
    }
    create_vec(nbM + ngM, cond_ws->tmp_nbgM, c_ptr);
    c_ptr += cond_ws->tmp_nbgM->memsize;
    c_tmp = c_ptr;
    create_vec(nuM + nxM, cond_ws->tmp_nuxM, c_ptr);
    c_ptr += cond_ws->tmp_nuxM->memsize;

    cond_ws->bs = N;

    cond_ws->memsize = d_cond_qp_ws_memsize(ocp_dim, cond_arg);

#if defined(RUNTIME_CHECKS)
    if (c_ptr > ((char*) mem) + cond_ws->memsize) {
        printf("\nCreate_cond_qp_ocp2dense: outsize memory bounds!\n\n");
        exit(1);
    }
#endif
}


void d_cond_qp_cond(struct d_ocp_qp* ocp_qp, struct d_dense_qp* dense_qp, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws) {

    d_cond_BAbt(ocp_qp, NULL, NULL, cond_arg, cond_ws);

    d_cond_RSQrq(ocp_qp, dense_qp->Hv, dense_qp->gz, cond_arg, cond_ws);

    d_cond_DCtd(ocp_qp, dense_qp->idxb, dense_qp->Ct, dense_qp->d, dense_qp->d_mask, dense_qp->idxs_rev, dense_qp->Z, dense_qp->gz, cond_arg, cond_ws);
}


void d_cond_qp_cond_lhs(struct d_ocp_qp* ocp_qp, struct d_dense_qp* dense_qp, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws) {

    d_cond_BAt(ocp_qp, NULL, cond_arg, cond_ws);

    d_cond_RSQ(ocp_qp, dense_qp->Hv, cond_arg, cond_ws);

    d_cond_DCt(ocp_qp, dense_qp->idxb, dense_qp->Ct, dense_qp->idxs_rev, dense_qp->Z, cond_arg, cond_ws);
}


void d_cond_qp_cond_rhs(struct d_ocp_qp* ocp_qp, struct d_dense_qp* dense_qp, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws) {

    d_cond_b(ocp_qp, NULL, cond_arg, cond_ws);

    d_cond_rq(ocp_qp, dense_qp->gz, cond_arg, cond_ws);

    d_cond_d(ocp_qp, dense_qp->d, dense_qp->d_mask, dense_qp->gz, cond_arg, cond_ws);
}


void d_cond_qp_expand_sol(struct d_ocp_qp* ocp_qp, struct d_dense_qp_sol* dense_qp_sol, struct d_ocp_qp_sol* ocp_qp_sol, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws) {

    d_expand_sol(ocp_qp, dense_qp_sol, ocp_qp_sol, cond_arg, cond_ws);
}


// TODO remove
void d_cond_qp_expand_primal_sol(struct d_ocp_qp* ocp_qp, struct d_dense_qp_sol* dense_qp_sol, struct d_ocp_qp_sol* ocp_qp_sol, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws) {

    d_expand_primal_sol(ocp_qp, dense_qp_sol, ocp_qp_sol, cond_arg, cond_ws);
}


/************************************************
 * update cond
 ************************************************/

void d_cond_qp_update(int* idxc, struct d_ocp_qp* ocp_qp, struct d_dense_qp* dense_qp, struct d_cond_qp_arg* cond_arg, struct d_cond_qp_ws* cond_ws) {

    d_update_cond_BAbt(idxc, ocp_qp, NULL, NULL, cond_arg, cond_ws);

    d_update_cond_RSQrq_N2nx3(idxc, ocp_qp, dense_qp->Hv, dense_qp->gz, cond_arg, cond_ws);

    d_update_cond_DCtd(idxc, ocp_qp, dense_qp->idxb, dense_qp->Ct, dense_qp->d, dense_qp->idxs_rev, dense_qp->Z, dense_qp->gz, cond_arg, cond_ws);
}
