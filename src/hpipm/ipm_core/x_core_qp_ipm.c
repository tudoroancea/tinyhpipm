hpipm_size_t MEMSIZE_CORE_QP_IPM(int nv, int ne, int nc) {

    hpipm_size_t size;

    int nv0 = nv;
    int ne0 = ne;
    int nc0 = nc;
    // if target avx
    // nv0 = ...

    size = 0;

    size += 3 * nv0 * sizeof(REAL);  // v_bkp dv res_g
    size += 3 * ne0 * sizeof(REAL);  // pi_bkp dpi res_b
    size += 10 * nc0 * sizeof(REAL);  // lam_bkp t_bkp dlam dt res_d res_m res_m_bkp t_inv Gamma gamma

    size = (size + 63) / 64 * 64;  // make multiple of cache line size

    return size;
}


void CREATE_CORE_QP_IPM(int nv, int ne, int nc, struct CORE_QP_IPM_WORKSPACE* workspace, void* mem) {

    workspace->nv = nv;
    workspace->ne = ne;
    workspace->nc = nc;

    int nv0 = nv;
    int ne0 = ne;
    int nc0 = nc;
    // if target avx NO!!!!
    // nv0 = ...

    REAL* d_ptr = (REAL*) mem;

    workspace->t_inv = d_ptr;  // t_inv
    d_ptr += nc0;

    workspace->v_bkp = d_ptr;  // v_bkp
    d_ptr += nv0;

    workspace->pi_bkp = d_ptr;  // pi_bkp
    d_ptr += ne0;

    workspace->lam_bkp = d_ptr;  // lam_bkp
    d_ptr += nc0;

    workspace->t_bkp = d_ptr;  // t_bkp
    d_ptr += nc0;

    workspace->dv = d_ptr;  // dv
    d_ptr += nv0;

    workspace->dpi = d_ptr;  // dpi
    d_ptr += ne0;

    workspace->dlam = d_ptr;  // dlam
    d_ptr += nc0;

    workspace->dt = d_ptr;  // dt
    d_ptr += nc0;

    workspace->res_g = d_ptr;  // res_g
    d_ptr += nv0;

    workspace->res_b = d_ptr;  // res_b
    d_ptr += ne0;

    workspace->res_d = d_ptr;  // res_d
    d_ptr += nc0;

    workspace->res_m = d_ptr;  // res_m
    d_ptr += nc0;

    workspace->res_m_bkp = d_ptr;  // res_m_bkp
    d_ptr += nc0;

    workspace->Gamma = d_ptr;  // Gamma
    d_ptr += nc0;

    workspace->gamma = d_ptr;  // gamma
    d_ptr += nc0;

    // by default no constr masking
    workspace->nc_mask = nc;

    if (nc > 0) {
        workspace->nc_inv = 1.0 / nc;
        workspace->nc_mask_inv = workspace->nc_inv;
    } else {
        workspace->nc_inv = 0.0;
        workspace->nc_mask_inv = 0.0;
    }


    workspace->lam_min = 0.0;
    workspace->t_min = 0.0;
    workspace->t_min_inv = 1e30;
    workspace->tau_min = 0.0;
    workspace->split_step = 0;
    workspace->t_lam_min = 2;


    workspace->memsize = MEMSIZE_CORE_QP_IPM(nv, ne, nc);


    char* c_ptr = (char*) d_ptr;


#if defined(RUNTIME_CHECKS)
    if (c_ptr > ((char*) mem) + workspace->memsize) {
        printf("\nCreate_core_qp_ipm: outsize memory bounds!\n\n");
        exit(1);
    }
#endif


    return;
}
