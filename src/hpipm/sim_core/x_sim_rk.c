hpipm_size_t SIM_RK_DATA_MEMSIZE(int ns) {

    hpipm_size_t size = 0;

    size += 1 * ns * ns * sizeof(REAL);  // A
    size += 2 * ns * sizeof(REAL);  // B C

    return size;
}


void SIM_RK_DATA_CREATE(int ns, struct SIM_RK_DATA* rk_data, void* mem) {

    rk_data->ns = ns;

    REAL* d_ptr = mem;

    //
    rk_data->A_rk = d_ptr;
    d_ptr += ns * ns;
    //
    rk_data->B_rk = d_ptr;
    d_ptr += ns;
    //
    rk_data->C_rk = d_ptr;
    d_ptr += ns;


    rk_data->memsize = SIM_RK_DATA_MEMSIZE(ns);


    char* c_ptr = (char*) d_ptr;


#if defined(RUNTIME_CHECKS)
    if (c_ptr > ((char*) mem) + rk_data->memsize) {
        printf("\nSIM_RK_DATA_CREATE: outsize memory bounds!\n\n");
        exit(1);
    }
#endif


    return;
}


void SIM_RK_DATA_INIT_DEFAULT(char* field, struct SIM_RK_DATA* rk_data) {
    if (hpipm_strcmp(field, "ERK4")) {
        if (rk_data->ns == 4) {
            int ns = 4;
            REAL* A = rk_data->A_rk;
            REAL* B = rk_data->B_rk;
            REAL* C = rk_data->C_rk;
            rk_data->expl = 1;
            // classic ERK4
            // A
            A[0 + ns * 0] = 0.0;
            A[0 + ns * 1] = 0.0;
            A[0 + ns * 2] = 0.0;
            A[0 + ns * 3] = 0.0;
            A[1 + ns * 0] = 0.5;
            A[1 + ns * 1] = 0.0;
            A[1 + ns * 2] = 0.0;
            A[1 + ns * 3] = 0.0;
            A[2 + ns * 0] = 0.0;
            A[2 + ns * 1] = 0.5;
            A[2 + ns * 2] = 0.0;
            A[2 + ns * 3] = 0.0;
            A[3 + ns * 0] = 0.0;
            A[3 + ns * 1] = 0.0;
            A[3 + ns * 2] = 1.0;
            A[3 + ns * 3] = 0.0;
            // b
            B[0] = 1.0 / 6.0;
            B[1] = 1.0 / 3.0;
            B[2] = 1.0 / 3.0;
            B[3] = 1.0 / 6.0;
            // c
            C[0] = 0.0;
            C[1] = 0.5;
            C[2] = 0.5;
            C[3] = 1.0;
        } else {
            printf("error: SIM_RK_DATA_INIT_DEFAULT: rk_data->ns=%d != 4. Exiting.\n", rk_data->ns);
            exit(1);
        }
    } else if (hpipm_strcmp(field, "ERK2")) {
        if (rk_data->ns == 2) {
            int ns = 2;
            REAL* A = rk_data->A_rk;
            REAL* B = rk_data->B_rk;
            REAL* C = rk_data->C_rk;
            rk_data->expl = 1;
            // explicit midpoint method
            // A
            A[0 + ns * 0] = 0.0;
            A[0 + ns * 1] = 0.0;
            A[1 + ns * 0] = 0.5;
            A[1 + ns * 1] = 0.0;
            // b
            B[0] = 0.0;
            B[1] = 1.0;
            // c
            C[0] = 0.0;
            C[1] = 0.5;
        } else {
            printf("error: SIM_RK_DATA_INIT_DEFAULT: rk_data->ns=%d != 4. Exiting.\n", rk_data->ns);
            exit(1);
        }
    } else {
        printf("error: SIM_RK_DATA_INIT_DEFAULT: wrong field name '%s'. Exiting.\n", field);
        exit(1);
    }
    return;
}


void SIM_RK_DATA_SET_ALL(int expl, REAL* A, REAL* B, REAL* C, struct SIM_RK_DATA* rk_data) {

    int ii, jj;

    rk_data->expl = expl;

    int ns = rk_data->ns;
    REAL* A_rk = rk_data->A_rk;
    REAL* B_rk = rk_data->B_rk;
    REAL* C_rk = rk_data->C_rk;

    for (jj = 0; jj < ns; jj++)
        for (ii = 0; ii < ns; ii++)
            A_rk[ii + ns * jj] = A[ii + ns * jj];
    for (ii = 0; ii < ns; ii++)
        B_rk[ii] = B[ii];
    for (ii = 0; ii < ns; ii++)
        C_rk[ii] = C[ii];

    return;
}
