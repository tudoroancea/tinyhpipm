hpipm_size_t DENSE_QP_DIM_STRSIZE() {

    hpipm_size_t size = 0;

    size += sizeof(struct DENSE_QP_DIM);

    return size;
}


hpipm_size_t DENSE_QP_DIM_MEMSIZE() {

    hpipm_size_t size = 0;

    size = (size + 8 - 1) / 8 * 8;

    return size;
}


void DENSE_QP_DIM_CREATE(struct DENSE_QP_DIM* dim, void* mem) {

    // zero memory (to avoid corrupted memory like e.g. NaN)
    //	hpipm_size_t memsize = DENSE_QP_DIM_MEMSIZE();
    //	hpipm_zero_memset(memsize, mem);

    dim->memsize = DENSE_QP_DIM_MEMSIZE();

    // initialize dims to zero by default

    dim->nv = 0;
    dim->ne = 0;
    dim->nb = 0;
    dim->ng = 0;
    dim->ns = 0;
    dim->nsb = 0;
    dim->nsg = 0;

    return;
}


void DENSE_QP_DIM_SET_ALL(int nv, int ne, int nb, int ng, int nsb, int nsg, struct DENSE_QP_DIM* dim) {

    dim->nv = nv;
    dim->ne = ne;
    dim->nb = nb;
    dim->ng = ng;
    dim->ns = nsb + nsg;
    dim->nsb = nsb;
    dim->nsg = nsg;

    return;
}


void DENSE_QP_DIM_SET(char* field_name, int value, struct DENSE_QP_DIM* dim) {
    if (hpipm_strcmp(field_name, "nv")) {
        DENSE_QP_DIM_SET_NV(value, dim);
    } else if (hpipm_strcmp(field_name, "ne")) {
        DENSE_QP_DIM_SET_NE(value, dim);
    } else if (hpipm_strcmp(field_name, "nb")) {
        DENSE_QP_DIM_SET_NB(value, dim);
    } else if (hpipm_strcmp(field_name, "ng")) {
        DENSE_QP_DIM_SET_NG(value, dim);
    } else if (hpipm_strcmp(field_name, "nsb")) {
        DENSE_QP_DIM_SET_NSB(value, dim);
    } else if (hpipm_strcmp(field_name, "nsg")) {
        DENSE_QP_DIM_SET_NSG(value, dim);
    } else if (hpipm_strcmp(field_name, "ns")) {
        DENSE_QP_DIM_SET_NS(value, dim);
    } else {
        printf("error: SET_OCP_QP_DIM: wrong field %s\n", field_name);
        exit(1);
    }
    return;
}


void DENSE_QP_DIM_SET_NV(int value, struct DENSE_QP_DIM* dim) {
    dim->nv = value;

    return;
}


void DENSE_QP_DIM_SET_NE(int value, struct DENSE_QP_DIM* dim) {
    dim->ne = value;

    return;
}


void DENSE_QP_DIM_SET_NB(int value, struct DENSE_QP_DIM* dim) {
    dim->nb = value;

    return;
}


void DENSE_QP_DIM_SET_NG(int value, struct DENSE_QP_DIM* dim) {
    dim->ng = value;

    return;
}


void DENSE_QP_DIM_SET_NSB(int value, struct DENSE_QP_DIM* dim) {
    dim->nsb = value;
    dim->ns = dim->nsb + dim->nsg;

    return;
}


void DENSE_QP_DIM_SET_NSG(int value, struct DENSE_QP_DIM* dim) {
    dim->nsg = value;
    dim->ns = dim->nsb + dim->nsg;

    return;
}


void DENSE_QP_DIM_SET_NS(int value, struct DENSE_QP_DIM* dim) {
    dim->ns = value;

    return;
}


void DENSE_QP_DIM_GET_NV(struct DENSE_QP_DIM* dim, int* value) {
    *value = dim->nv;

    return;
}


void DENSE_QP_DIM_GET_NE(struct DENSE_QP_DIM* dim, int* value) {
    *value = dim->ne;

    return;
}


void DENSE_QP_DIM_GET_NB(struct DENSE_QP_DIM* dim, int* value) {
    *value = dim->nb;

    return;
}


void DENSE_QP_DIM_GET_NG(struct DENSE_QP_DIM* dim, int* value) {
    *value = dim->ng;

    return;
}
