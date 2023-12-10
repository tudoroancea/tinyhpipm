hpipm_size_t DENSE_QCQP_DIM_STRSIZE() {

    hpipm_size_t size = 0;

    size += sizeof(struct DENSE_QCQP_DIM);

    return size;
}


hpipm_size_t DENSE_QCQP_DIM_MEMSIZE() {

    hpipm_size_t size = 0;

    size += 1 * sizeof(struct DENSE_QP_DIM);
    size += 1 * DENSE_QP_DIM_MEMSIZE();

    size = (size + 63) / 64 * 64;  // make multiple of typical cache line size
    size += 1 * 64;  // align once to typical cache line size

    return size;
}


void DENSE_QCQP_DIM_CREATE(struct DENSE_QCQP_DIM* dim, void* mem) {

    // zero memory (to avoid corrupted memory like e.g. NaN)
    hpipm_size_t memsize = DENSE_QCQP_DIM_MEMSIZE();
    hpipm_zero_memset(memsize, mem);

    // qp_dim struct
    struct DENSE_QP_DIM* dim_ptr = mem;

    dim->qp_dim = dim_ptr;
    dim_ptr += 1;

    // align to typical cache line size
    hpipm_size_t s_ptr = (hpipm_size_t) dim_ptr;
    s_ptr = (s_ptr + 63) / 64 * 64;

    // void
    char* c_ptr = (char*) s_ptr;

    DENSE_QP_DIM_CREATE(dim->qp_dim, c_ptr);
    c_ptr += dim->qp_dim->memsize;


    dim->memsize = DENSE_QCQP_DIM_MEMSIZE();

#if defined(RUNTIME_CHECKS)
    if (c_ptr > ((char*) mem) + dim->memsize) {
        printf("\nerror: DENSE_QCQP_DIM_CREATE: outside memory bounds!\n\n");
        exit(1);
    }
#endif

    // initialize dims to zero by default

    dim->nv = 0;
    dim->ne = 0;
    dim->nb = 0;
    dim->ng = 0;
    dim->nq = 0;
    dim->ns = 0;
    dim->nsb = 0;
    dim->nsg = 0;
    dim->nsq = 0;

    return;
}


void DENSE_QCQP_DIM_SET(char* field_name, int value, struct DENSE_QCQP_DIM* dim) {
    if (hpipm_strcmp(field_name, "nv")) {
        DENSE_QCQP_DIM_SET_NV(value, dim);
    } else if (hpipm_strcmp(field_name, "ne")) {
        DENSE_QCQP_DIM_SET_NE(value, dim);
    } else if (hpipm_strcmp(field_name, "nb")) {
        DENSE_QCQP_DIM_SET_NB(value, dim);
    } else if (hpipm_strcmp(field_name, "ng")) {
        DENSE_QCQP_DIM_SET_NG(value, dim);
    } else if (hpipm_strcmp(field_name, "nq")) {
        DENSE_QCQP_DIM_SET_NQ(value, dim);
    } else if (hpipm_strcmp(field_name, "nsb")) {
        DENSE_QCQP_DIM_SET_NSB(value, dim);
    } else if (hpipm_strcmp(field_name, "nsg")) {
        DENSE_QCQP_DIM_SET_NSG(value, dim);
    } else if (hpipm_strcmp(field_name, "nsq")) {
        DENSE_QCQP_DIM_SET_NSQ(value, dim);
    } else if (hpipm_strcmp(field_name, "ns")) {
        DENSE_QCQP_DIM_SET_NS(value, dim);
    } else {
        printf("error: SET_OCP_QCQP_DIM: wrong field %s\n", field_name);
        exit(1);
    }
    return;
}


void DENSE_QCQP_DIM_SET_NV(int value, struct DENSE_QCQP_DIM* dim) {
    dim->nv = value;

    DENSE_QP_DIM_SET_NV(dim->nv, dim->qp_dim);

    return;
}


void DENSE_QCQP_DIM_SET_NE(int value, struct DENSE_QCQP_DIM* dim) {
    dim->ne = value;

    DENSE_QP_DIM_SET_NE(dim->ne, dim->qp_dim);

    return;
}


void DENSE_QCQP_DIM_SET_NB(int value, struct DENSE_QCQP_DIM* dim) {
    dim->nb = value;

    DENSE_QP_DIM_SET_NB(dim->nb, dim->qp_dim);

    return;
}


void DENSE_QCQP_DIM_SET_NG(int value, struct DENSE_QCQP_DIM* dim) {
    dim->ng = value;

    DENSE_QP_DIM_SET_NG(dim->ng + dim->nq, dim->qp_dim);

    return;
}


void DENSE_QCQP_DIM_SET_NQ(int value, struct DENSE_QCQP_DIM* dim) {
    dim->nq = value;

    DENSE_QP_DIM_SET_NG(dim->ng + dim->nq, dim->qp_dim);

    return;
}


void DENSE_QCQP_DIM_SET_NSB(int value, struct DENSE_QCQP_DIM* dim) {
    dim->nsb = value;
    dim->ns = dim->nsb + dim->nsg + dim->nsq;

    DENSE_QP_DIM_SET_NSB(dim->nsb, dim->qp_dim);
    DENSE_QP_DIM_SET_NS(dim->ns, dim->qp_dim);

    return;
}


void DENSE_QCQP_DIM_SET_NSG(int value, struct DENSE_QCQP_DIM* dim) {
    dim->nsg = value;
    dim->ns = dim->nsb + dim->nsg + dim->nsq;

    DENSE_QP_DIM_SET_NSG(dim->nsg + dim->nsq, dim->qp_dim);
    DENSE_QP_DIM_SET_NS(dim->ns, dim->qp_dim);

    return;
}


void DENSE_QCQP_DIM_SET_NSQ(int value, struct DENSE_QCQP_DIM* dim) {
    dim->nsq = value;
    dim->ns = dim->nsb + dim->nsg + dim->nsq;

    DENSE_QP_DIM_SET_NSG(dim->nsg + dim->nsq, dim->qp_dim);
    DENSE_QP_DIM_SET_NS(dim->ns, dim->qp_dim);

    return;
}


void DENSE_QCQP_DIM_SET_NS(int value, struct DENSE_QCQP_DIM* dim) {
    dim->ns = value;

    DENSE_QP_DIM_SET_NS(dim->ns, dim->qp_dim);

    return;
}

void DENSE_QCQP_DIM_GET_NV(struct DENSE_QCQP_DIM* dim, int* value) {
    *value = dim->nv;

    return;
}


void DENSE_QCQP_DIM_GET_NE(struct DENSE_QCQP_DIM* dim, int* value) {
    *value = dim->ne;

    return;
}


void DENSE_QCQP_DIM_GET_NB(struct DENSE_QCQP_DIM* dim, int* value) {
    *value = dim->nb;

    return;
}


void DENSE_QCQP_DIM_GET_NG(struct DENSE_QCQP_DIM* dim, int* value) {
    *value = dim->ng;

    return;
}


void DENSE_QCQP_DIM_GET_NQ(struct DENSE_QCQP_DIM* dim, int* value) {
    *value = dim->nq;

    return;
}
