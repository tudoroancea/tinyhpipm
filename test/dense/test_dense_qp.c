#include <stdio.h>
#include <sys/_types/_size_t.h>

#include "tinyhpipm/common.h"
#include "tinyhpipm/dense/d_dense_qp.h"
#include "tinyhpipm/dense/d_dense_qp_dim.h"
#include "tinyhpipm/dense/d_dense_qp_ipm.h"
#include "tinyhpipm/dense/d_dense_qqp_sol.h"

int main() {
    /************************************************
     * dense qp dim
     ************************************************/

    size_t dim_size = d_dense_qp_dim_memsize();
    void* dim_mem = malloc(dim_size);

    struct d_dense_qp_dim dim;
    d_dense_qp_dim_create(&dim, dim_mem);

    // d_dense_qp_dim_set_all(nv, ne, nb, ng, nsb, nsg, &dim);
    d_dense_qp_dim_set_nv(nv, &dim);
    d_dense_qp_dim_set_ne(ne, &dim);
    d_dense_qp_dim_set_nb(nb, &dim);
    d_dense_qp_dim_set_ng(ng, &dim);
    d_dense_qp_dim_set_nsb(nsb, &dim);
    d_dense_qp_dim_set_nsg(nsg, &dim);

    //	d_dense_qp_dim_codegen("examples/c/data/test_d_dense_data.c", "w", &dim);
}
