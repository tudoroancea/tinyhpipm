#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_m_aux.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/dense/d_dense_qp_dim.h"
#include "hpipm/dense/s_dense_qp_dim.h"

void cvt_d2s_dense_qp_dim(struct d_dense_qp_dim* qpd, struct s_dense_qp_dim* qps) {

    qps->nv = qpd->nv;
    qps->ne = qpd->ne;
    qps->nb = qpd->nb;
    qps->ng = qpd->ng;
    qps->ns = qpd->ns;

    return;
}
