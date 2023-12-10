#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_m_aux.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/dense/d_dense_qp.h"
#include "hpipm/dense/d_dense_qp_dim.h"
#include "hpipm/dense/s_dense_qp.h"
#include "hpipm/dense/s_dense_qp_dim.h"


void cvt_d2s_dense_qp(struct d_dense_qp* qpd, struct s_dense_qp* qps) {

    int ii;

    int nv = qpd->dim->nv;
    int ne = qpd->dim->ne;
    int nb = qpd->dim->nb;
    int ng = qpd->dim->ng;
    int ns = qpd->dim->ns;

    blasfeo_cvt_d2s_mat(nv, nv, qpd->Hv, 0, 0, qps->Hv, 0, 0);
    blasfeo_cvt_d2s_vec(2 * ns, qpd->Z, 0, qps->Z, 0);
    blasfeo_cvt_d2s_vec(nv + 2 * ns, qpd->gz, 0, qps->gz, 0);
    blasfeo_cvt_d2s_mat(ne, nv, qpd->A, 0, 0, qps->A, 0, 0);
    blasfeo_cvt_d2s_vec(ne, qpd->b, 0, qps->b, 0);
    blasfeo_cvt_d2s_mat(nv, ng, qpd->Ct, 0, 0, qps->Ct, 0, 0);
    blasfeo_cvt_d2s_vec(2 * nb + 2 * ng + 2 * ns, qpd->d, 0, qps->d, 0);
    blasfeo_cvt_d2s_vec(2 * nb + 2 * ng + 2 * ns, qpd->m, 0, qps->m, 0);
    for (ii = 0; ii < nb; ii++) qps->idxb[ii] = qpd->idxb[ii];
    for (ii = 0; ii < nb + ng; ii++) qps->idxs_rev[ii] = qpd->idxs_rev[ii];

    return;
}
