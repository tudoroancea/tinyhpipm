#if defined(RUNTIME_CHECKS)
#include <stdlib.h>
#endif

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_m_aux.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/ocp/d_ocp_qp.h"
#include "hpipm/ocp/s_ocp_qp.h"


void m_cvt_d_ocp_qp_to_s_ocp_qp(struct d_ocp_qp* d_qp, struct s_ocp_qp* s_qp) {

    // TODO check that they have equal sizes !!!!!

    int N = d_qp->N;
    int* nx = d_qp->nx;
    int* nu = d_qp->nu;
    int* nb = d_qp->nb;
    int* ng = d_qp->ng;

    int ii, jj;

    for (ii = 0; ii < N; ii++) {
        m_cvt_d2blasfeo_smat(nu[ii] + nx[ii] + 1, nx[ii + 1], d_qp->BAbt + ii, 0, 0, s_qp->BAbt + ii, 0, 0);
        m_cvt_d2blasfeo_smat(nu[ii] + nx[ii] + 1, nu[ii] + nx[ii], d_qp->RSQrq + ii, 0, 0, s_qp->RSQrq + ii, 0, 0);
        m_cvt_d2blasfeo_smat(nu[ii] + nx[ii], ng[ii], d_qp->DCt + ii, 0, 0, s_qp->DCt + ii, 0, 0);
        m_cvt_d2blasfeo_svec(nx[ii + 1], d_qp->b + ii, 0, s_qp->b + ii, 0);
        m_cvt_d2blasfeo_svec(nu[ii] + nx[ii], d_qp->rq + ii, 0, s_qp->rq + ii, 0);
        m_cvt_d2blasfeo_svec(nb[ii], d_qp->d_lb + ii, 0, s_qp->d_lb + ii, 0);
        m_cvt_d2blasfeo_svec(nb[ii], d_qp->d_ub + ii, 0, s_qp->d_ub + ii, 0);
        m_cvt_d2blasfeo_svec(ng[ii], d_qp->d_lg + ii, 0, s_qp->d_lg + ii, 0);
        m_cvt_d2blasfeo_svec(ng[ii], d_qp->d_ug + ii, 0, s_qp->d_ug + ii, 0);
        for (jj = 0; jj < nb[ii]; jj++) s_qp->idxb[ii][jj] = d_qp->idxb[ii][jj];
    }
    ii = N;
    m_cvt_d2blasfeo_smat(nu[ii] + nx[ii] + 1, nu[ii] + nx[ii], d_qp->RSQrq + ii, 0, 0, s_qp->RSQrq + ii, 0, 0);
    m_cvt_d2blasfeo_smat(nu[ii] + nx[ii], ng[ii], d_qp->DCt + ii, 0, 0, s_qp->DCt + ii, 0, 0);
    m_cvt_d2blasfeo_svec(nu[ii] + nx[ii], d_qp->rq + ii, 0, s_qp->rq + ii, 0);
    m_cvt_d2blasfeo_svec(nb[ii], d_qp->d_lb + ii, 0, s_qp->d_lb + ii, 0);
    m_cvt_d2blasfeo_svec(nb[ii], d_qp->d_ub + ii, 0, s_qp->d_ub + ii, 0);
    m_cvt_d2blasfeo_svec(ng[ii], d_qp->d_lg + ii, 0, s_qp->d_lg + ii, 0);
    m_cvt_d2blasfeo_svec(ng[ii], d_qp->d_ug + ii, 0, s_qp->d_ug + ii, 0);
    for (jj = 0; jj < nb[ii]; jj++) s_qp->idxb[ii][jj] = d_qp->idxb[ii][jj];

    return;
}
