#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_d_aux.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/auxiliary/aux_mem.h"
#include "hpipm/auxiliary/aux_string.h"
#include "hpipm/dense/d_dense_qcqp.h"
#include "hpipm/dense/d_dense_qcqp_dim.h"
#include "hpipm/dense/d_dense_qcqp_sol.h"


#define CREATE_STRVEC blasfeo_create_dvec
#define UNPACK_VEC blasfeo_unpack_dvec
#define PACK_VEC blasfeo_pack_dvec
#define DENSE_QCQP d_dense_qcqp
#define DENSE_QCQP_DIM d_dense_qcqp_dim
#define DENSE_QCQP_SOL d_dense_qcqp_sol
#define REAL double
#define STRVEC blasfeo_dvec
#define SIZE_STRVEC blasfeo_memsize_dvec
#define VECCP_LIBSTR blasfeo_dveccp

#define DENSE_QCQP_SOL_STRSIZE d_dense_qcqp_sol_strsize
#define DENSE_QCQP_SOL_MEMSIZE d_dense_qcqp_sol_memsize
#define DENSE_QCQP_SOL_CREATE d_dense_qcqp_sol_create
#define DENSE_QCQP_SOL_GET d_dense_qcqp_sol_get
#define DENSE_QCQP_SOL_GET_V d_dense_qcqp_sol_get_v
#define DENSE_QCQP_SOL_GET_PI d_dense_qcqp_sol_get_pi
#define DENSE_QCQP_SOL_GET_LAM_LB d_dense_qcqp_sol_get_lam_lb
#define DENSE_QCQP_SOL_GET_LAM_UB d_dense_qcqp_sol_get_lam_ub
#define DENSE_QCQP_SOL_GET_LAM_LG d_dense_qcqp_sol_get_lam_lg
#define DENSE_QCQP_SOL_GET_LAM_UG d_dense_qcqp_sol_get_lam_ug
#define DENSE_QCQP_SOL_GET_LAM_UQ d_dense_qcqp_sol_get_lam_uq
#define DENSE_QCQP_SOL_SET_V d_dense_qcqp_sol_set_v


#include "x_dense_qcqp_sol.c"
