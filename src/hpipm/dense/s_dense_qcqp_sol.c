#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_s_aux.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/auxiliary/aux_mem.h"
#include "hpipm/auxiliary/aux_string.h"
#include "hpipm/dense/s_dense_qcqp.h"
#include "hpipm/dense/s_dense_qcqp_dim.h"
#include "hpipm/dense/s_dense_qcqp_sol.h"


#define CREATE_STRVEC blasfeo_create_svec
#define UNPACK_VEC blasfeo_unpack_svec
#define PACK_VEC blasfeo_pack_svec
#define DENSE_QCQP s_dense_qcqp
#define DENSE_QCQP_DIM s_dense_qcqp_dim
#define DENSE_QCQP_SOL s_dense_qcqp_sol
#define REAL float
#define STRVEC blasfeo_svec
#define SIZE_STRVEC blasfeo_memsize_svec
#define VECCP_LIBSTR blasfeo_sveccp

#define DENSE_QCQP_SOL_STRSIZE s_dense_qcqp_sol_strsize
#define DENSE_QCQP_SOL_MEMSIZE s_dense_qcqp_sol_memsize
#define DENSE_QCQP_SOL_CREATE s_dense_qcqp_sol_create
#define DENSE_QCQP_SOL_GET s_dense_qcqp_sol_get
#define DENSE_QCQP_SOL_GET_V s_dense_qcqp_sol_get_v
#define DENSE_QCQP_SOL_GET_PI s_dense_qcqp_sol_get_pi
#define DENSE_QCQP_SOL_GET_LAM_LB s_dense_qcqp_sol_get_lam_lb
#define DENSE_QCQP_SOL_GET_LAM_UB s_dense_qcqp_sol_get_lam_ub
#define DENSE_QCQP_SOL_GET_LAM_LG s_dense_qcqp_sol_get_lam_lg
#define DENSE_QCQP_SOL_GET_LAM_UG s_dense_qcqp_sol_get_lam_ug
#define DENSE_QCQP_SOL_GET_LAM_UQ s_dense_qcqp_sol_get_lam_uq
#define DENSE_QCQP_SOL_SET_V s_dense_qcqp_sol_set_v



#include "x_dense_qcqp_sol.c"



