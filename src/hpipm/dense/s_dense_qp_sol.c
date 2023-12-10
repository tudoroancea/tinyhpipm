#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_s_aux.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/auxiliary/aux_mem.h"
#include "hpipm/auxiliary/aux_string.h"
#include "hpipm/dense/s_dense_qp.h"
#include "hpipm/dense/s_dense_qp_dim.h"
#include "hpipm/dense/s_dense_qp_sol.h"



#define CREATE_STRVEC blasfeo_create_svec
#define UNPACK_VEC blasfeo_unpack_svec
#define PACK_VEC blasfeo_pack_svec
#define DENSE_QP s_dense_qp
#define DENSE_QP_DIM s_dense_qp_dim
#define DENSE_QP_SOL s_dense_qp_sol
#define REAL float
#define STRVEC blasfeo_svec
#define SIZE_STRVEC blasfeo_memsize_svec
#define VECCP_LIBSTR blasfeo_sveccp

#define DENSE_QP_SOL_STRSIZE s_dense_qp_sol_strsize
#define DENSE_QP_SOL_MEMSIZE s_dense_qp_sol_memsize
#define DENSE_QP_SOL_CREATE s_dense_qp_sol_create
#define DENSE_QP_SOL_GET_ALL s_dense_qp_sol_get_all
#define DENSE_QP_SOL_GET s_dense_qp_sol_get
#define DENSE_QP_SOL_GET_V s_dense_qp_sol_get_v
#define DENSE_QP_SOL_GET_PI s_dense_qp_sol_get_pi
#define DENSE_QP_SOL_GET_LAM_LB s_dense_qp_sol_get_lam_lb
#define DENSE_QP_SOL_GET_LAM_UB s_dense_qp_sol_get_lam_ub
#define DENSE_QP_SOL_GET_LAM_LG s_dense_qp_sol_get_lam_lg
#define DENSE_QP_SOL_GET_LAM_UG s_dense_qp_sol_get_lam_ug
#define DENSE_QP_SOL_GET_VALID_OBJ s_dense_qp_sol_get_valid_obj
#define DENSE_QP_SOL_GET_OBJ s_dense_qp_sol_get_obj
#define DENSE_QP_SOL_SET_V s_dense_qp_sol_set_v



#include "x_dense_qp_sol.c"


