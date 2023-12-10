#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_d_aux.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/auxiliary/aux_mem.h"
#include "hpipm/auxiliary/aux_string.h"
#include "hpipm/dense/d_dense_qp.h"
#include "hpipm/dense/d_dense_qp_dim.h"
#include "hpipm/dense/d_dense_qp_sol.h"


#define CREATE_STRVEC blasfeo_create_dvec
#define UNPACK_VEC blasfeo_unpack_dvec
#define PACK_VEC blasfeo_pack_dvec
#define DENSE_QP d_dense_qp
#define DENSE_QP_DIM d_dense_qp_dim
#define DENSE_QP_SOL d_dense_qp_sol
#define REAL double
#define STRVEC blasfeo_dvec
#define SIZE_STRVEC blasfeo_memsize_dvec
#define VECCP_LIBSTR blasfeo_dveccp

#define DENSE_QP_SOL_STRSIZE d_dense_qp_sol_strsize
#define DENSE_QP_SOL_MEMSIZE d_dense_qp_sol_memsize
#define DENSE_QP_SOL_CREATE d_dense_qp_sol_create
#define DENSE_QP_SOL_GET_ALL d_dense_qp_sol_get_all
#define DENSE_QP_SOL_GET d_dense_qp_sol_get
#define DENSE_QP_SOL_GET_V d_dense_qp_sol_get_v
#define DENSE_QP_SOL_GET_PI d_dense_qp_sol_get_pi
#define DENSE_QP_SOL_GET_LAM_LB d_dense_qp_sol_get_lam_lb
#define DENSE_QP_SOL_GET_LAM_UB d_dense_qp_sol_get_lam_ub
#define DENSE_QP_SOL_GET_LAM_LG d_dense_qp_sol_get_lam_lg
#define DENSE_QP_SOL_GET_LAM_UG d_dense_qp_sol_get_lam_ug
#define DENSE_QP_SOL_GET_VALID_OBJ d_dense_qp_sol_get_valid_obj
#define DENSE_QP_SOL_GET_OBJ d_dense_qp_sol_get_obj
#define DENSE_QP_SOL_SET_V d_dense_qp_sol_set_v


#include "x_dense_qp_sol.c"
