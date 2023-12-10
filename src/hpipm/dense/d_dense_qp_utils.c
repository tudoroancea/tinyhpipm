#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_d_aux.h"
#include "blasfeo/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/blasfeo_i_aux_ext_dep.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/dense/d_dense_qp.h"
#include "hpipm/dense/d_dense_qp_dim.h"
#include "hpipm/dense/d_dense_qp_ipm.h"
#include "hpipm/dense/d_dense_qp_res.h"
#include "hpipm/dense/d_dense_qp_sol.h"


#define DOUBLE_PRECISION


#define BLASFEO_PRINT_MAT blasfeo_print_dmat
#define BLASFEO_PRINT_TRAN_VEC blasfeo_print_tran_dvec
#define DENSE_QP d_dense_qp
#define DENSE_QP_IPM_ARG d_dense_qp_ipm_arg
#define DENSE_QP_DIM d_dense_qp_dim
#define DENSE_QP_IPM_ARG d_dense_qp_ipm_arg
#define DENSE_QP_RES d_dense_qp_res
#define DENSE_QP_SOL d_dense_qp_sol


#define DENSE_QP_DIM_PRINT d_dense_qp_dim_print
#define DENSE_QP_DIM_CODEGEN d_dense_qp_dim_codegen
#define DENSE_QP_PRINT d_dense_qp_print
#define DENSE_QP_CODEGEN d_dense_qp_codegen
#define DENSE_QP_SOL_PRINT d_dense_qp_sol_print
#define DENSE_QP_RES_PRINT d_dense_qp_res_print
#define DENSE_QP_IPM_ARG_PRINT d_dense_qp_ipm_arg_print
#define DENSE_QP_IPM_ARG_CODEGEN d_dense_qp_ipm_arg_codegen


#include "x_dense_qp_utils.c"
