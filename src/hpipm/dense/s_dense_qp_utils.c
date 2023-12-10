#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_s_aux.h"
#include "blasfeo/blasfeo_s_aux_ext_dep.h"
#include "blasfeo/blasfeo_i_aux_ext_dep.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/dense/s_dense_qp.h"
#include "hpipm/dense/s_dense_qp_dim.h"
#include "hpipm/dense/s_dense_qp_ipm.h"
#include "hpipm/dense/s_dense_qp_res.h"
#include "hpipm/dense/s_dense_qp_sol.h"


#define SINGLE_PRECISION



#define BLASFEO_PRINT_MAT blasfeo_print_smat
#define BLASFEO_PRINT_TRAN_VEC blasfeo_print_tran_svec
#define DENSE_QP s_dense_qp
#define DENSE_QP_IPM_ARG s_dense_qp_ipm_arg
#define DENSE_QP_DIM s_dense_qp_dim
#define DENSE_QP_IPM_ARG s_dense_qp_ipm_arg
#define DENSE_QP_RES s_dense_qp_res
#define DENSE_QP_SOL s_dense_qp_sol



#define DENSE_QP_DIM_PRINT s_dense_qp_dim_print
#define DENSE_QP_DIM_CODEGEN s_dense_qp_dim_codegen
#define DENSE_QP_PRINT s_dense_qp_print
#define DENSE_QP_CODEGEN s_dense_qp_codegen
#define DENSE_QP_SOL_PRINT s_dense_qp_sol_print
#define DENSE_QP_RES_PRINT s_dense_qp_res_print
#define DENSE_QP_IPM_ARG_PRINT s_dense_qp_ipm_arg_print
#define DENSE_QP_IPM_ARG_CODEGEN s_dense_qp_ipm_arg_codegen



#include "x_dense_qp_utils.c"

