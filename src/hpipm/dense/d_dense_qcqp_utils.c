#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_d_aux.h"
#include "blasfeo/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/blasfeo_i_aux_ext_dep.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/dense/d_dense_qcqp.h"
#include "hpipm/dense/d_dense_qcqp_dim.h"
#include "hpipm/dense/d_dense_qcqp_res.h"
#include "hpipm/dense/d_dense_qcqp_sol.h"
// #include "hpipm_d_dense_qcqp_ipm.h"


#define DOUBLE_PRECISION


#define BLASFEO_PRINT_MAT blasfeo_print_dmat
#define BLASFEO_PRINT_TRAN_MAT blasfeo_print_tran_dmat
#define BLASFEO_PRINT_TRAN_VEC blasfeo_print_tran_dvec
#define DENSE_QCQP d_dense_qcqp
#define DENSE_QCQP_DIM d_dense_qcqp_dim
#define DENSE_QCQP_IPM_ARG d_dense_qcqp_ipm_arg
#define DENSE_QCQP_RES d_dense_qcqp_res
#define DENSE_QCQP_SOL d_dense_qcqp_sol


#define DENSE_QCQP_DIM_PRINT d_dense_qcqp_dim_print
#define DENSE_QCQP_PRINT d_dense_qcqp_print
#define DENSE_QCQP_SOL_PRINT d_dense_qcqp_sol_print
#define DENSE_QCQP_RES_PRINT d_dense_qcqp_res_print


#include "x_dense_qcqp_utils.c"
