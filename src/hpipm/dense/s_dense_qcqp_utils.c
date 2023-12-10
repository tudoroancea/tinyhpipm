#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_s_aux.h"
#include "blasfeo/blasfeo_s_aux_ext_dep.h"
#include "blasfeo/blasfeo_i_aux_ext_dep.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/dense/s_dense_qcqp.h"
#include "hpipm/dense/s_dense_qcqp_dim.h"
#include "hpipm/dense/s_dense_qcqp_res.h"
#include "hpipm/dense/s_dense_qcqp_sol.h"
// #include "hpipm_s_dense_qcqp_ipm.h"



#define SINGLE_PRECISION



#define BLASFEO_PRINT_MAT blasfeo_print_smat
#define BLASFEO_PRINT_TRAN_MAT blasfeo_print_tran_smat
#define BLASFEO_PRINT_TRAN_VEC blasfeo_print_tran_svec
#define DENSE_QCQP s_dense_qcqp
#define DENSE_QCQP_DIM s_dense_qcqp_dim
#define DENSE_QCQP_IPM_ARG s_dense_qcqp_ipm_arg
#define DENSE_QCQP_RES s_dense_qcqp_res
#define DENSE_QCQP_SOL s_dense_qcqp_sol



#define DENSE_QCQP_DIM_PRINT s_dense_qcqp_dim_print
#define DENSE_QCQP_PRINT s_dense_qcqp_print
#define DENSE_QCQP_SOL_PRINT s_dense_qcqp_sol_print
#define DENSE_QCQP_RES_PRINT s_dense_qcqp_res_print



#include "x_dense_qcqp_utils.c"


