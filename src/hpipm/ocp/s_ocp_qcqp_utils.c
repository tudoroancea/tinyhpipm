#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_i_aux_ext_dep.h"
#include "blasfeo/blasfeo_s_aux.h"
#include "blasfeo/blasfeo_s_aux_ext_dep.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/ocp/s_ocp_qcqp_dim.h"
#include "hpipm/ocp/s_ocp_qcqp_ipm.h"
#include "hpipm/ocp/s_ocp_qcqp_sol.h"
#include "hpipm/ocp/s_ocp_qp.h"


#define SINGLE_PRECISION


#define BLASFEO_PRINT_MAT blasfeo_print_smat
#define BLASFEO_PRINT_TRAN_MAT blasfeo_print_tran_smat
#define BLASFEO_PRINT_TRAN_VEC blasfeo_print_tran_svec
#define BLASFEO_PRINT_EXP_TRAN_VEC blasfeo_print_exp_tran_svec
#define OCP_QCQP s_ocp_qcqp
#define OCP_QCQP_DIM s_ocp_qcqp_dim
#define OCP_QCQP_IPM_ARG s_ocp_qcqp_ipm_arg
#define OCP_QCQP_RES s_ocp_qcqp_res
#define OCP_QCQP_SOL s_ocp_qcqp_sol


#define OCP_QCQP_DIM_PRINT s_ocp_qcqp_dim_print
#define OCP_QCQP_DIM_CODEGEN s_ocp_qcqp_dim_codegen
#define OCP_QCQP_PRINT s_ocp_qcqp_print
#define OCP_QCQP_CODEGEN s_ocp_qcqp_codegen
#define OCP_QCQP_SOL_PRINT s_ocp_qcqp_sol_print
#define OCP_QCQP_IPM_ARG_CODEGEN s_ocp_qcqp_ipm_arg_codegen
#define OCP_QCQP_RES_PRINT s_ocp_qcqp_res_print


#include "x_ocp_qcqp_utils.c"
