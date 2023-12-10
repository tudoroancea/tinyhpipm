#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_i_aux_ext_dep.h"
#include "blasfeo/blasfeo_s_aux.h"
#include "blasfeo/blasfeo_s_aux_ext_dep.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/ocp/s_ocp_qp.h"
#include "hpipm/ocp/s_ocp_qp_dim.h"
#include "hpipm/ocp/s_ocp_qp_ipm.h"
#include "hpipm/ocp/s_ocp_qp_sol.h"


#define SINGLE_PRECISION


#define BLASFEO_PRINT_MAT blasfeo_print_smat
#define BLASFEO_PRINT_TRAN_VEC blasfeo_print_tran_svec
#define OCP_QP s_ocp_qp
#define OCP_QP_SOL s_ocp_qp_sol
#define OCP_QP_IPM_ARG s_ocp_qp_ipm_arg
#define OCP_QP_RES s_ocp_qp_res
#define OCP_QP_DIM s_ocp_qp_dim


#define OCP_QP_DIM_PRINT s_ocp_qp_dim_print
#define OCP_QP_DIM_CODEGEN s_ocp_qp_dim_codegen
#define OCP_QP_PRINT s_ocp_qp_print
#define OCP_QP_CODEGEN s_ocp_qp_codegen
#define OCP_QP_SOL_PRINT s_ocp_qp_sol_print
#define OCP_QP_IPM_ARG_PRINT s_ocp_qp_ipm_arg_print
#define OCP_QP_IPM_ARG_CODEGEN s_ocp_qp_ipm_arg_codegen
#define OCP_QP_RES_PRINT s_ocp_qp_res_print


#include "x_ocp_qp_utils.c"
