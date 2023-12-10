#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_d_aux.h"
#include "blasfeo/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/blasfeo_i_aux_ext_dep.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/ocp/d_ocp_qp.h"
#include "hpipm/ocp/d_ocp_qp_dim.h"
#include "hpipm/ocp/d_ocp_qp_ipm.h"
#include "hpipm/ocp/d_ocp_qp_sol.h"


#define DOUBLE_PRECISION


#define BLASFEO_PRINT_MAT blasfeo_print_dmat
#define BLASFEO_PRINT_TRAN_VEC blasfeo_print_tran_dvec
#define OCP_QP d_ocp_qp
#define OCP_QP_DIM d_ocp_qp_dim
#define OCP_QP_IPM_ARG d_ocp_qp_ipm_arg
#define OCP_QP_RES d_ocp_qp_res
#define OCP_QP_SOL d_ocp_qp_sol


#define OCP_QP_DIM_PRINT d_ocp_qp_dim_print
#define OCP_QP_DIM_CODEGEN d_ocp_qp_dim_codegen
#define OCP_QP_PRINT d_ocp_qp_print
#define OCP_QP_CODEGEN d_ocp_qp_codegen
#define OCP_QP_SOL_PRINT d_ocp_qp_sol_print
#define OCP_QP_IPM_ARG_PRINT d_ocp_qp_ipm_arg_print
#define OCP_QP_IPM_ARG_CODEGEN d_ocp_qp_ipm_arg_codegen
#define OCP_QP_RES_PRINT d_ocp_qp_res_print


#include "x_ocp_qp_utils.c"
