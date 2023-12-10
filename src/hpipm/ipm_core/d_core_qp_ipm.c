#if defined(RUNTIME_CHECKS)
#include <stdio.h>
#include <stdlib.h>
#endif

#include "hpipm/ipm_core/d_core_qp_ipm.h"
#include "hpipm/ipm_core/d_core_qp_ipm_aux.h"


#define CORE_QP_IPM_WORKSPACE d_core_qp_ipm_workspace
#define REAL double

#define MEMSIZE_CORE_QP_IPM d_memsize_core_qp_ipm
#define CREATE_CORE_QP_IPM d_create_core_qp_ipm


#include "x_core_qp_ipm.c"
