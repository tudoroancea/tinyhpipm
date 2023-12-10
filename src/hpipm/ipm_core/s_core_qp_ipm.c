#if defined(RUNTIME_CHECKS)
#include <stdio.h>
#include <stdlib.h>
#endif

#include "hpipm/ipm_core/s_core_qp_ipm.h"
#include "hpipm/ipm_core/s_core_qp_ipm_aux.h"


#define CORE_QP_IPM_WORKSPACE s_core_qp_ipm_workspace
#define REAL float

#define MEMSIZE_CORE_QP_IPM s_memsize_core_qp_ipm
#define CREATE_CORE_QP_IPM s_create_core_qp_ipm


#include "x_core_qp_ipm.c"
