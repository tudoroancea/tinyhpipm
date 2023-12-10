#include <stdlib.h>
#include <stdio.h>

#include "blasfeo/blasfeo_target.h"
#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_s_aux.h"
#include "blasfeo/blasfeo_s_blas.h"
#include "hpipm/sim_core/s_sim_rk.h"
#include "hpipm/sim_core/s_sim_erk.h"
#include "hpipm/auxiliary/aux_mem.h"



#define BLASFEO_AXPY blasfeo_saxpy
#define BLASFEO_VEC blasfeo_svec
#define REAL float
#define SIM_ERK_ARG s_sim_erk_arg
#define SIM_RK_DATA s_sim_rk_data
#define SIM_ERK_WS s_sim_erk_ws



#define SIM_ERK_ARG_MEMSIZE s_sim_erk_arg_memsize
#define SIM_ERK_ARG_CREATE s_sim_erk_arg_create
#define SIM_ERK_ARG_SET_ALL s_sim_erk_arg_set_all

#define SIM_ERK_WS_MEMSIZE s_sim_erk_ws_memsize
#define SIM_ERK_WS_CREATE s_sim_erk_ws_create
#define SIM_ERK_WS_SET_ALL s_sim_erk_ws_set_all
#define SIM_ERK_WS_SET_NF s_sim_erk_ws_set_nf
#define SIM_ERK_WS_SET_X s_sim_erk_ws_set_x
#define SIM_ERK_WS_SET_FS s_sim_erk_ws_set_fs
#define SIM_ERK_WS_GET_X s_sim_erk_ws_get_x
#define SIM_ERK_WS_SET_P s_sim_erk_ws_set_p
#define SIM_ERK_WS_SET_ODE s_sim_erk_ws_set_ode
#define SIM_ERK_WS_SET_VDE_FOR s_sim_erk_ws_set_vde_for
#define SIM_ERK_WS_SET_ODE_ARGS s_sim_erk_ws_set_ode_args
#define SIM_ERK_WS_GET_FS s_sim_erk_ws_get_fs
#define SIM_ERK_SOLVE s_sim_erk_solve



#include "x_sim_erk.c"
