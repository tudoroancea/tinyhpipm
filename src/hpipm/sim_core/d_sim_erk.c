#include <stdio.h>
#include <stdlib.h>

#include "blasfeo/blasfeo_common.h"
#include "blasfeo/blasfeo_d_aux.h"
#include "blasfeo/blasfeo_d_blas.h"
#include "blasfeo/blasfeo_target.h"
#include "hpipm/auxiliary/aux_mem.h"
#include "hpipm/sim_core/d_sim_erk.h"
#include "hpipm/sim_core/d_sim_rk.h"


#define BLASFEO_AXPY blasfeo_daxpy
#define BLASFEO_VEC blasfeo_dvec
#define REAL double
#define SIM_ERK_ARG d_sim_erk_arg
#define SIM_RK_DATA d_sim_rk_data
#define SIM_ERK_WS d_sim_erk_ws


#define SIM_ERK_ARG_MEMSIZE d_sim_erk_arg_memsize
#define SIM_ERK_ARG_CREATE d_sim_erk_arg_create
#define SIM_ERK_ARG_SET_ALL d_sim_erk_arg_set_all

#define SIM_ERK_WS_MEMSIZE d_sim_erk_ws_memsize
#define SIM_ERK_WS_CREATE d_sim_erk_ws_create
#define SIM_ERK_WS_SET_ALL d_sim_erk_ws_set_all
#define SIM_ERK_WS_SET_NF d_sim_erk_ws_set_nf
#define SIM_ERK_WS_SET_X d_sim_erk_ws_set_x
#define SIM_ERK_WS_SET_FS d_sim_erk_ws_set_fs
#define SIM_ERK_WS_GET_X d_sim_erk_ws_get_x
#define SIM_ERK_WS_SET_P d_sim_erk_ws_set_p
#define SIM_ERK_WS_SET_ODE d_sim_erk_ws_set_ode
#define SIM_ERK_WS_SET_VDE_FOR d_sim_erk_ws_set_vde_for
#define SIM_ERK_WS_SET_ODE_ARGS d_sim_erk_ws_set_ode_args
#define SIM_ERK_WS_GET_FS d_sim_erk_ws_get_fs
#define SIM_ERK_SOLVE d_sim_erk_solve


#include "x_sim_erk.c"
