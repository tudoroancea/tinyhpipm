#include <stdio.h>
#include <stdlib.h>

#include "hpipm/auxiliary/aux_string.h"
#include "hpipm/sim_core/s_sim_rk.h"


#define REAL float
#define SIM_RK_DATA s_sim_rk_data


#define SIM_RK_DATA_MEMSIZE s_sim_rk_data_memsize
#define SIM_RK_DATA_CREATE s_sim_rk_data_create
#define SIM_RK_DATA_INIT_DEFAULT s_sim_rk_data_init_default
#define SIM_RK_DATA_SET_ALL s_sim_rk_data_set_all


#include "x_sim_rk.c"
