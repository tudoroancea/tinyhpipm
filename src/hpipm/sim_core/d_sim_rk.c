#include <stdio.h>
#include <stdlib.h>

#include "hpipm/auxiliary/aux_string.h"
#include "hpipm/sim_core/d_sim_rk.h"

#define REAL double
#define SIM_RK_DATA d_sim_rk_data


#define SIM_RK_DATA_MEMSIZE d_sim_rk_data_memsize
#define SIM_RK_DATA_CREATE d_sim_rk_data_create
#define SIM_RK_DATA_INIT_DEFAULT d_sim_rk_data_init_default
#define SIM_RK_DATA_SET_ALL d_sim_rk_data_set_all


#include "x_sim_rk.c"
