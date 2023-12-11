#ifndef HPIPM_S_SIM_RK_H_
#define HPIPM_S_SIM_RK_H_

#include "hpipm/common.h"

#ifdef __cplusplus
extern "C" {
#endif

struct s_sim_rk_data {
    float* A_rk;  // A in butcher tableau
    float* B_rk;  // b in butcher tableau
    float* C_rk;  // c in butcher tableau
    int expl;  // erk vs irk
    int ns;  // number of stages
    hpipm_size_t memsize;
};


//
hpipm_size_t s_sim_rk_data_memsize(int ns);
//
void s_sim_rk_data_create(int ns, struct s_sim_rk_data* rk_data, void* memory);
//
void s_sim_rk_data_init_default(char* field, struct s_sim_rk_data* rk_data);
//
void s_sim_rk_data_set_all(int expl, float* A_rk, float* B_rk, float* C_rk, struct s_sim_rk_data* rk_data);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // HPIPM_S_SIM_RK_H_
