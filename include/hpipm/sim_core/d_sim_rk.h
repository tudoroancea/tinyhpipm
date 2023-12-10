#ifndef HPIPM_D_SIM_RK_H_
#define HPIPM_D_SIM_RK_H_

#include "hpipm/common.h"

#ifdef __cplusplus
extern "C" {
#endif

struct d_sim_rk_data {
    double* A_rk;  // A in butcher tableau
    double* B_rk;  // b in butcher tableau
    double* C_rk;  // c in butcher tableau
    int expl;  // erk vs irk
    int ns;  // number of stages
    hpipm_size_t memsize;
};


//
hpipm_size_t d_sim_rk_data_memsize(int ns);
//
void d_sim_rk_data_create(int ns, struct d_sim_rk_data* rk_data, void* memory);
//
void d_sim_rk_data_init_default(char* field, struct d_sim_rk_data* rk_data);
//
void d_sim_rk_data_set_all(int expl, double* A_rk, double* B_rk, double* C_rk, struct d_sim_rk_data* rk_data);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // HPIPM_D_SIM_RK_H_
