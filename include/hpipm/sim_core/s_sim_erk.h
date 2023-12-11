#ifndef HPIPM_S_SIM_ERK_H_
#define HPIPM_S_SIM_ERK_H_

#ifdef __cplusplus
extern "C" {
#endif

struct s_sim_erk_arg {
    struct s_sim_rk_data* rk_data;  // integrator data
    float h;  // step size
    int steps;  // number of steps
    //	int for_sens; // compute adjoint sensitivities
    //	int adj_sens; // compute adjoint sensitivities
    hpipm_size_t memsize;
};


struct s_sim_erk_ws {
    void (*ode)(int t, float* x, float* p, void* ode_args, float* xdot);  // function pointer to ode
    void (*vde_for)(int t, float* x, float* p, void* ode_args, float* xdot);  // function pointer to forward vde
    void (*vde_adj)(int t, float* adj_in, void* ode_args, float* adj_out);  // function pointer to adjoint vde
    void* ode_args;  // pointer to ode args
    struct s_sim_erk_arg* erk_arg;  // erk arg
    float* K;  // internal variables
    float* x_for;  // states and forward sensitivities
    float* x_traj;  // states at all steps
    float* l;  // adjoint sensitivities
    float* p;  // parameter
    float* x_tmp;  // temporary states and forward sensitivities
    float* adj_in;
    float* adj_tmp;
    int nx;  // number of states
    int np;  // number of parameters
    int nf;  // number of forward sensitivities
    int na;  // number of adjoint sensitivities
    int nf_max;  // max number of forward sensitivities
    int na_max;  // max number of adjoint sensitivities
    hpipm_size_t memsize;
};


//
hpipm_size_t s_sim_erk_arg_memsize();
//
void s_sim_erk_arg_create(struct s_sim_erk_arg* erk_arg, void* mem);
//
void s_sim_erk_arg_set_all(struct s_sim_rk_data* rk_data, float h, int steps, struct s_sim_erk_arg* erk_arg);

//
hpipm_size_t s_sim_erk_ws_memsize(struct s_sim_erk_arg* erk_arg, int nx, int np, int nf_max, int na_max);
//
void s_sim_erk_ws_create(struct s_sim_erk_arg* erk_arg, int nx, int np, int nf_max, int na_max, struct s_sim_erk_ws* work, void* memory);
//
void s_sim_erk_ws_set_all(int nf, int na, float* x, float* fs, float* bs, float* p, void (*ode)(int t, float* x, float* p, void* ode_args, float* xdot), void (*vde_for)(int t, float* x, float* p, void* ode_args, float* xdot), void (*vde_adj)(int t, float* adj_in, void* ode_args, float* adj_out), void* ode_args, struct s_sim_erk_ws* work);
// number of directions for forward sensitivities
void s_sim_erk_ws_set_nf(int* nf, struct s_sim_erk_ws* work);
// parameters (e.g. inputs)
void s_sim_erk_ws_set_p(float* p, struct s_sim_erk_ws* work);
// state
void s_sim_erk_ws_set_x(float* x, struct s_sim_erk_ws* work);
// forward sensitivities
void s_sim_erk_ws_set_fs(float* fs, struct s_sim_erk_ws* work);
// ode funtion
void s_sim_erk_ws_set_ode(void (*ode)(int t, float* x, float* p, void* ode_args, float* xdot), struct s_sim_erk_ws* work);
// forward vde function
void s_sim_erk_ws_set_vde_for(void (*ode)(int t, float* x, float* p, void* ode_args, float* xdot), struct s_sim_erk_ws* work);
// ode_args, passed straight to the ode/vde_for/vde_adj functions
void s_sim_erk_ws_set_ode_args(void* ode_args, struct s_sim_erk_ws* work);
// state
void s_sim_erk_ws_get_x(struct s_sim_erk_ws* work, float* x);
// forward sensitivities
void s_sim_erk_ws_get_fs(struct s_sim_erk_ws* work, float* fs);
//
void s_sim_erk_solve(struct s_sim_erk_arg* arg, struct s_sim_erk_ws* work);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // HPIPM_D_SIM_ERK_H_
