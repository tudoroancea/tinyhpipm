#ifndef HPIPM_COMMON_H_
#define HPIPM_COMMON_H_

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef size_t hpipm_size_t;

enum hpipm_mode {
    SPEED_ABS,  // focus on speed, absolute IPM formulation
    SPEED,  // focus on speed, relative IPM formulation
    BALANCE,  // balanced mode, relative IPM formulation
    ROBUST,  // focus on robustness, relative IPM formulation
};

enum hpipm_status {
    SUCCESS,  // found solution satisfying accuracy tolerance
    MAX_ITER,  // maximum iteration number reached
    MIN_STEP,  // minimum step length reached
    NAN_SOL,  // NaN in solution detected
    INCONS_EQ,  // unconsistent equality constraints
};

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // HPIPM_COMMON_H_
