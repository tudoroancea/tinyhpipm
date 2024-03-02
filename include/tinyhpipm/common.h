#ifndef HPIPM_COMMON_H_
#define HPIPM_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <time.h>

#define ALIGNED(VEC, BYTES) VEC __attribute__((aligned(BYTES)))

typedef size_t hpipm_size_t;
typedef int hpipm_int_t;
typedef double hpipm_float_t;

enum tinyhpipm_mode {
    SPEED_ABS,  // focus on speed, absolute IPM formulation
    SPEED,  // focus on speed, relative IPM formulation
    BALANCE,  // balanced mode, relative IPM formulation
    ROBUST,  // focus on robustness, relative IPM formulation
};

enum tinyhpipm_status {
    SUCCESS,  // found solution satisfying accuracy tolerance
    MAX_ITER,  // maximum iteration number reached
    MIN_STEP,  // minimum step length reached
    NAN_SOL,  // NaN in solution detected
    INCONS_EQ,  // unconsistent equality constraints
};

int hpipm_strcmp(char* str1, char* str2);
void hpipm_zero_memset(hpipm_size_t memsize, void* mem);

struct tinyhpipm_timer {
    struct timespec tic;
    struct timespec toc;
};

void tic(struct tinyhpipm_timer* t);

double toc(struct tinyhpipm_timer* t);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // HPIPM_COMMON_H_
