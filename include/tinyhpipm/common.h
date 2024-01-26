#ifndef HPIPM_COMMON_H_
#define HPIPM_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>

#define ALIGNED(VEC, BYTES) VEC __attribute__((aligned(BYTES)))

typedef size_t hpipm_size_t;


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

// #if (defined __APPLE__)

// #include <mach/mach_time.h>

// /** A structure for keeping internal timer data. */
// struct timer {
//     uint64_t tic;
//     uint64_t toc;
//     mach_timebase_info_data_t tinfo;
// };

// #else
/* Use POSIX clock_gettime() for timing on non-Windows machines. */
#include <time.h>

// #if __STDC_VERSION__ >= 199901L  // C99 Mode

// #include <sys/stat.h>
// #include <sys/time.h>

// struct timer {
//     struct timeval tic;
//     struct timeval toc;
// };

// #else  // ANSI C Mode

struct tinyhpipm_timer {
    struct timespec tic;
    struct timespec toc;
};

// #endif  // __STDC_VERSION__ >= 199901L
// #endif // __APPLE__

/** A function for measurement of the current time. */
void tic(struct tinyhpipm_timer* t);

/** A function which returns the elapsed time. */
double toc(struct tinyhpipm_timer* t);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // HPIPM_COMMON_H_
