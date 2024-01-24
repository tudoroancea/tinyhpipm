#include "hpipm/common.h"
#include <stdint.h>

int hpipm_strcmp(char* str1, char* str2) {
    int i = 0;
    while (str1[i] == str2[i] & str1[i] != '\0')
        i++;
    if ((str1[i] > str2[i]) | (str1[i] < str2[i]))
        return 0;
    else
        return 1;
}

void hpipm_zero_memset(hpipm_size_t memsize, void* mem) {
    hpipm_size_t ii;
    hpipm_size_t memsize_m8 = memsize / 8;  // sizeof(double) is 8
    hpipm_size_t memsize_r8 = memsize % 8;
    double* double_ptr = mem;
    ii = 0;
    if (memsize_m8 > 7) {
        for (; ii < memsize_m8 - 7; ii += 8) {
            double_ptr[ii + 0] = 0.0;
            double_ptr[ii + 1] = 0.0;
            double_ptr[ii + 2] = 0.0;
            double_ptr[ii + 3] = 0.0;
            double_ptr[ii + 4] = 0.0;
            double_ptr[ii + 5] = 0.0;
            double_ptr[ii + 6] = 0.0;
            double_ptr[ii + 7] = 0.0;
        }
    }
    for (; ii < memsize_m8; ii++) {
        double_ptr[ii] = 0.0;
    }
    char* char_ptr = (char*) (&double_ptr[ii]);
    for (ii = 0; ii < memsize_r8; ii++) {
        char_ptr[ii] = 0;
    }
}

// #if (defined __APPLE__)
// void tic(struct timer* t) {
//     /* read current clock cycles */
//     t->tic = mach_absolute_time();
// }

// double toc(struct timer* t) {
//     uint64_t duration; /* elapsed time in clock cycles*/

//     t->toc = mach_absolute_time();
//     duration = t->toc - t->tic;

//     /*conversion from clock cycles to nanoseconds*/
//     mach_timebase_info(&(t->tinfo));
//     duration *= t->tinfo.numer;
//     duration /= t->tinfo.denom;

//     return (double) duration / 1e9;
// }
// #else
// #if __STDC_VERSION__ >= 199901L  // C99 Mode

// /* read current time */
// void tic(struct timer* t) {
//     gettimeofday(&t->tic, 0);
// }

// /* return time passed since last call to tic on this timer */
// double toc(struct timer* t) {
//     struct timeval temp;

//     gettimeofday(&t->toc, 0);

//     if ((t->toc.tv_usec - t->tic.tv_usec) < 0) {
//         temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec - 1;
//         temp.tv_usec = 1000000 + t->toc.tv_usec - t->tic.tv_usec;
//     } else {
//         temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec;
//         temp.tv_usec = t->toc.tv_usec - t->tic.tv_usec;
//     }

//     return (double) temp.tv_sec + (double) temp.tv_usec / 1e6;
// }

// #else  // ANSI C Mode

/* read current time */
void tic(struct timer* t) {
    clock_gettime(CLOCK_MONOTONIC, &t->tic);
}


/* return time passed since last call to tic on this timer */
double toc(struct timer* t) {
    struct timespec temp;

    clock_gettime(CLOCK_MONOTONIC, &t->toc);

    if ((t->toc.tv_nsec - t->tic.tv_nsec) < 0) {
        temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec - 1;
        temp.tv_nsec = 1000000000 + t->toc.tv_nsec - t->tic.tv_nsec;
    } else {
        temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec;
        temp.tv_nsec = t->toc.tv_nsec - t->tic.tv_nsec;
    }

    return (double) temp.tv_sec + (double) temp.tv_nsec / 1e9;
}

// #endif  // __STDC_VERSION__ >= 199901L
// #endif
