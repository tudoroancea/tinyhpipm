#ifndef HPIPM_TIMING_H_
#define HPIPM_TIMING_H_

#include "blasfeo/blasfeo_timing.h"


#ifdef __cplusplus
extern "C" {
#endif


typedef blasfeo_timer hpipm_timer;


// create time point (start timer)
void hpipm_tic(hpipm_timer* t);
// return time elapsed since specified tic (stop timer)
double hpipm_toc(hpipm_timer* t);


#ifdef __cplusplus
}  // #extern "C"
#endif


#endif  // HPIPM_TIMING_H_
