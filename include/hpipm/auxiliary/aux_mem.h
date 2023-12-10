#ifndef HPIPM_AUX_MEM_H_
#define HPIPM_AUX_MEM_H_

#include "hpipm/common.h"

#ifdef __cplusplus
extern "C" {
#endif

void hpipm_zero_memset(hpipm_size_t memsize, void* mem);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // HPIPM_AUX_MEM_H_
