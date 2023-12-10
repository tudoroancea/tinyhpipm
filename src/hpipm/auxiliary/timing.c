#include "hpipm/auxiliary/timing.h"

void hpipm_tic(hpipm_timer* t) {
    blasfeo_tic(t);
    return;
}


double hpipm_toc(hpipm_timer* t) {
    return blasfeo_toc(t);
}
