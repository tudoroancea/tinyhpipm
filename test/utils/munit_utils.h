#include "../utils/munit.h"  // TODO: make munit include non-relative
#include "tinyhpipm/blas/struct.h"
#include <stdlib.h>


#define NULL_TEST \
    { NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL }

#define NULL_SUITE \
    { NULL, NULL, NULL, 0, MUNIT_SUITE_OPTION_NONE }

#define NULL_PARAM \
    { NULL, NULL }

 // Helper function to generate random double between -100 and 100
double random_double(void) {
    return munit_rand_double() * 200.0 - 100.0;
}

// Helper function to initialize a vector with random data
void init_random_vec(struct vec* v) {
    for (int i = 0; i < v->m; i++) {
        VECEL(v, i) = random_double();
    }
}
