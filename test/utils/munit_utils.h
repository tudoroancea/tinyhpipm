#include "../utils/munit.h"  // TODO: make munit include non-relative
#include "tinyhpipm/blas/struct.h"
#include <stdlib.h>

#define munit_assert_isnan(x) munit_assert_int(isnan((x)), ==, 1)

#define NULL_TEST \
    { NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL }

#define NULL_SUITE \
    { NULL, NULL, NULL, 0, MUNIT_SUITE_OPTION_NONE }

#define NULL_PARAM \
    { NULL, NULL }

/*
 * @brief Helper function to generate random double between -100 and 100.
 * @return A random double between -100 and 100.
 */
double random_double(void) {
    return munit_rand_double() * 200.0 - 100.0;
}

/*
 * @brief Helper function to initialize a vector with random data.
 * @param v The vector to initialize.
 */
void init_random_vec(struct vec* v) {
    for (int i = 0; i < v->m; i++) {
        VECEL(v, i) = random_double();
    }
}

/*
 * @brief Helper function to initialize a matrix with random data.
 * @param mat The matrix to initialize.
 */
void init_random_mat(struct mat* mat) {
    for (int i = 0; i < mat->m; i++) {
        for (int j = 0; j < mat->n; ++j) {
            MATEL(mat, i, j) = random_double();
        }
    }
}

/*
 * @brief Helper function to initialize a buffer of doubles with random numbers.
 * @param buff The buffer to initialize.
 * @param n The size of the buffer.
 */
void init_random_buffer(double* buff, int n) {
    for (int i = 0; i < n; ++i) {
        buff[i] = random_double();
    }
}
