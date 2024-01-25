// tests for struct.h
#include "../munit.h"  // TODO: make it non-relative
#include "tinyhpipm/blas/struct.h"
#include <stdio.h>
#include <stdlib.h>

// one suite per function, several tests per suite for different mat/vec sizes.
// after the first tests for matrix creation, we create setup and teardown functions
// to create and destroy mat and vec structs.
// suite naming: /blas/struct/<function name>
// test naming: /blas/struct/<function name>/<mat or vec size>

MunitResult test_memsize_mat(const MunitParameter params[], void* fixture) {
    munit_assert_int(memsize_mat(10, 4), ==, (12 * 4 + 4 * 4) * sizeof(double));
    munit_assert_int(memsize_mat(10, 5), ==, (12 * 8 + 4 * 4) * sizeof(double));
}

MunitResult test_memsize_vec(const MunitParameter params[], void* fixture) {
    munit_assert_int(memsize_vec(10), ==, 12 * sizeof(double));
    int n = 1;
    for (int i = 0; i < 10; i++) {
        n *= D_PS;
        munit_assert_int(memsize_vec(n), ==, n * sizeof(double));
    }
}

#define NULL_TEST \
    { NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL }
static MunitTest memsize_tests[] = {
        {"/memsize_vec", test_memsize_vec, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        {"/memsize_mat", test_memsize_mat, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        /* Mark the end of the array with an entry where the test function is NULL */
        NULL_TEST,
};

static const MunitSuite suite = {
        "/blas/struct/memsize",
        memsize_tests,
        NULL,
        1,
        MUNIT_SUITE_OPTION_NONE};

int main(int argc, char* argv[]) {
    return munit_suite_main(&suite, NULL, argc, argv);
}

// int main() {
//     struct vec u, v;
//     create_vec(10, &u, malloc(memsize_vec(10)));
//     create_vec(4, &v, malloc(memsize_vec(4)));
//     munit_assert_int(u.m, ==, 10);
//     munit_assert_int(u.pm, ==, 12);
//     munit_assert_int(v.m, ==, 4);
//     munit_assert_int(v.pm, ==, 4);
//     return 0;
// }
