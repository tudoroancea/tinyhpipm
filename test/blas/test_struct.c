#include "../utils/munit.h"  // TODO: make it non-relative
#include "tinyhpipm/blas/struct.h"
#include <stdio.h>
#include <stdlib.h>


#define NULL_TEST \
    { NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL }

#define NULL_SUITE \
    { NULL, NULL, NULL, 0, MUNIT_SUITE_OPTION_NONE }

/***************************************************************************************
 *  memsize
 ***************************************************************************************/

MunitResult test_memsize_mat(const MunitParameter params[], void* fixture) {
    munit_assert_int(memsize_mat(10, 4), ==, (12 * 4 + 4 * 4) * sizeof(double));
    munit_assert_int(memsize_mat(10, 5), ==, (12 * 8 + 4 * 4) * sizeof(double));
    return MUNIT_OK;
}

MunitResult test_memsize_vec(const MunitParameter params[], void* fixture) {
    munit_assert_int(memsize_vec(10), ==, 12 * sizeof(double));
    int n = 1;
    for (int i = 0; i < 10; i++) {
        n *= D_PS;
        munit_assert_int(memsize_vec(n), ==, n * sizeof(double));
    }
    return MUNIT_OK;
}

static MunitTest memsize_tests[] = {
        {"/memsize_vec", test_memsize_vec, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        {"/memsize_mat", test_memsize_mat, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        /* Mark the end of the array with an entry where the test function is NULL */
        NULL_TEST,
};

/***************************************************************************************
 *  create
 ***************************************************************************************/

MunitResult test_create_mat_1(const MunitParameter params[], void* fixture) {
    struct mat A;
    create_mat(10, 4, &A, munit_malloc(memsize_mat(10, 4)));
    munit_assert_int(A.m, ==, 10);
    munit_assert_int(A.n, ==, 4);
    munit_assert_int(A.pm, ==, 12);
    munit_assert_int(A.cn, ==, 4);
    free_mat(&A);
    return MUNIT_OK;
}

MunitResult test_create_mat_2(const MunitParameter params[], void* fixture) {
    struct mat A;
    create_mat(5, 11, &A, munit_malloc(memsize_mat(5, 11)));
    munit_assert_int(A.m, ==, 5);
    munit_assert_int(A.n, ==, 11);
    munit_assert_int(A.pm, ==, 8);
    munit_assert_int(A.cn, ==, 12);
    free_mat(&A);
    return MUNIT_OK;
}

MunitResult test_create_vec(const MunitParameter params[], void* fixture) {
    struct vec v;
    create_vec(10, &v, munit_malloc(memsize_vec(10)));
    munit_assert_int(v.m, ==, 10);
    munit_assert_int(v.pm, ==, 12);
    free_vec(&v);
    return MUNIT_OK;
}

static MunitTest create_tests[] = {
        {"/create_mat_1", test_create_mat_1, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        {"/create_mat_2", test_create_mat_2, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        {"/create_vec", test_create_vec, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        /* Mark the end of the array with an entry where the test function is NULL */
        NULL_TEST,
};

/***************************************************************************************
 *  indexing
 ***************************************************************************************/

MunitResult test_index_mat(const MunitParameter params[], void* fixture) {
    struct mat A;
    int m = 5, n = 4;
    double* mem = (double*) munit_malloc(memsize_mat(m, n));
    for (int i = 0; i < 8 * 4; i++) {
        mem[i] = 1.0;
    }
    create_mat(m, n, &A, (void*) mem);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; ++j) {
            munit_assert_double(MATEL(&A, i, j), ==, 1.0);
        }
    }
    return MUNIT_OK;
}

MunitResult test_index_vec(const MunitParameter params[], void* fixture) {
    struct vec v;
    int m = 5;
    double* mem = (double*) munit_malloc(memsize_vec(m));
    for (int i = 0; i < m; i++) {
        mem[i] = 1.0;
    }
    create_vec(m, &v, (void*) mem);
    for (int i = 0; i < m; i++) {
        munit_assert_double(VECEL(&v, i), ==, 1.0);
    }
    return MUNIT_OK;
}

static MunitTest indexing_tests[] = {
        {"/index_mat", test_index_mat, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        {"/index_vec", test_index_vec, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        /* Mark the end of the array with an entry where the test function is NULL */
        NULL_TEST,
};

/***************************************************************************************
 *  pack
 ***************************************************************************************/
// pack whole mat
MunitResult test_pack_mat_all_rows_all_cols(const MunitParameter params[], void* fixture) {
    struct mat sA;
    int m = 15, n = 12;
    create_mat(m, n, &sA, munit_malloc(memsize_mat(m, n)));
    double* A = (double*) munit_malloc(m * n);
    for (int k = 0; k < m * n; k++) {
        A[k] = (double) k;
    }
    pack_mat(m, n, A, m, &sA, 0, 0);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; ++j) {
            munit_assert_double(MATEL(&sA, i, j), ==, (double) (j * m + i));
        }
    }
    return MUNIT_OK;
}
// pack only some cols
MunitResult test_pack_mat_all_rows_some_cols(const MunitParameter params[], void* fixture) {

    struct mat sA;
    int m = 15, n = 12;
    create_mat(m, n, &sA, munit_malloc(memsize_mat(m, n)));
    double* A = (double*) munit_malloc(m * n);
    for (int k = 0; k < m * n; k++) {
        A[k] = (double) k;
    }
    pack_mat(m, n / 2, A, 2 * m, &sA, 0, 0);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; ++j) {
            if (j < n / 2) {
                munit_assert_double(MATEL(&sA, i, j), ==, (double) (2 * j * m + i));

            } else {
                munit_assert_double(MATEL(&sA, i, j), ==, 0.0);
            }
        }
    }
    return MUNIT_OK;
}
// pack only some rows
MunitResult test_pack_mat_some_rows_all_cols(const MunitParameter params[], void* fixture) {
    struct mat sA;
    int m = 15, n = 12;
    create_mat(m, n, &sA, munit_malloc(memsize_mat(m, n)));
    double* A = (double*) munit_malloc(m * n);
    for (int k = 0; k < m * n; k++) {
        A[k] = (double) k;
    }
    pack_mat(m / 2, n, A, m, &sA, 0, 0);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; ++j) {
            if (i < m / 2) {
                munit_assert_double(MATEL(&sA, i, j), ==, (double) (j * m + i));

            } else {
                munit_assert_double(MATEL(&sA, i, j), ==, 0.0);
            }
        }
    }
    return MUNIT_OK;
}
// pack only some rows and cols
MunitResult test_pack_mat_some_rows_some_cols(const MunitParameter params[], void* fixture) {
    struct mat sA;
    int m = 15, n = 12;
    create_mat(m, n, &sA, munit_malloc(memsize_mat(m, n)));
    double* A = (double*) munit_malloc(m * n);
    for (int k = 0; k < m * n; k++) {
        A[k] = (double) k;
    }
    pack_mat(m / 2, n / 2, A, 2 * m, &sA, 0, 0);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; ++j) {
            if (j < n / 2 && i < m / 2) {
                munit_assert_double(MATEL(&sA, i, j), ==, (double) (2 * j * m + i));

            } else {
                munit_assert_double(MATEL(&sA, i, j), ==, 0.0);
            }
        }
    }
    return MUNIT_OK;
}
// pack whole mat
MunitResult test_pack_tran_mat_all_rows_all_cols(const MunitParameter params[], void* fixture) {
    struct mat sA;
    int m = 15, n = 12;
    create_mat(m, n, &sA, munit_malloc(memsize_mat(m, n)));
    double* A = (double*) munit_malloc(m * n);  // A is to be understood of dims n x m
    for (int k = 0; k < m * n; k++) {
        A[k] = (double) k;
    }
    pack_tran_mat(n, m, A, n, &sA, 0, 0);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; ++j) {
            munit_assert_double(MATEL(&sA, i, j), ==, (double) (i * n + j));
        }
    }
    return MUNIT_OK;
}
MunitResult test_pack_l_mat_all_rows_all_cols(const MunitParameter params[], void* fixture) {
    struct mat sA;
    int m = 15, n = 12;
    create_mat(m, n, &sA, munit_malloc(memsize_mat(m, n)));
    double* A = (double*) munit_malloc(m * n);  // A is to be understood of dims m x n
    for (int k = 0; k < m * n; k++) {
        A[k] = (double) k;
    }
    pack_l_mat(m, n, A, m, &sA, 0, 0);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; ++j) {
            if (i >= j) {
                munit_assert_double(MATEL(&sA, i, j), ==, (double) (j * m + i));
            } else {
                munit_assert_double(MATEL(&sA, i, j), ==, 0.0);
            }
        }
    }
    return MUNIT_OK;
}
MunitResult test_pack_u_mat_all_rows_all_cols(const MunitParameter params[], void* fixture) {
    struct mat sA;
    int m = 15, n = 12;
    create_mat(m, n, &sA, munit_malloc(memsize_mat(m, n)));
    double* A = (double*) munit_malloc(m * n);  // A is to be understood of dims m x n
    for (int k = 0; k < m * n; k++) {
        A[k] = (double) k;
    }
    pack_u_mat(m, n, A, m, &sA, 0, 0);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; ++j) {
            if (i <= j) {
                munit_assert_double(MATEL(&sA, i, j), ==, (double) (j * m + i));
            } else {
                munit_assert_double(MATEL(&sA, i, j), ==, 0.0);
            }
        }
    }
    return MUNIT_OK;
}


MunitResult test_pack_vec_all_elts(const MunitParameter params[], void* fixture) {
    struct vec sv;
    int m = 15;
    create_vec(m, &sv, munit_malloc(memsize_vec(m)));
    double* v = (double*) munit_malloc(m);
    for (int k = 0; k < m; k++) {
        v[k] = (double) k;
    }
    pack_vec(m, v, 1, &sv, 0);
    for (int i = 0; i < m; i++) {
        munit_assert_double(VECEL(&sv, i), ==, (double) i);
    }
    return MUNIT_OK;
}
MunitResult test_pack_vec_some_elts(const MunitParameter params[], void* fixture) {
    struct vec sv;
    int m = 15;
    create_vec(m, &sv, munit_malloc(memsize_vec(m)));
    double* v = (double*) munit_malloc(m);
    for (int k = 0; k < m; k++) {
        v[k] = (double) k;
    }
    pack_vec(m / 2, v, 2, &sv, 0);
    for (int i = 0; i < m; i++) {
        if (i < m / 2) {
            munit_assert_double(VECEL(&sv, i), ==, (double) (2 * i));
        } else {
            munit_assert_double(VECEL(&sv, i), ==, 0.0);
        }
    }
    return MUNIT_OK;
}

static MunitTest packing_tests[] = {
        {"/pack_mat_all_rows_all_cols", test_pack_mat_all_rows_all_cols, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        {"/pack_mat_all_rows_some_cols", test_pack_mat_all_rows_some_cols, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        {"/pack_mat_some_rows_all_cols", test_pack_mat_some_rows_all_cols, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        {"/pack_mat_some_rows_some_cols", test_pack_mat_some_rows_some_cols, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        {"/pack_tran_mat_all_rows_all_cols", test_pack_tran_mat_all_rows_all_cols, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        {"/pack_u_mat_all_rows_all_cols", test_pack_u_mat_all_rows_all_cols, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        {"/pack_l_mat_all_rows_all_cols", test_pack_l_mat_all_rows_all_cols, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        {"/pack_vec_all_elts", test_pack_vec_all_elts, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        {"/pack_vec_some_elts", test_pack_vec_some_elts, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        /* Mark the end of the array with an entry where the test function is NULL */
        NULL_TEST,
};

/***************************************************************************************
 *  unpack
 ***************************************************************************************/

MunitResult test_unpack_mat(const MunitParameter params[], void* fiture) {
    struct mat sA;
    int m = 15, n = 12;
    create_mat(m, n, &sA, munit_malloc(memsize_mat(m, n)));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            MATEL(&sA, i, j) = (double) (j * m + i);
        }
    }
    double* A = (double*) munit_malloc(m * n);  // A is to be understood of dims m x n
    unpack_mat(m, n, &sA, 0, 0, A, m);
    for (int k = 0; k < m * n; k++) {
        munit_assert_double(A[k], ==, (double) k);
    }
    return MUNIT_OK;
}

MunitResult test_unpack_tran_mat(const MunitParameter params[], void* fiture) {
    struct mat sA;
    int m = 12, n = 12;
    create_mat(m, n, &sA, munit_malloc(memsize_mat(m, n)));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; ++j) {
            MATEL(&sA, i, j) = (double) (j * m + i);
        }
    }
    double* A = (double*) munit_malloc(m * n);  // A is to be understood of dims m x n
    unpack_tran_mat(m, n, &sA, 0, 0, A, m);
    struct vec sv;
    create_vec(m * n, &sv, munit_malloc(memsize_vec(m * n)));
    pack_vec(m * n, A, 1, &sv, 0);

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            munit_assert_double(A[i * n + j], ==, MATEL(&sA, i, j));
        }
    }
    return MUNIT_OK;
}

MunitResult test_unpack_vec(const MunitParameter params[], void* fiture) {
    struct vec sv;
    int m = 15;
    create_vec(m, &sv, munit_malloc(memsize_vec(m)));
    for (int i = 0; i < m; i++) {
        VECEL(&sv, i) = (double) i;
    }
    double* v = (double*) munit_malloc(m);
    unpack_vec(m, &sv, 0, v, 1);
    for (int k = 0; k < m; k++) {
        munit_assert_double(v[k], ==, (double) k);
    }
    return MUNIT_OK;
}

static MunitTest unpacking_tests[] = {
        {"/unpack_mat", test_unpack_mat, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        {"/unpack_tran_mat", test_unpack_tran_mat, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        {"/unpack_vec", test_unpack_vec, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        /* Mark the end of the array with an entry where the test function is NULL */
        NULL_TEST,
};


/***************************************************************************************
 *  main
 ***************************************************************************************/
static MunitSuite all_suites[] = {
        {"/memsize", memsize_tests, NULL, 1, MUNIT_SUITE_OPTION_NONE},
        {"/create", create_tests, NULL, 1, MUNIT_SUITE_OPTION_NONE},
        {"/indexing", indexing_tests, NULL, 1, MUNIT_SUITE_OPTION_NONE},
        {"/packing", packing_tests, NULL, 1, MUNIT_SUITE_OPTION_NONE},
        {"/unpacking", unpacking_tests, NULL, 1, MUNIT_SUITE_OPTION_NONE},
        NULL_SUITE,
};

static const MunitSuite giga_suite = {"/blas/struct", NULL, all_suites, 1, MUNIT_SUITE_OPTION_NONE};

int main(int argc, char* argv[]) {
    return munit_suite_main(&giga_suite, NULL, argc, argv);
}
