#include "../utils/munit.h"  // TODO: make munit include non-relative
#include "../utils/munit_utils.h"
#include "../utils/naive_blas.h"
#include "tinyhpipm/blas/blas1.h"
#include "tinyhpipm/blas/struct.h"
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/***************************************************************************************
 *  daxpy
 ***************************************************************************************/

static const int m = 8;
static double v1_data[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
static double v2_data[] = {8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0};
static double vnan_data[] = {NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN};
// static double vinf_data[] = {INFINITY, INFINITY, INFINITY, INFINITY, INFINITY, INFINITY, INFINITY, INFINITY};

MunitResult test_daxpy(const MunitParameter params[], void* fixture) {
    (void) params;
    (void) fixture;

    struct vec x, y, z;
    create_vec(m, &x, munit_malloc(memsize_vec(m)));
    create_vec(m, &y, munit_malloc(memsize_vec(m)));
    create_vec(m, &z, munit_malloc(memsize_vec(m)));
    pack_vec(m, v1_data, 1, &x, 0);
    pack_vec(m, v2_data, 1, &y, 0);

    daxpy(m, 0.0, &x, 0, &y, 0, &z, 0);
    for (int i = 0; i < m; i++) {
        munit_assert_double_equal(VECEL(&z, i), v2_data[i], 10);
    }

    daxpy(m, 1.0, &x, 0, &y, 0, &z, 0);
    for (int i = 0; i < m; i++) {
        munit_assert_double_equal(VECEL(&z, i), v2_data[i] + v1_data[i], 10);
    }

    daxpy(m, -1.0, &x, 0, &y, 0, &z, 0);
    for (int i = 0; i < m; i++) {
        munit_assert_double_equal(VECEL(&z, i), v2_data[i] - v1_data[i], 10);
    }
    return MUNIT_OK;
}

MunitResult test_daxpy_nan(const MunitParameter params[], void* fixture) {
    (void) params;
    (void) fixture;

    struct vec x, y, z;
    create_vec(m, &x, munit_malloc(memsize_vec(m)));
    create_vec(m, &y, munit_malloc(memsize_vec(m)));
    create_vec(m, &z, munit_malloc(memsize_vec(m)));
    pack_vec(m, v1_data, 1, &x, 0);
    pack_vec(m, vnan_data, 1, &y, 0);

    const double alpha = 0.0;
    daxpy(m, alpha, &x, 0, &y, 0, &z, 0);
    for (int i = 0; i < m; i++) {
        munit_assert_int(isnan(VECEL(&z, i)), ==, 1);
    }
    return MUNIT_OK;
}

static char* m_values[] = {"1", "10", "100", "1000", NULL};
static char* a_values[] = {"random", "min", "0.0", "1.0", "-1.0", NULL};
static MunitParameterEnum test_daxpy_random_params[] = {
        {"m", m_values},
        {"a", a_values},
        NULL_PARAM,
};

MunitResult test_daxpy_random(const MunitParameter params[], void* fixture) {
    (void) fixture;

    const char* m_str = munit_parameters_get(params, "m");
    const int m = atoi(m_str);
    if (m <= 0) {
        return MUNIT_SKIP;
    }
    const char* a_str = munit_parameters_get(params, "a");
    double a;
    if (strcmp(a_str, "random") == 0) {
        a = random_double();
    } else if (strcmp(a_str, "min") == 0) {
        a = DBL_MIN;
    } else {
        a = atof(a_str);
    }

    struct vec x, y, z;
    create_vec(m, &x, munit_malloc(memsize_vec(m)));
    create_vec(m, &y, munit_malloc(memsize_vec(m)));
    create_vec(m, &z, munit_malloc(memsize_vec(m)));

    // Initialize with random data
    init_random_vec(&x);
    init_random_vec(&y);

    // Run daxpy
    daxpy(m, a, &x, 0, &y, 0, &z, 0);

    // Compare results
    for (int i = 0; i < m; i++) {
        munit_assert_double_equal(VECEL(&z, i), a * VECEL(&x, i) + VECEL(&y, i), 10);
    }
    return MUNIT_OK;
}

// TODO: add tests with NaNs and Infs
static MunitTest daxpy_tests[] = {
        {"/normal", test_daxpy, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        {"/nan", test_daxpy_nan, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        {"/random", test_daxpy_random, NULL, NULL, MUNIT_TEST_OPTION_NONE, test_daxpy_random_params},
        NULL_TEST,
};

/***************************************************************************************
 *  dvecmul
 ***************************************************************************************/

MunitResult test_dvecmul(const MunitParameter params[], void* fixture) {
    (void) params;
    (void) fixture;
    struct vec x, y, z;
    create_vec(m, &x, munit_malloc(memsize_vec(m)));
    create_vec(m, &y, munit_malloc(memsize_vec(m)));
    create_vec(m, &z, munit_malloc(memsize_vec(m)));
    pack_vec(m, v1_data, 1, &x, 0);
    pack_vec(m, v2_data, 1, &y, 0);

    dvecmul(m, &x, 0, &y, 0, &z, 0);
    for (int i = 0; i < m; i++) {
        munit_assert_double_equal(VECEL(&z, i), v1_data[i] * v2_data[i], 10);
    }
    return MUNIT_OK;
}

MunitResult test_dvecmul_nan(const MunitParameter params[], void* fixture) {
    (void) params;
    (void) fixture;
    struct vec x, y, z;
    create_vec(m, &x, munit_malloc(memsize_vec(m)));
    create_vec(m, &y, munit_malloc(memsize_vec(m)));
    create_vec(m, &z, munit_malloc(memsize_vec(m)));
    pack_vec(m, v1_data, 1, &x, 0);
    pack_vec(m, vnan_data, 1, &y, 0);

    dvecmul(m, &x, 0, &y, 0, &z, 0);
    for (int i = 0; i < m; i++) {
        munit_assert_int(isnan(VECEL(&z, i)), ==, 1);
    }
    return MUNIT_OK;
}

static MunitTest dvecmul_tests[] = {
        {"/general", test_dvecmul, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        {"/nan", test_dvecmul_nan, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        /* Mark the end of the array with an entry where the test function is NULL */
        NULL_TEST,
};

/***************************************************************************************
 *  dvecmulacc
 ***************************************************************************************/

MunitResult test_dvecmulacc(const MunitParameter params[], void* fixture) {
    (void) params;
    (void) fixture;
    struct vec x, y, z;
    create_vec(m, &x, munit_malloc(memsize_vec(m)));
    create_vec(m, &y, munit_malloc(memsize_vec(m)));
    create_vec(m, &z, munit_malloc(memsize_vec(m)));
    pack_vec(m, v1_data, 1, &x, 0);
    pack_vec(m, v2_data, 1, &y, 0);
    for (int i = 0; i < m; i++) {
        VECEL(&z, i) = 0.0;
    }

    dvecmulacc(m, &x, 0, &y, 0, &z, 0);
    dvecmulacc(m, &x, 0, &y, 0, &x, 0);
    for (int i = 0; i < m; i++) {
        munit_assert_double_equal(VECEL(&z, i), v1_data[i] * v2_data[i], 10);
        munit_assert_double_equal(VECEL(&x, i), v1_data[i] * v2_data[i] + v1_data[i], 10);
    }
    return MUNIT_OK;
}

MunitResult test_dvecmulacc_nan(const MunitParameter params[], void* fixture) {
    (void) params;
    (void) fixture;
    struct vec x, y, z;
    create_vec(m, &x, munit_malloc(memsize_vec(m)));
    create_vec(m, &y, munit_malloc(memsize_vec(m)));
    create_vec(m, &z, munit_malloc(memsize_vec(m)));
    pack_vec(m, v1_data, 1, &x, 0);
    pack_vec(m, vnan_data, 1, &y, 0);
    for (int i = 0; i < m; i++) {
        VECEL(&z, i) = 0.0;
    }

    dvecmulacc(m, &x, 0, &y, 0, &z, 0);
    for (int i = 0; i < m; i++) {
        munit_assert_int(isnan(VECEL(&z, i)), ==, 1);
    }
    return MUNIT_OK;
}

static MunitTest dvecmulacc_tests[] = {
        {"/general", test_dvecmulacc, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        {"/nan", test_dvecmulacc_nan, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        /* Mark the end of the array with an entry where the test function is NULL */
        NULL_TEST,
};

/************************************************************************************
 *  dvecmuldot
 ***************************************************************************************/

MunitResult test_dvecmuldot(const MunitParameter params[], void* fixture) {
    (void) params;
    (void) fixture;
    struct vec x, y, z;
    create_vec(m, &x, munit_malloc(memsize_vec(m)));
    create_vec(m, &y, munit_malloc(memsize_vec(m)));
    create_vec(m, &z, munit_malloc(memsize_vec(m)));
    pack_vec(m, v1_data, 1, &x, 0);
    pack_vec(m, v2_data, 1, &y, 0);

    const double ddot = dvecmuldot(m, &x, 0, &y, 0, &z, 0);
    munit_assert_double_equal(ddot, v1_data[0] * v2_data[0] + v1_data[1] * v2_data[1] + v1_data[2] * v2_data[2] + v1_data[3] * v2_data[3] + v1_data[4] * v2_data[4] + v1_data[5] * v2_data[5] + v1_data[6] * v2_data[6] + v1_data[7] * v2_data[7], 10);
    for (int i = 0; i < m; i++) {
        munit_assert_double_equal(VECEL(&z, i), v1_data[i] * v2_data[i], 10);
    }
    return MUNIT_OK;
}

static MunitTest dvecmuldot_tests[] = {
        {"/general", test_dvecmuldot, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        /* Mark the end of the array with an entry where the test function is NULL */
        NULL_TEST,
};


/***************************************************************************************
 *  main
 ***************************************************************************************/
static MunitSuite all_suites[] = {
        {"/daxpy", daxpy_tests, NULL, 1, MUNIT_SUITE_OPTION_NONE},
        {"/dvecmul", dvecmul_tests, NULL, 1, MUNIT_SUITE_OPTION_NONE},
        {"/dvecmulacc", dvecmulacc_tests, NULL, 1, MUNIT_SUITE_OPTION_NONE},
        {"/dvecmuldot", dvecmuldot_tests, NULL, 1, MUNIT_SUITE_OPTION_NONE},
        NULL_SUITE,
};

static const MunitSuite file_suite = {"/blas/blas1", NULL, all_suites, 1, MUNIT_SUITE_OPTION_NONE};

int main(int argc, char* argv[]) {
    return munit_suite_main(&file_suite, NULL, argc, argv);
}
