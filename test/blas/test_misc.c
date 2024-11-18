#include "../utils/munit.h"  // TODO: make munit include non-relative
#include "../utils/munit_utils.h"  // TODO: make munit include non-relative
#include "tinyhpipm/blas/blas1.h"
#include "tinyhpipm/blas/misc.h"
#include "tinyhpipm/blas/print.h"
#include "tinyhpipm/blas/struct.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/***************************************************************************************
 *  get and set
 ***************************************************************************************/

MunitResult test_dvecex_sp(const MunitParameter params[], void* fixture) {
    (void) params;
    (void) fixture;
    return MUNIT_OK;
}
MunitResult test_drowin(const MunitParameter params[], void* fixture) {
    (void) params;
    (void) fixture;
    return MUNIT_OK;
}
MunitResult test_drowex(const MunitParameter params[], void* fixture) {
    (void) params;
    (void) fixture;
    return MUNIT_OK;
}
MunitResult test_dcolin(const MunitParameter params[], void* fixture) {
    (void) params;
    (void) fixture;
    return MUNIT_OK;
}
MunitResult test_dcolex(const MunitParameter params[], void* fixture) {
    (void) params;
    (void) fixture;
    return MUNIT_OK;
}
MunitResult test_ddiaex(const MunitParameter params[], void* fixture) {
    (void) params;
    (void) fixture;
    const int m = 5, n = 5;
    double* Ap = munit_newa(double, n* m);
    init_random_buffer(Ap, n * m);
    struct mat A;
    create_mat(m, n, &A, munit_malloc(memsize_mat(m, n)));
    pack_mat(m, n, Ap, m, &A, 0, 0);
    struct vec v, v2;
    create_vec(m, &v, munit_malloc(memsize_vec(m)));
    create_vec(m, &v2, munit_malloc(memsize_vec(m + 2)));
    daxpy(m + 1, -1.0, &v2, 0, &v2, 0, &v2, 0);  // convoluted way to set v = 0.0

    ddiaex(m, 1.0, &A, 0, 0, &v, 0);
    for (int i = 0; i < m; ++i) {
        munit_assert_double_equal(VECEL(&v, i), MATEL(&A, i, i), 10);
    }

    ddiaex(m, NAN, &A, 0, 0, &v, 0);
    for (int i = 0; i < m; ++i) {
        munit_assert_isnan(VECEL(&v, i));
    }

    ddiaex(m - 1, 2.0, &A, 0, 0, &v2, 1);
    munit_assert_double_equal(VECEL(&v2, 0), 0.0, 10);
    munit_assert_double_equal(VECEL(&v2, m), 0.0, 10);
    for (int i = 0; i < m - 1; ++i) {
        munit_assert_double_equal(VECEL(&v2, i + 1), 2.0 * MATEL(&A, i, i), 10);
    }

    return MUNIT_OK;
}
MunitResult test_dvecse(const MunitParameter params[], void* fixture) {
    (void) params;
    (void) fixture;
    return MUNIT_OK;
}
MunitResult test_dgese(const MunitParameter params[], void* fixture) {
    (void) params;
    (void) fixture;
    const int m = 5, n = 7;
    struct mat A;
    create_mat(m, n, &A, munit_malloc(memsize_mat(m, n)));

    double alpha = 1.0;
    dgese(m, n, alpha, &A, 0, 0);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            munit_assert_double(MATEL(&A, i, j), ==, alpha);
        }
    }

    alpha = 2.0;
    dgese(m - 1, n - 1, alpha, &A, 1, 1);
    for (int i = 1; i < m; ++i) {
        for (int j = 1; j < n; ++j) {
            munit_assert_double(MATEL(&A, i, j), ==, alpha);
        }
    }

    alpha = 3.0;
    dgese(m - 1, n - 1, alpha, &A, 0, 1);
    for (int i = 0; i < m - 1; ++i) {
        for (int j = 1; j < n; ++j) {
            munit_assert_double(MATEL(&A, i, j), ==, alpha);
        }
    }

    return MUNIT_OK;
}

static MunitTest get_set_tests[] = {
        {"/dvecex_sp", test_dvecex_sp, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        {"/drowin", test_drowin, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        {"/drowex", test_drowex, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        {"/dcolin", test_dcolin, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        {"/dcolex", test_dcolex, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        {"/ddiaex", test_ddiaex, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        {"/dvecse", test_dvecse, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        {"/dgese", test_dgese, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        NULL_TEST};


/***************************************************************************************
 *  copy
 ***************************************************************************************/

MunitResult test_dgecp(const MunitParameter params[], void* fixture) {
    (void) params;
    (void) fixture;
    const int m = 5, n = 7;
    struct mat A, B;
    create_mat(m, n, &A, munit_malloc(memsize_mat(m, n)));
    init_random_mat(&A);
    create_mat(m, n, &B, munit_malloc(memsize_mat(m, n)));
    dgese(m, n, 0.0, &B, 0, 0);

    dgecp(m, n, &A, 0, 0, &B, 0, 0);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            munit_assert_double(MATEL(&A, i, j), ==, MATEL(&B, i, j));
        }
    }

    dgese(m, n, 0.0, &B, 0, 0);
    dgecp(m - 1, n - 2, &A, 1, 1, &B, 0, 0);
    for (int i = 0; i < m - 1; ++i) {
        for (int j = 0; j < n - 2; ++j) {
            munit_assert_double(MATEL(&A, i + 1, j + 1), ==, MATEL(&B, i, j));
        }
    }

    return MUNIT_OK;
}

MunitResult test_dtrcp(const MunitParameter params[], void* fixture) {
    (void) params;
    (void) fixture;
    const int m = 5, n = 7;
    struct mat A, B;
    create_mat(m, n, &A, munit_malloc(memsize_mat(m, n)));
    init_random_mat(&A);
    create_mat(m, n, &B, munit_malloc(memsize_mat(m, n)));
    dgese(m, n, 0.0, &B, 0, 0);

    dtrcp_l(m, &A, 0, 0, &B, 0, 0);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < i; ++j) {
            munit_assert_double(MATEL(&A, i, j), ==, MATEL(&B, i, j));
        }
    }

    dgese(m, n, 0.0, &B, 0, 0);
    dtrcp_l(m - 1, &A, 1, 0, &B, 0, 0);
    for (int i = 0; i < m - 1; ++i) {
        for (int j = 0; j <= i; ++j) {
            munit_assert_double(MATEL(&A, i + 1, j), ==, MATEL(&B, i, j));
        }
    }

    return MUNIT_OK;
}

static MunitTest copy_tests[] = {
        {"/dgecp", test_dgecp, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        {"/dtrcp", test_dtrcp, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        NULL_TEST,
};

/***************************************************************************************
 * transpositions
 ***************************************************************************************/

MunitResult test_dgetr(const MunitParameter params[], void* fixture) {
    (void) params;
    (void) fixture;
    const int m = 5, n = 7;
    double *Ap = munit_newa(double, n* m), *Bp = munit_newa(double, n* m);
    init_random_buffer(Ap, n * m);
    init_random_buffer(Bp, n * m);
    struct mat A, B;
    create_mat(m, n, &A, munit_malloc(memsize_mat(m, n)));
    pack_mat(m, n, Ap, m, &A, 0, 0);
    create_mat(m, n, &B, munit_malloc(memsize_mat(m, n)));
    dgese(m, n, 0.0, &B, 0, 0);

    dgetr(m, n, &A, 0, 0, &B, 0, 0);
    for (int i = 0; i < m; ++i) {
        int j = 0;
        for (; j < m; ++j) {
            munit_assert_double_equal(MATEL(&B, i, j), MATEL(&A, j, i), 10);
        }
        for (; j < n; ++j) {
            munit_assert_double_equal(MATEL(&B, i, j), 0.0, 10);
        }
    }

    dgese(m, n, 0.0, &B, 0, 0);
    dgetr(m, n - 2, &A, 0, 0, &B, 0, 0);
    for (int i = 0; i < m; ++i) {
        int j = 0;
        for (; j < m; ++j) {
            munit_assert_double_equal(MATEL(&B, i, j), MATEL(&A, j, i), 10);
        }
        for (; j < n; ++j) {
            munit_assert_double_equal(MATEL(&B, i, j), 0.0, 10);
        }
    }

    return MUNIT_OK;
}
MunitResult test_dtrtr_l(const MunitParameter params[], void* fixture) {
    (void) params;
    (void) fixture;
    const int m = 5, n = 7;
    double *Ap = munit_newa(double, n* m), *Bp = munit_newa(double, n* m);
    init_random_buffer(Ap, n * m);
    init_random_buffer(Bp, n * m);
    struct mat A, B;
    create_mat(m, n, &A, munit_malloc(memsize_mat(m, n)));
    pack_mat(m, n, Ap, m, &A, 0, 0);
    create_mat(m, n, &B, munit_malloc(memsize_mat(m, n)));
    dgese(m, n, 0.0, &B, 0, 0);

    dtrtr_l(m, &A, 0, 0, &B, 0, 0);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j <= i; ++j) {
            if (j < i || j >= m) {
                munit_assert_double_equal(MATEL(&B, i, j), 0.0, 10);
            } else {
                munit_assert_double_equal(MATEL(&B, i, j), MATEL(&A, j, i), 10);
            }
        }
    }

    return MUNIT_OK;
}

static MunitTest transposition_tests[] = {
        {"/dgetr", test_dgetr, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        {"/dtrtr_l", test_dtrtr_l, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        NULL_TEST,
};


/***************************************************************************************
 * extended blas level 1 routines
 ***************************************************************************************/

MunitResult test_dgead(const MunitParameter params[], void* fixture) {
    (void) params;
    (void) fixture;
    const int m = 5, n = 7;
    double *Ap = munit_newa(double, n* m), *Bp = munit_newa(double, n* m);
    init_random_buffer(Ap, n * m);
    init_random_buffer(Bp, n * m);
    struct mat A, B;
    create_mat(m, n, &A, munit_malloc(memsize_mat(m, n)));
    pack_mat(m, n, Ap, m, &A, 0, 0);
    create_mat(m, n, &B, munit_malloc(memsize_mat(m, n)));
    pack_mat(m, n, Bp, m, &B, 0, 0);

    dgead(m, n, 1.0, &A, 0, 0, &B, 0, 0);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            munit_assert_double_equal(MATEL(&B, i, j), Bp[i + j * m] + Ap[i + j * m], 10);
        }
    }

    pack_mat(m, n, Bp, m, &B, 0, 0);
    dgead(m - 2, n - 1, 1.0, &A, 1, 0, &B, 2, 1);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i < 2 || j < 1) {
                munit_assert_double_equal(MATEL(&B, i, j), Bp[i + j * m], 10);
            } else {
                munit_assert_double_equal(MATEL(&B, i, j), Bp[i + j * m] + Ap[i - 1 + (j - 1) * m], 10);
            }
        }
    }

    return MUNIT_OK;
}


MunitResult test_ddiare(const MunitParameter params[], void* fixture) {
    (void) params;
    (void) fixture;
    const int m = 5, n = 5;
    double* Ap = munit_newa(double, n* m);
    init_random_buffer(Ap, n * m);
    struct mat A;
    create_mat(m, n, &A, munit_malloc(memsize_mat(m, n)));
    pack_mat(m, n, Ap, m, &A, 0, 0);

    ddiare(m, 1.0, &A, 0, 0);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                munit_assert_double_equal(MATEL(&A, i, j), Ap[i + j * m] + 1.0, 10);
            } else {
                munit_assert_double_equal(MATEL(&A, i, j), Ap[i + j * m], 10);
            }
        }
    }

    ddiare(m, NAN, &A, 0, 0);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                munit_assert_isnan(MATEL(&A, i, i));
            } else {
                munit_assert_double_equal(MATEL(&A, i, j), Ap[i + j * m], 10);
            }
        }
    }
    return MUNIT_OK;
}


MunitResult test_ddiaad(const MunitParameter params[], void* fixture) {
    (void) params;
    (void) fixture;
    const int m = 5, n = 5;
    double* Ap = munit_newa(double, n* m);
    init_random_buffer(Ap, n * m);
    struct mat A;
    create_mat(m, n, &A, munit_malloc(memsize_mat(m, n)));
    pack_mat(m, n, Ap, m, &A, 0, 0);
    const int d = 3;
    struct vec v;
    create_vec(d, &v, munit_malloc(memsize_vec(d)));
    init_random_vec(&v);
    const int idx[] = {0, 2, 4};

    ddiaad_sp(d, 1.0, &v, 0, idx, &A, 0, 0);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j && (i == idx[0] || i == idx[1] || i == idx[2])) {
                munit_assert_double_equal(MATEL(&A, i, j), Ap[i + j * m] + VECEL(&v, i / 2), 10);
            } else {
                munit_assert_double_equal(MATEL(&A, i, j), Ap[i + j * m], 10);
            }
        }
    }

    pack_mat(m, n, Ap, m, &A, 0, 0);
    ddiaad_sp(d, 0.0, &v, 0, idx, &A, 0, 0);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            munit_assert_double_equal(MATEL(&A, i, j), Ap[i + j * m], 10);
        }
    }

    pack_mat(m, n, Ap, m, &A, 0, 0);
    ddiaad_sp(d, NAN, &v, 0, idx, &A, 0, 0);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j && (i == idx[0] || i == idx[1] || i == idx[2])) {
                munit_assert_isnan(MATEL(&A, i, j));
            } else {
                munit_assert_double_equal(MATEL(&A, i, j), Ap[i + j * m], 10);
            }
        }
    }
    return MUNIT_OK;
}

static MunitTest extended_blas1_tests[] = {
        {"/dgead", test_dgead, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        {"/ddiare", test_ddiare, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        {"/ddiad", test_ddiaad, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
        NULL_TEST,
};


/***************************************************************************************
 * swap
 ***************************************************************************************/


static MunitTest swap_tests[] = {NULL_TEST};

/***************************************************************************************
 * norm
 ***************************************************************************************/

static MunitTest norm_tests[] = {NULL_TEST};


/***************************************************************************************
 *  main
 ***************************************************************************************/

static MunitSuite all_suites[] = {
        // {"/dgese", dgese_tests, NULL, 1, MUNIT_SUITE_OPTION_NONE},
        // {"/dgecp", dgecp_tests, NULL, 1, MUNIT_SUITE_OPTION_NONE},
        // {"/dtrcp", dtrcp_tests, NULL, 1, MUNIT_SUITE_OPTION_NONE},
        // {"/dgead", dgead_tests, NULL, 1, MUNIT_SUITE_OPTION_NONE},
        // {"/dgetr", dgetr_tests, NULL, 1, MUNIT_SUITE_OPTION_NONE},
        // {"/dtrtr", dtrtr_tests, NULL, 1, MUNIT_SUITE_OPTION_NONE},
        // {"/ddiare", ddiare_tests, NULL, 1, MUNIT_SUITE_OPTION_NONE},
        // {"/ddiaex", ddiaex_tests, NULL, 1, MUNIT_SUITE_OPTION_NONE},
        // {"/ddiaad", ddiaad_tests, NULL, 1, MUNIT_SUITE_OPTION_NONE},

        {"/get_set", get_set_tests, NULL, 1, MUNIT_SUITE_OPTION_NONE},
        {"/copy", copy_tests, NULL, 1, MUNIT_SUITE_OPTION_NONE},
        {"/transposition", transposition_tests, NULL, 1, MUNIT_SUITE_OPTION_NONE},
        {"/swap", swap_tests, NULL, 1, MUNIT_SUITE_OPTION_NONE},
        {"/extended_blas1", extended_blas1_tests, NULL, 1, MUNIT_SUITE_OPTION_NONE},
        {"/norm", norm_tests, NULL, 1, MUNIT_SUITE_OPTION_NONE},
        NULL_SUITE,
};

static const MunitSuite file_suite = {"/blas/misc", NULL, all_suites, 1, MUNIT_SUITE_OPTION_NONE};

int main(int argc, char* argv[]) {
    return munit_suite_main(&file_suite, NULL, argc, argv);
}
