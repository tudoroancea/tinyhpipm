#include "../utils/munit.h"  // TODO: make it non-relative
#include "tinyhpipm/blas.h"
#include "tinyhpipm/ocp/d_ocp_qp.h"
#include "tinyhpipm/ocp/d_ocp_qp_dim.h"
#include "tinyhpipm/ocp/d_ocp_qp_ipm.h"
#include "tinyhpipm/ocp/d_ocp_qp_sol.h"
#include <stdio.h>
#include <stdlib.h>

#define NULL_TEST \
    { NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL }

#define NULL_SUITE \
    { NULL, NULL, NULL, 0, MUNIT_SUITE_OPTION_NONE }
#define PRINT 0


/***************************************************************************************
 *  main
 ***************************************************************************************/
static MunitSuite all_suites[] = {
        NULL_SUITE,
};

static const MunitSuite giga_suite = {"/ocp/qp", NULL, all_suites, 1, MUNIT_SUITE_OPTION_NONE};

int main(int argc, char* argv[]) {
    return munit_suite_main(&giga_suite, NULL, argc, argv);
}
