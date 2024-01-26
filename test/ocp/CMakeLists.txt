# create_unit_test(test_ocp_qcqp)
# create_unit_test(test_ocp_qp)

add_executable(test_ocp_qcqp_safety_filter generic_ocp_qcqp_test.c data/safety_filter.c ${TEST_UTILS_HEADERS})
target_link_libraries(test_ocp_qcqp_safety_filter tinyhpipm_static test_utils m)
target_compile_definitions(test_ocp_qcqp_safety_filter PUBLIC TEST_NAME="/safety_filter")
add_test(NAME test_ocp_qcqp_safety_filter COMMAND test_ocp_qcqp_safety_filter)
