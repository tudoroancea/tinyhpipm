add_executable(test_struct test_struct.c ${MUNIT_HEADERS})
target_link_libraries(test_struct tinyhpipm_static munit)
add_test(NAME test_struct COMMAND test_struct)
