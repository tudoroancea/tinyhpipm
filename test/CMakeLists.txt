set(TEST_UTILS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/utils/munit.h)
add_library(test_utils utils/munit.h utils/munit.c utils/naive_blas.h utils/naive_blas.c)
target_link_libraries(test_utils m)

function(create_unit_test test_name)
    add_executable(${test_name} ${test_name}.c ${TEST_UTILS_HEADERS})
    target_link_libraries(${test_name} tinyhpipm_static test_utils m)
    add_test(NAME ${test_name} COMMAND ${test_name})
endfunction()

add_subdirectory(blas)
add_subdirectory(ocp)
