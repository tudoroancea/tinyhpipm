include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_d_cond test_d_cond.c d_tools.c)
add_executable(test_d_part_cond test_d_part_cond.c d_tools.c)
add_executable(test_d_dense test_d_dense.c)
add_executable(test_d_ocp test_d_ocp.c d_tools.c)

target_link_libraries(test_d_cond hpipm_dynamic blasfeo_dynamic m)
target_link_libraries(test_d_part_cond hpipm_dynamic blasfeo_dynamic m)
target_link_libraries(test_d_dense hpipm_dynamic blasfeo_dynamic m)
target_link_libraries(test_d_ocp hpipm_dynamic blasfeo_dynamic m)

add_test(NAME test_d_cond COMMAND test_d_cond)
add_test(NAME test_d_part_cond COMMAND test_d_part_cond)
add_test(NAME test_d_dense COMMAND test_d_dense)
add_test(NAME test_d_ocp COMMAND test_d_ocp)
