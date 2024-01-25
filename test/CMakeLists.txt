set(MUNIT_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/munit.h)
add_library(munit munit.h munit.c)

add_subdirectory(blas)
