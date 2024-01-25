cmake_minimum_required(VERSION 3.5)
project(tinyhpipm C ASM)

# #################################################################################
# general options ################################################################
# #################################################################################
set(CMAKE_C_STANDARD 11)

# set default build type and generator
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug, Release" FORCE)
endif()

# check that the compiler is supported
if(NOT CMAKE_C_COMPILER_ID MATCHES "GNU" AND NOT CMAKE_C_COMPILER_ID MATCHES "Clang")
	message(FATAL_ERROR "The compiler ${CMAKE_C_COMPILER_ID} is not supported, only GCC and Clang are supported")
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # create a compile_commands.json file

set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/install" CACHE STRING "Installation path" FORCE)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mcpu=native")

# #################################################################################
# create targets and install them #################################################
# #################################################################################
file(GLOB_RECURSE TINYHPIPM_HEADERS "include/tinyhpipm/*.h")
file(GLOB_RECURSE TINYHPIPM_SRC "src/*.c")
add_library(tinyhpipm_static STATIC)
target_sources(tinyhpipm_static PRIVATE ${TINYHPIPM_SRC} PUBLIC FILE_SET HEADERS BASE_DIRS ${PROJECT_SOURCE_DIR}/include FILES ${TINYHPIPM_HEADERS})
target_link_libraries(tinyhpipm_static m)
add_library(tinyhpipm_dynamic SHARED)
target_sources(tinyhpipm_dynamic PRIVATE ${TINYHPIPM_SRC} PUBLIC FILE_SET HEADERS BASE_DIRS ${PROJECT_SOURCE_DIR}/include FILES ${TINYHPIPM_HEADERS})
target_link_libraries(tinyhpipm_dynamic m)

install(TARGETS tinyhpipm_static tinyhpipm_dynamic
	EXPORT ${PROJECT_NAME}-config
	FILE_SET HEADERS
	LIBRARY DESTINATION lib)
install(EXPORT ${PROJECT_NAME}-config NAMESPACE ${PROJECT_NAME}:: DESTINATION lib/cmake)

# #################################################################################
# testing #########################################################################
# #################################################################################

# tests
# enable_testing()
# add_subdirectory(test/blasfeo)
# add_subdirectory(test/hpipm)
