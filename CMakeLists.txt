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

# set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -g -DDEBUG -Wall -Wextra")
# set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -DNDEBUG")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mcpu=native")

# #################################################################################
# create targets and install them #################################################
# #################################################################################
# add_library(blasfeo_static STATIC)
# target_sources(blasfeo_static PRIVATE ${BLASFEO_SRC} PUBLIC FILE_SET HEADERS BASE_DIRS ${PROJECT_SOURCE_DIR}/include FILES ${BLASFEO_HEADERS})
# target_link_libraries(blasfeo_static m)  # necessary on Linux
# add_library(blasfeo_dynamic SHARED)
# target_sources(blasfeo_dynamic PRIVATE ${BLASFEO_SRC} PUBLIC FILE_SET HEADERS BASE_DIRS ${PROJECT_SOURCE_DIR}/include FILES ${BLASFEO_HEADERS})
# target_link_libraries(blasfeo_dynamic m)  # necessary on Linux
file(GLOB_RECURSE TINYHPIPM_HEADERS "include/hpipm/*.h")
file(GLOB_RECURSE TINYHPIPM_SRC "src/*.c")
message(STATUS "TINYHPIPM_SRC: ${TINYHPIPM_SRC}")
add_library(hpipm_static STATIC)
target_sources(hpipm_static PRIVATE ${TINYHPIPM_SRC} PUBLIC FILE_SET HEADERS BASE_DIRS ${PROJECT_SOURCE_DIR}/include FILES ${TINYHPIPM_HEADERS})
target_link_libraries(hpipm_static m)
add_library(hpipm_dynamic SHARED)
target_sources(hpipm_dynamic PRIVATE ${TINYHPIPM_SRC} PUBLIC FILE_SET HEADERS BASE_DIRS ${PROJECT_SOURCE_DIR}/include FILES ${TINYHPIPM_HEADERS})
target_link_libraries(hpipm_dynamic m)

install(TARGETS hpipm_static hpipm_dynamic
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
