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

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # create a compile_commands.json file

set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/install" CACHE STRING "Installation path" FORCE)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC") 
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -g -DDEBUG -Wall -Wextra")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -DNDEBUG")


# #################################################################################
# target selection ###############################################################
# #################################################################################
set(TARGET "GENERIC" CACHE STRING "Set CPU architecture target")
set_property(CACHE TARGET PROPERTY STRINGS GENERIC X86_SSE3 X86_AVX X64_SSE3 X64_AVX X64_AVX_FMA X64_AVX2_FMA X64_AV512_FMA ARMV7A_NEON_VPFV3 ARMV8A_NEON2_VPFV4)

if(${TARGET} MATCHES X86_SSE3)
  set(TARGET_NEED_FEATURE_SSE3 1)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m32 -msse3")
elseif(${TARGET} MATCHES X86_AVX)
	set(TARGET_NEED_FEATURE_AVX 1)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m32 -mavx")
elseif(${TARGET} MATCHES X64_SSE3)
	set(TARGET_NEED_FEATURE_SSE3 1)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m64 -msse3")
elseif(${TARGET} MATCHES X64_AVX)
	set(TARGET_NEED_FEATURE_AVX 1)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m64 -mavx")
elseif(${TARGET} MATCHES X64_AVX_FMA)
	set(TARGET_NEED_FEATURE_AVX 1)
	set(TARGET_NEED_FEATURE_FMA 1)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m64 -mavx -mfma")
elseif(${TARGET} MATCHES X64_AVX2_FMA)
	set(TARGET_NEED_FEATURE_AVX2 1)
	set(TARGET_NEED_FEATURE_FMA 1)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m64 -mavx2 -mfma")
elseif(${TARGET} MATCHES X64_AV512_FMA)
	set(TARGET_NEED_FEATURE_AVX512 1)
	set(TARGET_NEED_FEATURE_FMA 1)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m64 -mavx512f -mavx512vl -mfma")
elseif(${TARGET} MATCHES ARMV7A_NEON_VPFV3)
	set(TARGET_NEED_FEATURE_VFPv3 1)
	set(TARGET_NEED_FEATURE_NEON 1)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -marm -mfloat-abi=hard -mfpu=neon-vfpv3")
elseif(${TARGET} MATCHES ARMV8A_NEON2_VPFV4)
	set(TARGET_NEED_FEATURE_VFPv4 1)
	set(TARGET_NEED_FEATURE_NEONv2 1)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=armv8-a+crc+crypto+simd")
endif()

# #################################################################################
# blasfeo ########################################################################
# #################################################################################
set(K_MAX_STACK 300 CACHE STRING "Maximum k value using stack memory")
add_compile_definitions(K_MAX_STACK=${K_MAX_STACK})
add_compile_definitions(LA_HIGH_PERFORMANCE)
add_compile_definitions(TARGET_${TARGET})
add_compile_definitions(EXT_DEP)

# Create the target.h file with the proper target definition
set(LA HIGH_PERFORMANCE)
set(MF PANELMAJ)
configure_file(${PROJECT_SOURCE_DIR}/blasfeo_target.h.in
	${CMAKE_CURRENT_SOURCE_DIR}/include/blasfeo/blasfeo_target.h @ONLY)

# compiler and assembly flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DTARGET_${TARGET}")
set(CMAKE_ASM_FLAGS "${CMAKE_ASM_FLAGS} -DTARGET_${TARGET}")

if(CMAKE_C_COMPILER_ID MATCHES "GNU" OR CMAKE_C_COMPILER_ID MATCHES "Clang")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_FLAGS_TARGET_${TARGET}}")
	set(CMAKE_ASM_FLAGS "${CMAKE_ASM_FLAGS} ${ASM_FLAGS_TARGET_${TARGET}}")
endif()

# blasfeo source files
file(GLOB AUX_COMMON_SRC
	${PROJECT_SOURCE_DIR}/src/blasfeo/auxiliary/blasfeo_processor_features.c
	${PROJECT_SOURCE_DIR}/src/blasfeo/auxiliary/blasfeo_stdlib.c
	${PROJECT_SOURCE_DIR}/src/blasfeo/auxiliary/memory.c
	${PROJECT_SOURCE_DIR}/src/blasfeo/auxiliary/d_aux_common.c
	${PROJECT_SOURCE_DIR}/src/blasfeo/auxiliary/s_aux_common.c
)

file(GLOB AUX_EXT_DEP_SRC
	${PROJECT_SOURCE_DIR}/src/blasfeo/auxiliary/d_aux_ext_dep_common.c
	${PROJECT_SOURCE_DIR}/src/blasfeo/auxiliary/s_aux_ext_dep_common.c
	${PROJECT_SOURCE_DIR}/src/blasfeo/auxiliary/d_aux_ext_dep.c
	${PROJECT_SOURCE_DIR}/src/blasfeo/auxiliary/s_aux_ext_dep.c
	${PROJECT_SOURCE_DIR}/src/blasfeo/auxiliary/v_aux_ext_dep_lib.c
	${PROJECT_SOURCE_DIR}/src/blasfeo/auxiliary/i_aux_ext_dep_lib.c
	${PROJECT_SOURCE_DIR}/src/blasfeo/auxiliary/timing.c
)

if(${TARGET} MATCHES X64_AVX512_FMA OR ${TARGET} MATCHES X64_AVX2_FMA OR ${TARGET} MATCHES X64_AVX)
	file(GLOB BLASFEO_HP_PM_SRC
		${PROJECT_SOURCE_DIR}/src/blasfeo/blasfeo_hp_pm/d_blas1_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/blasfeo_hp_pm/d_blas2_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/blasfeo_hp_pm/d_blas2_diag_lib.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/blasfeo_hp_pm/d_blas3_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/blasfeo_hp_pm/d_blas3_diag_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/blasfeo_hp_pm/d_lapack_lib4.c

		${PROJECT_SOURCE_DIR}/src/blasfeo/blasfeo_hp_pm/s_blas1_lib8.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/blasfeo_hp_pm/s_blas2_lib8.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/blasfeo_hp_pm/s_blas2_diag_lib.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/blasfeo_hp_pm/s_blas3_lib8.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/blasfeo_hp_pm/s_blas3_diag_lib8.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/blasfeo_hp_pm/s_lapack_lib8.c
	)

	file(GLOB AUX_HP_PM_SRC
		${PROJECT_SOURCE_DIR}/src/blasfeo/auxiliary/d_aux_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/auxiliary/s_aux_lib8.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/auxiliary/m_aux_lib48.c
	)
elseif(${TARGET} MATCHES X64_SSE3 OR ${TARGET} MATCHES X64_AVX_FMA OR ${TARGET} MATCHES ARMV8A_NEON2_VPFV4 OR ${TARGET} MATCHES ARMV7A_NEON_VPFV3 OR ${TARGET} MATCHES GENERIC)
	file(GLOB BLASFEO_HP_PM_SRC
		${PROJECT_SOURCE_DIR}/src/blasfeo/blasfeo_hp_pm/d_blas1_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/blasfeo_hp_pm/d_blas2_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/blasfeo_hp_pm/d_blas2_diag_lib.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/blasfeo_hp_pm/d_blas3_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/blasfeo_hp_pm/d_blas3_diag_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/blasfeo_hp_pm/d_lapack_lib4.c

		${PROJECT_SOURCE_DIR}/src/blasfeo/blasfeo_hp_pm/s_blas1_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/blasfeo_hp_pm/s_blas2_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/blasfeo_hp_pm/s_blas2_diag_lib.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/blasfeo_hp_pm/s_blas3_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/blasfeo_hp_pm/s_blas3_diag_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/blasfeo_hp_pm/s_lapack_lib4.c
	)

	file(GLOB AUX_HP_PM_SRC
		${PROJECT_SOURCE_DIR}/src/blasfeo/auxiliary/d_aux_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/auxiliary/s_aux_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/auxiliary/m_aux_lib44.c
	)
endif()

if(${TARGET} MATCHES X64_AVX2_FMA)
	file(GLOB KERNEL_SRC
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx2/kernel_dgemm_12x4_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx2/kernel_dgemm_8x8_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx2/kernel_dgemm_8x4_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx2/kernel_dgemm_4x4_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx2/kernel_dgemv_8_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx2/kernel_dgemv_4_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx2/kernel_dsymv_6_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx2/kernel_dger_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx2/kernel_dgetrf_pivot_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx2/kernel_dgebp_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx2/kernel_dgelqf_4_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx2/kernel_dgetr_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_dgeqrf_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_dgemm_diag_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_dgecp_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_dpack_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_dgetr_lib.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgemv_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dsymv_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dpack_buffer_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dger_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_ddot_lib.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_daxpy_lib.c

		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx2/kernel_sgemm_24x4_lib8.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx2/kernel_sgemm_16x4_lib8.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx2/kernel_sgemm_8x8_lib8.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx2/kernel_sgemm_8x4_lib8.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_sgemm_diag_lib8.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_sgemv_8_lib8.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_sgemv_4_lib8.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_sgecpsc_lib8.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_sgead_lib8.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_sgetr_lib8.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_spack_lib8.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgemm_8x4_lib8.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgemm_4x4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_spack_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sdot_lib.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_saxpy_lib.c

		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/sse3/kernel_align_x64.S

		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_d_aux_lib.S
	)
elseif(${TARGET} MATCHES X64_AVX)
	file(GLOB KERNEL_SRC
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_dgemm_12x4_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_dgemm_8x4_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_dgemm_4x4_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_dgemm_diag_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_dgemv_12_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_dgemv_8_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_dgemv_4_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_dsymv_6_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_dger_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_dgetrf_pivot_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_dgeqrf_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_dgebp_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_dgecp_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_dgetr_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_dpack_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_dgetr_lib.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgemv_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dsymv_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dpack_buffer_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dger_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_ddot_lib.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_daxpy_lib.c

		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_sgemm_16x4_lib8.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_sgemm_8x8_lib8.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_sgemm_8x4_lib8.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_sgemm_diag_lib8.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_sgemv_8_lib8.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_sgemv_4_lib8.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_sgecpsc_lib8.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_sgetr_lib8.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_sgead_lib8.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_spack_lib8.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgemm_8x4_lib8.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgemm_4x4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_spack_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sdot_lib.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_saxpy_lib.c

		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/sse3/kernel_align_x64.S

		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/avx/kernel_d_aux_lib.S
	)
elseif(${TARGET} MATCHES X64_SSE3)
	file(GLOB KERNEL_SRC
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/sse3/kernel_dgemm_4x4_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/sse3/kernel_dgemv_4_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgemm_4x4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgemm_diag_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgemv_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dsymv_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgetrf_pivot_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgeqrf_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgecp_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgetr_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dpack_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dpack_buffer_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dger_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_ddot_lib.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_daxpy_lib.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgetr_lib.c

		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/sse3/kernel_sgemm_4x4_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgemm_4x4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgemm_diag_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgemv_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_ssymv_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgetrf_pivot_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgecp_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgetr_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_spack_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sdot_lib.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_saxpy_lib.c

		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/sse3/kernel_align_x64.S

		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_d_aux_lib.c
	)
elseif(${TARGET} MATCHES X64_AVX_FMA)
	file(GLOB KERNEL_SRC
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/fma/kernel_dgemm_4x4_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgemm_4x4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgemm_diag_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgemv_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dsymv_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgetrf_pivot_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgeqrf_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgecp_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgetr_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dpack_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dpack_buffer_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dger_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_ddot_lib.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_daxpy_lib.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgetr_lib.c

		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgemm_4x4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgemm_diag_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgemv_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_ssymv_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgetrf_pivot_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgecp_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgetr_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_spack_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sdot_lib.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_saxpy_lib.c

		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/sse3/kernel_align_x64.S

		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_d_aux_lib.c
	)
elseif(${TARGET} MATCHES ARMV8A_NEON2_VPFV4)
	file(GLOB KERNEL_SRC
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/armv8a/kernel_dgemm_12x4_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/armv8a/kernel_dgemm_8x4_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/armv8a/kernel_dgemm_4x4_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/armv8a/kernel_dpack_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/armv8a/kernel_dgetrf_pivot_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/armv8a/kernel_dgemv_4_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/armv8a/kernel_dgetr_lib.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/armv8a/kernel_dger_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgemm_4x4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgemm_diag_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgemv_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dsymv_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgetrf_pivot_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgeqrf_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgecp_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgetr_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dpack_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dpack_buffer_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dger_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_ddot_lib.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_daxpy_lib.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgetr_lib.c

		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/armv8a/kernel_sgemm_16x4_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/armv8a/kernel_sgemm_12x4_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/armv8a/kernel_sgemm_8x8_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/armv8a/kernel_sgemm_8x4_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/armv8a/kernel_sgemm_4x4_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/armv8a/kernel_sgemv_4_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/armv8a/kernel_spack_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgemm_4x4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgemm_diag_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgemv_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_ssymv_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgetrf_pivot_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgecp_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgetr_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_spack_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sdot_lib.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_saxpy_lib.c

		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_align_generic.c

		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/armv8a/kernel_d_aux_lib.S
	)
elseif(${TARGET} MATCHES ARMV7A_NEON_VPFV3)
	file(GLOB KERNEL_SRC
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/armv7a/kernel_dgemm_4x4_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgemm_4x4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgemm_diag_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgemv_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dsymv_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgetrf_pivot_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgeqrf_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgecp_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgetr_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dpack_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dpack_buffer_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dger_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_ddot_lib.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_daxpy_lib.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgetr_lib.c

		# TODO: check if we should add ${PROJECT_SOURCE_DIR}/kernel/armv7a/kernel_sgemm_12x4_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/armv7a/kernel_sgemm_8x4_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/armv7a/kernel_sgemm_4x4_lib4.S
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgemm_4x4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgemm_diag_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgemv_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_ssymv_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgetrf_pivot_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgecp_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgetr_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_spack_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sdot_lib.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_saxpy_lib.c

		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_align_generic.c

		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_d_aux_lib.c
	)
elseif(${TARGET} MATCHES GENERIC)
	file(GLOB KERNEL_SRC
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgemm_4x4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgemm_diag_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgemv_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dsymv_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgecp_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgetr_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgetrf_pivot_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgeqrf_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dpack_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dpack_buffer_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dger_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_ddot_lib.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_daxpy_lib.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_dgetr_lib.c

		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgemm_4x4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgemm_diag_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgemv_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_ssymv_4_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgetrf_pivot_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgecp_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sgetr_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_spack_lib4.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_sdot_lib.c
		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_saxpy_lib.c

		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_align_generic.c

		${PROJECT_SOURCE_DIR}/src/blasfeo/kernel/generic/kernel_d_aux_lib.c
	)
endif()

set(BLASFEO_SRC ${AUX_COMMON_SRC} ${KERNEL_SRC} ${AUX_HP_PM_SRC} ${BLASFEO_HP_PM_SRC} ${AUX_EXT_DEP_SRC})
file(GLOB_RECURSE BLASFEO_HEADERS "include/blasfeo/*.h")

# #################################################################################
# HPIPM ##########################################################################
# #################################################################################
file(GLOB HPIPM_COND_SRC
	${PROJECT_SOURCE_DIR}/src/hpipm/cond/d_cond.c
	${PROJECT_SOURCE_DIR}/src/hpipm/cond/d_cond_aux.c
	${PROJECT_SOURCE_DIR}/src/hpipm/cond/d_part_cond.c
	${PROJECT_SOURCE_DIR}/src/hpipm/cond/d_cond_qcqp.c
	${PROJECT_SOURCE_DIR}/src/hpipm/cond/d_cast_qcqp.c
	${PROJECT_SOURCE_DIR}/src/hpipm/cond/d_part_cond_qcqp.c

	${PROJECT_SOURCE_DIR}/src/hpipm/cond/s_cond.c
	${PROJECT_SOURCE_DIR}/src/hpipm/cond/s_cond_aux.c
	${PROJECT_SOURCE_DIR}/src/hpipm/cond/s_part_cond.c
	${PROJECT_SOURCE_DIR}/src/hpipm/cond/s_cond_qcqp.c
	${PROJECT_SOURCE_DIR}/src/hpipm/cond/s_cast_qcqp.c
	${PROJECT_SOURCE_DIR}/src/hpipm/cond/s_part_cond_qcqp.c
)

file(GLOB HPIPM_DENSE_SRC
	${PROJECT_SOURCE_DIR}/src/hpipm/dense/d_dense_qp.c
	${PROJECT_SOURCE_DIR}/src/hpipm/dense/d_dense_qp_dim.c
	${PROJECT_SOURCE_DIR}/src/hpipm/dense/d_dense_qp_sol.c
	${PROJECT_SOURCE_DIR}/src/hpipm/dense/d_dense_qp_kkt.c
	${PROJECT_SOURCE_DIR}/src/hpipm/dense/d_dense_qp_ipm.c
	${PROJECT_SOURCE_DIR}/src/hpipm/dense/d_dense_qp_res.c
	${PROJECT_SOURCE_DIR}/src/hpipm/dense/d_dense_qp_utils.c
	${PROJECT_SOURCE_DIR}/src/hpipm/dense/d_dense_qcqp_dim.c
	${PROJECT_SOURCE_DIR}/src/hpipm/dense/d_dense_qcqp.c
	${PROJECT_SOURCE_DIR}/src/hpipm/dense/d_dense_qcqp_sol.c
	${PROJECT_SOURCE_DIR}/src/hpipm/dense/d_dense_qcqp_res.c
	${PROJECT_SOURCE_DIR}/src/hpipm/dense/d_dense_qcqp_ipm.c
	${PROJECT_SOURCE_DIR}/src/hpipm/dense/d_dense_qcqp_utils.c

	${PROJECT_SOURCE_DIR}/src/hpipm/dense/s_dense_qp.c
	${PROJECT_SOURCE_DIR}/src/hpipm/dense/s_dense_qp_dim.c
	${PROJECT_SOURCE_DIR}/src/hpipm/dense/s_dense_qp_sol.c
	${PROJECT_SOURCE_DIR}/src/hpipm/dense/s_dense_qp_kkt.c
	${PROJECT_SOURCE_DIR}/src/hpipm/dense/s_dense_qp_ipm.c
	${PROJECT_SOURCE_DIR}/src/hpipm/dense/s_dense_qp_res.c
	${PROJECT_SOURCE_DIR}/src/hpipm/dense/s_dense_qp_utils.c
	${PROJECT_SOURCE_DIR}/src/hpipm/dense/s_dense_qcqp_dim.c
	${PROJECT_SOURCE_DIR}/src/hpipm/dense/s_dense_qcqp.c
	${PROJECT_SOURCE_DIR}/src/hpipm/dense/s_dense_qcqp_sol.c
	${PROJECT_SOURCE_DIR}/src/hpipm/dense/s_dense_qcqp_res.c
	${PROJECT_SOURCE_DIR}/src/hpipm/dense/s_dense_qcqp_ipm.c
	${PROJECT_SOURCE_DIR}/src/hpipm/dense/s_dense_qcqp_utils.c
)

if(${TARGET} MATCHES X64_AVX_FMA OR ${TARGET} MATCHES X64_AVX2_FMA OR ${TARGET} MATCHES X64_AVX512_FMA)
	file(GLOB HPIPM_IPM_CORE_SRC
		${PROJECT_SOURCE_DIR}/src/hpipm/ipm_core/d_core_qp_ipm_aux_avx.c
		${PROJECT_SOURCE_DIR}/src/hpipm/ipm_core/d_core_qp_ipm.c
		${PROJECT_SOURCE_DIR}/src/hpipm/ipm_core/s_core_qp_ipm_aux_avx.c
		${PROJECT_SOURCE_DIR}/src/hpipm/ipm_core/s_core_qp_ipm.c
	)
elseif(${TARGET} MATCHES ARMV8A_NEON2_VPFV4)
	file(GLOB HPIPM_IPM_CORE_SRC
		${PROJECT_SOURCE_DIR}/src/hpipm/ipm_core/d_core_qp_ipm_aux_neon.c
		# ${PROJECT_SOURCE_DIR}/src/hpipm/ipm_core/d_core_qp_ipm_aux.c
		${PROJECT_SOURCE_DIR}/src/hpipm/ipm_core/d_core_qp_ipm.c
		${PROJECT_SOURCE_DIR}/src/hpipm/ipm_core/s_core_qp_ipm_aux.c
		${PROJECT_SOURCE_DIR}/src/hpipm/ipm_core/s_core_qp_ipm.c
	)
else()
	file(GLOB HPIPM_IPM_CORE_SRC
		${PROJECT_SOURCE_DIR}/src/hpipm/ipm_core/d_core_qp_ipm_aux.c
		${PROJECT_SOURCE_DIR}/src/hpipm/ipm_core/d_core_qp_ipm.c
		${PROJECT_SOURCE_DIR}/src/hpipm/ipm_core/s_core_qp_ipm_aux.c
		${PROJECT_SOURCE_DIR}/src/hpipm/ipm_core/s_core_qp_ipm.c
	)
endif()

file(GLOB HPIPM_OCP_SRC
	${PROJECT_SOURCE_DIR}/src/hpipm/ocp/d_ocp_qp.c
	${PROJECT_SOURCE_DIR}/src/hpipm/ocp/d_ocp_qp_dim.c
	${PROJECT_SOURCE_DIR}/src/hpipm/ocp/d_ocp_qp_sol.c
	${PROJECT_SOURCE_DIR}/src/hpipm/ocp/d_ocp_qp_kkt.c
	${PROJECT_SOURCE_DIR}/src/hpipm/ocp/d_ocp_qp_ipm.c
	${PROJECT_SOURCE_DIR}/src/hpipm/ocp/d_ocp_qp_res.c
	${PROJECT_SOURCE_DIR}/src/hpipm/ocp/d_ocp_qp_utils.c
	${PROJECT_SOURCE_DIR}/src/hpipm/ocp/d_ocp_qp_red.c
	${PROJECT_SOURCE_DIR}/src/hpipm/ocp/d_ocp_qcqp.c
	${PROJECT_SOURCE_DIR}/src/hpipm/ocp/d_ocp_qcqp_dim.c
	${PROJECT_SOURCE_DIR}/src/hpipm/ocp/d_ocp_qcqp_sol.c
	${PROJECT_SOURCE_DIR}/src/hpipm/ocp/d_ocp_qcqp_ipm.c
	${PROJECT_SOURCE_DIR}/src/hpipm/ocp/d_ocp_qcqp_res.c
	${PROJECT_SOURCE_DIR}/src/hpipm/ocp/d_ocp_qcqp_utils.c

	${PROJECT_SOURCE_DIR}/src/hpipm/ocp/s_ocp_qp.c
	${PROJECT_SOURCE_DIR}/src/hpipm/ocp/s_ocp_qp_dim.c
	${PROJECT_SOURCE_DIR}/src/hpipm/ocp/s_ocp_qp_sol.c
	${PROJECT_SOURCE_DIR}/src/hpipm/ocp/s_ocp_qp_kkt.c
	${PROJECT_SOURCE_DIR}/src/hpipm/ocp/s_ocp_qp_ipm.c
	${PROJECT_SOURCE_DIR}/src/hpipm/ocp/s_ocp_qp_res.c
	${PROJECT_SOURCE_DIR}/src/hpipm/ocp/s_ocp_qp_utils.c
	${PROJECT_SOURCE_DIR}/src/hpipm/ocp/s_ocp_qp_red.c
	${PROJECT_SOURCE_DIR}/src/hpipm/ocp/s_ocp_qcqp.c
	${PROJECT_SOURCE_DIR}/src/hpipm/ocp/s_ocp_qcqp_dim.c
	${PROJECT_SOURCE_DIR}/src/hpipm/ocp/s_ocp_qcqp_sol.c
	${PROJECT_SOURCE_DIR}/src/hpipm/ocp/s_ocp_qcqp_ipm.c
	${PROJECT_SOURCE_DIR}/src/hpipm/ocp/s_ocp_qcqp_res.c
	${PROJECT_SOURCE_DIR}/src/hpipm/ocp/s_ocp_qcqp_utils.c

	# ${PROJECT_SOURCE_DIR}/ocp_qp/m_ocp_qp.c
	# ${PROJECT_SOURCE_DIR}/ocp_qp/m_ocp_qp_kkt.c
	# ${PROJECT_SOURCE_DIR}/ocp_qp/m_ocp_qp_ipm_hard.c
)

file(GLOB HPIPM_AUXILIARY_SRC
	${PROJECT_SOURCE_DIR}/src/hpipm/auxiliary/aux_string.c
	${PROJECT_SOURCE_DIR}/src/hpipm/auxiliary/aux_mem.c
	${PROJECT_SOURCE_DIR}/src/hpipm/auxiliary/timing.c
)

file(GLOB HPIPM_SIM_SRC
	${PROJECT_SOURCE_DIR}/src/hpipm/sim_core/d_sim_rk.c
	${PROJECT_SOURCE_DIR}/src/hpipm/sim_core/d_sim_erk.c

	${PROJECT_SOURCE_DIR}/src/hpipm/sim_core/s_sim_rk.c
	${PROJECT_SOURCE_DIR}/src/hpipm/sim_core/s_sim_erk.c
)

set(HPIPM_SRC ${HPIPM_COND_SRC} ${HPIPM_DENSE_SRC} ${HPIPM_IPM_CORE_SRC} ${HPIPM_OCP_SRC} ${TREE_HPIPM_OCP_SRC} ${HPIPM_AUXILIARY_SRC} ${HPIPM_SIM_SRC})
file(GLOB_RECURSE HPIPM_HEADERS "include/hpipm/*.h")

# #################################################################################
# create targets and install them #################################################
# #################################################################################
add_library(blasfeo_static STATIC)
target_sources(blasfeo_static PRIVATE ${BLASFEO_SRC} PUBLIC FILE_SET HEADERS BASE_DIRS ${PROJECT_SOURCE_DIR}/include FILES ${BLASFEO_HEADERS})
target_link_libraries(blasfeo_static m)  # necessary on Linux
add_library(blasfeo_dynamic SHARED)
target_sources(blasfeo_dynamic PRIVATE ${BLASFEO_SRC} PUBLIC FILE_SET HEADERS BASE_DIRS ${PROJECT_SOURCE_DIR}/include FILES ${BLASFEO_HEADERS})
target_link_libraries(blasfeo_dynamic m)  # necessary on Linux

add_library(hpipm_static STATIC)
target_sources(hpipm_static PRIVATE ${HPIPM_SRC} PUBLIC FILE_SET HEADERS BASE_DIRS ${PROJECT_SOURCE_DIR}/include FILES ${HPIPM_HEADERS})
target_link_libraries(hpipm_static blasfeo_static)
add_library(hpipm_dynamic SHARED)
target_sources(hpipm_dynamic PRIVATE ${HPIPM_SRC} PUBLIC FILE_SET HEADERS BASE_DIRS ${PROJECT_SOURCE_DIR}/include FILES ${HPIPM_HEADERS})
target_link_libraries(hpipm_dynamic blasfeo_dynamic)

install(TARGETS blasfeo_static blasfeo_dynamic hpipm_static hpipm_dynamic
	EXPORT ${PROJECT_NAME}-config
	FILE_SET HEADERS
	LIBRARY DESTINATION lib)
install(EXPORT ${PROJECT_NAME}-config NAMESPACE ${PROJECT_NAME}:: DESTINATION lib/cmake)

# #################################################################################
# testing #########################################################################
# #################################################################################

# tests
enable_testing()
add_subdirectory(test/blasfeo)
add_subdirectory(test/hpipm)

# #################################################################################
# new #############################################################################
# #################################################################################

# add_library(new_double STATIC)
# target_sources(new_double)
# what do we do now?
