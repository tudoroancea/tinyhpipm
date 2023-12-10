#include <stdio.h>
#include <stdlib.h>

#include <blasfeo/blasfeo_block_size.h>
#include <blasfeo/blasfeo_common.h>
#include <blasfeo/blasfeo_d_aux.h>
#include <blasfeo/blasfeo_d_aux_ext_dep.h>
#include <blasfeo/blasfeo_stdlib.h>


#define CREATE_MAT blasfeo_create_dmat
#define CREATE_VEC blasfeo_create_dvec
#define MEMSIZE_MAT blasfeo_memsize_dmat
#define MEMSIZE_VEC blasfeo_memsize_dvec
#define REAL double
#define MAT blasfeo_dmat
#define MATEL BLASFEO_DMATEL
#define VEC blasfeo_dvec
#define VECEL BLASFEO_DVECEL

#define ALLOCATE_MAT blasfeo_allocate_dmat
#define ALLOCATE_VEC blasfeo_allocate_dvec
#define FREE_MAT blasfeo_free_dmat
#define FREE_VEC blasfeo_free_dvec
#define PRINT_MAT blasfeo_print_dmat
#define PRINT_TRAN_MAT blasfeo_print_tran_dmat
#define PRINT_VEC blasfeo_print_dvec
#define PRINT_TRAN_VEC blasfeo_print_tran_dvec
#define PRINT_TO_FILE_MAT blasfeo_print_to_file_dmat
#define PRINT_TO_FILE_EXP_MAT blasfeo_print_to_file_exp_dmat
#define PRINT_TO_FILE_VEC blasfeo_print_to_file_dvec
#define PRINT_TO_FILE_TRAN_VEC blasfeo_print_to_file_tran_dvec
#define PRINT_TO_STRING_MAT blasfeo_print_to_string_dmat
#define PRINT_TO_STRING_VEC blasfeo_print_to_string_dvec
#define PRINT_TO_STRING_TRAN_VEC blasfeo_print_to_string_tran_dvec
#define PRINT_EXP_MAT blasfeo_print_exp_dmat
#define PRINT_EXP_TRAN_MAT blasfeo_print_exp_tran_dmat
#define PRINT_EXP_VEC blasfeo_print_exp_dvec
#define PRINT_EXP_TRAN_VEC blasfeo_print_exp_tran_dvec

#include "x_aux_ext_dep.c"
