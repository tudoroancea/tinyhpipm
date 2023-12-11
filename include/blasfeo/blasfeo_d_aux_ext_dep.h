/*
 * auxiliary algebra operation external dependancies header
 *
 * include/blasfeo_d_aux_ext_dep.h
 *
 * - dynamic memory allocation
 * - print
 *
 */

#ifndef BLASFEO_D_AUX_EXT_DEP_H_
#define BLASFEO_D_AUX_EXT_DEP_H_


#include <stdio.h>


#include "blasfeo_common.h"


#ifdef __cplusplus
extern "C" {
#endif


#ifdef EXT_DEP

/* column-major matrices */

// dynamically allocate row*col doubles of memory and set accordingly a pointer to double; set allocated memory to zero
void d_zeros(double** pA, int row, int col);
// dynamically allocate row*col doubles of memory aligned to 64-byte boundaries and set accordingly a pointer to double; set allocated memory to zero
void d_zeros_align(double** pA, int row, int col);
// dynamically allocate size bytes of memory aligned to 64-byte boundaries and set accordingly a pointer to double; set allocated memory to zero
void d_zeros_align_bytes(double** pA, int size);
// free the memory allocated by d_zeros
void d_free(double* pA);
// free the memory allocated by d_zeros_align or d_zeros_align_bytes
void d_free_align(double* pA);
// print a column-major matrix
void d_print_mat(int m, int n, double* A, int lda);
// print the transposed of a column-major matrix
void d_print_tran_mat(int row, int col, double* A, int lda);
// print to file a column-major matrix
void d_print_to_file_mat(FILE* file, int row, int col, double* A, int lda);
// print to file a column-major matrix in exponential format
void d_print_to_file_exp_mat(FILE* file, int row, int col, double* A, int lda);
// print to string a column-major matrix
void d_print_to_string_mat(char** buf_out, int row, int col, double* A, int lda);
// print to file the transposed of a column-major matrix
void d_print_tran_to_file_mat(FILE* file, int row, int col, double* A, int lda);
// print to file the transposed of a column-major matrix in exponential format
void d_print_tran_to_file_exp_mat(FILE* file, int row, int col, double* A, int lda);
// print in exponential notation a column-major matrix
void d_print_exp_mat(int m, int n, double* A, int lda);
// print in exponential notation the transposed of a column-major matrix
void d_print_exp_tran_mat(int row, int col, double* A, int lda);

/* strmat and strvec */

// create a strmat for a matrix of size m*n by dynamically allocating memory
void blasfeo_allocate_dmat(int m, int n, struct blasfeo_dmat* sA);
// create a strvec for a vector of size m by dynamically allocating memory
void blasfeo_allocate_dvec(int m, struct blasfeo_dvec* sa);
// free the memory allocated by blasfeo_allocate_dmat
void blasfeo_free_dmat(struct blasfeo_dmat* sA);
// free the memory allocated by blasfeo_allocate_dvec
void blasfeo_free_dvec(struct blasfeo_dvec* sa);
// print a strmat
void blasfeo_print_dmat(int m, int n, struct blasfeo_dmat* sA, int ai, int aj);
// print in exponential notation a strmat
void blasfeo_print_exp_dmat(int m, int n, struct blasfeo_dmat* sA, int ai, int aj);
// print to file a strmat
void blasfeo_print_to_file_dmat(FILE* file, int m, int n, struct blasfeo_dmat* sA, int ai, int aj);
// print to file a strmat in exponential format
void blasfeo_print_to_file_exp_dmat(FILE* file, int m, int n, struct blasfeo_dmat* sA, int ai, int aj);
// print to string a strmat
void blasfeo_print_to_string_dmat(char** buf_out, int m, int n, struct blasfeo_dmat* sA, int ai, int aj);
// print the transposed of a strmat
void blasfeo_print_tran_dmat(int m, int n, struct blasfeo_dmat* sA, int ai, int aj);
// print a strvec
void blasfeo_print_dvec(int m, struct blasfeo_dvec* sa, int ai);
// print in exponential notation a strvec
void blasfeo_print_exp_dvec(int m, struct blasfeo_dvec* sa, int ai);
// print to file a strvec
void blasfeo_print_to_file_dvec(FILE* file, int m, struct blasfeo_dvec* sa, int ai);
// print to string a strvec
void blasfeo_print_to_string_dvec(char** buf_out, int m, struct blasfeo_dvec* sa, int ai);
// print the transposed of a strvec
void blasfeo_print_tran_dvec(int m, struct blasfeo_dvec* sa, int ai);
// print in exponential notation the transposed of a strvec
void blasfeo_print_exp_tran_dvec(int m, struct blasfeo_dvec* sa, int ai);
// print to file the transposed of a strvec
void blasfeo_print_to_file_tran_dvec(FILE* file, int m, struct blasfeo_dvec* sa, int ai);
// print to string the transposed of a strvec
void blasfeo_print_to_string_tran_dvec(char** buf_out, int m, struct blasfeo_dvec* sa, int ai);

#endif  // EXT_DEP


#ifdef __cplusplus
}
#endif

#endif  // BLASFEO_D_AUX_EXT_DEP_H_
