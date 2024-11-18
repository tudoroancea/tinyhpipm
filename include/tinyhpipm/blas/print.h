#ifndef HPIPM_BLAS_print_H
#define HPIPM_BLAS_print_H

#include "tinyhpipm/blas/struct.h"
#include <stdio.h>

// void print_mat_ref(int m, int n,  mat_ref* sA, int ai, int aj);
// void allocate_mat_ref(int m, int n,  mat_ref* sA);
// void allocate_vec_ref(int m,  vec_ref* sa);
// void free_mat_ref( mat_ref* sA);
// void free_vec_ref( vec_ref* sa);
// void print_exp_mat_ref(int m, int n,  mat_ref* sA, int ai, int aj);
// void print_to_file_mat_ref(FILE* file, int m, int n,  mat_ref* sA, int ai, int aj);
// void print_to_file_exp_mat_ref(FILE* file, int m, int n,  mat_ref* sA, int ai, int aj);
// void print_to_string_mat_ref(char** buf_out, int m, int n,  mat_ref* sA, int ai, int aj);

// void print_vec(int m, struct vec* sa, int ai);
// void print_exp_vec(int m, struct vec* sa, int ai);
// void print_to_file_vec(FILE* file, int m, struct vec* sa, int ai);
// void print_to_string_vec(char** buf_out, int m, struct vec* sa, int ai);
// void print_tran_vec(int m, struct vec* sa, int ai);
// void print_exp_tran_vec(int m, struct vec* sa, int ai);
// void print_to_file_tran_vec(FILE* file, int m, struct vec* sa, int ai);
// void print_to_string_tran_vec(char** buf_out, int m, struct vec* sa, int ai);

void print_mat(int m, int n, struct mat* sA, int ai, int aj);
// print in exponential notation a strmat
void print_exp_mat(int m, int n, struct mat* sA, int ai, int aj);
// print to file a strmat
void print_to_file_mat(FILE* file, int m, int n, struct mat* sA, int ai, int aj);
// print to file a strmat in exponential format
void print_to_file_exp_mat(FILE* file, int m, int n, struct mat* sA, int ai, int aj);
// print to string a strmat
void print_to_string_mat(char** buf_out, int m, int n, struct mat* sA, int ai, int aj);
// print the transposed of a strmat
void print_tran_mat(int m, int n, struct mat* sA, int ai, int aj);
// print a strvec
void print_vec(int m, struct vec* sa, int ai);
// print in exponential notation a strvec
void print_exp_vec(int m, struct vec* sa, int ai);
// print to file a strvec
void print_to_file_vec(FILE* file, int m, struct vec* sa, int ai);
// print to string a strvec
void print_to_string_vec(char** buf_out, int m, struct vec* sa, int ai);
// print the transposed of a strvec
void print_tran_vec(int m, struct vec* sa, int ai);
// print in exponential notation the transposed of a strvec
void print_exp_tran_vec(int m, struct vec* sa, int ai);
// print to file the transposed of a strvec
void print_to_file_tran_vec(FILE* file, int m, struct vec* sa, int ai);
// print to string the transposed of a strvec
void print_to_string_tran_vec(char** buf_out, int m, struct vec* sa, int ai);

// print an integer matrix in colmajor format
void int_print_mat(int row, int col, int* A, int lda);

#endif  // HPIPM_BLAS_print_H
