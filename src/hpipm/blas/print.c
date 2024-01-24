
#include "hpipm/blas/print.h"
#include "hpipm/blas/struct.h"

// print a matrix structure
void print_mat(int m, int n, struct mat* sA, int ai, int aj) {
    int ii, jj;
    for (ii = 0; ii < m; ii++) {
        for (jj = 0; jj < n; jj++) {
            printf("%9.5f ", MATEL(sA, ai + ii, aj + jj));
        }
        printf("\n");
    }
    printf("\n");
}


// print the transposed of a matrix structure
void print_tran_mat(int m, int n, struct mat* sA, int ai, int aj) {
    int ii, jj;
    for (jj = 0; jj < n; jj++) {
        for (ii = 0; ii < m; ii++) {
            printf("%9.5f ", MATEL(sA, ai + ii, aj + jj));
        }
        printf("\n");
    }
    printf("\n");
}


// print a vector structure
void print_vec(int m, struct vec* sA, int ai) {
    int ii;
    for (ii = 0; ii < m; ii++) {
        printf("%9.5f\n", VECEL(sA, ai + ii));
    }
    printf("\n");
}


// print the transposed of a vector structure
void print_tran_vec(int m, struct vec* sA, int ai) {
    int ii;
    for (ii = 0; ii < m; ii++) {
        printf("%9.5f ", VECEL(sA, ai + ii));
    }
    printf("\n\n");
}


// print a matrix structure
void print_to_file_mat(FILE* file, int m, int n, struct mat* sA, int ai, int aj) {
    int ii, jj;
    for (ii = 0; ii < m; ii++) {
        for (jj = 0; jj < n; jj++) {
            fprintf(file, "%9.5f ", MATEL(sA, ai + ii, aj + jj));
        }
        fprintf(file, "\n");
    }
    fprintf(file, "\n");
}


// print a matrix structure
void print_to_file_exp_mat(FILE* file, int m, int n, struct mat* sA, int ai, int aj) {
    for (int ii = 0; ii < m; ii++) {
        for (int jj = 0; jj < n; jj++) {
            fprintf(file, "%9.5e ", MATEL(sA, ai + ii, aj + jj));
        }
        fprintf(file, "\n");
    }
    fprintf(file, "\n");
}


// print a matrix structure
void print_to_string_mat(char** buf_out, int m, int n, struct mat* sA, int ai, int aj) {
    int ii, jj;
    for (ii = 0; ii < m; ii++) {
        for (jj = 0; jj < n; jj++) {
            *buf_out += sprintf(*buf_out, "%9.5f ", MATEL(sA, ai + ii, aj + jj));
        }
        *buf_out += sprintf(*buf_out, "\n");
    }
    *buf_out += sprintf(*buf_out, "\n");
}


// print a vector structure
void print_to_file_vec(FILE* file, int m, struct vec* sA, int ai) {
    int ii;
    for (ii = 0; ii < m; ii++) {
        fprintf(file, "%9.5f\n", VECEL(sA, ai + ii));
    }
    fprintf(file, "\n");
}


// print a vector structure
void print_to_string_vec(char** buf_out, int m, struct vec* sA, int ai) {
    int ii;
    for (ii = 0; ii < m; ii++) {
        *buf_out += sprintf(*buf_out, "%9.5f\n", VECEL(sA, ai + ii));
    }
    *buf_out += sprintf(*buf_out, "\n");
}


// print the transposed of a vector structure
void print_to_file_tran_vec(FILE* file, int m, struct vec* sA, int ai) {
    int ii;
    for (ii = 0; ii < m; ii++) {
        fprintf(file, "%9.5f ", VECEL(sA, ai + ii));
    }
    fprintf(file, "\n\n");
}


// print the transposed of a vector structure
void print_to_string_tran_vec(char** buf_out, int m, struct vec* sA, int ai) {
    int ii;
    for (ii = 0; ii < m; ii++) {
        *buf_out += sprintf(*buf_out, "%9.5f ", VECEL(sA, ai + ii));
    }
    *buf_out += sprintf(*buf_out, "\n\n");
}


// print a matrix structure
void print_exp_mat(int m, int n, struct mat* sA, int ai, int aj) {
    int ii, jj;
    for (ii = 0; ii < m; ii++) {
        for (jj = 0; jj < n; jj++) {
            printf("%9.5e ", MATEL(sA, ai + ii, aj + jj));
        }
        printf("\n");
    }
    printf("\n");
}


// print the transposed of a matrix structure
void print_exp_tran_mat(int m, int n, struct mat* sA, int ai, int aj) {
    int ii, jj;
    for (jj = 0; jj < n; jj++) {
        for (ii = 0; ii < m; ii++) {
            printf("%9.5e ", MATEL(sA, ai + ii, aj + jj));
        }
        printf("\n");
    }
    printf("\n");
}


// print a vector structure
void print_exp_vec(int m, struct vec* sA, int ai) {
    int ii;
    for (ii = 0; ii < m; ii++) {
        printf("%9.5e\n", VECEL(sA, ai + ii));
    }
    printf("\n");
}


// print the transposed of a vector structure
void print_exp_tran_vec(int m, struct vec* sA, int ai) {
    int ii;
    for (ii = 0; ii < m; ii++) {
        printf("%9.5e ", VECEL(sA, ai + ii));
    }
    printf("\n\n");
}

void int_print_mat(int row, int col, int* A, int lda) {
    int i, j;
    for (i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            printf("%d ", A[i + lda * j]);
        }
        printf("\n");
    }
    printf("\n");
}