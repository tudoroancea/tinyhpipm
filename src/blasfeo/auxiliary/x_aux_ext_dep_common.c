


/* creates a zero matrix */
void ZEROS(REAL** pA, int row, int col) {
    blasfeo_malloc((void**) pA, (row * col) * sizeof(REAL));
    REAL* A = *pA;
    int i;
    for (i = 0; i < row * col; i++) A[i] = 0.0;
}


/* creates a zero matrix aligned to a cache line */
void ZEROS_ALIGN(REAL** pA, int row, int col) {
    blasfeo_malloc_align((void**) pA, (row * col) * sizeof(REAL));
    REAL* A = *pA;
    int i;
    for (i = 0; i < row * col; i++) A[i] = 0.0;
}


/* frees matrix */
void FREE(REAL* pA) {
    blasfeo_free(pA);
}


/* frees aligned matrix */
void FREE_ALIGN(REAL* pA) {
    blasfeo_free_align(pA);
}


/* prints a matrix in column-major format */
void PRINT_MAT(int m, int n, REAL* A, int lda) {
    int i, j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            printf("%9.5f ", A[i + lda * j]);
        }
        printf("\n");
    }
    printf("\n");
    return;
}


/* prints the transposed of a matrix in column-major format */
void PRINT_TRAN_MAT(int row, int col, REAL* A, int lda) {
    int i, j;
    for (j = 0; j < col; j++) {
        for (i = 0; i < row; i++) {
            printf("%9.5f ", A[i + lda * j]);
        }
        printf("\n");
    }
    printf("\n");
}


/* prints a matrix in column-major format */
void PRINT_TO_FILE_MAT(FILE* file, int row, int col, REAL* A, int lda) {
    int i, j;
    for (i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            fprintf(file, "%9.5f ", A[i + lda * j]);
        }
        fprintf(file, "\n");
    }
    fprintf(file, "\n");
}

/* prints a matrix in column-major format */
void PRINT_TO_FILE_EXP_MAT(FILE* file, int row, int col, REAL* A, int lda) {
    int i, j;
    for (i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            fprintf(file, "%9.5e ", A[i + lda * j]);
        }
        fprintf(file, "\n");
    }
    fprintf(file, "\n");
}


/* prints a matrix in column-major format */
void PRINT_TO_STRING_MAT(char** buf_out, int row, int col, REAL* A, int lda) {
    int i, j;
    for (i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            *buf_out += sprintf(*buf_out, "%9.5f ", A[i + lda * j]);
        }
        *buf_out += sprintf(*buf_out, "\n");
    }
    *buf_out += sprintf(*buf_out, "\n");
    return;
}


/* prints the transposed of a matrix in column-major format */
void PRINT_TO_FILE_TRAN_MAT(FILE* file, int row, int col, REAL* A, int lda) {
    int i, j;
    for (j = 0; j < col; j++) {
        for (i = 0; i < row; i++) {
            fprintf(file, "%9.5f ", A[i + lda * j]);
        }
        fprintf(file, "\n");
    }
    fprintf(file, "\n");
}


/* prints a matrix in column-major format (exponential notation) */
void PRINT_EXP_MAT(int m, int n, REAL* A, int lda) {
    int i, j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            printf("%e\t", A[i + lda * j]);
        }
        printf("\n");
    }
    printf("\n");
}


/* prints the transposed of a matrix in column-major format (exponential notation) */
void PRINT_EXP_TRAN_MAT(int row, int col, REAL* A, int lda) {
    int i, j;
    for (j = 0; j < col; j++) {
        for (i = 0; i < row; i++) {
            printf("%e\t", A[i + lda * j]);
        }
        printf("\n");
    }
    printf("\n");
}
