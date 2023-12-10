// create a matrix structure for a matrix of size m*n by dynamically allocating the memory
void ALLOCATE_MAT(int m, int n, struct MAT* sA) {
    size_t size = MEMSIZE_MAT(m, n);
    void* mem;
    blasfeo_malloc_align(&mem, size);
    CREATE_MAT(m, n, sA, mem);
    return;
}


// free memory of a matrix structure
void FREE_MAT(struct MAT* sA) {
    blasfeo_free_align(sA->mem);
    return;
}


// create a vector structure for a vector of size m by dynamically allocating the memory
void ALLOCATE_VEC(int m, struct VEC* sa) {
    size_t size = MEMSIZE_VEC(m);
    void* mem;
    blasfeo_malloc_align(&mem, size);
    CREATE_VEC(m, sa, mem);
    return;
}


// free memory of a matrix structure
void FREE_VEC(struct VEC* sa) {
    blasfeo_free_align(sa->mem);
    return;
}


// print a matrix structure
void PRINT_MAT(int m, int n, struct MAT* sA, int ai, int aj) {
    int ii, jj;
    for (ii = 0; ii < m; ii++) {
        for (jj = 0; jj < n; jj++) {
            printf("%9.5f ", MATEL(sA, ai + ii, aj + jj));
        }
        printf("\n");
    }
    printf("\n");
    return;
}


// print the transposed of a matrix structure
void PRINT_TRAN_MAT(int m, int n, struct MAT* sA, int ai, int aj) {
    int ii, jj;
    for (jj = 0; jj < n; jj++) {
        for (ii = 0; ii < m; ii++) {
            printf("%9.5f ", MATEL(sA, ai + ii, aj + jj));
        }
        printf("\n");
    }
    printf("\n");
    return;
}


// print a vector structure
void PRINT_VEC(int m, struct VEC* sa, int ai) {
    int ii;
    for (ii = 0; ii < m; ii++) {
        printf("%9.5f\n", VECEL(sa, ai + ii));
    }
    printf("\n");
    return;
}


// print the transposed of a vector structure
void PRINT_TRAN_VEC(int m, struct VEC* sa, int ai) {
    int ii;
    for (ii = 0; ii < m; ii++) {
        printf("%9.5f ", VECEL(sa, ai + ii));
    }
    printf("\n\n");
    return;
}


// print a matrix structure
void PRINT_TO_FILE_MAT(FILE* file, int m, int n, struct MAT* sA, int ai, int aj) {
    int ii, jj;
    for (ii = 0; ii < m; ii++) {
        for (jj = 0; jj < n; jj++) {
            fprintf(file, "%9.5f ", MATEL(sA, ai + ii, aj + jj));
        }
        fprintf(file, "\n");
    }
    fprintf(file, "\n");
    return;
}


// print a matrix structure
void PRINT_TO_FILE_EXP_MAT(FILE* file, int m, int n, struct MAT* sA, int ai, int aj) {
    int ii, jj;
    for (ii = 0; ii < m; ii++) {
        for (jj = 0; jj < n; jj++) {
            fprintf(file, "%9.5e ", MATEL(sA, ai + ii, aj + jj));
        }
        fprintf(file, "\n");
    }
    fprintf(file, "\n");
    return;
}


// print a matrix structure
void PRINT_TO_STRING_MAT(char** buf_out, int m, int n, struct MAT* sA, int ai, int aj) {
    int ii, jj;
    for (ii = 0; ii < m; ii++) {
        for (jj = 0; jj < n; jj++) {
            *buf_out += sprintf(*buf_out, "%9.5f ", MATEL(sA, ai + ii, aj + jj));
        }
        *buf_out += sprintf(*buf_out, "\n");
    }
    *buf_out += sprintf(*buf_out, "\n");
    return;
}


// print a vector structure
void PRINT_TO_FILE_VEC(FILE* file, int m, struct VEC* sa, int ai) {
    int ii;
    for (ii = 0; ii < m; ii++) {
        fprintf(file, "%9.5f\n", VECEL(sa, ai + ii));
    }
    fprintf(file, "\n");
    return;
}


// print a vector structure
void PRINT_TO_STRING_VEC(char** buf_out, int m, struct VEC* sa, int ai) {
    int ii;
    for (ii = 0; ii < m; ii++) {
        *buf_out += sprintf(*buf_out, "%9.5f\n", VECEL(sa, ai + ii));
    }
    *buf_out += sprintf(*buf_out, "\n");
    return;
}


// print the transposed of a vector structure
void PRINT_TO_FILE_TRAN_VEC(FILE* file, int m, struct VEC* sa, int ai) {
    int ii;
    for (ii = 0; ii < m; ii++) {
        fprintf(file, "%9.5f ", VECEL(sa, ai + ii));
    }
    fprintf(file, "\n\n");
    return;
}


// print the transposed of a vector structure
void PRINT_TO_STRING_TRAN_VEC(char** buf_out, int m, struct VEC* sa, int ai) {
    int ii;
    for (ii = 0; ii < m; ii++) {
        *buf_out += sprintf(*buf_out, "%9.5f ", VECEL(sa, ai + ii));
    }
    *buf_out += sprintf(*buf_out, "\n\n");
    return;
}


// print a matrix structure
void PRINT_EXP_MAT(int m, int n, struct MAT* sA, int ai, int aj) {
    int ii, jj;
    for (ii = 0; ii < m; ii++) {
        for (jj = 0; jj < n; jj++) {
            printf("%9.5e ", MATEL(sA, ai + ii, aj + jj));
        }
        printf("\n");
    }
    printf("\n");
    return;
}


// print the transposed of a matrix structure
void PRINT_EXP_TRAN_MAT(int m, int n, struct MAT* sA, int ai, int aj) {
    int ii, jj;
    for (jj = 0; jj < n; jj++) {
        for (ii = 0; ii < m; ii++) {
            printf("%9.5e ", MATEL(sA, ai + ii, aj + jj));
        }
        printf("\n");
    }
    printf("\n");
    return;
}


// print a vector structure
void PRINT_EXP_VEC(int m, struct VEC* sa, int ai) {
    int ii;
    for (ii = 0; ii < m; ii++) {
        printf("%9.5e\n", VECEL(sa, ai + ii));
    }
    printf("\n");
    return;
}


// print the transposed of a vector structure
void PRINT_EXP_TRAN_VEC(int m, struct VEC* sa, int ai) {
    int ii;
    for (ii = 0; ii < m; ii++) {
        printf("%9.5e ", VECEL(sa, ai + ii));
    }
    printf("\n\n");
    return;
}
