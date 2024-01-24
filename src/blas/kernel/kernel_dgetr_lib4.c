// transposed of general matrices, read along panels, write across panels
void kernel_dgetr_4_lib4(int tri, int kmax, int kna, double alpha, double* A, double* C, int sdc) {

    if (tri == 1) {
        // A is lower triangular, C is upper triangular
        // kmax+1 4-wide + end 3x3 triangle

        kmax += 1;
    }

    const int bs = 4;

    int k;

    k = 0;

    if (kmax < kna)
        goto cleanup_loop;

    if (kna > 0) {
        for (; k < kna; k++) {
            C[0 + bs * 0] = alpha * A[0 + bs * 0];
            C[0 + bs * 1] = alpha * A[1 + bs * 0];
            C[0 + bs * 2] = alpha * A[2 + bs * 0];
            C[0 + bs * 3] = alpha * A[3 + bs * 0];

            C += 1;
            A += bs;
        }
        C += bs * (sdc - 1);
    }

    for (; k < kmax - 3; k += 4) {
        C[0 + bs * 0] = alpha * A[0 + bs * 0];
        C[0 + bs * 1] = alpha * A[1 + bs * 0];
        C[0 + bs * 2] = alpha * A[2 + bs * 0];
        C[0 + bs * 3] = alpha * A[3 + bs * 0];

        C[1 + bs * 0] = alpha * A[0 + bs * 1];
        C[1 + bs * 1] = alpha * A[1 + bs * 1];
        C[1 + bs * 2] = alpha * A[2 + bs * 1];
        C[1 + bs * 3] = alpha * A[3 + bs * 1];

        C[2 + bs * 0] = alpha * A[0 + bs * 2];
        C[2 + bs * 1] = alpha * A[1 + bs * 2];
        C[2 + bs * 2] = alpha * A[2 + bs * 2];
        C[2 + bs * 3] = alpha * A[3 + bs * 2];

        C[3 + bs * 0] = alpha * A[0 + bs * 3];
        C[3 + bs * 1] = alpha * A[1 + bs * 3];
        C[3 + bs * 2] = alpha * A[2 + bs * 3];
        C[3 + bs * 3] = alpha * A[3 + bs * 3];

        C += bs * sdc;
        A += bs * bs;
    }

cleanup_loop:

    for (; k < kmax; k++) {
        C[0 + bs * 0] = alpha * A[0 + bs * 0];
        C[0 + bs * 1] = alpha * A[1 + bs * 0];
        C[0 + bs * 2] = alpha * A[2 + bs * 0];
        C[0 + bs * 3] = alpha * A[3 + bs * 0];

        C += 1;
        A += bs;
    }

    if (tri == 1) {
        // end 3x3 triangle
        kna = (bs - (bs - kna + kmax) % bs) % bs;

        if (kna == 1) {
            C[0 + bs * 1] = alpha * A[1 + bs * 0];
            C[0 + bs * 2] = alpha * A[2 + bs * 0];
            C[0 + bs * 3] = alpha * A[3 + bs * 0];
            C[1 + bs * (sdc + 1)] = alpha * A[2 + bs * 1];
            C[1 + bs * (sdc + 2)] = alpha * A[3 + bs * 1];
            C[2 + bs * (sdc + 2)] = alpha * A[3 + bs * 2];
        } else if (kna == 2) {
            C[0 + bs * 1] = alpha * A[1 + bs * 0];
            C[0 + bs * 2] = alpha * A[2 + bs * 0];
            C[0 + bs * 3] = alpha * A[3 + bs * 0];
            C[1 + bs * 2] = alpha * A[2 + bs * 1];
            C[1 + bs * 3] = alpha * A[3 + bs * 1];
            C[2 + bs * (sdc + 2)] = alpha * A[3 + bs * 2];
        } else {
            C[0 + bs * 1] = alpha * A[1 + bs * 0];
            C[0 + bs * 2] = alpha * A[2 + bs * 0];
            C[0 + bs * 3] = alpha * A[3 + bs * 0];
            C[1 + bs * 2] = alpha * A[2 + bs * 1];
            C[1 + bs * 3] = alpha * A[3 + bs * 1];
            C[2 + bs * 3] = alpha * A[3 + bs * 2];
        }
    }
}


// transposed of general matrices, read along panels, write across panels
void kernel_dgetr_3_lib4(int tri, int kmax, int kna, double alpha, double* A, double* C, int sdc) {

    if (tri == 1) {
        // A is lower triangular, C is upper triangular
        // kmax+1 3-wide + end 2x2 triangle

        kmax += 1;
    }

    const int bs = 4;

    int k;

    k = 0;

    if (kmax < kna)
        goto cleanup_loop;

    if (kna > 0) {
        for (; k < kna; k++) {
            C[0 + bs * 0] = alpha * A[0 + bs * 0];
            C[0 + bs * 1] = alpha * A[1 + bs * 0];
            C[0 + bs * 2] = alpha * A[2 + bs * 0];

            C += 1;
            A += bs;
        }
        C += bs * (sdc - 1);
    }

    for (; k < kmax - 3; k += 4) {
        C[0 + bs * 0] = alpha * A[0 + bs * 0];
        C[0 + bs * 1] = alpha * A[1 + bs * 0];
        C[0 + bs * 2] = alpha * A[2 + bs * 0];

        C[1 + bs * 0] = alpha * A[0 + bs * 1];
        C[1 + bs * 1] = alpha * A[1 + bs * 1];
        C[1 + bs * 2] = alpha * A[2 + bs * 1];

        C[2 + bs * 0] = alpha * A[0 + bs * 2];
        C[2 + bs * 1] = alpha * A[1 + bs * 2];
        C[2 + bs * 2] = alpha * A[2 + bs * 2];

        C[3 + bs * 0] = alpha * A[0 + bs * 3];
        C[3 + bs * 1] = alpha * A[1 + bs * 3];
        C[3 + bs * 2] = alpha * A[2 + bs * 3];

        C += bs * sdc;
        A += bs * bs;
    }

cleanup_loop:

    for (; k < kmax; k++) {
        C[0 + bs * 0] = alpha * A[0 + bs * 0];
        C[0 + bs * 1] = alpha * A[1 + bs * 0];
        C[0 + bs * 2] = alpha * A[2 + bs * 0];

        C += 1;
        A += bs;
    }

    if (tri == 1) {
        // end 2x2 triangle
        kna = (bs - (bs - kna + kmax) % bs) % bs;

        if (kna == 1) {
            C[0 + bs * 1] = alpha * A[1 + bs * 0];
            C[0 + bs * 2] = alpha * A[2 + bs * 0];
            C[1 + bs * (sdc + 1)] = alpha * A[2 + bs * 1];
        } else {
            C[0 + bs * 1] = alpha * A[1 + bs * 0];
            C[0 + bs * 2] = alpha * A[2 + bs * 0];
            C[1 + bs * 2] = alpha * A[2 + bs * 1];
        }
    }
}


// transposed of general matrices, read along panels, write across panels
void kernel_dgetr_2_lib4(int tri, int kmax, int kna, double alpha, double* A, double* C, int sdc) {

    if (tri == 1) {
        // A is lower triangular, C is upper triangular
        // kmax+1 2-wide + end 1x1 triangle

        kmax += 1;
    }

    const int bs = 4;

    int k;

    k = 0;

    if (kmax < kna)
        goto cleanup_loop;

    if (kna > 0) {
        for (; k < kna; k++) {
            C[0 + bs * 0] = alpha * A[0 + bs * 0];
            C[0 + bs * 1] = alpha * A[1 + bs * 0];

            C += 1;
            A += bs;
        }
        C += bs * (sdc - 1);
    }

    for (; k < kmax - 3; k += 4) {
        C[0 + bs * 0] = alpha * A[0 + bs * 0];
        C[0 + bs * 1] = alpha * A[1 + bs * 0];

        C[1 + bs * 0] = alpha * A[0 + bs * 1];
        C[1 + bs * 1] = alpha * A[1 + bs * 1];

        C[2 + bs * 0] = alpha * A[0 + bs * 2];
        C[2 + bs * 1] = alpha * A[1 + bs * 2];

        C[3 + bs * 0] = alpha * A[0 + bs * 3];
        C[3 + bs * 1] = alpha * A[1 + bs * 3];

        C += bs * sdc;
        A += bs * bs;
    }

cleanup_loop:

    for (; k < kmax; k++) {
        C[0 + bs * 0] = alpha * A[0 + bs * 0];
        C[0 + bs * 1] = alpha * A[1 + bs * 0];

        C += 1;
        A += bs;
    }

    if (tri == 1) {
        // end 1x1 triangle
        C[0 + bs * 1] = alpha * A[1 + bs * 0];
    }
}


// transposed of general matrices, read along panels, write across panels
void kernel_dgetr_1_lib4(int tri, int kmax, int kna, double alpha, double* A, double* C, int sdc) {

    if (tri == 1) {
        // A is lower triangular, C is upper triangular
        // kmax+1 1-wide

        kmax += 1;
    }

    const int bs = 4;

    int k;

    k = 0;

    if (kmax < kna)
        goto cleanup_loop;

    if (kna > 0) {
        for (; k < kna; k++) {
            C[0 + bs * 0] = alpha * A[0 + bs * 0];

            C += 1;
            A += bs;
        }
        C += bs * (sdc - 1);
    }

    for (; k < kmax - 3; k += 4) {
        C[0 + bs * 0] = alpha * A[0 + bs * 0];

        C[1 + bs * 0] = alpha * A[0 + bs * 1];

        C[2 + bs * 0] = alpha * A[0 + bs * 2];

        C[3 + bs * 0] = alpha * A[0 + bs * 3];

        C += bs * sdc;
        A += bs * bs;
    }

cleanup_loop:

    for (; k < kmax; k++) {
        C[0 + bs * 0] = alpha * A[0 + bs * 0];

        C += 1;
        A += bs;
    }
}


// transposed of general matrices, read across panels, write along panels
void kernel_dgetr_4_0_lib4(int kmax, double* A, int sda, double* B) {
    const int ps = 4;
    int k;
    for (k = 0; k < kmax - 3; k += 4) {
        //
        B[0 + ps * 0] = A[0 + ps * 0];
        B[0 + ps * 1] = A[1 + ps * 0];
        B[0 + ps * 2] = A[2 + ps * 0];
        B[0 + ps * 3] = A[3 + ps * 0];
        //
        B[1 + ps * 0] = A[0 + ps * 1];
        B[1 + ps * 1] = A[1 + ps * 1];
        B[1 + ps * 2] = A[2 + ps * 1];
        B[1 + ps * 3] = A[3 + ps * 1];
        //
        B[2 + ps * 0] = A[0 + ps * 2];
        B[2 + ps * 1] = A[1 + ps * 2];
        B[2 + ps * 2] = A[2 + ps * 2];
        B[2 + ps * 3] = A[3 + ps * 2];
        //
        B[3 + ps * 0] = A[0 + ps * 3];
        B[3 + ps * 1] = A[1 + ps * 3];
        B[3 + ps * 2] = A[2 + ps * 3];
        B[3 + ps * 3] = A[3 + ps * 3];

        A += ps * sda;
        B += ps * ps;
    }
    for (; k < kmax; k++) {
        //
        B[0 + ps * 0] = A[0 + ps * 0];
        B[1 + ps * 0] = A[0 + ps * 1];
        B[2 + ps * 0] = A[0 + ps * 2];
        B[3 + ps * 0] = A[0 + ps * 3];

        A += 1;
        B += ps;
    }
}
