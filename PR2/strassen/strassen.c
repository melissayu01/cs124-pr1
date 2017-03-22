
// MATRIX HELPERS

// me(m, stride, i, j)
//    Return a pointer to matrix element `m[i][j]` -- the element
//    at row `i` and column `j`. The matrix has dimensions
//    `n_row` and `n_col`. Requires: `i < n_row && j < n_col`
static inline int* me(int* m, size_t stride, size_t i, size_t j) {
    assert(j < stride);
    return &m[i * stride + j];
}

// pad(m, sz, i, j)
//    Return a pointer to a new matrix statically padded to the
//    smallest size `padded_sz` such that `padded_sz >= thresh * 1 >> k`
//    for any non-negative integer `k`. Padding is applied on
//    bottom and right columns.
int* pad(int* m, size_t sz, size_t thresh) {
    // find size of statically padded matrix
    size_t padded_sz = thresh;
    while (padded_sz < sz)
        padded_sz << 1;

    // create new matrix with size padded_sz and fill it
    padded_m = calloc(sizeof(int) * padded_sz * padded_sz);
    for (size_t i = 0; i < padded_sz; ++i)
        for (size_t j = 0; j < sz; ++j)
            *me(padded_m, padded_sz, i, j) = *me(m, sz, i, j);
}

// add_submat(c, a, b, sz, stride, i, j, sub_flag)
//    adds submatrices `a` and `b`, each of dimension `sz x sz`
//    and places the result in the submatrix of `c` with the
//    same dimensions, where `a`, `b`, `c` are `sz x sz` matrices.
void add_submat(int* c, int* a, int* b, size_t sz, size_t stride, int sub_flag) {
    for (size_t i = 0; i < sz; ++i)
        for (size_t j = 0; j < sz; ++j)
            *me(c, stride, i, j) = *me(a, stride, i, j) +
                math.pow(-1, sub_flag) * *me(b, stride, i, j);
}

// add(c, a, b, sz, i, j, sub_flag)
//    adds matrices `a` and `b`, each of dimension `sz x sz`
//    and places the result in the matrix `c` with the
//    same dimensions.
void add(int* c, int* a, int* b, size_t sz, int sub_flag) {
    for (size_t i = 0; i < sz; ++i)
        for (size_t j = 0; j < sz; ++j)
            *me(c, sz, i, j) = *me(a, sz, i, j) + math.pow(-1, sub_flag) * *me(b, sz, i, j);
}


// MATRIX MULTIPLICATION

// base_matrix_multiply(c, sz, a, b)
//    `a`, `b`, and `c` are square matrices with dimension `sz`.
//    Computes the matrix product `a x b` and stores it in `c`.
void base_matrix_multiply(int* c, size_t sz, size_t stride, int* a, int* b) {
    // clear submatrix of `c`
    for (size_t i = 0; i < sz; ++i)
        for (size_t j = 0; j < sz; ++j)
            *me(c, stride, i, j) = 0;

    // compute product and update `c`
    for (size_t i = 0; i < sz; ++i)
        for (size_t k = 0; k < sz; ++k)
            for (size_t j = 0; j < sz; ++j)
                *me(c, stride, i, j) += *me(a, stride, i, k) * *me(b, stride, k, j);
}

// rec_strassen_mult(c, a, b, s_sz, sz, thresh)
//    `a` and `b` are square submatrices with dimension `s_sz`, taken from
//    larger matrices with dimension `sz`. Recursively computes the matrix
//    product `a x b` using Strassen's algorithm and stores this result to
//    the submatrix of `c` with the same dimensions. Switches to using n^3 base
//    algorithm when `sz < thresh`.
double* rec_strassen_mult(int* c, int* a, int* b, size_t sz, size_t stride, size_t thresh) {
    // allocate scratch matrices
    int* t = malloc(sizeof(int) * sz * sz / 2);

    if (s_sz <= thresh) // base algorithm
        base_matrix_multiply(c, sz, stride, a, b);
    else {
        // get submatrices
        size_t subsz = sz/2;
        int* a11 = me(a, stride, 0, 0);
        int* a12 = me(a, stride, 0, subsz);
        int* a21 = me(a, stride, subsz, 0);
        int* a22 = me(a, stride, subsz, subsz);

        int* b11 = me(b, stride, 0, 0);
        int* b12 = me(b, stride, 0, subsz);
        int* b21 = me(b, stride, subsz, 0);
        int* b22 = me(b, stride, subsz, subsz);

        int* c11 = me(c, stride, 0, 0);
        int* c12 = me(c, stride, 0, subsz);
        int* c21 = me(c, stride, subsz, 0);
        int* c22 = me(c, stride, subsz, subsz);

        // compute result
        add_submat(c12, a21, a11, sz, sz, 1);
        add_submat(c21, b11, b12, s_sz, sz, 0);
        rec_strassen_mult(c22, c12, c21, s_sz, sz, thresh); // M6

        add_submat(c12, a12, a22, s_sz, sz, 1);
        add_submat(c21, b21, b22, s_sz, sz, 0);
        rec_strassen_mult(c11, c12, c21, s_sz, sz, thresh); // M7

        add_submat(c12, a11, a22, s_sz, sz, 0);
        add_submat(c21, b11, b22, s_sz, sz, 0);
        rec_strassen_mult(t1, c12, c21, s_sz, sz, thresh); // M1



    }

    // free scratch matrices
    free(t);
}

double* strassen(double* a, double* b, size_t sz, size_t cutoff) {
    // pa = pad(a, sz, cutoff)
    // pb = pad(b, sz, cutoff)
    // pc = rec_mult(double* a, double* b, size_t sz, size_t cutoff)

    // c = malloc()
    // for (size_t i = 0; i < sz; ++i)
    // for (size_t j = 0; j < sz; ++j)
    // *me(c, padded_sz, i, j) = *me(pc, sz, i, j)
}
