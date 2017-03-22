#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>
#include <getopt.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>

// MATRIX HELPER FUNCTIONS

// me(m, stride, i, j)
//    Return a pointer to submatrix element `m[i][j]` -- the element
//    at row `i` and column `j`, where the stride (number of columns)
//    in the original matrix is `stride`. Requires: `i < sz && j < sz`
static inline int* me(int* m, size_t stride, size_t i, size_t j) {
    return &m[i * stride + j];
}

// print_mat(m, sz, stride)
//    Prints square matrix `m` with dimension `sz`.
void print_mat(int* m, size_t sz, size_t stride) {
    for (size_t i = 0; i < sz; ++i) {
        for (size_t j = 0; j < sz; ++j)
            printf("\t%d", *me(m, stride, i, j));
        printf("\n");
    }
    printf("\n");
}

// equals(a, b, sz)
//    Returns 1 if matrices `a` and `b` with size `sz` are equal element-wise.
int equals(int* a, int* b, size_t sz) {
    for (size_t i = 0; i < sz; ++i)
        for (size_t j = 0; j < sz; ++j)
            if (*me(a, sz, i, j) != *me(b, sz, i, j)) {
                printf("%d != %d\n", *me(a, sz, i, j), *me(b, sz, i, j));
                return 0;
            }
    return 1;
}

// zero_mat(m, sz, stride)
//    Initializes all entries in submatrix at `m` with `sz` and `stride` to 0.
void zero_mat(int* m, size_t sz, size_t stride) {
    for (size_t i = 0; i < sz; ++i)
        for (size_t j = 0; j < sz; ++j)
            *me(m, stride, i, j) = 0;
}

// pad(pm, m, sz, thresh)
//    Allocates memory for and statically pads matrix `m` to the
//    smallest size `padded_sz` such that `padded_sz >= thresh * 1 >> k`
//    for any non-negative integer `k`. Padding is applied on the
//    bottom and right columns. Returns padded size.
size_t pad(int** pm, int* m, size_t sz, size_t thresh) {
    // find size of statically padded matrix
    size_t padded_sz = thresh ? thresh : sz;
    while (padded_sz < sz)
        padded_sz = padded_sz << 1;

    // create new matrix with size padded_sz and fill it
    *pm = calloc(padded_sz * padded_sz, sizeof(int));
    for (size_t i = 0; i < sz; ++i)
        for (size_t j = 0; j < sz; ++j)
            *me(*pm, padded_sz, i, j) = *me(m, sz, i, j);

    return padded_sz;
}

// add_submat(c, a, b, sz, stride, sub_flag)
//    Adds submatrices `a` and `b`, each of dimension `sz x sz`
//    and places the result in the submatrix of `c` with the
//    same dimensions, where the stride for all submatrices is `stride`.
//    Subtracts `a` and `b` if `sub_flag == 1`.
void add_submat(int* c, int* a, int* b, size_t sz, size_t stride, int sub_flag) {
    for (size_t i = 0; i < sz; ++i)
        for (size_t j = 0; j < sz; ++j)
            *me(c, stride, i, j) = *me(a, stride, i, j) +
                (sub_flag ? -1 : 1) * *me(b, stride, i, j);
}


// MATRIX MULTIPLICATION

// base_matrix_multiply(c, a, b, sz, stride)
//    `a`, `b`, and `c` are square submatrices with dimension `sz` and
//    stride `stride`. Computes the matrix product `a x b` and stores it in `c`.
void base_matrix_multiply(int* c, int* a, int* b, size_t sz, size_t stride) {
    // clear submatrix of `c`
    zero_mat(c, sz, stride);

    // compute product and update `c`
    for (size_t i = 0; i < sz; ++i)
        for (size_t k = 0; k < sz; ++k)
            for (size_t j = 0; j < sz; ++j)
                *me(c, stride, i, j) += *me(a, stride, i, k) * *me(b, stride, k, j);
}

// rec_strassen_mult(c, a, b, sz, stride, thresh)
//    `a` and `b` are square submatrices with dimension `sz`, taken from
//    larger matrices with dimension `stride`. Recursively computes the matrix
//    product `a x b` using Strassen's algorithm and stores this result to
//    the submatrix of `c` with the same dimensions. Switches to using n^3 base
//    algorithm when `sz < thresh`.
void rec_strassen_mult(int* c, int* a, int* b, size_t sz, size_t stride, size_t thresh) {
    size_t subsz = (size_t) sz / 2;

    if (sz <= thresh) // base algorithm
        base_matrix_multiply(c, a, b, sz, stride);
    else if (sz == 1) {
        c[0] = a[0] * b[0];
    } else {
        // allocate scratch matrices
        int* t = (int*) calloc(subsz * stride, sizeof(int));

        // get submatrices
        int* t1 = me(t, stride, 0, 0);
        int* t2 = me(t, stride, 0, subsz);

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
        add_submat(c12, a21, a11, subsz, stride, 1);
        add_submat(c21, b11, b12, subsz, stride, 0);
        rec_strassen_mult(c22, c12, c21, subsz, stride, thresh); // M6

        add_submat(c12, a12, a22, subsz, stride, 1);
        add_submat(c21, b21, b22, subsz, stride, 0);
        rec_strassen_mult(c11, c12, c21, subsz, stride, thresh); // M7

        add_submat(c12, a11, a22, subsz, stride, 0);
        add_submat(c21, b11, b22, subsz, stride, 0);
        rec_strassen_mult(t1, c12, c21, subsz, stride, thresh); // M1

        add_submat(c11, t1, c11, subsz, stride, 0);
        add_submat(c22, t1, c22, subsz, stride, 0);
        add_submat(t2, a21, a22, subsz, stride, 0);
        rec_strassen_mult(c21, t2, b11, subsz, stride, thresh); // M2

        add_submat(c22, c22, c21, subsz, stride, 1);
        add_submat(t1, b21, b11, subsz, stride, 1);
        rec_strassen_mult(t2, a22, t1, subsz, stride, thresh); // M4

        add_submat(c21, c21, t2, subsz, stride, 0);
        add_submat(c11, c11, t2, subsz, stride, 0);
        add_submat(t1, b12, b22, subsz, stride, 1);
        rec_strassen_mult(c12, a11, t1, subsz, stride, thresh); // M3

        add_submat(c22, c22, c12, subsz, stride, 0);
        add_submat(t2, a11, a12, subsz, stride, 0);
        rec_strassen_mult(t1, t2, b22, subsz, stride, thresh); // M5

        add_submat(c12, c12, t1, subsz, stride, 0);
        add_submat(c11, c11, t1, subsz, stride, 1);

        // free scratch matrices
        free(t);
    }
}

// strassen_matrix_multiply(c, a, b, sz, thresh)
//    `a` and `b` are square matrices with dimension `sz`. Computes the matrix
//    product `a x b` using Strassen's algorithm and stores this result to
//    the matrix `c` with the same dimensions. Switches to using n^3 base
//    algorithm when `sz < thresh`.
void strassen_matrix_multiply(int* c, int* a, int* b, size_t sz, size_t thresh) {
    int *pa, *pb, *pc;
    size_t psz;

    // pad input
    psz = pad(&pa, a, sz, thresh);
    psz = pad(&pb, b, sz, thresh);

    // allocate padded output matrix
    pc = (int*) malloc(sizeof(int) * psz * psz);

    // clear output matrices
    zero_mat(pc, psz, psz);
    zero_mat(c, sz, sz);

    // call recursive multiply
    rec_strassen_mult(pc, pa, pb, psz, psz, thresh);

    // remove padding and store to c
    for (size_t i = 0; i < sz; ++i)
        for (size_t j = 0; j < sz; ++j)
            *me(c, sz, i, j) = *me(pc, psz, i, j);

    // free padded matrices
    free(pa);
    free(pb);
    free(pc);
}


int main(int argc, char* argv[]) {
    srand(time(NULL));

    int max_int = 2;
    size_t sz = 1000;
    size_t thresh = 100;

    // setup
    assert(sz > 0);
    assert(sz < (size_t) sqrt(SIZE_MAX / sizeof(int)));
    struct timeval s_time0, s_time1, b_time0, b_time1;

    // allocate matrices
    int* a = (int*) malloc(sizeof(int) * sz * sz);
    int* b = (int*) malloc(sizeof(int) * sz * sz);
    int* c_base = (int*) malloc(sizeof(int) * sz * sz);
    int* c_strassen = (int*) malloc(sizeof(int) * sz * sz);

    // fill in source matrices
    for (size_t i = 0; i < sz; ++i)
	   for (size_t j = 0; j < sz; ++j)
            *me(a, sz, i, j) = rand() % max_int;

    for (size_t i = 0; i < sz; ++i)
	   for (size_t j = 0; j < sz; ++j)
	      *me(b, sz, i, j) = rand() % max_int;

    // strassen matrix multiply
    printf("computing strassen...\n");
    gettimeofday(&s_time0, NULL);
    strassen_matrix_multiply(c_strassen, a, b, sz, thresh);
    gettimeofday(&s_time1, NULL);

    // base matrix multiply
    printf("computing base...\n");
    gettimeofday(&b_time0, NULL);
    base_matrix_multiply(c_base, a, b, sz, sz);
    gettimeofday(&b_time1, NULL);

    // print results
    int is_equal = equals(c_base, c_strassen, sz);
    printf("matrices are equal: %d\n", is_equal);

    // compute times, print times and ratio
    timersub(&s_time1, &s_time0, &s_time1);
    timersub(&b_time1, &b_time0, &b_time1);
    printf("base multiply time %ld.%06lds\n", b_time1.tv_sec, b_time1.tv_usec);
    double time_ratio = (s_time1.tv_sec + s_time1.tv_usec * 0.000001)
        / (b_time1.tv_sec + b_time1.tv_usec * 0.000001);
    printf("strassen multiply time %ld.%06lds (%gx)\n",
           s_time1.tv_sec, s_time1.tv_usec, time_ratio);

    free(a);
    free(b);
    free(c_base);
    free(c_strassen);
}
