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

// RANDOM

// xrandom()
//    Return a pseudo-random number in the range [0, XRAND_MAX].
//    We use our own generator to ensure values computed on different
//    OSes will follow the same sequence.
#define XRAND_MAX 10 // 0x7FFFFFFF
static uint64_t xrandom_seed;
unsigned xrandom() {
    xrandom_seed = xrandom_seed * 6364136223846793005U + 1U;
    return (xrandom_seed >> 32) & XRAND_MAX;
}


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
//    Allocates memory for and fills matrix statically padded to the
//    smallest size `padded_sz` such that `padded_sz >= thresh * 1 >> k`
//    for any non-negative integer `k`. Padding is applied on
//    bottom and right columns. Returns padded size.
size_t pad(int** pm, int* m, size_t sz, size_t thresh) {
    // find size of statically padded matrix
    size_t padded_sz = thresh;
    while (padded_sz < sz)
        padded_sz = padded_sz << 1;

    // create new matrix with size padded_sz and fill it
    *pm = calloc(padded_sz * padded_sz, sizeof(int));
    for (size_t i = 0; i < padded_sz; ++i)
        for (size_t j = 0; j < sz; ++j)
            *me(*pm, padded_sz, i, j) = *me(m, sz, i, j);

    return padded_sz;
}

// add_submat(c, a, b, sz, stride, i, j, sub_flag)
//    adds submatrices `a` and `b`, each of dimension `sz x sz`
//    and places the result in the submatrix of `c` with the
//    same dimensions, where `a`, `b`, `c` are `sz x sz` matrices.
void add_submat(int* c, int* a, int* b, size_t sz, size_t stride, int sub_flag) {
    for (size_t i = 0; i < sz; ++i)
        for (size_t j = 0; j < sz; ++j)
            *me(c, stride, i, j) = *me(a, stride, i, j) +
                pow(-1, sub_flag) * *me(b, stride, i, j);
}

// add(c, a, b, sz, i, j, sub_flag)
//    adds matrices `a` and `b`, each of dimension `sz x sz`
//    and places the result in the matrix `c` with the
//    same dimensions.
void add(int* c, int* a, int* b, size_t sz, int sub_flag) {
    for (size_t i = 0; i < sz; ++i)
        for (size_t j = 0; j < sz; ++j)
            *me(c, sz, i, j) = *me(a, sz, i, j) + pow(-1, sub_flag) * *me(b, sz, i, j);
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
void rec_strassen_mult(int* c, int* a, int* b, size_t sz, size_t stride, size_t thresh) {
    size_t subsz = sz / 2;

    // allocate scratch matrices
    int* t = (int*) malloc(sizeof(int) * subsz * sz);

    if (sz <= thresh) // base algorithm
        base_matrix_multiply(c, sz, stride, a, b);
    else {
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
        add_submat(t2, a21, a11, subsz, stride, 0);
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
    }

    // free scratch matrices
    free(t);
}

void strassen_matrix_multiply(int* c, int* a, int* b, size_t sz, size_t thresh) {
    int *pa, *pb, *pc;

    // pad input
    size_t psz = pad(&pa, a, sz, thresh);
    size_t pszb = pad(&pb, b, sz, thresh);
    assert(psz == pszb);

    // allocate output matrices
    pc = (int*) malloc(sizeof(int) * psz * psz);
    rec_strassen_mult(pc, pa, pb, psz, psz, thresh);

    // remove padding
    for (size_t i = 0; i < sz; ++i)
        for (size_t j = 0; j < sz; ++j)
            *me(c, sz, i, j) = *me(pc, psz, i, j);

    // free padded matrices
    free(pa);
    free(pb);
    free(pc);
}

void print_mat(int* m, size_t sz) {
    for (size_t i = 0; i < sz; ++i) {
        for (size_t j = 0; j < sz; ++j)
            printf(" %d", *me(m, sz, i, j));
        printf("\n");
    }
}

int main(int argc, char* argv[]) {
    size_t sz = 2;
    size_t thresh = 1;

    // allocate matrices
    int* a = (int*) malloc(sizeof(int) * sz * sz);
    int* b = (int*) malloc(sizeof(int) * sz * sz);
    int* c = (int*) malloc(sizeof(int) * sz * sz);

    // fill in source matrices
    for (size_t i = 0; i < sz; ++i)
	   for (size_t j = 0; j < sz; ++j)
            *me(a, sz, i, j) = xrandom();

    for (size_t i = 0; i < sz; ++i)
	   for (size_t j = 0; j < sz; ++j)
	      *me(b, sz, i, j) = xrandom();

    print_mat(a, sz);
    print_mat(b, sz);

    // compute `c = a x b`
    printf("computing base...\n");
    base_matrix_multiply(c, sz, sz, a, b);
    print_mat(c, sz);

    // clear c
    for (size_t i = 0; i < sz; ++i)
        for (size_t j = 0; j < sz; ++j)
            *me(c, sz, i, j) = 0;

    printf("computing strassen...\n");
    strassen_matrix_multiply(c, a, b, sz, thresh);
    print_mat(c, sz);

    free(a);
    free(b);
    free(c);
}
