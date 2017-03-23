#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>

// MATRIX HELPER FUNCTIONS

// me(m, stride, i, j)
//    Return a pointer to submatrix element `m[i][j]` -- the element
//    at row `i` and column `j`, where the stride (number of columns)
//    in the original matrix is `stride`. Requires: `i < sz && j < sz`
static inline int* me(int* m, size_t stride, size_t i, size_t j) {
    return &m[i * stride + j];
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
