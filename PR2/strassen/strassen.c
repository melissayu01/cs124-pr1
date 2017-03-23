#include "strassen.h"
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <getopt.h>
#include <limits.h>

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
        add_submat(t1, a11, a21, subsz, stride, 1);
        add_submat(t2, b22, b12, subsz, stride, 1);
        rec_strassen_mult(c21, t1, t2, subsz, stride, thresh); // M6

        add_submat(t1, a21, a22, subsz, stride, 0);
        add_submat(t2, b12, b11, subsz, stride, 1);
        rec_strassen_mult(c22, t1, t2, subsz, stride, thresh); // M7

        add_submat(t1, t1, a11, subsz, stride, 1);
        add_submat(t2, b22, t2, subsz, stride, 1);
        rec_strassen_mult(c11, t1, t2, subsz, stride, thresh); // M1

        add_submat(t1, a12, t1, subsz, stride, 1);
        rec_strassen_mult(c12, t1, b22, subsz, stride, thresh); // M2

        add_submat(c12, c22, c12, subsz, stride, 0);
        rec_strassen_mult(t1, a11, b11, subsz, stride, thresh); // M4

        add_submat(c11, c11, t1, subsz, stride, 0);
        add_submat(c12, c11, c12, subsz, stride, 0);
        add_submat(c11, c11, c21, subsz, stride, 0);
        add_submat(t2, t2, b21, subsz, stride, 1);
        rec_strassen_mult(c21, a22, t2, subsz, stride, thresh); // M3

        add_submat(c21, c11, c21, subsz, stride, 1);
        add_submat(c22, c11, c22, subsz, stride, 0);
        rec_strassen_mult(c11, a12, b21, subsz, stride, thresh); // M5

        add_submat(c11, t1, c11, subsz, stride, 0);

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
    pc = (int*) calloc(psz * psz, sizeof(int));

    // call recursive multiply
    rec_strassen_mult(pc, pa, pb, psz, psz, thresh);

    // clear output matrix
    zero_mat(c, sz, sz);

    // remove padding and store to c
    for (size_t i = 0; i < sz; ++i)
        for (size_t j = 0; j < sz; ++j)
            *me(c, sz, i, j) = *me(pc, psz, i, j);

    // free padded matrices
    free(pa);
    free(pb);
    free(pc);
}


// I/O AND TESTING

// file_to_mat(fname, a, b, sz)
//    Reads 2 square matrices, each with dimension `sz`, sequentially into
//    the memory at `a` and `b`. (Does not allocate memory for `a` and `b`.)
void file_to_mat(char* fname, int* a, int* b, size_t sz) {
    FILE* fp = fopen(fname, "r");
    if (!fp) {
        perror(fname);
        exit(EXIT_FAILURE);
    }
    for (size_t i = 0; i < pow(sz, 2); ++i)
        fscanf(fp, "%d", &a[i]);
    for (size_t i = 0; i < pow(sz, 2); ++i)
        fscanf(fp, "%d", &b[i]);

    fclose(fp);
}

// rand_mat_to_file(fname, sz, max_val)
//    Writes 2 square matrices, each with dimension `sz`, sequentially into
//    file `fname`. All entries in matrices are random ints in the range [min_val, max_val).
void rand_mat_to_file(char* fname, size_t sz, int min_val, int max_val) {
    srand(time(NULL));
    FILE* fp = fopen(fname, "w");
    if (!fp) {
        perror(fname);
        exit(EXIT_FAILURE);
    }

    int range = max_val - min_val;
    for (size_t i = 0; i < 2 * pow(sz, 2); ++i)
        fprintf(fp, "%d\n", rand() % range + min_val);
    fclose(fp);
}

// usage()
//    Prints usage instructions to stderr.
static inline void usage(void) {
    fprintf(stderr, "Usage: ./strassen <verbosity> <dimension> <inputfile>\n");
}

typedef struct matrix_statistics {
    int corner[4];
    int diagonal_sum;
    int* diagonal;
} matrix_statistics;

// compute_statistics(m, sz)
//    Compute and return some statistics about matrix `m`.
matrix_statistics compute_statistics(int* m, size_t sz) {
    matrix_statistics mstat;
    mstat.corner[0] = *me(m, sz, 0, 0);
    mstat.corner[1] = *me(m, sz, 0, sz-1);
    mstat.corner[2] = *me(m, sz, sz-1, 0);
    mstat.corner[3] = *me(m, sz, sz-1, sz-1);
    mstat.diagonal = calloc(sz, sizeof(int));
    mstat.diagonal_sum = 0;
    for (size_t i = 0; i < sz; ++i) {
        mstat.diagonal[i] = *me(m, sz, i, i);
        mstat.diagonal_sum += mstat.diagonal[i];
    }
    return mstat;
}

// run_strassen(vb, c, a, b, sz, thresh)
//    Run strassen's algorithm.
struct timeval run_strassen(int vb, int* c, int* a, int* b, size_t sz, size_t thresh) {
    struct timeval time0, time1, time2, time3;

    // compute `c = a x b`
    if (vb > 0)
        printf("computing strassen...\n");
    gettimeofday(&time0, NULL);
    strassen_matrix_multiply(c, a, b, sz, thresh);
    gettimeofday(&time1, NULL);
    timersub(&time1, &time0, &time1);
    matrix_statistics strassen_mstat = compute_statistics(c, sz);

    // compute times, print times and ratio
    if (vb > 0) {
        zero_mat(c, sz, sz);

        printf("computing base...\n");
        gettimeofday(&time2, NULL);
        base_matrix_multiply(c, a, b, sz, sz);
        gettimeofday(&time3, NULL);
        timersub(&time3, &time2, &time3);
        matrix_statistics base_mstat = compute_statistics(c, sz);

        printf("base multiply time %ld.%06ds\n", time3.tv_sec, time3.tv_usec);
        double time_ratio = (time1.tv_sec + time1.tv_usec * 0.000001)
            / (time3.tv_sec + time3.tv_usec * 0.000001);
        printf("strassen multiply time %ld.%06ds (%gx)\n",
               time1.tv_sec, time1.tv_usec, time_ratio);

        // print statistics and differences
        for (int i = 0; i < 4; ++i)
            printf("corner statistic %d: base %d\n"
                  "                    strassen %d (%d%% off)\n",
                  i, base_mstat.corner[i], strassen_mstat.corner[i],
                  100 * abs(strassen_mstat.corner[i] - base_mstat.corner[i]) / base_mstat.corner[i]);
        printf("diagonal sum statistic: base %d\n"
          "                        strassen %d (%d%% off)\n",
          base_mstat.diagonal_sum, strassen_mstat.diagonal_sum,
          100 * abs(strassen_mstat.diagonal_sum - base_mstat.diagonal_sum) / base_mstat.diagonal_sum);

        free(base_mstat.diagonal);
    }

    // print diagonal
    if (vb == 0)
        for (size_t i = 0; i < sz; ++i)
            printf("%d\n", strassen_mstat.diagonal[i]);
    free(strassen_mstat.diagonal);

    fflush(stdout);
    return time1;
}

// compute_optimal_thresh()
//    Compute optimal thresholds for varying matrix sizes.
void compute_optimal_thresh(void) {
    // parameters
    char* fname = "mat.txt";
    int vb = -1;
    int min = -1;
    int max = 2;
    int num_trials = 5;

    int *a, *b, *c;
    double tm;
    struct timeval t;

    for (size_t sz = 550; sz <= 1024; sz += 75) {
        size_t best_thresh = -1;
        double best_tm = INT_MAX;

        for (size_t thresh = 0; thresh <= sz; thresh += 1) {
            tm = 0;
            for (int i = 0; i < num_trials; ++i) {
                // generate random matrix for testing
                rand_mat_to_file(fname, sz, min, max);

                // allocate matrices
                a = (int*) calloc(sz * sz, sizeof(int));
                b = (int*) calloc(sz * sz, sizeof(int));
                c = (int*) calloc(sz * sz, sizeof(int));

                // fill in source matrices
                file_to_mat(fname, a, b, sz);

                // run strassen's algorithm
                t = run_strassen(vb, c, a, b, sz, thresh);
                tm += t.tv_sec + t.tv_usec * 0.000001;

                // free matrices
                free(a);
                free(b);
                free(c);
            }
            if (tm / num_trials < best_tm) {
                best_tm = tm / num_trials;
                best_thresh = thresh;
            }
        }
        printf("sz=%zu, thresh=%zu, time=%f\n", sz, best_thresh, best_tm);
    }
}

int main(int argc, char* argv[]) {
    compute_optimal_thresh();
    // // parameters
    // size_t thresh = 125;
    //
    // // read options
    // int vb = strtoul(argv[1], NULL, 0);
    // size_t sz = strtoul(argv[2], NULL, 0);
    // char* fname = argv[3];
    //
    // if (!(argc == 4 && sz > 0 && sz < (size_t) sqrt(SIZE_MAX / sizeof(int))
    //     && (vb == 0 || vb == 1) && fname != NULL)) {
    //     usage();
    //     exit(EXIT_FAILURE);
    // }
    //
    // // TODO: REMOVE BEFORE SUBMITTING -- generate random matrix for testing
    // rand_mat_to_file(fname, sz, -1, 2);
    //
    // // allocate matrices
    // int* a = (int*) calloc(sz * sz, sizeof(int));
    // int* b = (int*) calloc(sz * sz, sizeof(int));
    // int* c = (int*) calloc(sz * sz, sizeof(int));
    //
    // // fill in source matrices
    // file_to_mat(fname, a, b, sz);
    //
    // // run strassen's algorithm
    // run_strassen(vb, c, a, b, sz, thresh);
    //
    // free(a);
    // free(b);
    // free(c);
}
