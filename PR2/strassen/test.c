#include "strassen.h"
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <getopt.h>
#include <limits.h>

int main(int argc, char* argv[]) {
    // parameters
    char* fname = "mat.txt";
    int vb = -1;
    int min = -1;
    int max = 2;

    int *a, *b, *c;
    double tm;
    struct timeval t;

    for (size_t sz = 100; sz <= 1000; sz += 100) {
        size_t best_thresh = -1;
        double best_tm = INT_MAX;

        for (size_t thresh = 0; thresh <= sz; thresh += 10) {
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
            tm = t.tv_sec + t.tv_usec * 0.000001;
            if (tm < best_tm) {
                best_tm = tm;
                best_thresh = thresh;
            }

            // free matrices
            free(a);
            free(b);
            free(c);
        }
        printf("sz=%zu, thresh=%zu, time=%f\n", sz, best_thresh, best_tm);
    }
}
