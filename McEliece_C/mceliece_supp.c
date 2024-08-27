#include <math.h>
#include <stdlib.h>

#include "galois.h"
#include "sodium.h"

int galois_pow(int a, int b, int w) {
    int result = 1;
    while (b > 0) {
        result = galois_single_multiply(result, a, w);
        b--;
    }
    return result;
}

int galois_eval_poly(int a, int degree, int poly[], int w) {
    // degree is degree of poly: due to type decay this information is lost after passing
    int result = 0;
    for (int bit = 0; bit <= degree; bit++) {
            int temp = poly[bit];
            temp = galois_single_multiply(temp, galois_pow(a, bit, w), w);
            result = result ^ temp;
        }
    return result;
}

bool has_zeroes(int poly[], int degree, int w) {
    for (int i = 0; i < (1 << w); i++) {
        if (galois_eval_poly(i, degree, poly, w) == 0) {
            return true;
        }
    }
    return false;
}

int* random_perm(int n) {
    int* perm = malloc(sizeof(int) * n);
    for (int i = 0; i < n; i++) {
        perm[i] = i;
    }
    for (int i = n-1; i >= 0; i--) {
        int j = randombytes_uniform(i+1);
        int temp = perm[i];
        perm[i] = perm[j];
        perm[j] = temp;
    }
    return perm;
}

void multiply_matrix(int m1[][], int m2[][], int d1, int d2, int d3, int result[d1][d3]) {
    // m1 has dimensions d1 x d2, and m2 has dimensions d2 x d3
    for (int i = 0; i < d1; i++) {
        for (int j = 0; j < d2; j++) {
            result[i][j] = 0;
            for (int k = 0; k < d3; k++) {
                result[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }
}

// the following are specific algorithms for row/col operations on matrices containing elements in F_2^m interpreted as m-bit vectors

void swap_row(int i, int r1, int j, int r2, int d1, int d2, int H[d1][d2]) {
    // note we are swapping row m*i + r1 and row m*j + r2
    for (int k = 0; k < d2; k++) {
        if ((H[i][k] << r1) & 1 && ~(H[j][k] >> r2) & 1)  {
            H[i][k] &= ~(1 << r1);
            H[j][k] |= 1 << r2;
        }
        else if (~(H[i][k] << r1) & 1 && (H[j][k] >> r2) & 1) {
            H[j][k] &= ~(1 << r2);
            H[i][k] |= 1 << r1;
        }
    }
}

void add_row(int i, int r1, int j, int r2, int d1, int d2, int H[d1][d2]) {
    // note we are adding row m*i + r1 to row m*j + r2
    for (int k = 0; k < d2; k++) {
        H[j][k] ^= ((H[i][k] & ~(1 << r1)) << (r2-r1));
    }
}