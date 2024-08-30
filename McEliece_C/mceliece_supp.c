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

void generate_error(int n, int w, int error[n]) {
    // create n-bit error vector with Hamming weight w
    // first, create one with w 1s at the beginning
    memset(error, 0, sizeof(int) * n);
    for (int i = 0; i < w; i++) {
        error[i] = 1;
    }
    // now shuffle error using the Fisher-Yates shuffle
    for (int i = n-1; i > 0; i--) {
        int j = randombytes_uniform(i+1);
        int temp = error[i];
        error[i] = error[j];
        error[j] = temp;
    }
}

void multiply_matrix(const int d1, const int d2, const int d3, int m1[d1][d2], int m2[d2][d3], int result[d1][d3], const int m) {
    // m1 has dimensions d1 x d2, and m2 has dimensions d2 x d3
    for (int i = 0; i < d1; i++) {
        for (int j = 0; j < d3; j++) {
            result[i][j] = 0;
            for (int k = 0; k < d2; k++) {
                result[i][j] ^= galois_single_multiply(m1[i][k], m2[k][j],m);
            }
        }
    }
}

void swap_col(const int i, const int j, const int d1, const int d2, int H[d1][d2]) {
    for (int row = 0; row < d1; row++) {
        const int temp = H[row][i];
        H[row][i] = H[row][j];
        H[row][j] = temp;
    }

}

int* add_poly(const int d1, const int d2, const int poly1[d1], const int poly2[d2]) {
    // adds two polys with coefficients in F_2^m (actually for any m, not just the fixed one here)
    int max;
    if (d1 > d2) {
        max = d1;
    }
    else {
        max = d2;
    }
    int* result = malloc(sizeof(int) * max);
    for (int i = 0; i < max; i++) {
        result[i] = poly1[i]^poly2[i];
    }
    return result;
}

void reduce_mod_poly(const int d1, const int d2, const int poly2[d1+1], const int poly2[d2+1] int result[d2+1]) {
    // reduces poly1 mod poly2 within F_2^m
}

void invert_poly(const int degree, const int poly[degree+1], int result[degree+1], const int w) {

}

// the following are specific algorithms for row/col operations on matrices containing elements in F_2^m interpreted as m-bit vectors

void swap_row(int i, int r1, int j, int r2, int d1, int d2, int H[d1][d2]) {
    // note we are swapping row m*i + r1 and row m*j + r2
    for (int k = 0; k < d2; k++) {
        if ((H[i][k] >> r1) & 1 && ~(H[j][k] >> r2) & 1)  {
            H[i][k] &= ~(1 << r1);
            H[j][k] |= 1 << r2;
        }
        else if (~(H[i][k] >> r1) & 1 && (H[j][k] >> r2) & 1) {
            H[j][k] &= ~(1 << r2);
            H[i][k] |= 1 << r1;
        }
    }
}

void add_row(int i, int r1, int j, int r2, int d1, int d2, int H[d1][d2]) {
    // note we are adding row m*i + r1 to row m*j + r2
    if (r2 > r1) {
        for (int col = 0; col < d2; col++) {
            H[j][col] ^= ((H[i][col] & (1 << r1)) << (r2-r1));
        }
    }
    else {
        for (int col = 0; col < d2; col++) {
            H[j][col] ^= ((H[i][col] & (1 << r1)) >> (r1-r2));
        }
    }
}

void row_reduce(const int d1, const int d2, int H[d1][d2], const int m) {
    // this row reduces a matrix without doing any column operations, might have some cols of zeros afterwards
    int count = 0;
    bool increase = false;
    int col = 0;
    while (count < m*d1) {
        if (increase) {
            count++;
            increase = false;
        }
        for (int k = count; k < m*d1; k++) {
            int row = k / m;
            int r = k % m;
                if ((H[row][col] >> r & 1) == 1) {
                    increase = true;
                    // move up row with a leading one as high as possible
                    swap_row(row, r, count / m, count % m, d1, d2, H);
                    // add row in question to all rows below if they have a 1 in the leading col
                    for (int j = count+1; j < m*d1; j++) {
                        const int entry = j / m;
                        const int bit = j % m;
                        if (((H[entry][col] >> bit) & 1) == 1) {
                            add_row(count / m, count % m, entry, bit, d1, d2, H);
                        }
                    }
                    // add row in question to all rows above if they have a 1 in the leading col
                    for (int j = 0; j < count; j++) {
                        const int entry = j / m;
                        const int bit = j % m;
                        if (((H[entry][col] >> bit) & 1) == 1) {
                            add_row(count / m, count % m, entry, bit, d1, d2, H);
                        }
                    }
                }
        }
        col++;
    }

}