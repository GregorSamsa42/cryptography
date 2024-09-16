// supplementary functions for the McEliece cryptosystem

#include <stdlib.h>
#include "sodium.h"
#include <string.h>
#include <stdbool.h>
#include "../galois/galois.h"

#ifndef MCELIECE_SUPP_H
#define MCELIECE_SUPP_H

int galois_pow(const int a, const int b, const int w);

// polynomial manipulation methods

int galois_eval_poly(const int a, const int d, const int poly[d+1], const int w);

bool has_zeroes(const int d, const int poly[d+1], const int w);

void add_poly(const int d1, const int d2, const int poly1[d1+1], const int poly2[d2+1], int* result);

void mult_by_var_poly(const int d, const int n, const int poly1[d+1], int* result);

void mult_poly(const int d1, const int d2, const int poly1[d1+1], const int poly2[d2+1], int* result, const int w);

bool fully_zero(const int d, const int poly[d+1]);

int poly_degree(const int d, const int poly[d+1]);

void reduce_mod_poly(const int d1, int d2, int poly1[d1+1], const int poly2[d2+1], int* quotient, const int w);

void EEA_standard(const int d1, const int d2, const int poly1[d1+1], const int poly2[d2+1], const int w, int* factor1, int* factor2, int* gcd);

void EEA_patterson(const int d1, const int d2, const int poly1[d1+1], const int goppa[d2+1], int* A, int* B, const int w);

bool invert_poly(const int d, int poly[d+1], const int dm, const int modulus[dm+1], int result[dm], const int w);

void poly_pow_two_mod(const int d1, const int d2, int poly1[d1+1], const int poly2[d2+1], const int exponent, int* result, const int w);

void sqrt_poly(const int d1, const int d2, int poly1[d1+1], const int poly2[d2+1], int* result, const int w);

// misc

void random_perm(const int n, int* perm);

void generate_error(const int n, const int w, int error[n]);

// matrix methods

void multiply_matrix(const int d1, const int d2, const int d3, int m1[d1][d2], int m2[d2][d3], int result[d1][d3], const int m);

void swap_col(const int i, const int j, const int d1, const int d2, int H[d1][d2]);

void swap_row(const int i, const int r1, const int j, const int r2, int d1, const int d2, int H[d1][d2]);

void add_row(const int i, const int r1, const int j, const int r2, int d1, const int d2, int H[d1][d2]);

void row_reduce(const int d1, const int d2, int H[d1][d2], const int m) ;

bool fully_reduced_parity(const int d1, const int d2, const int H[d1][d2], const int m, int* failed_col);


#endif //MCELIECE_SUPP_H
