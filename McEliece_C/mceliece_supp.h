//
// Created by linux on 09.09.24.
//
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include "../galois/galois.h"

#ifndef MCELIECE_SUPP_H
#define MCELIECE_SUPP_H

int galois_pow(const int a, const int b, const int w);

void swap_bytes(unsigned char* a, unsigned char* b, const int bit_a, const int bit_b);

void XOR_bits_in_place(unsigned char* a, const uint16_t b, const int bit_a, const int bit_b);

// polynomial manipulation methods

int galois_eval_poly(const int a, const int d, const uint16_t poly[d+1], const int w);

bool has_zeroes(const int d, const uint16_t poly[d+1], const int w);

void add_poly(const int d1, const int d2, const uint16_t poly1[d1+1], const uint16_t poly2[d2+1], uint16_t* result);

void mult_by_var_poly(const int d, const int n, const uint16_t poly1[d+1], uint16_t* result);

void mult_poly(const int d1, const int d2, const uint16_t poly1[d1+1], const uint16_t poly2[d2+1], uint16_t* result, const int w);

bool fully_zero(const int d, const uint16_t poly[d+1]);

int poly_degree(const int d, const uint16_t poly[d+1]);

void reduce_mod_poly(const int d1, int d2, uint16_t poly1[d1+1], const uint16_t poly2[d2+1], uint16_t* quotient, const int w);

void EEA_standard(const int d1, const int d2, const uint16_t poly1[d1+1], const uint16_t poly2[d2+1], const int w, uint16_t* factor1, uint16_t* factor2, uint16_t* gcd);

void EEA_patterson(const int d1, const int d2, const uint16_t poly1[d1+1], const uint16_t goppa[d2+1], uint16_t* A, uint16_t* B, const int w);

bool invert_poly(const int d, uint16_t poly[d+1], const int dm, const uint16_t modulus[dm+1], uint16_t result[dm], const int w);

void poly_pow_two_mod(const int d1, const int d2, uint16_t poly1[d1+1], const uint16_t poly2[d2+1], const int exponent, uint16_t* result, const int w);

void sqrt_poly(const int d1, const int d2, uint16_t poly1[d1+1], const uint16_t poly2[d2+1], uint16_t* result, const int w);

// misc

void random_perm(const int n, uint16_t* perm);

void generate_error(const int n, const int w, unsigned char error[n]);

// matrix methods

void multiply_matrix(const int d1, const int d2, const int d3, uint16_t m1[d1][d2], uint16_t m2[d2][d3], uint16_t result[d1][d3], const int m);

void swap_col(const int i, const int j, const int d1, const int d2, uint16_t H[d1][d2]);

void swap_row(const int i, const int r1, const int j, const int r2, int d1, const int d2, uint16_t H[d1][d2]);

void add_row(const int i, const int r1, const int j, const int r2, int d1, const int d2, uint16_t H[d1][d2]);

void row_reduce(const int d1, const int d2, uint16_t H[d1][d2], const int m) ;

bool fully_reduced_parity(const int d1, const int d2, const uint16_t H[d1][d2], const int m, int* failed_col);


#endif //MCELIECE_SUPP_H
