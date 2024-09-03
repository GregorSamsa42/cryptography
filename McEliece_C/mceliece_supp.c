
#include <stdlib.h>
#include "sodium.h"
#include <stdbool.h>

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

void random_perm(int n, int* perm) {
    // perm needs to have allocated n ints of memory
    for (int i = 0; i < n; i++) {
        perm[i] = i;
    }
    for (int i = n-1; i >= 0; i--) {
        int j = randombytes_uniform(i+1);
        int temp = perm[i];
        perm[i] = perm[j];
        perm[j] = temp;
    }
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

// polynomial manipulation methods: note we always give the degree of a polynomial: the methods then ignore any coefficients above

void add_poly(const int d1, const int d2, const int poly1[d1+1], const int poly2[d2+1], int* result) {
    // adds two polys with coefficients in F_2^m (actually for any m, not just the fixed one here)
    // result should have as many integers allocated as the maximum of d1+1 and d2+1
    int max;
    int min;
    bool first_larger;
    if (d1 > d2) {
        max = d1;
        min = d2;
        first_larger = true;
    }
    else {
        max = d2;
        min = d1;
        first_larger = false;
    }
    for (int i = 0; i <= min; i++) {
        result[i] = poly1[i]^poly2[i];
    }
    if (first_larger) {
        for (int i = min+1; i <= max; i++) {
            result[i] = poly1[i];
        }
    }
    else {
        for (int i = min+1; i <= max; i++) {
            result[i] = poly2[i];
        }
    }
}

void mult_by_var_poly(const int d, const int n, const int poly1[d+1], int* result) {
    // multiplies p(x) by x^^n
    // result must have enough space allocated already (but this is guaranteed by mult_poly)
    memset(result, 0, sizeof(int) * (d+n+1));
    for (int i = n; i <= d+n; i++) {
        result[i] = poly1[i-n];
    }
}

void mult_poly(const int d1, const int d2, const int poly1[d1+1], const int poly2[d2+1], int* result, int w) {
    // result needs d1+d2 * sizeof(int) allocated space
    memset(result, 0, sizeof(int) * (d1+d2+1));
    for (int j = 0; j <= d2; j++) {
        int* temp = malloc(sizeof(int) * (d1+d2+1));
        mult_by_var_poly(d1, j, poly1, temp);
        for (int k = 0; k <= d1+d2; k++) {
            temp[k] = galois_single_multiply(temp[k], poly2[j], w);
        }
        add_poly(d1+d2, d1+d2, result, temp, result);
        free(temp);
    }
}

int poly_degree(const int d, const int poly[d+1]) {
    // outputs the actual degree of a polynomial with d coefficients
    // if it is zero, output -1
    for (int i = d; i >= 0; i--) {
        if (poly[i] != 0) {
            return i;
        }
    }
    return -1;
}

void reduce_mod_poly(const int d1, int d2, int poly1[d1+1], const int poly2[d2+1], int* quotient, const int w) {

    // reduces poly1 mod poly2 within F_2^w in place and returns quotient
    // poly2 must be non-zero
    // quotient must have d1-d2+1 allocated integers
    d2 = poly_degree(d2, poly2);
    if (d2 == -1) {
        quotient = NULL;
    }
    memset(quotient, 0, sizeof(int) * (d1-d2+1));
    for (int i = d1; i >= d2; i--) {
        // create temp poly that is poly2 scaled appropriately, and subtract it from poly1
        int* temp = malloc(sizeof(int) * (i+1));
        mult_by_var_poly(d1, i-d2, poly2, temp);
        // note poly2[d2] is guaranteed to be non-zero here
        const int factor = galois_single_divide(poly1[i],poly2[d2],w);
        for (int k = 0; k <= i; k++) {
            temp[k] = galois_single_multiply(temp[k], factor,w);
        }
        // add the appropriate poly to the quotient
        // temp_quotient will have degree at most i-d2 but give it more to prevent garbage values
        int temp_quotient[d1-d2+1];
        memset(temp_quotient, 0, sizeof(int) * (d1-d2+1));
        temp_quotient[i-d2] = factor;
        add_poly(i-d2, d1-d2, temp_quotient, quotient, quotient);
        // subtract temp off of poly1
        add_poly(d1, d1, poly1, temp, poly1);
        free(temp);
    }
}

void EEA_standard(int d1, int d2, const int poly1[d1+1], const int poly2[d2+1], const int w, int* factor1, int* factor2, int* gcd) {
    // initialise variables in Euclid's algorithm. Want r - q*g = rem = u*poly1 + v*poly2.
    int* r = malloc(sizeof(int) * (d1+1));
    for (int i = 0; i <= d1; i++) {
        r[i] = poly1[i];
    }
    int* g = malloc(sizeof(int) * (d2+1));
    for (int i = 0; i <= d2; i++) {
        g[i] = poly2[i];
    }
    // note that for the Bezout coefficients we know they have degrees bounded by max(d1,d2)/2
    int* u = malloc(sizeof(int) * (d1+d2+1));
    memset(u, 0, sizeof(int) * (d1+d2+1));
    int* u1 = malloc(sizeof(int) * (d1+d2+1));
    memset(u1, 0, sizeof(int) * (d1+d2+1));
    u1[0] = 1;
    int* v = malloc(sizeof(int) * (d1+d2+1));
    memset(v, 0, sizeof(int) * (d1+d2+1));
    v[0] = 1;
    int* v1 = malloc(sizeof(int) * (d1+d2+1));
    memset(v1, 0, sizeof(int) * (d1+d2+1));
    int* u2 = malloc(sizeof(int) * (d1+d2+1));
    memset(u2, 0, sizeof(int) * (d1+d2+1));
    int* v2 = malloc(sizeof(int) * (d1+d2+1));
    memset(v2, 0, sizeof(int) * (d1+d2+1));
    int remdeg = 1;
    // terminate when the remainder is zero (i.e. degree -1), then can read off Bezout coefficients from u and v
    while (remdeg > -1) {
        int* q = malloc(sizeof(int) * (d1+d2+1));
        reduce_mod_poly(d1, d2, r, g, q, w);
        int* temp = r;
        r = g;
        g = temp;
        d1 = poly_degree(d1, r);
        d2 = poly_degree(d2, g);
        remdeg = d2;
        // copy over arrays, don't just change pointers
        for (int i = 0; i <= d1+d2; i++) {
            u2[i] = u1[i];
        }
        for (int i = 0; i <= d1+d2; i++) {
            v2[i] = v1[i];
        }
        for (int i = 0; i <= d1+d2; i++) {
            u1[i] = u[i];
        }
        for (int i = 0; i <= d1+d2; i++) {
            v1[i] = v[i];
        }
        int* q_times_v = malloc(sizeof(int) * (d1+d2+1));
        mult_poly(d1, d2, q, v, q_times_v,w);
        int* q_times_u = malloc(sizeof(int) * (d1+d2+1));
        mult_poly(d1, d2, q, u, q_times_u,w);
        add_poly(d1+d2,d1+d2,u2,q_times_u,u);
        add_poly(d1+d2,d1+d2,v2,q_times_v,v);
        free(q_times_v);
        free(q_times_u);
        free(q);
    }
}
    void EEA_patterson(int d1, int d2, const int poly1[d1+1], const int goppa[d2+1], const int w, int* A, int* B) {
    // initialise variables in Euclid's algorithm. Want r - q*g = rem = u*poly1 + v*poly2.
    // same as before, except the break condition and the outputs are different
    const int d2_fixed = d2;
    int* r = malloc(sizeof(int) * (d1+1));
    for (int i = 0; i <= d1; i++) {
        r[i] = poly1[i];
    }
    int* g = malloc(sizeof(int) * (d2+1));
    for (int i = 0; i <= d2; i++) {
        g[i] = goppa[i];
    }
    // note that for the Bezout coefficients we know they have degrees bounded by max(d1,d2)/2
    int* u = malloc(sizeof(int) * (d1+d2+1));
    memset(u, 0, sizeof(int) * (d1+d2+1));
    int* u1 = malloc(sizeof(int) * (d1+d2+1));
    memset(u1, 0, sizeof(int) * (d1+d2+1));
    u1[0] = 1;
    int* v = malloc(sizeof(int) * (d1+d2+1));
    memset(v, 0, sizeof(int) * (d1+d2+1));
    v[0] = 1;
    int* v1 = malloc(sizeof(int) * (d1+d2+1));
    memset(v1, 0, sizeof(int) * (d1+d2+1));
    int* u2 = malloc(sizeof(int) * (d1+d2+1));
    memset(u2, 0, sizeof(int) * (d1+d2+1));
    int* v2 = malloc(sizeof(int) * (d1+d2+1));
    memset(v2, 0, sizeof(int) * (d1+d2+1));
    int remdeg = 1;
    // terminate when the remainder is zero (i.e. degree -1), then can read off Bezout coefficients from u and v
    while (remdeg > d2_fixed/2 || poly_degree(d1+d2,u1) > (d2_fixed-1)/2) {
        int* q = malloc(sizeof(int) * (d1+d2+1));
        reduce_mod_poly(d1, d2, r, g, q, w);
        int* temp = r;
        r = g;
        g = temp;
        d1 = poly_degree(d1, r);
        d2 = poly_degree(d2, g);
        remdeg = d2;
        // copy over arrays, don't just change pointers
        for (int i = 0; i <= d1+d2; i++) {
            u2[i] = u1[i];
        }
        for (int i = 0; i <= d1+d2; i++) {
            v2[i] = v1[i];
        }
        for (int i = 0; i <= d1+d2; i++) {
            u1[i] = u[i];
        }
        for (int i = 0; i <= d1+d2; i++) {
            v1[i] = v[i];
        }
        int* q_times_v = malloc(sizeof(int) * (d1+d2+1));
        mult_poly(d1, d2, q, v, q_times_v,w);
        int* q_times_u = malloc(sizeof(int) * (d1+d2+1));
        mult_poly(d1, d2, q, u, q_times_u,w);
        add_poly(d1+d2,d1+d2,u2,q_times_u,u);
        add_poly(d1+d2,d1+d2,v2,q_times_v,v);
        free(q_times_v);
        free(q_times_u);
        free(q);
    }

    for (int i = 0; i <= d2_fixed; i++) {
        A[i] = r[i];
    }
    for (int i = 0; i <= d2_fixed; i++) {
        B[i] = u1[i];
    }
    free(u);
    free(u1);
    free(u2);
    free(v);
    free(v1);
    free(v2);
    free(r);
    free(g);
}


bool invert_poly(const int d, const int poly[d+1], const int dm, const int modulus[dm+1], int result[dm+1], const int w) {
    // inverts poly modulo the polynomial modulus, within F_2^w
    // if invertible, returns true, else returns false.
    memset(result, 0, sizeof(int) * (dm+1));
    int factor[dm+1];
    int gcd[dm+1];
    EEA_standard(d,dm,poly,modulus,w, result, factor, gcd);
    const int degree = poly_degree(dm, gcd);
    if (degree == 0) {
        for (int i = 0; i <= dm; i++) {
            result[i] = galois_single_divide(result[i], gcd[0], w);
        }
        return true;
    }
    else {
        return false;
    }

}

void poly_pow_mod(const int d1, const int d2, const int poly1[d1+1], const int poly2[d2+1], const int exponent, int* result, const int w) {
    // only works for positive exponents. computes poly1^exponent mod poly2.
    int q[d1-d2+1];
    memset(result, 0, sizeof(int)*d2);
    result[0] = 1;
    // poly1_var is reduced poly1, as don't want to change poly1 which is passed by reference/is a pointer
    int poly1_var[d1+1];
    for (int i = 0; i <= d1; i++) {
        poly1_var[i] = poly1[i];
    }
    reduce_mod_poly(d1, d2, poly1_var, poly2, q, w);
    for (int i = 0; i < exponent; i++) {
        int temp[d1+d2+1];
        mult_poly(d1, d2, poly1_var, result, temp, w);
        reduce_mod_poly(d1, d2, temp, poly2, q, w);
        for (int j = 0; j < d2; j++) {
            result[j] = temp[j];
        }
    }

}

// the following are specific algorithms for row/col operations on matrices containing elements in F_2^m interpreted as m-bit vectors

void swap_row(const int i, const int r1, const int j, const int r2, int d1, const int d2, int H[d1][d2]) {
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

void add_row(const int i, const int r1, const int j, const int r2, int d1, const int d2, int H[d1][d2]) {
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