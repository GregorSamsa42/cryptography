
#include "mceliece_supp.h"
#include "sodium.h"



void XOR_bits_in_place(unsigned char* a, const int b, const int bit_a, const int bit_b) {
    // XORs the bit_a-th bit of a with the bit_b-th bit of b
    if (bit_a >= bit_b) {
        *a ^= (unsigned char)((b & (1 << bit_b)) << (bit_a - bit_b));
    }
    else {
        *a ^= (unsigned char)((b & (1 << bit_b)) >> (bit_b - bit_a));
    }
}

void swap_bytes(unsigned char* a, unsigned char* b, const int bit_a, const int bit_b) {
    if (((*a >> bit_a) & 1) != ((*b >> bit_b) & 1)) {
        *a ^= (1 << bit_a);
        *b ^= (1 << bit_b);
    }
}

bool fully_zero(const int d, const int poly[d+1]) {
    // tests if a polynomial is zero
    int sum = 0;
    for (int i = 0; i < d; i++) {
        sum |= poly[i];
    }
    if (sum == 0) {
        return true;
        }
    else {
        return false;
    }

}

int galois_pow(const int a, const int b, const int w) {
    int result = 1;
    for (int i = 0; i < b; i++) {
        result = galois_single_multiply(result, a, w);
    }
    return result;
}


void random_perm(const int n, int perm[n]) {
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

void generate_error(const int n, const int w, unsigned char error[n]) {
    // create n-char error vector with Hamming weight w
    // must have w <= n/8, then error is saved as chars
    // first, create one with w 1s at the beginning
    // error needs n allocated ints
    memset(error, 0, sizeof(char) * n);
    for (int i = 0; i < w/8; i++) {
        error[i] = 255;
    }
    error[w/8] = (1 << (w%8)) - 1;
    // now shuffle error using the Fisher-Yates shuffle
    for (int i = 8*n-1; i > 0; i--) {
        const int j = randombytes_uniform(i+1);
        swap_bytes(&error[i/8], &error[j/8], i % 8, j % 8);
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


int galois_eval_poly(const int a, const int d, const int poly[d+1], const int w) {
    // degree is degree of poly: due to type decay this information is lost after passing
    int result = 0;
    for (int bit = 0; bit <= d; bit++) {
        int temp = poly[bit];
        temp = galois_single_multiply(temp, galois_pow(a, bit, w), w);
        result = result ^ temp;
    }
    return result;
}

bool has_zeroes(const int d, const int poly[d+1], const int w) {
    for (int i = 0; i < (1 << w); i++) {
        if (galois_eval_poly(i, d, poly, w) == 0) {
            return true;
        }
    }
    return false;
}

void add_poly(const int d1, const int d2, const int poly1[d1+1], const int poly2[d2+1], int* result) {
    // adds two polys with coefficients in F_2^m (actually for any m, not just the fixed one here)
    // we will only ever call this method with d1 > d2, so assume this is the case for simplicity.
    // result should have d1+1 integers allocated
    // this is particularly convenient as can let poly1 = result (not so for other methods)

    for (int i = 0; i <= d2; i++) {
        result[i] = poly1[i]^poly2[i];
    }
    for (int i = d2+1; i <= d1; i++) {
        result[i] = poly1[i];
    }
}

void mult_by_var_poly(const int d, const int n, const int poly1[d+1], int* result) {
    // multiplies p(x) by x^^n
    // result must have d+n+1 integers allocated already (but this is guaranteed by mult_poly)
    memset(result, 0, sizeof(int) * (d+n+1));
    for (int i = n; i <= d+n; i++) {
        result[i] = poly1[i-n];
    }
}

void mult_poly(const int d1, const int d2, const int poly1[d1+1], const int poly2[d2+1], int* result, const int w) {
    // result needs d1+d2+1 * sizeof(int) allocated space
    for (int j = 0; j <= d2; j++) {
        int* temp = malloc(sizeof(int) * (d1+d2+1));
        memset(temp, 0, sizeof(int) * (d1+d2+1));
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
    if (d2 == -1 || d2 > d1) {
        quotient = NULL;
        return;
    }
    for (int i = d1; i >= d2; i--) {
        // create temp poly that is poly2 scaled appropriately, and subtract it from poly1
        int temp[i+1];
        mult_by_var_poly(d2, i-d2, poly2, temp);
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
        add_poly(d1-d2, i-d2, quotient, temp_quotient, quotient);
        // subtract temp off of poly1
        add_poly(d1,i, poly1, temp, poly1);
    }
}

void EEA_standard(const int d1, const int d2, const int poly1[d1+1], const int poly2[d2+1], const int w, int* factor1, int* factor2, int* gcd) {
    // initialise variables in Euclid's algorithm. Want r - q*g = rem = u*poly1 + v*poly2.
    // again for simplicity wlog d1 <= d2. If not code will not work but this is OK - we are only using it in this way
    // in practice just reduce_mod the first poly before inputting
    // factor1,factor2,gcd require d2+1 allocated integers
    int* r = malloc(sizeof(int) * (d2+1));
    memset(r, 0, sizeof(int) * (d2+1));
    for (int i = 0; i <= d1; i++) {
        r[i] = poly1[i];
    }
    int* g = malloc(sizeof(int) * (d2+1));
    memset(g, 0, sizeof(int) * (d2+1));
    for (int i = 0; i <= d2; i++) {
        g[i] = poly2[i];
    }
    // note that for the Bezout coefficients we know they have degrees bounded by d2 in every step.
    int* u = malloc(sizeof(int) * (d2+1));
    memset(u, 0, sizeof(int) * (d2+1));
    int* u1 = malloc(sizeof(int) * (d2+1));
    memset(u1, 0, sizeof(int) * (d2+1));
    u1[0] = 1;
    int* v = malloc(sizeof(int) * (d2+1));
    memset(v, 0, sizeof(int) * (d2+1));
    v[0] = 1;
    int* v1 = malloc(sizeof(int) * (d2+1));
    memset(v1, 0, sizeof(int) * (d2+1));
    int* u2 = malloc(sizeof(int) * (d2+1));
    memset(u2, 0, sizeof(int) * (d2+1));
    int* v2 = malloc(sizeof(int) * (d2+1));
    memset(v2, 0, sizeof(int) * (d2+1));
    int remdeg = 1;
    // terminate when the remainder is zero (i.e. degree -1), then can read off Bezout coefficients from u and v
    while (remdeg > -1) {
        int* q = malloc(sizeof(int) * (d2+1));
        memset(q, 0, sizeof(int) * (d2+1));
        reduce_mod_poly(d2, d2, r, g, q, w);
        int* temp = malloc(sizeof(int) * (d2+1));
        memset(temp, 0, sizeof(int) * (d2+1));
        for (int j = 0; j <= d2; j++) {
            temp[j] = r[j];
        }
        for (int j = 0; j <= d2; j++) {
            r[j] = g[j];
        }
        for (int j = 0; j <= d2; j++) {
            g[j] = temp[j];
        }
        remdeg = poly_degree(d2, g);
        // copy over arrays, don't just change pointers
        for (int i = 0; i <= d2; i++) {
            u2[i] = u1[i];
        }
        for (int i = 0; i <= d2; i++) {
            v2[i] = v1[i];
        }
        for (int i = 0; i <= d2; i++) {
            u1[i] = u[i];
        }
        for (int i = 0; i <= d2; i++) {
            v1[i] = v[i];
        }
        // must allocate more room for q_times_v because of how mult_poly works
        int* q_times_v = malloc(sizeof(int) * (2*d2+1));
        memset(q_times_v, 0, sizeof(int) * (2*d2+1));
        mult_poly(d2, d2, q, v, q_times_v,w);
        int* q_times_u = malloc(sizeof(int) * (2*d2+1));
        memset(q_times_u, 0, sizeof(int) * (2*d2+1));
        mult_poly(d2, d2, q, u, q_times_u,w);
        // but we know that in reality q_times_u and q_times_v have degree at most d2.
        add_poly(d2,d2,u2,q_times_u,u);
        add_poly(d2,d2,v2,q_times_v,v);
        free(q);
        free(q_times_v);
        free(q_times_u);
        free(temp);
    }
    for (int i = 0; i <= d2; i++) {
        gcd[i] = r[i];
    }
    for (int i = 0; i <= d2; i++) {
        factor1[i] = u1[i];
    }
    for (int i = 0; i <= d2; i++) {
        factor2[i] = v1[i];
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
    void EEA_patterson(const int d1, const int d2, const int poly1[d1+1], const int goppa[d2+1], int* A, int* B, const int w) {
    // initialise variables in Euclid's algorithm. Want r - q*g = rem = u*poly1 + v*poly2.
    // same as before, except the break condition and the outputs are different
    // factor1,factor2,gcd require d2+1 allocated integers
    int* r = malloc(sizeof(int) * (d2+1));
    memset(r, 0, sizeof(int) * (d2+1));
    for (int i = 0; i <= d1; i++) {
        r[i] = poly1[i];
    }
    int* g = malloc(sizeof(int) * (d2+1));
    memset(g, 0, sizeof(int) * (d2+1));
    for (int i = 0; i <= d2; i++) {
        g[i] = goppa[i];
    }
    // note that for the Bezout coefficients we know they have degrees bounded by d2 in every step.
    int* u = malloc(sizeof(int) * (d2+1));
    memset(u, 0, sizeof(int) * (d2+1));
    int* u1 = malloc(sizeof(int) * (d2+1));
    memset(u1, 0, sizeof(int) * (d2+1));
    u1[0] = 1;
    int* v = malloc(sizeof(int) * (d2+1));
    memset(v, 0, sizeof(int) * (d2+1));
    v[0] = 1;
    int* v1 = malloc(sizeof(int) * (d2+1));
    memset(v1, 0, sizeof(int) * (d2+1));
    int* u2 = malloc(sizeof(int) * (d2+1));
    memset(u2, 0, sizeof(int) * (d2+1));
    int* v2 = malloc(sizeof(int) * (d2+1));
    memset(v2, 0, sizeof(int) * (d2+1));
    int remdeg = d2;
    // terminate when the remainder is zero (i.e. degree -1), then can read off Bezout coefficients from u and v
    while (remdeg > d2/2 || poly_degree(d2,u1) > (d2-1)/2) {
        int* q = malloc(sizeof(int) * (d2+1));
        memset(q, 0, sizeof(int) * (d2+1));
        reduce_mod_poly(d2, d2, r, g, q, w);
        int* temp = malloc(sizeof(int) * (d2+1));
        memset(temp, 0, sizeof(int) * (d2+1));
        for (int j = 0; j <= d2; j++) {
            temp[j] = r[j];
        }
        for (int j = 0; j <= d2; j++) {
            r[j] = g[j];
        }
        for (int j = 0; j <= d2; j++) {
            g[j] = temp[j];
        }
        remdeg = poly_degree(d2, g);
        // copy over arrays, don't just change pointers
        for (int i = 0; i <= d2; i++) {
            u2[i] = u1[i];
        }
        for (int i = 0; i <= d2; i++) {
            v2[i] = v1[i];
        }
        for (int i = 0; i <= d2; i++) {
            u1[i] = u[i];
        }
        for (int i = 0; i <= d2; i++) {
            v1[i] = v[i];
        }
        // must allocate more room for q_times_v because of how mult_poly works
        int* q_times_v = malloc(sizeof(int) * (2*d2+1));
        memset(q_times_v, 0, sizeof(int) * (2*d2+1));
        mult_poly(d2, d2, q, v, q_times_v,w);
        int* q_times_u = malloc(sizeof(int) * (2*d2+1));
        memset(q_times_u, 0, sizeof(int) * (2*d2+1));
        mult_poly(d2, d2, q, u, q_times_u,w);
        // but we know that in reality q_times_u and q_times_v have degree at most d2.
        add_poly(d2,d2,u2,q_times_u,u);
        add_poly(d2,d2,v2,q_times_v,v);
        free(q);
        free(q_times_v);
        free(q_times_u);
        free(temp);
    }

    for (int i = 0; i <= d2; i++) {
        A[i] = g[i];
    }
    for (int i = 0; i <= d2; i++) {
        B[i] = u[i];
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


bool invert_poly(const int d, int poly[d+1], const int dm, const int modulus[dm+1], int result[dm], const int w) {
    // inverts poly modulo the polynomial modulus, within F_2^w
    // if invertible, returns true, else returns false.
    // memset(result, 0, sizeof(int) * (dm+1));
    int factor[dm+1];
    int gcd[dm+1];
    memset(gcd, 0, sizeof(int) * (dm+1));
    memset(factor, 0, sizeof(int) * (dm+1));
    // following line isn't really needed in our implementation but might be needed in theory
    // reduce_mod_poly(d,dm,poly,modulus,factor,w);
    EEA_standard(d,dm,poly,modulus,w, result, factor, gcd);
    const int degree = poly_degree(dm, gcd);
    if (degree == 0) {
        for (int i = 0; i < dm; i++) {
            result[i] = galois_single_divide(result[i], gcd[0], w);
        }
        return true;
    }
    else {
        return false;
    }

}

void poly_pow_two_mod(const int d1, const int d2, int poly1[d1+1], const int poly2[d2+1], const int exponent, int* result, const int w) {
    // only works for positive exponents. computes poly1^(2^exponent) mod poly2.
    // result must have d2 integers allocated
    // note poly1 is changed here, but this does not pose a problem for our implementations
    int* q = malloc(sizeof(int) * (d1+d2+1));
    memset(q, 0, sizeof(int) * (d1+d2+1));
    memset(result, 0, sizeof(int) * d2);
    reduce_mod_poly(d1, d2, poly1, poly2, q, w);
    for (int i = 0; i <= d1; i++) {
        result[i] = poly1[i];
    }
    for (int i = 0; i < exponent; i++) {
        int temp[2*d2-1];
        memset(temp, 0, sizeof(int) * (2*d2-1));
        mult_poly(d2-1, d2-1, result, result, temp, w);
        reduce_mod_poly(2*d2-2, d2, temp, poly2, q, w);
        for (int j = 0; j < d2; j++) {
            result[j] = temp[j];
        }
    }
    free(q);

}

void sqrt_poly(const int d1, const int d2, int poly1[d1+1], const int poly2[d2+1], int* result, const int w) {
    // this works if poly2 is irreducible over F_2^w.
    // Have a field with 2^(mt) elements, so Frobenius is the identity
    // result must have d2 integers allocated
    // again note poly1 is not constant.
    int t = poly_degree(d2, poly2);
    poly_pow_two_mod(d1,t,poly1,poly2,w*t-1, result, w);
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

bool fully_reduced_parity(const int d1, const int d2, const int H[d1][d2], const int m, int* failed_col) {
    bool fully_reduced = true;
    for (int col = 0; col < d1*m; col++) {
        if ((H[col / m][col] >> (col % m) & 1) != 1) {
            fully_reduced = false;
            *failed_col = col;
            break;
        }
    }
    return fully_reduced;
}