#include "sodium.h"
#include <stdio.h>
#include "../galois/galois.h"
#include <string.h>

#include "mceliece_supp.c"

// TODO: goppa polys of degree t not prime

int* randomgoppa(int m, int t, int goppa[]) {
    // this picks a random polynomial of degree t that has no zeros.
    // Pick t prime, this makes the polynomial guaranteed to be irreducible (Rabin's algorithm for irreducibility).

    int* quotient = malloc(sizeof(int) * (t+1));
    int X[t+1];
    memset(X, 0, sizeof(X));
    X[1] = 1;
    // now X is the identity poly
    int test[t];
    test[0] = 1; // this makes sure test is not 0 in the first run
    goppa[t] = 1;
    // goppa is not necessarily irreducible! But according to Rabin, f of degree t with no zeros is irreducible if and only if it divides X^{2^(mt)}-X
    while (has_zeroes(goppa, t, m) || !fully_zero(t-1,test)) {
        for (int i = 0; i < t; i++) {
            goppa[i] = randombytes_uniform(1 << m);
        }
        // note we compute X^{2^(mt)}-X over and over in every run. This is because it is vastly more efficient to immediately compute it mod goppa.
        poly_pow_two_mod(t,t,X,goppa,m*t,test, m);
        add_poly(t,t,X,test, test);
    }
    free(quotient);
    for (int i = 0; i < t; i++) {
     //    printf("%d ", goppa[i]);
    }
    return goppa;
}

void keygen(const int m, const int t, int Q[t][(1 << m)-m*t], int goppa[t+1], int private_perm[1 << m]) {
    // generates public and private keys. Q is the public key, and private_perm and private
    memset(goppa, 0, sizeof(int) * (t+1));
    randomgoppa(m,t, goppa);
    // generate random perm of all field elements
    random_perm(1 << m, private_perm);
    // create parity check matrix H = HG * Hhat
    int Hhat[t][1 << m];
    for (int i = 0; i < t; i++) {
        for (int j = 0; j <= (1 << m); j++) {
            Hhat[i][j] = galois_single_divide(galois_pow(private_perm[j],i,m), galois_eval_poly(private_perm[j],t,goppa, m),m);
        }
    }
    int HG[t][t];
    memset(HG, 0, sizeof(HG));
    for (int i = 0; i < t; i++) {
        for (int j = 0; j <= i; j++) {
            HG[i][j] = goppa[t-i+j];
        }
    }
    int H[t][1 << m];
    multiply_matrix(t, t, 1<<m, HG,Hhat,H, m);
    // row reduce H while interpreted as a matrix of 0s and 1s
    // have H[m*i+r][j] = (H[i][j] << r) % 1
    row_reduce(t,1<<m, H, m);
    // want parity to be of the form (I_mt | Q) where Q is mt x (n-mt) over F_2 and t x (n-mt) over F_2^m
    // this doesn't necessarily happen immediately
    // swap columns, swapping the variables in the permutation as well
    for (int col = 0; col < (1 << m); col++) {
        for (int row = 0; row < t; row++) {
            if (H[row][col] == 1) {
                break;
            }
        }
        for (int next_col = col+1; next_col < t; next_col++) {
            for (int next_row = 0; next_row < t; next_row++) {
                if (H[next_row][next_col] == 1) {
                    swap_col(col, next_col, t, 1<<m, H);
                    const int temp = private_perm[col];
                    private_perm[col] = private_perm[next_col];
                    private_perm[next_col] = temp;
                    goto end;
                }
            }
        }
        end:
    }
    row_reduce(t,1<<m, H, m);
    // H is now in systematic form. Return Q as the public key.
    for (int i = 0; i < t; i++) {
        for (int j = 0; j < (1 << m)-m*t; j++) {
            Q[i][j] = H[i][j+m*t];
        }
    }

}
void encrypt(int m, int t, int msg[(1 << m) - m*t], int Q[t][(1 << m)-m*t], int codeword[(1 << m)]) {
    // want to use matrix multiplication but this is inefficient...
    // Q*msg would give the top m*t bits of the codeword (the redundancy), interpreted as t field elements of F_2^m
    memset(codeword, 0, sizeof(int) * (1 << m));
    // implement quick XOR matrix-vector multiplication
    for (int i = 0; i < t; i++) {
        for (int k = 0; k < m; k++) {
            for (int j = 0; j < (1 << m)-m*t; j++) {
                if (msg[j] == 1) {
                    codeword[i*m+k] ^= (Q[i][j] >> k) & 1;
                }
            }
        }
    }
    for (int i = t; i < (1<<m); i++) {
        codeword[i] = msg[i];
    }
    // now codeword is the result of Q*msg, and must add error
    int error[1 << m];
    generate_error(1<<m, t, error);
    for (int i = 0; i < (1 << m); i++) {
        codeword[i] ^= error[i];
    }
}

void errorlocator(int m, int t, int codeword[1 << m], int goppa[t+1], int private_perm[1 << m], int sigma[t+1]) {
    // Apply Patterson's algorithm to derive the error locator polynomial over F_2^m
    // first compute the syndrome polynomial s
    int s[t+1];
    memset(s, 0, sizeof(int) * (t+1));
    for (int i = 0; i < (1 << m); i++) {
        int linear[t+1];
        memset(linear, 0, sizeof(int) * (t+1));
        if (codeword[i] == 1) {
            linear[1] = 1;
            linear[0] = private_perm[i];
        }
        else {
            continue;
        }
        int temp[t+1];
        if (invert_poly(1,linear,t,goppa,temp,m)) {
            add_poly(t,t,s,temp,s);
        }
    }
    // if it is zero, then there are no errors.
    if (fully_zero(t, s)) {
        for (int i = 0; i <= t; i++) {
            sigma[i] = 0;
        }
        return;
    }
    int T[t+1];
    invert_poly(t,s,t,goppa,T,m);
    // if T = s^{-1} is x, then the error locator poly is x
    T[1] ^= 1;
    // T has been set to s^{-1}+x
    if (fully_zero(t, T)) {
        for (int i = 0; i <= t; i++) {
            sigma[i] = 0;
        }
        sigma[1] = 1;
        return;
    }
    // now find square root r of T(x).
    int r[t+1];
    sqrt_poly(t,t,T,goppa,r,m);
    // the polys in following have smaller degree in theory, but we give them degree t to be sure
    int A[t+1];
    int B[t+1];
    EEA_patterson(t,t,r,goppa,A,B,m);
    int A2[t+1];
    mult_poly(t,t,A,A,A2,m);
    int B2[t+1];
    mult_poly(t,t,B,B,B2,m);
    int xB2[t+1];
    mult_by_var_poly(t,1,B2,xB2);
    // the error locating poly sigma is A^2+xB^2
    add_poly(t,t,A2,xB2,sigma);
}

void error_from_errorlocator(int m, int t, int sigma[t], int error[1 << m]) {
    // given the error locator poly sigma, convert it into the error by plugging in each value consecutively
    // this is very time intensive!!!
    memset(error, 0, sizeof(int) * (1 << m));
    for (int i = 0; i < (1 << m) ; i++) {
        if (galois_eval_poly(i, t, sigma, m) == 0) {
            error[i] = 1;
        }
    }
}

int main() {
    if (sodium_init() < 0) {
        /* panic! the library couldn't be initialized; it is not safe to use */
    }
    //
    const int m = 8;
    const int t = 11; // segfaults for t > 11, why ??
    int Q[t][(1 << m)-m*t]; int goppa[t+1]; int private_perm[1 << m];
    // generate a random message of 2^m-mt bits to be sent
    int msg[(1 << m) - m*t];
    for (int i = 0; i < (1<<m)-m*t; i++) {
        msg[i] = randombytes_uniform(2);
    }
    keygen(m,t,Q,goppa,private_perm);
    int codeword[1 << m];
    encrypt(m,t,msg,Q,codeword);
    int sigma[t+1];
    errorlocator(m,t,codeword,goppa,private_perm, sigma);
    for (int i = 0; i < t; i++) {
        printf("%d ",sigma[i]);
    }
    int poly1[12] = {1,34,12,44,2,144,11,241,2,4,1,4};
    int poly2[2] = {2,4};

    return 0;
}
