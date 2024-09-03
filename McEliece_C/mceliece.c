#include "sodium.h"
#include <stdio.h>
#include "../galois/galois.h"
#include <string.h>

#include "mceliece_supp.c"

int* randomgoppa(int m, int t, int goppa[]) {
    while (has_zeroes(goppa, t, m)) {
        for (int i = 0; i <= t; i++) {
            goppa[i] = randombytes_uniform(1 << m);
        }
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
        for (int j = 0; j < (1 << m)-m*t; j++) {
            if (msg[j] == 1) {
                codeword[i] ^= Q[i][j];
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
        linear[1] = galois_shift_inverse(codeword[i],m);
        linear[0] = galois_single_divide(private_perm[i],codeword[i],m);
        int temp[t+1];
        invert_poly(t,linear,t,goppa,temp,m);
        add_poly(t,t,s,temp,s);
    }
    // if it is zero, then there are no errors.
    int sum = 0;
    for (int i = 0; i < (1 << m); i++) {
        sum |= s[i];
    }
    if (sum == 0) {
        for (int i = 0; i <= t; i++) {
            sigma[i] = 0;
        }
        return;
    }
    int T[t+1];
    invert_poly(t,s,t,goppa,T,m);
    // if T is x, then the error locator poly is x
    T[1] ^= 1;
    sum = 0;
    for (int i = 0; i < (1 << m); i++) {
        sum |= T[i];
    }
    if (sum == 0) {
        for (int i = 0; i <= t; i++) {
            sigma[i] = 0;
        }
        sigma[1] = 1;
        return;
    }
    T[1] ^= 1;
    // now find square root r of T(z)+z. Q is a field with 2^(mt) elements, so Frobenius is the identity
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
    const int t = 15;
    int Q[t][(1 << m)-m*t]; int goppa[t+1]; int private_perm[1 << m];
    // generate a random message of 2^m-mt bits to be sent
    int msg[(1 << m) - m*t];
    for (int i = 0; i < (1<<m)-m*t; i++) {
        msg[i] = randombytes_uniform(2);
    }
    keygen(m,t,Q,goppa,private_perm);
    int codeword[1 << m];
    encrypt(m,t,msg,Q,codeword);
    int poly1[1] = {142};
    int poly2[2] = {2,1};
    int result[5] = {0};
    poly_pow_mod(0,1,poly1,poly2,256,result,8);
    for (int i = 0; i < 1; i++) {
        printf("result[%d] = %d\n",i,result[i]);
    }

    return 0;
}
