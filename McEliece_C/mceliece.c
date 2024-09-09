#include "sodium.h"
#include <stdio.h>
#include "../galois/galois.h"
#include <string.h>

#include "mceliece_supp.h"

// TODO: goppa polys of degree t not prime

int* randomgoppa(const int m, const int t, int goppa[t+1]) {
    // this picks a random polynomial of degree t that has no zeros.
    // Pick t prime, this makes the polynomial guaranteed to be irreducible (Rabin's algorithm for irreducibility).

    int* quotient = malloc(sizeof(int) * (t+1));
    int X[2];
    X[0] = 0;
    X[1] = 1;
    // now X is the identity poly
    int test[t];
    test[0] = 1; // this makes sure test is not 0 in the first run
    goppa[t] = 1;
    // goppa is not necessarily irreducible! But according to Rabin, f of degree t with no zeros is irreducible if and only if it divides X^{2^(mt)}-X
    while (has_zeroes(t, goppa, m) || !fully_zero(t-1,test)) {
        for (int i = 0; i < t; i++) {
            goppa[i] = randombytes_uniform(1 << m);
        }
        memset(test, 0, sizeof(test));
        // note we compute X^{2^(mt)}-X over and over in every run. This is because it is vastly more efficient to immediately compute it mod goppa.
        poly_pow_two_mod(1,t,X,goppa,m*t,test, m);
        add_poly(t,1,test,X, test);
    }
    free(quotient);
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
    // this doesn't necessarily happen immediately (we have only row reduced, not swapped cols)
    // swap columns, swapping the variables in the permutation as well
    int failed_col = 0;
    while (!fully_reduced_parity(t,1<<m,H,m,&failed_col)) {
        printf("\n---H half reduced---\n");
        for (int row=0; row<t; row++)
        {
            for(int columns=0; columns<(1 << m); columns++)
            {
                printf("%d     ", H[row][columns]);
            }
            printf("\n");
        }
        printf("---H half reduced---\n");
        for (int col = failed_col+1; col < (1<<m); col++) {
                if(H[failed_col/m][col] > (1<<((failed_col) % m))) {
                    swap_col(col, failed_col, t, 1<<m, H);
                    const int temp = private_perm[col];
                    private_perm[col] = private_perm[failed_col];
                    private_perm[failed_col] = temp;
                    goto end;
                }
        }
        for (int row = (failed_col / m)+1; row < t; row++) {
            for (int col = failed_col+1; col < (1<<m); col++) {
                if(H[row][col] > 0) {
                    swap_col(col, failed_col, t, 1<<m, H);
                    const int temp = private_perm[col];
                    private_perm[col] = private_perm[failed_col];
                    private_perm[failed_col] = temp;
                    goto end;
                }
            }
        }
        end:
        row_reduce(t,1<<m, H, m);
    }
    // H is now in systematic form. Return Q as the public key.
    for (int i = 0; i < t; i++) {
        for (int j = 0; j < (1 << m)-m*t; j++) {
            Q[i][j] = H[i][j+m*t];
        }
    }
}
void encrypt(const int m, const int t, int msg[(1 << m) - m*t], int Q[t][(1 << m)-m*t], int codeword[(1 << m)]) {
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
    for (int i = m*t; i < (1<<m); i++) {
        codeword[i] = msg[i-m*t];
    }
    // now codeword is the result of Q*msg, and must add error
    int error[1 << m];
    generate_error(1<<m, t, error);
    for (int i = 0; i < (1 << m); i++) {
        codeword[i] ^= error[i];
    }
}

void errorlocator(const int m, const int t, int codeword[1 << m], int goppa[t+1], int private_perm[1 << m], int sigma[t+1]) {
    // Apply Patterson's algorithm to derive the error locator polynomial over F_2^m
    // first compute the syndrome polynomial s
    int s[t];
    memset(s, 0, sizeof(int) * (t+1));
    for (int i = 0; i < (1 << m); i++) {
        int linear[2];
        if (codeword[i] == 1) {
            linear[1] = 1;
            linear[0] = private_perm[i];
        }
        else {
            continue;
        }
        int temp[t+1];
        if (invert_poly(1,linear,t,goppa,temp,m)) {
            add_poly(t-1,t-1,s,temp,s);
        }
    }
    // if it is zero, then there are no errors.
    if (fully_zero(t, s)) {
        for (int i = 0; i <= t; i++) {
            sigma[i] = 0;
        }
        return;
    }
    printf("\n Goppa poly:");
    for (int i = 0; i <= t; i++) {
        printf("%d ",goppa[i]);
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
    memset(A2, 0, sizeof(int) * (t+1));
    // it is guaranteed that A has degree at most t/2 and B at most (t-1)/2.
    mult_poly(t/2,t/2,A,A,A2,m);
    int B2[t+1];
    memset(B2, 0, sizeof(int) * (t+1));
    mult_poly((t-1)/2,(t-1)/2,B,B,B2,m);
    int xB2[t+1];
    memset(xB2, 0, sizeof(int) * (t+1));
    mult_by_var_poly(t-1,1,B2,xB2);

    // the error locating poly sigma is A^2+xB^2
    add_poly(t,t,A2,xB2,sigma);
}

void error_from_errorlocator(const int m, const int t, int sigma[t], int error[1 << m], int private_perm[1 << m]) {
    // given the error locator poly sigma, convert it into the error by plugging in each value consecutively
    // this is very time intensive!!!
    memset(error, 0, sizeof(int) * (1 << m));
    printf("\nLoc.Error:");
    for (int i = 0; i < (1 << m) ; i++) {
        if (galois_eval_poly(private_perm[i], t, sigma, m) == 0) {
            error[i] = 1;
        }
    }
    for (int i = 0; i < (1 << m) ; i++) {
        printf("%d", error[i]);
    }
}

void decrypt(const int m, const int t, int codeword[1 << m], int goppa[t+1], int private_perm[1 << m], int cleartext[(1 << m) - m*t]) {
    int sigma[t+1];
    errorlocator(m, t, codeword, goppa, private_perm, sigma);
    printf("\n Error locator poly:");
    for (int i = 0; i <= t; i++) {
        printf("%d ",sigma[i]);
    }
    int error[1 << m];
    error_from_errorlocator(m,t,sigma,error, private_perm);
    for (int i = 0; i < (1 << m)-m*t; i++) {
        cleartext[i] = codeword[i+m*t]^error[i+m*t];
    }
}

int main(int argc, char** argv) {
    if (sodium_init() < 0) {
        /* panic! the library couldn't be initialized; it is not safe to use */
    }
    //
    const int m = 8;
    const int t = 17; // pick prime
    int Q[t][(1 <<m)-m*t]; int goppa[t+1]; int private_perm[1 << m];
    // generate a random message of 2^m-mt bits to be sent
    int msg[(1 << m) - m*t];
    for (int i = 0; i < (1<<m)-m*t; i++) {
        msg[i] = randombytes_uniform(2);
    }
    printf("  Message:");
    for (int i = 0; i < (1 << m)-m*t; i++) {
        printf("%d",msg[i]);
    }
    keygen(m,t,Q,goppa,private_perm);
    int codeword[1 << m];
    encrypt(m,t,msg,Q,codeword);
    printf("\n Codeword:");
    for (int i = 0; i < (1 << m); i++) {
        printf("%d",codeword[i]);
    }
    int cleartext[(1 << m) - m*t];
    decrypt(m,t,codeword,goppa,private_perm, cleartext);
    printf("\nDecrypted:");
    for (int i = 0; i < (1 << m)-m*t; i++) {
        printf("%d",cleartext[i]);
    }

    bool same = true;
    for (int i = 0; i < (1 << m)-m*t; i++) {
        if (cleartext[i] != msg[i]) {
            same = false;
            printf("\n%d", i);
        }
    }
    if (same) {
        printf("\nmessage and its decryption are the SAME");
    }


    return 0;
}
