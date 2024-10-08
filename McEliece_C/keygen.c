// Key pair generation for the McEliece cryptosystem
// saves a private key in the format:
// t+1 integers corresponding to the coefficients of the Goppa polynomial, starting with the constant coeff.
// 2^m integers corresponding to the randomly selected private permutation
// saves a public key in the form:
// consecutive rows of the matrix Q of size t * (2^m-m*t) from left to right, each entry is an integer
// Author: Georgi Kocharyan


#include "sodium.h"
#include <stdio.h>
#include "../galois/galois.h"
#include <string.h>

#include "mceliece_supp.h"

int* randomgoppa(const int m, const int t, int goppa[t+1]) {
    // this picks a random polynomial of degree t that has no zeros.
    // Pick t prime, this makes the polynomial guaranteed to be irreducible (Rabin's algorithm for irreducibility).

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
        add_poly(t-1,1,test,X, test);
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
        for (int j = 0; j < (1 << m); j++) {
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
        // printf("\n---H half reduced---\n");
        // for (int row=0; row<t; row++)
        // {
        //     for(int columns=0; columns<(1 << m); columns++)
        //     {
        //         printf("%d     ", H[row][columns]);
        //     }
        //     printf("\n");
        // }
        // printf("---H half reduced---\n");
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
    FILE* fp = fopen("mceliece_public.key","wb");
    for (int i = 0; i < t; i++) {
        for (int j = 0; j < (1 << m)-m*t; j++) {
            Q[i][j] = H[i][j+m*t];
        }
        fwrite(Q[i], sizeof(int), ((1 << m)-m*t), fp);
    }
    fclose(fp);
}
#define m 11
#define t 47 // pick prime
int Q[t][(1 <<m)-m*t];
int main() {
    if (sodium_init() < 0) {
        /* panic! the library couldn't be initialized; it is not safe to use */
        return 1;
    }
    //
    int* goppa = malloc(sizeof(int) * (t+1));
    int* private_perm = malloc(sizeof(int) * (1 << m));
    keygen(m,t,Q,goppa,private_perm);
    FILE* fp = fopen("mceliece_secret.key","wb");
    if (fp == NULL) {
        perror("Error opening file for writing");;
    }
    size_t written = fwrite(goppa, sizeof(int), t+1, fp);
    if (written != (t+1)) {
        perror("Error writing to file");
        fclose(fp);
    }
    written = fwrite(private_perm, sizeof(int), 1 << m, fp);
    if (written != (1 << m)) {
        perror("Error writing to file");
        fclose(fp);
    }
    fclose(fp);
    free(goppa);
    free(private_perm);
}