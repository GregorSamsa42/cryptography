#include "sodium.h"
#include <stdio.h>
#include "../galois/galois.h"
#include <math.h>
#include <stdbool.h>
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

void keygen(const int m, const int t) {
    int goppa[t+1];
    memset(goppa, 0, sizeof(goppa));
    randomgoppa(m,t, goppa);
    // generate random perm of all field elements
    int* list_field_perm = random_perm(1 << m);
    // create parity check matrix H = HG * Hhat
    int Hhat[t][1 << m];
    for (int i = 0; i < t; i++) {
        for (int j = 0; j <= (1 << m); j++) {
            Hhat[i][j] = galois_single_divide(galois_pow(list_field_perm[j],i,m), galois_eval_poly(list_field_perm[j],t,goppa, m),m);
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
    printf("okay %b\n", H[0][0]);
    printf("%b\n", H[0][1]);
    printf("%b\n", H[1][0]);
    printf("%b\n", H[1][1]);
    row_reduce(t,1<<m, H, m);
    printf("okay %b\n", H[0][0]);
    printf("%b\n", H[0][1]);
    printf("%b\n", H[1][0]);
    printf("%b\n", H[1][1]);
    // want parity to be of the form (I_mt | Q) where Q is mt x (n-mt)
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
                    int temp = list_field_perm[col];
                    list_field_perm[col] = list_field_perm[next_col];
                    list_field_perm[next_col] = temp;
                    goto end;
                }
            }
        }
        end:
    }
    row_reduce(t,1<<m, H, m);
    printf("%d\n", H[0][0]);
    printf("%d\n", H[0][1]);
    printf("%d\n", H[1][0]);
    printf("%d\n", H[1][1]);
}

int main() {
    if (sodium_init() < 0) {
        /* panic! the library couldn't be initialized; it is not safe to use */
    }
    //
    // int H[2][3] = {0b01, 0b11, 0b10, 0b10, 0b00, 0b11};
    // printf("%d\n", H[0][0]);
    // printf("%d\n", H[0][1]);
    // printf("%d\n", H[1][0]);
    // printf("%d\n", H[1][1]);
    // row_reduce(2,3, H, 2);
    // printf("%d\n", H[0][0]);
    // printf("%d\n", H[0][1]);
    // printf("%d\n", H[1][0]);
    // printf("%d\n", H[1][1]);
    int m = 8;
    int t = 15;
    // generate a random message of 2^m-mt bits to be sent
    bool msg[(1 << m) - m*t];
    keygen(8, 15);
    return 0;
}
