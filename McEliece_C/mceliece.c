#include "sodium.h"
#include <stdio.h>
#include "../galois/galois.h"
#include <math.h>
#include <stdbool.h>
#include "mceliece_supp.c"

auto randomgoppa(int m, int t, int goppa[]) {
    while (has_zeroes(goppa, t, m)) {
        for (int i = 0; i <= t; i++) {
            goppa[i] = randombytes_uniform(1 << m);
        }
    }
    return goppa;
}

void keygen(int m, int t) {
    int goppa[t+1];
    randomgoppa(m,t, goppa);
    // generate random perm of all field elements
    int* list_field_perm = random_perm(1 << m);
    // create parity check matrix H = HG * Hhat
    int HG[t][t];
    for (int i = 0; i < t; i++) {
        for (int j = 0; j <= i; j++) {
            HG[i][j] = goppa[t-i-j];
        }
    }
    int Hhat[t][1 << m];
    for (int i = 0; i < t; i++) {
        for (int j = 0; j <= (1 << m); j++) {
            Hhat[i][j] = galois_single_divide(galois_pow(list_field_perm[j],i,m), galois_eval_poly(list_field_perm[j],t,goppa, m),m);
        }
    }
    int H[t][1 << m];
    multiply_matrix(HG,Hhat,t,t,1<<m,H);
    // row reduce H while interpreted as a matrix of 0s and 1s
    // have H[m*i+r][j] = (H[i][j] << r) % 1
}

int main()
{
    if (sodium_init() < 0) {
        /* panic! the library couldn't be initialized; it is not safe to use */
    }
    const int m = 8;
    const int t = 15;
    // generate a random message of 2^m-mt bits to be sent
    bool msg[(1 << m) - m*t];
    for (int i = 0; i < m; i++) {
        msg[i] = randombytes_uniform(2);
    }
    printf("The message is %d\n", msg[1]);

    return 0;
}
