//
// Created by linux on 09.09.24.
//

#include "mceliece_supp.h"

void encrypt(const int m, const int t, int msg[(1 << m) - m*t], int Q[t][(1 << m)-m*t], int codeword[(1 << m)]) {
    // want to use matrix multiplication but this is inefficient...
    // Q*msg would give the top m*t bits of the codeword (the redundancy), interpreted as t field elements of F_2^m
    // implement quick XOR matrix-vector multiplication

    // must parse the input file.

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
int main(int argc, char *argv[]) {

}