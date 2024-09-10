//
// Created by linux on 09.09.24.
//
#include <errno.h>
#include "mceliece_supp.h"



void encrypt(const int m, const int t, unsigned char** msg, int Q[t][(1 << m)-m*t], unsigned char codeword[(1 << (m-3))], int* shift) {
    // want to use matrix multiplication but this is inefficient...
    // Q*msg would give the top m*t bits of the codeword (the redundancy), interpreted as t field elements of F_2^m
    // implement quick XOR matrix-vector multiplication
    // a char is 8 bits, so need a char of length 2^{m-3} for a bit string of length 2^m

    // must parse the input buffer: this is done by reading the chars starting from the passed pointer, shifted by 0 <= shift < 8 bits.

    for (int i = 0; i < t; i++) {
        for (int k = 0; k < m; k++) {
            for (int j = *shift; j < (1 << m)-m*t; j++) {
                if ((**msg >> (8-(j % 8))) == 1) {
                    XOR_bits(&codeword[(i*m+k)/8],Q[i][j], (i*m+k) % 8, k);
                }
                if (j % 8 == 7) {
                    (*msg)++;
                }
            }
        }
    }
    *shift = ((*shift + ((1 << m) - m*t)) % 8);
    for (int i = m*t; i < (1<<m); i++) {
        XOR_bits(&codeword[i/8],*msg[(i-m*t)/8], (i) % 8, (i-m*t) % 8);
    }
    // now codeword is the result of Q*msg, and must add error
    unsigned char error[1 << (m-3)];
    generate_error(1<<(m-3), t, error);
    for (int i = 0; i < (1 << (m-3)); i++) {
        codeword[i] ^= error[i];
    }
}
int main(int argc, char *argv[]) {
 // argv[1] is the file to be encrypted, argv[2] is the public key. If argv[3] is supplied, this is the name of the output.
    if (argc != 3 && argc != 4) {
        printf("There should be two or three arguments.\n");
        return 1;
    }
    FILE *fp_pubkey = fopen(argv[2], "rb");
    if (fp_pubkey == NULL) {
        printf("Value of errno: %d\n", errno);
        perror("Error opening public key file:");
        return 2;
    }

    const int m = 6; // > 3
    const int t = 7; // prime
    // pubkey has length t(2^m-m*t) integers

    // length of cleartext is unknown

    FILE *fp_cleartext = fopen(argv[1], "rb");
    if (fp_cleartext == NULL) {
        printf("Value of errno: %d\n", errno);
        perror("Error opening file to be encrypted:");
        return 3;
    }
    FILE* fp_encrypted;
    if (argc == 4) {
        fp_encrypted = fopen(argv[3], "wb");
    }
    else {
        fp_encrypted = fopen("../encrypted_file", "wb");
    }
    // read the cleartext into buffer and find out its size
    fseek(fp_cleartext, 0, SEEK_END);
    const long file_size = ftell(fp_cleartext);
    fseek(fp_cleartext, 0, SEEK_SET);
    // allocate enough space for the buffer to store the file plus padding
    // file_size is the number of chars. Add 2^m-mt bits to guarantee there are enough zeros at the end
    // must store how much padding we required for decryption!
    unsigned char* buffer = malloc((file_size+ (1 << m - m*t)/8) * sizeof(char));
    memset(buffer, 0, (file_size+ (1 << m - m*t)/8) * sizeof(char));
    fread(buffer, file_size, 1, fp_cleartext);
    fclose(fp_cleartext);
    int padding = (8*file_size) % ((1 << m) - m*t);
    // read the public key into a matrix
    int Q[t][(1<<m) - m*t];
    for (int i = 0; i < t; i++) {
        fseek(fp_pubkey, i*((1<<m)-m*t), SEEK_SET);
        fread(Q[i], sizeof(int), (1<<m)-m*t, fp_pubkey);
    }
    fclose(fp_pubkey);
    // the first integer in fp_encrypted denotes the size of the padding
    fwrite(&padding, sizeof(int), 1, fp_encrypted);
    fclose(fp_encrypted);
    // now reopen as append to add the encrypted data.
    if (argc == 4) {
        fp_encrypted = fopen(argv[3], "ab");
    }
    else {
        fp_encrypted = fopen("../encrypted_file", "ab");
    }
    int shift = 0;
    // split the file to be encrypted into chunks of (2^m-mt) bits
    for (int i = 0; i < 8*file_size/((1<<m) - m*t)+1; i++) {
        unsigned char codeword[(1<<(m-3))];
        memset(codeword, 0, sizeof(codeword));
        unsigned char* ptr = &buffer[0];
        encrypt(m, t, &ptr, Q, codeword, &shift);
        fwrite(codeword, sizeof(char), (1<<(m-3)), fp_encrypted);
    }
    fclose(fp_encrypted);
    free(buffer);
}