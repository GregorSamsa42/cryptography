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

    unsigned char* temp = *msg;
    for (int i = 0; i < t; i++) {
        for (int k = 0; k < m; k++) {
            temp = *msg;
            for (int j = *shift; j < (1 << m)-m*t+*shift; j++) {
                // maybe 8-j % 8 ??
                if (((*temp >> ((j % 8))) & 1) == 1) {
                    XOR_bits_in_place(&codeword[(i*m+k)/8],Q[i][j-*shift], (i*m+k) % 8, k);
                }
                if (j % 8 == 7) {
                    temp++;
                }
            }
        }
    }
    for (int i = m*t; i < (1<<m); i++) {
        XOR_bits_in_place(&codeword[i/8],(*msg)[(i-m*t+*shift)/8], (i) % 8, (i-m*t+*shift) % 8);
    }
    *msg = temp;
    *shift = ((*shift + ((1 << m) - m*t)) % 8);
    // now codeword is the result of Q*msg, and must add error
    unsigned char error[1 << (m-3)];
    memset(error,0,sizeof(error));
    generate_error(1<<(m-3), t, error);
    for (int i = 0; i < (1 << (m-3)); i++) {
        codeword[i] ^= error[i];
    }
    printf("Codeword: \n");
    for (int i = 0; i < (1 << (m-3)); i++) {
        printf("%d ", codeword[i]);
    }
    printf("Error: \n");
for (int i = 0; i < (1 << (m-3)); i++) {
    printf("%d ", error[i]);
}
}
int main(int argc, char *argv[]) {
 // argv[1] is the file to be encrypted, argv[2] is the public key. If argv[3] is supplied, this is the name of the output.
    if (argc != 3 && argc != 4) {
        printf("There should be two or three arguments.\n");
        return 1;
    }

    const int m = 11; // > 3
    const int t = 50; // prime
    // pubkey has length t(2^m-m*t) integers

    // length of cleartext is unknown

    FILE* fp_cleartext = fopen(argv[1], "rb");
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
    const u_int16_t padding = ((1 << m) - m*t) - ((8*file_size) % ((1 << m) - m*t));
    // read the public key into a matrix
    FILE *fp_pubkey = fopen(argv[2], "rb");
    if (fp_pubkey == NULL) {
        printf("Value of errno: %d\n", errno);
        perror("Error opening public key file:");
        return 2;
    }
    int Q[t][(1<<m) - m*t];
    for (int i = 0; i < t; i++) {
        fseek(fp_pubkey, sizeof(int)*i*((1<<m)-m*t), SEEK_SET);
        fread(Q[i], sizeof(int), (1<<m)-m*t, fp_pubkey);
    }
    fclose(fp_pubkey);
    // the first 16bit integer in fp_encrypted denotes the padding, and overwrites file if it exists by reading
    fwrite(&padding, sizeof(u_int16_t), 1, fp_encrypted);
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
    // can compute the amount of required chunks easily from the file size and from m,t
    unsigned char* ptr = &buffer[0];
    for (int i = 0; i < 8*file_size/((1<<m) - m*t)+1; i++) {
        unsigned char codeword[(1<<(m-3))];
        memset(codeword, 0, sizeof(codeword));
        encrypt(m, t, &ptr, Q, codeword, &shift);
        fwrite(codeword, sizeof(char), (1<<(m-3)), fp_encrypted);
    }
    fclose(fp_encrypted);
    free(buffer);
}