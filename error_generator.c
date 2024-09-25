#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <sys/types.h>
// This file prepares a timed attack on the Mceliece cryptosystem running on a SOM. It takes in a ciphertext
// and generates many error files which can be sent as input, and by measuring their execution time we can deduce the
// cleartext.
// Idea in paper "Side Channels in the McEliece PKC" by Falko Strenzke, Erik Tews, H. Gregor Molter, Raphael Overbeck,
// Abdulhadi Shoufan

void XOR_bits_in_place(unsigned char* a, const int b, const int bit_a, const int bit_b) {
    // XORs the bit_a-th bit of a with the bit_b-th bit of b
    if (bit_a >= bit_b) {
        *a ^= (unsigned char)((b & (1 << bit_b)) << (bit_a - bit_b));
    }
    else {
        *a ^= (unsigned char)((b & (1 << bit_b)) >> (bit_b - bit_a));
    }
}

int main(int argc, char **argv) {
    int size = 256;
    if (argc != 3) {
        printf("There should be three arguments.\n");
        return 1;
    }

    FILE *fp = fopen(argv[1], "rb");
    if (fp == NULL) {
        printf("Value of errno: %d\n", errno);
        perror("Error opening file to be decrypted:");
        return 3;
    }
    // read the file into buffer: it must have size 258 chars
    u_int16_t padding;
    fread(&padding, sizeof(u_int16_t), 1, fp);
    unsigned char* buffer = malloc(((size) * sizeof(char)));
    fread(buffer, sizeof(char), size, fp);
    fclose(fp);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < 8; j++) {
            XOR_bits_in_place(&buffer[i], 1, j, 0);
            FILE* fp_encrypted;
            char buf[20];
            snprintf(buf, sizeof(buf), "%s%s%d", argv[2],"test", 8*i+j);
            fp_encrypted = fopen(buf, "wb");
            fwrite(&padding, sizeof(u_int16_t), 1, fp_encrypted);
            fwrite(buffer, sizeof(char), 256, fp_encrypted);
            XOR_bits_in_place(&buffer[i], 1, j, 0);
        }
    }
    free(buffer);
}