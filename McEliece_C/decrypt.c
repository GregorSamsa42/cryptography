//
// Created by linux on 11.09.24.
//
#include <errno.h>
#include "mceliece_supp.h"

void errorlocator(const int m, const int t, unsigned char codeword[1 << (m-3)], int goppa[t+1], int private_perm[1 << m], int sigma[t+1]) {
    // Apply Patterson's algorithm to derive the error locator polynomial over F_2^m
    // first compute the syndrome polynomial s
    int s[t];
    memset(s, 0, sizeof(int) * (t+1));
    for (int i = 0; i < (1 << m); i++) {
        int linear[2];
        // i-th bit of codeword is (i%8)-th bit of (i/8)-th char in codeword.
        if (((codeword[i/8] >> (i % 8)) & 1) == 1) {
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
        sigma[0] = 1;
        for (int i = 1; i <= t; i++) {
            sigma[i] = 0;
        }
        return;
    }
    // printf("\n \n Goppa poly:");
    // for (int i = 0; i <= t; i++) {
    //     printf("%d ",goppa[i]);
    // }
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

void error_from_errorlocator(const int m, const int t, int sigma[t], unsigned char error[1 << (m-3)], int private_perm[1 << m]) {
    // given the error locator poly sigma, convert it into the error by plugging in each value consecutively
    // this is very time intensive!!!
    memset(error, 0, sizeof(unsigned char) * (1 << (m-3)));

    for (int i = 0; i < (1 << m) ; i++) {
        if (galois_eval_poly(private_perm[i], t, sigma, m) == 0) {
            // set corresponding bit of error to 1
            error[i/8] |= (1 << (i % 8));
        }
    }
    // printf("\nLoc.Error:");
    // for (int i = 0; i < (1 << (m-3)) ; i++) {
    //     printf("%d ", error[i]);
    // }
}

void decrypt(const int m, const int t, unsigned char* codeword, int goppa[t+1], int private_perm[1 << m], unsigned char** cleartext, int* shift) {
    int sigma[t+1];
    errorlocator(m, t, codeword, goppa, private_perm, sigma);
    // printf("\n Codeword:");
    // for (int i = 0; i < (1 << (m-3)); i++) {
    //     printf("%d ",codeword[i]);
    // }
    // printf("\n Error locator poly:");
    // for (int i = 0; i <= t; i++) {
    //     printf("%d ",sigma[i]);
    // }
    unsigned char error[1 << (m-3)];
    error_from_errorlocator(m,t,sigma,error, private_perm);
    // remove error from the codeword and move pointers appropriately
    for (int i = 0; i < (1 << m)-m*t; i++) {
        // append the last 2^m-mt bits of the output of the decrypted codeword to cleartext
        XOR_bits_in_place(*cleartext, codeword[(m*t+i)/8]^error[(m*t+i)/8], (i+*shift)%8, (m*t+i) % 8);
        if ((i+*shift) % 8 == 7) {
            (*cleartext)++;
        }
    }
    *shift = ((*shift + ((1 << m) - m*t)) % 8);
}

int main(int argc, char *argv[]) {
    // argv[1] is file to be decrypted, argv[2] the private key, and argv[3] the output if given
    if (argc != 3 && argc != 4) {
        printf("There should be two or three arguments.\n");
        return 1;
    }


    const int m = 8; // > 3
    const int t = 11; // prime
    // pubkey has length t(2^m-m*t) integers

    // length of cleartext is unknown

    FILE *fp_ciphertext = fopen(argv[1], "rb");
    if (fp_ciphertext == NULL) {
        printf("Value of errno: %d\n", errno);
        perror("Error opening file to be decrypted:");
        return 3;
    }
    FILE* fp_decrypted;
    if (argc == 4) {
        fp_decrypted = fopen(argv[3], "wb");
    }
    else {
        fp_decrypted = fopen("../decrypted_file", "wb");
    }
    // read the ciphertext into buffer and find out its size
    // note the first two chars of the ciphertext are a 16-bit integer denoting the size of the padding
    fseek(fp_ciphertext, 0, SEEK_END);
    const long file_size = ftell(fp_ciphertext);
    fseek(fp_ciphertext, 0, SEEK_SET);
    u_int16_t padding;
    fread(&padding, sizeof(u_int16_t), 1, fp_ciphertext);
    fseek(fp_ciphertext, 2, SEEK_SET);
    unsigned char* buffer = malloc(((file_size-2) * sizeof(char)));
    memset(buffer, 0, ((file_size-2) * sizeof(char)));
    fread(buffer, file_size, 1, fp_ciphertext);
    fclose(fp_ciphertext);

    // read the private key and split it into its two components goppa poly and private perm
    int goppa[t+1]; int private_perm[1 << m];
    FILE *fp_secretkey = fopen(argv[2], "rb");
    if (fp_secretkey == NULL) {
        printf("Value of errno: %d\n", errno);
        perror("Error opening secret key file:");
        return 2;
    }
    fread(goppa, sizeof(int), t+1, fp_secretkey);
    fseek(fp_secretkey, (t+1)*sizeof(int), SEEK_SET);
    fread(private_perm, 1 <<m, sizeof(int), fp_secretkey);
    fclose(fp_secretkey);
    // printf("Goppa poly:");
    // for (int i = 0; i <= t; i++) {
    //     printf("%d ",goppa[i]);
    // }
    // printf("\n Perm:");
    // for (int i = 0; i < (1<<m); i++) {
    //     printf("%d ",private_perm[i]);
    // }
    // initialise cleartext buffer with the maximal possible length it could attain
    unsigned char cleartext[(((file_size-2)/(1 << (m-3)))*((1<<m)-m*t))/8];
    unsigned char* ptr_text = &cleartext[0];
    unsigned char* ptr_buf = &buffer[0];
    int shift = 0; // start from the 7-th i.e. last bit of a char
    memset(cleartext, 0, sizeof(cleartext));
    // now decrypt the ciphertext 2^{m-3} chars at a time, which is guaranteed to divide the total amount.
    for (int i = 0; i < (file_size-2)/(1 << (m-3)); i++) {
        decrypt(m, t, ptr_buf, goppa, private_perm, &ptr_text, &shift);
        ptr_buf = ptr_buf + (1 << (m-3));
    }
    // write to file, making sure to disregard the padding
    fwrite(cleartext, sizeof(unsigned char), ((file_size-2)/(1 << (m-3)))*((1<<m)-m*t)/8-((padding+7)/8), fp_decrypted);
    free(buffer);
    fclose(fp_decrypted);
}