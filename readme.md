This project depends on the finite field arithmetic package available at
https://web.eecs.utk.edu/~jplank/plank/papers/CS-07-593/

Place galois.c and galois.h in a folder called galois in the root directory. Note libsodium is required as a dependency as well.

./keygen generates a public and private key pair
./encrypt file mceliece_public.key output encrypts a file using a public key.
./decrypt file mceliece_secret.key output decrypts a file using a private key.
