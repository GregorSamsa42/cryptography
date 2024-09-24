This project depends on the finite field arithmetic package available at
https://web.eecs.utk.edu/~jplank/plank/papers/CS-07-593/

Place galois.c and galois.h in a folder called galois in the root directory. Note libsodium is required as a dependency as well.

./keygen generates a public and private key pair
./encrypt file mceliece_public.key output encrypts a file using a public key.
./decrypt file mceliece_secret.key output decrypts a file using a private key.


--- 


zum cross-compilen für das Variscite DART-MX8M-PLUS verwende:

aarch64-linux-gnu-gcc -static ../McEliece_C/keygen.c ../McEliece_C/mceliece_supp.c ../galois/galois.c -L../lib -lsodium -o keygen -O2

bzw mit encrypt/decrypt statt keygen.

zum Verbinden nutze sudo picocom -b 115200 /dev/ttyUSBx, je nachdem wie das dev heißt (sudo dmesg gibt das aus)
username: root
password: pw1234


Zum debuggen: clang-15 -fsanitize=address,undefined -g ../galois/galois.c  mceliece_supp.c keygen.c /usr/lib/x86_64-linux-gnu/libsodium.a
./a.out

Die version mit Trigger cross-compilet wie folgt:

aarch64-linux-gnu-gcc -static ../McEliece_C/decrypt.c ../McEliece_C/mceliece_supp.c ../galois/galois.c ../lib/libgpiod-2.1.3/lib/chip.c ../lib/libgpiod-2.1.3/lib/line-config.c ../lib/libgpiod-2.1.3/lib/line-settings.c ../lib/libgpiod-2.1.3/lib/request-config.c ../lib/libgpiod-2.1.3/lib/internal.c ../lib/libgpiod-2.1.3/lib/chip-info.c ../lib/libgpiod-2.1.3/lib/line-info.c ../lib/libgpiod-2.1.3/lib/line-request.c ../lib/libgpiod-2.1.3/lib/info-event.c ../lib/libgpiod-2.1.3/lib/edge-event.c -I../lib/libgpiod-2.1.3/include -L../lib -lsodium -o decrypt
