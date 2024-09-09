cmake_minimum_required(VERSION 3.29)
project(cryptography C)

set(CMAKE_C_STANDARD 11)

include_directories(galois)

add_executable(encrypt
        galois/galois.c
        galois/galois.h
        McEliece_C/encrypt.c
        McEliece_C/mceliece_supp.h
        McEliece_C/mceliece_supp.c)

add_executable(keygen
        galois/galois.c
        galois/galois.h
        McEliece_C/keygen.c
            McEliece_C/mceliece_supp.h
            McEliece_C/mceliece_supp.c)


target_link_libraries(keygen libsodium.a -lm)

