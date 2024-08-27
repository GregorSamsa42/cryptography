cmake_minimum_required(VERSION 3.29)
project(cryptography C)

set(CMAKE_C_STANDARD 11)

include_directories(galois)

add_executable(cryptography
    galois/galois.c
    galois/galois.h
    McEliece_C/mceliece.c)

target_link_libraries(cryptography libsodium.a -lm)


