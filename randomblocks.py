import os
import serial
from timeit import default_timer as timer
ser = serial.Serial(
    port='/dev/ttyUSB0',
    baudrate=115200
)

file_number = 1000
ser.timeout = None


for i in range(file_number):
    input=b'./decrypt_trigger random_blocks/random_block' + bytes(str(i),'ascii') + b' mceliece_secret.key'
    ser.write(input + b'\n')
    print(input)
    start = timer()
    ser.read_until(b'@')
    end = timer()
    print("Decryption took " + str(end - start) + " seconds.")
