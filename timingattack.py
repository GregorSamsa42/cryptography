from timeit import default_timer as timer
import serial
import numpy
import time

size = 256
t = 47

orginput = input
# configuration for Variscite DART-MX8M-PLUS
ser = serial.Serial(
    port='/dev/ttyUSB0',
    baudrate=115200
)

ser.timeout = None

times = [0]*(size*8)

for i in range(size*8):
    input=b'./decrypt_trigger errors/test' + bytes(str(i),'ascii') + b' mceliece_secret.key'
    ser.write(input + b'\n')
    print(input)
    start = timer()
    ser.read_until(b'@')
    end = timer()
    print("Decryption took " + str(end - start) + " seconds.")
    times[i] = end-start

error_guess = numpy.argsort(times).tolist()[0:t]
print(error_guess)
