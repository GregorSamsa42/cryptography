from timeit import default_timer as timer
import serial
# import time

orginput = input
# configuration for Variscite DART-MX8M-PLUS
ser = serial.Serial(
    port='/dev/ttyUSB0',
    baudrate=115200
)

ser.timeout = None

input=b'./keygen'
# send the command to the device
ser.write(input + b'\n')
# wait until we receive the prompt again (which contains an @)
start = timer()
ser.read_until(b'@')
end = timer()
print("Key generation took " + str(end - start) + " seconds.")

# send over file to be decrypted (a block of 258 bytes)
##if 0:
##    ser.write(b'echo -e ' + b'\x27' + b'hello world\x22\x2B' + b'\x27' + b' > test1' + b'\n')
##
##    with open("ciphertext", "rb") as fp:
##        input = fp.read(258);
##        ser.write(b'echo -e ' + b'\x27' + input + b'\x27' + b' > test' + b'\n')
##        print(b'echo -e ' + b'\x27' + input + b'\x27' + b' >> test' + b'\n')
##
##
##for i,b in enumerate(range(90,116)):
##    input = bytes([b])
##    ser.write(b'echo -e ' + b'\x27' + input + b'\x27' + b' > test' + bytes(str(i), "ascii") + b'\n')
##    time.sleep(0.1)
##    ser.write(b"cat test"+ bytes(str(i), "ascii") + b"\n")
##    ser.timeout = 0.1
##    print(ser.read(256))
##    
##    orginput()

input=b'./encrypt readme.md mceliece_public.key'
# send the command to the device
ser.write(input + b'\n')
# wait until we receive the prompt again (which contains an @)
start = timer()
ser.read_until(b'@')
end = timer()
print("Encryption took " + str(end - start) + " seconds.")

input=b'./decrypt encrypted_file_new mceliece_secret.key'
ser.write(input + b'\n')
start = timer()
ser.read_until(b'@')
end = timer()
print("Decryption took " + str(end - start) + " seconds.")
