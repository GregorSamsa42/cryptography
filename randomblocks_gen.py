import os

fileSizeInBytes = 258 # 2^11/8 + 2
file_number = 1000

for i in range(file_number):
    with open('random_blocks/random_block'+str(i), 'wb') as output:
        output.write(os.urandom(fileSizeInBytes))
