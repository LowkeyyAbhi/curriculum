import numpy as np

str1 = "GATTACAAGTCC"
str2 = "TTACAGTCA"

match = 2
gap = -1
mismatch = -3

matrix = np.zeros((len(str1)+1, len(str2)+1))

# print(matrix)

#Step 1: Initialisation
for i in range(len(str1)+1):
    matrix[i][0] = i*gap

for j in range(len(str2)+1):
    matrix[0][j] = j*gap

print(matrix)