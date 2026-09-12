import numpy as np

str1 = "GATTACAAGTCC" # DNA sequence 1
str2 = "TTACAGTCA" #DNA sequence 2

match = 3 #match score as provided in assignment pdf
mismatch = -2 #mismatch penalty as provided in assignment pdf
gap = -2 #gap penalty score as provided in assignment pdf

matrix = np.zeros((len(str1)+1, len(str2)+1)) # Create a matrix through numpy library

#Step 1: Initialisation of the matrix

# Set first row and first column as 0
for i in range(len(str1)+1):
    matrix[i][0] = 0 #Across first column

for j in range(len(str2)+1):
    matrix[0][j] = 0 # Across first row

#Step 2: Matrix filling according to the scheme

for i in range(1, len(str1)+1):
    for j in range(1, len(str2)+1):
        #Assigning match sore if sequence matches
        if (str1[i-1] == str2[j-1]):
            diagonal_score = match 
        else:
            diagonal_score = mismatch #Assigning mismatch penalty if sequence mismatch

        matrix[i][j] = max(0, matrix[i-1][j-1]+diagonal_score, matrix[i][j-1]+gap, matrix[i-1][j]+gap)

# Matrix filling is done by taking the maximum of these four values

# Step 3: Finding starting point for Traceback step

#Finding the highest score in matrix and cells where they occur
max_value = np.max(matrix)
max_cells = np.argwhere(matrix == max_value)

ti, tj = max_cells[0] # First maximum cell for tracebcak indices

#Step 4: Traceback

#Initialised two empty strings which will be filled using aligned sequences and third string to represent vertical bar
aligned1 = ""
aligned2 = ""
aligned = ""

while ti>0 and tj>0 and matrix[ti][tj] != 0:

    if str1[ti-1] == str2[tj-1]:
        diagonal_score = match #In case of match
    else:
        diagonal_score = mismatch #In case of mismatch

    # Diagonal movement
    if (matrix[ti][tj] == matrix[ti-1][tj-1]+diagonal_score):

        if (str1[ti-1] == str2[tj-1]):
            aligned = "|"+aligned
        else:
            aligned = " "+aligned

        #Adds character to both the sequences
        aligned1 = str1[ti-1]+aligned1
        aligned2 = str2[tj-1]+aligned2

        ti = ti-1
        tj = tj-1

    # Upward movement
    elif matrix[ti][tj] == matrix[ti-1][tj]+gap:

        #Adds character to first sequence and gap to second
        aligned1 = str1[ti-1]+aligned1
        aligned2 = "-"+aligned2
        aligned = " "+aligned

        ti = ti-1

    # Left movement
    elif matrix[ti][tj] == matrix[ti][tj-1]+gap:

        #Adds gap to first and character to second sequence
        aligned1 = "-"+aligned1
        aligned2 = str2[tj-1]+aligned2
        aligned = " "+aligned

        tj = tj-1


#Print the complete locala alignment matrix
print("a. Complete Local Alignment matrix: \n")
print(matrix)

#Print the maximum score
print("\nb. Maximum Score: ", int(max_value))
print("Cells where maximum occur: ", end="")
for i in max_cells:
    print(i, end="")

#Print the optimal local alignment
print("\n\nd. Optimal Local Alignment: ")
print(aligned1)
print(aligned)
print(aligned2)

#Alignment score
print("\ne. Optimal Local Alignment Score: ", int(max_value))