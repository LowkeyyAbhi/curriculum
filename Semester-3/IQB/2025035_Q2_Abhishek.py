import numpy as np

str1 = "ATGAT"  #DNA sequence 1
str2 = "ATT" #DNA sequence 2 

# str1 = "GCTTAGC"
# str2 = "GCATTGC"

match = 2 #match score as provided in assignment pdf
gap = -2 #gap penalty score as provided in assignment pdf
mismatch = -1 #mismatch penalty as provided in assignment pdf

matrix = np.zeros((len(str1)+1, len(str2)+1))  #Create a matrix through numpy library in python

# print(matrix)

#Step 1: Initialisation the matrix by assigning gap penalty as part of the
# initialisaton process

for i in range(len(str1)+1): 
    matrix[i][0] = i*gap  # Across first column

for j in range(len(str2)+1):
    matrix[0][j] = j*gap # Across first row

#Step 2: Matrix Filling according to scoring scheme
for i in range(1, len(str1)+1):
    for j in range(1, len(str2)+1):
        if str1[i-1] == str2[j-1]: # Assigning match score if two sequence match
            diagonal_score = match
        else:
            diagonal_score = mismatch #Assigning mismatch penalty if sequences didn't match

        matrix[i][j] = max(matrix[i-1][j-1]+diagonal_score, matrix[i][j-1]+gap, matrix[i-1][j]+gap)
#Matrix is filled by taking the maximum of the three

#Step 3: Traceback

#Initialised two empty strings which will be filled using aligned sequences and third string to represnt vertical bar
aligned1 = ""
aligned_2 = ""
aligned = ""

#Since traceback step begins from bottom-most we assign ti and tj as 
ti = len(str1) #Length of sequence 1
tj = len(str2) #Length of sequence 2

while ti>0 or tj>0: # Or of two condition so that both sequence reaches their end
    if ti>0 and tj>0:

        #Checking match or mismatch
        if str1[ti-1] == str2[tj-1]:
            diagonal_score = match
        else:
            diagonal_score = mismatch

    #Diagonal movement
    if (ti>0 and tj>0 and matrix[ti][tj] == matrix[ti-1][tj-1]+diagonal_score):

        #Adds character to both the sequences
        aligned1 = str1[ti-1]+aligned1
        aligned_2 = str2[tj-1]+aligned_2

        if (str1[ti-1] == str2[tj-1]):
            aligned = "|"+aligned
        else:
            aligned = " "+aligned

        ti = ti-1
        tj = tj-1


    #Upward movement 
    elif (ti>0 and matrix[ti][tj] == matrix[ti-1][tj]+gap):

        #Adds character to first sequence and gap to second
        aligned1 = str1[ti-1]+aligned1
        aligned_2 = "-"+aligned_2
        aligned = " "+aligned

        ti = ti-1

    #Left movement
    elif (tj>0 and matrix[ti][tj] == matrix[ti][tj-1]+gap):

        #Adds gap to first and character to second
        aligned1 = "-"+aligned1
        aligned_2 = str2[tj-1]+aligned_2
        aligned = " "+aligned

        tj = tj-1

#Printing the complete global alignment matrix
print("b. Complete Global Alignment matrix: \n")
print(matrix)

#Printing the optimal alignment

print("\nc. Optimal Global Alignment: ")
print(aligned1)
print(aligned)
print(aligned_2)

#Alignment Score
print("\nc. Optimal Alignment Score: ", int(matrix[len(str1)][len(str2)])) #Last entry in the aligned matrix
