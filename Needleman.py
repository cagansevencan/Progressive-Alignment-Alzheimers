my_matrix = []
#Fill out the 0th row
my_matrix.append([1,3,5,7])
# Fill out the 1st row
my_matrix.append([2,3,4,5])
# Fill out the 2st row
my_matrix.append([5,2,20,3])

#Helper function to print out matrices
def print_matrix(mat):
    #Loop over all rows
    for i in range(0, len(mat)):
        print("[", end ='')
        # Loop over each column in row i
        for j in range(0,len(mat[i])):
            # Print out the value in row i, column j
            print(mat[i][j], end='')
            # Only add a tab if we are not in the last column
            if j != len(mat[i]) -1:
                print("\t", end = '')
        print("]")

print_matrix(my_matrix)


# Use these values to calculate scores
gap_penalty = -1
match_award = 1
mismatch_penalty = -1




with open("APOE Variant 1.txt",  'r') as data1:
    d1 = data1.read()

with open("APOELocal.txt",  'r') as data2:
    localA = data2.read()

#print(d1)

with open("APOE Variant 5.txt",  'r') as data5:
    d5 = data5.read()

#print(d2)

with open("APOE Variant 3.txt",  'r') as data3:
    d3 = data3.read()


with open("ABCA7 - Brown Rat.txt",  'r') as data7:
    bRat = data7.read()

with open("ABCA7 - House Mouse V1.txt",  'r') as data8:
    hM1 = data8.read()

with open("ABCA7 Human.txt",  'r') as data9:
    humanA = data9.read()



# Make a score matrix with these two sequences
seq1 = "ATTACA"
seq2 = "ATGCT"
seq3 = "ATTGAC"

# A function for making a matrix of zeroes
def zeros(rows, cols):
    # Define an empty list
    retval = []
    # Set up the rows of the matrix
    for x in range(rows):
        # For each row, add an empty list
        retval.append([])
        # Set up the columns in each row
        for y in range(cols):
            # Add a zero to each column in each row
            retval[-1].append(0)
    # Return the matrix of zeros
    return retval

# A function for determining the score between any two bases in alignment
def match_score(alpha, beta):
    if alpha == beta:
        return match_award
    elif alpha == '-' or beta == '-':
        return gap_penalty
    else:
        return mismatch_penalty

# The function that actually fills out a matrix of scores
def needleman_wunsch(seq1, seq2):
    # length of two sequences
    n = len(seq1)
    m = len(seq2)

    # Generate matrix of zeros to store scores
    score = zeros(m + 1, n + 1)

    # 1. Fill out first column
    for i in range(0, m+1):
        score[i][0] = gap_penalty * i

    # 2. Fill out first row
    for j in range(0, n+1):
        score[0][j] = gap_penalty * j


    # 3. Fill out all other values in the score matrix
    for i in range(1, m+1):
        for j in range(1, n+1):
        # Calculate the score by checking the top, left, and diagonal cells
            delete = score[i-1][j] + gap_penalty
            insert = score[i][j-1] + gap_penalty
            diag = score[i-1][j-1] + match_score(seq1[j-1], seq2[i-1])
        # Record the maximum score from the three possible scores calculated above
            score[i][j] = max(diag, delete, insert)

    # Back-trace from right of the score matrix and compute the alignment
    # Variables to store alignment
    align1 = ""
    align2 = ""

    # Starting from the bottom right cell in the matrix
    i = m
    j = n

    # Using i and j to keep track of the position in the matrix
    while i > 0 and j > 0:  # making sure we are inside the borders
        currentScore = score[i][j]
        diagonalScore = score[i-1][j-1]
        leftScore = score[i][j-1]
        aboveScore = score[i-1][j]
        # Checking to see which cell the current score was coming from,
        # next update i and j accordingly to that cell.
        if currentScore == diagonalScore + match_score(seq1[j-1], seq2[i-1]):
            align1 += seq1[j-1]
            align2 += seq2[i-1]
            i -=1
            j -= 1
        elif currentScore == leftScore + gap_penalty:
            align1 += seq1[j-1]
            align2 += '-'
            j -= 1
        elif currentScore == aboveScore + gap_penalty:
            align1 += '-'
            align2 += seq2[i-1]
            i -= 1
    # Finish tracing up until top left
    while j > 0:
        align1 += seq1[j-1]
        align2 += '-'
        j -= 1
    while i > 0:
        align1 += '-'
        align2 += seq2[i-1]
        i -= 1
    # Since we started backwards we need to reverse our sequences
    align1 = align1[::-1]
    align2 = align2[::-1]

    sym = ''
    ident = 0
    seq_score =0
    seq_size = len(align1)
    for i in range(seq_size):
        a1 = align1[i]
        a2 = align2[i]
        if a1 == a2:
            sym += a1
            ident +=1
            seq_score += match_score(a1, a2)

        else:
            seq_score += match_score(a1,a2)
            sym += ' '

    ident = ident / seq_size * 100

    print('Identity = %2.1f percent' % ident)
    print('Score = %d\n'% seq_score)
    print(align1)
    print(sym)
    print(align2)

    #print_matrix(score)

    return(align1, align2)



def positions_of_diff(aligned_str):
    diff_position = []
    str_length = len(aligned_str)
    for x in range(str_length):
        if aligned_str[x] == '-':
            diff_position.append(x+1)

    return diff_position



    # Test out the needleman_wunsch() function
print_matrix(needleman_wunsch(bRat, humanA))