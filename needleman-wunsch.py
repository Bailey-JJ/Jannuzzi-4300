'''
Module: needleman-wunsch.py
Author: Bailey Jannuzzi
Project: Needleman-Wunsch Algorithm for Global Alignment
Description:
    Simple implementation of Needleman-Wunsch Algorithm. 
    Used for global alignment of two sequences.
    Declares inputs, outputs, and scoring variables. 
    Script initializes matrix, fills matrix using scoring values, defines a traceback and scoring function, then runs traceback procedure to determine the best alignment. 
    Output from terminal is the filled matrix and final alignment with score attached.
   
    There can be multiple optimal alignments, this script will only produce one possible best alignment option.

    Functions are defined first, with type hints, then main() is written and run.
'''

# Printing matrix to check it
def print_matrix(matrix: list[list[int]], seq_01: str = '', seq_02: str = '') -> None:
    '''
    Function takes matrix values and formats into rows and columns.
    Sets top header and left labels as well.
    '''
    print('        ', end = '')
    for letter in '-' + seq_01:
        print(f'{letter:4}', end = '  ')
    print()

    print('       ----------------------', end = '\n')
    for i, row in enumerate(matrix):
        if i == 0:
            letter = '-'
        else:
            letter = seq_02[i-1]

        print(f'{letter:4}', end = '|')

        for _ in row:
            print(f'{_:4}', end = '  ')
        print()
    print()



# Initializaing the matrix
def init_matrix(seq_01: str = '', seq_02: str = '', gap: int = -2) -> list[list[int]]:
    '''
    For this I used the class example sequences, so I just kept sequence 2 as the row sequence. 
    But, length of sequence doesn't affect which spot it can go in.
    '''
    rows: int = len(seq_02) + 1
    cols: int = len(seq_01) + 1

    matrix: list[list[int]] = [[0 for j in range(cols)] for i in range(rows)]

    for i in range(1, rows):
        matrix[i][0] = i * gap
    for j in range(1, cols):
        matrix[0][j] = j * gap

    return matrix



# Filling matrix
def fill_matrix(matrix: list[list[int]], seq_01: str = '', seq_02: str = '', match: int = 1, mismatch: int = -1, gap: int = -2) -> None:
    '''
    Calculates each cell's score.
    '''
    rows: int = len(seq_02) + 1
    cols: int = len(seq_01) + 1

    for i in range(1, rows):
        for j in range(1, cols):
            # If the letters match, add 1 to the score, then move diagonally
            if seq_02[i-1] == seq_01[j-1]:
                diag = matrix[i-1][j-1] + match
            # If the letters dont match, subtract 1 from the score, then move diagonally
            else:
                diag = matrix[i-1][j-1] + mismatch
            
            # Calculate up and left moves
            up = matrix[i-1][j] + gap
            left = matrix[i][j-1] + gap

            # Pick the most positive score to put into the cell
            matrix[i][j] = max(diag, up, left)



# Defining traceback and scoring function
def trace(matrix: list[list[int]], seq_01: str = '', seq_02: str = '', match: int = 1, mismatch: int = -1, gap: int = -2) -> str:
    '''
    Details traceback procedure and has a scoring function so that end alignment output can include score.
    [add more thoughts later]
    '''
    # Function to score an alignment
    def score_alignment(aligned_seq_01: str = '', aligned_seq_02: str = '') -> int:
        '''
        Takes the two completed alignments and recalculates their score based on what was needed to get to that score.
        '''
        total: int = 0
        for x, y in zip(aligned_seq_01, aligned_seq_02):
            if x == y and x != '-':
                total += match
            elif x == '-' or y == '-':
                total += gap
            else:
                total += mismatch
        return total

    # Declaring whatever variables are needed for this section
    i: int = (len(seq_02) + 1) - 1
    j: int = (len(seq_01) + 1) - 1
    aligned_seq_01: str = ''
    aligned_seq_02: str = ''

    while i > 0 or j > 0:
        scoring: int = matrix[i][j]

        if i > 0 and j > 0:
            if seq_02[i-1] == seq_01[j-1]:
                # Basically, if the letters match, check if cell's score is from a match. If they match, add the letter to the alignments.
                if scoring == matrix[i-1][j-1] + match:
                    aligned_seq_01 = seq_01[j-1] + aligned_seq_01
                    aligned_seq_02 = seq_02[i-1] + aligned_seq_02
                    # Changing the i and j values, moves to check the next diagonal value
                    i -= 1
                    j -= 1
            
            # Otherwise, check if the letters don't match and if cell's score is from a mismatch. If its a mismatch, add the two different letters to the alignments.
            elif scoring == matrix[i-1][j-1] + mismatch:
                aligned_seq_01 = seq_01[j-1] + aligned_seq_01
                aligned_seq_02 = seq_02[i-1] + aligned_seq_02
                i -= 1
                j -= 1

        # Checking if the cell's core came from an up move
        if i > 0 and scoring == matrix[i-1][j] + gap:
            aligned_seq_01 = '-' + aligned_seq_01
            aligned_seq_02 = seq_02[i-1] + aligned_seq_02
            # I just want to update i so it moves up in the matrix
            i -= 1

        # Checking if the cell's score came from a left move
        if j > 0 and scoring == matrix[i][j-1] + gap:
            aligned_seq_01 = seq_01[j-1] + aligned_seq_01
            aligned_seq_02 = '-' + aligned_seq_02
            # Updates only j to move left in the matrix
            j -= 1

    return f'Best Alignment:\n{aligned_seq_01}\n{aligned_seq_02}\nAlignment Score: {score_alignment(aligned_seq_01, aligned_seq_02)}'


# Defining main function
def main(seq_01: str = 'ATG', seq_02: str = 'GGAATGG') -> None:
    '''
    Example sequences from class slides are set as the default values for the two sequences.
    '''

    # Initializing Matrix
    matrix: list[list[int]] = init_matrix(seq_01, seq_02)
    print_matrix(matrix, seq_01, seq_02)


    # Filling Matrix
    fill_matrix(matrix, seq_01, seq_02)


    # Printing matrix to check results
    print_matrix(matrix, seq_01, seq_02)


    # Best Alignment and Alignment Score
    print(trace(matrix, seq_01, seq_02))


if __name__ == "__main__":
    main()
