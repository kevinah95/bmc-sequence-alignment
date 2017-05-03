import numpy
seq1="SEND"
seq2="AND"
#TODO cambiar
seq1 = list(seq1)
seq2 = list(seq2)

match_award      = 1
mismatch_penalty = -1
gap_penalty      = -2

def match_score(alpha, beta):
    if alpha == beta:
        return match_award
    else:
        return mismatch_penalty


def needleman_wunsch(seq1, seq2):
    m, n = len(seq2), len(seq1)  # length of two sequences
    score_matrix = numpy.zeros((m+1, n+1))      # the DP table
    # Calculate DP table
    for i in range(0, m + 1):
        score_matrix[i][0] = gap_penalty * i
    for j in range(0, n + 1):
        score_matrix[0][j] = gap_penalty * j
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diagonal = score_matrix[i - 1][j - 1] + match_score(seq2[i-1], seq1[j-1])
            up = score_matrix[i - 1][j] + gap_penalty
            left = score_matrix[i][j - 1] + gap_penalty
            score_matrix[i][j] = max(diagonal,up,left)
    print(score_matrix)
    # Traceback and compute the alignment 
    align1, align2 = '', ''
    i,j = m,n # start from the bottom right cell
    while i > 0 and j > 0: # end toching the top or the left edge
        score_current = score_matrix[i][j]
        score_diagonal = score_matrix[i-1][j-1]
        score_up = score_matrix[i][j-1]
        score_left = score_matrix[i-1][j]

        if score_current == score_diagonal + match_score(seq2[i-1], seq1[j-1]):
            align1 += seq2[i-1]
            align2 += seq1[j-1]
            i -= 1
            j -= 1
        elif score_current == score_left + gap_penalty:
            align1 += seq2[i-1]
            align2 += '-'
            i -= 1
        elif score_current == score_up + gap_penalty:
            align1 += '-'
            align2 += seq1[j-1]
            j -= 1

    # Finish tracing up to the top left cell
    while i > 0:
        align1 += seq2[i-1]
        align2 += '-'
        i -= 1
    while j > 0:
        align1 += '-'
        align2 += seq1[j-1]
        j -= 1
    print(score_matrix)
    print(align1,align2)

needleman_wunsch(seq1,seq2)
