import numpy as np
from Bio.Seq import Seq
#------------------------------
penalty = {'MATCH': 1, 'MISMATCH': -1, 'GAP': -2}
arrows = np.array([[0, 0, 0, 0]])


def match_score(alpha, beta):
    if alpha == beta:
        return penalty['MATCH']
    else:
        return penalty['MISMATCH']


def finalize(align1, align2, score):
    align1 = align1[::-1]  # reverse sequence 1
    align2 = align2[::-1]  # reverse sequence 2

    i, j = 0, 0

    # calcuate identity, score and aligned sequeces
    symbol = ''
    found = 0
    for i in range(0, len(align1)):
        # if two AAs are the same, then output the letter
        if align1[i] == align2[i]:
            symbol += '|'

        # if they are not identical and none of them is gap
        elif align1[i] != align2[i] and align1[i] != '-' and align2[i] != '-':
            symbol += ' '
            found = 0

        # if one of them is a gap, output a space
        elif align1[i] == '-' or align2[i] == '-':
            symbol += ' '

    print ('Score =', score)
    print (align1)
    print (symbol)
    print (align2)


def insideBand(i, j, k):
    return (-k <= i - j <= k)


def kband(seq_alpha_col, seq_beta_row, k_value, p_penalty, score_only):
    if not seq_alpha_col or not seq_beta_row:
        print("Alguna de las secuencias está vacía.")
        return
    global penalty
    penalty = p_penalty
    m, n = len(seq_beta_row), len(seq_alpha_col)  # length of two sequences
    if m != n:
        print("LAs secuencias no tienen el mismo largo.")
        return
    #---------
    pointer_mat = []  # pointer matrix

    if not score_only:
        # Custom datatype 3-tuple for ('D','U','L')
        dt = np.dtype([('diagonal', np.str, 1),
                       ('up', np.str, 1), ('left', np.str, 1)])

        pointer_mat = np.zeros((m + 1, n + 1), dtype=dt)

    #---------
    score_matrix = np.zeros((m + 1, n + 1), dtype=int)      # the DP table
    # Calculate DP table
    for i in range(0, k_value):
        score_matrix[i][0] = penalty['GAP'] * i
        if not score_only:
            pointer_mat[i][0]['up'] = 'U'
    for j in range(0, k_value):
        score_matrix[0][j] = penalty['GAP'] * j
        if not score_only:
            pointer_mat[0][j]['left'] = 'L'
    for i in range(1, m + 1):
        for h in range(-k_value, k_value):
            j = i + h
            if 1 <= j <= n:
                diagonal = score_matrix[i - 1][j - 1] + \
                    match_score(seq_beta_row[i - 1], seq_alpha_col[j - 1])
                up = score_matrix[i - 1][j] + penalty['GAP']
                left = score_matrix[i][j - 1] + penalty['GAP']
                score_matrix[i][j] = diagonal
                max_pointer = diagonal
                if insideBand(i - 1, j, k_value):
                    max_pointer = max(diagonal,up)
                    score_matrix[i][j] = max_pointer
                if insideBand(i, j - 1, k_value):
                    max_pointer = max(diagonal,left)
                    score_matrix[i][j] = max_pointer
                
                if not score_only:
                    if (diagonal == max_pointer):
                        pointer_mat[i][j]['diagonal'] = 'D'
                    if (up == max_pointer):
                        pointer_mat[i][j]['up'] = 'U'
                    if (left == max_pointer):
                        pointer_mat[i][j]['left'] = 'L'
    if score_only:
        print(score_matrix)
    # print(pt_mat)

    # Traceback and compute the alignment
    align1, align2 = '', ''
    i, j = m, n  # start from the bottom right cell
    a = np.array([[0, 0, 0, 0]])

    first = True
    global arrows
    score = score_matrix[i][j]
    while i > 0 and j > 0:  # end toching the top or the left edge
        score_current = score_matrix[i][j]
        score_diagonal = score_matrix[i - 1][j - 1]
        score_up = score_matrix[i][j - 1]
        score_left = score_matrix[i - 1][j]
        if not score_only:
            a[:, :2] = j, i
        if score_current == score_diagonal + match_score(seq_beta_row[i - 1], seq_alpha_col[j - 1]):
            align1 += seq_beta_row[i - 1]
            align2 += seq_alpha_col[j - 1]
            i -= 1
            j -= 1
            score = score + score_diagonal
        elif score_current == score_left + penalty['GAP']:
            align1 += seq_beta_row[i - 1]
            align2 += '-'
            i -= 1
            score = score + score_left
        elif score_current == score_up + penalty['GAP']:
            align1 += '-'
            align2 += seq_alpha_col[j - 1]
            j -= 1
            score = score + score_up
        if not score_only:
            a[:, 2:] = j, i
            if first:
                arrows = np.copy(a)
                first = False
            else:
                arrows = np.concatenate((arrows, a), axis=0)
    if not score_only:
        a[:, :2] = a[:, 2:]
    # print(arrows)
    # print(a)

    # Finish tracing up to the top left cell
    while i > 0:
        align1 += seq_beta_row[i - 1]
        align2 += '-'
        i -= 1
        if not score_only:
            a[:, 2:] = j, i
            arrows = np.concatenate((arrows, a), axis=0)
    while j > 0:
        align1 += '-'
        align2 += seq_alpha_col[j - 1]
        j -= 1
        if not score_only:
            a[:, 2:] = j, i
            arrows = np.concatenate((arrows, a), axis=0)

    finalize(align1, align2, score)
    return score_matrix, pointer_mat, arrows


if __name__ == '__main__':
    alpha = Seq("MADRE")
    beta = Seq("PADRE")
    kband(alpha, beta,2, penalty, score_only=True)
