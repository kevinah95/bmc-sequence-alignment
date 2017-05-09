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


def needleman_wunsch(seq_alpha_col, seq_beta_row, score_only):
    if not seq_alpha_col or not seq_beta_row:
        print("Alguna de las secuencias está vacía.")
        return
    m, n = len(seq_beta_row), len(seq_alpha_col)  # length of two sequences
    #---------
    pointer_mat = [] # pointer matrix
    
    if not score_only:
        # Custom datatype 3-tuple for ('D','U','L') 
        dt = np.dtype([('diagonal', np.str, 1),
                   ('up', np.str, 1), ('left', np.str, 1)])
        
        pointer_mat = np.zeros((m + 1, n + 1), dtype=dt)
        
    #---------
    score_matrix = np.zeros((m + 1, n + 1), dtype=int)      # the DP table
    # Calculate DP table
    for i in range(0, m + 1):
        if not score_only:
            pointer_mat[i][0]['up'] = 'U'
    for j in range(0, n + 1):
        if not score_only:
            pointer_mat[0][j]['left'] = 'L'
    max_score = 0
    max_pos   = None    # The row and columbn of the highest score in matrix.
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diagonal = score_matrix[i - 1][j - 1] + \
                match_score(seq_beta_row[i - 1], seq_alpha_col[j - 1])
            up = score_matrix[i - 1][j] + penalty['GAP']
            left = score_matrix[i][j - 1] + penalty['GAP']
            max_pointer = max(0,diagonal, up, left)
            if max_pointer > max_score:
                max_score = max_pointer
                max_pos   = (i, j)
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
    #print(pt_mat)

    # Traceback and compute the alignment
    align1, align2 = '', ''
    i, j = m, n  # start from the bottom right cell
    a = np.array([[0, 0, 0, 0]])

    first = True
    global arrows
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
        elif score_current == score_left + penalty['GAP']:
            align1 += seq_beta_row[i - 1]
            align2 += '-'
            i -= 1
        elif score_current == score_up + penalty['GAP']:
            align1 += '-'
            align2 += seq_alpha_col[j - 1]
            j -= 1
        if not score_only:
            a[:, 2:] = j, i
            if first:
                arrows = np.copy(a)
                first = False
            else:
                arrows = np.concatenate((arrows, a), axis=0)
    if not score_only:
        a[:, :2] = a[:, 2:]
    #print(arrows)
    #print(a)

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
    
    
    return score_matrix, pointer_mat, arrows
if __name__ == '__main__':
    alpha = Seq("SEND")
    beta = Seq("AND")
    needleman_wunsch(alpha,beta,score_only=True)