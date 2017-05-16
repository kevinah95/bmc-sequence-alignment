import numpy as np
from Bio.Seq import Seq
#------------------------------
penalty = {'MATCH': 3, 'MISMATCH': -3, 'GAP': -2}
arrows = np.array([[0, 0, 0, 0]])


def match_score(alpha, beta):
    if alpha == beta:
        return penalty['MATCH']
    else:
        return penalty['MISMATCH']


#traceback
def finalize(align1, align2,score):
    align1 = align1[::-1]    #reverse sequence 1
    align2 = align2[::-1]    #reverse sequence 2
    
    
    #calcuate identity, score and aligned sequeces
    sym = ''
    for i in range(len(align1)):
        a1 = align1[i]
        a2 = align2[i]
        if a1 == a2:                
            sym += '|'
        elif a1 != a2 and a1 != '-' and a2 != '-': 
            sym += ' '
        elif a1 == '-' or a2 == '-':          
            sym += ' '
    print('Score =', score)
    print(align1)
    print(sym)
    print(align2)


def next_move(score_matrix, i, j):
    diag = score_matrix[i - 1][j - 1]
    up   = score_matrix[i - 1][j]
    left = score_matrix[i][j - 1]
    if diag >= up and diag >= left:     # Tie goes to the DIAG move.
        return 1 if diag != 0 else 0    # 1 signals a DIAG move. 0 signals the end.
    elif up > diag and up >= left:      # Tie goes to UP move.
        return 2 if up != 0 else 0      # UP move or end.
    elif left > diag and left > up:
        return 3 if left != 0 else 0    # LEFT move or end.
    else:
        # Execution should not reach here.
        raise ValueError('invalid move during traceback')


def smith_waterman(seq_alpha_col, seq_beta_row,p_penalty, score_only):
    """
    Llena las matrices de acuerdo al algoritmo Smith-Waterman.

    Esta función realiza las operaciones necesarias para lograr el 
    alineamiento local (Smith-Waterman) entre dos secuencias. 
    Un alineamiento local encuentra la mejor concordancia entre 
    los caracteres de dos subsecuencias de las secuencias principales.

    Al realizar alineamientos, se puede especificar el match score,
    el mismath score y el gap penalty. El match score indica la compatibilidad 
    entre un alineamiento entre dos caracteres en las secuencias.
    Los caracteres altamente compatibles deben recibir la puntuación del match, 
    y los que no sean compatibles deben recibir la puntuación de mismatch. 
    Los gaps deben ser negativos.

    Parameters
    ----------
    seq_alpha_col : Bio.Seq.Seq
        Primer secuencia a ser comparada
    
    seq_beta_row : Bio.Seq.Seq
        Segunda secuencia a ser comparada
    
    p_penalty : dict(str -> int)
        Diccionario que contiene los valores de MATCH, MISMATCH y GAP

    score_only : boolean
        Cuando es True solamente muestra el alineamiento y el valor del score, 
        así se utiliza menos memoria y es más rápido.
        Cuando es False guarda la matriz y la traza que recorre el alineamiento
        en la carpeta `/output`.

    Returns
    -------
    np.array, np.array, np.array
        Tres arreglos que contienen los puntajes, los punteros y la matriz del mejor alineamiento

    """
    if not seq_alpha_col or not seq_beta_row:
        print("Alguna de las secuencias está vacía.")
        return
    
    m, n = len(seq_beta_row), len(seq_alpha_col)  # length of two sequences
    global penalty
    penalty = p_penalty
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
            if not (max_pointer == 0):
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
    a = np.array([[0, 0, 0, 0]])

    first = True
    global arrows
    END, DIAG, UP, LEFT = range(4)
    i,j        = max_pos
    print(max_pos)
    score=score_matrix[i][j]
    move         = next_move(score_matrix, i, j)
    while move != END:  # end toching the top or the left edge
        #-------------
        score_diagonal = score_matrix[i - 1][j - 1]
        score_up = score_matrix[i][j - 1]
        score_left = score_matrix[i - 1][j]
        if not score_only:
            a[:, :2] = j, i
        if move == DIAG:
            align1 += seq_beta_row[i - 1]
            align2 += seq_alpha_col[j - 1]
            i -= 1
            j -= 1
            score=score+penalty['MATCH']
        elif move == UP:
            align1 += seq_beta_row[i - 1]
            align2 += '-'
            i -= 1
            score=score+penalty['GAP']
        else:
            align1 += '-'
            align2 += seq_alpha_col[j - 1]
            j -= 1
            score=score+penalty['GAP']

        if not score_only:
            a[:, 2:] = j, i
            if first:
                arrows = np.copy(a)
                first = False
            else:
                arrows = np.concatenate((arrows, a), axis=0)
        move = next_move(score_matrix, i, j)

    
    finalize(align1,align2,score)
    return score_matrix, pointer_mat, arrows
if __name__ == '__main__':
    alpha = Seq("TGTTACGG")
    beta = Seq("GGTTGACTA")
    penalty = {'MATCH': 3, 'MISMATCH': -3, 'GAP': -2}
    smith_waterman(alpha,beta,penalty,score_only=True)