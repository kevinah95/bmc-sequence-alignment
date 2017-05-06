import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from Bio.Seq import Seq
#------------------------------
seq1 = Seq("GATCCA")
seq2 = Seq("GTGCCT")

match_award = 1
mismatch_penalty = -1
gap_penalty = -2
arrows = None


def match_score(alpha, beta):
    if alpha == beta:
        return match_award
    else:
        return mismatch_penalty


def needleman_wunsch(seq1, seq2):
    m, n = len(seq2), len(seq1)  # length of two sequences
    #---------
    dt = np.dtype([('diagonal', np.str, 1),
                   ('up', np.str, 1), ('left', np.str, 1)])
    """
    Pointer Matrix size: (m + 1) x (n + 1)
    Datatype: 3 tuple to save direction of the path
    
    For more info: https://docs.scipy.org/doc/numpy/user/basics.rec.html
    """
    pt_mat = np.zeros((m + 1, n + 1), dtype=dt)

    #---------
    score_matrix = np.zeros((m + 1, n + 1), dtype=int)      # the DP table
    # Calculate DynamicProgramming table
    for i in range(0, m + 1):
        score_matrix[i][0] = gap_penalty * i
        pt_mat[i][0]['up'] = 'U'
    for j in range(0, n + 1):
        score_matrix[0][j] = gap_penalty * j
        pt_mat[0][j]['left'] = 'L'
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diagonal = score_matrix[i - 1][j - 1] + \
                match_score(seq2[i - 1], seq1[j - 1])
            up = score_matrix[i - 1][j] + gap_penalty
            left = score_matrix[i][j - 1] + gap_penalty
            max_pointer = max(diagonal, up, left)
            score_matrix[i][j] = max_pointer
            if (diagonal == max_pointer):
                pt_mat[i][j]['diagonal'] = 'D'
            if (up == max_pointer):
                pt_mat[i][j]['up'] = 'U'
            if (left == max_pointer):
                pt_mat[i][j]['left'] = 'L'

    # print(score_matrix)
    # print(pt_mat)

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
        #print("current ",i,j,end='')
        a[:, :2] = j, i
        if score_current == score_diagonal + match_score(seq2[i - 1], seq1[j - 1]):
            align1 += seq2[i - 1]
            align2 += seq1[j - 1]
            i -= 1
            j -= 1
        elif score_current == score_left + gap_penalty:
            align1 += seq2[i - 1]
            align2 += '-'
            i -= 1
        elif score_current == score_up + gap_penalty:
            align1 += '-'
            align2 += seq1[j - 1]
            j -= 1
        a[:, 2:] = j, i
        if first:
            # First loop copy a
            arrows = np.copy(a)
            first = False
        else:
            # Concatenate origin-target
            arrows = np.concatenate((arrows, a), axis=0)
    a[:, :2] = a[:, 2:]

    # Finish tracing up to the top left cell
    while i > 0:
        align1 += seq2[i - 1]
        align2 += '-'
        i -= 1
        a[:, 2:] = j, i
        arrows = np.concatenate((arrows, a), axis=0)
    while j > 0:
        align1 += '-'
        align2 += seq1[j - 1]
        j -= 1
        a[:, 2:] = j, i
        arrows = np.concatenate((arrows, a), axis=0)
    return score_matrix, pt_mat


#-------------------------------
plt.rcParams["figure.figsize"] = 4, 5
param = {"grid.linewidth": 1.6,
         "grid.color": "lightgray",
         "axes.linewidth": 1.6,
         "axes.edgecolor": "lightgray"}
plt.rcParams.update(param)

# Data
headh = seq1
headv = seq2

v, pt_mat = needleman_wunsch(seq1, seq2)

# Plot
fig, ax = plt.subplots()
ax.set_xlim(-1.5, v.shape[1] - .5)
ax.set_ylim(-1.5, v.shape[0] - .5)
ax.invert_yaxis()
for i in range(v.shape[0]):
    for j in range(v.shape[1]):
        ax.text(j, i, v[i, j], ha="center", va="center")
for i, l in enumerate(headh):
    ax.text(i + 1, -1, l, ha="center", va="center", fontweight="semibold")
for i, l in enumerate(headv):
    ax.text(-1, i + 1, l, ha="center", va="center", fontweight="semibold")

ax.xaxis.set_minor_locator(ticker.FixedLocator(
    np.arange(-1.5, v.shape[1] - .5, 1)))
ax.yaxis.set_minor_locator(ticker.FixedLocator(
    np.arange(-1.5, v.shape[1] - .5, 1)))
plt.tick_params(axis='both', which='both', bottom='off', top='off',
                left="off", right="off", labelbottom='off', labelleft='off')
ax.grid(True, which='minor')


arrowprops = dict(facecolor='blue', alpha=0.5, lw=0,
                  shrink=0.2, width=2, headwidth=7, headlength=7)

# all paths
for i in range(1, pt_mat.shape[0]):
    for j in range(1, pt_mat.shape[1]):
        if(pt_mat[i][j]['left'] != ''):
            ax.annotate("", xy=(j - 1, i),
                        xytext=(j, i), arrowprops=arrowprops)
        if(pt_mat[i][j]['diagonal'] != ''):
            ax.annotate("", xy=(j - 1, i - 1),
                        xytext=(j, i), arrowprops=arrowprops)
        if(pt_mat[i][j]['up'] != ''):
            ax.annotate("", xy=(j, i - 1),
                        xytext=(j, i), arrowprops=arrowprops)

# optimal path
arrowprops.update(facecolor='crimson')
for i in range(arrows.shape[0]):
    ax.annotate("", xy=arrows[i, 2:],
                xytext=arrows[i, :2], arrowprops=arrowprops)
plt.show()
