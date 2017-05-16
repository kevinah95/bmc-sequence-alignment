import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import datetime
from Bio.Seq import Seq
if __name__ == '__main__':
    from glocal import glocal
else:
    from .glocal import glocal

#-------------------------------
def plot_gl(seq_alpha_col,seq_beta_row,p_penalty):

    if not seq_alpha_col or not seq_beta_row:
        print("Alguna de las secuencias está vacía.")
        return

    
    plt.rcParams["figure.figsize"] = 18, 18
    param = {"grid.linewidth": 1.6,
            "grid.color": "lightgray",
            "axes.linewidth": 1.6,
            "axes.edgecolor": "lightgray",
            "font.size": 8}
    plt.rcParams.update(param)

    # Data
    headh = seq_alpha_col
    headv = seq_beta_row

    score_matrix, pt_mat, arrows = glocal(seq_alpha_col,seq_beta_row,p_penalty,score_only=False)

    # Plot
    fig, ax = plt.subplots()
    ax.set_xlim(-1.5, score_matrix.shape[1] - .5)
    ax.set_ylim(-1.5, score_matrix.shape[0] - .5)
    ax.invert_yaxis()
    for i in range(score_matrix.shape[0]):
        for j in range(score_matrix.shape[1]):
            ax.text(j, i, score_matrix[i, j], ha="center", va="center")
    for i, l in enumerate(headh):
        ax.text(i + 1, -1, l, ha="center", va="center", fontweight="semibold")
    for i, l in enumerate(headv):
        ax.text(-1, i + 1, l, ha="center", va="center", fontweight="semibold")

    ax.xaxis.set_minor_locator(ticker.FixedLocator(
        np.arange(-1.5, score_matrix.shape[1] - .5, 1)))
    ax.yaxis.set_minor_locator(ticker.FixedLocator(
        np.arange(-1.5, score_matrix.shape[1] - .5, 1)))
    plt.tick_params(axis='both', which='both', bottom='off', top='off',
                    left="off", right="off", labelbottom='off', labelleft='off')
                

    #-----------ax.set_aspect('auto')
    ax.grid(True, which='minor')



    arrowprops = dict(facecolor='blue', alpha=0.5, lw=0,
                    shrink=0.2, width=2, headwidth=7, headlength=7)

    # all path
    for i in range(1,pt_mat.shape[0]):
        for j in range(1,pt_mat.shape[1]):
            if(pt_mat[i][j]['left'] != ''):
                ax.annotate("", xy=(j-1,i),
                            xytext=(j,i), arrowprops=arrowprops)
            if(pt_mat[i][j]['diagonal'] != ''):
                ax.annotate("", xy=(j-1,i-1),
                            xytext=(j,i), arrowprops=arrowprops)
            if(pt_mat[i][j]['up'] != ''):
                ax.annotate("", xy=(j,i-1),
                            xytext=(j,i), arrowprops=arrowprops)

    # optimal path
    arrowprops.update(facecolor='crimson')
    for i in range(arrows.shape[0]):
        ax.annotate("", xy=arrows[i, 2:],  # origin
                    xytext=arrows[i, :2], arrowprops=arrowprops)
    #------------
    
    plt.gca().set_aspect('auto')
    time = '{:%Y-%m-%d_%H-%M-%S}'.format(datetime.datetime.now())
    plt.savefig("output/glocal/output-gl_"+time+".pdf", dpi=600)
    #plt.show()

if __name__ == '__main__':
    alpha = Seq("ACTCA")
    beta = Seq("TTCAT")
    penalty = {'MATCH': 1, 'MISMATCH': -1, 'GAP': -2}
    plot_gl(alpha,beta,penalty)