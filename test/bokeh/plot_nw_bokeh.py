from bokeh.plotting import figure, show, output_file
from bokeh.models import Range1d, LabelSet, Label, FixedTicker
from bokeh.models import TickFormatter, Arrow, NormalHead, VeeHead, Text
from bokeh.core.properties import Dict, Int, String
import numpy as np
from Bio.Seq import Seq
if __name__ == '__main__':
    from needleman_wunsch import needleman_wunsch
else:
    from .needleman_wunsch import needleman_wunsch


def plot_nw(seq_alpha_col, seq_beta_row, penalty):
    if not seq_alpha_col or not seq_beta_row:
        print("Alguna de las secuencias está vacía.")
        return

    text_props = {
        "angle": 0,
        "text_color": "black",
        "text_align": "center",
        "text_baseline": "middle"
    }

    m, n = len(seq_alpha_col), len(seq_beta_row)

    p = figure(title='Dist. of 10th Grade Students at Lee High',
               x_axis_location="above")
    p.y_range = Range1d(n + 1, -1.5)

    #-----
    score_matrix, pt_mat, arrows = needleman_wunsch(
        seq_alpha_col, seq_beta_row, penalty, score_only=False)

    for i in range(score_matrix.shape[0]):
        for j in range(score_matrix.shape[1]):
            p.rect(j, i, 0.9, 0.9,
                   fill_alpha=0.6, color="#CAB2D6")
            glyph = Text(x=j, y=i, text=[
                         str(score_matrix[i, j])], text_font_size="9pt", **text_props)
            p.add_glyph(glyph)

    '''for i in range(m + 1):
        for j in range(n + 1):
            p.rect(i, j, 0.9, 0.9,
                   fill_alpha=0.6, color="#CAB2D6")
            p.add_layout(Arrow(end=VeeHead(line_width=4, size=4, line_alpha=.5),
                               x_start=i + 0.3, y_start=j + 0.3, x_end=i + 1 - 0.3, y_end=j + 1 - 0.3, line_width=4, line_cap="round", line_alpha=.5))'''
    # headv
    for i, l in enumerate(seq_beta_row):
        glyph = Text(x=-1, y=i + 1, text=[l],
                     text_font_size="9pt", **text_props)
        p.add_glyph(glyph)

    # headh
    for i, l in enumerate(seq_alpha_col):
        glyph = Text(x=i + 1, y=-1, text=[l],
                     text_font_size="9pt", **text_props)
        p.add_glyph(glyph)

    # all path
    '''padding = 0.3
    for i in range(1, pt_mat.shape[0]):
        for j in range(1, pt_mat.shape[1]):
            if(pt_mat[i][j]['left'] != ''):
                p.add_layout(Arrow(end=VeeHead(line_width=4, size=4, line_alpha=.5),
                                   x_start=(j - 1) + padding, y_start=i, x_end=j - padding, y_end=i, line_width=4, line_cap="round", line_alpha=.5))
            if(pt_mat[i][j]['diagonal'] != ''):
                p.add_layout(Arrow(end=VeeHead(line_width=4, size=4, line_alpha=.5),
                                   x_start=(j - 1) + padding, y_start=(i - 1) + padding, x_end=j - padding, y_end=i - padding, line_width=4, line_cap="round", line_alpha=.5))
            if(pt_mat[i][j]['up'] != ''):
                p.add_layout(Arrow(end=VeeHead(line_width=4, size=4, line_alpha=.5),
                                   x_start=j, y_start=(i - 1) + padding, x_end=j, y_end=i - padding, line_width=4, line_cap="round", line_alpha=.5))

    for i in range(arrows.shape[0]):
        p.add_layout(Arrow(end=VeeHead(line_width=4, size=4, line_alpha=.5),
                           x_start=arrows[i, 2] + padding, y_start=arrows[i, 3] + padding, x_end=arrows[i, 0] + padding, y_end=arrows[i, 1] + padding, line_width=4, line_cap="round", line_color="#DC143C", line_alpha=.5))
'''
    show(p)


if __name__ == '__main__':
    alpha = Seq("MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG\
KKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTP\
AVHASLDKFLASVSTVLTSKYR")
    beta = Seq("MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK\
VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFG\
KEFTPPVQAAYQKVVAGVANALAHKYH")
    penalty = {'MATCH': 1, 'MISMATCH': -1, 'GAP': -2}
    plot_nw(alpha, beta, penalty)
