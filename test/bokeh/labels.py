from bokeh.plotting import figure, show, output_file
from bokeh.models import ColumnDataSource, Range1d, LabelSet, Label, FixedTicker
from bokeh.models import TickFormatter, Arrow, NormalHead, VeeHead,Text
from bokeh.core.properties import Dict, Int, String
import numpy as np
from Bio.Seq import Seq


matrix = np.zeros((2, 2))

source = ColumnDataSource(data=dict(x=[str(x) for x in range(6)],
                                    y=[str(y) for y in range(6)],
                                    text_v=["100" for y in range(6)]))
seq1 = Seq("ABCDEF")
p = figure(title='Dist. of 10th Grade Students at Lee High',
           x_axis_location="above")
p.y_range = Range1d(len(seq1), -1.5)



for i in range(len(seq1)+1):
    for j in range(len(seq1)+1):
        p.rect(i, j, 0.9, 0.9,
            fill_alpha=0.6, color="#CAB2D6")
        p.add_layout(Arrow(end=VeeHead(line_width=4, size=4, line_alpha=.5),
                   x_start=i + 0.3, y_start=j + 0.3, x_end=i + 1 - 0.3, y_end=j + 1 - 0.3, line_width=4, line_cap="round", line_alpha=.5))



text_props = {
    
    "angle": 0,
    "text_color": "black",
    "text_align": "center",
    "text_baseline": "middle"
}

for i, l in enumerate(seq1):
    glyph = Text(x=-1,y=i + 1, text=[l], text_font_size="9pt", **text_props)
    p.add_glyph(glyph)
for i, l in enumerate(seq1):
    glyph = Text(x=i + 1,y=-1, text=[l], text_font_size="9pt", **text_props)
    p.add_glyph(glyph)

show(p)
