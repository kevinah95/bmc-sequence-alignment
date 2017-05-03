from __future__ import print_function
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import numpy as np

ax = plt.axes([0.1, 0.3, 0.5, 0.5])
ax.xaxis.tick_top()
seq1 = Seq(" CASA")
seq2 = Seq(" MASA")

ax.pcolormesh(np.array([[1, 2], [3, 4]]))
plt.yticks(range(len(seq1)), list(seq1))
plt.gca().invert_yaxis()
plt.xticks(range(len(seq2)), list(seq2))
plt.ylabel("My y-label")
plt.xlabel("Check saved figures for their bboxes $^\searrow$")

plt.show()