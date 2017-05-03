import numpy as np
import matplotlib.pyplot as plt
from Bio.Seq import Seq

seq1 = Seq("CASA")
seq2 = Seq("MASA")
fig, ax = plt.subplots()
ax2 = ax.twiny()
ax3 = ax.twinx()

l_s1 = len(seq1)+1
l_s2 = len(seq2)+1

min_val_1, max_val_1 = 0, l_s1
ind_array_x = np.arange(min_val_1 + 0.5, max_val_1 + 0.5, 1.0)
#----------
min_val_2, max_val_2 = 0, l_s2
ind_array_y = np.arange(min_val_2 + 0.5, max_val_2 + 0.5, 1.0)
#----------
x, y = np.meshgrid(ind_array_x, ind_array_y, indexing='ij')

dummy_y = 0
for i, (x_val, y_val) in enumerate(zip(x.flatten(), y.flatten())):
    c = 'x' if i % 2 else 'o'
    an1 = ax.text(x_val, y_val, c, va='center', ha='center')
    #down
    ax.annotate("",
                xy=(x_val, y_val), xycoords='data',
                xytext=(x_val + 0.0, y_val - 0.5), textcoords='data',
                arrowprops=dict(arrowstyle="<-",
                                connectionstyle="arc3"),)
    #up
    ax.annotate("",
                xy=(x_val, y_val), xycoords='data',
                xytext=(x_val + 0.0, y_val + 0.5), textcoords='data',
                arrowprops=dict(arrowstyle="<-",
                                connectionstyle="arc3"),)
    #right
    ax.annotate("",
                xy=(x_val, y_val), xycoords='data',
                xytext=(x_val + 0.5, y_val - 0.0), textcoords='data',
                arrowprops=dict(arrowstyle="<-",
                                connectionstyle="arc3"),)
    #left
    ax.annotate("",
                xy=(x_val, y_val), xycoords='data',
                xytext=(x_val - 0.5, y_val - 0.0), textcoords='data',
                arrowprops=dict(arrowstyle="<-",
                                connectionstyle="arc3"),)
    dummy_x = i % l_s1
    print("(",dummy_x,",",dummy_y,")",end="")
    if(dummy_x == l_s1-1):
        dummy_y+=1

ax.set_xlim(min_val_1, max_val_1)
ax.set_ylim(min_val_2, max_val_2)
#--------------
ax.set_xticks(np.arange(max_val_1))
ax.set_xticklabels([])
ax.set_yticks(np.arange(max_val_2))
ax.set_yticklabels([])
ax.xaxis.tick_top()
ax.grid()
#-----------
ax2.set_xlim(ax.get_xlim())
ax3.set_ylim(ax.get_ylim())
seq1 = " "+seq1
seq2 = " "+seq2
ax2.set_xticks(np.arange(0.5,len(seq1),1))
ax2.set_xticklabels(list(seq1))
#-----------
ax3.set_yticks(np.arange(0.5,len(seq2),1))
ax3.set_yticklabels(list(seq2))
ax3.yaxis.tick_left()
plt.gca().invert_yaxis()

plt.draw()
plt.show()
