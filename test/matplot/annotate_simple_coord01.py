import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots()

min_val, max_val = 0, 10
ind_array = np.arange(min_val + 0.5, max_val + 0.5, 1.0)
print(np.zeros((5, 5)))
print(ind_array)
x, y = np.meshgrid(ind_array, ind_array)

print(np.meshgrid(ind_array, ind_array))

for i, (x_val, y_val) in enumerate(zip(x.flatten(), y.flatten())):
    c = 'x' if i%2 else 'o' 
    an1 = ax.text(x_val, y_val, c, va='center', ha='center')
    ax.annotate("",
            xy=(x_val, y_val), xycoords='data',
            xytext=(x_val+0.5, y_val-0.5), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3"), 
            )

ax.set_xlim(min_val, max_val)
ax.set_ylim(min_val, max_val)
ax.set_xticks(np.arange(max_val))
ax.set_yticks(np.arange(max_val))
ax.grid()

plt.draw()
plt.show()