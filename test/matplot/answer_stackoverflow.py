import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
np.random.seed(5)
plt.rcParams["figure.figsize"] = 4, 5
param = {"grid.linewidth": 1.6,
         "grid.color": "lightgray",
         "axes.linewidth": 1.6,
         "axes.edgecolor": "lightgray"}
plt.rcParams.update(param)

# Data
headh = list("GATCCA")
headv = list("GTGCCT")

v = np.zeros((7, 7), dtype=int)
v[1:, 1:] = np.random.randint(-2, 7, size=(6, 6))
arrows = np.random.randint(0, v.shape[1], size=(14, 4))
print(arrows)
opt = np.array([(0, 1), (1, 0), (1, 1)])
arrows[:, 2:] = arrows[:, :2] + opt[np.random.randint(0, 3, size=14)]
print(arrows)
arrowsb = np.random.randint(0, v.shape[1], size=(7, 4))
optb = np.array([(0, 1), (1, 0), (1, 1)])
arrowsb[:, 2:] = arrowsb[:, :2] + optb[np.random.randint(0, 3, size=7)]

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


arrowprops = dict(facecolor='crimson', alpha=0.5, lw=0,
                  shrink=0.2, width=2, headwidth=7, headlength=7)
for i in range(arrows.shape[0]):
    ax.annotate("", xy=arrows[i, 2:],#origin
                xytext=arrows[i, :2], arrowprops=arrowprops)
arrowprops.update(facecolor='blue')
for i in range(arrowsb.shape[0]):
    ax.annotate("", xy=arrowsb[i, 2:],
                xytext=arrowsb[i, :2], arrowprops=arrowprops)
plt.show()
