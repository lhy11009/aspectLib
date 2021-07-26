import numpy as np
from matplotlib import pyplot as plt

A=np.loadtxt("contour.xyz")
print(A)

fig, ax = plt.subplots()
ax.plot(A[:,0], A[:, 1], '.')
fig.savefig('temp.png')

