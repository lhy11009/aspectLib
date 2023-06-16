# -*- coding: utf-8 -*-
r"""(one line description)
generate simple plots
"""

import os # edit
import numpy as np # edit
from matplotlib import pyplot as plt

cores = np.array([1020, 1428, 1700])
times = np.array([8.94e+03, 1.06e+04, 1.54e+04])
stokes_dofs = np.array([27436497, 46431861, 84356115]) 
x0 = 1.0
unit_time0 = cores[0] * times[0] / stokes_dofs[0]
unit_times = cores * times / stokes_dofs
xs = cores 
ys = unit_time0 / unit_times
_color = "tab:blue"
_label = "Stokes DOFs = 27.4 Ma, 46.4 Ma, 84.3 Ma"
x_label = "cores"
y_label = "parallel efficiency"
file_name = "stampede2_affinity_test_06082023.png"


fig, ax= plt.subplots(tight_layout=True, figsize=(5, 5))  # plot of wallclock
ax.plot(xs, ys, ".-", color=_color, label=_label)
ax.grid()
ax.set_xlabel(x_label)
ax.set_ylabel(y_label)
ax.legend()
fig.savefig(file_name)
print("file saved: ", file_name)