# -*- coding: utf-8 -*-
r"""(one line description)
generate simple plots
"""

import os # edit
import numpy as np # edit
from matplotlib import pyplot as plt

xs = []
ys = []
_color = "tab:blue"
_label = ""
x_label = 
y_label =
file_name = "foo.png"

fig, ax= plt.subplots(tight_layout=True, figsize=(5, 5))  # plot of wallclock
ax.plot(xs, ys, ".", color=_color, label=_label)
ax.grid()
ax.set_xlabel(x_label)
ax.set_ylabel(y_label)
ax.legend()
fig.savefig(file_name)
print("file saved: ", file_name)