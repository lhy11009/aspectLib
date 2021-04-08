# -*- coding: utf-8 -*-
r"""plot the shape of an ellipse slab

This exports: 

  - a figure:
        initial_slab.png

This depends on:

  - cpp executable: 
        os.path.join(ASPECT_LAB_DIR, 'cpp_scripts', 'ellipse_distance.o')

Examples of usage:

  - default usage:

        python PlotEllipeSlab.py plot 2.0 1.0 0.1
            here a = 2.0, b = 1.0, d = 0.1

descriptions
""" 
import numpy as np
import sys, os, argparse
# import json, re
# import pathlib
import subprocess
import numpy as np
# from matplotlib import cm
from matplotlib import pyplot as plt

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')
# executables this one depends on
# compute distance to the ellipse surface
cpp_edist = os.path.join(ASPECT_LAB_DIR, 'cpp_scripts', 'ellipse_distance.o')


def SurfEllipse(a, b, d):
    '''
    shape of slab
    Inputs:
        a - semi-axis, x
        b - semi-axis, y
        d - depth of the points from the surface
    Returns:
        surfs - surface
        bots - bottom
    '''
    n = 100
    surfs = np.zeros((2, n))
    # bots = np.zeros((2, n))
    internals_x = []
    internals_y = []


    # ellipse surface
    surfs[0, :] = np.linspace(0.0, a, n)
    surfs[1, :] = b * (1 - surfs[0, :]**2.0 / a**2.0 )**0.5
    
    # ellipse internal
    for xp in np.linspace(0, a, 200):
        for yp in np.linspace(0, b, 100):
            if (xp**2.0 / a**2.0 + yp**2.0 / b**2.0 > 1.0):
                # check point is in the ellipse
                continue
            outputs = subprocess.check_output([cpp_edist, str(a), str(b), str(xp), str(yp)])
            distance = float(outputs.split()[1])
            if distance < d:
                # check point is in the slab
                internals_x.append(xp)
                internals_y.append(yp)

    # generate output
    internals = np.array([internals_x, internals_y])
    
    return surfs, internals


def main():
    '''
    main function of this module
    Inputs:
        sys.arg[1](str):
            commend
        sys.arg[2, :](str):
            options
    '''
    _commend = sys.argv[1]
    # parse options
    # parser = argparse.ArgumentParser(description='Parse parameters')
    # parser.add_argument('-i', '--inputs', type=str,
    #                     default='',
    #                    help='Some inputs')
    _options = []
    # try:
    #    _options = sys.argv[2: ]
    # except IndexError:
    #    pass
    # arg = parser.parse_args(_options)

    # commands
    if _commend == 'plot':
        a = float(sys.argv[2])
        b = float(sys.argv[3])
        d = float(sys.argv[4])
        # get surface of the ellipse
        surfs, internals = SurfEllipse(a, b, d)
        # plot
        fig, ax = plt.subplots()
        ax.plot(surfs[0, :], surfs[1, :], 'b')
        ax.plot(internals[0, :], internals[1, :], 'k.')
        filename = "initial_slab.png"
        ax.axis('scaled')
        fig.savefig(filename)

# run script
if __name__ == '__main__':
    main()