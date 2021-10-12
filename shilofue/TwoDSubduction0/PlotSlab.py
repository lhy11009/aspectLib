# -*- coding: utf-8 -*-
r"""Plot Slab Morphology

This exports: 

  -

This depends on:

  -  

Examples of usage:

  - default usage:

        python -m 

descriptions
""" 
import numpy as np
import sys, os, argparse
# import json, re
# import pathlib
# import subprocess
import numpy as np
# from matplotlib import cm
from matplotlib import pyplot as plt

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')

def Usage():
    print("\
(One liner description\n\
\n\
Examples of usage: \n\
\n\
  - default usage: \n\
\n\
        python -m \
        ")

class SLABPLOT():
    '''
    Plot slab morphology
    Inputs:
        -
    Returns:
        -
    '''
    def __init__(self):
        pass

    def ReadData(self, filename):
        data = np.loadtxt(filename)
        xs = data[:, 0]
        ys = data[:, 1]
        rs = (xs**2.0 + ys**2.0)**0.5
        thetas = np.arccos(xs/rs)
        indexes = np.argsort(thetas)
        self.rs = rs[indexes]
        self.thetas = thetas[indexes]
    
    def PlotSlabSurface(self, fileout):
        fig, ax = plt.subplots()
        ax.plot(self.rs * np.cos(self.thetas), self.rs * np.sin(self.thetas), label='slab surface')
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        plt.savefig(fileout)
        print("plot slab surface: ", fileout)
    
    def FindPoints(self, rf , depth, n_points, i_find=None, found=False):
        '''
        get the points at depth
        I use a iteration with a criteria of finding the right amount of points
        Inputs:
            ro: outer radis
            depth: depth of trench
            n_points: right amount of points, n_point[0] < number of found points < n_point[1]
        '''
        max_depth = 15e3  # max depth deviation from surface
        assert(type(n_points) == list and len(n_points) == 2)
        if found == True:
            # found before, only search previous results
            i_find_old = i_find.copy()  # initiate
            search_range = i_find_old
        else:
            # not found, search all indexes
            i_find_old = None
            search_range = range(self.rs.size)
        i_find = []
        # searching
        for i in search_range:
            r = self.rs[i]
            if np.abs(r-rf) < depth:
                i_find.append(i)
        i_find = np.array(i_find)
        if i_find.size > n_points[1]:
            # range too big
            theta, _ = self.FindPoints(rf, depth/2.0, n_points, i_find, True)
            return theta, True
        elif i_find.size < n_points[0]:
            # range too small
            if found == True:
                # already found in previous loop
                return self.thetas[i_find_old], True
            elif found == False and depth < max_depth:
                # not found, max limit is not reached
                theta, _found = self.FindPoints(rf, depth*2.0, n_points, None, False)
                return theta, _found
            else:
                # not found, max limit is reached
                return None, False
            pass
        else:
            # find and right range
            # print(self.thetas[i_find]) # Debug
            return self.thetas[i_find], True


def slab_morph(file_path):
    '''
    Plot slab morphology
    Inputs:
        file_path(str):
            slab surface from vtk output
    '''
    outputs = {}
    SlabPlot = SLABPLOT()
    SlabPlot.ReadData(file_path)
    SlabPlot.PlotSlabSurface(os.path.join(ASPECT_LAB_DIR, 'output', 'slab.png'))
    ro = np.max(SlabPlot.rs)
    ri = 2890e3 # todo: modify with prm
    thetas, found = SlabPlot.FindPoints(ro, 0.1e3, [0, 5])
    if found == False:
        raise ValueError("slab_morph: trench hinge cannot be found")
    trench_theta = np.mean(thetas)
    outputs['trench'] = {'theta': trench_theta}  # append to output
    morp_rs = []  # initiate
    morp_thetas = []
    num = 100
    for i in range(num):
        rf = ro - 1.0 * i / (num-1) * (ro - ri)
        # morphology
        thetas, found = SlabPlot.FindPoints(rf, 0.1e3, [4, 10])
        if found == False:
            slab_depth = ro - rf
            break
        morp_rs.append(rf)
        morp_thetas.append(thetas)
    outputs['morp'] = {'radius': morp_rs, 'theta': morp_thetas}  # append to output
    return outputs


def vtk_and_slab_morph():
    '''
    run vtk and read in slab morph
    '''



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
    parser = argparse.ArgumentParser(description='Parse parameters')
    parser.add_argument('-i', '--inputs', type=str,
                        default='',
                        help='Some inputs')
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    # commands
    if (_commend in ['-h', '--help']):
        # example:
        Usage()
    elif _commend == 'morph':
        # slab_morphology, inputs is the vtk output file:
        slab_morph(arg.inputs)
    elif _commend == 'morph_case':
        # slab_morphology, input is the case name
        # example:
        file_path = os.path.join(arg.inputs, 'vtk_outputs', 'contour_slab%d.txt' % arg.snapshot)
        slab_morph(arg.inputs)
    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()