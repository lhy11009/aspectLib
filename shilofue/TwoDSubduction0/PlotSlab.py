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
from shilofue.PlotVisit import PrepareVTKOptions, RunVTKScripts, VISIT_OPTIONS
import shilofue.Utilities as Utilities

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
    plan to be deprecated an move into cpp
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


def vtk_and_slab_morph(case_dir, pvtu_step, **kwargs):
    '''
    run vtk and read in slab morph
    Inputs:
        case_dir (str): case directory
        pvtu_step (int): time step
        kwargs:
            new: remove old output file
    '''
    print("pvtu_step: %s\n" % str(pvtu_step))
    output_dir = os.path.join(case_dir, 'vtk_outputs')
    output_file = os.path.join(output_dir, 'slab_morph.txt')
    vtk_option_path, _time, step = PrepareVTKOptions(case_dir, 'TwoDSubduction_SlabAnalysis', vtk_step=pvtu_step)
    contour_file = RunVTKScripts('TwoDSubduction_SlabAnalysis', vtk_option_path)
    slab_morph_outputs = slab_morph(contour_file)
    print("Trench: %s" % slab_morph_outputs['trench']['theta'])
    # header
    file_header = "# %s\n# %s\n# %s\n# %s\n# %s\n" % \
    ("pvtu_step", "step", "time (yr)", "trench (rad)", "minimum_radius (km)")
    trench_theta = slab_morph_outputs['trench']['theta']
    minimum_radius = slab_morph_outputs['morp']['radius'][-1]
    outputs = "%-12s%-12d%-14.4e%-14.4e%-14.4e\n" % (pvtu_step, step, _time, trench_theta, minimum_radius)
    # remove old file
    is_new = kwargs.get('new', False)
    if is_new and os.path.isfile(output_file):
        os.remove(output_file)
    # output data
    if not os.path.isfile(output_file):
        with open(output_file, 'w') as fout:
            fout.write(file_header)
            fout.write(outputs)
        print('Create output: %s' % output_file)
    else:
        with open(output_file, 'a') as fout:
            fout.write(outputs)
        print('Update output: %s' % output_file)


def vtk_and_slab_morph_case(case_dir):
    '''
    run vtk and get outputs for every snapshots
    '''
    # get all available snapshots
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    available_pvtu_snapshots = Utilities.string2list(Visit_Options.options['ALL_AVAILABLE_GRAPHICAL_SNAPSHOTS'])
    available_pvtu_steps = [i - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']) for i in available_pvtu_snapshots]
    # get slab morphology
    for pvtu_step in available_pvtu_steps:
        if pvtu_step == 0:
            # start new file with the 0th step
            vtk_and_slab_morph(case_dir, pvtu_step, new=True)
        else:
            vtk_and_slab_morph(case_dir, pvtu_step)
    pass
    

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
    parser.add_argument('-s', '--step', type=str,
                        default='0',
                        help='Timestep')
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
    elif _commend == 'morph_step':
        # slab_morphology, input is the case name
        # example:
        vtk_and_slab_morph(arg.inputs, arg.step)
    elif _commend == 'morph_case':
        # slab_morphology, input is the case name
        # example:
        vtk_and_slab_morph_case(arg.inputs)
    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()