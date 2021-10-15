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
from shilofue.Plot import LINEARPLOT
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

class SLABPLOT(LINEARPLOT):
    '''
    Plot slab morphology
    Inputs:
        -
    Returns:
        -
    '''
    def PlotSlabSurface(self, fileout):
        fig, ax = plt.subplots()
        ax.plot(self.rs * np.cos(self.thetas), self.rs * np.sin(self.thetas), label='slab surface')
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        plt.savefig(fileout)
        print("plot slab surface: ", fileout)


def slab_morph(inputs):
    '''
    Plot slab morphology
    Inputs:
        inputs (str): the outputs from another executable (slab analysis)
    '''
    outputs = {}
    outputs['trench_theta'] = float(Utilities.re_read_variable_from_string(inputs, 'Trench theta', ':'))
    outputs['slab_depth'] = float(Utilities.re_read_variable_from_string(inputs, 'Slab depth', ':'))
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
    _stdout = RunVTKScripts('TwoDSubduction_SlabAnalysis', vtk_option_path)
    slab_outputs = slab_morph(_stdout)
    # header
    file_header = "# %s\n# %s\n# %s\n# %s\n# %s\n" % \
    ("1: pvtu_step", "2: step", "3: time (yr)", "4: trench (rad)", "5: slab depth (m)")
    # output string
    outputs = "%-12s%-12d%-14.4e%-14.4e%-14.4e\n" % (pvtu_step, step, _time, slab_outputs['trench_theta'], slab_outputs['slab_depth'])
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
    # todo get where previous session ends
    SlabPlot = SLABPLOT('slab')
    slab_morph_file = os.path.join(case_dir, 'vtk_outputs', 'slab_morph.txt')
    SlabPlot.ReadHeader(slab_morph_file)
    SlabPlot.ReadData(slab_morph_file)
    col_pvtu_step = SlabPlot.header['pvtu_step']['col']
    last_pvtu_step = int(SlabPlot.data[-1, col_pvtu_step])
    # get slab morphology
    for pvtu_step in available_pvtu_steps:
        if pvtu_step <= last_pvtu_step:
            # skip existing steps
            continue
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
    elif _commend == 'morph_step':
        # slab_morphology, input is the case name
        # example:
        vtk_and_slab_morph(arg.inputs, int(arg.step))
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