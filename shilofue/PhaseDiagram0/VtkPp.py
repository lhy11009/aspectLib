# -*- coding: utf-8 -*-
r"""(one line description)

This exports: 

  -

This depends on:

  -  

descriptions

    This script handles analysis using the vtk package for this project
"""
#### 3rd parties 
import numpy as np
import sys, os, argparse
# import json, re
# import pathlib
# import subprocess
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from numpy import linalg as LA 
import multiprocessing
#### self
from shilofue.PlotVisit import PrepareVTKOptions, RunVTKScripts, PARALLEL_WRAPPER_FOR_VTK
from shilofue.FOO0.PlotVisit import VISIT_OPTIONS
from shilofue.ParsePrm import ReadPrmFile
from shilofue.Plot import LINEARPLOT
import shilofue.VtkPp as VtkPp
# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')
# import utilities in subdirectiory
sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities


def Usage():
    print("\
(One liner description\n\
\n\
Examples of usage: \n\
")


class VTKP(VtkPp.VTKP):
    '''
    Class inherited from a parental class
    Attributes:
    '''
    def __init__(self, **kwargs):
        VtkPp.VTKP.__init__(self)  # initiation of parental class


def WrapperForStep(case_dir, vtu_step, **kwargs):
    '''
    Wrapper for using PVTK class to analyze one step
    Inputs:
        case_dir (str): case directory
        vtu_step (int): step in vtu outputs
    '''
    filein = os.path.join(case_dir, "output", "solution", "solution-%05d.pvtu" % vtu_step)
    assert(os.path.isfile(filein))
    vtk_option_path, _time, step = PrepareVTKOptions(VISIT_OPTIONS, case_dir, 'foo',\
    vtu_step=vtu_step, include_step_in_filename=True, generate_horiz=True)
    VtkP = VTKP()
    VtkP.ReadFile(filein)
    field_names = ['T', 'density']
    VtkP.ConstructPolyData(field_names, include_cell_center=True)
    # generate outputs
    outputs = "%-12s%-12d%-14.4e\n"\
    % (vtu_step, step, _time)
    print(outputs) # debug
    return vtu_step, outputs


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
    parser.add_argument('-o', '--outputs', type=str,
                        default='',
                        help='Some outputs')
    parser.add_argument('-s', '--step', type=int,
                        default=0,
                        help='step')
    parser.add_argument('-vs', '--vtu_step', type=int,
                        default=0,
                        help='vtu_step')
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
    elif _commend == 'foo':
        # use the wrapper
        WrapperForStep(arg.inputs, int(arg.vtu_step), rewrite=1)
    else:
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()