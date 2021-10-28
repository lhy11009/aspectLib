# -*- coding: utf-8 -*-
r"""(one line description)

This exports:

  -

This depends on:

  -

Examples of usage:

  - default usage: plot case running results

        Lib_TwoDSubduction0_PlotCase plot_case -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear32/eba1_MRf12_iter20

descriptions
"""
import numpy as np
import sys, os, argparse, re
# import json, re
# import pathlib
# import subprocess
import numpy as np
# from matplotlib import cm
from matplotlib import pyplot as plt
from shilofue.ParsePrm import ReadPrmFile
import shilofue.PlotVisit as PlotVisit
import shilofue.PlotCase as PlotCase
import shilofue.PlotRunTime as PlotRunTime
import shilofue.PlotStatistics as PlotStatistics

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')

sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities

def Usage():
    print("\
This scripts generate plots for a single case in TwoDSubduction project\n\
\n\
Examples of usage: \n\
\n\
  - default usage: plot case running results, -t option deals with a time range, default is a whole range\n\
\n\
        Lib_TwoDSubduction0_PlotCase plot_case -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear32/eba1_MRf12_iter20\
          -t 0.0 -t1 0.5e6\
        ")


def TwoDSubduction_PlotCaseRun(case_path, **kwargs):
    '''
    Plot case run result
    Inputs:
        case_path(str): path to the case
        kwargs:
            time_range
    Returns:
        -
    '''
    # get case parameters
    prm_path = os.path.join(case_path, 'case.prm')

    # plot visit
    Visit_Options = PlotVisit.VISIT_OPTIONS(case_path)
    Visit_Options.Interpret(last_steps=3)  # interpret scripts, plot the last 3 steps
    odir = os.path.join(case_path, 'visit_scripts')
    if not os.path.isdir(odir):
        os.mkdir(odir)
    ofile = os.path.join(odir, 'slab.py')
    visit_script = os.path.join(ASPECT_LAB_DIR, 'visit_scripts', 'TwoDSubduction', 'slab.py')
    Visit_Options.read_contents(visit_script)
    Visit_Options.substitute()
    ofile_path = Visit_Options.save(ofile, relative=True)
    print("Visualizing using visit")
    PlotVisit.RunScripts(ofile_path)  # run scripts


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
    parser.add_argument('-t', '--time', type=float,
                        default=None,
                        help='Time')
    parser.add_argument('-t1', '--time1', type=float,
                        default=None,
                        help='Time1')
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    # commands
    if _commend == 'plot_case':
        # example:
        if arg.time != None and arg.time1 != None:
            PlotCase.PlotCaseRun(arg.inputs, time_range=[arg.time, arg.time1])
            TwoDSubduction_PlotCaseRun(arg.inputs, time_range=[arg.time, arg.time1])
        else:
            PlotCase.PlotCaseRun(arg.inputs)
            TwoDSubduction_PlotCaseRun(arg.inputs)
    elif (_commend in ['-h', '--help']):
        # example:
        Usage()

# run script
if __name__ == '__main__':
    main()
