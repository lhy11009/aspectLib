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
import shilofue.Utilities as Utilities
import shilofue.PlotRunTime as PlotRunTime
import shilofue.PlotStatistics as PlotStatistics

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')

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


def PlotCaseRun(case_path, **kwargs):
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
    inputs = ReadPrmFile(prm_path)
    # get solver scheme
    solver_scheme = inputs.get('Nonlinear solver scheme', 'single Advection, single Stokes')

    log_file = os.path.join(case_path, 'output', 'log.txt')
    assert(os.access(log_file, os.R_OK))
    statistic_file = os.path.join(case_path, 'output', 'statistics')
    assert(os.access(statistic_file, os.R_OK))

    # statistic
    print('Ploting statistic results')
    fig_path = os.path.join(case_path, 'img', 'Statistic.png')
    PlotStatistics.PlotFigure(statistic_file, fig_path)

    # run time
    print('Ploting run time')
    fig_path = os.path.join(case_path, 'img', 'run_time.png')
    if not os.path.isdir(os.path.dirname(fig_path)):
        os.mkdir(os.path.dirname(fig_path))
    # time range
    try:
        time_range = kwargs['time_range']
        fig_output_path, step_range = PlotRunTime.PlotFigure(log_file, fig_path, fix_restart=True, time_range=time_range)
    except KeyError:
        fig_output_path = PlotRunTime.PlotFigure(log_file, fig_path, fix_restart=True)
        step_range = None

    # Newton history
    # determine whether newton is used
    print('solver_scheme, ', solver_scheme) # debug
    fig_path = os.path.join(case_path, 'img', 'newton_solver_history.png')
    # match object in solver scheme, so we only plot for these options
    match_obj = re.search('Newton', solver_scheme)
    match_obj1 = re.search('iterated defect correction Stokes', solver_scheme)
    if match_obj or match_obj1:
        print("Plotting newton solver history")
        PlotRunTime.PlotNewtonSolverHistory(log_file, fig_path, step_range=step_range)
    else:
        print("Skipping newton solver history")


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
            PlotCaseRun(arg.inputs, time_range=[arg.time, arg.time1])
        else:
            PlotCaseRun(arg.inputs)
    elif (_commend in ['-h', '--help']):
        # example:
        Usage()

# run script
if __name__ == '__main__':
    main()
