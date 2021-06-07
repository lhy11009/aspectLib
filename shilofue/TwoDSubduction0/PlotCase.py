# -*- coding: utf-8 -*-
r"""(one line description)

This exports: 

  -

This depends on:

  -  

Examples of usage:

  - default usage: plot case running results

        python -m shilofue.TwoDSubduction0.PlotCase plot_case -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear32/eba1_MRf12_iter20

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
import shilofue.PlotRunTime as PlotRunTime
import shilofue.PlotStatistics as PlotStatistics

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')

def PlotCaseRun(case_path):
    '''
    Plot case run result
    Inputs:
        case_path(str): path to the case
    Returns:
        -
    '''
    log_file = os.path.join(case_path, 'output', 'log.txt')
    assert(os.access(log_file, os.R_OK))
    statistic_file = os.path.join(case_path, 'output', 'statistics')
    assert(os.access(statistic_file, os.R_OK))
    
    # statistic
    fig_path = os.path.join(case_path, 'img', 'Statistic.png')
    PlotStatistics.PlotFigure(statistic_file, fig_path)

    # run time
    fig_path = os.path.join(case_path, 'img', 'run_time.png')
    if not os.path.isdir(os.path.dirname(fig_path)):
        os.mkdir(os.path.dirname(fig_path))
    PlotRunTime.PlotFigure(log_file, fig_path, fix_restart=True)

    # Newton history
    fig_path = os.path.join(case_path, 'img', 'newton_solver_history.png')
    PlotRunTime.PlotNewtonSolverHistory(log_file, fig_path) 


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
    if _commend == 'plot_case':
        # example:
        PlotCaseRun(arg.inputs)

# run script
if __name__ == '__main__':
    main()