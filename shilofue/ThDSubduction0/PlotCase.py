# -*- coding: utf-8 -*-
r"""(one line description)

This exports:

  -

This depends on:

  -

Examples of usage:

  - default usage: plot case running results

        Lib_FOO0_PlotCase plot_case -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear32/eba1_MRf12_iter20

descriptions
    First replate FOO with the name of the project
"""
import numpy as np
import sys, os, argparse, re
import json, re
# import pathlib
# import subprocess
import numpy as np
# from matplotlib import cm
from matplotlib import pyplot as plt
from shilofue.ParsePrm import ReadPrmFile
from shilofue.PlotVisit import RunScripts
import shilofue.PlotCombine as PlotCombine
from shilofue.ThDSubduction0.PlotVisit import VISIT_OPTIONS, PREPARE_RESULT_OPTIONS
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
  - default usage: plot case running results, -t option deals with a time range, default is a whole range.\
-ti option deals with a time interval and is useful for animation.\n\
\n\
        Lib_ThDSubduction0_PlotCase plot_case -i  `pwd`\
 -t 0.0 -t1 0.5e6 -ti 0.1e6\n\
\n\
  - plot cases in a directory (loop), same options as before:\n\
        Lib_ThDSubduction0_PlotCase  plot_case_in_dir -i `pwd`\n\
\n\
  - prepare result of a single step by combining figures:\n\
        Lib_ThDSubduction0_PlotCase prepare_result_step -i `pwd`\n\
\n\
  - plot and then prepare result of a single step by combining figures:\n\
        Lib_ThDSubduction0_PlotCase plot_prepare_result_step -i `pwd`\n\
\n\
  - generate animation: \n\
        Lib_ThDSubduction0_PlotCase animate_case -i `pwd`\n\
\n\
  - generate animation for cases in a directory: \n\
        Lib_ThDSubduction0_PlotCase animate_case_in_dir -i `pwd`\n\
\n\
        ")


def PlotCaseRun(case_path, **kwargs):
    '''
    Plot case run result
    Inputs:
        case_path(str): path to the case
        kwargs:
            time_range
            step(int): if this is given as an int, only plot this step
    Returns:
        -
    '''
    print("PlotCaseRun in ThDSubduction: operating")
    time_interval = kwargs.get('time_interval', None)
    # get case parameters
    prm_path = os.path.join(case_path, 'output', 'original.prm')
    # plot with paraview
    # initiate class object
    Paraview_Options = VISIT_OPTIONS(case_path)
    # call function
    Paraview_Options.Interpret(time_interval=time_interval)
    # ofile = os.path.join('visit_scripts', 'slab_sph.py')
    ofile = os.path.join(case_path, 'paraview_scripts', 'slab.py')
    paraview_script = os.path.join(ASPECT_LAB_DIR, 'paraview_scripts',"ThDSubduction", 'slab.py')
    paraview_base_script = os.path.join(ASPECT_LAB_DIR, 'paraview_scripts', 'base.py')  # base.py : base file
    Paraview_Options.read_contents(paraview_base_script, paraview_script)  # this part combines two scripts
    Paraview_Options.substitute()  # substitute keys in these combined file with values determined by Interpret() function
    ofile_path = Paraview_Options.save(ofile, relative=False)  # save the altered script


class PLOTTER(PlotCase.PLOTTER):
    '''
    A class for preparing results
    '''
    def __init__(self, module, PlotCaseRuns, **kwargs):
        '''
        Initiation
        module (class) - class to use for generating results at each step
        kwargs(dict):

        '''
        PlotCase.PLOTTER.__init__(self, module, PlotCaseRuns)   


def PrScriptToUse(case_path, default):
    '''
    Return the script to use for preparing results
    This is mainly used for generating figures for animation
    Inputs:
    '''
    pr_script = default
    return pr_script


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
    parser.add_argument('-s', '--step', type=int,
                        default=0,
                        help='step')
    parser.add_argument('-r', '--rewrite', type=int,
                        default=0,
                        help='If rewrite previous result')
    parser.add_argument('-ti', '--time_interval', type=float,
                        default=None,
                        help='Time interval, affecting the time steps to visualize')
    
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    default = os.path.join(ASPECT_LAB_DIR, "files", "ThDSubduction", "figure_step_template.json")

    # commands
    # rearrange the time_range entry
    if arg.time == None or arg.time1 == None:
        time_range = None
    else:
        assert(type(arg.time) == float and type(arg.time1) == float)
        time_range = [arg.time, arg.time1]
    if _commend == 'plot_case':
        # PlotCase.PlotCaseCombined([PlotCase.PlotCaseRun, PlotCaseRun], arg.inputs, time_range=time_range, time_interval=arg.time_interval)
        PlotCase.PlotCaseCombined([PlotCaseRun], arg.inputs, time_range=time_range, time_interval=arg.time_interval)
    elif _commend == 'plot_case_in_dir':
        PlotCase.PlotCaseCombinedDir([PlotCase.PlotCaseRun, PlotCaseRun], arg.inputs, time_range=time_range)
        pass
    elif _commend == 'prepare_result_step':
        pr_script = PrScriptToUse(arg.inputs, default)
        Plotter = PLOTTER(PREPARE_RESULT_OPTIONS, [PlotCase.PlotCaseRun, PlotCaseRun])
        Plotter.PrepareResultStep(arg.inputs, pr_script, arg.step)
    elif _commend == 'plot_prepare_result_step':
        pr_script = PrScriptToUse(arg.inputs, default)
        Plotter = PLOTTER(PREPARE_RESULT_OPTIONS, [PlotCase.PlotCaseRun, PlotCaseRun])
        Plotter.PlotPrepareResultStep(arg.inputs, pr_script, arg.step)
    elif _commend == 'animate_case':
        pr_script = PrScriptToUse(arg.inputs, default)
        Plotter = PLOTTER(PREPARE_RESULT_OPTIONS, [PlotCase.PlotCaseRun, PlotCaseRun])
        PlotCase.AnimateCaseResults(Plotter.PlotPrepareResultStep, arg.inputs, pr_script, time_interval=arg.time_interval)
    elif _commend == 'animate_case_in_dir':
        pr_script = PrScriptToUse(arg.inputs, default)
        Plotter = PLOTTER(PREPARE_RESULT_OPTIONS, [PlotCase.PlotCaseRun, PlotCaseRun])
        PlotCase.AnimateCombinedDir(Plotter.PlotPrepareResultStep, arg.inputs, pr_script)
    elif (_commend in ['-h', '--help']):
        # example:
        Usage()
    else:
        raise ValueError("Invalid command: use -h for help information")

# run script
if __name__ == '__main__':
    main()
