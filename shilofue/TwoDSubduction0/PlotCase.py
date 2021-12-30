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
from shilofue.PlotVisit import RunScripts
import shilofue.PlotCombine as PlotCombine
from shilofue.TwoDSubduction0.PlotVisit import VISIT_OPTIONS, PREPARE_RESULT_OPTIONS
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
 -t 0.0 -t1 0.5e6\n\
\n\
  - plot cases in a directory (loop), same options as before:\n\
        Lib_TwoDSubduction0_PlotCase  plot_case_in_dir -i ~/ASPECT_PROJECT/TwoDSubduction/EBA_CDPT\n\
\n\
  - prepare result of a single step by combining figures:\n\
        Lib_TwoDSubduction0_PlotCase prepare_result_step -i ~/ASPECT_PROJECT/TwoDSubduction/test_peierls1/peierls\
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
    prm_path = os.path.join(case_path, 'output', 'original.prm')

    # plot visit
    Visit_Options = VISIT_OPTIONS(case_path)
    Visit_Options.Interpret(last_steps=3)  # interpret scripts, plot the last 3 steps
    odir = os.path.join(case_path, 'visit_scripts')
    if not os.path.isdir(odir):
        os.mkdir(odir)
    if Visit_Options.get_geometry() == 'chunk':
        py_script = 'slab_sph.py'
    elif Visit_Options.get_geometry() == 'box':
        py_script = 'slab_cart.py'
    else:
        raise ValueError('%s: no option related to geometry : %s' % (Utilities.func_name(), Visit_Options.get_geometry()))
    ofile = os.path.join(odir, py_script)
    visit_script = os.path.join(ASPECT_LAB_DIR, 'visit_scripts', 'TwoDSubduction', py_script)
    visit_script_base = os.path.join(ASPECT_LAB_DIR, 'visit_scripts', 'base.py')
    Visit_Options.read_contents(visit_script_base, visit_script)
    Visit_Options.substitute()
    ofile_path = Visit_Options.save(ofile, relative=True)
    print("Visualizing using visit")
    RunScripts(ofile_path)  # run scripts


def PrepareResultStep(case_path, step):
    '''
    Prepare results
    Inputs:
        case_path(str): path to the case
        step(int): step in computation
    '''
    # Generate json file
    pr_script = os.path.join(ASPECT_LAB_DIR, "files", "TwoDSubduction", "figure_step_template.json")
    odir = os.path.join(case_path, 'json_files')
    ofile = os.path.join(odir, 'figure_step.json')
    assert(os.path.isfile(pr_script))
    Prepare_Result = PREPARE_RESULT_OPTIONS(case_path)
    Prepare_Result.Interpret(step=step)
    Prepare_Result.read_contents(pr_script)
    Prepare_Result.substitute()
    ofile_path = Prepare_Result.save(ofile, relative=True)
    PlotCombine.PrepareResults(ofile_path)  # call function with the generated json file


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
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    # commands
    if _commend == 'plot_case':
        PlotCase.PlotCaseCombined([PlotCase.PlotCaseRun, PlotCaseRun], arg.inputs, [arg.time, arg.time1])
    elif _commend == 'plot_case_in_dir':
        PlotCase.PlotCaseCombinedDir([PlotCase.PlotCaseRun, PlotCaseRun], arg.inputs, [arg.time, arg.time1])
        pass
    elif _commend == 'prepare_result_step':
        PrepareResultStep(arg.inputs, arg.step)
    elif (_commend in ['-h', '--help']):
        # example:
        Usage()

# run script
if __name__ == '__main__':
    main()
