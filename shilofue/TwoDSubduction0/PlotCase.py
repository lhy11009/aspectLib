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
import json, re
# import pathlib
# import subprocess
import numpy as np
# from matplotlib import cm
from matplotlib import pyplot as plt
from shilofue.ParsePrm import ReadPrmFile
from shilofue.PlotVisit import RunScripts
import shilofue.PlotCombine as PlotCombine
from shilofue.TwoDSubduction0.PlotVisit import VISIT_OPTIONS, PREPARE_RESULT_OPTIONS
from shilofue.TwoDSubduction0.PlotSlab import vtk_and_slab_morph_case, SLABPLOT
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
        Lib_TwoDSubduction0_PlotCase plot_case -i  ~/ASPECT_PROJECT/TwoDSubduction/EBA_CDPT1/eba_cdpt_SA80.0_OA40.0\
 -t 0.0 -t1 0.5e6\n\
\n\
  - plot cases in a directory (loop), same options as before:\n\
        Lib_TwoDSubduction0_PlotCase  plot_case_in_dir -i ~/ASPECT_PROJECT/TwoDSubduction/EBA_CDPT\n\
\n\
  - prepare result of a single step by combining figures:\n\
        Lib_TwoDSubduction0_PlotCase prepare_result_step -i ~/ASPECT_PROJECT/TwoDSubduction/test_peierls1/peierls\n\
\n\
  - plot and then prepare result of a single step by combining figures:\n\
        Lib_TwoDSubduction0_PlotCase plot_prepare_result_step -i ~/ASPECT_PROJECT/TwoDSubduction/test_peierls1/peierls\n\
\n\
  - generate animation: \n\
        Lib_TwoDSubduction0_PlotCase animate_case -i ~/ASPECT_PROJECT/TwoDSubduction/test_peierls1/peierls\n\
\n\
  - generate animation for cases in a directory: \n\
        Lib_TwoDSubduction0_PlotCase animate_case_in_dir -i `pwd`\n\
\n\
  - generate slab_morph.txt: \n\
    (what this does is looping throw all visualizing steps, so it takes time)\n\
        Lib_TwoDSubduction0_PlotCase morph_case -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear34/eba_low_tol_newton_shift_CFL0.8 \n\
\n\
  - plot trench movement: \n\
    (note you have to have a slab_morph.txt generated) \n\
        Lib_TwoDSubduction0_PlotCase plot_morph -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear34/eba_low_tol_newton_shift_CFL0.8 \n\
\n\
  - generate slab_morph.txt for cases in a directory: \n\
    Lib_TwoDSubduction0_PlotCase morph_case_in_dir -i ~/ASPECT_PROJECT/TwoDSubduction/EBA_CDPT1 \n\
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
    print("PlotCaseRun in TwoDSubduction0: operating")
    # get case parameters
    prm_path = os.path.join(case_path, 'output', 'original.prm')

    # plot visit
    Visit_Options = VISIT_OPTIONS(case_path)
    # provide steps to plot and interpret
    step = kwargs.get('step', None)
    if type(step) == int:
        Visit_Options.Interpret(steps=[step])  # only plot a single step
    else:
        Visit_Options.Interpret(last_step=3)  # by default, plot the last 3 steps
    odir = os.path.join(case_path, 'visit_scripts')
    if not os.path.isdir(odir):
        os.mkdir(odir)
    py_script = 'slab.py'
    ofile = os.path.join(odir, py_script)
    visit_script = os.path.join(ASPECT_LAB_DIR, 'visit_scripts', 'TwoDSubduction', py_script)
    visit_script_base = os.path.join(ASPECT_LAB_DIR, 'visit_scripts', 'base.py')
    Visit_Options.read_contents(visit_script_base, visit_script)  # combine these two scripts
    Visit_Options.substitute()
    ofile_path = Visit_Options.save(ofile, relative=True)
    print("Visualizing using visit")
    RunScripts(ofile_path)  # run scripts

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


def PrepareResultStep(case_path, step):
    '''
    Prepare results
    Inputs:
        case_path(str): path to the case
        step(int): step in computation
    '''
    # Generate json file
    odir = os.path.join(case_path, 'json_files')
    ofile = os.path.join(odir, 'figure_step.json')
    Prepare_Result = PREPARE_RESULT_OPTIONS(case_path)
    Prepare_Result.Interpret(step=step)
    # script to use
    if Prepare_Result.get_geometry() == "chunk":
        pr_script = os.path.join(ASPECT_LAB_DIR, "files", "TwoDSubduction", "figure_step_template.json")
    elif Prepare_Result.get_geometry() == "box":
        pr_script = os.path.join(ASPECT_LAB_DIR, "files", "TwoDSubduction", "figure_step_template_box.json")
    assert(os.path.isfile(pr_script))
    Prepare_Result.read_contents(pr_script)
    Prepare_Result.substitute()
    ofile_path = Prepare_Result.save(ofile, relative=True)
    PlotCombine.PrepareResults(ofile_path)  # call function with the generated json file


def PlotPrepareResultStep(case_path, step):
    '''
    First plot and then prepare results
    Inputs:
        case_path(str): path to the case
        step(int): step in computation
    '''
    # prepare result
    PlotCaseRun(case_path, step=step)
    PrepareResultStep(case_path, step)


def PrScriptToUse(case_path, template_for_chunk, template_for_box):
    '''
    Return model geometry
    Inputs:
    '''
    Prepare_Result = PREPARE_RESULT_OPTIONS(case_path)
    if Prepare_Result.get_geometry() == "chunk":
        pr_script = template_for_chunk
    elif Prepare_Result.get_geometry() == "box":
        pr_script = template_for_box
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
    
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    default_chunk = os.path.join(ASPECT_LAB_DIR, "files", "TwoDSubduction", "figure_step_template.json")
    default_box = os.path.join(ASPECT_LAB_DIR, "files", "TwoDSubduction", "figure_step_template_box.json")


    # commands
    # rearrange the time_range entry
    if arg.time == None or arg.time1 == None:
        time_range = None
    else:
        assert(type(arg.time) == float and type(arg.time1) == float)
        time_range = [arg.time, arg.time1]
    if _commend == 'plot_case':
        PlotCase.PlotCaseCombined([PlotCase.PlotCaseRun, PlotCaseRun], arg.inputs, time_range=time_range)
    elif _commend == 'plot_case_in_dir':
        PlotCase.PlotCaseCombinedDir([PlotCase.PlotCaseRun, PlotCaseRun], arg.inputs, time_range=time_range)
        pass
    elif _commend == 'prepare_result_step':
        pr_script = PrScriptToUse(arg.inputs, default_chunk, default_box)
        Plotter = PLOTTER(PREPARE_RESULT_OPTIONS, [PlotCase.PlotCaseRun, PlotCaseRun])
        Plotter.PrepareResultStep(arg.inputs, pr_script, arg.step)
    elif _commend == 'plot_prepare_result_step':
        pr_script = PrScriptToUse(arg.inputs, default_chunk, default_box)
        Plotter = PLOTTER(PREPARE_RESULT_OPTIONS, [PlotCase.PlotCaseRun, PlotCaseRun])
        Plotter.PlotPrepareResultStep(arg.inputs, pr_script, arg.step)
    elif _commend == 'animate_case':
        pr_script = PrScriptToUse(arg.inputs, default_chunk, default_box)
        Plotter = PLOTTER(PREPARE_RESULT_OPTIONS, [PlotCase.PlotCaseRun, PlotCaseRun])
        PlotCase.AnimateCaseResults(Plotter.PlotPrepareResultStep, arg.inputs, pr_script)
    elif _commend == 'animate_case_in_dir':
        pr_script = PrScriptToUse(arg.inputs, default_chunk, default_box)
        Plotter = PLOTTER(PREPARE_RESULT_OPTIONS, [PlotCase.PlotCaseRun, PlotCaseRun])
        PlotCase.AnimateCombinedDir(Plotter.PlotPrepareResultStep, arg.inputs, pr_script)
    elif _commend == 'morph_case':
        # slab_morphology, input is the case name
        vtk_and_slab_morph_case(arg.inputs, rewrite=arg.rewrite)
    elif _commend == 'morph_case_in_dir':
        # slab_morphology for cases in directory, input is the case name
        PlotCase.PlotCaseCombinedDir([vtk_and_slab_morph_case], arg.inputs, rewrite=arg.rewrite)
    elif _commend == 'plot_morph':
        # plot slab morphology
        SlabPlot = SLABPLOT('slab')
        prm_file = os.path.join(arg.inputs, 'output', 'original.prm')
        SlabPlot.ReadPrm(prm_file)
        SlabPlot.PlotMorph(arg.inputs)
    elif (_commend in ['-h', '--help']):
        # example:
        Usage()
    else:
        raise ValueError("Invalid command: use -h for help information")

# run script
if __name__ == '__main__':
    main()
