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
from shilofue.TwoDSubduction0.PlotSlab import vtk_and_slab_morph_case
from shilofue.TwoDSubduction0.VtkPp import SlabMorphologyCase, PlotWedgeTCase, WedgeTCase, SLABPLOT, PlotTrenchThermalState
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
\n\
  - combine a few steps (visualization, morphology, wedge T, etc): \n\
    Lib_TwoDSubduction0_PlotCase plot_case_composite -i /mnt/lochy0/ASPECT_DATA/TwoDSubduction/EBA_CDPT4/eba_cdpt_SA80.0_OA40.0_CpEcl \n\
        ")

def PlotCaseRun(case_path, **kwargs):
    '''
    Plot case run result
    Inputs:
        case_path(str): path to the case
        kwargs:
            time_range
            step(int): if this is given as an int, only plot this step
            visualization (str): visualization software, visit or paraview.
            last_step: number of last steps to plot
    Returns:
        -
    '''
    run_visual = kwargs.get('run_visual', 0)
    step = kwargs.get('step', None)
    time_interval = kwargs.get('time_interval', None)
    visualization = kwargs.get('visualization', 'visit')
    plot_axis = kwargs.get('plot_axis', False)
    last_step = kwargs.get('last_step', 3)
    max_velocity = kwargs.get('max_velocity', -1.0)
    plot_types = kwargs.get("plot_types", ["upper_mantle"])
    rotation_plus = kwargs.get("rotation_plus", 0.0)
    assert(visualization in ["paraview", "visit", "pygmt"])
    print("PlotCaseRun in TwoDSubduction0: operating")
    # get case parameters
    prm_path = os.path.join(case_path, 'output', 'original.prm')

    # steps to plot: here I use the keys in kwargs to allow different
    # options: by steps, a single step, or the last step
    if type(step) == int:
        kwargs["steps"] = [step]
    elif type(step) == list:
        kwargs["steps"] = step
    else:
        kwargs["last_step"] = last_step

    # Inititiate the class and intepret the options
    # Note that all the options defined by kwargs is passed to the interpret function
    Visit_Options = VISIT_OPTIONS(case_path)
    Visit_Options.Interpret(**kwargs)

    # generate scripts base on the method of plotting
    if visualization == 'visit':
        odir = os.path.join(case_path, 'visit_scripts')
        if not os.path.isdir(odir):
            os.mkdir(odir)
        print("Generating visit scripts")
        py_script = 'slab.py'
        ofile = os.path.join(odir, py_script)
        visit_script = os.path.join(ASPECT_LAB_DIR, 'visit_scripts', 'TwoDSubduction', py_script)
        visit_script_base = os.path.join(ASPECT_LAB_DIR, 'visit_scripts', 'base.py')
        Visit_Options.read_contents(visit_script_base, visit_script)  # combine these two scripts
        Visit_Options.substitute()
    elif visualization == 'paraview':
        odir = os.path.join(case_path, 'paraview_scripts')
        if not os.path.isdir(odir):
            os.mkdir(odir)
        print("Generating paraview scripts")
        py_script = 'slab.py'
        ofile = os.path.join(odir, py_script)
        paraview_script = os.path.join(ASPECT_LAB_DIR, 'paraview_scripts', 'TwoDSubduction', py_script)
        paraview_script_base = os.path.join(ASPECT_LAB_DIR, 'paraview_scripts', 'base.py')
        Visit_Options.read_contents(paraview_script_base, paraview_script)  # combine these two scripts
        Visit_Options.substitute()
    elif visualization == 'pygmt':
        odir = os.path.join(case_path, 'pygmt_scripts')
        if not os.path.isdir(odir):
            os.mkdir(odir)
        print("Generating pygmt scripts")
        py_script = 'make_lateral_flow_fig.py'
        ofile = os.path.join(odir, py_script)
        pygmt_script = os.path.join(ASPECT_LAB_DIR, 'pygmt_scripts', 'TwoDSubduction', py_script)
        pygmt_script_base = os.path.join(ASPECT_LAB_DIR, 'pygmt_scripts', 'aspect_plotting_util.py')
        Visit_Options.read_contents(pygmt_script_base, pygmt_script)  # combine these two scripts
        Visit_Options.substitute()

    ofile_path = Visit_Options.save(ofile, relative=True)
    if run_visual == 1:
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
    parser.add_argument('-ti', '--time_interval', type=float,
                        default=None,
                        help='Time interval, affecting the time steps to visualize')
    parser.add_argument('-s', '--step', type=int,
                        default=0,
                        help='step')
    parser.add_argument('-r', '--rewrite', type=int,
                        default=0,
                        help='If rewrite previous result')
    parser.add_argument('-rv', '--run_visualization', type=int,
                        default=0,
                        help='if visualization programs run with the script we generate to get figures')
    
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    # this two are the json files for the order of options to do in imageio
    # following the options defined in this two files, the results would be a combination of result for
    # one single computation step.
    # default_chunk = os.path.join(ASPECT_LAB_DIR, "files", "TwoDSubduction", "figure_step_template.json")
    default_chunk = os.path.join(ASPECT_LAB_DIR, "files", "TwoDSubduction", "figure_step_template_chunk_01032024.json")
    default_box = os.path.join(ASPECT_LAB_DIR, "files", "TwoDSubduction", "figure_step_template_box.json")


    # commands
    # rearrange the time_range entry
    if arg.time == None or arg.time1 == None:
        time_range = None
    else:
        assert(type(arg.time) == float and type(arg.time1) == float)
        time_range = [arg.time, arg.time1]
    if _commend == 'plot_case':
        PlotCase.PlotCaseCombined([PlotCase.PlotCaseRun, PlotCaseRun], arg.inputs, time_range=time_range, run_visual=arg.run_visualization,\
        time_interval=arg.time_interval, visualization="paraview", last_step=1)
        # PlotCase.PlotCaseCombined([PlotCaseRun], arg.inputs, time_range=time_range, run_visual=arg.run_visualization, time_interval=arg.time_interval)
    elif _commend == 'plot_case_in_dir':
        PlotCase.PlotCaseCombinedDir([PlotCase.PlotCaseRun, PlotCaseRun], arg.inputs, time_range=time_range, run_visual=arg.run_visualization, time_interval=arg.time_interval)
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
        Plotter = PLOTTER(PREPARE_RESULT_OPTIONS, []) # note we don't want to replot things here
        PlotCase.AnimateCaseResults(Plotter.PlotPrepareResultStep, arg.inputs, pr_script, time_interval=arg.time_interval, duration=0.2)
    elif _commend == 'animate_case_in_dir':
        pr_script = PrScriptToUse(arg.inputs, default_chunk, default_box)
        Plotter = PLOTTER(PREPARE_RESULT_OPTIONS, [PlotCaseRun])
        PlotCase.AnimateCombinedDir(Plotter.PlotPrepareResultStep, arg.inputs, pr_script)
    elif _commend == 'morph_case':
        # slab_morphology, input is the case name
        SlabPlot = SLABPLOT('slab')
        PlotCase.PlotCaseCombined([vtk_and_slab_morph_case, SlabPlot.PlotMorph], arg.inputs, rewrite=arg.rewrite)
    elif _commend == 'morph_case_in_dir':
        # slab_morphology for cases in directory, input is the case name
        SlabPlot = SLABPLOT('slab')
        PlotCase.PlotCaseCombinedDir([vtk_and_slab_morph_case, SlabPlot.PlotMorph], arg.inputs, rewrite=arg.rewrite)
    elif _commend == 'plot_morph':
        # plot slab morphology
        SlabPlot = SLABPLOT('slab')
        SlabPlot.PlotMorph(arg.inputs)
    elif _commend == 'plot_case_composite':
        # todo_composite
        # plot case composite
        # plot the linear plot and visualization
        time_interval_visual = 0.5e6
        PlotCase.PlotCaseCombined([PlotCase.PlotCaseRun, PlotCaseRun],\
                                  arg.inputs,\
                                  time_range=time_range, \
                                  run_visual=arg.run_visualization,\
                                  time_interval=time_interval_visual)
        # plot the slab morphology
        time_interval_morph = 0.5e6 
        SlabMorphologyCase(arg.inputs, rewrite=1, time_interval=time_interval_morph)
        SlabPlot = SLABPLOT('slab')
        SlabPlot.PlotMorph(arg.inputs)
        # plot the wedge T
        time_interval_wedgt_T = 0.5e6 
        WedgeTCase(arg.inputs, time_interval=time_interval_wedgt_T)
        PlotWedgeTCase(arg.inputs, time_interval=time_interval_wedgt_T)
        # plot the thermal state
        # this has to use the same interval as the plot morphology 
        PlotTrenchThermalState(arg.inputs, time_interval=time_interval_morph, silent=True)
    elif (_commend in ['-h', '--help']):
        # example:
        Usage()
    else:
        raise ValueError("Invalid command: use -h for help information")

# run script
if __name__ == '__main__':
    main()
