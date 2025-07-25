# -*- coding: utf-8 -*-
r"""(one line description)

This exports:

  -

This depends on:

  -

Examples of usage:

  - default usage: plot case running results

         plot_case -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear32/eba1_MRf12_iter20

descriptions
"""
import numpy as np
import sys, os, argparse, re
import json
# import pathlib
# import subprocess
import numpy as np
# from matplotlib import cm
from matplotlib import pyplot as plt
from matplotlib import gridspec
from shilofue.ParsePrm import ReadPrmFile
import shilofue.PlotVisit as PlotVisit
import shilofue.PlotRunTime as PlotRunTime
import shilofue.PlotStatistics as PlotStatistics
import shilofue.PlotCombine as PlotCombine
import imageio

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')

sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities

HaMaGeoLib_DIR = "/home/lochy/ASPECT_PROJECT/HaMaGeoLib"
if os.path.abspath(HaMaGeoLib_DIR) not in sys.path:
    sys.path.append(os.path.abspath(HaMaGeoLib_DIR))
from hamageolib.research.haoyuan_2d_subduction.legacy_tools import PlotNewtonSolverHistory

def Usage():
    print("\
This scripts generate plots for a single case in TwoDSubduction project\n\
\n\
Examples of usage: \n\
\n\
  - default usage: plot case running results, -t option deals with a time range, default is a whole range\n\
\n\
        Lib_plot_case plot_case -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear32/eba1_MRf12_iter20\
          -t 0.0 -t1 0.5e6\
        ")


class PLOTTER():
    '''
    A class for preparing results
    '''
    def __init__(self, module, PlotCaseRuns, **kwargs):
        '''
        Initiation
        module (class) - class to use for generating results at each step
        PlotCaseRuns (list of functions) - a list of functions for generating result at a given step.
        kwargs(dict):
        '''
        self.module = module
        self.PlotCaseRuns = PlotCaseRuns
        pass

    def GenerateJson(self, case_path, pr_script, step):
        '''
        Interpret the required json file for plotting
        Inputs:
            case_path(str): path to the case
            step(int): step in computation
        '''
        # Generate json file
        odir = os.path.join(case_path, 'json_files')
        ofile = os.path.join(odir, 'figure_step.json')
        Prepare_Result = self.module(case_path)
        Prepare_Result.Interpret(step=step)
        Utilities.my_assert(os.path.isfile(pr_script), AssertionError, "%s doesn't exist." % pr_script) # check script to use
        Prepare_Result.read_contents(pr_script)
        Prepare_Result.substitute()
        ofile_path = Prepare_Result.save(ofile, relative=True)
        return ofile_path
    
    def PrepareResultStep(self, case_path, pr_script, step, **kwargs):
        '''
        Prepare results
        Inputs:
            case_path(str): path to the case
            step(int): step in computation
            kwargs(dict):
                update: update result if the final result is presented
        '''
        update = kwargs.get('update', True)
        generate_json = kwargs.get('generate_json', True)
        # Generate json file for combining result
        odir = os.path.join(case_path, 'json_files')
        ofile = os.path.join(odir, 'figure_step.json')
        if generate_json:
            Prepare_Result = self.module(case_path)
            Prepare_Result.Interpret(step=step)
            assert(os.path.isfile(pr_script)) # check script to use
            Prepare_Result.read_contents(pr_script)
            Prepare_Result.substitute()
            ofile_path = Prepare_Result.save(ofile, relative=True)  # save json file
            print("PrepareResultStep: json file generated %s" % ofile_path)
        else:
            # skip this step
            ofile_path = ofile
            assert(os.path.isfile(ofile_path))
            print("PrepareResultStep: use previous json file %s" % ofile_path)
        with open(ofile_path, 'r') as fin:
            contents = json.load(fin)
        file_to_expect = contents['figures'][-1]["save path"]
        if update or not os.path.isfile(file_to_expect): 
            PlotCombine.PrepareResults(ofile_path)  # call function with the generated json file
        else:
            print("%s: result at step %d is found (%s), skip" % (Utilities.func_name(), step, file_to_expect))

    def PlotPrepareResultStep(self, case_path, pr_script, step, **kwargs):
        '''
        First plot and then prepare results
        Inputs:
            case_path(str): path to the case
            step(int): step in computation
            kwargs(dict):
                update: update result if the final result is presented
        Returns:
            file_to_expect(str): file to generate, this is the file at the bottom of the list
        '''
        update = kwargs.get('update', True)  # if update on results
        ofile_path = self.GenerateJson(case_path, pr_script, step) # generate json file first
        with open(ofile_path, 'r') as fin:
            contents = json.load(fin)
        file_to_expect = contents['figures'][-1]["save path"]  # find the file we expect
        if update or not os.path.isfile(file_to_expect): 
            # plot results, here step means visualization step
            for PlotCaseRun0 in self.PlotCaseRuns:
                # we can take multiple functions from different sources
                PlotCaseRun0(case_path, step=step)
            # combine results, as the json file is generated, skip here.
            self.PrepareResultStep(case_path, step, pr_script, update=update, generate_json=False)
        else:
            print("%s: result at step %d is found (%s), skip" % (Utilities.func_name(), step, file_to_expect))
        return file_to_expect

# todo_legacy: move plotting functions to new package
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
    print("PlotCaseRun: operating")
    # get case parameters
    prm_path = os.path.join(case_path, 'output', 'original.prm')
    if not os.path.isfile(prm_path):
        print('original.prm is missing, it\'s likely that there is no outputs. Skip.')
        return 1
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
    time_range = kwargs.get('time_range', None)
    fig = plt.figure(tight_layout=True, figsize=(6, 18))
    gs = gridspec.GridSpec(3,1)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[2, 0])
    if time_range is not None:
        try:
            fig_output_path, step_range = PlotRunTime.PlotFigure(log_file, fig_path, axis=ax1, ax2=ax2, ax3=ax3,\
                                                             savefig=False, fix_restart=True, time_range=time_range) 
            fig.savefig(fig_path)
            print("New figure: %s" % fig_path)
        except Exception:
            print("PlotRunTime is unsuccessful\n")
    else:
        try:
            fig_output_path = PlotRunTime.PlotFigure(log_file, fig_path, axis=ax1, ax2=ax2, ax3=ax3,\
                                                             savefig=False, fix_restart=True, save_temp_file_local=True) 
        
            fig.savefig(fig_path)
            step_range = None
            print("New figure: %s" % fig_path)
        except Exception:
            step_range = None
            print("PlotRunTime is unsuccessful\n")

    # Newton history
    # determine whether newton is used
    print('solver_scheme, ', solver_scheme) # debug
    fig_path = os.path.join(case_path, 'img', 'newton_solver_history.png')
    # match object in solver scheme, so we only plot for these options
    match_obj = re.search('Newton', solver_scheme)
    match_obj1 = re.search('iterated defect correction Stokes', solver_scheme)
    if match_obj or match_obj1:
        print("Plotting newton solver history")
        log_path = os.path.join(RESULT_DIR, 'run_time_output_newton')
        # PlotNewtonSolverHistory(log_path, fig_path, step_range=step_range)
    else:
        print("Skipping newton solver history")
    return 0


def PlotCaseCombined(modules, inputs, **kwargs):
    '''
    combine several modules in plotting
    inputs:
        modules (list of function)
        inputs (str): path to a case
        **kwargs:
            time_range ([None, None] or [double, double]): range of time to plot
    '''
    assert(type(modules) == list)
    # check time_range
    time_range = kwargs.get('time_range', None)
    if time_range is not None:
        # Check there is a meaningful time range
        assert(type(time_range) == list and\
         type(time_range[0]) == float and\
         type(time_range[1]) == float)
    # call functions to plot
    for module in modules:
        module(inputs, **kwargs)


def PlotCaseCombinedDir(modules, pdir, **kwargs):
    '''
    combine several modules in plotting
    inputs:
        modules (list of function)
        pdir (str): path to a directory
        **kwargs:
            time_range ([None, None] or [double, double]): range of time to plot
    '''
    assert(os.path.isdir(pdir))
    for _name in os.listdir(pdir):
        _path = os.path.join(pdir, _name)
        if os.path.isdir(_path):
            case_prm = os.path.join(_path, 'case.prm')
            if os.path.isfile(case_prm):
                print("\nFind case %s" % _path)
                PlotCaseCombined(modules, _path, **kwargs)


def AnimateCaseResults(PrepareS, case_path, pr_script, **kwargs):
    '''
    create animation
    this function will be called upon with the "animate_case" option
    of each project.
    Note that his would plot, at each step, the Statistic & Newton solver files again
    by default, so it's better to not do that.
    Inputs:
        PrepareS: a function to plot case result for a single step
            I will pass "Plotter.PlotPrepareResultStep" for each
            project to this function.
        case_path(str): path to the case
        kwargs(dict):
            step_range(list of 2): a list of steps to animate
            name(str): name of the animation
            duration: time for each frame
    '''
    # initiation, note that we use the interfaces in VISIT_OPTIONS
    # class to figure out the steps to plot
    name = kwargs.get('name', 'ani')  # name of the animation
    duration = kwargs.get('duration', 0.2)  # time for each frame
    time_interval = kwargs.get('time_interval', None)
    VisitOptions = PlotVisit.VISIT_OPTIONS(case_path)
    VisitOptions.Interpret(time_interval=time_interval)
    steps = VisitOptions.options['GRAPHICAL_STEPS']
    # plot results
    filenames = []
    for step in steps:
        # prepare results, if the figure is already generated, skip
        filename = PrepareS(case_path, pr_script, step, update=False)
        filenames.append(filename)
    # use the imageio module to generate gif
    o_dir = os.path.dirname(filenames[-1])  # just use this directory as output
    o_path = os.path.join(o_dir, "%s.gif" % name)
    print("%s: saving file %s" % (Utilities.func_name(), o_path))
    frames = []
    for filename in filenames:
        frames.append(imageio.imread(filename))
    imageio.mimsave(o_path, frames, 'GIF', duration=duration)


def AnimateCombinedDir(PrepareS, dir, pr_script, **kwargs):
    '''
    combine animation for cases in a directory
    inputs:
        PrepareS
        case_path(str): path to the case
        kwargs(dict):
            step_range(list of 2): a list of steps to animate
            name(str): name of the animation
    '''
    assert(os.path.isdir(dir))
    name = kwargs.get('name', 'ani')  # name of the animation
    for subdir, dirs, _ in os.walk(dir):
        for _dir in dirs:
            _path = os.path.join(subdir, _dir)
            case_prm = os.path.join(_path, 'case.prm')
            if os.path.isfile(case_prm):
                print("\nFind case %s" % _path)
                try:
                    step_range = kwargs['step_range']
                except KeyError:
                    AnimateCaseResults(PrepareS, _path, pr_script, name=name)
                else:
                    AnimateCaseResults(PrepareS, _path, pr_script, step_range=step_range, name=name)



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
