# -*- coding: utf-8 -*-
r"""(one line description)

This exports: 

  -

This depends on:

  -  

Examples of usage:

  - default usage:

        python -m shilofue.PlotRunTime plot 
        -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear30/eba_re/output/log.txt 
        -o /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear30/eba_re/img/run_time.png

descriptions
""" 
import numpy as np
import sys, os, argparse
# import json, re
# import pathlib
import subprocess
import numpy as np
import shilofue.Plot as Plot
from shilofue.Utilities import UNITCONVERT
# from matplotlib import cm
from matplotlib import pyplot as plt

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')


def PlotFigure(log_path, fig_path):
    '''
    Read runtime info from log file and then plot
    Inputs:
        log_path(str) - path to log file
    Returns:
        RunTime.png
    '''
    # read log file
    temp_path = os.path.join(RESULT_DIR, 'run_time_output')
    os.system("awk -f %s/bash_scripts/awk_states/parse_block_output %s > %s" % (ASPECT_LAB_DIR, log_path, temp_path))

    # plot statistics ouput #####
    Plotter = Plot.LINEARPLOT('RunTime', {})
    Plotter.ReadHeader(temp_path)
    Plotter.ReadData(temp_path)
    col_time = Plotter.header['Time']['col']
    col_step = Plotter.header['Time_step_number']['col']
    col_wallclock = Plotter.header['Wall_Clock']['col']
    times = Plotter.data[:, col_time]
    steps = Plotter.data[:, col_step]
    wallclocks = Plotter.data[:, col_wallclock]
    # line 1: time
    fig, ax1 = plt.subplots(figsize=(5, 5)) 
    color = 'tab:blue'
    ax1.plot(steps, times, '-', color=color, label='Time') 
    ax1.set_ylabel('Time [yr]', color=color) 
    ax1.set_xlabel('Step') 
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.set_title('Run Time')
    # line 2: wall clock
    ax2 = ax1.twinx()
    color = 'tab:red'
    ax2.plot(steps, wallclocks, '-', color=color, label='Wall Clock [s]') 
    ax2.set_ylabel('Wall Clock [s]', color=color) 
    ax2.set_xlabel('Step') 
    ax2.tick_params(axis='y', labelcolor=color)
    # save figure
    fig.tight_layout()
    plt.savefig(fig_path) 


def PlotNewtonSolver(log_path, fig_path_base, **kwargs):
    '''
    Read runtime info from log file and then plot
    Inputs:
        log_path(str) - path to log file
    Returns:
        RunTime.png
    '''
    # read log file
    temp_path = os.path.join(RESULT_DIR, 'run_time_output_newton')
    os.system("awk -f %s/bash_scripts/awk_states/parse_block_newton %s > %s" % (ASPECT_LAB_DIR, log_path, temp_path))
    Plotter = Plot.LINEARPLOT('RunTime', {})
    Plotter.ReadHeader(temp_path)
    Plotter.ReadData(temp_path)

    step = kwargs.get('step', 0)  # step to plot
    col_step = Plotter.header['Time_step_number']['col']
    mask_step = (Plotter.data[:, col_step] == step)
    data = Plotter.data[mask_step, :]

    col_number_of_iteration = Plotter.header['Index_of_nonlinear_iteration']['col']
    number_of_iterations = data[:, col_number_of_iteration]
    col_residual = Plotter.header['Relative_nonlinear_residual']['col']
    residuals = data[:, col_residual]
    col_scaling_factor = Plotter.header["Newton_Derivative_Scaling_Factor"]['col']
    scaling_factors = data[:, col_scaling_factor]
    # line1: residual
    fig, ax1 = plt.subplots(figsize=(5, 5)) 
    color = 'tab:blue'
    ax1.semilogy(number_of_iterations, residuals, '-', color=color, label='Residuals') 
    ax1.set_ylabel('Relative non-linear residual', color=color) 
    ax1.set_xlabel('Number of non-linear iterations') 
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.set_title('Newton Solver Output For Step %d' % (step))
    # line 2: scaling factor
    ax2 = ax1.twinx()
    color = 'tab:red'
    ax2.plot(number_of_iterations, scaling_factors, '--', color=color, label='Scaling factors') 
    ax2.set_ylabel('Newton derivative scaling factor', color=color) 
    ax2.set_xlabel('Number of non-linear iterations') 
    ax2.tick_params(axis='y', labelcolor=color)

    # save figure 
    fig.tight_layout()
    fig_path_base0 = fig_path_base.rpartition('.')[0]
    fig_path_type = fig_path_base.rpartition('.')[2]
    fig_path = "%s_s%06d.%s" % (fig_path_base0, step, fig_path_type)
    plt.savefig(fig_path)
    return fig_path



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
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    # commands
    if _commend == 'plot':
        # example:
        PlotFigure(arg.inputs, arg.outputs)
    
    if _commend == 'plot_case':
        # example:
        log_file = os.path.join(arg.inputs, 'output', 'log.txt')
        assert(log_file)
        fig_path = os.path.join(arg.inputs, 'img', 'run_time.png')
        if not os.path.isdir(os.path.dirname(fig_path)):
            os.mkdir(os.path.dirname(fig_path))
        PlotFigure(log_file, fig_path)
    
    if _commend == "plot_newton_solver_step":
        # plot newton solver output
        o_path = PlotNewtonSolver(arg.inputs, arg.outputs) 

# run script
if __name__ == '__main__':
    main()