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
    temp_path = os.path.join(ASPECT_LAB_DIR, 'run_time_output')
    os.system("awk -f %s/bash_scripts/awk_states/parse_block_output %s > %s" % (ASPECT_LAB_DIR, log_path, temp_path))

    # Init the UnitConvert class
    UnitConvert = UNITCONVERT()
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

# run script
if __name__ == '__main__':
    main()