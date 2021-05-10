# -*- coding: utf-8 -*-
r"""Plot statistics output

This exports: 

  -

This depends on:

  -  

Examples of usage:

  - default usage: plot statistic output

        python -m shilofue.PlotStatistics plot
        -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear30/eba_re_mesh/output/statistics
        -o /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear30/eba_re_mesh/img/Statistics.png

  - plot a case
        
        python -m shilofue.PlotStatistics plot_case -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear30/eba_re_mesh

descriptions
""" 
import numpy as np
import sys, os, argparse
import json
# import pathlib
# import subprocess
import numpy as np
# from matplotlib import cm
from matplotlib import pyplot as plt
from shilofue.Plot import LINEARPLOT
from shilofue.Utilities import UNITCONVERT

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')


class STATISTICS_PLOT(LINEARPLOT):
    '''
    Class for plotting depth average file.
    This is an inheritage of the LINEARPLOT class

    Attributes:
    Args:
    '''
    def __init__(self, _name, **kwargs):
        LINEARPLOT.__init__(self, _name, kwargs)  # call init from base function
    
    def GetStep(self, time):
        '''
        Inputs:
            time(double)
        get step corresponding to a value of model time
        '''
        # get data
        col_t = self.header['Time']['col']
        col_step = self.header['Time_step_number']['col']
        times = self.data[:, col_t]
        steps = self.data[:, col_step]

        # get step
        idx = np.argmin(abs(times - time))
        step = int(steps[idx])
        return step
    
    def GetTime(self, step):
        '''
        future
        Inputs:
            step(int)
        get time to a value of model step
        '''
        time = 0.0
        return time


def PlotFigure(file_path, fig_path):
    '''
    descriptions
    Inputs:
        - file_path(str): path of a statistic file of aspect
        - figure_path(str): path of the output figure
    Returns:
        -
    '''
    # Init the UnitConvert class
    UnitConvert = UNITCONVERT()
    json_file = os.path.join(ASPECT_LAB_DIR, 'shilofue', 'json_files', 'post_process.json')
    with open(json_file, 'r') as fin:
        pdict = json.load(fin)
    plot_options = pdict.get('Statistics', {})
    Plotter = STATISTICS_PLOT('Statistics', unit_convert=UnitConvert, options=plot_options)
    fig_generated_path = Plotter(file_path, fileout=fig_path)  # plot figure
    print("New figure: %s" % fig_generated_path)
    return fig_generated_path
    
    pass


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
    
    elif _commend == 'plot_case':
        # example:
        statistic_file = os.path.join(arg.inputs, 'output', 'statistics')
        fig_path = os.path.join(arg.inputs, 'img', 'Statistic.png')
        PlotFigure(statistic_file, fig_path)

# run script
if __name__ == '__main__':
    main()