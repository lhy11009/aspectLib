# -*- coding: utf-8 -*-
r"""(one line description)

This exports:

  -

This depends on:

  -

Examples of usage:

  - default usage:

        python -m

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
import shilofue.PlotCase as PlotCase
import shilofue.PlotDepthAverage as PlotDepthAverage
from matplotlib import pyplot as plt

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')

def Usage():
    print("\
(One liner description\n\
\n\
Examples of usage: \n\
\n\
  - default usage: \n\
\n\
    analyze result: \n\
        python -m shilofue.Others.latent_heat_benchmark plot_case -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/latent_heat_issue/latent-heat_pseudo_1d\n\
        ")

def AnalyzeResult(case_dir, **kwargs):
    '''
    Analyze case result
    Inputs:
        case_dir(str): case directory
    Returns:
        -
    '''
    PlotCase.PlotCaseRun(case_dir)
    depth_average_path = os.path.join(case_dir, 'output', 'depth_average.txt')
    if os.access(depth_average_path, os.R_OK):
        # work with depth_average file
        DepthAverage = PlotDepthAverage.DEPTH_AVERAGE_PLOT('DepthAverage')
        DepthAverage.ReadHeader(depth_average_path)
        DepthAverage.ReadData(depth_average_path)
        DepthAverage.SplitTimeStep()
        data_length = len(DepthAverage.time_step_times)
        # plot the initial time step
        fig_path_base = os.path.join(case_dir, 'img', 'DepthAverage.png')
        PlotDepthAverage.PlotDaFigure(depth_average_path, fig_path_base)
        # extract data by step
        query_depth = 5e5
        # query_temperatures = []
        band_widths = []
        max_Ts = []
        limit_Ts = []
        times = []
        max_time = kwargs.get('max_time', None)
        for i in range(data_length):
            time = DepthAverage.time_step_times[i]
            if time > max_time:
                break
            times.append(time)
            temperatures, _= DepthAverage.ExportDataByTime(time, ['depth', 'temperature'])
            max_T = np.max(temperatures[:, 1])  # max temperature
            max_Ts.append(max_T)
            min_T = np.min(temperatures[:, 1])  # max temperature
            limit_T = max_T - (max_T - min_T) / np.e
            limit_Ts.append(limit_T)
            mask_band = (temperatures[:, 1] > limit_T) # bandwidth of the "hot band"
            band_width = 0.0
            for i in range(1, temperatures.shape[0]):
                if (temperatures[i, 1] > limit_T):
                    band_width += temperatures[i, 0] - temperatures[i-1, 0]
            band_widths.append(band_width)
            # query_temperature = np.interp(query_depth, temperatures[:, 0], temperatures[:, 1])  # look for temperature in the middle
            # query_temperatures.append(query_temperature)
        fig, ax = plt.subplots()
        ax.plot(times, band_widths, label='bandwidth (m)')  # plot bandwidth
        ax.set_ylabel('Band Width (m)', color='tab:blue')
        ax.set_xlabel('Time (yr)')
        ax1 = ax.twinx()
        ax1.plot(times, max_Ts, 'r--', label='max T (K)')
        ax1.set_ylabel('max T (K)', color='tab:red')
        if max_time is None:
            fig_basename = 'HotBandWidth.png'
        else:
            # append max time information
            fig_basename = 'HotBandWidth_t%.4e.png' % max_time
        fig_path = os.path.join(case_dir, 'img', fig_basename)
        fig.savefig(fig_path)
        print('saved figure: ', fig_path)

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
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    # commands
    if (_commend in ['-h', '--help']):
        # example:
        Usage()
    elif _commend == 'plot_case':
        # example:
        AnalyzeResult(arg.inputs, max_time=2e9)
    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()
