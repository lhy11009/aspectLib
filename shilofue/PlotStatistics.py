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

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')

sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities

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
        Inputs:
            step(int)
        get time to a value of model step
        '''
        col_t = self.header['Time']['col']
        col_step = self.header['Time_step_number']['col']
        times = self.data[:, col_t]
        steps = self.data[:, col_step]
        assert(len(times) == len(steps))
        time = 0.0
        found = False
        for i in range(len(steps)):  # search for step
            if step == int(steps[i]):
                time = times[i]
                found = True
        Utilities.my_assert(found, ValueError, "step %d is not a valid step" % step)
        return time

    def GetLastStep(self):
        '''
        get step and time of the last time step
        Return:
            last step, model time of the last step
        '''
        # get data
        col_t = self.header['Time']['col']
        col_step = self.header['Time_step_number']['col']
        times = self.data[:, col_t]
        steps = self.data[:, col_step]
        idx = np.argmax(steps)
        return int(steps[idx]), times[idx]

    # todo_combine
    def PlotNumberOfCells(self, **kwargs):
        '''
        plot the number of cells
        '''
        ax = kwargs.get('axis', None)
        label = kwargs.get('label', None)
        color = kwargs.get('color', None)
        if ax == None:
            raise ValueError("Not implemented")
        col_t = self.header['Time']['col']
        unit_t = self.header['Time']['unit']
        col_step = self.header['Time_step_number']['col']
        col_noc = self.header['Number_of_mesh_cells']['col']
        times = self.data[:, col_t]
        steps = self.data[:, col_step]
        nocs = self.data[:, col_noc]
        if self.UnitConvert is not None:
            to_myr = self.UnitConvert(unit_t, 'myr')
        else:
            raise ValueError("a UNITCONVERT class must be given")
        ax.plot(times * to_myr, nocs, label=label, color=color)
        ax.set_xlabel('Time (myr)')
        ax.set_ylabel('Number of mesh cells')
        pass
    
    # todo_combine
    def PlotNumberOfNonlinearIterations(self, **kwargs):
        '''
        plot the number of cells
        '''
        ax = kwargs.get('axis', None)
        label = kwargs.get('label', None)
        color = kwargs.get('color', None)
        if ax == None:
            raise ValueError("Not implemented")
        col_t = self.header['Time']['col']
        unit_t = self.header['Time']['unit']
        col_step = self.header['Time_step_number']['col']
        col_noni = self.header['Number_of_nonlinear_iterations']['col']
        times = self.data[:, col_t]
        steps = self.data[:, col_step]
        nonis = self.data[:, col_noni]
        if self.UnitConvert is not None:
            to_myr = self.UnitConvert(unit_t, 'myr')
        else:
            raise ValueError("a UNITCONVERT class must be given")
        ax.plot(times * to_myr, nonis, '.', label=label, color=color)
        ax.set_xlabel('Time (myr)')
        ax.set_ylabel('Number of nonlinear iterations')
        pass
    
    # todo_combine
    def PlotDegreeOfFreedom(self, **kwargs):
        '''
        plot the number of cells
        '''
        ax = kwargs.get('axis', None)
        label = kwargs.get('label', None)
        label_all = kwargs.get('label_all', False) # figure out labels
        if label_all:
            label_stokes = "(stokes)"
            label_temperature = "(temperature)"
            label_composition = "(composition)"
            if label == None:
                label_total = "(total)"
            else:
                label_total = label + " (total)"
        else:
            label_stokes = None
            label_temperature = None
            label_composition = None
            label_total = label
        color = kwargs.get('color', None)
        if ax == None:
            raise ValueError("Not implemented")
        col_t = self.header['Time']['col']
        unit_t = self.header['Time']['unit']
        col_step = self.header['Time_step_number']['col']
        col_dof_stokes = self.header['Number_of_Stokes_degrees_of_freedom']['col']
        col_dof_temperature = self.header['Number_of_temperature_degrees_of_freedom']['col']
        col_dof_composition = self.header['Number_of_degrees_of_freedom_for_all_compositions']['col']
        times = self.data[:, col_t]
        steps = self.data[:, col_step]
        dofs_stokes = self.data[:, col_dof_stokes]
        dofs_temperature = self.data[:, col_dof_temperature]
        dofs_composition = self.data[:, col_dof_composition]
        dofs_total = dofs_stokes + dofs_temperature + dofs_composition
        if self.UnitConvert is not None:
            to_myr = self.UnitConvert(unit_t, 'myr')
        else:
            raise ValueError("a UNITCONVERT class must be given")
        ax.plot(times * to_myr, dofs_total, '-', color=color, label=label_total)
        ax.plot(times * to_myr, dofs_stokes, ':', color=color, label=label_stokes)
        ax.plot(times * to_myr, dofs_temperature, '--', color=color, label=label_temperature)
        ax.plot(times * to_myr, dofs_composition, '-.', color=color, label=label_composition)
        ax.set_xlabel('Time (myr)')
        ax.set_ylabel('Number of degree of freedoms')
        pass
    
    # todo_combine
    def PlotTemperature(self, **kwargs):
        '''
        plot the number of cells
        '''
        ax = kwargs.get('axis', None)
        label = kwargs.get('label', None)
        label_all = kwargs.get('label_all', False) # figure out labels
        if label_all:
            label_average = "(average temperature)"
            label_minimum = "(minimum temperature)"
            label_maximum = "(maximum temperature)"
            if label == None:
                label_average = "(average temperature)"
            else:
                label_average = label + "(average temperature)"
        else:
            label_average = None
            label_minimum = None
            label_maximum = None
            label_average = label
        color = kwargs.get('color', None)
        if ax == None:
            raise ValueError("Not implemented")
        col_t = self.header['Time']['col']
        unit_t = self.header['Time']['unit']
        col_step = self.header['Time_step_number']['col']
        col_min_T = self.header['Minimal_temperature']['col']
        col_avg_T = self.header['Average_temperature']['col']
        col_max_T = self.header['Maximal_temperature']['col']
        times = self.data[:, col_t]
        steps = self.data[:, col_step]
        min_Ts = self.data[:, col_min_T]
        avg_Ts = self.data[:, col_avg_T]
        max_Ts = self.data[:, col_max_T]
        if self.UnitConvert is not None:
            to_myr = self.UnitConvert(unit_t, 'myr')
        else:
            raise ValueError("a UNITCONVERT class must be given")
        ax.plot(times * to_myr, avg_Ts, '-', color=color, label=label_average)
        ax.plot(times * to_myr, min_Ts, '-.', color=color, label=label_minimum)
        ax.plot(times * to_myr, max_Ts, '--', color=color, label=label_maximum)
        ax.set_xlabel('Time (myr)')
        ax.set_ylabel('Temperature (K)')
        pass
    
    # todo_combine
    def PlotVelocity(self, **kwargs):
        '''
        plot the velocity outputs
        '''
        ax = kwargs.get('axis', None)
        label = kwargs.get('label', None)
        label_all = kwargs.get('label_all', False) # figure out labels
        if label_all:
            label_maximum = "(maximum velocity)"
            if label == None:
                label_rms = "(rms velocity)"
            else:
                label_rms = label + "(rms velocity)"
        else:
            label_maximum = None
            label_rms = label
        color = kwargs.get('color', None)
        if ax == None:
            raise ValueError("Not implemented")
        col_t = self.header['Time']['col']
        unit_t = self.header['Time']['unit']
        col_step = self.header['Time_step_number']['col']
        col_rms_V = self.header['RMS_velocity']['col']
        unit_V = self.header['RMS_velocity']['unit']
        col_max_V = self.header['Max._velocity']['col']
        times = self.data[:, col_t]
        steps = self.data[:, col_step]
        rms_Vs = self.data[:, col_rms_V]
        max_Vs = self.data[:, col_max_V]
        if self.UnitConvert is not None:
            to_myr = self.UnitConvert(unit_t, 'myr')
        else:
            raise ValueError("a UNITCONVERT class must be given")
        ax.plot(times * to_myr, rms_Vs, '-', color=color, label=label_rms)
        ax.plot(times * to_myr, max_Vs, '--', color=color, label=label_maximum)
        ax.set_xlabel('Time (myr)')
        ax.set_ylabel('Velocity (%s)' % unit_V)
        pass

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
    UnitConvert = Utilities.UNITCONVERT()
    json_file = os.path.join(ASPECT_LAB_DIR, 'shilofue', 'json_files', 'post_process.json')
    with open(json_file, 'r') as fin:
        pdict = json.load(fin)
    plot_options = pdict.get('Statistics', {})
    Plotter = STATISTICS_PLOT('Statistics', unit_convert=UnitConvert, options=plot_options)
    fig_dir = os.path.dirname(fig_path)
    if not os.path.isdir(fig_dir):
        os.mkdir(fig_dir)
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