# -*- coding: utf-8 -*-
r"""plot resutls from depth-average file

This exports: 

  -

This depends on:

  -  

Examples of usage:

  - default usage:

        python -m shilofue.PlotDepthAverage plot_by_time -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear30/eba_intial_T/output/depth_average.txt 
    
  - plot case result:
        
        python -m shilofue.PlotDepthAverage plot_case -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear30/eba_re_mesh

        going to fix output name and time

descriptions
""" 
import numpy as np
import sys, os, argparse
import json 
# re
# import pathlib
# import subprocess
import numpy as np
# from matplotlib import cm
from matplotlib import pyplot as plt
import shilofue.Plot as Plot
from shilofue.Utilities import ReadHeader2, my_assert, UNITCONVERT

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')


class DEPTH_AVERAGE_PLOT(Plot.LINEARPLOT):
    '''
    Class for plotting depth average file.
    This is an inheritage of the LINEARPLOT class

    Attributes:
    Args:
    '''
    def __init__(self, _name, **kwargs):
        Plot.LINEARPLOT.__init__(self, _name, kwargs)  # call init from base function
        self.time_step_length = None
        # both these two arrays have the length of total time steps
        # the first records the time for each time step
        # the second points to the actual step within data
        self.time_step_times = None
        self.time_step_indexes = None

    def __call__(self, _filename, **kwargs):
        '''
        Read and plot
        Attributes:
            _filename(string):
                filename for data file
        Returns:
            _fileout(string or list):
                filename for output figure
        future:
            add in option for unit
        '''
        _fileout = kwargs.get('fileout', _filename + '.pdf')
        _time = kwargs.get('time', 'last')  # default is 'last' which means the last step
        self.ReadHeader(_filename)  # inteprate header information
        self.ReadData(_filename)  # read data
        self.ManageUnits()  # mange unit to output
        self.SplitTimeStep()  # split time step data
        # work out the name of output files
        _fname_base = _fileout.rpartition('.')[0]
        _fname_type = _fileout.rpartition('.')[2]
        _fname_list = []  # initialize a list for return
        if type(_time) in [float, int]:
            _time_list = [_time]
        elif type(_time) in [list, np.ndarray]:
            _time_list = _time
        else:
            raise TypeError('type of time needs to be in [float, int, list, np.ndarrayy]')
        for _t in _time_list:
            if type(_t) not in [float, int]:
                raise TypeError('type of values in time needs to be in [float, int, list, np.ndarrayy]')
            _data_list = self.ManageData(_t)  # manage output data
            _t_in_myr = _t * self.UnitConvert(self.header['time']['unit'], 'myr')
            _fname = "%s_t%.8e.%s" % (_fname_base, _t_in_myr, _fname_type)
            _figure_title = "Detph Average, t = %.10f myr" % _t_in_myr
            _fname = self.PlotCombine(_data_list, _fname, title=_figure_title)
            _fname_list.append(_fname)
        if len(_fname_list) == 0:
            # if there is only one name, just return this name
            return _fname_list[0]
        else:
            return _fname_list
        
    def ReadDataStep(self, _filename, **kwargs):
        '''
        Read data of a time step, currently only read the first time step.
        Attributes:
            _filename(string):
                filename for data file
        Returns:
            _datalist:
                a list of data for
        future:
            add in option for unit
        '''
        _fileout = kwargs.get('fileout', _filename + '.pdf')
        self.ReadHeader(_filename)  # inteprate header information
        self.ReadData(_filename)  # read data
        self.ManageUnits()  # mange unit to output
        self.SplitTimeStep()  # split time step data
        
        _t = 0.0  # t is the time of the 0th tep
        if type(_t) not in [float, int]:
            raise TypeError('type of values in time needs to be in [float, int, list, np.ndarrayy]')
        _data_list = self.ManageData(_t)  # manage output data

        data_type = kwargs.get('datatype', None)
        if data_type is None:
            return _data_list
        else:
            _data_list_o = []
            for _type in data_type:
                col = self.header[_type]['col']
                _data_list_o.append(_data_list[col])
            return _data_list_o


    def ReadHeader(self, _filename):
        '''
        Read header information from file.
        overload base function, use ReadHeader2
        function in utilities.py
        Args:
            _filename(str):
                filename for data file
        '''
        assert(os.access(_filename, os.R_OK))
        with open(_filename, 'r') as fin:
            _texts = fin.readlines()  # read the text of the file header
        self.header = ReadHeader2(_texts)

    def SplitTimeStep(self):
        '''
        split time steps, since the data is a big chunck
        '''
        _time_step_times = []  # initialize
        _time_step_indexes = []
        _col_time = self.header['time']['col']
        _col_depth = self.header['depth']['col']
        _times = self.data[:, _col_time]
        _depths = self.data[:, _col_depth]
        # get the length of a single time step
        for i in range(1, _depths.size):
            if _depths[i] < _depths[i-1]:
                self.time_step_length = i
                break
            elif i == _depths.size - 1:
                # as the exiting value from python is simply _depths.size - 1
                self.time_step_length = i + 1
        # make a ndarray of different value of time
        _step_times = [_times[_idx] for _idx in range(0, _times.size, self.time_step_length)]
        i = 0  # first sub list for first step
        _time_step_times.append(_step_times[0])
        _time_step_indexes.append([0])
        # loop to group data at the same step
        for j in range(1, len(_step_times)):
            _time = _step_times[j]
            if abs(_time - _step_times[j-1]) > 1e-16:
                _time_step_indexes.append([])
                _time_step_times.append(_time)
                i += 1
            _time_step_indexes[i].append(j)
        # both these two arrays have the length of total time steps
        # the first records the time for each time step
        # the second points to the actual step within data
        self.time_step_times = np.array(_time_step_times)
        self.time_step_indexes = _time_step_indexes
    
    def ManageData(self, _time):
        '''
        manage data, get new data for this class
        Returns:
            _data_list(list):
                list of data for ploting
            _time(float):
                time of plotting
        '''
        _data_list = []
        _time_step = np.argmin(abs(self.time_step_times - _time))  # time_step
        _index0 = self.time_step_indexes[_time_step][-1] * self.time_step_length
        if _time_step == len(self.time_step_times) - 1:
            # this is the last step
            _index1 = self.data.shape[0]
        else:
            _index1 = self.time_step_indexes[_time_step + 1][0] * self.time_step_length
        for i in range(self.data.shape[1]):
            _data_list.append(self.data[_index0 : _index1, i])
        # get the super adiabatic temperature
        _col_temperature = self.header['temperature']['col']
        _col_adiabatic_temperature = self.header['adiabatic_temperature']['col']
        _super_adiabatic_temperature = self.data[_index0 : _index1, _col_temperature] - self.data[_index0 : _index1, _col_adiabatic_temperature]
        _data_list.append(_super_adiabatic_temperature)
        self.header['super_adiabatic_temperature'] = {}
        self.header['super_adiabatic_temperature']['col'] = self.header['total_col']
        self.header['super_adiabatic_temperature']['unit'] = 'K'
        self.header['total_col'] += 1
        return _data_list
    
    def ManageUnits(self):
        '''
        manage units, get units for data.
        This is due to the bad form of the header of this file
        '''
        self.header['time']['unit'] = 'yr'
        self.header['depth']['unit'] = 'm'
        self.header['temperature']['unit'] = 'K'
        self.header['adiabatic_temperature']['unit'] = 'K'
        self.header['viscosity']['unit'] = 'Pa s'
        self.header['velocity_magnitude']['unit'] = 'm/yr'
        if self.dim == 2:
            self.header['vertical_heat_flux']['unit'] = 'mw/m'
        elif self.dim == 3:
            self.header['vertical_heat_flux']['unit'] = 'mw/m^2'


def PlotDaFigure(depth_average_path, fig_path_base):
    '''
    plot figure
    '''
    assert(os.access(depth_average_path, os.R_OK))
    # read that
    DepthAverage = DEPTH_AVERAGE_PLOT('DepthAverage')
    DepthAverage.ReadHeader(depth_average_path)
    DepthAverage.ReadData(depth_average_path)
    # manage data
    DepthAverage.SplitTimeStep()
    time_step = 0
    try:
        i0 = DepthAverage.time_step_indexes[time_step][-1] * DepthAverage.time_step_length
        if time_step == len(DepthAverage.time_step_times) - 1:
            # this is the last step
            i1 = DepthAverage.data.shape[0]
        else:
            i1 = DepthAverage.time_step_indexes[time_step + 1][0] * DepthAverage.time_step_length
    except IndexError:
        print("PlotDaFigure: File may not contain any depth average output, abort")
        return
    data = DepthAverage.data[i0:i1, :]
    # get depth
    col_depth = DepthAverage.header['depth']['col']
    depths = data[:, col_depth]
    # get pressure
    col_P = DepthAverage.header['adiabatic_pressure']['col']
    pressures = data[:, col_P]
    # get temperature
    col_T = DepthAverage.header['temperature']['col']
    temperatures = data[:, col_T]
    col_Tad = DepthAverage.header['adiabatic_temperature']['col']
    adiabat = data[:, col_Tad]
    # get viscosity
    col_eta = DepthAverage.header['viscosity']['col']
    eta = data[:, col_eta]
    # heat_flux
    col_hf = DepthAverage.header['vertical_heat_flux']['col']
    hf = data[:, col_hf]
    # heat_flux
    col_mf = DepthAverage.header['vertical_mass_flux']['col']
    mf = data[:, col_mf]

    # plot
    fig, axs = plt.subplots(2, 2, figsize=(10, 10))
    color = 'tab:blue'
    axs[0, 0].plot(pressures/1e9, depths/1e3, color=color, label='pressure')
    axs[0, 0].invert_yaxis()
    axs[0, 0].set_ylabel('Depth [km]') 
    axs[0, 0].set_xlabel('Pressure [GPa]', color=color) 
    # axs[0].invert_yaxis()
    # ax2: temperature
    color = 'tab:red'
    ax2 = axs[0, 0].twiny()
    ax2.plot(temperatures, depths/1e3, color=color, label='temperature')
    ax2.plot(adiabat, depths/1e3, '--', color=color, label='adiabat')
    ax2.set_xlabel('Temperature [K]', color=color) 
    ax2.legend()
    # second: viscosity
    axs[0, 1].semilogx(eta, depths/1e3, 'c', label='Viscosity')
    axs[0, 1].set_xlim([1e18,1e25])
    axs[0, 1].invert_yaxis()
    axs[0, 1].grid()
    axs[0, 1].set_ylabel('Depth [km]') 
    axs[0, 1].set_xlabel('Viscosity [Pa*s]')
    axs[0, 1].legend()
    # third: heat flux
    axs[1, 0].plot(hf, depths/1e3, color='r', label='vertical heat flux')
    axs[1, 0].invert_yaxis()
    axs[1, 0].grid()
    axs[1, 0].set_ylabel('Depth [km]') 
    axs[1, 0].set_xlabel('Heat Flux [W / m2]') 
    # fourth: mass flux
    yr_to_s = 365 * 24 * 3600
    axs[1, 1].plot(mf * yr_to_s, depths/1e3, color='b', label='vertical mass flux')
    axs[1, 1].invert_yaxis()
    axs[1, 1].grid()
    axs[1, 1].set_ylabel('Depth [km]') 
    axs[1, 1].set_xlabel('Mass Flux [kg / m2 / yr]') 
    # layout:w
    fig.tight_layout()
    # check directory
    fig_dir = os.path.dirname(fig_path_base)
    if not os.path.isdir(fig_dir):
        os.mkdir(fig_dir)
    fig_path_base0 = fig_path_base.rpartition('.')[0]
    fig_path_type = fig_path_base.rpartition('.')[2]
    fig_path = "%s.%s" % (fig_path_base0, fig_path_type)
    plt.savefig(fig_path)
    print("New figure: %s" % fig_path)


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
                        default='./DepthAverage.png',
                        help='output file(png)')
    parser.add_argument('-t', '--time', type=str,
                        default=0.0,
                        help='model time(s or yr)')
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    # commands
    if _commend == 'plot_by_time':
        # example:
        # use a json file
        json_file = os.path.join(shilofue_DIR, 'json_files', 'post_process.json')
        assert(os.access(json_file, os.R_OK))
        with open(json_file, 'r') as fin:
            json_options = json.load(fin)
        plot_options = json_options.get('DepthAverage', {})
    
        # Init the UnitConvert class
        UnitConvert = UNITCONVERT()

        # Initiate class for depth average plot
        DepthAverage = DEPTH_AVERAGE_PLOT('DepthAverage', unit_convert=UnitConvert, options=plot_options)

        # plot figure
        assert(os.access(arg.inputs, os.R_OK))
        if not os.path.isdir(os.path.dirname(arg.outputs)):
            os.mkdir(os.path.dirname(arg.outputs))
        DepthAverage(arg.inputs, fileout=arg.outputs, time=arg.time)

    # commands
    if _commend == 'plot_case':
        depth_average_path = os.path.join(arg.inputs, 'output', 'depth_average.txt')
        fig_path_base = os.path.join(arg.inputs, 'img', 'DepthAverage.png')
        PlotDaFigure(depth_average_path, fig_path_base)

# run script
if __name__ == '__main__':
    main()