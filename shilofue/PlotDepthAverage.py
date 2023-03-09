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
from scipy.interpolate import interp1d
# from matplotlib import cm
from matplotlib import pyplot as plt
import shilofue.Plot as Plot
from shilofue.Utilities import ReadHeader2, my_assert, UNITCONVERT

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')
# import utilities in subdirectiory
sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities


def Usage():
    print("\
Plot and process depth_average output from aspect\n\
\n\
Examples of usage: \n\
  - plot by giving a time: \n\
        python -m shilofue.PlotDepthAverage plot_by_time\n\
        -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear30/eba_intial_T/output/depth_average.txt\n\
\n\
  - plot figure: \n\
\n\
        python -m shilofue.PlotDepthAverage plot_case -i ~/ASPECT_PROJECT/LatentHeatBK/LatentHeat3/lh_test -t 1e9 \
\n\
  - export data: \n\
\n\
        python -m shilofue.PlotDepthAverage export_case \n\
        -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear32/eba1_MRf12_iter20_DET660 \n\
        -o .test/depth_average_export\
        ")


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
                filename for data file~/ASPECT_PROJECT/TwoDSubduction/non_linear32/eba1_MRf12_iter20_DET660/output$
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
    
    def ExportDataByTime(self, time, names):
        '''
        Export data as ndarray by time and names
        '''
        assert(type(names)==list)
        time_step = np.argmin(abs(self.time_step_times - time))  # time_step
        i0 = self.time_step_indexes[time_step][-1] * self.time_step_length
        if time_step == len(self.time_step_times) - 1:
            # this is the last step
            i1 = self.data.shape[0]
        else:
            i1 = self.time_step_indexes[time_step + 1][0] * self.time_step_length
        odata = self.export("", names, rows=[i for i in range(i0, i1)], include_size=True, data_only=True)
        return odata, self.time_step_times[time_step]


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


    def Import(self, _filename):
        '''
        Combine a few functions to read data, header, as well as split
        data to steps
        '''
        self.ReadHeader(_filename)
        self.ReadData(_filename)
        self.SplitTimeStep()


    def SplitTimeStep(self):
        '''
        split time steps, since the data is a big chunck
        '''
        time_step_times = []  # initialize
        time_step_indexes = []
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
        time_step_times.append(_step_times[0])
        time_step_indexes.append([0])
        # loop to group data at the same step
        for j in range(1, len(_step_times)):
            _time = _step_times[j]
            if abs(_time - _step_times[j-1]) > 1e-16:
                time_step_indexes.append([])
                time_step_times.append(_time)
                i += 1
            time_step_indexes[i].append(j)
        # both these two arrays have the length of total time steps
        # the first records the time for each time step
        # the second points to the actual step within data
        self.time_step_times = np.array(time_step_times)
        self.time_step_indexes = time_step_indexes
    
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

    def GetInterpolateFunc(self, time, field_name):
        names = ["depth", field_name]
        odata, _ = self.ExportDataByTime(time, names)
        _func = interp1d(odata[:, 0], odata[:, 1], assume_sorted=True, fill_value="extrapolate")
        return _func
    
    def GetIntegrateArray(self, time, field_name, dim, geometry, geometry_length, geometry_width = None, **kwargs):
        '''
        Returns:
            integretions - an array of added volume * field
            segmentations - an array of volume * field
        Note the data at depth 0 has to be fixed in some way (as depth 0 is not in the depth_average files)
        '''
        # initiating
        Utilities.my_assert(geometry in ["cartesian", "spherical"], ValueError,\
            "geometry must be either cartesian or spherical")
        Ro = kwargs.get('Ro', 6371e3)
        # get raw data
        names = ["depth", field_name]
        odata, _ = self.ExportDataByTime(time, names)
        depths = odata[:, 0]
        vals = odata[:, 1]
        # compute integretion
        segmentations = np.zeros(self.time_step_length)
        integretions = np.zeros(self.time_step_length)
        if geometry == "cartesian" and dim == 2:
            # compute the value for the shallowes depth,
            # note this only takes the first value in
            # the data
            integretions[0] = depths[0] * geometry_length * vals[0]
            segmentations[0] = integretions[0]
        elif geometry == "spherical" and dim == 2:
            integretions[0] = geometry_length/2.0 * (Ro**2.0 - (Ro - depths[0])**2.0) * vals[0]
            segmentations[0] = integretions[0]
        else:
            raise NotImplementedError()
        for i in range(1, self.time_step_length):
            if geometry == "cartesian" and dim == 2:
                volume = (depths[i] - depths[i-1]) * geometry_length
            elif geometry == "spherical" and dim == 2:
                r_last = Ro - depths[i-1]
                r = Ro - depths[i]
                volume = geometry_length/2.0 * (r_last**2.0 - r**2.0) # geometry_length in theta
            else:
                raise NotImplementedError()
            # compute average between the current point and the last point
            integretions[i] = integretions[i-1] + (vals[i-1] + vals[i])/2.0 * volume 
            segmentations[i] = (vals[i-1] + vals[i])/2.0 * volume 
        return integretions, segmentations

def PlotDaFigure(depth_average_path, fig_path_base, **kwargs):
    '''
    plot figure
    Inputs:
        kwargs:
            time_step - time_step to plot the figure, default is 0
    '''
    time = kwargs.get('time', 0.0)
    assert(os.access(depth_average_path, os.R_OK))
    # read that
    DepthAverage = DEPTH_AVERAGE_PLOT('DepthAverage')
    DepthAverage.ReadHeader(depth_average_path)
    DepthAverage.ReadData(depth_average_path)
    
    # manage data
    DepthAverage.SplitTimeStep()
    names = ['depth', 'adiabatic_pressure', 'temperature', 'adiabatic_temperature', 'viscosity', 'vertical_heat_flux', 'vertical_mass_flux', 'adiabatic_density']
    # plot the outputs of the log value of the viscosity if data presents
    plot_log_viscosity = DepthAverage.Has("log_viscosity")
    if plot_log_viscosity:
        names.append("log_viscosity")
    data, exact_time = DepthAverage.ExportDataByTime(time, names)

    # get depth
    depths = data[:, 0]
    # get pressure
    pressures = data[:, 1]
    # get temperature
    temperatures = data[:, 2]
    adiabat = data[:, 3]
    # get viscosity
    eta = data[:, 4]
    # heat_flux
    hf = data[:, 5]
    # heat_flux
    mf = data[:, 6]
    # density
    densities = data[:, 7]
    if plot_log_viscosity:
        log_eta = data[:, 8]
    else:
        log_eta = None

    # plot
    fig, axs = plt.subplots(3, 2, figsize=(10, 10))
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
    if plot_log_viscosity:
        axs[0, 1].semilogx(10**log_eta, depths/1e3, 'c--', label='Viscosity (from log)')
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
    # fifth: density
    axs[2, 0].plot(densities, depths/1e3, color='b', label='adiabatic density')
    axs[2, 0].invert_yaxis()
    axs[2, 0].grid()
    axs[2, 0].set_ylabel('Depth [km]') 
    axs[2, 0].set_xlabel('Density [kg / m3]') 
    # title
    fig.suptitle("Time = %.4e" % exact_time)
    # layout:w
    fig.tight_layout()
    # check directory
    fig_dir = os.path.dirname(fig_path_base)
    if not os.path.isdir(fig_dir):
        os.mkdir(fig_dir)
    fig_path_base0 = fig_path_base.rpartition('.')[0]
    fig_path_type = fig_path_base.rpartition('.')[2]
    fig_path = "%s_t%.4e.%s" % (fig_path_base0, exact_time, fig_path_type)
    plt.savefig(fig_path)
    print("New figure: %s" % fig_path)


def ExportData(depth_average_path, output_dir, **kwargs):
    '''
    Export data of a step to separate file
    Inputs:
        kwargs:
            time_step - time_step to plot the figure, default is 0
            fix_time_step - fix time step if it's beyong the limit, default is False.
    Returns:
        odata (ndarray): selected data for outputing
        output_path (str): path of file generated
    '''
    time_step = kwargs.get('time_step', 0)
    fix_time_step = kwargs.get('fix_time_step', False)
    assert(os.access(depth_average_path, os.R_OK))
    # read that
    DepthAverage = DEPTH_AVERAGE_PLOT('DepthAverage')
    DepthAverage.ReadHeader(depth_average_path)
    DepthAverage.ReadData(depth_average_path)
    # manage data
    DepthAverage.SplitTimeStep()
    if fix_time_step and time_step > len(DepthAverage.time_step_times) - 1:
        time_step = len(DepthAverage.time_step_times) - 2
    try:
        i0 = DepthAverage.time_step_indexes[time_step][-1] * DepthAverage.time_step_length
        if time_step == len(DepthAverage.time_step_times) - 1:
            # this is the last step
            i1 = DepthAverage.data.shape[0]
        else:
            i1 = DepthAverage.time_step_indexes[time_step + 1][0] * DepthAverage.time_step_length
    except IndexError:
        print("PlotDaFigure: File (%s) may not contain any depth average output, abort" % depth_average_path)
        return
    names = kwargs.get('names', ['depth', 'temperature', 'adiabatic_density'])
    output_path = os.path.join(output_dir, 'depth_average_output_s%d' % time_step)
    odata = DepthAverage.export(output_path, names, rows=[i for i in range(i0, i1)], include_size=True)
    return odata, output_path


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
    parser.add_argument('-t', '--time', type=float,
                        default=0.0,
                        help='model time(s or yr)')
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
    if (_commend in ['-h', '--help']):
        # example:
        Usage()
    elif _commend == 'plot_by_time':
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
        PlotDaFigure(depth_average_path, fig_path_base, time=arg.time)
    
    if _commend == 'export_case':
        depth_average_path = os.path.join(arg.inputs, 'output', 'depth_average.txt')
        ExportData(depth_average_path, arg.outputs, time_step=arg.step)

# run script
if __name__ == '__main__':
    main()