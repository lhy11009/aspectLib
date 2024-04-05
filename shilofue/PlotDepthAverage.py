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
from matplotlib import gridspec
import shilofue.Plot as Plot
import shilofue.PostHefesto as PHefesto


# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
JSON_FILE_DIR = os.path.join(ASPECT_LAB_DIR, "files", "json_examples")
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
        self.header = Utilities.ReadHeader2(_texts)


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


def PlotDaFigureByName(depth_average_path, **kwargs):
    '''
    plot figure
    Inputs:
        kwargs:
            time_step - time_step to plot the figure, default is 0
    '''
    time = kwargs.get('time', 0.0)
    assert(os.access(depth_average_path, os.R_OK))
    axT = kwargs.get('axT', None)
    axP = kwargs.get('axP', None)
    ax_density = kwargs.get('ax_density', None)
    # read that
    DepthAverage = DEPTH_AVERAGE_PLOT('DepthAverage')
    DepthAverage.ReadHeader(depth_average_path)
    DepthAverage.ReadData(depth_average_path)
    
    # manage data
    DepthAverage.SplitTimeStep()
    names = ['depth', 'adiabatic_pressure', 'temperature', 'adiabatic_temperature', 'viscosity', 'vertical_heat_flux', 'vertical_mass_flux', 'adiabatic_density']
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

    # plot temperature
    if axT is not None:
        axT.plot(temperatures, depths/1e3, "-", label="T", color="tab:red")
        axT.set_ylabel("Depth [km]")
        axT.set_xlabel("Temperature [K]")
        axT.legend()

    # plot pressure 
    if axP is not None:
        axP.plot(pressures / 1e9, depths/1e3, "-", label="Pressure", color="tab:green")
        axP.set_ylabel("Depth [km]")
        axP.set_xlabel("Pressure [GPa]")
        axT.legend()

    # plot density 
    if ax_density is not None:
        ax_density.plot(densities, depths/1e3, "-", label="Density", color="tab:blue")
        ax_density.set_ylabel("Depth [km]")
        ax_density.set_xlabel("Density [kg/m^3]")
        ax_density.legend()


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


def ComputeBuoyancy(da_file0, da_file1, **kwargs):
    '''
    Inputs:
        da_file0 (str): reference profile
        da_file1 (str): another profile
        kwargs:
            odir (str): output directory
            axT (matplotlib axis): axis to plot the temperature
            ax_density_ratio (matplotlib axis): axis to plot the density ratio
            ax_buoy (matplotlib axis): axis to plot the buoyancy
    '''
    assert(os.path.isfile(da_file0))
    assert(os.path.isfile(da_file1))
    axT = kwargs.get('axT', None)
    ax_density_ratio = kwargs.get('ax_density_ratio', None)
    ax_buoy = kwargs.get('ax_buoy', None)
    odir = kwargs.get("odir", RESULT_DIR)
    max_depth = 2890e3
    n_depth = 2891
    g = 9.8
    interval = max_depth / (n_depth - 1)

    # read data
    odata0, _ = ExportData(da_file0, RESULT_DIR, names=['depth', 'temperature', 'adiabatic_density'])
    odata1, _ = ExportData(da_file1, RESULT_DIR, names=['depth', 'temperature', 'adiabatic_density'])
    depths_0 = odata0[:, 0]
    Ts_0 = odata0[:, 1]
    densities_0 = odata0[:, 2]
    depths_1 = odata1[:, 0]
    Ts_1 = odata1[:, 1]
    densities_1 = odata1[:, 2]

    # interpolate data
    DensityFunc0 = interp1d(depths_0, densities_0, assume_sorted=True, fill_value="extrapolate")
    TFunc0 = interp1d(depths_0, Ts_0, assume_sorted=True, fill_value="extrapolate")
    DensityFunc1 = interp1d(depths_1, densities_1, assume_sorted=True, fill_value="extrapolate")
    TFunc1 = interp1d(depths_1, Ts_1, assume_sorted=True, fill_value="extrapolate")

    # compute buoyancy and buoyancy number
    # the buoyancy number is computed with buoyancy / density1
    # density1 is chosen instead of density0 to simulate the buoyancy ratio of density0
    depths = np.linspace(0, max_depth, n_depth)
    buoyancies = np.zeros(n_depth)
    density_ratios = np.zeros(n_depth)
    adiabatic_diff_densities = np.zeros(n_depth)
    alpha = 3.1e-5  # thermal expansivity: contant value
    for i in range(n_depth):
        depth = depths[i]
        density0 = DensityFunc0(depth)
        density1 = DensityFunc1(depth)
        T0 = TFunc0(depth)
        T1 = TFunc1(depth)
        adiabatic_diff_density = density1 - density0
        diff_density = adiabatic_diff_density + alpha*(T0 - T1)*density1
        buoyancy = -diff_density * g
        buoyancies[i] = buoyancy
        density_ratios[i] = diff_density / density1
    
    # plot temperature
    if axT is not None:
        axT.plot(Ts_0, depths_0/1e3, "-", label="T0", color="tab:red")
        axT.plot(Ts_1, depths_1/1e3, "-", label="T1", color="tab:blue")
        axT.set_ylabel("Depth [km]")
        axT.set_xlabel("Temperature [K]")
        axT.legend()

    # plot density ratio
    if ax_density_ratio is not None:
        ax_density_ratio.plot(density_ratios, depths/1e3, label="density ratio", color="tab:red")
        ax_density_ratio.set_ylabel("Depth [km]")
        ax_density_ratio.set_xlabel("Density Ratio")
        ax_density_ratio.legend()
    
    # plot buoyancy
    if ax_buoy is not None:
        ax_buoy.plot(buoyancies, depths/1e3, label="buoyancy", color="tab:blue")
        ax_buoy.set_ylabel("Depth [km]")
        ax_buoy.set_xlabel("Buoyancy [N/m^3]")
        ax_buoy.legend()

    # return variables
    return depths, buoyancies, density_ratios


def PlotBuoyancy(da_file0, da_file1, **kwargs):
    '''
    Plot Buoyancy from two profiles
    Inputs:
        da_file0 (str): reference profile
        da_file1 (str): another profile
        kwargs:
            odir: output directory
    '''
    odir = kwargs.get("odir", RESULT_DIR)
    fig_path = os.path.join(odir, "buoyancy.png")
    fig = plt.figure(tight_layout=True, figsize=(6, 15))
    gs = gridspec.GridSpec(3, 1)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[2, 0])
    if os.path.isfile(fig_path):
        os.remove(fig_path)

    # plot figures
    depths, buoyancies, buoyancy_ratios = ComputeBuoyancy(da_file0, da_file1, odir=odir, ax1=ax1, ax2=ax2, ax3=ax3)
    ax1.grid()
    ax2.grid()
    ax1.invert_yaxis()
    ax2.invert_yaxis()
    # save figures
    fig_path = os.path.join(fig_path)
    if os.path.isfile(fig_path):
        # remove previous files
        os.remove(fig_path)
    fig.savefig(fig_path)
    assert(os.path.isfile(fig_path))
    print("%s: New figue %s" % (Utilities.func_name(), fig_path))


def CompareHefestoBuoyancy(**kwargs):
    '''
    Compare the buoyancy from two ASPECT depth_average file
    to the buoyancy from two Hefesto files
    '''
    TwoDSubduction_DIR = os.environ['TwoDSubduction_DIR']
    file_type = kwargs.get("file_type", "png")
    source_dir = os.path.join(ASPECT_LAB_DIR, "tests", "integration", 'fixtures', 'post_hefesto')
    fort56_file0 = os.path.join(source_dir, "fort.56.0")
    fort56_file1 = os.path.join(source_dir, "fort.56.1")
    # fort56_file1 = "/home/lochy/Softwares/HeFESTo/HeFESToRepository/output_2.7749999983454376/fort.56"
    # fort56_file1 = "/home/lochy/Softwares/HeFESTo/HeFESToRepository/output_2.0249999999999124/fort.56"
    assert(os.path.isfile(fort56_file0) and os.path.isfile(fort56_file1))
    da_file0 = os.path.join(TwoDSubduction_DIR, "HeFesto_Compare/adb_T_1650.1_fix_width_refine/output/depth_average.txt")
    da_file1 = os.path.join(TwoDSubduction_DIR, "HeFesto_Compare/adb_T_676.23_fix_width_refine/output/depth_average.txt")
    # da_file1 = "/home/lochy/ASPECT_PROJECT/TwoDSubduction_DIR/HeFesto_Compare/adb_T_1921.7_fix_width_refine/output/depth_average.txt"
    # da_file1 = "/home/lochy/ASPECT_PROJECT/TwoDSubduction_DIR/HeFesto_Compare/adb_T_1059.4_fixT_2/output/depth_average.txt"
    assert(os.path.isfile(da_file0) and os.path.isfile(da_file1))

    # first compare outputs from ASPECT and HeFesto
    fig_path = os.path.join(RESULT_DIR, "HeFesto_DA_compare.png")
    fig = plt.figure(tight_layout=True, figsize=(6, 18))
    gs = gridspec.GridSpec(3, 1)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[2, 0])
    # plot depth average profile
    PlotDaFigureByName(da_file0, axT=ax1, ax_density=ax3, axP=ax2)
    # plot HeFesto Profile
    PHefesto.PlotHeFestoProfile(fort56_file0, axT=ax1, ax_density=ax3, axP=ax2)
    ax1.invert_yaxis()
    ax1.grid()
    ax2.invert_yaxis()
    ax2.grid()
    ax3.invert_yaxis()
    ax3.grid()
    fig.savefig(fig_path)
    assert(os.path.isfile(fig_path))
    print("%s: new figure %s" % (Utilities.func_name(), fig_path))
    ax1.set_ylim([1000, 0])
    ax2.set_ylim([1000, 0])
    ax3.set_ylim([1000, 0])
    fig_path = os.path.join(RESULT_DIR, "HeFesto_DA_compare_1000.%s" % file_type)
    fig.savefig(fig_path)
    assert(os.path.isfile(fig_path))
    print("%s: new figure %s" % (Utilities.func_name(), fig_path))


    # first compare outputs of cold internal from ASPECT and HeFesto
    fig_path = os.path.join(RESULT_DIR, "HeFesto_DA_compare_cold.%s" % file_type)
    fig = plt.figure(tight_layout=True, figsize=(6, 18))
    gs = gridspec.GridSpec(3, 1)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[2, 0])
    # plot depth average profile
    PlotDaFigureByName(da_file1, axT=ax1, ax_density=ax3, axP=ax2)
    # plot HeFesto Profile
    PHefesto.PlotHeFestoProfile(fort56_file1, axT=ax1, ax_density=ax3, axP=ax2)
    ax1.invert_yaxis()
    ax1.grid()
    ax2.invert_yaxis()
    ax2.grid()
    ax3.invert_yaxis()
    ax3.grid()
    fig.savefig(fig_path)
    assert(os.path.isfile(fig_path))
    print("%s: new figure %s" % (Utilities.func_name(), fig_path))
    ax1.set_ylim([1000, 0])
    ax2.set_ylim([1000, 0])
    ax3.set_ylim([1000, 0])
    fig_path = os.path.join(RESULT_DIR, "HeFesto_DA_compare_cold_1000.%s" % file_type)
    fig.savefig(fig_path)
    assert(os.path.isfile(fig_path))
    print("%s: new figure %s" % (Utilities.func_name(), fig_path))
    
    
    # then plot the buoyancy forces 
    fig_path = os.path.join(RESULT_DIR, "buoyancy_hefesto_compare.%s" % file_type)
    # initial plot 
    fig = plt.figure(tight_layout=True, figsize=(6, 15))
    gs = gridspec.GridSpec(3, 1)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[2, 0])
    if os.path.isfile(fig_path):
        os.remove(fig_path)
    # plot figures with Hefesto
    _, buoyancies, buoyancy_ratios = PHefesto.ComputeBuoyancy(fort56_file0, fort56_file1, axT=ax1,\
                                                              ax_density_ratio=ax2, ax_buoy=ax3)
    # plot figures with ASPECT profiles
    _, buoyancies, buoyancy_ratios = ComputeBuoyancy(da_file0, da_file1, axT=ax1, ax_density_ratio=ax2, ax_buoy=ax3)
    # save figures
    ax1.set_ylim([0.0, 1000.0])
    ax1.invert_yaxis()
    ax1.grid()
    ax2.set_ylim([0.0, 1000.0])
    ax2.invert_yaxis()
    ax2.grid()
    ax3.set_ylim([0.0, 1000.0])
    ax3.invert_yaxis()
    ax3.grid()
    # save
    fig.savefig(fig_path)
    assert(os.path.isfile(fig_path))
    print("%s: new figure %s" % (Utilities.func_name(), fig_path))



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
    parser.add_argument('-i1', '--inputs1', type=str,
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
    parser.add_argument('-ft', '--file_type', type=str,
                        default="png",
                        help='file output type')
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
        json_file = os.path.join(JSON_FILE_DIR, 'post_process.json')
        assert(os.access(json_file, os.R_OK))
        with open(json_file, 'r') as fin:
            json_options = json.load(fin)
        plot_options = json_options.get('DepthAverage', {})
    
        # Init the UnitConvert class
        UnitConvert = Utilities.UNITCONVERT()

        # Initiate class for depth average plot
        DepthAverage = DEPTH_AVERAGE_PLOT('DepthAverage', unit_convert=UnitConvert, options=plot_options)

        # plot figure
        assert(os.access(arg.inputs, os.R_OK))
        if not os.path.isdir(os.path.dirname(arg.outputs)):
            os.mkdir(os.path.dirname(arg.outputs))
        DepthAverage(arg.inputs, fileout=arg.outputs, time=arg.time)

    # commands
    elif _commend == 'plot_case':
        depth_average_path = os.path.join(arg.inputs, 'output', 'depth_average.txt')
        fig_path_base = os.path.join(arg.inputs, 'img', 'DepthAverage.png')
        PlotDaFigure(depth_average_path, fig_path_base, time=arg.time)
    
    elif _commend == 'export_case':
        depth_average_path = os.path.join(arg.inputs, 'output', 'depth_average.txt')
        ExportData(depth_average_path, arg.outputs, time_step=arg.step)
    
    elif _commend == 'plot_buoyancy':
        PlotBuoyancy(arg.inputs, arg.inputs1)
    
    elif _commend == 'compare_hefesto_buoyancy':
        CompareHefestoBuoyancy(file_type=arg.file_type)

# run script
if __name__ == '__main__':
    main()