import sys
import os
import json
import re
import shilofue.json
import numpy as np
from importlib import resources
from shilofue.Utilities import JsonOptions
from shilofue.Utilities import ReadHeader
from shilofue.Utilities import ReadHeader2
from matplotlib import pyplot as plt

class LINEARPLOT():
    '''
    LINEARPLOT():
    class for LINEARPLOT

    '''

    def __init__(self, _name, **kwargs):
        '''
        _name(str):
            name of the plotting
        kwargs:
            unit_convert(fun):
                a unit_convert function, default is None
        '''
        self.options = JsonOptions(_name)
        self.UnitConvert = kwargs.get('unit_convert', None)
        self.dim = kwargs.get('dim', 2)  # dimension
        assert(self.dim in [1, 2, 3])  # dimension must be 1, 2, 3
    
    def __call__(self, _filename, **kwargs):
        '''
        Read and plot
        Attributes:
            _filename(string):
                filename for data file
        Returns:
            _fileout(string):
                filename for output figure
        '''
        # canvas = kwargs.get('canvas', np.array([1, 1]))
        # ptype = kwargs['ptype']
        # Read plot options from a json file
        _jsonfile = kwargs.get('json', None)
        _fileout = kwargs.get('fileout', _filename + '.pdf')
        if _jsonfile is not None:
            with open(_jsonfile, 'r') as fin:
                _configs = json.load(fin)
        else:
            # default is to read from 'Statistics.json' in shilofue.json
            with resources.open_text(shilofue.json, 'Statistics.json') as fin:
                _configs = json.load(fin)
        assert(os.access(_filename, os.R_OK))
        with open(_filename, 'r') as fin:
            _texts = fin.readlines()  # read the text of the file header
        self.header = ReadHeader(_texts)  # inteprate header information
        assert(os.access(_filename, os.R_OK))  # read in data
        self.data = np.genfromtxt(_filename, comments='#')
        _data_list = []
        for i in range(self.data.shape[1]):
            _data_list.append(self.data[:, i])
        _fileout = self.PlotCombine(_data_list, _fileout, _configs)
        return _fileout

    def PlotCombine(self, _data_list, _fileout, _configs):
        '''
        Combine all plottings
        Arguments:
            _data_list(list<ndarray>):
                list of data, each member is a set of data for
                some variable
            _fileout(str):
                name of the output file
            _configs(dict):
                canvas(ndarray):
                    layout of canvas
                types(list):
                    types of ploting, should have the 
                    same size with canvas
        Returns:
            _filename(string):
                name of the plotting created
        Raises:
            AssertionError:
                if size of canvas doesn't match number of plottings
        '''
        # plot configuration
        assert(type(_data_list) is list)
        _canvas = _configs.get('canvas', [1, 1])
        assert(type(_canvas) is list and len(_canvas) == 2)
        _types = _configs['types']  # types of plotting
        assert(type(_types) is list and
            _canvas[0] * _canvas[1] == len(_types))  # size of canvas match size of _types
        _size = _configs.get('size', (12, 12))  # size of the plot
        # plot
        fig, axs = plt.subplots(_canvas[0], _canvas[1], figsize=_size)
        for i in range(len(_types)):
            # find the right axis
            if _canvas[0] == 1 and _canvas[1] == 1:
                # if canvas is single, the axs is an axis
                _ax = axs
            elif _canvas[0] == 1 or _canvas[1] == 1:
                # if canvas is a number, then axs is 1-d array
                _ax = axs[i]
            else:
                # if canvas is two dimensional, axs is also two dimensional
                _id1 = i // _canvas[1]
                _id2 = i % _canvas[1]
                _ax = axs[_id1, _id2]
            # find options for plot
            _type = _types[i]
            if type(_type) == str:
                _opt = self.options[_type]  # get the plot options
                self.Plot(_data_list, _ax, _opt)
            elif type(_type) == list:
                # when _type is a list, plot multiple lines in a subplot
                for _subtype in _type:
                    assert(type(_subtype) == str)
                    _opt = self.options[_subtype]  # get the plot options
                    self.Plot(_data_list, _ax, _opt)
        fig.tight_layout()
        fig.savefig(_fileout)
        plt.close(fig)
        return _fileout

    def Plot(self, _data_list, _ax, _opt):
        '''
        Plot(self, fig, ax, _opt):
        Plot a line plotting from data
        Arguments:
            _data_list(list<ndarray>):
                list of data, each member is a set of data for
                some variable
            _ax(axis):
                an axis object from matplotlib
            _opt(dict):
                dictionary for plot options
        Returns:
            _ax(axis):
                the same as input
        '''
        # Options for plot
        # xname and yname are labels for x and y axis,
        # color is the color for ploting
        # label is the label for legend
        # line is the type of line for plotting
        _xname = _opt.get('xname', 'Time')
        _yname = _opt.get('yname', 'Number_of_mesh_cells')
        _color = _opt.get('color', 'r')
        _label = _opt.get('label', 'default')
        _line = _opt.get('line', '-')
        _invert_x = _opt.get('invert_x', 0)  # invert x axis
        _invert_y = _opt.get('invert_y', 0)
        _log_x = _opt.get('log_x', 0)  # plot x as log
        _log_y = _opt.get('log_y', 0)  # plot y as log
        # then column is determined
        # by header information
        _colx = self.header[_xname]['col']
        _coly = self.header[_yname]['col']
        _unitx = self.header[_xname]['unit']
        _unity = self.header[_yname]['unit']
        _x = _data_list[_colx]
        _y = _data_list[_coly]
        # get the unit for plot
        _unit_x_plot = _opt.get('xunit', _unitx)
        _unit_y_plot = _opt.get('yunit', _unity)
        # convert the unit
        if self.UnitConvert is not None and _unit_x_plot != _unitx:
            x_convert_ratio = self.UnitConvert(_unitx, _unit_x_plot)
        else:
            x_convert_ratio = 1.0
        if self.UnitConvert is not None and _unit_y_plot != _unity:
            y_convert_ratio = self.UnitConvert(_unity, _unit_y_plot)
        else:
            y_convert_ratio = 1.0
        # construct label from xname and yname
        # put in unit
        if _unit_x_plot is not None:
            _xlabel = re.sub("_", " ", _xname) + ' (' + _unit_x_plot + ')'
        else:
            _xlabel = re.sub("_", " ", _xname)
        if _unit_y_plot is not None:
            _ylabel = re.sub("_", " ", _yname) + ' (' + _unit_y_plot + ')'
        else:
            _ylabel = re.sub("_", " ", _yname)
        if _log_x and _log_y:
            _ax.loglog(_x*x_convert_ratio, _y*y_convert_ratio, _line, color=_color, label=_label)
        elif _log_x:
            _ax.semilogx(_x*x_convert_ratio, _y*y_convert_ratio, _line, color=_color, label=_label)
        elif _log_y:
            _ax.semilogy(_x*x_convert_ratio, _y*y_convert_ratio, _line, color=_color, label=_label)
        else:
            _ax.plot(_x*x_convert_ratio, _y*y_convert_ratio, _line, color=_color, label=_label)
        _ax.set(xlabel=_xlabel, ylabel=_ylabel)
        if _invert_x and ~_ax.xaxis_inverted():
            _ax.invert_xaxis()
        if _invert_y and ~_ax.yaxis_inverted():
            _ax.invert_yaxis()
        _ax.legend()


class DEPTH_AVERAGE_PLOT(LINEARPLOT):
    '''
    Class for plotting depth average file.
    This is inheritage of the LINEARPLOT class

    Attributes:
        todo
    Args:
        __init__():
            initiation
        todo
    '''

    def __call__(self, _filename, **kwargs):
        '''
        call function of this class
        Read and plot

        Args:
            _filename(string):
                filename for data file
        Returns:
            _fileout(string):
                filename for output figure
        '''
        _fileout = kwargs.get('fileout', _filename + '.pdf')
        _jsonfile = kwargs.get('json', None)
        if _jsonfile is not None:
            with open(_jsonfile, 'r') as fin:
                _configs = json.load(fin)
        else:
            # default is to read from 'DepthAverage.json' in shilofue.json
            with resources.open_text(shilofue.json, 'DepthAverage.json') as fin:
                _configs = json.load(fin)
        assert(os.access(_filename, os.R_OK))
        with open(_filename, 'r') as fin:
            _texts = fin.readlines()  # read the text of the file header
        self.header = ReadHeader2(_texts)  # inteprate header information
        assert(os.access(_filename, os.R_OK))  # read in data
        self.data = np.genfromtxt(_filename, comments='#')
        _data_list = self.manage_data()
        self.manage_units()
        _fileout = self.PlotCombine(_data_list, _fileout, _configs)
        return _fileout
    
    def manage_data(self):
        '''
        manage data, get new data for this class
        Returns:
            _data_list(list):
                list of data for ploting
        '''
        _data_list = []
        for i in range(self.data.shape[1]):
            _data_list.append(self.data[:, i])
        # get the super adiabatic temperature
        _col_temperature = self.header['temperature']['col']
        _col_adiabatic_temperature = self.header['adiabatic_temperature']['col']
        _super_adiabatic_temperature = self.data[:, _col_temperature] - self.data[:, _col_adiabatic_temperature]
        _data_list.append(_super_adiabatic_temperature)
        self.header['super_adiabatic_temperature'] = {}
        self.header['super_adiabatic_temperature']['col'] = self.header['total_col']
        self.header['super_adiabatic_temperature']['unit'] = None
        self.header['total_col'] += 1
        return _data_list
    
    def manage_units(self):
        '''
        manage units, get units for data.
        This is due to the bad form of the header of this file
        '''
        self.header['depth']['unit'] = 'm'
        self.header['temperature']['unit'] = 'K'
        self.header['super_adiabatic_temperature']['unit'] = 'K'
        self.header['viscosity']['unit'] = 'Pa s'
        self.header['velocity_magnitude']['unit'] = 'm/yr'
        if self.dim == 2:
            self.header['vertical_heat_flux']['unit'] = 'mw/m'
        elif self.dim == 3:
            self.header['vertical_heat_flux']['unit'] = 'mw/m^2'
