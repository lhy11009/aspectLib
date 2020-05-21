import numpy as np
import sys
import os
import json
import re
from matplotlib import pyplot as plt
from importlib import resources
import shilofue.json

def PlotOptions(prefix):
    '''
    PlotOptions(prefix) 
    Read format options for plotting from json files
    Args:
        prefix(string):
            a prefix for options, all json files with this prefix will be imported
    Returns:
        options(dict):
            a dictionary for options. Every match with the prefix will generate a 
            key in this dictionary, while the value is read from a json file.
            For example, if the prefix is 'Statistics' and one json file is
            'Statistics_Number_of_Cells.json', then there will be a key called
            'Number_of_Cells' in this dictionary.
    '''
    _options = {}  # options is a dictionary
    for _filename in resources.contents(shilofue.json):
        # this resources.contents return a interable from a sub-package,
        # entries in this interable are filenames
        if re.match("^" + prefix + '_', _filename):
            # the following two lines eliminate the words before the first '_'
            # with that '_' as well as the '.json' in the end.
            # so that if filename is 'Statistics_Number_of_Cells.json',
            # _name is 'Number_of_Cells'
            _name = re.sub("^.*?_", '', _filename) # when 
            _name = re.sub(r".json", '', _name)
            with resources.open_text(shilofue.json, _filename) as fin:
                _options[_name] = json.load(fin)  # values are entries in this file
    assert(_options is not {})  # assert not vacant
    return _options


class UNITCONVERT():
    '''
    UNITCONVERT():
    class for unit convert

    '''

    class OptionNotValid(Exception):
        '''
        OptionNotValid(Exception)
        Error type

        '''
        pass

    def __init__(self, _filename):
        with resources.open_text(shilofue.json, _filename) as fin:
            _data = json.load(fin)
        self.Converts = _data['Converts']

    def __call__(self, _units):
        for option in self.Converts:
            if _units == option["units"]:
                return option["factor"]
        raise self.OptionNotValid("Error when getting the correct unit\
                             change option")


class STATISTICS():
    '''
    STATISTICS():
    class for statistic data

    '''

    def __init__(self):
        self.options = PlotOptions('Statistics')
    
    def __call__(self, _filename, **kwargs):
        '''
        Read and plot
        '''
        # canvas = kwargs.get('canvas', np.array([1, 1]))
        # ptype = kwargs['ptype']
        # Read plot options from a json file
        _jsonfile = kwargs.get('json', None)
        if _jsonfile is not None:
            with open(_jsonfile, 'r') as fin:
                _configs = json.load(fin)
        else:
            # default is to read from 'Statistics.json' in shilofue.json
            with resources.open_text(shilofue.json, 'Statistics.json') as fin:
                _configs = json.load(fin)
        self.ReadSTSHead(_filename)
        self.ReadSTS(_filename)
        # fileout = self.PlotSTSCombine(names=ptype, canvas=canvas)
        _fileout = self.PlotSTSCombine(_configs)
        return _fileout

    def ReadSTSHead(self, filename):
        '''
        Read header information
        '''
        self.header = {}
        assert(os.access(filename, os.R_OK))
        fin = open(filename, 'r')
        for line in fin:
            # derive column, unit and key from header
            if re.match("^#", line) is None:
                # only look at the first few lines starting with "#"
                break
            line = re.sub("^# ", '', line)
            match_obj = re.match("[0-9]*", line)
            col = int(match_obj[0]) - 1
            line = re.sub("^.*?: ", '', line)
            match_obj = re.search("[(].*[)]", line)
            if match_obj:
                unit = match_obj.group(0)
                unit = unit[1: -1]
            else:
                unit = None
            line = re.sub(" [(].*[)]", '', line)
            line = line.strip("\n")
            key = re.sub(" ", "_", line)
            # save information as key, options
            # options include col and unit ####
            info = {}
            info['col'] = col
            info['unit'] = unit
            self.header[key] = info

    def ReadSTS(self, filename):
        '''
        ReadSTS(self, filename):
        Read statistic from file
        '''
        # This file starts with lines of comment
        assert(os.access(filename, os.R_OK))
        self.data = np.genfromtxt(filename, comments='#')
    
    def PlotSTSCombine(self, _configs):
        '''
        Combine all plottings
        Arguments:
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
        _canvas = _configs.get('canvas', [1, 1])
        assert(type(_canvas) is list and len(_canvas) == 2)
        _types = _configs['types']  # types of plotting
        assert(type(_types) is list and
            _canvas[0] * _canvas[1] == len(_types))  # size of canvas match size of _types
        # plot
        fig, axs = plt.subplots(_canvas[0], _canvas[1])
        for i in range(len(_types)):
            _type = _types[i]
            _opt = self.options[_type]  # get the plot options
            try:
                # fix the bug where there is only one plot
                self.PlotSTS(axs[i], _opt)
            except TypeError:
                self.PlotSTS(axs, _opt)
        fig.tight_layout()
        _filename = './Statistics.pdf'
        fig.savefig(_filename)
        plt.close(fig)
        return _filename

    def PlotSTS(self, _ax, _opt):
        '''
        PlotSTS(self, fig, ax, _opt):
        Plot a statistic from data
        Arguments:
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
        # then column is determined
        # by header information
        _colx = self.header[_xname]['col']
        _coly = self.header[_yname]['col']
        _unitx = self.header[_xname]['unit']
        _unity = self.header[_yname]['unit']
        _x = self.data[:, _colx]
        _y = self.data[:, _coly]
        _ax.plot(_x, _y, _line, color=_color, label=_label)
        # construct label from xname and yname
        # put in unit
        if _unitx is not None:
            _xlabel = re.sub("_", " ", _xname) + ' (' + _unitx + ')'
        else:
            _xlabel = re.sub("_", " ", _xname)
        if _unity is not None:
            _ylabel = re.sub("_", " ", _yname) + ' (' + _unity + ')'
        else:
            _ylabel = re.sub("_", " ", _yname)
        _ax.set(xlabel=_xlabel, ylabel=_ylabel)
        _ax.legend()

