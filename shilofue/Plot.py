import sys
import os
import json
import re
import shilofue.json
import numpy as np
from importlib import resources
from shilofue.Utilities import JsonOptions
from shilofue.Utilities import ReadHeader
from matplotlib import pyplot as plt

class STATISTICS():
    '''
    STATISTICS():
    class for statistic data

    '''

    def __init__(self):
        self.options = JsonOptions('Statistics')
    
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
        assert(os.access(_filename, os.R_OK))
        with open(_filename, 'r') as fin:
            _texts = fin.readlines()  # read the text of the file header
        self.header = ReadHeader(_texts)  # inteprate header information
        self.ReadSTS(_filename)
        # fileout = self.PlotSTSCombine(names=ptype, canvas=canvas)
        _fileout = self.PlotSTSCombine(_configs)
        return _fileout


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

