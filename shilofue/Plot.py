import numpy as np
import sys
import os
import json
import re
from matplotlib import pyplot as plt
from importlib_resources import files
import shilofue.data

jsonDir = 'json'  # json files folder


def JsonDump(_dict, _name):
    '''
    JsonDump(_dict, _name):
    Dump configure in json

    '''
    if not os.path.isdir(jsonDir):
        os.mkdir(jsonDir)
    filename = os.path.join(jsonDir, "%s.json" % _name)
    with open(filename, 'w') as fout:
        json.dump(_dict, fout)


def JsonRead(_filename):
    '''
    JsonRead(_filename):
    Read configuration from json

    '''
    filename = os.path.join(jsonDir, _filename)
    assert(os.access(filename, os.R_OK))
    # source = files(shilofue.data).joinpath(_filename)
    # with as_file(source) as f_json:
        # data = json.load(f_json)
    fin = open(filename, 'r')
    data = json.load(fin)
    return data


def PlotOptions(prefix):
    '''
    PlotOptions(prefix):
    Format options

    '''
    Options = {}
    for dirname, dirnames, filenames in os.walk(jsonDir):
        for filename in filenames:
            if re.match("^" + prefix, filename):
                _name = re.sub("^.*?_", '', filename)
                _name = re.sub(r".json", '', _name)
                Options[_name] = JsonRead(filename)
    return Options


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

    def __init__(self, filename):
        data = JsonRead(filename)
        self.Converts = data['Converts']

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
        self.Options = PlotOptions('Statistics')

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
            option = {}
            option['col'] = col
            option['unit'] = unit
            self.header[key] = option

    def ReadSTS(self, filename):
        '''
        ReadSTS(self, filename):
        Read statistic from file
        '''
        # This file starts with lines of comment
        assert(os.access(filename, os.R_OK))
        self.data = np.genfromtxt(filename, comments='#')

    def PlotSTS(self, fig, ax, _opt):
        '''
        PlotSTS(self, fig, ax, _opt):
        Plot statistic from data
        '''
        # xname and yname are from json file, then column is determined
        # by header information ####
        xname = _opt.get('xname', 'Time')
        yname = _opt.get('yname', 'Number_of_mesh_cells')
        color = _opt.get('color', 'r')
        Label = _opt.get('label', 'default')
        Line = _opt.get('line', '-')
        print("xname:", xname)
        print("yname:", yname)  # debug
        colx = self.header[xname]['col']
        coly = self.header[yname]['col']
        unitx = self.header[xname]['unit']
        unity = self.header[yname]['unit']
        x = self.data[:, colx]
        y = self.data[:, coly]
        ax.plot(x, y, '%s%s' % (Line, color), label=Label)
        # construct label from xname and yname ####
        if unitx is not None:
            xLabel = re.sub("_", " ", xname) + ' [' + unitx + ']'
        else:
            xLabel = re.sub("_", " ", xname)
        if unity is not None:
            yLabel = re.sub("_", " ", yname) + ' [' + unity + ']'
        else:
            yLabel = re.sub("_", " ", yname)
        ax.set(xlabel=xLabel, ylabel=yLabel)
        ax.legend()
        fig.tight_layout()

    def PlotSTSCombine(self, **kwargs):
        '''
        Plot statistic from data
        '''
        canvas = kwargs.get('canvas', np.array([1, 1]))
        ptype = kwargs['ptype']
        print("ptype: ", ptype)  # debug
        fig, axs = plt.subplots(canvas[0], canvas[1])
        for i in range(len(ptype)):
            _type = ptype[i]
            _opt = self.Options[_type]
            try:
                self.PlotSTS(fig, axs[i], _opt)
            except TypeError:
                self.PlotSTS(fig, axs, _opt)
        filename = './Statistics.pdf'
        fig.savefig(filename)
        plt.close(fig)
        return filename

    def __call__(self, filename, **kwargs):
        '''
        Read and plot
        '''
        canvas = kwargs.get('canvas', np.array([1, 1]))
        ptype = kwargs['ptype']
        self.ReadSTSHead(filename)
        self.ReadSTS(filename)
        print(self.header)
        fileout = self.PlotSTSCombine(ptype=ptype, canvas=canvas)
        return fileout
