import numpy as np
import sys
import os
import json
import re
from matplotlib import pyplot as plt


jsonDir = 'json'  # json files folder


################################################################################
# Check OS status of file, end program if failed
def OS_Check(filename, status):
    if status is 'r':
        if not os.access(filename, os.R_OK):
            print("Read access not permitted on %s" % filename)
            sys.exit(1)
    elif status is 'w':
        if not os.access(filename, os.W_OK):
            print("Write access not permitted on %s" % filename)
            sys.exit(1)
    elif status is 'a':
        if not os.access(filename, os.W_OK):
            print("Write at the end of the file access not permitted on %s" % filename)
            sys.exit(1)
    elif status is 'f':
        if not os.access(filename, os.F_OK):
            print("File not Found on %s" % filename)
            sys.exit(1)
    else:
        print("Wrong input for parameter 'status'(second) in Function OS_Check")
        sys.exit(1)


################################################################################
# Read variable value and assign default if not found ##########
def Config(_kwargs, _name, _default):
    try:
        value = _kwargs[_name]
    except KeyError:
        value = _default
    return value


################################################################################
# Dump configure in json
def JsonDump(_dict, _name):
    if not os.path.isdir(jsonDir):
        os.mkdir(jsonDir)
    filename = os.path.join(jsonDir, "%s.json" % _name)
    with open(filename, 'w') as fout:
        json.dump(_dict, fout)


################################################################################
# Read configuration from json
def JsonRead(_filename):
    filename = os.path.join(jsonDir, _filename)
    OS_Check(filename, 'r')
    fin = open(filename, 'r')
    data = json.load(fin)
    return data


################################################################################
# Format options
def PlotOptions(prefix):
    Options = {}
    for dirname, dirnames, filenames in os.walk(jsonDir):
        for filename in filenames:
            if re.match("^" + prefix, filename):
                _name = re.sub("^.*?_", '', filename)
                _name = re.sub(r".json", '', _name)
                Options[_name] = JsonRead(filename)
    return Options


################################################################################
# class for unit convert
class UNITCONVERT():
    # Error type
    class OptionNotValid(Exception):
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


################################################################################
# class for statistic data
class STATISTICS():
    def __init__(self):
        self.Options = PlotOptions('Statistics')

    # Read header information
    def ReadSTSHead(self, filename):
        self.header = {}
        OS_Check(filename, 'r')
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

    # Read statistic from file ##########
    def ReadSTS(self, filename):
        # This file starts with lines of comment
        OS_Check(filename, 'r')
        self.data = np.genfromtxt(filename, comments='#')

    # Plot statistic from data ##########
    def PlotSTS(self, fig, ax, _opt):
        # xname and yname are from json file, then column is determined
        # by header information ####
        xname = Config(_opt, 'xname', 'Time')
        yname = Config(_opt, 'yname', 'Number_of_mesh_cells')
        print("xname:", xname)
        print("yname:", yname)  # debug
        colx = self.header[xname]['col']
        coly = self.header[yname]['col']
        unitx = self.header[xname]['unit']
        unity = self.header[yname]['unit']
        x = self.data[:, colx]
        y = self.data[:, coly]
        color = Config(_opt, 'color', 'r')
        Label = Config(_opt, 'label', 'default')
        Line = Config(_opt, 'line', '-')
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

    # Plot statistic from data ##########
    def PlotSTSCombine(self, **kwargs):
        canvas = Config(kwargs, 'canvas', np.array([1, 1]))
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

    # Read and plot
    def __call__(self, filename, **kwargs):
        canvas = Config(kwargs, 'canvas', np.array([1, 1]))
        ptype = kwargs['ptype']
        self.ReadSTSHead(filename)
        self.ReadSTS(filename)
        print(self.header)
        fileout = self.PlotSTSCombine(ptype=ptype, canvas=canvas)
        return fileout
