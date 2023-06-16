import sys
import os
import json
import re
import shilofue.json_files
import warnings
import argparse
import subprocess
import pathlib
import numpy as np
from importlib import resources
from shilofue.ParsePrm import ParseFromDealiiInput
from matplotlib import pyplot as plt


# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities


class LINEARPLOT():
    '''
    LINEARPLOT():
    class for LINEARPLOT

    '''
    def __init__(self, _name, options={}):
        '''
        options contains the instructions to plot figures, if none is given,
        then default options with be read in through json file.
        _name(str):
            name of the plotting
        options:
            unit_convert(fun):
                a unit_convert function, default is None
        '''
        self.name = _name
        _json_dir = options.get('json_dir', None)
        self.options = Utilities.JsonOptions(_name, _json_dir)
        self.UnitConvert = options.get('unit_convert', None)
        self.dim = options.get('dim', 2)  # dimension
        assert(self.dim in [1, 2, 3])  # dimension must be 1, 2, 3

        # reset the options with a option in the options
        with resources.open_text(shilofue.json_files, 'post_process.json') as fin:
            all_options = json.load(fin)
        self.options = all_options.get(self.name, {})
        try:
            options = options['options']
        except KeyError:
            pass
        else:
            self.options.update(options)

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
        _fileout = kwargs.get('fileout', _filename + '.pdf')

        # inteprate header information
        self.ReadHeader(_filename)

        # Read data
        # catch possible vacant file
        state=self.ReadData(_filename)
        if state == 1:
            return None

        # manage output data
        _data_list = self.ManageData()
        _fileout = self.PlotCombine(_data_list, _fileout)
        return _fileout

    class DataNotFoundWarning(UserWarning):
        # handle the circumstance that one data field is not found,
        # but continue the plot
        pass

    def ReadHeader(self, _filename):
        '''
        Read header information from file.
        An example of string is:
        '# 1: Time (years)'
        Args:
            _filename(str):
                filename for data file
        '''
        if not os.access(_filename, os.R_OK):  # read in data
            raise FileExistsError("%s cannot be read." % _filename)
        with open(_filename, 'r') as fin:
            _texts = fin.readlines()  # read the text of the file header
        try:
            self.header = Utilities.ReadHeader(_texts)
        except:
            raise Exception('Header for file %s cannot be read' % _filename)

    def SortHeader(self):
        '''
        sort header
        '''
        names = []
        cols = []
        units = []
        for key, value in self.header.items():
            if key == 'total_col':
                continue
            names.append(key)
            cols.append(int(value['col']))
            units.append(value['unit'])
        names = np.array(names)
        cols = np.array(cols)
        units = np.array(units)
        sort_indexes = np.argsort(cols)
        names = names[sort_indexes]
        cols = cols[sort_indexes]
        units = units[sort_indexes]
        return cols, names, units  ## debug

    def ReadData(self, _filename):
        '''
        Read Data
        Attributes:
            _filename(string):
                filename for data file
        '''
        if not os.access(_filename, os.R_OK):  # read in data
            raise FileExistsError("%s cannot be read." % _filename)
        # self.data = np.genfromtxt(_filename, comments='#')

        # import data via numpy buid in method
        # catch warning of empty file and return 1
        with warnings.catch_warnings(record=True) as w:
            self.data = np.genfromtxt(_filename, comments='#', filling_values=0.0)
            if (len(w) > 0):
                assert(issubclass(w[-1].category, UserWarning))
                assert('Empty input file' in str(w[-1].message))
                warnings.warn('ReadData: %s, abort' % str(w[-1].message))
                return 1

        if len(self.data.shape) == 1:
            # only one row, expand it too 2-d array
            self.data = np.array([self.data])
        return 0
    
    def ReadPrm(self, prm_file):
        '''
        Read prm file
        '''
        # read prm file
        assert(os.path.isfile(prm_file))
        with open(prm_file, 'r') as fin:
            self.prm = ParseFromDealiiInput(fin)


    def ManageData(self):
        '''
        manage data, get new data for this class
        for the base class, this method simply takes the combination
        of self.data
        Returns:
            _data_list(list):
                list of data for ploting
        '''
        _data_list = []
        for i in range(self.data.shape[1]):
            _data_list.append(self.data[:, i])
        return _data_list


    def Has(self, field_name):
        '''
        check if data has entry of field
        Inputs:
            field_name: name of field
        Returns:
            true: if field is present
            false: if field is not present
        '''
        has_field = (field_name in self.header)
        return has_field


    def HasData(self):
        '''
        Return true if there is data
        Returns:
            True: there is data
            False: there is no data
        '''
        if (self.data.size == 0):
            return False
        else:
            return True

    def PlotCombine(self, _data_list, _fileout, **kwargs):
        '''
        Combine all plottings
        Arguments:
            _data_list(list<ndarray>):
                list of data, each member is a set of data for
                some variable
            _fileout(str):
                name of the output file
            **kwargs:
                title(str):
                    title of the ploting
        Returns:
            _filename(string):
                name of the plotting created
        Raises:
            AssertionError:
                if size of canvas doesn't match number of plottings
        '''
        # plot configuration
        assert(type(_data_list) is list)
        _canvas = self.options.get('canvas', [1, 1])
        assert(type(_canvas) is list and len(_canvas) == 2)
        _types = self.options.get('types', [])  # types of plotting
        assert(type(_types) is list and
            _canvas[0] * _canvas[1] >= len(_types))  # size of canvas match size of _types
        _size = self.options.get('size', (5, 5))  # size of the plot
        _title = kwargs.get('title', None)  # get title

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
        # xname and yname are labels for x and y ax
        _xname = _opt.get('xname', 'Time')
        _yname = _opt.get('yname', 'Number_of_mesh_cells')
        _color = _opt.get('color', 'r')
        _label = _opt.get('label', None)
        _line = _opt.get('line', '-')
        _invert_x = _opt.get('invert_x', 0)  # invert x axis
        _invert_y = _opt.get('invert_y', 0)
        _log_x = _opt.get('log_x', 0)  # plot x as log
        _log_y = _opt.get('log_y', 0)  # plot y as log
        # then column is determined
        # by header information
        _colx = self.header[_xname]['col']
        try:
            _coly = self.header[_yname]['col']
        except KeyError:
            # there is no such field in the file
            # give out a warning, continue with the next plot
            warnings.warn('The field %s doesn\'t exist. We will keep ploting,\
but you will get a blank one for this field name' % _yname,
                           self.DataNotFoundWarning)
            return
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

        # legend
        if _label is not None:
            _ax.legend()


    def export(self, output_path, names, **kwargs):
        '''
        export data to file
        kwargs:
            rows: assign rows to output, default is None
        '''
        Utilities.my_assert(type(names) == list, TypeError, "%s: names must be a list" % Utilities.func_name())
        cols = []
        for _name in names:
            cols.append(int(self.header[_name]['col']))
        rows = kwargs.get("rows", None)
        if rows != None:
            my_assert((type(rows) == list), TypeError, "%s: rows must be a list" % Utilities.func_name())
            odata = self.data[np.ix_(rows, cols)]
        else:
            odata = self.data[:, cols]
        include_size=kwargs.get('include_size', False)
        data_only = kwargs.get('data_only', False)  # only return data, but not file
        if not data_only:
            with open(output_path, 'w') as fout:
                i = 1
                # header
                for _name in names:
                    fout.write("# %d: %s\n" % (i, _name))
                    i += 1
                if include_size:
                    fout.write("%d %d\n" % (odata.shape[0], odata.shape[1]))
                # data
                np.savetxt(fout, odata, fmt='%-20.8e')
            print('Export data to %s' % output_path)
            print('\tData layout: ', odata.shape)
        return odata


class STATISTICS_PLOT_OLD(LINEARPLOT):
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
        future
        Inputs:
            step(int)
        get time to a value of model step
        '''
        time = 0.0
        return time



class NEWTON_SOLVER_PLOT(LINEARPLOT):
    '''
    Class for plotting depth average file.
    This is an inheritage of the LINEARPLOT class

    Attributes:
    Args:
    '''
    def __init__(self, _name, **kwargs):
        LINEARPLOT.__init__(self, _name, kwargs)  # call init from base function

        # initiate
        # if step is None, then plot for all steps, otherwise plot for a step
        self.step = None

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
        _fileout = kwargs.get('fileout', _filename + '.pdf')
        self.ReadHeader(_filename)  # inteprate header information

        # catch possible vacant file
        state = self.ReadData(_filename)  # read data
        if state == 1:
            return None

        # manage output data
        _data_list = self.ManageData()

        # figure out name of file
        if self.step is None:
            # plot for all
            _fname = _fileout
        else:
            _fname_base = _fileout.rpartition('.')[0]
            _fname_type = _fileout.rpartition('.')[2]
            _fname = "%s_s%07d.%s" % (_fname_base, self.step, _fname_type)

        # plot
        _fileout = self.PlotCombine(_data_list, _fname)
        return _fileout

    def GetStep(self, step):
        '''
        Get time step
        Inputs:
            step(int): index of time step
        '''
        self.step = step

    def ManageData(self):
        '''
        manage data, get new data for this class
        Returns:
            _data_list(list):
                list of data for ploting
        '''
        if self.step is None:
            _data_list = self.ManageDataAll()
        else:
            _data_list = self.ManageDataStep()
        return _data_list

    def ManageDataStep(self):
        '''
        manage data for a single step
        Returns:
            _data_list(list):
                list of data for ploting
        '''
        _data_list = []

        # get list of step
        col_step = self.header['Time_step_number']['col']
        mask_step = (self.data[:, col_step] == self.step)
        for i in range(self.data.shape[1]):
            _data_list.append(self.data[mask_step, i])
        return _data_list

    def ManageDataAll(self):
        '''
        manage data for a single step
        Returns:
            _data_list(list):
                list of data for ploting
        '''
        _data_list = []
        for i in range(self.data.shape[1]):
            _data_list.append(self.data[:, i])

        # get number of nonlinear iteration
        nni = np.array([i for i in range(self.data.shape[0])])
        _data_list.append(nni)
        # mend header
        self.header['Number_of_nonlinear_iteration'] = {}
        self.header['Number_of_nonlinear_iteration']['col'] = self.header['total_col']
        self.header['Number_of_nonlinear_iteration']['unit'] = None
        self.header['total_col'] += 1
        return _data_list


class MACHINE_TIME_PLOT(LINEARPLOT):
    '''
    Class for plotting machine time file.
    This is an inheritage of the LINEARPLOT class

    Attributes:
    Args:
    '''
    def __init__(self, _name, **kwargs):
        LINEARPLOT.__init__(self, _name, kwargs)  # call init from base function

    def ManageData(self):
        '''
        manage data, get new data for this class
        for the base class, this method simply takes the combination
        of self.data
        Returns:
            _data_list(list):
                list of data for ploting
        '''
        _data_list = []
        for i in range(self.data.shape[1]):
            _data_list.append(self.data[:, i])

        # compute core time
        col_mt = self.header['Machine_time']['col']
        col_cpu = self.header['CPU_number']['col']
        core_time = self.data[:, col_mt] * self.data[:, col_cpu]

        # append core time
        _data_list.append(core_time)
        self.header['Core_time'] = {'col': self.header['total_col'], 'unit': None}
        self.header['total_col'] += 1

        return _data_list

    def GetStepMT(self, filename, step):
        '''
        get the total core hours spent before reaching a specific step
        Inputs:
            step(int): a specified step
        Returns:
            machine_time(float): machine time spent before reaching this step
        '''
        # Read header and data
        self.ReadHeader(filename)
        state=self.ReadData(filename)
        if state == 1:
            # empty file, return None
            return None

        # get steps and machine time
        col_step = self.header['Time_step_number']['col']
        col_mt = self.header['Machine_time']['col']
        steps = self.data[:, col_step]
        machine_times = self.data[:, col_mt]

        # check step is in range
        my_assert(type(step) == int, TypeError, "GetStepMT: step mush be an int value")
        my_assert(step <= np.max(steps), ValueError, "GetStepMT: step given is bigger than maximum step")

        # interp for step
        machine_time_at_step = np.interp(step, steps, machine_times)

        # compute total time
        col_cpu = self.header['CPU_number']['col']
        number_of_cpu = self.data[0, col_cpu]
        return machine_time_at_step * number_of_cpu, number_of_cpu


def ProjectPlot(case_dirs, _file_type, **kwargs):
    '''
    Plot figures for all cases in this project
    Inputs:
        kwargs:
            update(True or False): if True, update existing figures
    '''
    from shilofue.PlotDepthAverage import DEPTH_AVERAGE_PLOT

    update = kwargs.get('update', False)
    pdict = kwargs.get('pdict', {})
    # Init the UnitConvert class
    UnitConvert = Utilities.UNITCONVERT()

    # plot statistics ouput
    plot_options = pdict.get('Statistics', {})
    Statistics = STATISTICS_PLOT_OLD('Statistics', unit_convert=UnitConvert, options=plot_options)
    # depth average output
    plot_options = pdict.get('DepthAverage', {})
    DepthAverage = DEPTH_AVERAGE_PLOT('DepthAverage', unit_convert=UnitConvert, options=plot_options)
    # newton solver output
    # This is a little bit confusing to myself.
    # But yes, we need two of them to do different type of plot.
    plot_options = pdict.get('NewtonSolverStep', {})
    NewtonSolverStep = NEWTON_SOLVER_PLOT('NewtonSolverStep', options=plot_options)
    plot_options = pdict.get('NewtonSolver', {})
    NewtonSolver = NEWTON_SOLVER_PLOT('NewtonSolver', options=plot_options)
    # machine time
    plot_options = pdict.get('MachineTime', {})
    MachineTime = MACHINE_TIME_PLOT('MachineTime', options=plot_options)

    # loop for cases and post process
    for _case_dir in case_dirs:
        # cases
        _case_output_dir = os.path.join(_case_dir, 'output')
        _case_img_dir = os.path.join(_case_dir, 'img')
        if not os.path.isdir(_case_img_dir):
            # make img folder if not exists
            os.mkdir(_case_img_dir)

        # plot statistic
        _statistic_file = os.path.join(_case_output_dir, 'statistics')
        _ofile = os.path.join(_case_img_dir, 'Statistics.'+ _file_type)
        # compare the dates of files, determine whether to plot
        is_plot = False
        if os.path.isfile(_statistic_file):
            if not os.path.isfile(_ofile) or \
               os.stat(_statistic_file)[8] > os.stat(_ofile)[8]:
                is_plot = True
        # plot
        if is_plot:
            try:
                Statistics(_statistic_file, fileout=_ofile)
            except Exception as e:
                raise Exception("Plot statistics file failed for %s, please chech file content.\
One option is to delete incorrect file before running again" % _statistic_file) from e
            else:
                print('Plot has been generated: ', _ofile)  # screen output

        # plot depth-average
        _depth_average_file = os.path.join(_case_output_dir, 'depth_average.txt')
        _time = 0.0
        _ofile_route = os.path.join(_case_img_dir, 'DepthAverage.%s' % _file_type)
        _ofile = os.path.join(_case_img_dir, 'DepthAverage_t%.8e.%s' % (_time, _file_type))  # ofile has the exact time
        if os.path.isfile(_depth_average_file) and (not os.path.isfile(_ofile) or update is True):
            # check for ofile here is not precist, not intuitive. future: change the implementation
            try:
                _ofile_exact = DepthAverage(_depth_average_file, fileout=_ofile_route, time=_time)
            except Exception as e:
                raise Exception("Plot DepthAverage file failed for %s, please chech file content.\
One option is to delete incorrect file before running again" % _depth_average_file) from e
            else:
                if _ofile_exact is not None:
                    # output when there is file generated
                    print('Plot has been generated: ', _ofile_exact)  # screen output

        # add solver output
        # plot newton solver output
        _solver_file = os.path.join(_case_output_dir, 'solver_output')
        _ofile_route = os.path.join(_case_img_dir, 'NewtonSolverStep.%s' % _file_type)
        # plot step0
        _step = 0
        NewtonSolverStep.GetStep(_step)
        _ofile = os.path.join(_case_img_dir, 'NewtonSolverStep_s%07d.%s' % (_step, _file_type))
        if os.path.isfile(_solver_file) and (not os.path.isfile(_ofile) or update is True):
            # check for ofile here is not precist, not intuitive. future: change the implementation
            try:
                _ofile_exact = NewtonSolverStep(_solver_file, fileout=_ofile_route)
            except Exception as e:
                raise Exception("Plot NewtonSolver file failed for %s, please chech file content.\
One option is to delete incorrect file before running again" % _solver_file) from e
            else:
                if _ofile_exact is not None:
                    # output when there is file generated
                    print('Plot has been generated: ', _ofile_exact)  # screen output
        # plot step 1
        _step = 1
        NewtonSolverStep.GetStep(_step)
        _ofile = os.path.join(_case_img_dir, 'NewtonSolverStep_s%07d.%s' % (_step, _file_type))
        if os.path.isfile(_solver_file) and (not os.path.isfile(_ofile) or update is True):
            # check for ofile here is not precist, not intuitive. future: change the implementation
            try:
                _ofile_exact = NewtonSolverStep(_solver_file, fileout=_ofile_route)
            except Exception as e:
                raise Exception("Plot NewtonSolver file failed for %s, please chech file content.\
One option is to delete incorrect file before running again" % _solver_file) from e
            else:
                if _ofile_exact is not None:
                    # output when there is file generated
                    print('Plot has been generated: ', _ofile_exact)  # screen output
        # plot whole history
        _ofile = os.path.join(_case_img_dir, 'NewtonSolver.%s' % _file_type)
        if os.path.isfile(_solver_file) and (not os.path.isfile(_ofile) or update is True):
            # check for ofile here is not precist, not intuitive. future: change the implementation
            try:
                _ofile_exact = NewtonSolver(_solver_file, fileout=_ofile)
            except Exception as e:
                raise Exception("Plot NewtonSolver file failed for %s, please chech file content.\
One option is to delete incorrect file before running again" % _solver_file) from e
            else:
                if _ofile_exact is not None:
                    # output when there is file generated
                    print('Plot has been generated: ', _ofile_exact)  # screen output

        # plot machine_time
        _machine_time_file = os.path.join(_case_output_dir, 'machine_time')
        _time = 0.0
        _ofile = os.path.join(_case_img_dir, 'MachineTime.%s' % _file_type)  # ofile has the exact time
        # compare the dates of files, determine whether to plot
        is_plot = False
        if os.path.isfile(_machine_time_file):
            if not os.path.isfile(_ofile) or \
               os.stat(_machine_time_file)[8] > os.stat(_ofile)[8]:
                is_plot = True
        # plot
        if is_plot:
            # check for ofile here is not precist, not intuitive. future: change the implementation
            try:
                _ofile_exact = MachineTime(_machine_time_file, fileout=_ofile)
            except Exception as e:
                raise Exception("Plot MachineTime file failed for %s, please chech file content.\
One option is to delete incorrect file before running again" % _machine_time_file) from e
            else:
                if _ofile_exact is not None:
                    # output when there is file generated
                    print('Plot has been generated: ', _ofile_exact)  # screen output

    pass


# a main function
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
    parser.add_argument('-j', '--json_file', type=str,
                        default='./config_case.json',
                        help='Filename for json file')
    parser.add_argument('-i', '--inputs', type=str,
                        default='',
                        help='A directory that contains the input')
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    # commands

    if _commend == 'analyze_affinity_test_results':
        # example:
        # python -m shilofue.Plot analyze_affinity_test_results
        # -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/rene_affinity_test/results/spherical_shell_expensive_solver/peloton-ii-32tasks-core-openmpi-4.0.1
        analyze_affinity_test_results(arg.inputs)

# run script
if __name__ == '__main__':
    main()
