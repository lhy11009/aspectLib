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
import sys, os, argparse, json, re, warnings
import numpy as np
from matplotlib import pyplot as plt
current_dir = os.path.dirname(__file__)
JSON_FILE_DIR = os.path.join(current_dir, 'json_examples')
sys.path.append(os.path.join(current_dir, 'utilities', 'python_scripts'))
import Utilities
class sPlot_LINEARPLOT():
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
        self.options = Utilities.JsonOptions(_name, JSON_FILE_DIR)
        self.UnitConvert = options.get('unit_convert', None)
        self.dim = options.get('dim', 2)  # dimension
        assert(self.dim in [1, 2, 3])  # dimension must be 1, 2, 3

        # reset the options with a option in the options
        json_file_path = os.path.join(JSON_FILE_DIR, "post_process_std.json")
        with open(json_file_path, 'r') as fin:
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
            Utilities.my_assert((type(rows) == list), TypeError, "%s: rows must be a list" % Utilities.func_name())
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


class sPlot_LINEARPLOT():
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
        self.options = Utilities.JsonOptions(_name, JSON_FILE_DIR)
        self.UnitConvert = options.get('unit_convert', None)
        self.dim = options.get('dim', 2)  # dimension
        assert(self.dim in [1, 2, 3])  # dimension must be 1, 2, 3

        # reset the options with a option in the options
        json_file_path = os.path.join(JSON_FILE_DIR, "post_process_std.json")
        with open(json_file_path, 'r') as fin:
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
            Utilities.my_assert((type(rows) == list), TypeError, "%s: rows must be a list" % Utilities.func_name())
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


class sPlot_LINEARPLOT():
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
        self.options = Utilities.JsonOptions(_name, JSON_FILE_DIR)
        self.UnitConvert = options.get('unit_convert', None)
        self.dim = options.get('dim', 2)  # dimension
        assert(self.dim in [1, 2, 3])  # dimension must be 1, 2, 3

        # reset the options with a option in the options
        json_file_path = os.path.join(JSON_FILE_DIR, "post_process_std.json")
        with open(json_file_path, 'r') as fin:
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
            Utilities.my_assert((type(rows) == list), TypeError, "%s: rows must be a list" % Utilities.func_name())
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


class STATISTICS_PLOT(sPlot_LINEARPLOT):
    '''
    Class for plotting depth average file.
    This is an inheritage of the LINEARPLOT class

    Attributes:
    Args:
    '''
    def __init__(self, _name, **kwargs):
        sPlot_LINEARPLOT.__init__(self, _name, kwargs)  # call init from base function
    
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
    
    def PlotNumberOfIterationsStokessolver(self, **kwargs):
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
        col_iter_stokes = self.header['Iterations_for_Stokes_solver']['col']
        times = self.data[:, col_t]
        steps = self.data[:, col_step]
        iter_stokes = self.data[:, col_iter_stokes]
        if self.UnitConvert is not None:
            to_myr = self.UnitConvert(unit_t, 'myr')
        else:
            raise ValueError("a UNITCONVERT class must be given")
        ax.plot(times * to_myr, iter_stokes, '.', label=label, color=color)
        ax.set_xlabel('Time (myr)')
        ax.set_ylabel('Number of Iterations for Stokes solver')
    
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
    json_file = os.path.join(JSON_FILE_DIR, 'post_process.json')
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