import json
import re
import os
import inspect
import shilofue.json_files
import numpy as np
from importlib import resources
from pathlib import Path

class CODESUB():
    '''
    This is a class to substitute existing code with values & parameters
    Attributes:
        contents (str)
        options(disc)
    '''
    def __init__(self):
        contents=''
        options={}

    def read_contents(self, _path):
        '''
        read contents from a file
        '''
        my_assert(os.access(_path, os.R_OK), FileNotFoundError, "%s: %s cannot be opened" % (func_name(), _path))
        with open(_path, 'r') as fin:
            self.contents = fin.read()

    def read_options(self, _path):
        '''
        read options from a json file
        '''
        my_assert(os.access(_path, os.R_OK), FileNotFoundError, "%s: %s cannot be opened" % (func_name(), _path))
        with open(_path, 'r') as fin:
            self.options = json.load(fin)

    def substitute(self):
        '''
        substitute keys with values
        '''
        for key, value in self.options.items():
            self.contents = re.sub(key, str(value), self.contents)

    def save(self, _path):
        '''
        save contents to a new file
        '''
        # look for directory
        dir_path = os.path.dirname(_path)
        if not os.path.isdir(dir_path):
            os.mkdir(dir_path)
        # save file
        with open(_path, 'w') as fout:
            fout.write(self.contents)
        return _path

def JsonOptions(prefix, _dir=None):
    '''
    JsonOptions(prefix)
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
    if _dir == None:
        # default option
        for _filename in resources.contents(shilofue.json_files):
            # this resources.contents return a interable from a sub-package,
            # entries in this interable are filenames
            if re.match("^" + prefix + '_', _filename):
                # the following two lines eliminate the words before the first '_'
                # with that '_' as well as the '.json' in the end.
                # so that if filename is 'Statistics_Number_of_Cells.json',
                # _name is 'Number_of_Cells'
                _name = _filename.split('_', maxsplit=1)[1]
                _name = _name.rsplit(".", maxsplit=1)[0]
                with resources.open_text(shilofue.json_files, _filename) as fin:
                    _options[_name] = json.load(fin)  # values are entries in this file
    else:
        pathlist = Path(_dir).rglob('%s_*.json' % prefix)
        for path in pathlist:
            _filename = str(path)
            _base_name = os.path.basename(_filename)
            _name = _base_name.split('_', maxsplit=1)[1]
            _name = _name.rsplit(".", maxsplit=1)[0]
            with open(_filename, 'r') as fin:
                _options[_name] = json.load(fin)  # values are entries in this file
    assert(_options is not {})  # assert not vacant
    return _options


def ReadHeader(_texts):
    '''
    Read header information from file.
    An example of string is:
    '# 1: Time (years)'
    Args:
        _texts(list<string>): a list of string, each is a line of header
    Returns:
        _header(dict): header information
            key(dict): infomation under this key:
                'col':
                    column in file
                'unit':
                    unit of this variable
    '''
    _header = {}
    _header['total_col'] = 0  # total columes in file
    for _line in _texts:
        # derive column, unit and key from header
        if re.match("^#", _line) is None:
            # only look at the first few lines starting with "#"
            break
        _line = re.sub("^# ", '', _line)  # eliminate '#' in front
        # match the first number, which is the columne in file
        _match_obj = re.match("[0-9]*", _line)
        _col = int(_match_obj[0]) - 1
        _line = re.sub("^.*?: ", '', _line) # eliminate words before ':'
        # match string as '(foo)', thus get the unit
        _match_obj = re.search("[(].*[)]", _line)
        if _match_obj:
            # if matched, get the unit
            _unit = _match_obj.group(0)
            _unit = _unit[1: -1] # _unit[0] is ' ', thus not included
        else:
            _unit = None
        # eliminate '(foo)', what left is the name
        _line = re.sub(" [(].*[)]", '', _line)
        _line = _line.strip("\n")
        _key = re.sub(" ", "_", _line)  # substitute ' ' with '_'
        # save information as key, options
        # options include col and unit
        _info = {}
        _info['col'] = _col
        _header['total_col'] = max(_header['total_col'], _col+1)
        _info['unit'] = _unit
        _header[_key] = _info
    return _header


def ReadHeader2(_texts):
    '''
    Read header information from file.
    An example of string is:
    '# Time Depth Viscosity Temperature'
    Args:
        _texts(list<string>): a list of string, each is a line of header, for
        this method, there should be only one line
    Returns:
        _header(dict): header information
            key(dict): infomation under this key:
                'col':
                    column in file
                'unit':
                    unit of this variable
    '''
    _header = {}
    _header['total_col'] = 0  # total columes in file
    # derive column, unit and key from header
    _line = _texts[0]
    _line = re.sub("^# ", '', _line)  # eliminate '#' in front
    _line_segments = re.split('( |\t|\n)+', _line)  # split by words
    for _key in _line_segments:
        if _key not in ['', ' ', '\t', '\n']:
            _header[_key] = {} # initial dictionary for this key
            _header[_key]['col'] = _header['total_col']  # assign column number
            _header[_key]['unit'] = None
            _header['total_col'] += 1
    return _header


def my_assert(_condition, _errortype, _message):
    '''
    an assert function for runtime use
    Inputs:
        _condition(True or False):
            the condition to assert
        _errortype(some error):
            the type of error to raise
        _message(str):
            the message to raise
    '''
    if _condition is False:
        raise _errortype(_message)


class WarningTypes():
    '''
    A class of self defined warning types
    '''
    class FileHasNoContentWarning(Warning):
        pass


class UNITCONVERT():
    '''
    UNITCONVERT():
    class for unit convert

    Attributes:
        units(dict): magnitudes of units, values are [int, str].
            These are magnitude of this unit and the base unit for
            this unit.
        alias(dict): alias of units, values are string

    Functions:
        __init__:
            initiation
        __call__:
            Call function, return the convertion ratio
    '''


    def __init__(self, **kwargs):
        '''
        initiation
        Args:
            kwargs:
                'filename': filename of the file to load. If not exiting, use
                the default file  'UnitConvert.json' in shilofue.json instead
        '''
        _filename = kwargs.get('filename', None)
        if _filename is None:
            with resources.open_text(shilofue.json_files, 'UnitConvert.json') as fin:
                _data = json.load(fin)
        else:
            with open(_filename, 'r') as fin:
                _data = json.load(fin)
        # a dictionary, this is the magnitudes of units, every value is an list
        # of an int and a string
        self.units = _data['units']
        # a dictionary, this is the alias of units, every value is a string
        self.alias = _data['alias']

    def __call__(self, _unit_from, _unit_to):
        '''
        Call function, return the convertion from _unit_from to _unit_to so that
        Args:
            _unit_from(str):
                unit to convert from
            _unit_to(str):
                unit to convert to
        Raises:
            KeyError:
                when unit is neither a recorded unit or an alias for a unit
        Assertions:
            assert that the base of these two are the same
        returns:
            _convert_ratio(float):
                a ratio of magnitude of two units
        '''
        # units are either a recorded unit or an alias for a recorded unit
        try:
            _unit_from_magnitude = self.units[_unit_from][0]
            _unit_from_base = self.units[_unit_from][1]
        except KeyError:
            try:
                _unit_from_alias = self.alias[_unit_from]
            except KeyError:
                raise KeyError('unit_from(i.e. %s) is neither a recorded unit or an alias for a unit' % _unit_from)
            try:
                _unit_from_magnitude = self.units[_unit_from_alias][0]
            except KeyError:
                raise KeyError('The alias %s is not registered' % _unit_from_alias)
            _unit_from_base = self.units[_unit_from_alias][1]
       # same thing for unit_to
        try:
            _unit_to_magnitude = self.units[_unit_to][0]
            _unit_to_base = self.units[_unit_to][1]
        except KeyError:
            try:
                _unit_to_alias = self.alias[_unit_to]
            except KeyError:
                raise KeyError('unit_to(i.e %s) is neither a recorded unit or an alias for a unit' % _unit_to)
            try:
                _unit_to_magnitude = self.units[_unit_to_alias][0]
            except KeyError:
                raise KeyError('The alias %s is not registered' % _unit_to_alias)
            _unit_to_base = self.units[_unit_to_alias][1]
        # convert unit
        assert(_unit_from_base == _unit_to_base)  # assert that the base of these two are the same
        _convert_ratio = _unit_from_magnitude / _unit_to_magnitude  # ratio of the two magnitude
        return _convert_ratio


# re functions
def re_neat_word(_pattern):
    '''
    Eliminate ' ', '\\t', '\\n' in front of and at the back of _pattern
    future: add test function
    Inputs:
        _pattern(str): input
    Outputs:
        _neat(str): result
    '''
    _neat = re.sub('^[ \t]*', '', _pattern)
    _neat = re.sub('[ \t\n]*$', '', _neat)  # strip value
    return _neat


def re_count_indent(_pattern):
    '''
    conunt indentation at the start of a string
    future: use system ts for tab space
    Inputs:
        _pattern(str): input
    Outputs:
        _indent(int): indentation at the start of the string
    '''
    assert(type(_pattern) is str)
    _indent = 0
    for i in range(len(_pattern)):
        if _pattern[i] == ' ':
            _indent += 1
        elif _pattern[i] == '\t':
            _indent += 4
        else:
            break
    return _indent


def ggr2cart(lat,lon,r):
    # transform spherical lat,lon,r geographical coordinates
    # to global cartesian xyz coordinates
    #
    # input:  lat,lon,r in radians, meters
    # output: x,y,z in meters 3 x M
    sla = np.sin(lat)
    cla = np.cos(lat)
    slo = np.sin(lon)
    clo = np.cos(lon)

    x = r * cla * clo
    y = r * cla * slo
    z = r * sla

    return x,y,z

def ggr2cart2(lon, r):
    # transform spherical lon, r geographical coordinates
    # to global cartesian xy coordinates
    #
    # input:  phi, r in radians, meters
    # output: x,y in meters 2 x M
    slo = np.sin(lon)
    clo = np.cos(lon)

    x = r * clo
    y = r * slo

    return x, y


def cart2sph(x,y,z):
    """
    A function that transfers cartisen geometry to spherical geometry
    """
    r = (x**2+y**2+z**2)**0.5
    th = np.arctan2((x**2+y**2)**0.5,z)
    ph = np.arctan2(y,x)
    return r, th, ph


def cart2sph2(x,y):
    """
    A function that transfers cartisen geometry to spherical geometry in 2d
    """
    r = (x**2+y**2)**0.5
    ph = np.arctan2(y,x)
    return r, ph


def Make2dArray(x):
    """
    A function that returns 2d nparray, this could be useful to generate uniform output
    """
    if type(x) in [int, float]:
        array_ = Make2dArray([x])
    elif type(x) is list:
        array_ = np.array(x)
        array_ = array_.reshape((1, -1))
    elif type(x) is np.array and x.ndix == 1:
        array_ = x.reshape((1, -1))
    else:
        raise TypeError("Make2dArray, x must be int, float, list or 1d ndarray")
    return array_


def WriteFileHeader(ofile, header):
    """
    A function that writes statistic-like file header
    Inputs:
        ofile(str): file to output
        header(dict): dictionary of header to output
    """
    # assert file not existing
    with open(ofile, 'w') as fout:
        output = ''
        for key, value in header.items():
            col = value['col']
            unit = value.get('unit', None)
            output += '# %d: %s' % (col+1, key)
            if unit is not None:
                output += ' (%s)' % unit
            output += '\n'
        fout.write(output)


def PhaseFunction(x):
    """
    function defined in MaterialModel/Utilities in aspect
    """
    return 0.5 * (1 + np.tanh(x))


def AveragePhaseFunctionInputs(x1, x2):
    """
    todo
    Average phase function value
    """
    my_assert(x1.shape == x2.shape, ValueError, "Inputs(x1 and x2) need to be arrays of the same shape")
    average = np.zeros(x1.shape)

    # logical expressions for regions
    pin_point = 2.0
    limit_point = -2.0
    is_above = np.logical_and(x1 > pin_point, x2 > pin_point)
    is_upper_left = np.logical_and(x1 < pin_point, x2 > pin_point)
    is_lower_right = np.logical_and(x1 > pin_point, x2 < pin_point)
    is_middle = np.logical_and(x1 < pin_point, x2 < pin_point)
    # compute values
    average[is_above] = np.minimum(x1[is_above], x2[is_above])
    average[is_upper_left] = x1[is_upper_left]
    average[is_lower_right] = x2[is_lower_right]
    average[is_middle] = pin_point - np.sqrt((pin_point - x1[is_middle])**2.0 + (pin_point - x2[is_middle])**2.0)
    return average


def touch(path):
    """
    touch a file as in unix system
    """
    with open(path, 'a'):
        os.utime(path, None)


def func_name():
    """
    return the name of the calling function
    """
    return inspect.currentframe().f_back.f_code.co_name


def dump_message(fout, message):
    """
    dump_message to a log file
    notes: I haven't used this since I am not sure whether it's good to open up a log file at the beginning of a script.
    Inputs:
        fout: an object that denines an output stream (e.g. sys.stderr)
    """
    outputs = "%s: %s" % (inspect.currentframe().f_back.f_code.co_name, message)
    fout.write(outputs)


def get_name_and_extention(_path):
    """
    return name and extention
    """
    _name = _path.rsplit(".", maxsplit=1)[0]
    _extension = _path.rsplit(".", maxsplit=1)[1]
    return _name, _extension
