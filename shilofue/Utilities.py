import json
import re
import os
import shilofue.json
import numpy as np
from importlib import resources


def JsonOptions(prefix):
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
            with resources.open_text(shilofue.json, 'UnitConvert.json') as fin:
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
                raise KeyError('unit_from: %s is neither a recorded unit or an alias for a unit' % _unit_from)
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
                raise KeyError('unit_to is neither a recorded unit or an alias for a unit')
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
    todo: add test function
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
    todo: use system ts for tab space
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