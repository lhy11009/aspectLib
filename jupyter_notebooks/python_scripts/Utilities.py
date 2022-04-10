import json
import re
import os
import inspect
import numpy as np
from importlib import resources
from pathlib import Path
from PIL import Image, ImageFilter, ImageFont, ImageDraw


r'''
Alias for longer syntax in python
'''
def func_name():
    """
    return the name of the calling function
    """
    return inspect.currentframe().f_back.f_code.co_name


r'''
Functions that bridging unix operations
'''
def touch(path):
    """
    touch a file as in unix system
    """
    with open(path, 'a'):
        os.utime(path, None)


def var_subs(_str, **kwargs):
    """
    substitute variable as if in bash
    _str should looks like:
        "${HOME}/foo", in this case, HOME will be substitued with the
        corresponding environmental variable.
    """
    assert(type(_str)==str)
    parts = _str.split('${')
    result = ''  # initiation
    i = 0  # count
    for part in parts:
        if i >= 1:
            sub_parts = part.split('}', maxsplit=1)  # extract the name of the variable
            assert(len(sub_parts) == 2)
            part_substituted = os.environ[sub_parts[0]] + sub_parts[1]  # perform substitution
        else:
            part_substituted = part
        result += part_substituted
        i += 1
    return result
    

r'''
functions for Execption handling
'''

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


r'''
Functions for string options
'''
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


def re_read_variable_from_string(inputs, pattern, splitter):
    '''
    read variable from a input string, typically some output from another executable
    Inputs:
        inputs(str)
    '''
    found = False
    input_array = inputs.split('\n')
    for line in input_array:
        if re.match(pattern, line):
            output = line.split(splitter)[1]
            output = re_neat_word(output)
            found = True
    if not found:
        raise ValueError("%s: pattern(%s) not found" % (func_name(), pattern))
    return output


def get_name_and_extention(_path):
    """
    return name and extention
    """
    _name = _path.rsplit(".", maxsplit=1)[0]
    _extension = _path.rsplit(".", maxsplit=1)[1]
    return _name, _extension


def string2list(inputs, _type=int):
    """
    convert list input to string outputs
    """
    inputs = inputs.strip(']')
    inputs = inputs.strip('[')
    inputs = inputs.strip(' ')
    outputs_str = inputs.split(',')
    if _type == int:
        outputs = [int(i) for i in outputs_str]
    elif _type == float:
        outputs = [float(i) for i in outputs_str]
    else:
        raise ValueError("%s: invalid option of _type (int or float)" % func_name())
    return outputs


r'''
Functions for read files
'''
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


r'''
Functions for writing files
'''
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


def dump_message(fout, message):
    """
    dump_message to a log file
    notes: I haven't used this since I am not sure whether it's good to open up a log file at the beginning of a script.
    Inputs:
        fout: an object that denines an output stream (e.g. sys.stderr)
    """
    outputs = "%s: %s" % (inspect.currentframe().f_back.f_code.co_name, message)
    fout.write(outputs)


'''
Functions for substributing contents of file
'''
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

    def read_contents(self, *paths):
        '''
        read contents from a file
        '''
        self.contents=''
        i = 0 # conunt
        for _path in paths:
            my_assert(os.access(_path, os.R_OK), FileNotFoundError, "%s: %s cannot be opened" % (func_name(), _path))
            with open(_path, 'r') as fin:
                if i > 0:
                    self.contents += '\n\n'
                self.contents += fin.read()
            i += 1

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


r'''
Functions for coordinate transform
'''
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


r'''
Functions for additional opionts in numpy
'''
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



r'''
Functions for phase functions
'''
def PhaseFunction(x):
    """
    function defined in MaterialModel/Utilities in aspect
    """
    return 0.5 * (1 + np.tanh(x))


def AveragePhaseFunctionInputs(x1, x2):
    """
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


r"""
Classes and functions related to parsing json file as options
"""
class JSON_OPT():
    """
    This class defines json options and offer interfaces to use
    and document them
    Args:
        keys (list of list of str): each member contains a list of keys from the top level
        descriptions (list of str): each member is a description of the option of the file
        types (list of type): each member is the type of the variable
        values (list of variable (type)): each member is the value of a variable
        defaults (list of variable (type)): each member is the default value of a variable
        nicks (list of variable (type)): each member is the nickname of a variable
        start (int): number of entries in the parental class
    """
    def __init__(self):
        """
        Initiation
        """
        self.keys = []
        self.descriptions = []
        self.values = [] 
        self.types = []
        self.defaults = []
        self.nicks = []
        self.start = 0


    def read_json(self, _path):
        """
        Read in json option from file path
        Inputs:
            _path (str): path of the json file
        """
        assert(os.access(_path, os.R_OK))
        print("    Read options from json file: %s" % _path)
        with open(_path, 'r') as fin:
            options = json.load(fin)
        self.import_options(options)

    def import_options(self, options):
        """
        Import options from a dictionary and then check the values
        Inputs:
            options (dict): dictionary of options
        """
        assert(type(options)==dict)
        for i in range(len(self.keys)):
            if type(self.types[i]) == tuple:
                # this is a feature
                try:
                    feature_options = read_dict_recursive(options, self.keys[i])
                except KeyError:
                    pass
                else:
                    # note that self.types[i][1] is the class of sub-opitions
                    assert(type(feature_options) == list)
                    for feature_option in feature_options:
                        self.values[i].append(\
                            create_option_object(self.types[i][1], feature_option)\
                                )
            else:
                # this is a variable
                try:
                    value = read_dict_recursive(options, self.keys[i])
                    self.values[i] = value
                    my_assert(type(value) == self.types[i], TypeError,\
                    "%s: type of the default (%s) is not %s" % (func_name(), str(type(value)), str(self.types[i]))) # assert the type of value
                except KeyError:
                    pass
        self.check()

    def check(self):
        """
        check to see if these values make sense
        """
        pass

    def add_key(self, description, _type, keys, default_value, **kwargs):
        """
        Add an option (key, describtion)
        Inputs:
            description (str)
            _type (type)
            keys (list): a list of keys from the top level
            kwargs (dict):
                nick (str): nickname
        """
        for key in keys:
            assert(type(key) == str)
        nick = kwargs.get('nick', keys[-1])
        self.keys.append(keys)
        self.descriptions.append(description)
        self.types.append(_type)
        my_assert(_type == type(default_value), TypeError,\
        "%s: type of the default (%s) is not %s" % (func_name(), str(type(default_value)), str(_type))) # assert the type of default
        self.defaults.append(default_value)
        self.values.append(default_value)
        self.nicks.append(nick)
        pass

    def add_features(self, description, keys, SUB_OPT, **kwargs):
        """
        add a feature, which shows as a list in the json file.
        Inputs:
            description (str)
            keys (list): a list of keys from the top level
            SUB_OPT: a sub-class of JSON_OPT for storing options of
                individual features
            kwargs (dict):
                nick (str): nickname
        """
        nick = kwargs.get('nick', keys[-1])
        self.keys.append(keys)
        self.values.append([])  # initiate a vacant list
        self.descriptions.append(description)
        self.types.append(("feature", SUB_OPT))
        self.defaults.append(None) # this is not used
        self.nicks.append(nick)

    def get_value(self, keys):
        """
        Get the value through a list of keys
        Inputs:
            keys(list): a list of keys from the top level
        """
        pass

    def print(self, indent=0):
        """
        Print the keys and related values stored
        Inputs:
            indent(int): indentation at the front
        Returns:
            _str(str): string output
        """
        _str = indent*' ' + "All options stored in this object:\n"
        for i in range(len(self.keys)):
            _str += '\n' + (indent+4)*' ' + str(self.keys[i]) + '\n'
            _str += (indent+4)*' ' + "Description: %s" % self.descriptions[i] + '\n'
            _str += (indent+4)*' ' + "Value: %s" % str(self.values[i]) + '\n'
        return _str

    def document(self, indent=0, **kwargs):
        """
        Print the documentation of this class.
        Inputs:
            indent(int): indentation at the front
            kwargs (dict):
                start : from where to start. This is a option to print keys in parental
                    classes and daughter classes separately
        Returns:
            _str(str): string output
        """
        _start = kwargs.get('start', 0)
        _str = '\n' + indent*' ' + '(Note the first number means the index, \
while the second (in brackets) is the index relative to the parental class)'
        for i in range(len(self.keys)):
            if i > _start:
                _str += '\n'+ (indent+4)*' ' + '%d (%d):' % (i, i - _start) + '\n'
            else:
                _str += '\n'+ (indent+4)*' ' + '%d :' % i + '\n'
            _str += (indent+4)*' ' + str(self.keys[i]) + '\n'
            _str += (indent+8)*' ' + "Nickname: %s" % self.nicks[i] + '\n'
            _str += (indent+8)*' ' + "Description: %s" % self.descriptions[i] + '\n'
            _str += (indent+8)*' ' + "Default value: %s" % str(self.defaults[i]) + '\n'
            if type(self.types[i]) is tuple:
                # this is a feature
                sub_opt = self.types[i][1]()
                _str += sub_opt.document(indent+4)
        return _str
    
    def document_str(self):
        '''
        Wrapper for document, call this to properly layout documentation for daughter and parental classes
        Return:
            document string (str)
        '''
        return self.document(start=self.start)

    def number_of_keys(self):
        """
        Returns:
            n_keys (int): number of keys defined
        """
        return len(self.keys)
    
    def __call__(self, func_name):
        """
        Interface to some function.
        This should return the exact outputs that function need
        Implement in daughter classes
        """
        pass


def create_option_object(OPT, options):
    """
    create an object with the intended options    
    Inputs:
        OPT (class inherits JSON_OPT)
        options (dict)
    Return:
        obt (an instantiation of OPT) - object containing the right options
    """
    assert(type(options)==dict)
    obt = OPT()
    obt.import_options(options)
    return obt


def show_all_options(_path):
    '''
    Show all options in a json file as list of keys
    This aims to aid the process of interpreting and coding
    the options of this kind of file.
    Inputs:
        _path (str): path of the input json file
    Returns:
        -
    '''
    assert(os.access(_path, os.R_OK))
    with open(_path, 'r') as fin:
        json_contents = json.load(fin)
    all_options = show_all_options_dict(json_contents)
    return all_options


def show_all_options_dict(_dict):
    '''
    Show all options in a dict type as list of keys
    This aims to aid the process of interpreting and coding
    the options of this kind of file.
    Inputs:
        _dict (dict): input dictionary contains options
    Returns:
        all_options (list): all options put in a list with each entry
        being a list of keys.
    '''   
    all_options = []
    for key, value in _dict.items():
        if re.match("^_", key) is None:  # skip keys start with '_'
            all_options.append([key])
            if type(value) == dict:
                sub_options = show_all_options_dict(value)  # recursive call function with nested options
                for _option in sub_options:
                    all_options.append([key] + _option)  # append to the upper level key
    return all_options


def read_dict_recursive(_dict, list_of_keys):
    '''
    Read value by a list of keys
    i.e. ['visit', 'slab'], and 'visit' is a subdictionary within
    the given _dict variable.
    Inputs:
        _dict (dict): input dictionary
        list_of_keys (str): list of keys to look for
    '''
    if len(list_of_keys) > 1:
        sub_dict = _dict[list_of_keys[0]]
        sublist_of_keys = [list_of_keys[i] for i in range(1, len(list_of_keys))]
        value = read_dict_recursive(sub_dict, sublist_of_keys)
    else:
        value = _dict[list_of_keys[0]]
    return value


def write_dict_recursive(_dict, list_of_keys, value):
    '''
    Write value in a dict by a list of keys
    i.e. ['visit', 'slab'], and 'visit' is a subdictionary within
    the given _dict variable.
    Inputs:
        _dict (dict): input dictionary
        list_of_keys (str): list of keys to look for
    '''
    if len(list_of_keys) > 1:
        sub_dict = _dict[list_of_keys[0]]
        sublist_of_keys = [list_of_keys[i] for i in range(1, len(list_of_keys))]
        sub_dict = write_dict_recursive(sub_dict, sublist_of_keys, value)
    else:
        _dict[list_of_keys[0]] = value
    return _dict


r"""
Classes and functions related to post-process images
"""

class IMAGE_OPT(JSON_OPT):
    '''
    class to hold options for individual images
    '''
    def __init__(self):
        """
        Initiation
        """
        JSON_OPT.__init__(self)
        self.add_key("Path to the file", str, ["path"], "foo.png", nick='image_path')
        self.add_key("Operation to do", str, ["operation"], "new", nick='operation')
        self.add_key("Whether to create a mask (0 or 1)", int, ["mask"], 0, nick='mask')
        self.add_key("A scale to resize", float, ["resize"], 1.0, nick='resize')
        self.add_key("A position to put on the new figure", list, ["position"], [0, 0], nick='position')
        self.add_key("A path to save the figure", str, ["save path"], "", nick='save_path')
        self.add_key("Method to use", str, ["method"], "on_first_figure", nick='method')
        self.add_key("Whether this is a intermediate result to be removed",\
             int, ["temp"], 0, nick='is_temp')
        self.add_key("Text to use, only work with the \"text\" operation", str, ["text"], "", nick="text")
        self.add_key("Size of the Text to use, only work with the \"text\" operation", int, ["font size"], 40, nick="font size")
    
    def check(self):
        """
        check values
        """
        # my_assert(os.path.isfile(var_subs(self.values[0])), FileExistsError,\
        #    "%s: file %s doesn't exist" % (func_name(), self.values[0]))
        assert(self.values[1] in ["new", "paste", "crop", "text"])
        if self.values[1] == "paste":
            assert(len(self.values[4]) == 2)  # assert positon has the correct length
        elif self.values[1] == "crop":
            assert(len(self.values[4]) == 4)
        elif self.values[1] == "text":
            assert(type(self.values[8]) == str and self.values[8] != "")  # check text entry
            assert(type(self.values[9]) == int and self.values[9] > 0)
        
        assert(self.values[2] in [0, 1])  # mask is a bool value
        for i in self.values[4]:
            assert(type(i) == int)
        assert(self.values[7] in [0, 1])  # is_temp is a bool value
        if self.values[5] == "":
            assert(self.values[7] == 0)  # if not saving, no mark on intermediate result
    
    def to_pillow_run(self):
        '''
        interface to the PillowRun function
        '''
        im_path = var_subs(self.values[0])
        operation = self.values[1]
        resize = self.values[3]
        position = self.values[4]
        mask = self.values[2]
        method = self.values[6]
        save = var_subs(self.values[5])
        is_temp = self.values[7]
        text = self.values[8]
        font_size = self.values[9]
        return im_path, operation, resize, position, mask, method, save, is_temp, text, font_size


class PILLOW_OPT(JSON_OPT):
    '''
    class to hold options for using pillow to operation on images
    '''
    def __init__(self):
        """
        Initiation
        """
        JSON_OPT.__init__(self)
        self.add_features("Figures to operate with", ["figures"], IMAGE_OPT)
        pass

    def to_pillow_run(self):
        '''
        interface to the PillowRun function
        '''
        im_paths = []
        operations = []
        resizes = []
        positions = []
        masks = []
        methods = []
        saves = []
        is_temps = []
        texts = []
        font_sizes = []
        for feature in self.values[0]:
            im_path, operation, resize, position, mask, method, save, is_temp, text, font_size =\
                feature.to_pillow_run()
            im_paths.append(im_path)
            operations.append(operation)
            resizes.append(resize)
            positions.append(position)
            masks.append(mask)
            methods.append(method)
            saves.append(save)
            is_temps.append(is_temp)
            texts.append(text)
            font_sizes.append(font_size)
        return im_paths, operations, resizes, positions, masks, methods, saves, is_temps, texts, font_sizes


def ImageMerge(im_paths, **kwargs):
    '''
    Merge two or more plots
    Inputs:
        im_paths (list): list of path of image files
        positions (list of tuple): position on the new figure
        kwargs (dict):
            method: which method to use
                on_first_figure: take the first figure and lay everything on that
                use_new_one: take a new blank figure and lay everything on that
                masks: transparency mask
    '''
    length = len(im_paths)
    for im_path in im_paths:
        assert(os.path.isfile(im_path))  # assert file paths
    method = kwargs.get('method', 'on_first_figure')  # method
    masks = kwargs.get('mask', [0 for i in range(length)]) # transparency
    resizes = kwargs.get('resize', [1.0 for i in range(length)])
    positions = kwargs.get('position', [(0, 0) for i in range(length)])

    new_image = ImageMerge0(im_paths, method, masks, resizes, positions)
    return new_image


def ImageMerge0(im_paths, method, masks, resizes, positions):
    '''
    Merge two or more plots
    '''
    length = len(im_paths)
    assert(len(masks) == length)
    assert(len(resizes) == length)
    image0 = Image.open(im_paths[0])
    # get new size
    if method == 'on_first_figure':
        new_size = (image0.size[0], image0.size[1])
    elif method == 'use_new_one':
        new_size = [image0.size[0], image0.size[1]]
        for i in range(1, length):
            _image =  Image.open(im_paths[i])
            size = [ int(v * resizes[i]) for v in _image.size ]
        new_size[0] = max(positions[i][0]+size[0], new_size[0])        
        new_size[1] = max(positions[i][1]+size[1], new_size[1])        
    else:
        raise ValueError('method must be either \"on first figure\" or \"use_new_one\"')
    new_image = Image.new('RGB',new_size,(250,250,250))
    new_image.paste(image0, [0, 0])  # paste first figure
    for i in range(1, length):
        _image =  Image.open(im_paths[i])
        size = [ int(v * resizes[i]) for v in _image.size ]  # resize
        _image = _image.resize(size)
        if masks[i]:
            new_image.paste(_image, positions[i], mask=_image)
        else:
            new_image.paste(_image, positions[i])
    return new_image


def PillowRun(im_paths, operations, resizes, positions, masks, methods, saves, is_temps, texts, font_sizes):
    '''
    Image operation with pillow
    Inputs:
        im_paths (list): list of path of image files
        operations (list of str): operations to perform
        resizes (list of float): scale to resize to
        positions (list of tuple): position on the new figure
        masks: transparency mask
        methods: which method to use
            on_first_figure: take the first figure and lay everything on that
            use_new_one: take a new blank figure and lay everything on that
        saves (list of str): paths to save figure
        is_temps(list of int): whether to remove intermediate results
        texts (list of str): texts to append
        font_sizes (list of int): font sizes
    '''
    i = 0
    for im_path in im_paths:
        if operations[i] == 'new':
            last = im_path
            new_image = Image.open(im_path)
        elif operations[i] == 'paste':
            new_image = ImageMerge0([last, im_path], methods[i], (0, masks[i]), (1.0, resizes[i]), ((0.0, 0.0), positions[i]))
        elif operations[i] == 'crop':
            new_image = Image.open(im_path)
            new_image = new_image.crop(positions[i])
        elif operations[i] == 'text':
            fnt0 = ImageFont.truetype("Pillow/Tests/fonts/FreeMono.ttf", font_sizes[i])  # get a font
            new_image = Image.open(im_path)
            d = ImageDraw.Draw(new_image)
            d.text(positions[i], texts[i], font=fnt0, fill=(0, 0, 0))  # anchor option doesn't workf
        else:
            raise ValueError('operation must be in [new, paste, crop]')
        if saves[i] != '':
            if is_temps[i] == 0:
                print("%s: save figure %s" % (func_name(), saves[i]))
            if not os.path.isdir(os.path.dirname(saves[i])):
                os.mkdir(os.path.dirname(saves[i]))
            new_image.save(saves[i]) # save image
            last = saves[i]
        i += 1
    i = 0
    for save_path in saves:
        if is_temps[i] == 1:
            os.remove(save_path)  # remove intermediate results
        i += 1
