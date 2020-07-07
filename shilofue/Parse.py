import re
import os
import json
from shilofue.Utilities import my_assert, re_neat_word

'''
For now, my strategy is first defining a method to parse inputs for every key word,
then defining a class to parse different type of inputs.
Todo:
    a. put them in a bigger class
    b. add method for functions
'''


class COMPOSITION():
    """
    store value like
      'background:4.0e-6|1.5e-6, spcrust:0.0, spharz:4.0e-6, opcrust:4.0e-6, opharz:4.0e-6 '
    or parse value back
    """
    def __init__(self, line):
        # parse the format:
        # key1: val1|val2, key2: val3|val4
        # to a dictionary data where
        # data[key1] = [val1, val2]
        self.data = {}
        parts = line.split(',')
        for part in parts:
            key_str = part.split(':')[0]
            key = re_neat_word(key_str)
            values_str = part.split(':')[1].split('|')
            # convert string to float
            values = [float(re_neat_word(val)) for val in values_str]
            self.data[key] = values

    def parse_back(self):
        """
        def parse_back(self)

        parse data back to a string
        """
        line = ''
        j = 0
        for key, values in self.data.items():
            # construct the format:
            # key1: val1|val2, key2: val3|val4
            if j > 0:
                part_of_line = ', ' + key + ':'
            else:
                part_of_line = key + ':'
            i = 0
            for val in values:
                if i == 0:
                    part_of_line += '%.4e' % val
                else:
                    part_of_line += '|' + '%.4e' % val
                i += 1
            line += part_of_line
            j += 1
        return line


class CASE():
    '''
    class for a case
    Attributes:
        names(list):
            list for name of variables to change
        values(list):
            list of value of variables to change
        idict(dict):
            dictionary of parameters
        config(dict):
            dictionary of configuration for the parameters
    '''
    def __init__(self, _idict, **kwargs):
        '''
        initiate from a dictionary
        Inputs:
            case_name(str): case name of the model
            _idict(dict):
                dictionary import from a base file
            kwargs:
                config: (dict) - a dictionary that contains the configuration
                test: (dict) - a dictionary that contains the configuration to test
        '''
        my_assert(type(_idict)==dict, TypeError, "First entry mush be a dictionary")
        self.case_name = ''
        self.idict = _idict
        self.config = kwargs.get('config', {})
        my_assert(type(self.config)==dict, TypeError, 'Config must be a dictionary')
        self.test = kwargs.get('test', {})
        my_assert(type(self.test)==dict, TypeError, 'Test must be a dictionary')
        
    def UpdatePrmDict(self, _names, _values):
        '''
        Update the dictionary of prm file,
        call function ChangeDiscValues()
        Inputs:
            _names(list): list of names
            _values(list): list of values
        '''
        ChangeDiscValues(self.idict, _names, _values)  # change values in idict accordingly

    def Intepret(self, **kwargs):
        '''
        Intepret configuration,
        to be defined in children class
        '''
        pass


    def CaseName(self):
        '''
        Generate case name from self.config
        '''
        _case_name = ''
        for key, value in sorted(self.config.items(), key=lambda item: item[0]):
            _pattern = PatternFromStr(key)
            _pattern_value = PatternFromValue(value)
            _case_name += (_pattern + _pattern_value)
        if self.test != {}:
            _case_name += 'test'
            for key, value in sorted(self.test.items(), key=lambda item: item[0]):
                _pattern = PatternFromStr(key)
                _pattern_value = PatternFromValue(value)
                _case_name += (_pattern + _pattern_value)
        return _case_name

    def __call__(self, **kwargs):
        '''
        Create a .prm file
        inputs:
            kwargs:
                method: (str) - method of generate files
                dirname: (str) - output directory, in use with 'auto' method
                basename: (str) - base for case name, in use with 'auto' method
                filename: (str) - output file, in use with 'manual' method
        '''
        # assign file name with a method defined
        _method = kwargs.get('method', 'auto')
        if _method == 'auto':
            _dirname = kwargs.get('dirname', '.')
            _basename = kwargs.get('basename', '')
            _extra = kwargs.get('extra', {})  # extra dictionary, pass to intepret
            _operations = kwargs.get('operations', [])  # operations to do, pass to intepret
            # First intepret the configurations and update prm
            my_assert(self.config != None, ValueError,
                      'With the \'auto\' method, the config must exist')
            self.Intepret(extra=_extra, operations=_operations)
            # Next generate a case name
            self.case_name = _basename + self.CaseName()
            # After that, make a directory with case name
            _case_dir = os.path.join(_dirname, self.case_name)
            my_assert(os.path.isdir(_case_dir) is False, ValueError, 'The script doesn\'t support updating a pr-exiting group')
            os.mkdir(_case_dir)
            # write configs to _json
            _json_outputs = {'config': self.config, 'test': self.test, 'extra': _extra} # todo
            _json_ofile = os.path.join(_case_dir, 'config.json')
            with open(_json_ofile, 'w') as fout:
                json.dump(_json_outputs, fout)
            # At last, export a .prm file
            _filename = os.path.join(_case_dir, 'case.prm')
            with open(_filename, 'w') as fout:
                ParseToDealiiInput(fout, self.idict)
            pass
        elif _method == 'manual':
            # export a .prm file
            _filename = kwargs.get('filename', None)
            with open(_filename, 'w') as fout:
                ParseToDealiiInput(fout, self.idict)
            pass
        return self.case_name


def PatternFromStr(_str):
    '''
    Generate a pattern from a _str
    Inputs:
        _str(str): input string
    Returns:
        _pattern(str): pattern
    '''
    _pattern=''
    _parts = _str.split('_')
    for _part in _parts:
        # connect upper case of every first character
        _pattern += _part[0].upper()
    return _pattern


def PatternFromValue(_value):
    '''
    Generate a pattern from a _value
    Inputs:
        _str(str): input string
    Returns:
        _pattern(str or float or int): pattern
    '''
    if type(_value) is str:
        _pattern = ''
        # connect the lower case of every first character
        _parts = _value.split(' ')
        for _part in _parts:
            _pattern += _part[0].lower()
    elif type(_value) is float:
        _pattern='%.3e' % _value
    elif type(_value) is int:
        _pattern='%d' % _value
    else:
        raise(TypeError('value must be str or float or int'))
    return _pattern


class GROUP_CASE():
    '''
    Class for a group of cases
    Attributes:
        cases(list<class CASE>):
            a list of cases
    '''
    def __init__(self, CASE_CLASS, _inputs, _json_inputs):
        '''
        initiate from a dictionary
        Attribute:
            case_names(list): list of casa names
        Inputs:
            _inputs(dict): inputs from a base input file
            _json_inputs(dict): inputs from a json file containing configurations
        '''
        self.cases = []  # initiate a list to save parsed cases
        self.case_names = []
        self.configs = _json_inputs
        self.config_tests = GetGroupCaseFromDict1(_json_inputs)  # Get a list for names and a list for parameters from a dictionary read from a json file
        for _config_test in self.config_tests:
            _config = _config_test['config']
            _test = _config_test.get('test', {})
            self.cases.append(CASE_CLASS(_inputs, config=_config, test=_test))

    def get_case_names(self):
        '''
        return case names, todo
        Return:
            _case_names(list)
        '''
        pass
    
    def __call__(self, _odir, **kwargs):
        '''
        Inputs:
            _odir(str):
                name of the target directory
            kwargs:
                extra(dict): extra configurations
                operation(dict): operations to do
                basename(str): base name for cases
        '''
        _extra = kwargs.get('extra', {})
        _operations = kwargs.get('operations', [])
        _base_name = kwargs.get('basename', '')
        # write configs to _json
        _json_outputs = self.configs
        _json_outputs['extra'] = _extra 
        _json_ofile = os.path.join(_odir, 'config.json')
        with open(_json_ofile, 'w') as fout:
            json.dump(_json_outputs, fout)
        for _case in self.cases:
            _case_name = _case(dirname=_odir, extra=_extra, operations=_operations, basename=_base_name)
            self.case_names.append(_case_name)
        return self.case_names


def ParseFromDealiiInput(fin):
    """
    ParseFromDealiiInput(fin)

    Parse Dealii input file to a python dictionary
    """
    inputs = {}
    line = fin.readline()
    while line is not "":
        # Inputs formats are
        # comment: "# some comment"
        # start and end of new section:
        # "subsection name" and "end"
        # set variables values:
        # 'set key = val'
        if re.match('^(\t| )*#', line):
            # Skip comment lines, mark by '#' in file
            pass
        elif re.match('^.*set', line):
            # Parse key and value
            # from format in file as 'set key = val'
            # to a dictionary inputs
            # inputs[key] = val
            temp = re.sub('^(\t| )*set ', '', line, count=1)
            temp = temp.split('=', maxsplit=1)
            key = temp[0]
            key = key.strip(' ')
            value = temp[1]
            value = re.sub('^ *', '', value)
            value = re.sub(' *(#.*)?\n$', '', value)
            while value[-1] == '\\':
                # Deal with entries that extent to
                # multiple lines
                line = fin.readline()
                line = re.sub(' *(#.*)?\n$', '', line)
                value = value + '\n' + line
            inputs[key] = value
        elif re.match('^.*subsection', line):
            # Start a new subsection
            # Initialize new dictionary and interatively call function,
            key = re.sub('^.*subsection ', '', line)
            key = key.strip('\n')
            try:
                # Fix the bug where a subsection emerges
                # multiple times
                inputs[key]
            except KeyError:
                inputs[key] = ParseFromDealiiInput(fin)
            else:
                temp = ParseFromDealiiInput(fin)
                inputs[key].update(temp.items())
        elif re.match('^.*end', line):
            # Terminate and return, marked by 'end' in file
            return inputs
        line = fin.readline()
    return inputs


def ParseToDealiiInput(fout, outputs, layer=0):
    """
    def ParseToDealiiInput(fout, outputs, layer=0)

    Parse a python dictionary into a Dealii input file
    """
    indent = ' ' * 4 * layer  # Indentation of output
    for key, value in outputs.items():
        # output format is
        # same as input but no comment
        if type(value) is str:
            fout.write(indent + 'set %s = %s\n' % (key, value))
        elif type(value) is dict:
            if layer == 0:
                fout.write('\n')
            fout.write(indent + 'subsection %s\n' % key)
            layer1 = layer + 1
            ParseToDealiiInput(fout, value, layer1)
            fout.write(indent + 'end\n')
            if layer == 0:
                fout.write('\n')
        else:
            raise ValueError('Value in dict must be str')
    return


def GetGroupCaseFromDict(_idict):
    '''
    Get a list for names and a list for parameters from a dictionary read from a json file
    Inputs:
    _idict(dict):
        input dictionary
    Returns:
        _names(list):
            list of names, each member is a list itself
        _parameters(list):
            list of parameters, each member is a list itself
    '''
    my_assert(type(_idict) == dict, TypeError, 'Input is not a dictionary')
    _parameters = []  # initialize a array to load parameters
    _names = []  # initialize a array to load names of parameters
    for key, value in sorted(_idict.items(), key=lambda item: item[0]):
        if type(value) is dict:
            # in a top hierachy, append the name and call recursively
            _sub_names, _sub_parameters = GetGroupCaseFromDict(value)
            for i in range(len(_sub_names)):
                # concatenate names and append
                _names.append([key] + _sub_names[i])
            _parameters += _sub_parameters  # append parameters
        elif type(value) is list:
            _names.append([key])  # concatenate names
            _parameters.append(value)
        elif type(value) in [int, float, str]:
            _names.append([key])  # concatenate names
            _parameters.append([value])
        else:
            raise TypeError('%s is not int, float or str' % str(type(value)))
    return _names, _parameters


def GetGroupCaseFromDict1(_idict):
    '''
    second method(simpler one) of getting a list for names and a list for parameters from a dictionary read from a json file
    Inputs:
        _idict(dict):
            input dictionary
    returns:
        _config_tests(dict of dict):
            list of configuretions and tests for cases
    '''
    my_assert(type(_idict) == dict, TypeError, 'Input is not a dictionary')
    _config_tests=[]
    _configs = _idict['config']
    _tests = _idict.get('test', {})
    # get total number of cases
    _total = 1
    _totals = [1]
    for key, value in sorted(_configs.items(), key=lambda item: item[0]):
        if type(value) == list:
            _total *= len(value)
        _totals.append(_total)
    for key, value in sorted(_tests.items(), key=lambda item: item[0]):
        if type(value) == list:
            _total *= len(value)
        _totals.append(_total)
    # get case configuration and test
    for j in range(_total):
        # loop for every case
        _config_test = {'config': {}, 'test': {}}  # initiate dict
        i = 0  # current index of parameters
        for key, value in sorted(_configs.items(), key=lambda item: item[0]):
            # loop for configurations
            # derive index number by mod
            _ind = j
            _ind = int(_ind // _totals[i])
            try:
                _ind = _ind % len(value)
                # indexing by _ind and append value to key
                _config_test['config'][key] = value[_ind]
            except TypeError:
                # only one value
                _config_test['config'][key] = value
            i += 1
        for key, value in sorted(_tests.items(), key=lambda item: item[0]):
            # loop for tests
            _ind = j
            _ind = int(_ind // _totals[i])
            try:
                _ind = _ind % len(value)
                # indexing by _ind and append value to key
                _config_test['test'][key] = value[_ind]
            except TypeError:
                # only one value
                _config_test['test'][key] = value
            i += 1
        _config_tests.append(_config_test)
    return _config_tests


def ExpandNamesParameters(_names, _parameters):
    '''
    Inputs:
        _names(list):
            list of names, each member is a list itself
        _parameters(list):
            list of parameters, each member is a list itself
    Returns:
        _cases_config(list<dict>):
            a list of dictionaries. One dictionary is a config for a list file
    '''
    my_assert(type(_names) == list, TypeError, 'First Entry is not a list')
    my_assert(type(_parameters) == list, TypeError, 'Second Entry is not a list')
    my_assert(len(_names) == len(_parameters), ValueError, 'Length of first and second entry is not equal')
    _total = 1
    for _sub_parameters in _parameters:
        # take the value of total of all lengths multiplied
        _total *= len(_sub_parameters)
    # initialize this list of dictionaries for configurations
    _cases_config = []
    for i in range(_total):
        _cases_config.append({'names': [], 'values': []})
    # fill in all entries
    for j in range(len(_cases_config)):
        _cases_config[j]['names'] = _names.copy()
        for i in range(len(_names)):
            _ind = j  # get the index in _parameters[i]
            for k in range(len(_names)-1, i, -1):
                _ind = int(_ind // len(_parameters[k]))
            _ind = _ind % len(_parameters[i])
            _cases_config[j]['values'].append(_parameters[i][_ind])
    return _cases_config


def ChangeDiscValues(_idict, _names, _values):
    '''
    Change values in a complex dictionary with names and values
    Inputs:
        _idict(dict):
            Dictionary of parameters from a .prm file
        _names(list):
            list of parameters, each member is a list it self,
            which contains the path the the variable
        _values(list):
            list of values of variables
    '''
    my_assert(type(_idict) == dict, TypeError, 'First Entry needs to be a dict')
    my_assert(type(_names) == list, TypeError, 'Second Entry needs to be a list')
    my_assert(type(_values) == list, TypeError, 'Third Entry needs to be a list')
    my_assert(len(_names) == len(_values), ValueError, 'Length of second and third entries must match')
    for i in range(len(_names)):
        _name = _names[i]
        _sub_dict = _idict
        for _key in _name[0: len(_name)-1]:
            _sub_dict = _sub_dict[_key]
        _sub_dict[_name[-1]] = _values[i]
    
    
def AutoMarkdownCase(_case_name, _idict, **kwargs):
    '''
    generate a automatic markdown file for a case
    Inputs:
        _idict(dict): dictionary that holds configurations for this case
        _case_name(str): case name
        kwargs:
            md: file name of markdown file
    '''
    _md = kwargs.get('md', 'auto.md')
    _dirname = kwargs.get('dirname', '.')
    _md_file = os.path.join(_dirname, _md)
    _contents = '# Case %s\n\n' % _case_name  # an string to hold content of file
    _contents += '## Overview\n\n'
    _config = _idict['config']  # append configuration part
    _contents += '### The case is configured with:\n\n'
    for key, value in sorted(_config.items(), key=lambda item: item[0]):
        _contents += '%s: %s\n\n' % (key, str(value))
    _test = _idict['test']  # append test part
    _contents += '### The case is tested with:\n\n'
    for key, value in sorted(_test.items(), key=lambda item: item[0]):
        _contents += '%s: %s\n\n' % (key, str(value))
    _extra = _idict['extra']  # append extra setting part
    _contents += '### The case is genearated with extra settings:\n\n'
    for key, value in sorted(_extra.items(), key=lambda item: item[0]):
        _contents += '%s: %s\n\n' % (key, str(value))
    with open(_md_file, 'w') as fout:
        # output file
        fout.write(_contents)
    pass


def AutoMarkdownGroup(_group_name, _idict, **kwargs):
    '''
    generate a automatic markdown file for a group
    Inputs:
        _idict(dict): dictionary that holds configurations for this case
        _case_name(str): case name
        kwargs:
            md: file name of markdown file
    '''
    _md = kwargs.get('md', 'auto.md')
    _dirname = kwargs.get('dirname', '.')
    _md_file = os.path.join(_dirname, _md)
    _contents = '# Group %s\n\n' % _group_name  # an string to hold content of file
    _contents += '## Overview\n\n'
    _config = _idict.get('config', {})  # append configuration part
    _contents += '### The group is configured with:\n\n'
    for key, value in sorted(_config.items(), key=lambda item: item[0]):
        _contents += '%s: %s\n\n' % (key, str(value))
    _test = _idict.get('test', {})  # append test part
    _contents += '### The group is tested with:\n\n'
    for key, value in sorted(_test.items(), key=lambda item: item[0]):
        _contents += '%s: %s\n\n' % (key, str(value))
    _extra = _idict.get('extra', {})  # append extra setting part
    _contents += '### The group is genearated with extra settings:\n\n'
    for key, value in sorted(_extra.items(), key=lambda item: item[0]):
        _contents += '%s: %s\n\n' % (key, str(value))
    with open(_md_file, 'w') as fout:
        # output file
        fout.write(_contents)
    pass
        
        
def UpdateProjectMd(_project_dict, _project_dir):
    '''
    Update auto mkd file for all cases in this project
    '''
    for key, value in _project_dict.items():
        if key == 'cases':
            _dir = _project_dir
        else:
            # groups
            _dir = os.path.join(_project_dir, key)
            _group_json = os.path.join(_dir, 'config.json')
            with open(_group_json, 'r') as fin:
                _group_config_dict = json.load(fin)
            AutoMarkdownGroup(key, _group_config_dict, dirname=_dir)
        for _case in _project_dict[key]:
            # cases
            _case_dir = os.path.join(_dir, _case)
            _case_json = os.path.join(_case_dir, 'config.json')
            with open(_case_json, 'r') as fin:
                _case_config_dict = json.load(fin)
            AutoMarkdownCase(_case, _case_config_dict, dirname=_case_dir)


def UpdateProjectJson(_dir, **kwargs):
    '''
    update groups and files information of a project
    export to a json file
    todo
    Inputs:
        _dir(str): directory of a project
        kwargs:
            json: json_file
    Write:
        the file by kwargs['json']
    Return:
        the dictionary generated
    '''
    _json = kwargs.get('json', 'project.json')
    _cases = []
    _groups = []
    # walk through the directory and look for groups and cases
    for _subname in os.listdir(_dir):
        # loop in _dir
        _fullsubname = os.path.join(_dir, _subname)
        if os.path.isdir(_fullsubname):
            if 'config.json' in os.listdir(_fullsubname):
                if 'case.prm' in os.listdir(_fullsubname):
                    _cases.append(_subname)
                else:
                    _groups.append(_subname)
    # construct a output dictionary
    _odict = {"cases": _cases}
    for _group in _groups:
        _group_dir = os.path.join(_dir, _group)
        _sub_cases = []  # an array to hold names of sub-cases
        for _subname in os.listdir(_group_dir):
            _fullsubname = os.path.join(_group_dir, _subname)
            if os.path.isdir(_fullsubname):
                if 'case.prm' in os.listdir(_fullsubname):
                    _sub_cases.append(_subname)
        _odict[_group] = _sub_cases
    # output to a json file
    _json_file = os.path.join(_dir, _json)
    with open(_json_file, 'w') as fout:
        json.dump(_odict, fout)
    return _odict