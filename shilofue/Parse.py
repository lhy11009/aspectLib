import re
from shilofue.Utilities import my_assert

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
            values_str = part.split(':')[1].split('|')
            # convert string to float
            values = [float(val) for val in values_str]
            self.data[part.split(':')[0]] = values

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


class CASE():
    '''
    class for a case
    Attributes:
    '''
    def __init__(self, _config):
        '''
        initiate from a dictionary
        Inputs:
            _config(dict):
                a dictionary that has two members: names and values
                'names'(list):
                    a sequence of keys
                'values'(list):
                    a sequence of values
        '''
        self.names = _config('names')
        self.values = _config('values')

    def Parse(self, _ofile):
        '''
        parse case to a .prm file
        '''
        pass


class GROUP_CASE():
    '''
    Class for a group of cases
    Attributes:
        cases(list<class CASE>):
            a list of cases
    '''
    def __init__(self, _idict):
        '''
        initiate from a dictionary
        '''
        self.cases = []  # initiate a list to save parsed cases
        _names, _parameters = GetGroupCaseFromDict(_idict)  # Get a list for names and a list for parameters from a dictionary read from a json file
        _cases_config = ExpandNamesParameters(_names, _parameters)
        for _case_config in _cases_config:
            # initiate a new case with __init__ of clase CASE and append
            self.cases.append(CASE(_case_config))
    
    def Parse(self, _target_dir):
        '''
        Inputs:
            _target_dir(str):
                name of the target directory
        '''
        for _case in self.cases:
            _ofile = 'foo.prm'
            _case.Parse(_ofile)


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
