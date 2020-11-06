import re
import os
import shutil
import json
import numpy as np
from shilofue.Utilities import my_assert, re_neat_word, WriteFileHeader

'''
For now, my strategy is first defining a method to parse inputs for every key word,
then defining a class to parse different type of inputs.
Future:
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
        particle_data(list):
            list of particle coordinates
    '''
    # future: add interface for extra
    def __init__(self, _idict, **kwargs):
        '''
        initiate from a dictionary
        Inputs:
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
        # configurations
        my_assert(type(self.config)==dict, TypeError, 'Config must be a dictionary')
        self.extra = kwargs.get('extra', {})
        my_assert(type(self.extra)==dict, TypeError, 'extra must be a dictionary')
        self.test = kwargs.get('test', {})
        my_assert(type(self.test)==dict, TypeError, 'Test must be a dictionary')
        # list of particle coordinates
        self.particle_data = None
        
    def UpdatePrmDict(self, _names, _values):
        '''
        Update the dictionary of prm file,
        call function ChangeDiscValues()
        Inputs:
            _names(list): list of names
            _values(list): list of values
        '''
        ChangeDiscValues(self.idict, _names, _values)  # change values in idict accordingly

    def Intepret(self, parse_operations, **kwargs):
        '''
        Intepret configuration by calling function from the parse_operations class
        kwargs:
            extra(dict) - extra configuration
        '''
        _extra = kwargs.get('extra', {})
        _operations = kwargs.get('operations', [])
        _config = { **self.config, **self.test, **_extra }  # append _extra to self.config
        for operation in parse_operations.ALL_OPERATIONS:
            getattr(parse_operations, operation)(self.idict, _config)
        pass
    
    def process_particle_data(self):
        '''
        process the coordinates of particle, doing nothing here.
        Future class that inherit this one need to reload this method.
        '''
        pass

    def CaseName(self):
        '''
        Generate case name from self.config
        '''
        _case_name = ''
        for key, value in sorted(self.config.items(), key=lambda item: item[0]):
            if re.match("^_", key):
                # skip keys start with '_'
                continue
            _pattern = PatternFromStr(key)
            _pattern_value = PatternFromValue(value)
            _case_name += (_pattern + _pattern_value)
        if self.test != {}:
            _case_name += 'test'
            for key, value in sorted(self.test.items(), key=lambda item: item[0]):
                if re.match("^_", key):
                    # skip keys start with '_'
                    continue
                _pattern = PatternFromStr(key)
                _pattern_value = PatternFromValue(value)
                _case_name += (_pattern + _pattern_value)
        return _case_name
    
    def output_particle_ascii(self, fout):
        '''
        Output to a ascii file that contains Particle coordinates, containing the coordinates of each particle
        '''
        # header information
        _header = '# Ascii file for particle coordinates\n'
        # output particle file
        fout.write(_header)
        for i in range(self.particle_data.shape[0]):
            _string = ''
            for j in range(self.particle_data.shape[1]):
                if j == self.particle_data.shape[1] - 1:
                    _string += '%.4e\n' % self.particle_data[i, j]
                else:
                    _string += '%.4e ' % self.particle_data[i, j]
            fout.write(_string)
        pass

    def __call__(self, parse_operations, **kwargs):
        '''
        Create a .prm file
        inputs:
            parse_operations(class): operations to do
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
            _extra_files = kwargs.get('extra_file', {})  # extra dictionary, pass to intepret
            # Process particle data
            self.process_particle_data()
            # First intepret the configurations and update prm
            my_assert(self.config != None, ValueError,
                      'With the \'auto\' method, the config must exist')
            self.Intepret(parse_operations, extra=_extra)
            
            # Next generate a case name
            self.case_name = _basename + self.CaseName()
            
            # After that, make a directory with case name
            _case_dir = os.path.join(_dirname, self.case_name)
            # By default, we don't update
            update_ = kwargs.get('update', 0)
            if not update_:
                my_assert(os.path.isdir(_case_dir) is False, ValueError, 'Going to update a pr-exiting case, but update is not included in the option')
            if not os.path.isdir(_case_dir):
                os.mkdir(_case_dir)            

            # write configs to _json
            _json_outputs = {'basename': _basename, 'config': self.config, 'test': self.test, 'extra': _extra, 'extra_file': _extra_files}
            _json_ofile = os.path.join(_case_dir, 'config.json')
            with open(_json_ofile, 'w') as fout:
                json.dump(_json_outputs, fout)

            # At last, export a .prm file
            _filename = os.path.join(_case_dir, 'case.prm')
            with open(_filename, 'w') as fout:
                ParseToDealiiInput(fout, self.idict)

            # output particle data to an ascii file
            if self.particle_data is not None:
                _filename = os.path.join(_case_dir, 'particle.dat')
                with open(_filename, 'w') as fout:
                    self.output_particle_ascii(fout)

            # also copy the extra files
            if type(_extra_files) is str:
                _extra_files = [_extra_files]
            for _extra_file in _extra_files:
                shutil.copy2(_extra_file, _case_dir)

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
        return case names
        Return:
            _case_names(list)
        '''
        pass
    
    def __call__(self, parse_operations, _odir, **kwargs):
        '''
        Inputs:
            parse_operations(class): operations to do
            _odir(str):
                name of the target directory
            kwargs:
                extra(dict): extra configurations
                operation(dict): operations to do
                basename(str): base name for cases
        '''
        _extra = kwargs.get('extra', {})
        _base_name = kwargs.get('basename', '')
        # write configs to _json
        _json_outputs = self.configs
        _json_outputs['extra'] = _extra 
        _json_ofile = os.path.join(_odir, 'config.json')
        with open(_json_ofile, 'w') as fout:
            json.dump(_json_outputs, fout)

        # create cases in this group
        update_ = kwargs.get('update', 0)
        for _case in self.cases:
            _case_name = _case(parse_operations, dirname=_odir, extra=_extra, basename=_base_name, update=update_)
            self.case_names.append(_case_name)
        return self.case_names


class PARSE_OPERATIONS():
    """
    put global parse operations in a single class
    Attributes:
        ALL_OPERATIONS(list):
            all avalable operations
    """
    def __init__(self):
        """
        Initiation
        """
        self.ALL_OPERATIONS = ['MeshRefinement', 'Termination', 'Solver', 'MaterialModel']
        pass

    def MeshRefinement(self, Inputs, _config):
        '''
        change mesh refinement
        '''
        try:
            # Initial global refinement
            _initial_global_refinement = int(_config['initial_global_refinement'])
        except KeyError:
            pass
        else:
            Inputs['Mesh refinement']['Initial global refinement'] = str(_initial_global_refinement)
    
        try:
            # initial_adaptive_refinement
            _initial_adaptive_refinement = int(_config['initial_adaptive_refinement'])
        except KeyError:
            pass
        else:
            Inputs['Mesh refinement']['Initial adaptive refinement'] = str(_initial_adaptive_refinement)
        
        try:
            # refinement_fraction
            _refinement_fraction = float(_config['refinement_fraction'])
        except KeyError:
            pass
        else:
            Inputs['Mesh refinement']['Refinement fraction'] = str(_refinement_fraction)
        _complement_refine_coarse = _config.get('complement_refine_coarse', 0)
        if _complement_refine_coarse == 1:
            # make the sum of Refinement fraction and Coarsening fraction to 1.0
            try:
                Inputs['Mesh refinement']['Coarsening fraction'] = str(1.0 - float(Inputs['Mesh refinement']['Refinement fraction']))
            except KeyError:
                # there is no such settings in the inputs
                pass
        else:
            # import coarsening_fraction from file
            try:
                _coarsening_fraction = float(_config['coarsening_fraction'])
            except KeyError:
                pass
            else:
                Inputs['Mesh refinement']['Coarsening fraction'] = str(_coarsening_fraction)
    
        try:
            # only use minimum_refinement_function
            _only_refinement_function = int(_config['only_refinement_function'])
        except KeyError:
            pass
        else:
            if _only_refinement_function == 1:
                # only use minimum refinement function
                Inputs['Mesh refinement']['Strategy'] = 'minimum refinement function'
            if _only_refinement_function == 2:
                # include compositional gradient
                Inputs['Mesh refinement']['Strategy'] = 'minimum refinement function, composition approximate gradient'
            if _only_refinement_function == 3:
                # include viscosity
                Inputs['Mesh refinement']['Strategy'] = 'minimum refinement function, composition approximate gradient, viscosity'
            if _only_refinement_function == 4:
                # include strain rate
                Inputs['Mesh refinement']['Strategy'] = 'minimum refinement function, composition approximate gradient, strain rate'
        try:
            # longitude repetitions for chunk geometry
            _longitude_repetitions = int(_config['longitude_repetitions'])
        except KeyError:
            pass
        else:
            Inputs['Geometry model']['Chunk']['Longitude repetitions'] = str(_longitude_repetitions)

    def Termination(self, Inputs, _config):
        """
        Different termination strategy
        """
        try:
            # Initial global refinement
            _only_one_step = int(_config['only_one_step'])
        except KeyError:
            pass
        else:
            if _only_one_step == 1:
                Inputs['Termination criteria']['Termination criteria'] = 'end step'
                Inputs['Termination criteria']['End step'] = '1'
    
    def Solver(self, Inputs, _config):
        """
        solver parameters
        """
        # change the CFL number
        try:
            CFL_number = _config['CFL']
        except KeyError:
            pass
        else:
            Inputs['CFL number'] = str(CFL_number)
        
        # change the Max nonlinear iterations number
        try:
            max_nonlinear_iterations = _config['max_nonlinear_iterations']
        except KeyError:
            pass
        else:
            Inputs['Max nonlinear iterations'] = str(max_nonlinear_iterations)
        
        # change the nonlinear solver tolerance
        try:
            nonlinear_solver_tolerance = _config['nonlinear_solver_tolerance']
        except KeyError:
            pass
        else:
            Inputs['Nonlinear solver tolerance'] = str(nonlinear_solver_tolerance)

        # change Newton solver configuration
        try:
            newton_solver = Inputs['Solver parameters']['Newton solver parameters']
        except KeyError:
            pass
        else:
            try:
                # change "Nonlinear Newton solver switch tolerance"
                newton_solver_switch = _config['newton_solver_switch']
            except KeyError:
                pass
            else:
                newton_solver['Nonlinear Newton solver switch tolerance'] = str(newton_solver_switch)
            try:
                # change "Max pre-Newton nonlinear iterations"
                max_pre_newton_iteration = _config['max_pre_newton_iteration']
            except KeyError:
                pass
            else:
                newton_solver['Max pre-Newton nonlinear iterations'] = str(max_pre_newton_iteration)
            try:
                # change "Maximum linear Stokes solver tolerance"
                max_linear_tolerance = _config['max_linear_tolerance']
            except KeyError:
                pass
            else:
                newton_solver['Maximum linear Stokes solver tolerance'] = str(max_linear_tolerance)
    
    def MaterialModel(self, Inputs, _config):
        """
        properties in the material model
        only support changing parameters inside a material model
        """

        # Get model configurations from a prm file
        try: 
            model_name = Inputs['Material model']['Model name']
        except KeyError:
            return

        if model_name == 'visco plastic': 
            model = Inputs['Material model']['Visco Plastic']
            # change lower limit
            try:
                minimum_viscosity = _config['minimum_viscosity']
            except KeyError:
                pass
            else:
                model['Minimum viscosity'] = str(minimum_viscosity)
            # change upper limit
            try:
                maximum_viscosity = _config['maximum_viscosity']
            except KeyError:
                pass
            else:
                model['Maximum viscosity'] = str(maximum_viscosity)
            # parse back
            Inputs['Material model']['Visco Plastic'] = model


class BASH_OPTIONS():
    """
    parse .prm file to a option file that bash can easily read
    Attributes:
        _case_dir(str): path of this case
        _output_dir(str): path of the output
        _visit_file(str): path of the visit file
        odict(dict): dictionary of key and value to output
    """
    def __init__(self, case_dir):
        """
        Initiation
        Args:
            case_dir(str): directory of case
        """
        # check directory
        self._case_dir = case_dir
        my_assert(os.path.isdir(self._case_dir), FileNotFoundError,
                  'BASH_OPTIONS.__init__: case directory - %s doesn\'t exist' % self._case_dir)
        self._output_dir = os.path.join(case_dir, 'output')
        my_assert(os.path.isdir(self._output_dir), FileNotFoundError,
                  'BASH_OPTIONS.__init__: case output directory - %s doesn\'t exist' % self._output_dir)
        self._visit_file = os.path.join(self._output_dir, 'solution.visit')
        my_assert(os.access(self._visit_file, os.R_OK), FileNotFoundError,
                  'BASH_OPTIONS.__init__: case visit file - %s cannot be read' % self._visit_file)
        # output dir
        self._output_dir = os.path.join(case_dir, 'output')
        if not os.path.isdir(self._output_dir):
            os.mkdir(self._output_dir)
        # img dir
        self._img_dir = os.path.join(case_dir, 'img')
        if not os.path.isdir(self._img_dir):
            os.mkdir(self._img_dir)
        
        # get inputs from .prm file
        prm_file = os.path.join(self._case_dir, 'case.prm')
        my_assert(os.access(prm_file, os.R_OK), FileNotFoundError,
                  'BASH_OPTIONS.__init__: case prm file - %s cannot be read' % prm_file)
        with open(prm_file, 'r') as fin:
            self.idict = ParseFromDealiiInput(fin)

        # initiate a dictionary
        self.odict = {}
    
    def Interpret(self, kwargs):
        """
        Interpret the inputs, to be reloaded in children
        """
        pass

    def __call__(self, ofile, kwargs):
        """
        Call function
        Args:
            ofile(str): path of output
        """
        # interpret
        self.Interpret(kwargs)

        # open ofile for output
        # write outputs by keys and values
        with open(ofile, 'w') as fout:
            for key, value in self.odict.items():
                fout.write("%s       %s\n" % (key, value))
        pass


class VISIT_OPTIONS(BASH_OPTIONS):
    """
    todo
    parse .prm file to a option file that bash can easily read
    """
    def Interpret(self, kwargs={}):
        """
        Interpret the inputs, to be reloaded in children
        """
        # call function from parent
        BASH_OPTIONS.Interpret(self, kwargs)
        
        # visit file
        self.odict["VISIT_FILE"] = self._visit_file

        # particle file
        particle_file = os.path.join(self._output_dir, 'particles.visit')
        if os.access(particle_file, os.R_OK):
            self.odict["VISIT_PARTICLE_FILE"] = particle_file

        # directory to output data 
        self.odict["DATA_OUTPUT_DIR"] = self._output_dir

        # directory to output images
        if not os.path.isdir(self._img_dir):
            os.mkdir(self._img_dir)
        self.odict["IMG_OUTPUT_DIR"] = self._img_dir

        # own implementations
        # initial adaptive refinement
        self.odict['INITIAL_ADAPTIVE_REFINEMENT'] = self.idict['Mesh refinement'].get('Initial adaptive refinement', '6')
        # SNAPSHOT index
        self.odict['SINGLE_SNAPSHOT'] = kwargs.get('SINGLE_SNAPSHOT', 0)
        self.odict['MULTIPLE_SNAPSHOTS'] = kwargs.get('MULTIPLE_SNAPSHOTS', [])


class VISIT_XYZ():
    """
    Read .xyz file exported from visit and do analysis
    Attributes:
        data(nparray): data
        header(dict): headers in file
        output_data(nparray): data to output
        output_header(dict): output headers in file
        self.column_indexes(dict): record index of col in self.data
    """
    def __init__(self):
        """
        initiation
        """
        self.data = []
        self.output_data = np.array([])
        self.header = {}
        self.output_header = {}
        self.column_indexes = {}
        pass

    def ReadData(self, filein):
        """
        read date form file
        Args:
            filein(str): file input
        """
        # assert file
        my_assert(os.access(filein, os.R_OK), FileNotFoundError,
                  'VISIT_XYZ.__init__: visit xyz file - %s cannot be read' % filein)

        # construct cols
        cols = []
        i = 0
        for key, value in self.header.items():
            cols.append(value['col'])
            # record column in self.data
            self.column_indexes[key] = i
            i += 1

        # read data
        self.data = np.loadtxt(filein, usecols=(cols), skiprows=2)

    def Analyze(self, kwargs):
        """
        analyze data
        """
        pass

    def Output(self, ofile):
        """
        output data
        Inputs:
            kwargs(dict): options
                ofile(str): file to output
        """
        # create a new file and write header
        if not os.path.isfile(ofile):
            WriteFileHeader(ofile, self.output_header)

        with open(ofile, 'a') as fout:
            np.savetxt(fout, self.output_data, fmt='%-20.8e')

    def __call__(self, filein, **kwargs):
        """
        call function
        Args:
            filein(str): file input
            kwarg(dict): options
                cols(tuple): columns to import
        """
        # options
        default_header = {
            'x': {'col': 1 },
            'y': {'col': 2 },
            'id': {'col': 3 }
        }
        self.header = kwargs.get("header", default_header)

        # read data 
        self.ReadData(filein)

        # analyze
        self.Analyze(kwargs)

        # output if file given 
        ofile = kwargs.get('ofile', None)
        if ofile is not None:
            self.Output(ofile)


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
        elif re.match('^(\t| )*set', line):
            # Parse key and value
            # from format in file as 'set key = val'
            # to a dictionary inputs
            # inputs[key] = val
            temp = re.sub('^(\t| )*set ', '', line, count=1)
            temp = temp.split('=', maxsplit=1)
            key = temp[0]
            key = re.sub('(\t| )*$', '', key)
            # key = key.strip(' ')
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
        # if re.match("sub_group", key):
        #   continue
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

    # header of the file
    _contents = ''
    if not os.path.isfile(_md_file):
        _contents = '# Group %s\n\n' % _group_name  # an string to hold content of file
        _contents += '## Overview\n\n'

    # configures of group    
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
    with open(_md_file, 'a') as fout:
        # output file
        fout.write(_contents)
    pass
        
        
def UpdateProjectMd(_project_dict, _project_dir):
    '''
    Update auto mkd file for all cases in this project
    '''
    # deal with cases and groups
    for key, value in _project_dict.items():
        if key == 'cases':
            # cases, handled afterwards
            _dir = _project_dir
        else:
            # groups, first handle the group level
            _dir = os.path.join(_project_dir, key)
            _group_json = os.path.join(_dir, 'config.json')
            with open(_group_json, 'r') as fin:
                _group_config_dict = json.load(fin)
            AutoMarkdownGroup(key, _group_config_dict, dirname=_dir)
        for _case in _project_dict[key]:
            # cases, as there are the value related to the key 'cases'
            # or related to the name of a group
            _case_dir = os.path.join(_dir, _case)
            _case_json = os.path.join(_case_dir, 'config.json')
            with open(_case_json, 'r') as fin:
                _case_config_dict = json.load(fin)
            AutoMarkdownCase(_case, _case_config_dict, dirname=_case_dir)
    

def UpdateProjectJson(_dir, **kwargs):
    '''
    update groups and files information of a project
    export to a json file
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