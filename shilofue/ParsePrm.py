# -*- coding: utf-8 -*-
r"""Parse prm file and apply operations
future:
replace functions in Parse.py

This exports:

  -

This depends on:

  -

Examples of usage:

  - default usage:

    generate .prm file with fast zero step:
        python -m shilofue.ParsePrm fast_zero_step -i case.prm -o case_test.prm

descriptions:
    copy and pasted all the function in the origin PARSE_OPERATION class
"""
from multiprocessing.sharedctypes import Value
import numpy as np
import sys, os, argparse
import re
# import pathlib
# import subprocess
import numpy as np
import json
import copy

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')

sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities

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
            key = Utilities.re_neat_word(key_str)
            values_str = part.split(':')[1].split('|')
            # convert string to float
            values = [float(Utilities.re_neat_word(val)) for val in values_str]
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


def ParseFromDealiiInput(fin):
    """
    ParseFromDealiiInput(fin)

    Parse Dealii input file to a python dictionary
    """
    inputs = {}
    line = fin.readline()
    while line != "":
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
                print('%s is already presented, going to update.' % key)
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
            raise ValueError('Value in dict must be str, get:\n key: '\
            + key + "\n type of value: " + str(type(value)) + "\n value: " + str(value))
    return


def ReadPrmFile(_path):
    """
    simply a wrapper for the parse function
    Inputs:
        _path: path for a prm file
    """
    Utilities.my_assert(os.access(_path, os.R_OK), FileNotFoundError, "%s: prm file %s doesn't exist" % (Utilities.func_name(), _path))
    with open(_path, 'r') as fin:
        inputs = ParseFromDealiiInput(fin)
    return inputs


def WritePrmFile(_path, outputs):
    """
    simply a wrapper for the output function
    Inputs:
        _path: path for a prm file
        outputs: a dictionary contains the settings
    """
    with open(_path, 'w') as fout:
        ParseToDealiiInput(fout, outputs)


def MeshRefinement(Inputs, _config):
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
        # Time steps between mesh refinement
        _steps_between_refinement = int(_config['steps_between_refinement'])
    except KeyError:
        pass
    else:
        Inputs['Mesh refinement']['Time steps between mesh refinement'] = str(_steps_between_refinement)

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


def Termination(Inputs, _config):
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


def Solver(Inputs, _config):
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

    # change Stokes solver configuration
    try:
        stokes_solver = Inputs['Solver parameters']['Stokes solver parameters']
    except KeyError:
        pass
    else:
        try:
            # change linear tolerance
            linear_stokes_tolerance = _config['linear_stokes_tolerance']
        except KeyError:
            pass
        else:
            stokes_solver['Linear solver tolerance'] = str(linear_stokes_tolerance)
        try:
            # change cheap stokes steps
            cheap_stokes_steps = _config['cheap_stokes_steps']
        except KeyError:
            pass
        else:
            stokes_solver['Number of cheap Stokes solver steps '] = str(cheap_stokes_steps)

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
        try:
            # change "line search iteration"
            max_line_search = int(_config['max_line_search'])
        except KeyError:
            pass
        else:
            newton_solver['Max Newton line search iterations'] = str(max_line_search)


def GetStokesSolverTypeFromPrm(_path):
    '''
    read the type of the stokes solver from a file
    '''
    inputs = ReadPrmFile(_path)
    stokes_solver_type = GetStokesSolverType(inputs)
    return stokes_solver_type


def GetStokesSolverType(inputs):
    '''
    read the type of the stokes solver from a dictionary
    '''
    try:
        stokes_solver_type = inputs["Solver parameters"]["Stokes solver parameters"]["Stokes solver type"]
    except KeyError:
        stokes_solver_type = "block AMG"
    return stokes_solver_type


def MaterialModel(Inputs, _config):
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


def UpperMantleRheologyViscoPlastic(Inputs):
    '''
    parse upper mante rheology
    '''
    visco_plastic = Inputs["Material model"]['Visco Plastic']
    grain_size = float(visco_plastic["Grain size"])
    prefactors_for_diffusion_creep = COMPOSITION(visco_plastic["Prefactors for diffusion creep"])
    grain_size_exponents_for_diffusion_creep  = COMPOSITION(visco_plastic["Grain size exponents for diffusion creep"])
    activation_energies_for_diffusion_creep = COMPOSITION(visco_plastic["Activation energies for diffusion creep"])
    activation_volumes_for_diffusion_creep  = COMPOSITION(visco_plastic["Activation volumes for diffusion creep"])
    prefactors_for_dislocation_creep = COMPOSITION(visco_plastic["Prefactors for dislocation creep"])
    activation_energies_for_dislocation_creep = COMPOSITION(visco_plastic["Activation energies for dislocation creep"])
    activation_volumes_for_dislocation_creep  = COMPOSITION(visco_plastic["Activation volumes for dislocation creep"])
    stress_exponents_for_dislocation_creep = COMPOSITION(visco_plastic["Stress exponents for dislocation creep"])
    # call GetLowerMantleRheology to derive parameters for lower mantle flow law
    backgroud_upper_mantle_diffusion = {}
    backgroud_upper_mantle_diffusion['A'] = prefactors_for_diffusion_creep.data['background'][0]
    backgroud_upper_mantle_diffusion['d'] = grain_size
    backgroud_upper_mantle_diffusion['n'] = 1.0
    backgroud_upper_mantle_diffusion['m'] = grain_size_exponents_for_diffusion_creep.data['background'][0]
    backgroud_upper_mantle_diffusion['E'] = activation_energies_for_diffusion_creep.data['background'][0]
    backgroud_upper_mantle_diffusion['V'] = activation_volumes_for_diffusion_creep.data['background'][0]
    backgroud_upper_mantle_dislocation = {}
    backgroud_upper_mantle_dislocation['A'] = prefactors_for_dislocation_creep.data['background'][0]
    backgroud_upper_mantle_dislocation['d'] = grain_size
    backgroud_upper_mantle_dislocation['n'] = stress_exponents_for_dislocation_creep.data['background'][0]
    backgroud_upper_mantle_dislocation['m'] = 0.0
    backgroud_upper_mantle_dislocation['E'] = activation_energies_for_dislocation_creep.data['background'][0]
    backgroud_upper_mantle_dislocation['V'] = activation_volumes_for_dislocation_creep.data['background'][0]
    return backgroud_upper_mantle_diffusion, backgroud_upper_mantle_dislocation


def FastZeroStep(Inputs, output_first_step = False):
    '''
    Generate a prm file to run only the 0th step, and real fast
    Inputs:
        output_first_step: if true, then output the 1st step as well.
    '''
    Inputs['Nonlinear solver scheme'] = 'no Advection, no Stokes'
    if output_first_step:
        time_between_graphical_output = None
        try:
            time_between_graphical_output = Inputs['Postprocess']['Visualization']['Time between graphical output']
        except KeyError:
            raise KeyError("input file has to have \'Time between graphical output\' if the option of output_first_step is True")
        Inputs['End time'] = time_between_graphical_output
    else:
        Inputs['End time'] = '0' # end time is 0
        # don't solve it


def TestInitalSteps(Inputs, n_outputs, output_interval):
    '''
    options for generating a case to test the initial steps
    '''
    Inputs['End time'] = str(output_interval * n_outputs)
    # output timing information at every step
    Inputs["Timing output frequency"] = "0"
    # fix post-process section
    pp_dict = Inputs["Postprocess"]
    if "Depth average" in pp_dict:
        # output Depth average at every step
        pp_dict["Depth average"]["Time between graphical output"] = "0"
    if "Visualization" in pp_dict:
        pp_dict["Visualization"]["Time between graphical output"] = str(output_interval)
    if "Particles" in pp_dict:
        pp_dict["Particles"]["Time between data output"] = str(output_interval)
    # fix checkpointing - chekcpoitn eveyr step
    Inputs["Checkpointing"] = {
        "Steps between checkpoint": "1"
    }


class WBFeatureNotFoundError(Exception):
    pass


def FindWBFeatures(Inputs_wb, key):
    '''
    find index of feature in a world builder inputs by its key
    Inputs:
        Inputs_wb (dict): World buider dictionary
        key (str): name of the feature
    '''
    assert(type(Inputs_wb) == dict)
    Features = Inputs_wb['features']
    i = 0
    for feature in Features:
        if feature['name'] == key:
            break
        i += 1
        if i == len(Features):  # not found
            raise WBFeatureNotFoundError("%s: There is no feature named %s" % (Utilities.func_name(), key))
    return i


def RemoveWBFeatures(Inputs_wb, i):
    '''
    remove a feature in World builder dictionary with an index i
    Inputs:
        Inputs_wb (dict): World buider dictionary
        i (int): index of the feature
    '''
    assert(type(Inputs_wb) == dict)
    Outputs_wb = Inputs_wb.copy()
    try:
        Features = Inputs_wb['features']
    except KeyError:
        raise WBFeatureNotFoundError()
    Features.pop(i)
    Outputs_wb['features'] = Features
    return Outputs_wb


def DumpToJson(filein, fileout):
    '''
    Read in an input file and dump the content to a json file
    '''
    assert(os.access(filein, os.R_OK))
    with open(filein, 'r') as fin:
        i_dict = ParseFromDealiiInput(fin)
    with open(fileout, 'w') as fout:
        json.dump(i_dict, fout, indent=2)

#######
# functions and classes for parsing slurm files
#######

def ParseFromSlurmBatchFile(fin):
    '''
    read options from a slurm batch file.
    Note I allow multiple " " in front of each line, this may or may not be a good option.
    '''
    inputs = {}
    inputs["header"] = []
    inputs["config"] = {}
    inputs["load"] = []
    inputs["unload"] = []
    inputs["command"] = []
    inputs["others"] = []
    line = fin.readline()
    while line != "":
        if re.match('^(\t| )*#!', line):
            line1 = re.sub('(\t| )*\n$', '', line)  # eliminate the \n at the end
            inputs["header"].append(line1)
        elif re.match('^(\t| )*#(\t| )*SBATCH', line):
            line1 = re.sub('^(\t| )*#(\t| )*SBATCH ', '', line, count=1)
            key, value = Utilities.ReadDashOptions(line1)
            inputs["config"][key] = value
        elif re.match('(\t| )*module(\t| )*load', line):
            line1 = re.sub('(\t| )*module(\t| )*load(\t| )*', '', line, count=1)
            value = re.sub('(\t| )*\n$', '', line1) 
            inputs["load"].append(value)
        elif re.match('(\t| )*module(\t| )*unload', line):
            line1 = re.sub('(\t| )*module(\t| )*unload(\t| )*', '', line, count=1)
            value = re.sub('(\t| )*\n$', '', line1) 
            inputs["unload"].append(value)
        elif re.match('(\t| )*srun', line):
            temp = re.sub('(\t| )*srun(\t| )*', '', line, count=1)
            inputs["command"].append("srun")
            line1 = re.sub('(\t| )*\n$', '',temp)  # eliminate the \n at the end
            for temp in line1.split(' '):
                if not re.match("^(\t| )*$", temp):
                    # not ' '
                    temp1 = re.sub('^(\t| )*', '', temp) # eliminate ' '
                    value = re.sub('^(\t| )$', '', temp1)
                    inputs["command"].append(value)
        elif re.match('(\t| )*ibrun', line):
            temp = re.sub('(\t| )*ibrun(\t| )*', '', line, count=1)
            inputs["command"].append("ibrun")
            line1 = re.sub('(\t| )*\n$', '',temp)  # eliminate the \n at the end
            for temp in line1.split(' '):
                if not re.match("^(\t| )*$", temp):
                    # not ' '
                    temp1 = re.sub('^(\t| )*', '', temp) # eliminate ' '
                    value = re.sub('^(\t| )$', '', temp1)
                    inputs["command"].append(value)
        else:
            inputs["others"].append(re.sub('(\t| )*\n$', '', line))
            pass
        line = fin.readline()
    return inputs


def ParseToSlurmBatchFile(fout, outputs, **kwargs):
    '''
    export options to a slurm batch file.
    '''
    contents = ""
    # write header
    for header in outputs["header"]:
        contents += (header + "\n")
    # write config
    for key, value in outputs["config"].items():
        if re.match('^--', key):
            contents += ("#SBATCH" + " " + key + "=" + value + "\n")
        elif re.match('^-', key):
            contents += ("#SBATCH" + " " + key + " " + value + "\n")
        else:
            raise ValueError("The format of key (%s) is incorrect" % key)
    # load and unload 
    for module in outputs["unload"]:
        contents += ("module unload %s\n" % module)
    for module in outputs["load"]:
        contents += ("module load %s\n" % module)
    # others
    for line in outputs['others']:
        contents += (line + '\n')
    # command
    is_first = True
    for component in outputs["command"]:
        if is_first:
            is_first = False
        else:
            contents += " "
        contents += component
        if component == "mpirun":
            # append the cpu options after mpirun
            contents += " "
            contents += ("-n " + outputs['config']['-n'])
    fout.write(contents)
    pass


class SLURM_OPERATOR():
    '''
    Class for modifying slurm operation files
    Attributes:
        i_dict (dict): input options
        o_dict (dict): output options
    '''
    def __init__(self, slurm_base_file_path):
        assert(os.path.isfile(slurm_base_file_path))
        with open(slurm_base_file_path, 'r') as fin:
            self.i_dict = ParseFromSlurmBatchFile(fin)
        self.check() # call the check function to check the contents of the file
        self.o_dict = {}
    
    def check(self):
        '''
        check the options in inputs
        '''
        assert('-N' in self.i_dict['config'])
        assert('-n' in self.i_dict['config'])
        assert('--threads-per-core' in self.i_dict['config'])
        assert('--tasks-per-node' in self.i_dict['config'])
        assert('--partition' in self.i_dict['config'])

    def SetAffinity(self, nnode, nthread, nthreads_per_cpu, **kwargs):
        '''
        set options for affinities
        '''
        partition = kwargs.get('partition', None)
        use_mpirun = kwargs.get('use_mpirun', False)
        bind_to = kwargs.get('bind_to', None)
        self.o_dict = copy.deepcopy(self.i_dict)
        self.o_dict['config']['-N'] = str(int(nnode))
        self.o_dict['config']['-n'] = str(nthread)
        self.o_dict['config']['--threads-per-core'] = str(nthreads_per_cpu)
        self.o_dict['config']['--tasks-per-node'] = str(int(nthread//nnode))
        if partition is not None:
             self.o_dict['config']['--partition'] = partition 
        if use_mpirun:
            self.o_dict['command'][0] = 'mpirun'
        if bind_to != None:
            assert(bind_to in ["socket", 'core'])
            if self.o_dict['command'][0] == "mpirun":
                self.o_dict['command'].insert(1, "--bind-to " + bind_to)
            if self.o_dict['command'][0] == "srun":
                self.o_dict['command'].insert(1, "--cpu-bind=" + bind_to + "s")

            

    def SetName(self, _name):
        '''
        set the name of the job
        '''
        self.o_dict['config']['--job-name'] = _name

    def SetCommand(self, build_directory, prm_file):
        '''
        Set the command to use
        '''
        if build_directory != "":
            self.o_dict['command'][-2] = "${ASPECT_SOURCE_DIR}/build_%s/aspect" % build_directory
        else:
            self.o_dict['command'][-2] = "${ASPECT_SOURCE_DIR}/build/aspect"
        self.o_dict['command'][-1] = prm_file

    def __call__(self, slurm_file_path):
        with open(slurm_file_path, 'w') as fout:
            ParseToSlurmBatchFile(fout, self.o_dict)
        print("Slurm file created: ", slurm_file_path)

    pass


class SLURM_OPT(Utilities.JSON_OPT):
    def __init__(self):
        '''
        Initiation, first perform parental class's initiation,
        then perform daughter class's initiation.
        '''
        Utilities.JSON_OPT.__init__(self)
        self.add_key("Slurm file (inputs)", str, ["slurm file"], "slurm.sh", nick='slurm_base_file')
        self.add_key("Openmpi version", str, ["openmpi version"], "", nick='openmpi')
        self.add_key("build directory", str, ["build directory"], "", nick="build_directory")
        self.add_key("Flag", str, ["flag"], "", nick="flag")
        self.add_key("Tasks per node", int, ["tasks per node"], 32, nick='tasks_per_node')
        self.add_key("cpus", int, ["cpus"], 1, nick="cpus")
        self.add_key("Threads per cpu", int, ["threads per cpu"], 1, nick='threads_per_cpu')
        self.add_key("List of nodes", list, ["node list"], [], nick="nodelist")
        self.add_key("Path to the prm file", str, ["prm path"], "./case.prm", nick="prm_path")
        self.add_key("Output directory", str, ["output directory"], ".", nick="output_directory")
        self.add_key("Base directory", str, ["base directory"], ".", nick="base_directory")

    def check(self):
        slurm_base_path = self.values[0]
        os.path.isfile(slurm_base_path)
        prm_path = self.values[8]
        os.path.isfile(prm_path)
        base_directory = self.values[10]
        assert(os.path.isdir(base_directory))
        output_directory = self.values[9]
        if not os.path.isdir(output_directory):
            os.mkdir(output_directory)

    def get_base_path(self):
        '''
        get the path to the base file (i.e job_p-billen.sh)
        '''
        base_directory = self.values[10]
        slurm_base_file = self.values[0]
        slurm_base_path = os.path.join(base_directory, slurm_base_file)
        return slurm_base_path

    def to_set_affinity(self):
        tasks_per_node = self.values[4]
        cpus = self.values[5]
        threads_per_cpu = self.values[6]
        threads = int(cpus * threads_per_cpu)
        nnode = int(np.ceil(threads / tasks_per_node))
        return nnode, threads, threads_per_cpu

    def to_set_command(self):
        build_directory = self.values[2]
        prm_path = self.values[8]
        return build_directory, prm_path

    def get_job_name(self):
        '''
        use the basename of the output directory as the job name
        '''
        output_directory = self.values[9]
        job_name = os.path.basename(output_directory)
        return job_name
        pass

    def get_output_path(self):
        '''
        get the path of the output file (i.e. job_p-billen.sh)
        '''
        slurm_base_path = self.values[0]
        output_directory = self.values[9]
        output_path = os.path.join(output_directory, os.path.basename(slurm_base_path))
        return output_path

    def fix_base_dir(self, base_directory):
        self.values[10] = base_directory
    
    def fix_output_dir(self, output_directory):
        self.values[9] = output_directory
    
    def get_output_dir(self):
        output_directory = self.values[9]
        return output_directory

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
    if _commend == 'fast_zero_step':
        # Generate a prm file to run only the 0th step, and real fast
        _path = arg.inputs
        out_path = arg.outputs
        inputs = ReadPrmFile(_path)
        FastZeroStep(inputs)
        WritePrmFile(out_path, inputs)
    
    if _commend == 'dump_to_json':
        DumpToJson(arg.inputs, arg.outputs)

# run script
if __name__ == '__main__':
    main()
