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
        python shilofue/ParsePrm.py fast_zero_step -i case.prm -o case_test.prm

descriptions:
    copy and pasted all the function in the origin PARSE_OPERATION class
""" 
import numpy as np
import sys, os, argparse
import re
# import pathlib
# import subprocess
import numpy as np
# from matplotlib import cm
# from matplotlib import pyplot as plt
from shilofue.Parse import COMPOSITION

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')


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


def ReadPrmFile(_path):
    """
    simply a wrapper for the parse function
    Inputs:
        _path: path for a prm file
    """
    assert(os.access(_path, os.R_OK))
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
    backgroud_upper_mantle_dislocation['m'] = 1.0
    backgroud_upper_mantle_dislocation['E'] = activation_energies_for_dislocation_creep.data['background'][0] 
    backgroud_upper_mantle_dislocation['V'] = activation_volumes_for_dislocation_creep.data['background'][0]
    return backgroud_upper_mantle_diffusion, backgroud_upper_mantle_dislocation


def FastZeroStep(Inputs):
    '''
    Generate a prm file to run only the 0th step, and real fast
    todo
    '''
    Inputs['End time'] = '0' # end time is 0
    # limit the number of non-linear iterations to only 1
    Inputs['Max nonlinear iterations'] = '1'
    Inputs['Max nonlinear iterations in pre-refinement'] = '0'
    # make the linear tolerance really big
    Inputs['Solver parameters']['Stokes solver parameters'] = \
    {'Number of cheap Stokes solver steps': '0', 'Linear solver tolerance': '0.9999'}


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


# run script
if __name__ == '__main__':
    main()