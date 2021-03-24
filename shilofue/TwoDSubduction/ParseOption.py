# -*- coding: utf-8 -*-
r"""(one line description)
future
This exports: 

  -

This depends on:

  -  

Examples of usage:

  - default usage:

        python -m 

descriptions
""" 
import numpy as np
import sys, os, argparse
# import json, re
# import pathlib
# import subprocess
import numpy as np
# from matplotlib import cm
from matplotlib import pyplot as plt

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')

def SomeFunction(foo):
    '''
    descriptions
    Inputs:
        -
    Returns:
        -
    '''
    pass

def LowerMantle(Inputs, _config):
    """
    calculate flow law parameters
    """
    _type = _config.get('model_type', 0)
    my_assert(type(_type) == int, TypeError, "Type of input \'model_type\' must be int")
    if _type == 0:
        # when phase transition only happens on mantle composition
        LowerMantle0(Inputs, _config)
    elif _type == 1:
        # when phase transition only happens on all compositions
        # There is a eclogite transition of crustal layer
        LowerMantle1(Inputs, _config)
    elif _type == 2:
        # setup lower mantle rheology for phase transition in a CDPT model
        # There is a eclogite transition of crustal layer
        LowerMantle2(Inputs, _config)
    
def LowerMantle0(Inputs, _config):
    """
    calculate flow law parameters, when phase transition only happens on mantle composition
    """
    # parse from input
    jump = _config['upper_lower_viscosity']
    T = _config['T660']
    P = _config['P660']
    V1 = _config['LowerV']
    visco_plastic = Inputs["Material model"]['Visco Plastic']
    prefactors_for_diffusion_creep = Parse.COMPOSITION(visco_plastic["Prefactors for diffusion creep"])
    grain_size = float(visco_plastic["Grain size"])
    grain_size_exponents_for_diffusion_creep  = Parse.COMPOSITION(visco_plastic["Grain size exponents for diffusion creep"])
    activation_energies_for_diffusion_creep = Parse.COMPOSITION(visco_plastic["Activation energies for diffusion creep"])
    activation_volumes_for_diffusion_creep  = Parse.COMPOSITION(visco_plastic["Activation volumes for diffusion creep"])
    # call GetLowerMantleRheology to derive parameters for lower mantle flow law 
    backgroud_upper_mantle_diffusion = {}
    backgroud_upper_mantle_diffusion['A'] = prefactors_for_diffusion_creep.data['background'][0] 
    backgroud_upper_mantle_diffusion['d'] = grain_size
    backgroud_upper_mantle_diffusion['n'] = 1.0 
    backgroud_upper_mantle_diffusion['m'] = grain_size_exponents_for_diffusion_creep.data['background'][0] 
    backgroud_upper_mantle_diffusion['E'] = activation_energies_for_diffusion_creep.data['background'][0] 
    backgroud_upper_mantle_diffusion['V'] = activation_volumes_for_diffusion_creep.data['background'][0] 
    backgroud_lower_mantle_diffusion = Rheology.GetLowerMantleRheology(backgroud_upper_mantle_diffusion, jump, T, P, V1=V1, strategy='A')
    # future: add in choice of phases
    prefactors_for_diffusion_creep.data['background'] = [backgroud_upper_mantle_diffusion['A'], backgroud_lower_mantle_diffusion['A']]
    grain_size_exponents_for_diffusion_creep.data['background'] = [backgroud_upper_mantle_diffusion['m'], backgroud_lower_mantle_diffusion['m']]
    activation_energies_for_diffusion_creep.data['background'] = [backgroud_upper_mantle_diffusion['E'], backgroud_lower_mantle_diffusion['E']]
    activation_volumes_for_diffusion_creep.data['background'] = [backgroud_upper_mantle_diffusion['V'], backgroud_lower_mantle_diffusion['V']] 
    # parse back
    visco_plastic["Prefactors for diffusion creep"] = prefactors_for_diffusion_creep.parse_back()
    visco_plastic["Grain size exponents for diffusion creep"] = grain_size_exponents_for_diffusion_creep.parse_back()
    visco_plastic["Activation energies for diffusion creep"] = activation_energies_for_diffusion_creep.parse_back()
    visco_plastic["Activation volumes for diffusion creep"] = activation_volumes_for_diffusion_creep.parse_back()
    return Inputs

def LowerMantle1(Inputs, _config):
    """
    calculate flow law parameters, when phase transition only happens on all compositions
    There is a eclogite transition of crustal layer
    """
    # parse from input
    jump = _config['upper_lower_viscosity']
    T = _config['T660']
    P = _config['P660']
    V1 = _config['LowerV']
    visco_plastic = Inputs["Material model"]['Visco Plastic']
    prefactors_for_diffusion_creep = Parse.COMPOSITION(visco_plastic["Prefactors for diffusion creep"])
    grain_size = float(visco_plastic["Grain size"])
    grain_size_exponents_for_diffusion_creep  = Parse.COMPOSITION(visco_plastic["Grain size exponents for diffusion creep"])
    activation_energies_for_diffusion_creep = Parse.COMPOSITION(visco_plastic["Activation energies for diffusion creep"])
    activation_volumes_for_diffusion_creep  = Parse.COMPOSITION(visco_plastic["Activation volumes for diffusion creep"])
    # call GetLowerMantleRheology to derive parameters for lower mantle flow law 
    backgroud_upper_mantle_diffusion = {}
    backgroud_upper_mantle_diffusion['A'] = prefactors_for_diffusion_creep.data['background'][0] 
    backgroud_upper_mantle_diffusion['d'] = grain_size
    backgroud_upper_mantle_diffusion['n'] = 1.0 
    backgroud_upper_mantle_diffusion['m'] = grain_size_exponents_for_diffusion_creep.data['background'][0] 
    backgroud_upper_mantle_diffusion['E'] = activation_energies_for_diffusion_creep.data['background'][0] 
    backgroud_upper_mantle_diffusion['V'] = activation_volumes_for_diffusion_creep.data['background'][0] 
    backgroud_lower_mantle_diffusion = Rheology.GetLowerMantleRheology(backgroud_upper_mantle_diffusion, jump, T, P, V1=V1, strategy='A')
    # future: add in choice of phases
    prefactors_for_diffusion_creep.data['background'] = [backgroud_upper_mantle_diffusion['A'], backgroud_lower_mantle_diffusion['A']]
    grain_size_exponents_for_diffusion_creep.data['background'] = [backgroud_upper_mantle_diffusion['m'], backgroud_lower_mantle_diffusion['m']]
    activation_energies_for_diffusion_creep.data['background'] = [backgroud_upper_mantle_diffusion['E'], backgroud_lower_mantle_diffusion['E']]
    activation_volumes_for_diffusion_creep.data['background'] = [backgroud_upper_mantle_diffusion['V'], backgroud_lower_mantle_diffusion['V']] 
    # spcrust
    prefactors_for_diffusion_creep.data['spcrust'][2] = backgroud_lower_mantle_diffusion['A']
    grain_size_exponents_for_diffusion_creep.data['spcrust'][2] = backgroud_lower_mantle_diffusion['m']
    activation_energies_for_diffusion_creep.data['spcrust'][2] = backgroud_lower_mantle_diffusion['E']
    activation_volumes_for_diffusion_creep.data['spcrust'][2] = backgroud_lower_mantle_diffusion['V']
    # harz
    prefactors_for_diffusion_creep.data['spharz'][1] = backgroud_lower_mantle_diffusion['A']
    grain_size_exponents_for_diffusion_creep.data['spharz'][1] = backgroud_lower_mantle_diffusion['m']
    activation_energies_for_diffusion_creep.data['spharz'][1] = backgroud_lower_mantle_diffusion['E']
    activation_volumes_for_diffusion_creep.data['spharz'][1] = backgroud_lower_mantle_diffusion['V']
    # parse back
    visco_plastic["Prefactors for diffusion creep"] = prefactors_for_diffusion_creep.parse_back()
    visco_plastic["Grain size exponents for diffusion creep"] = grain_size_exponents_for_diffusion_creep.parse_back()
    visco_plastic["Activation energies for diffusion creep"] = activation_energies_for_diffusion_creep.parse_back()
    visco_plastic["Activation volumes for diffusion creep"] = activation_volumes_for_diffusion_creep.parse_back()
    return Inputs

def LowerMantle2(Inputs, _config):
    """
    calculate flow law parameters of the lower mantle for CDPT model
    There is a eclogite transition of crustal layer
    """
    index_660 = 4
    index_all = 8
    index_660_crust = 2
    index_all_crust = 4
    # parse from input
    jump = _config['upper_lower_viscosity']
    T = _config['T660']
    P = _config['P660']
    V1 = _config['LowerV']
    visco_plastic = Inputs["Material model"]['Visco Plastic']
    prefactors_for_diffusion_creep = Parse.COMPOSITION(visco_plastic["Prefactors for diffusion creep"])
    grain_size = float(visco_plastic["Grain size"])
    grain_size_exponents_for_diffusion_creep  = Parse.COMPOSITION(visco_plastic["Grain size exponents for diffusion creep"])
    activation_energies_for_diffusion_creep = Parse.COMPOSITION(visco_plastic["Activation energies for diffusion creep"])
    activation_volumes_for_diffusion_creep  = Parse.COMPOSITION(visco_plastic["Activation volumes for diffusion creep"])
    # call GetLowerMantleRheology to derive parameters for lower mantle flow law 
    backgroud_upper_mantle_diffusion = {}
    backgroud_upper_mantle_diffusion['A'] = prefactors_for_diffusion_creep.data['background'][0] 
    backgroud_upper_mantle_diffusion['d'] = grain_size
    backgroud_upper_mantle_diffusion['n'] = 1.0 
    backgroud_upper_mantle_diffusion['m'] = grain_size_exponents_for_diffusion_creep.data['background'][0] 
    backgroud_upper_mantle_diffusion['E'] = activation_energies_for_diffusion_creep.data['background'][0] 
    backgroud_upper_mantle_diffusion['V'] = activation_volumes_for_diffusion_creep.data['background'][0] 
    backgroud_lower_mantle_diffusion = Rheology.GetLowerMantleRheology(backgroud_upper_mantle_diffusion, jump, T, P, V1=V1, strategy='A')
    # background (pyrolite)
    prefactors_for_diffusion_creep.data['background'] = []
    grain_size_exponents_for_diffusion_creep.data['background'] = []
    activation_energies_for_diffusion_creep.data['background'] = []
    activation_volumes_for_diffusion_creep.data['background'] = []
    for i in range(index_660):
        prefactors_for_diffusion_creep.data['background'].append(backgroud_upper_mantle_diffusion['A'])
        grain_size_exponents_for_diffusion_creep.data['background'].append(backgroud_upper_mantle_diffusion['m'])
        activation_energies_for_diffusion_creep.data['background'].append(backgroud_upper_mantle_diffusion['E'])
        activation_volumes_for_diffusion_creep.data['background'].append(backgroud_upper_mantle_diffusion['V'])
    for i in range(index_660, index_all):
        prefactors_for_diffusion_creep.data['background'].append(backgroud_lower_mantle_diffusion['A'])
        grain_size_exponents_for_diffusion_creep.data['background'].append(backgroud_lower_mantle_diffusion['m'])
        activation_energies_for_diffusion_creep.data['background'].append(backgroud_lower_mantle_diffusion['E'])
        activation_volumes_for_diffusion_creep.data['background'].append(backgroud_lower_mantle_diffusion['V'])
    # spcrust
    prefactors_for_diffusion_creep.data['spcrust'] = [prefactors_for_diffusion_creep.data['spcrust'][0]]
    grain_size_exponents_for_diffusion_creep.data['spcrust'] = [grain_size_exponents_for_diffusion_creep.data['spcrust'][0]]
    activation_energies_for_diffusion_creep.data['spcrust'] = [activation_energies_for_diffusion_creep.data['spcrust'][0]]
    activation_volumes_for_diffusion_creep.data['spcrust'] = [activation_volumes_for_diffusion_creep.data['spcrust'][0]]
    for i in range(1, index_660_crust):
        # index 0 is the shear zone, it shouldn't be changed here
        prefactors_for_diffusion_creep.data['spcrust'].append(backgroud_upper_mantle_diffusion['A'])
        grain_size_exponents_for_diffusion_creep.data['spcrust'].append(backgroud_upper_mantle_diffusion['m'])
        activation_energies_for_diffusion_creep.data['spcrust'].append(backgroud_upper_mantle_diffusion['E'])
        activation_volumes_for_diffusion_creep.data['spcrust'].append(backgroud_upper_mantle_diffusion['V'])
    for i in range(index_660_crust, index_all_crust):
        prefactors_for_diffusion_creep.data['spcrust'].append(backgroud_lower_mantle_diffusion['A'])
        grain_size_exponents_for_diffusion_creep.data['spcrust'].append(backgroud_lower_mantle_diffusion['m'])
        activation_energies_for_diffusion_creep.data['spcrust'].append(backgroud_lower_mantle_diffusion['E'])
        activation_volumes_for_diffusion_creep.data['spcrust'].append(backgroud_lower_mantle_diffusion['V'])
    # harz
    prefactors_for_diffusion_creep.data['spharz'] = []
    grain_size_exponents_for_diffusion_creep.data['spharz'] = []
    activation_energies_for_diffusion_creep.data['spharz'] = []
    activation_volumes_for_diffusion_creep.data['spharz'] = []
    for i in range(index_660):
        prefactors_for_diffusion_creep.data['spharz'].append(backgroud_upper_mantle_diffusion['A'])
        grain_size_exponents_for_diffusion_creep.data['spharz'].append(backgroud_upper_mantle_diffusion['m'])
        activation_energies_for_diffusion_creep.data['spharz'].append(backgroud_upper_mantle_diffusion['E'])
        activation_volumes_for_diffusion_creep.data['spharz'].append(backgroud_upper_mantle_diffusion['V'])
    for i in range(index_660, index_all):
        prefactors_for_diffusion_creep.data['spharz'].append(backgroud_lower_mantle_diffusion['A'])
        grain_size_exponents_for_diffusion_creep.data['spharz'].append(backgroud_lower_mantle_diffusion['m'])
        activation_energies_for_diffusion_creep.data['spharz'].append(backgroud_lower_mantle_diffusion['E'])
        activation_volumes_for_diffusion_creep.data['spharz'].append(backgroud_lower_mantle_diffusion['V'])
  
    # parse back
    visco_plastic["Prefactors for diffusion creep"] = prefactors_for_diffusion_creep.parse_back()
    visco_plastic["Grain size exponents for diffusion creep"] = grain_size_exponents_for_diffusion_creep.parse_back()
    visco_plastic["Activation energies for diffusion creep"] = activation_energies_for_diffusion_creep.parse_back()
    visco_plastic["Activation volumes for diffusion creep"] = activation_volumes_for_diffusion_creep.parse_back()
    return Inputs
    
    
def Particle(Inputs, _config):
    """
    Define a particle section in the .prm file
    """
    # First, check the option in json file
    if_particle_in_slab = _config.get("if_particle_in_slab", 0)
    if if_particle_in_slab == 0:
        return
    # List of postprocessors
    list_postproc = Inputs['Postprocess']['List of postprocessors']
    if re.match('particles', list_postproc) is None:
        list_postproc += ', particles'
    Inputs['Postprocess']['List of postprocessors'] = list_postproc
    # initiate a dictionary for particle settings
    p_dict = {}
    p_dict['Particle generator name'] = 'ascii file'
    p_dict['Data output format'] = 'vtu'
    p_dict['List of particle properties'] = 'initial position, pT path'  # particle properties
    # initiate a new dictionary for the generator subsection
    p_generator_dict = {}
    p_generator_dict['Ascii file'] = {'Data directory': './', 'Data file name': 'particle.dat'}
    p_dict['Generator'] = p_generator_dict
    # assign it to Inputs
    Inputs['Postprocess']['Particles'] = p_dict

def Phases(Inputs, _config):
    '''
    Define phase transition variables
    '''
    visco_plastic = Inputs["Material model"]['Visco Plastic']
    try:
        eclogite_transition = visco_plastic["Eclogite transition"]
    except KeyError:
        pass
    else:
        # temperature for eclogite transition 
        try:
            temperature_for_eclogite_transition = _config["temperature_for_eclogite_transition"]
        except KeyError:
            pass
        else:
            eclogite_transition["Temperature for eclogite transition"] = str(temperature_for_eclogite_transition)
        # temperature width for eclogite transition 
        try:
            temperature_width_for_eclogite_transition = _config["temperature_width_for_eclogite_transition"]
        except KeyError:
            pass
        else:
            eclogite_transition["Temperature width for eclogite transition"] = str(temperature_width_for_eclogite_transition)
       
        visco_plastic["Eclogite transition"] = eclogite_transition
    try:
        # decoupled eclogite transition
        decoupling_eclogite_transition = _config["decoupling_eclogite_transition"]
    except KeyError:
        pass
    else:
        if decoupling_eclogite_transition == 1:
            visco_plastic["Decoupling eclogite viscosity"] = "true"
        else:
            visco_plastic["Decoupling eclogite viscosity"] = "false"
        
    Inputs["Material model"]['Visco Plastic'] = visco_plastic
    return Inputs

def InitialTemperature(Inputs, _config):
    """
    Initial temperature
    """
    # get the initial temperature module 
    initial_temperature = Inputs["Initial temperature model"]
    
    # subsection Subduction 2d temperature
    try:
        thermal_boundary_width_factor_in = _config["slab_thermal_boundary_width_factor_in"]
        initial_temperature["Subduction 2d temperature"]["Thermal boundary width factor in"] = str(thermal_boundary_width_factor_in)
    except KeyError:
        pass
    try:
        thermal_boundary_width_factor_out = _config["slab_thermal_boundary_width_factor_out"]
        initial_temperature["Subduction 2d temperature"]["Thermal boundary width factor out"] = str(thermal_boundary_width_factor_out)
    except KeyError:
        pass
    
    Inputs["Initial temperature model"] = initial_temperature


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
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    # commands
    if _commend == 'foo':
        # example:
        SomeFunction('foo')

# run script
if __name__ == '__main__':
    main()