
r"""Create cases for TwoDSubduction project

This exports: 

  -

This depends on:

  -  

Examples of usage:

  - print parameterization of lower mantle rheology:

        python -m shilofue.TwoDSubduction0.CreateCases
lower_mantle -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear30/eba/case1.prm 
-j files/TwoDSubduction/eba_lower_mantle.json
-o /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear30/eba_intial_T/case_o.prm

descriptions
""" 
import numpy as np
import sys, os, argparse
import json #, re
# import pathlib
# import subprocess
import numpy as np
# from matplotlib import cm
from matplotlib import pyplot as plt
import shilofue.Parse as Parse
import shilofue.ParsePrm as ParsePrm
from shilofue.Rheology import GetLowerMantleRheology
from shilofue.Utilities import my_assert, ggr2cart2, cart2sph2

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']


class MY_PARSE_OPERATIONS(Parse.PARSE_OPERATIONS):
    """
    put parse operations in a single class
    inherit from the general class
    Attributes:
        ALL_OPERATIONS(list):
            all avalable operations
    """
    def __init__(self):
        """
        Initiation
        """
        Parse.PARSE_OPERATIONS.__init__(self)
        self.ALL_OPERATIONS += ["LowerMantle", "Particle", "Phases", "InitialTemperature"]
    
    def LowerMantle(self, Inputs, _config):
        """
        calculate flow law parameters
        """
        _type = _config.get('model_type', 0)
        my_assert(type(_type) == int, TypeError, "Type of input \'model_type\' must be int")
        if _type == 0:
            # when phase transition only happens on mantle composition
            LowerMantle0(Inputs, _config)
        if _type == 1:
            # when phase transition only happens on all compositions
            # There is a eclogite transition of crustal layer
            LowerMantle1(Inputs, _config)


    def Particle(self, Inputs, _config):
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
    
    def Phases(self, Inputs, _config):
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
    
    def InitialTemperature(self, Inputs, _config):
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
    backgroud_lower_mantle_diffusion = GetLowerMantleRheology(backgroud_upper_mantle_diffusion, jump, T, P, V1=V1, strategy='A')
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
    return Inputs, backgroud_lower_mantle_diffusion


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
    backgroud_lower_mantle_diffusion = GetLowerMantleRheology(backgroud_upper_mantle_diffusion, jump, T, P, V1=V1, strategy='A')
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
    return Inputs, backgroud_lower_mantle_diffusion


def LowerMantle2(Inputs, _config):
    """
    calculate flow law parameters, when phase transition are described by the CDPT model
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
    backgroud_lower_mantle_diffusion = GetLowerMantleRheology(backgroud_upper_mantle_diffusion, jump, T, P, V1=V1, strategy='A')
    # write them back
    # background composition
    upper_mantle_index = 4 
    total_index = 8
    for i in range(upper_mantle_index, total_index):
        prefactors_for_diffusion_creep.data['background'][i] = backgroud_lower_mantle_diffusion['A']
        grain_size_exponents_for_diffusion_creep.data['background'][i] = backgroud_lower_mantle_diffusion['m']
        activation_energies_for_diffusion_creep.data['background'][i] = backgroud_lower_mantle_diffusion['E']
        activation_volumes_for_diffusion_creep.data['background'][i] = backgroud_lower_mantle_diffusion['V']
    # spcrust composition
    upper_mantle_index = 2
    total_index = 4
    for i in range(upper_mantle_index, total_index):
        prefactors_for_diffusion_creep.data['spcrust'][i] = backgroud_lower_mantle_diffusion['A']
        grain_size_exponents_for_diffusion_creep.data['spcrust'][i] = backgroud_lower_mantle_diffusion['m']
        activation_energies_for_diffusion_creep.data['spcrust'][i] = backgroud_lower_mantle_diffusion['E']
        activation_volumes_for_diffusion_creep.data['spcrust'][i] = backgroud_lower_mantle_diffusion['V']
    # spharz composition
    upper_mantle_index = 4
    total_index = 8
    for i in range(upper_mantle_index, total_index):
        prefactors_for_diffusion_creep.data['spharz'][i] = backgroud_lower_mantle_diffusion['A']
        grain_size_exponents_for_diffusion_creep.data['spharz'][i] = backgroud_lower_mantle_diffusion['m']
        activation_energies_for_diffusion_creep.data['spharz'][i] = backgroud_lower_mantle_diffusion['E']
        activation_volumes_for_diffusion_creep.data['spharz'][i] = backgroud_lower_mantle_diffusion['V']
    # parse back
    visco_plastic["Prefactors for diffusion creep"] = prefactors_for_diffusion_creep.parse_back()
    visco_plastic["Grain size exponents for diffusion creep"] = grain_size_exponents_for_diffusion_creep.parse_back()
    visco_plastic["Activation energies for diffusion creep"] = activation_energies_for_diffusion_creep.parse_back()
    visco_plastic["Activation volumes for diffusion creep"] = activation_volumes_for_diffusion_creep.parse_back()
    return Inputs, backgroud_lower_mantle_diffusion



class MYCASE(Parse.CASE):
    '''
    Inherit from class CASE in Parse.py
    '''
    def process_particle_data(self):
        '''
        process the coordinates of particle, doing nothing here.
        Reload here to add particles in the crust
        '''
        # all configurations
        _config = { **self.config, **self.test, **self.extra }
        # import values
        slab_phi_c = _config.get("slab_phi_c", 0.628319)
        R0 = _config.get("R0", 6.371e6)
        Rc = _config.get("Rc", 4.0e5)
        slab_to = _config.get("slab_to", 2.0e5)
        depth_particle_in_slab = _config.get("depth_particle_in_slab", 100.0)
        number_particle_in_slab = _config.get("number_particle_in_slab", 1000)
        my_assert(type(number_particle_in_slab)==int, TypeError, "number_particle_in_slab must be an int value")
        # initiate particle_data
        self.particle_data = np.zeros((number_particle_in_slab, 2))
        # get phi value at the tip of initial slab
        phi_st = slab_phi_c + (2 * Rc * slab_to - slab_to**2.0)**0.5 / R0
        for i in range(number_particle_in_slab):
            # get particle coordinates
            # First, angle is divided uniformly.
            # Then, radius is computed accordingly.
            phi = i * phi_st / number_particle_in_slab
            if phi < slab_phi_c:
                r = R0 - depth_particle_in_slab
            else:
                r = R0 + (Rc**2.0 - R0**2.0 * (phi - slab_phi_c)**2.0)**0.5  - Rc - depth_particle_in_slab
            # apply transform to cartisian coordinates
            x, y = ggr2cart2(phi, r)
            # assign value in particle_data
            self.particle_data[i, 0] = x
            self.particle_data[i, 1] = y


def SomeFunction(foo):
    '''
    descriptions
    Inputs:
        -
    Returns:
        -
    '''
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
                        help='Some inputs(prm file)')
    parser.add_argument('-j', '--json', type=str,
                        default='',
                        help='json file inputs')
    parser.add_argument('-o', '--outputs', type=str,
                        default='case_o.prm',
                        help='Some outputs(prm file)')
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    # commands
    if _commend == 'lower_mantle':
        # derive parameterization of lower mantle from a prm file and configurations
        # check file existence
        assert(os.access(arg.inputs, os.R_OK))
        assert(os.access(arg.json, os.R_OK))
        
        # read parameters and configuration
        with open(arg.inputs, 'r') as fin:
            _inputs = ParsePrm.ParseFromDealiiInput(fin)
        with open(arg.json, 'r') as fin:
            _config = json.load(fin)
        outputs, lower_mantle_diffusion = LowerMantle2(_inputs, _config)
        # screen output
        print(lower_mantle_diffusion)
        # file output
        with open(arg.outputs, 'w') as fout:
            ParsePrm.ParseToDealiiInput(fout, outputs)

# run script
if __name__ == '__main__':
    main()