
r"""Create cases for TwoDSubduction project

This exports: 

  -

This depends on:

  -  

Examples of usage:

  - print parameterization of lower mantle rheology:

        python -m shilofue.TwoDSubduction0.Parse
lower_mantle -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear30/eba/case1.prm 
-j files/TwoDSubduction/eba_lower_mantle.json
-o /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear30/eba_intial_T/case_o.prm

  - setup rheology in prm file
        python -m shilofue.TwoDSubduction0.Parse rheology_cdpt
-i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear31/eba_test_wet_mod/case.prm
-j /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear31/eba_test_wet_mod/rheology.json 
-o /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear31/eba_test_wet_mod/case_o.prm
    
    in case the lower mantle rheology is read in from json file as well:
        -b 1

  - init new case:
        python -m shilofue.TwoDSubduction0.Parse create -j ~/ASPECT_PROJECT/TwoDSubduction/non_linear32/init.json

descriptions
    available configuretions of case generating
    available options of of case generating:
        max_refinement: change the maximum level of refinemment, appendix: '_MRf%d'
""" 
import numpy as np
import sys, os, argparse
import json #, re
# import pathlib
# import subprocess
import numpy as np
# from matplotlib import cm
from matplotlib import pyplot as plt
from shutil import rmtree, copy2, copytree
from pathlib import Path 
import shilofue.Parse as Parse
import shilofue.ParsePrm as ParsePrm
from shilofue.Rheology import GetLowerMantleRheology, CreepRheologyInAspectViscoPlastic, ComputeComposite
from shilofue.Cases import CASE

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']

# import utilities
sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities


def Usage():
    print("\
This scripts parses inputs and outputs for cases in TwoDSubduction project\n\
\n\
Examples of usage: \n\\n\
  - print parameterization of lower mantle rheology:\n\
\n\
        python -m shilofue.TwoDSubduction0.Parse\n\
lower_mantle -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear30/eba/case1.prm \n\
-j files/TwoDSubduction/eba_lower_mantle.json\n\
-o /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear30/eba_intial_T/case_o.prm\n\
\n\
  - setup rheology in prm file\n\
\n\
        python -m shilofue.TwoDSubduction0.Parse rheology_cdpt\n\
-i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear31/eba_test_wet_mod/case.prm\n\
-j /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear31/eba_test_wet_mod/rheology.json \n\
-o /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear31/eba_test_wet_mod/case_o.prm\n\
\n\
    in case the lower mantle rheology is read in from json file as well:\n\
        -b 1\n\
\n\
  - init new case:\n\
\n\
        python -m shilofue.TwoDSubduction0.Parse create -j ~/ASPECT_PROJECT/TwoDSubduction/non_linear32/init.json\n\
\n\
  - convert a case in cartesan geometry to one in spherical(chunk) geometry and vice versa\n\
\n\
        Lib_TwoDSubduction0_Parse convert_sph_cart -i ~/ASPECT_PROJECT/TwoDSubduction/wb_sd_issue/wb_cart_cdd50\n\
        -o ~/ASPECT_PROJECT/TwoDSubduction/wb_sd_issue/wb_sph_cdd50\n\
\n\
        default is to transform from cartesian to spherical\n\
        ")


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
        Utilities.my_assert(type(_type) == int, TypeError, "Type of input \'model_type\' must be int")
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


def RheologyCDPT(Inputs, _config, **kwargs):
    """
    calculate flow law parameters, when phase transition only happens on mantle composition
    """
    # creep
    strain_rate = _config['strain_rate']
    diff = _config['diffusion_creep']
    disl = _config['dislocation_creep']
    # shear zone
    shear_zone = _config['shear_zone']
    spcrust_diff = shear_zone['viscosity']
    spcrust_disl = 1e31
    # lower mantle
    read_lower_mantle_full = kwargs.get('lower_mantle_full', 0)
    if read_lower_mantle_full == 1:
        diff_lm = _config['diffusion_lm']
    else:
        lower_mantle = _config['lower_mantle']
        strategy = lower_mantle.get('strategy', 'composite')  # see comments in json file
        jump = lower_mantle['upper_lower_viscosity']
        T = lower_mantle['T660']
        P = lower_mantle['P660']
        V1 = lower_mantle['LowerV']
        eta_diff_660 = CreepRheologyInAspectViscoPlastic(diff, strain_rate, P, T)
        eta_disl_660 = CreepRheologyInAspectViscoPlastic(disl, strain_rate, P, T)
        eta_660 = ComputeComposite(eta_diff_660, eta_disl_660)
        diff_lm = GetLowerMantleRheology(diff, jump, T, P, V1=V1, strategy=strategy, eta=eta_660)
    disl_lm = {'A': 5e-32, 'E': 0.0, 'V': 0.0, 'n': 1.0, 'm':0.0}
    # fix input file
    visco_plastic = Inputs["Material model"]['Visco Plastic TwoD']
    prefactors_for_diffusion_creep = Parse.COMPOSITION()
    grain_size_exponents_for_diffusion_creep = Parse.COMPOSITION()
    activation_energies_for_diffusion_creep = Parse.COMPOSITION()
    activation_volumes_for_diffusion_creep = Parse.COMPOSITION()
    # diffusion creep
    # background composition
    background_up = 4
    background_low = 4
    prefactors_for_diffusion_creep.data['background'] = [diff['A'] for i in range(background_up)] + [diff_lm['A'] for i in range(background_low)]
    grain_size_exponents_for_diffusion_creep.data['background'] = [diff['m'] for i in range(background_up)] + [diff_lm['m'] for i in range(background_low)]
    activation_energies_for_diffusion_creep.data['background'] = [diff['E'] for i in range(background_up)] + [diff_lm['E'] for i in range(background_low)]
    activation_volumes_for_diffusion_creep.data['background'] = [diff['V'] for i in range(background_up)] + [diff_lm['V'] for i in range(background_low)]
    # spcrust
    spcrust_sz = 1
    spcrust_up = 1
    spcrust_low = 2
    prefactors_for_diffusion_creep.data['spcrust'] = [0.5/spcrust_diff for i in range(spcrust_sz)]\
    + [diff['A'] for i in range(spcrust_up)] + [diff_lm['A'] for i in range(spcrust_low)]
    grain_size_exponents_for_diffusion_creep.data['spcrust'] = [0.0 for i in range(spcrust_sz)]\
    + [diff['m'] for i in range(spcrust_up)] + [diff_lm['m'] for i in range(spcrust_low)]
    activation_energies_for_diffusion_creep.data['spcrust'] = [0.0 for i in range(spcrust_sz)]\
    + [diff['E'] for i in range(spcrust_up)] + [diff_lm['E'] for i in range(spcrust_low)]
    activation_volumes_for_diffusion_creep.data['spcrust'] = [0.0 for i in range(spcrust_sz)]\
    + [diff['V'] for i in range(spcrust_up)] + [diff_lm['V'] for i in range(spcrust_low)]
    # spharz
    spharz_up = 4
    spharz_low = 4
    prefactors_for_diffusion_creep.data['spharz'] = [diff['A'] for i in range(spharz_up)] + [diff_lm['A'] for i in range(spharz_low)]
    grain_size_exponents_for_diffusion_creep.data['spharz'] = [diff['m'] for i in range(spharz_up)] + [diff_lm['m'] for i in range(spharz_low)]
    activation_energies_for_diffusion_creep.data['spharz'] = [diff['E'] for i in range(spharz_up)] + [diff_lm['E'] for i in range(spharz_low)]
    activation_volumes_for_diffusion_creep.data['spharz'] = [diff['V'] for i in range(spharz_up)] + [diff_lm['V'] for i in range(spharz_low)]
    # opcrust
    opcrust_up = 1
    opcrust_low = 0
    prefactors_for_diffusion_creep.data['opcrust'] = [diff['A'] for i in range(opcrust_up)] + [diff_lm['A'] for i in range(opcrust_low)]
    grain_size_exponents_for_diffusion_creep.data['opcrust'] = [diff['m'] for i in range(opcrust_up)] + [diff_lm['m'] for i in range(opcrust_low)]
    activation_energies_for_diffusion_creep.data['opcrust'] = [diff['E'] for i in range(opcrust_up)] + [diff_lm['E'] for i in range(opcrust_low)]
    activation_volumes_for_diffusion_creep.data['opcrust'] = [diff['V'] for i in range(opcrust_up)] + [diff_lm['V'] for i in range(opcrust_low)]
    # spharz
    opharz_up = 1
    opharz_low = 0
    prefactors_for_diffusion_creep.data['opharz'] = [diff['A'] for i in range(opharz_up)] + [diff_lm['A'] for i in range(opharz_low)]
    grain_size_exponents_for_diffusion_creep.data['opharz'] = [diff['m'] for i in range(opharz_up)] + [diff_lm['m'] for i in range(opharz_low)]
    activation_energies_for_diffusion_creep.data['opharz'] = [diff['E'] for i in range(opharz_up)] + [diff_lm['E'] for i in range(opharz_low)]
    activation_volumes_for_diffusion_creep.data['opharz'] = [diff['V'] for i in range(opharz_up)] + [diff_lm['V'] for i in range(opharz_low)]
    # dislocation creep
    prefactors_for_dislocation_creep = Parse.COMPOSITION()
    stress_exponents_for_dislocation_creep = Parse.COMPOSITION()
    activation_energies_for_dislocation_creep = Parse.COMPOSITION()
    activation_volumes_for_dislocation_creep = Parse.COMPOSITION()
    # background composition
    background_up = 4
    background_low = 4
    prefactors_for_dislocation_creep.data['background'] = [disl['A'] for i in range(background_up)] + [disl_lm['A'] for i in range(background_low)]
    stress_exponents_for_dislocation_creep.data['background'] = [disl['n'] for i in range(background_up)] + [disl_lm['n'] for i in range(background_low)]
    activation_energies_for_dislocation_creep.data['background'] = [disl['E'] for i in range(background_up)] + [disl_lm['E'] for i in range(background_low)]
    activation_volumes_for_dislocation_creep.data['background'] = [disl['V'] for i in range(background_up)] + [disl_lm['V'] for i in range(background_low)]
    # spcrust
    spcrust_sz = 1
    spcrust_up = 1
    spcrust_low = 2
    prefactors_for_dislocation_creep.data['spcrust'] = [0.5/spcrust_disl for i in range(spcrust_sz)]\
    + [disl['A'] for i in range(spcrust_up)] + [disl_lm['A'] for i in range(spcrust_low)]
    stress_exponents_for_dislocation_creep.data['spcrust'] = [1.0 for i in range(spcrust_sz)]\
    + [disl['n'] for i in range(spcrust_up)] + [disl_lm['n'] for i in range(spcrust_low)]
    activation_energies_for_dislocation_creep.data['spcrust'] = [0.0 for i in range(spcrust_sz)]\
    + [disl['E'] for i in range(spcrust_up)] + [disl_lm['E'] for i in range(spcrust_low)]
    activation_volumes_for_dislocation_creep.data['spcrust'] = [0.0 for i in range(spcrust_sz)]\
    + [disl['V'] for i in range(spcrust_up)] + [disl_lm['V'] for i in range(spcrust_low)]
    # spharz
    spharz_up = 4
    spharz_low = 4
    prefactors_for_dislocation_creep.data['spharz'] = [disl['A'] for i in range(spharz_up)] + [disl_lm['A'] for i in range(spharz_low)]
    stress_exponents_for_dislocation_creep.data['spharz'] = [disl['n'] for i in range(spharz_up)] + [disl_lm['n'] for i in range(spharz_low)]
    activation_energies_for_dislocation_creep.data['spharz'] = [disl['E'] for i in range(spharz_up)] + [disl_lm['E'] for i in range(spharz_low)]
    activation_volumes_for_dislocation_creep.data['spharz'] = [disl['V'] for i in range(spharz_up)] + [disl_lm['V'] for i in range(spharz_low)]
    # opcrust
    opcrust_up = 1
    opcrust_low = 0
    prefactors_for_dislocation_creep.data['opcrust'] = [disl['A'] for i in range(opcrust_up)] + [disl_lm['A'] for i in range(opcrust_low)]
    stress_exponents_for_dislocation_creep.data['opcrust'] = [disl['n'] for i in range(opcrust_up)] + [disl_lm['n'] for i in range(opcrust_low)]
    activation_energies_for_dislocation_creep.data['opcrust'] = [disl['E'] for i in range(opcrust_up)] + [disl_lm['E'] for i in range(opcrust_low)]
    activation_volumes_for_dislocation_creep.data['opcrust'] = [disl['V'] for i in range(opcrust_up)] + [disl_lm['V'] for i in range(opcrust_low)]
    # spharz
    opharz_up = 1
    opharz_low = 0
    prefactors_for_dislocation_creep.data['opharz'] = [disl['A'] for i in range(opharz_up)] + [disl_lm['A'] for i in range(opharz_low)]
    stress_exponents_for_dislocation_creep.data['opharz'] = [disl['n'] for i in range(opharz_up)] + [disl_lm['n'] for i in range(opharz_low)]
    activation_energies_for_dislocation_creep.data['opharz'] = [disl['E'] for i in range(opharz_up)] + [disl_lm['E'] for i in range(opharz_low)]
    activation_volumes_for_dislocation_creep.data['opharz'] = [disl['V'] for i in range(opharz_up)] + [disl_lm['V'] for i in range(opharz_low)]

    visco_plastic["Reference strain rate"] = str(strain_rate)
    visco_plastic["Prefactors for diffusion creep"] = prefactors_for_diffusion_creep.parse_back()
    visco_plastic["Grain size exponents for diffusion creep"] = grain_size_exponents_for_diffusion_creep.parse_back()
    visco_plastic["Activation energies for diffusion creep"] = activation_energies_for_diffusion_creep.parse_back()
    visco_plastic["Activation volumes for diffusion creep"] = activation_volumes_for_diffusion_creep.parse_back()
    visco_plastic["Prefactors for dislocation creep"] = prefactors_for_dislocation_creep.parse_back()
    visco_plastic["Stress exponents for dislocation creep"] = stress_exponents_for_dislocation_creep.parse_back()
    visco_plastic["Activation energies for dislocation creep"] = activation_energies_for_dislocation_creep.parse_back()
    visco_plastic["Activation volumes for dislocation creep"] = activation_volumes_for_dislocation_creep.parse_back()
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


def MaxRefinement(inputs, value):
    """
    options for max_refinement: change the maximum level of refinemment
    """
    isosurfaces_input = '10, %s, spcrust: 0.5 | 1.0; 9, %s, \
spharz: 0.5 | 1.0;  9, %s, opcrust: 0.5 | 1.0; 8, %s, \
opharz: 0.5 | 1.0; 9, %s, Temperature: 270.0 | 1173.0' % (value, value, value, value, value)
    inputs['Mesh refinement']['Initial adaptive refinement'] = str(value - 1 - int(inputs['Mesh refinement']['Initial global refinement']))
    inputs['Mesh refinement']['Isosurfaces']['Isosurfaces'] = isosurfaces_input
    return inputs, '_MRf%d' % value  # second entry is an appendix to case name
    pass


def DecouplingEclogiteTransiton(inputs, value):
    """
    options for decoupling_eclogite_transition: whether to decouple phase transiton from viscosity change
    """
    if value == 'CET':
        inputs['Material model']['Visco Plastic']['Decoupling eclogite viscosity'] = 'false'
    elif value == 'DET':
        inputs['Material model']['Visco Plastic']['Decoupling eclogite viscosity'] = 'true'
    return inputs, '_%s' % value  # second entry is an appendix to case name


def Restart(inputs, value):
    """
    calculate flow law parameters, when phase transition only happens on mantle composition
    """
    inputs['Resume computation'] = 'true'
    return inputs

        

def CreateNew(config, **kwargs):
    """        
        create cases under a directory
        read json file and prm file
        location of prm file is given by the json file
    """
    # read parameters
    paths = config['path']
    _root = paths['root']
    prm_path = paths['prm']
    extra_paths = paths['extra']
    # create case, using the interface defined in Cases.py.
    case_name = paths['base']
    newCase = CASE(case_name, prm_path)
    newCase.configure(RheologyCDPT, config)  # rheology
    for extra_path in extra_paths:
        newCase.add_extra_file(extra_path)  # add an extra file
    # hold, then only return
    hold = kwargs.get('hold', 0)
    if hold == 1:
        pass
    else:
        newCase.create(_root, fast_first_step=1)
    return newCase


def CreateNewWithOptions(config, options, **kwargs):
    """        
        create cases under a directory
        read json file and prm file
        location of prm file is given by the json file
    """
    # create case, using the interface defined in Cases.py.
    paths = config['path']
    _root = paths['root']
    key = options['key']
    values = options['values']
    newCases = []
    # I include options in order to generate multiple cases at a time
    for i in range(len(values)):
        newCase = CreateNew(config, hold=1)
        if key == 'max_refinement':
            newCase.configure(MaxRefinement, values[i], rename=1)
        if key == 'decoupling_eclogite_transiton':
            newCase.configure(DecouplingEclogiteTransiton, values[i], rename=1)
        newCases.append(newCase)
    # Create cases
    for newCase in newCases:
        newCase.create(_root, fast_first_step=1)


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
        Utilities.my_assert(type(number_particle_in_slab)==int, TypeError, "number_particle_in_slab must be an int value")
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
            x, y = Utilities.ggr2cart2(phi, r)
            # assign value in particle_data
            self.particle_data[i, 0] = x
            self.particle_data[i, 1] = y


def trans_point_in_wb(i_coord, std_coord, trans=0, **kwargs):
    '''
    todo
    Perform transform of coordinates of a point in world builder file.
    Inputs:
        i_coord (list of 2): coordinates of a point
        std_coord (list of 2): standard coordinates of a point, Part of this is copied into the result.
        trans (int): 0 - sph to cart; 1 - cart to sph
        kwargs (dict):
            Ro (float): outer radius
    '''
    Ro = kwargs.get('Ro', 6371e3)
    o_coord = []  # initiate
    # convert 1st dimention
    if trans == 0:
        o_coord.append(Ro*np.deg2rad(i_coord[0])) # length convertion, not coordinate transformation
    elif trans == 1:
        o_coord.append(np.rad2deg(i_coord[0]/Ro))
    else:
        raise ValueError("%s: wrong value for trans." % Utilities.func_name())
    o_coord.append(std_coord[1])  # copy 2nd dimention
    return o_coord


def trans_coord_in_wb(i_coords, std_coords, trans=0, **kwargs):
    '''
    todo
    Perform transform of coordinates in world builder file.
    This calls the trans_point_in_wb to transform each point
    This is intended to be called in the 'sph_cart_copy' funciton.
    The operation being done here is length convertion, not coordinate transformation
    Inputs:
        i_coords (list of float): input coordinate
        std_coords (list of float): 
            standard coordinate. Part of this is copied into the result.
            That is basically the 'redundant' demension.
        trans (int): 0 - sph to cart; 1 - cart to sph
        kwargs (dict):
            Ro (float): outer radius
    '''
    assert(type(i_coords) == list and type(std_coords) == list\
    and len(i_coords) == len(std_coords))
    o_coords = []  # initiate
    # convert first dimension, note one object contains multiple coordinate pairs
    Ro = kwargs.get('Ro', 6371e3)
    for i in range(len(i_coords)):
        i_coord = i_coords[i]
        std_coord = std_coords[i]
        o_coord = trans_point_in_wb(i_coord, std_coord, trans, Ro=Ro)
        o_coords.append(o_coord.copy())
    # copy second dimension
    # return
    return o_coords


def find_trans_coord_in_wb_feature(i_feature, std_feature, trans=0, **kwargs):
    '''
    todo
    Perform transform of coordinates in world builder features.
    This function would look for coordinates within features and perform the transformation
    Note that there are still coordinates in sections with these features that are unchanged,
    be sure to change those by hand.
    Inputs:
        i_feature (str): feature from input
        std_feature (std): a standard output feature
        trans (int): 0 - sph to cart; 1 - cart to sph
        kwargs (dict):
            Ro (float): outer radius
    '''
    Ro = kwargs.get('Ro', 6371e3)
    o_feature = i_feature.copy()
    coord_keys = ['coordinates', 'ridge coordinates']  # keys for components of coordinates
    point_keys = ['dip point']
    for key in coord_keys:
        if key in i_feature:
            i_coords = i_feature[key]
            std_coords = std_feature[key]
            o_coords = trans_coord_in_wb(i_coords, std_coords, trans, Ro=Ro)
            o_feature[key] = o_coords
    for key in point_keys:
        if key in i_feature:
            i_coord = i_feature[key]
            std_coord = std_feature[key]
            o_coord = trans_point_in_wb(i_coord, std_coord, trans, Ro=Ro)
            o_feature[key] = o_coord
    # recursively look for subfeatures
    for key, value in i_feature.items():
        if type(value) == dict:
            # in case of dict, call function recursively
            o_subfeature = find_trans_coord_in_wb_feature(value, std_feature[key], trans, Ro=Ro)
            o_feature[key] = o_subfeature
    return o_feature


def sph_cart_convert_case(case_dir, case_odir, trans=0):
    '''
    convert a case in sph geometry to cartesan and vice versa
    Inputs:
        case_dir (str): 
            directory to look for 'case.prm' and 'case.wb' as inputs, note that this directory must
            contain these two files.
        case_odir (str):
            directory to hold new case files. This will be made if not exist yet.
        trans(int):
            0 - sph to cart; 1 - cart to sph
    '''
    # todo
    assert(os.path.isdir(case_dir))
    if not os.path.isdir(case_odir):
        os.mkdir(case_odir)
    # call sph_cart_convert function
    prm_in_path = os.path.join(case_dir, 'case.prm')
    wb_in_path = os.path.join(case_dir, 'case.wb')
    outputs, outputs_wb = sph_cart_convert(prm_in_path, wb_in_path, trans)
    prm_o_path = os.path.join(case_odir, 'case.prm')
    wb_o_path = os.path.join(case_odir, 'case.wb')
    ParsePrm.WritePrmFile(prm_o_path, outputs)
    with open(wb_o_path, 'w') as fout:
        json.dump(outputs_wb, fout, indent=2)
    print("Note that there are still coordinates in sections with these features that are unchanged,\
    make sure to change those by hand.")


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
    parser.add_argument('-p', '--options', type=str,
                        default=None,
                        help='json file inputs for additional options')
    parser.add_argument('-o', '--outputs', type=str,
                        default='case_o.prm',
                        help='Some outputs(prm file)')
    parser.add_argument('-b', '--bool', type=int,
                        default=0,
                        help='a bool option')
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
    
    elif _commend == 'rheology_cdpt':
        # derive parameterization of lower mantle from a prm file and configurations
        # check file existence
        assert(os.access(arg.inputs, os.R_OK))
        assert(os.access(arg.json, os.R_OK))
        
        # read parameters and configuration
        with open(arg.inputs, 'r') as fin:
            _inputs = ParsePrm.ParseFromDealiiInput(fin)
        with open(arg.json, 'r') as fin:
            _config = json.load(fin)
        outputs = RheologyCDPT(_inputs, _config, lower_mantle_full=arg.bool)
        # screen output
        # file output
        with open(arg.outputs, 'w') as fout:
            ParsePrm.ParseToDealiiInput(fout, outputs)
    
    elif _commend == 'create':
        # create cases under a directory
        # read json file and prm file
        # location of prm file is given by the json file
        assert(os.access(arg.json, os.R_OK))
        with open(arg.json, 'r') as fin:
            _config = json.load(fin)
        CreateNew(_config)
        
    elif _commend == 'create_with_options':
        # create cases under a directory
        # read json file and prm file
        # location of prm file is given by the json file
        assert(os.access(arg.json, os.R_OK))
        with open(arg.json, 'r') as fin:
            _config = json.load(fin)
        # also read a 2nd json for options
        assert(os.access(arg.options, os.R_OK))
        with open(arg.options, 'r') as fin:
            options = json.load(fin)
        CreateNewWithOptions(_config, options)
    
    elif _commend == 'convert_sph_cart':
        #
        # todo
        sph_cart_convert_case(arg.inputs, arg.outputs, 1)
    elif (_commend in ['-h', '--help']):
        # example:
        Usage()

# run script
if __name__ == '__main__':
    main()