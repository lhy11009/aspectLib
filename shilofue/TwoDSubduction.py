import os
import sys
import json
import argparse
import re
import warnings
import pdb
import numpy as np
import shilofue.Parse as Parse
import shilofue.ParsePrm as ParsePrm
import shilofue.Doc as Doc
import shilofue.Plot as Plot
import shilofue.Rheology as Rheology
from numpy import linalg as LA
from matplotlib import pyplot as plt
from shilofue.Utilities import my_assert, ggr2cart2, cart2sph2, Make2dArray, UNITCONVERT


# global varibles
# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']

# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')

project = "TwoDSubduction"


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
            self.LowerMantle0(Inputs, _config)
        elif _type == 1:
            # when phase transition only happens on all compositions
            # There is a eclogite transition of crustal layer
            self.LowerMantle1(Inputs, _config)
        elif _type == 2:
            # setup lower mantle rheology for phase transition in a CDPT model
            # There is a eclogite transition of crustal layer
            self.LowerMantle2(Inputs, _config)

    def LowerMantle0(self, Inputs, _config):
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
    
    def LowerMantle1(self, Inputs, _config):
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
    
    def LowerMantle2(self, Inputs, _config):
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


class VISIT_OPTIONS(Parse.VISIT_OPTIONS):
    """
    inherite the VISIT_OPTIONS clase from Parse.py
    in order to add in additional settings
    """
    def Interpret(self, kwargs={}):
        """
        Interpret the inputs, to be reloaded in children
        """
        # call function from parent
        Parse.VISIT_OPTIONS.Interpret(self, kwargs)

        # default settings
        self.odict['IF_PLOT_SLAB'] = 'False'
        self.odict['IF_EXPORT_SLAB_MORPH'] = 'False'
        particle_output_dir = os.path.join(self._output_dir, "slab_morphs")
        self.odict["PARTICLE_OUTPUT_DIR"] = particle_output_dir

        # optional settings
        for key, value in kwargs.items():
            # slab
            if key == 'slab':
                self.odict['IF_PLOT_SLAB'] = 'True'
                self.odict['PLOT_SLAB_STEPS'] = value.get('steps', [0])
                self.odict['IF_DEFORM_MECHANISM'] = value.get('deform_mechanism', 0)
            # export particles for slab morph
            elif key == 'slab_morph':
                self.odict['IF_EXPORT_SLAB_MORPH'] = 'True'
                # check directory
                if not os.path.isdir(particle_output_dir):
                    os.mkdir(particle_output_dir)
        


class VISIT_XYZ(Parse.VISIT_XYZ):
    """
    Read .xyz file exported from visit and do analysis
    Attributes:
        max_depth: max_depth of slab

    """
    def __init__(self):
        """
        initiation
        """
        # call parental function
        Parse.VISIT_XYZ.__init__(self)

    def Analyze(self, kwargs):
        """
        analyze data
        Args:
            kwargs(dict): options
                radius(float): radius of the earth
                depth_ranges(float): ranges of depth to compute dip angle
        """
        # get options
        col_x = self.column_indexes['x']
        col_y = self.column_indexes['y']
        col_id = self.column_indexes['id']
        radius = kwargs.get("radius", 6371e3)

        # sort by id
        self.data = self.data[self.data[:, col_id].argsort()]

        # transfer to sph
        x = self.data[:, col_x]
        y = self.data[:, col_y]
        r, ph = cart2sph2(x, y)
        depth = radius - r

        # get maximum depth
        max_depth = np.max(depth)

        # get trench position
        # a depth for trench, in m
        trench_depth = kwargs.get('trench_depth', 1e3)
        mask_slab = (depth > trench_depth)
        trench_position = ph[mask_slab][0]

        # get length of slab
        # both length and dip has n-1 component because they are computed on the intervals between 2 points
        length = ((r[0: -1] - r[1:])**2.0 + r[0: -1]**2.0*(ph[0: -1] - ph[1:])**2.0)**0.5
        slab_length = LA.norm(length, 1)

        # get dip angle
        dip = SlabDip(r[0: -1], ph[0: -1], r[1:], ph[1:])

        # get average curvature in depth range
        depth_ranges = kwargs.get('depth_ranges', [[0.0, 6371e3]])
        my_assert(type(depth_ranges) is list, TypeError, "VISIT_XYZ.Analyze: depth_ranges must be a list")
        dips_in_ranges = np.zeros(len(depth_ranges))
        limit = 1e-6
        for i in range(len(depth_ranges)):
            depth_range = depth_ranges[i]
            mask_range = (dip > depth_range[0]) * (dip < depth_range[1])
            total = dip[mask_range].dot(length[mask_range])
            weight = LA.norm(length[mask_range], 1)
            if weight < limit:
                # i.e. max_depth < depth_range[1]
                dips_in_ranges[i] = 0.0
            else:
                dips_in_ranges[i] = total / weight
    
        # construct header
        # append time if present
        try:
            _time = kwargs['time']
        except KeyError:
            self.output_header = {
             'Maximum depth': {'col': 0, 'unit': self.header['x']['unit']},
             'Trench position': {'col': 1, 'unit': 'rad'},
             'Slab length': {'col': 2, 'unit': self.header['x']['unit']}
            }
            total_cols = 3
        else:
            _time_unit = kwargs.get('time_unit', 'yr')
            self.output_header = {
             'Time': {'col': 0, 'unit': _time_unit},
             'Maximum depth': {'col': 1, 'unit': self.header['x']['unit']},
             'Trench position': {'col': 2, 'unit': 'rad'},
             'Slab length': {'col': 3, 'unit': self.header['x']['unit']}
            }
            total_cols = 4
        for i in range(len(depth_ranges)):
            key = 'Dip angle %d_%d (rad)' % (int(depth_ranges[i][0]), int(depth_ranges[i][1]))
            self.output_header[key] = {'col': i+total_cols}

        # manage output
        # append time if present
        try:
            _time = kwargs['time']
        except KeyError:
            output_data_temp = [max_depth, trench_position, slab_length]
        else:
            output_data_temp = [_time, max_depth, trench_position, slab_length]
        for i in range(len(depth_ranges)):
            output_data_temp.append(dips_in_ranges[i])
        self.output_data = Make2dArray(output_data_temp)


class SLAB_MORPH_PLOT(Plot.LINEARPLOT):
    '''
    Class for plotting depth average file.
    This is an inheritage of the LINEARPLOT class

    Attributes:
    Args:
    '''
    def __init__(self, _name, **kwargs):
        Plot.LINEARPLOT.__init__(self, _name, kwargs)  # call init from base function
    
    def ManageData(self):
        '''
        manage data, get new data for this class
        for the base class, this method simply takes the combination
        of self.data
        Returns:
            _data_list(list):
                list of data for ploting
        '''
        _data_list = Plot.LINEARPLOT.ManageData(self)
        # add subduction rate
        _col_depth = self.header['Maximum_depth']['col']
        _unit_depth = self.header['Maximum_depth']['unit']
        _col_time = self.header['Time']['col']
        _unit_time = self.header['Time']['unit']
        _depths = self.data[:, _col_depth]
        _times = self.data[:, _col_time]
        # get derivative
        _size = _depths.size
        _depths_i = np.zeros(_size + 1)
        _depths_i[0: _size] = _depths
        _depths_i[1: _size+1] += _depths 
        _times_i = np.zeros(_size + 1)
        _times_i[0: _size] = _times
        _times_i[1: _size+1] += _times 
        _rate = np.diff(_depths_i) / np.diff(_times_i)
        _data_list.append(_rate)
        self.header['Subduction_rate'] = {}
        self.header['Subduction_rate']['col'] = self.header['total_col']
        self.header['Subduction_rate']['unit'] = '%s/%s' % (_unit_depth, _unit_time)
        self.header['total_col'] += 1
        return _data_list


def SlabDip(r0, ph0, r1, ph1):
    """
    compute the dip angle between 2 adjacent point
    """
    alpha = np.arctan2((r0 - r1), (r1 * (ph1 - ph0)))
    return alpha


def SlabMorph(case_dir, kwargs={}):
    """
    Slab morphology
    Inputs:
        case_dir(str): directory of case
        kwargs(dict): options
    """
    case_output_dir = os.path.join(case_dir, 'output')
    case_morph_dir = os.path.join(case_output_dir, 'slab_morphs')

    # Initiation
    Visit_Xyz = VISIT_XYZ()
    
    # a header for interpreting file format
    # note that 'col' starts form 0
    header = {
        'x': {'col': 1, 'unit': 'm'},
        'y': {'col': 2, 'unit': 'm' },
        'id': {'col': 4}
    }

    # depth range
    # this is for computing dip angles with different ranges
    depth_ranges = kwargs.get('depth_ranges', [[0, 100e3], [100e3, 400e3], [400e3, 6371e3]])
    my_assert(type(depth_ranges) == list, TypeError, "depth_ranges mush be a list")

    # remove older results 
    ofile = os.path.join(case_output_dir, 'slab_morph') 
    if os.path.isfile(ofile):
        os.remove(ofile)
    
    #   loop for every snap and call function
    snaps, times, _= Parse.GetSnapsSteps(case_dir, 'particle')

    for i in snaps:
        visit_xyz_file = os.path.join(case_morph_dir, 'visit_particles_%06d.xyz' % i)
        Visit_Xyz(visit_xyz_file, header=header, ofile=ofile, depth_ranges=depth_ranges, time=times[i])


def ProjectPlot(case_dirs, _file_type, **kwargs):
    '''
    Plot figures for all cases in this project
    Inputs:
        kwargs:
            update(True or False): if True, update existing figures
    '''
    # future: not used, may remove
    update = kwargs.get('update', False)
    pdict = kwargs.get('pdict', {})
    # Init the UnitConvert class
    UnitConvert = UNITCONVERT()

    # plot statistics ouput
    plot_options = pdict.get('slab_morph', {})
    Slab_morph_plot = SLAB_MORPH_PLOT('slab_morph', unit_convert=UnitConvert, options=plot_options)

    # loop for cases and post process
    for case_dir in case_dirs:
        # first generate slab_morph output
        img_dir = os.path.join(case_dir, 'img')
        if not os.path.isdir(img_dir):
            os.mkdir(img_dir)
        
        # conditional plot for slab morph
        # future extra options
        # with open(arg.json_file, 'r') as fin:
        #     dict_in = json.load(fin)
        #     extra_options = dict_in.get('slab_morph', {})
        particle_file = os.path.join(case_dir, 'output', 'particles.visit')
        ofile = os.path.join(img_dir, 'slab_morph.png')
        is_plot = False
        if os.path.isfile(particle_file):
            if not os.path.isfile(ofile) or \
               os.stat(particle_file)[8] > os.stat(ofile)[8]:
                is_plot = True
        if is_plot:
            # process slab morph with extra options
            extra_options = {}
            try:
                SlabMorph(case_dir, extra_options)
            except FileNotFoundError:
                warnings.warn('process_slab_morph: file existence requirements are not met')
            # then plot the slab morph figure
            filein = os.path.join(case_dir, 'output', 'slab_morph')
            # Get options
            # plot
            if os.path.isfile(filein):
                Slab_morph_plot(filein, fileout=ofile)


def PlotTestResults(source_dir, **kwargs):
        # initialize
        _case_img_dir = kwargs.get('output_dir', './test_results')
        if not os.path.isdir(_case_img_dir):
            os.mkdir(_case_img_dir)
        _file_type = kwargs.get('type', 'pdf')
        # convert unit 
        UnitConvert = UNITCONVERT()

        # case TwoDSubduction_pyrolite_density_1_0: density depth-average plot
        _case_output_dir = os.path.join(source_dir, 'output-TwoDSubduction_pyrolite_density_1_0')
        assert(os.path.isdir(_case_output_dir))
        # read data
        DepthAverage = Plot.DEPTH_AVERAGE_PLOT('DepthAverage', unit_convert=UnitConvert)
        _depth_average_file = os.path.join(_case_output_dir, 'depth_average.txt')
        assert(os.access(_depth_average_file, os.R_OK))
        data_list = DepthAverage.ReadDataStep(_depth_average_file, datatype=["depth", "adiabatic_density", "viscosity"])
        # plot
        fig, axs = plt.subplots(1, 2, figsize=(10, 5))
        axs[0].plot(data_list[0]/1e3, data_list[1], '-b', label='density')
        axs[0].grid()
        axs[0].set_ylabel('Density [kg/m3]')
        axs[0].set_xlabel('Depth [km]')
        axs[0].set_ylim((3100, 4800))
        axs[0].legend()
        # data from Chust 17
        Chust17_data = np.array([
          [0.37078107742533817, 3.1425201813692216], 
          [0.950936590564412, 3.2249689538700497],
          [4.289074088639609, 3.3645200063366243],
          [13.733855181574159, 3.6325618236181016],
          [13.913762104113376, 3.7473406677161245],
          [15.29053749408711, 3.7902418219342437],
          [18.542044934634184, 3.850878430998938],
          [19.224893623385242, 3.908204043696113],
          [22.7711869758023, 3.9795712587613234],
          [23.26083512704819, 4.008218112771951],
          [23.325863060256637, 4.227045194967481],
          [29.223242999526967, 4.430894804301371],
          [54.81555663135366, 4.7545978736862855]
        ])
        axs[1].plot(Chust17_data[:, 0], Chust17_data[:,1]*1000, '-k', label='Chust17_data')
        axs[1].set_xlabel('Pressure [GPa]')
        axs[1].set_ylim((3100, 4800))
        axs[1].grid()
        fig.tight_layout()
        # save
        _time = 0.0
        _ofile_route = os.path.join(_case_img_dir, 'Pyrolite_density_1_0.%s' % _file_type)
        plt.savefig(_ofile_route)
        print('Plot has been generated: ', _ofile_route)  # screen output


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
    parser = argparse.ArgumentParser(description='TwoDSubdunction Project')
    parser.add_argument('-b', '--base_file', type=str,
                        default='./files/TwoDSubduction/base.prm',
                        help='Filename for base file')
    parser.add_argument('-U', '--use_basename_as_base_file', type=int,
                        default=1,
                        help='Whether we use basename as base file')
    parser.add_argument('-j', '--json_file', type=str,
                        default='./config_case.json',
                        help='Filename for json file')
    parser.add_argument('-o', '--output_dir', type=str,
                        default='../TwoDSubduction/',
                        help='Directory for output')
    parser.add_argument('-e', '--operations_file', type=str,
                        default=None,
                        help='A file that has a list of operations, if not given, do all the available operations')
    parser.add_argument('-i', '--input_dir', type=str,
                        default=shilofue_DIR,
                        help='A directory that contains the input')
    parser.add_argument('-s', '--step', type=int,
                        default=0,
                        help='timestep')
    parser.add_argument('-ex', '--extension', type=str,
                        default='png',
                        help='extension for output')
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    # execute commend
    if _commend == 'create_group':
        # create a group
        # example usage:
        #    python -m shilofue.TwoDSubduction create_group -j config_group.json 2>&1 > .temp
        # create a group of cases
        # read files
        # read configuration
        with open(arg.json_file, 'r') as fin:
            _config = json.load(fin)
        _base_name = _config.get('basename', '')
        # read base file
        if arg.use_basename_as_base_file == 1:
            _filename = './files/TwoDSubduction/%s.prm' % _base_name
        else:
            _filename = arg.base_file
        with open(_filename, 'r') as fin:
            _inputs = ParsePrm.ParseFromDealiiInput(fin)
        if not os.path.isdir(arg.output_dir):
            os.mkdir(arg.output_dir)
            print('Now we create a group of cases:')  # screen output
        else:
            print('Now we update a group of cases:')  # screen output

        # create a directory under the name of the group
        _group_name = _config.get('name', 'foo')
        _odir = os.path.join(arg.output_dir, _group_name)
        # By default, we don't update
        update_ = _config.get('update', 0)
        if not update_:
            my_assert(os.path.isdir(_odir) is False, ValueError, 'Going to update a pr-exiting group, but update is not included in the option')
        if not os.path.isdir(_odir):
            os.mkdir(_odir)
        
        # initialte a class instance
        MyGroup = Parse.GROUP_CASE(MYCASE, _inputs, _config)
        # call __call__ function to generate
        _extra = _config.get('extra', {})
        # add an entry for parse_operations
        parse_operations = MY_PARSE_OPERATIONS()
        _case_names = MyGroup(parse_operations, _odir, extra=_extra, basename=_base_name, update=update_)
        # generate auto.md
        # check if there is alread a preexisting group
        Parse.AutoMarkdownGroup(_group_name, _config, dirname=_odir)
        for _case_name in _case_names:
            _case_dir = os.path.join(_odir, _case_name)
            _case_json_file = os.path.join(_case_dir, 'config.json')
            with open(_case_json_file, 'r') as fin:
                _case_config = json.load(fin)
            Parse.AutoMarkdownCase(_case_name, _case_config, dirname=_case_dir)
        print(_group_name)
        for _case_name in _case_names:
            # ouptut to screen
            print(_case_name)

    elif _commend == 'create':
        print('Now we create a single case:')  # screen output
        # create a case
        # read files
        # read configuration
        with open(arg.json_file, 'r') as fin:
            _config = json.load(fin)
        _base_name = _config.get('basename', '')
        # read base file
        if arg.use_basename_as_base_file == 1:
            _filename = './files/TwoDSubduction/%s.prm' % _base_name
        else:
            _filename = arg.base_file
        with open(_filename, 'r') as fin:
            _inputs = ParsePrm.ParseFromDealiiInput(fin)
        if not os.path.isdir(arg.output_dir):
            os.mkdir(arg.output_dir)
        # Initial a case
        MyCase = MYCASE(_inputs, config=_config['config'], test=_config['test'])
        # call __call__ function to generate
        _extra = _config.get('extra', {})
        # also add extra files
        _extra_files = _config.get('extra_file', {})
        # add an entry for parse_operations
        parse_operations = MY_PARSE_OPERATIONS()
        _case_name = MyCase(parse_operations, dirname=arg.output_dir, extra=_config['extra'], basename=_base_name, extra_file=_extra_files)
        # generate markdown file
        _case_dir = os.path.join(arg.output_dir, _case_name)
        Parse.AutoMarkdownCase(_case_name, _config, dirname=_case_dir)
        print(_case_name)
        # check this group exist
        my_assert(os.path.isdir(arg.output_dir), FileExistsError, "%s doesn't exist" % arg.output_dir)
        # initial class instance, future
        # MyCase = MYCASE(_inputs, config=_config['config'], test=_config['test'])
        # call function to return case names
        # check that these cases exit

        pass

    elif _commend == 'update_docs':
        # update the contents of the mkdocs
        # example usage:
        #   python -m shilofue.TwoDSubduction update_docs -o /home/lochy/ASPECT_PROJECT/TwoDSubduction -j post_process.json
        _project_dir = arg.output_dir
        _project_dict = Parse.UpdateProjectJson(_project_dir)  # update project json file
        
        # load options for post_process
        # load the project level configuration as default
        project_pp_json = os.path.join(ASPECT_LAB_DIR, 'files', project, 'post_process.json')
        with open(project_pp_json, 'r') as fin:
            pdict = json.load(fin)
        # load explicitly defined parameters
        with open(arg.json_file, 'r') as fin:
            pdict1 = json.load(fin)
        pdict.update(pdict1)
        
        # append analysis
        analysis_file = os.path.join(ASPECT_LAB_DIR, 'analysis.json')
        if os.path.isfile(analysis_file):
            with open(analysis_file, 'r') as fin:
                analysis_dict = json.load(fin)
        else:
            analysis_dict = {}
        
        # update docs
        docs_dict = pdict.get('docs', {})
        imgs = docs_dict.get('imgs', [])
        Doc.UpdateProjectDoc(_project_dict, _project_dir, images=imgs, analysis=analysis_dict)

    elif _commend == 'update':
        # update a case
        # example usage:
        #   python -m shilofue.TwoDSubduction update -o /home/lochy/ASPECT_PROJECT/TwoDSubduction -j post_process.json
        _project_dir = arg.output_dir
        _project_dict = Parse.UpdateProjectJson(_project_dir)  # update project json file
        
        # load options for post_process
        # load the project level configuration as default
        project_pp_json = os.path.join(ASPECT_LAB_DIR, 'files', project, 'post_process.json')
        with open(project_pp_json, 'r') as fin:
            pdict = json.load(fin)
        # load explicitly defined parameters
        with open(arg.json_file, 'r') as fin:
            pdict1 = json.load(fin)
        # update every entry in pdict1
        for key, value in pdict1.items():
            if type(value) == dict:
                try:
                    _ = pdict[key]
                    pdict[key].update(value)
                except KeyError:
                    pdict[key] = value
            else:
                pdict[key] = value

        # update auto.md file for every case
        Parse.UpdateProjectMd(_project_dict, _project_dir)

        # plot figures for every case
        # get sub cases
        pp_source_dirs = pdict.get('dirs', [])
        _format = pdict.get('py_format', 'png')
        for pp_source_dir_base in pp_source_dirs:
            pp_source_dir = os.path.join(_project_dir, pp_source_dir_base)
            pp_case_dirs = Parse.GetSubCases(pp_source_dir)
            Plot.ProjectPlot(pp_case_dirs, _format, update=False, pdict=pdict)
            # deal with project defined plots
            ProjectPlot(pp_case_dirs, _format, update=False, pdict=pdict)

    elif _commend == 'plot_newton_solver_step':
        # Plot one step from Newton solver
        # use -i option as input and -o option as output dir
        # example usage:
        #   python -m shilofue.TwoDSubduction plot_newton_solver_step -i tests/integration/fixtures/test-plot/newton_solver -o .test -s 1 --extension pdf
        filein = arg.input_dir
        output_dir = arg.output_dir
        step = arg.step
        extension = arg.extension
        ofile_route = os.path.join(output_dir, 'NewtonSolverStep.%s' % extension)
        # plot newton solver output
        NewtonSolverStep = Plot.NEWTON_SOLVER_PLOT('NewtonSolverStep')
        # plot step0
        NewtonSolverStep.GetStep(step)
        NewtonSolverStep(filein, fileout=ofile_route)
        pass
    
    elif _commend == 'plot_newton_solver':
        # plot the whole history outputs from Newton solver
        # use -i option as input and -o option as output dir
        # example usages:
        #   python -m shilofue.TwoDSubduction plot_newton_solver -i tests/integration/fixtures/test-plot/newton_solver -o .test
        filein = arg.input_dir
        output_dir = arg.output_dir
        step = arg.step
        ofile_route = os.path.join(output_dir, 'NewtonSolver.pdf')
        # plot newton solver output
        NewtonSolverStep = Plot.NEWTON_SOLVER_PLOT('NewtonSolver')
        # plot step0
        NewtonSolverStep(filein, fileout=ofile_route)
        pass
    
    elif _commend == 'plot_machine_time':
        # plot the machine time output
        # use -i option as input and -o option as output dir
        # example usages:
        #   python -m shilofue.TwoDSubduction plot_machine_time -i tests/integration/fixtures/test-plot/machine_time -o .test
        filein = arg.input_dir
        output_dir = arg.output_dir
        ofile = os.path.join(output_dir, 'MachineTime.pdf')
        # plot newton solver output
        MachineTime = Plot.MACHINE_TIME_PLOT('MachineTime')
        # plot step0
        MachineTime(filein, fileout=ofile)
        pass
    
    elif _commend == 'plot_slab_morph':
        # plot the slab morph output
        # use -i option as input and -o option as output dir
        # example usages:
        #   python -m shilofue.TwoDSubduction plot_slab_morph 
        #       -i /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear26/cr80w5ULV3.000e+01/output/slab_morph 
        #       -o /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear26/cr80w5ULV3.000e+01/img
        filein = arg.input_dir
        output_dir = arg.output_dir
        ofile = os.path.join(output_dir, 'slab_morph.png')
        # Init the UnitConvert class
        UnitConvert = UNITCONVERT()
        # Get options
        project_pp_json = os.path.join(ASPECT_LAB_DIR, 'files', 'TwoDSubduction', 'post_process.json')
        with open(project_pp_json, 'r') as fin:
            pdict = json.load(fin)
        plot_options = pdict.get('slab_morph', {})
        Slab_morph_plot = SLAB_MORPH_PLOT('slab_morph', unit_convert=UnitConvert, options=plot_options)
        # plot
        Slab_morph_plot(filein, fileout=ofile)

    elif _commend == 'process_slab_morph':
        # process slab morphology from visit particle output
        # generate a file that can be used for plot
        # example usages:
        # python -m shilofue.TwoDSubduction process_slab_morph -i 
        #   /home/lochy/ASPECT_PROJECT/TwoDSubduction/non_linear26/cr80w5ULV3.000e+01 -j post_process.json
        case_dir = arg.input_dir
        # process slab morph with extra options
        with open(arg.json_file, 'r') as fin:
            dict_in = json.load(fin)
            extra_options = dict_in.get('slab_morph', {})
        try:
            SlabMorph(case_dir, extra_options)
        except FileNotFoundError:
            warnings.warn('process_slab_morph: file existence requirements are not met')
    
            
    elif _commend == 'plot_slab_morph_case':
        # plot the slab morph output for a case
        # first generate slab_morph output
        case_dir = arg.input_dir
        # process slab morph with extra options
        with open(arg.json_file, 'r') as fin:
            dict_in = json.load(fin)
            extra_options = dict_in.get('slab_morph', {})
        try:
            SlabMorph(case_dir, extra_options)
        except FileNotFoundError:
            warnings.warn('process_slab_morph: file existence requirements are not met')
        # then plot the slab morph figure
        filein = os.path.join(case_dir, 'output', 'slab_morph')
        output_dir = os.path.join(case_dir, 'img')
        ofile = os.path.join(output_dir, 'slab_morph.png')
        # Init the UnitConvert class
        UnitConvert = UNITCONVERT()
        # Get options
        project_pp_json = os.path.join(ASPECT_LAB_DIR, 'files', 'TwoDSubduction', 'post_process.json')
        with open(project_pp_json, 'r') as fin:
            pdict = json.load(fin)
        plot_options = pdict.get('slab_morph', {})
        Slab_morph_plot = SLAB_MORPH_PLOT('slab_morph', unit_convert=UnitConvert, options=plot_options)
        # plot
        Slab_morph_plot(filein, fileout=ofile)


    elif _commend == 'plot':
        # future
        # plot something
        pass

    elif _commend == 'visit_options':
        # output bash options to a file that could be
        # read by bash script
        # initiate class object
        case_dir = arg.input_dir

        Visit_Options = VISIT_OPTIONS(case_dir)

        # load extra options
        if arg.json_file == './config_case.json':
            # no json file is giving
            extra_options = {}
        else:
            with open(arg.json_file, 'r') as fin:
                dict_in = json.load(fin)
                extra_options = dict_in.get('visit', {})

        # call function
        ofile = os.path.join(ASPECT_LAB_DIR, 'visit_keys_values')
        Visit_Options(ofile, extra_options)
        pass
    
    elif _commend == 'plot_test_results':
        # plot the result of tests
        # example:
        # python -m shilofue.TwoDSubduction plot_test_results -i 
        #  /home/lochy/softwares/aspect/build_TwoDSubduction/tests/ -o $TwoDSubduction_DIR/test_results
        source_dir = arg.input_dir
        # todo
        PlotTestResults(source_dir, output_dir=arg.output_dir)
    
    else:
        raise ValueError('Commend %s is not available.' % _commend)


# run script
if __name__ == '__main__':
    main()