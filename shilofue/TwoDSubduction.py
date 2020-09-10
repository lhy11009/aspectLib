import os
import sys
import json
import argparse
import re
import numpy as np
from shilofue.Parse import COMPOSITION
from shilofue.Parse import ParseFromDealiiInput
from shilofue.Parse import ParseToDealiiInput
from shilofue.Parse import CASE, GROUP_CASE, PARSE_OPERATIONS,\
                           UpdateProjectMd, UpdateProjectJson, AutoMarkdownCase, AutoMarkdownGroup
from shilofue.Rheology import GetLowerMantleRheology
from shilofue.Utilities import my_assert, ggr2cart2
from shilofue.Doc import UpdateProjectDoc
from shilofue.Plot import ProjectPlot


class MY_PARSE_OPERATIONS(PARSE_OPERATIONS):
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
        PARSE_OPERATIONS.__init__(self)
        self.ALL_OPERATIONS += ["LowerMantle", "Particle"]

    def LowerMantle(self, Inputs, _config):
        """
        calculate flow law parameters
        """
        # parse from input
        jump = _config['upper_lower_viscosity']
        T = _config['T660']
        P = _config['P660']
        V1 = _config['LowerV']
        visco_plastic = Inputs["Material model"]['Visco Plastic']
        prefactors_for_diffusion_creep = COMPOSITION(visco_plastic["Prefactors for diffusion creep"])
        grain_size = float(visco_plastic["Grain size"])
        grain_size_exponents_for_diffusion_creep  = COMPOSITION(visco_plastic["Grain size exponents for diffusion creep"])
        activation_energies_for_diffusion_creep = COMPOSITION(visco_plastic["Activation energies for diffusion creep"])
        activation_volumes_for_diffusion_creep  = COMPOSITION(visco_plastic["Activation volumes for diffusion creep"])
        # call GetLowerMantleRheology to derive parameters for lower mantle flow law 
        backgroud_upper_mantle_diffusion = {}
        backgroud_upper_mantle_diffusion['A'] = prefactors_for_diffusion_creep.data['background'][0] 
        backgroud_upper_mantle_diffusion['d'] = grain_size
        backgroud_upper_mantle_diffusion['n'] = 1.0 
        backgroud_upper_mantle_diffusion['m'] = grain_size_exponents_for_diffusion_creep.data['background'][0] 
        backgroud_upper_mantle_diffusion['E'] = activation_energies_for_diffusion_creep.data['background'][0] 
        backgroud_upper_mantle_diffusion['V'] = activation_volumes_for_diffusion_creep.data['background'][0] 
        backgroud_lower_mantle_diffusion = GetLowerMantleRheology(backgroud_upper_mantle_diffusion, jump, T, P, V1=V1, strategy='A')
        # todo_future: add in choice of phases
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


class MYCASE(CASE):
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

    
def main():
    '''
    main function of this module
    Inputs:
        sys.arg[1](str):
            commend
        sys.arg[2, :](str):
            options
    '''
    # parse commend
    _available_commends = ['create', 'create_group', 'plot', 'update']  # only these commends are available now
    _commend = sys.argv[1]
    if _commend not in _available_commends:
        raise ValueError('Commend %s is not available.' % _commend)
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
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)
    # execute commend
    if _commend == 'create_group':
        print('Now we create a group of cases:')  # screen output
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
            _inputs = ParseFromDealiiInput(fin)
        if not os.path.isdir(arg.output_dir):
            os.mkdir(arg.output_dir)
        # create a directory under the name of the group
        _group_name = _config.get('name', 'foo')
        _odir = os.path.join(arg.output_dir, _group_name)
        my_assert(not os.path.isdir(_odir), ValueError, "The script doesn't support updating a pr-exiting group")
        os.mkdir(_odir)
        # initialte a class instance
        MyGroup = GROUP_CASE(MYCASE, _inputs, _config)
        # call __call__ function to generate
        _extra = _config.get('extra', {})
        # add an entry for parse_operations
        parse_operations = MY_PARSE_OPERATIONS()
        _case_names = MyGroup(parse_operations, _odir, extra=_extra, basename=_base_name)
        # generate auto.md
        AutoMarkdownCase(_group_name, _config, dirname=_odir)
        for _case_name in _case_names:
            _case_dir = os.path.join(_odir, _case_name)
            _case_json_file = os.path.join(_case_dir, 'config.json')
            with open(_case_json_file, 'r') as fin:
                _case_config = json.load(fin)
            AutoMarkdownCase(_case_name, _case_config, dirname=_case_dir)
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
            _inputs = ParseFromDealiiInput(fin)
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
        AutoMarkdownCase(_case_name, _config, dirname=_case_dir)
        print(_case_name)
    
    elif _commend == 'query':
        # for now, only out put the cases in this group
        print('Now we query into a group')
        _config_file = os.path.join(arg.output_dir, 'config.json')
        # check this group exist
        my_assert(os.path.isdir(arg.output_dir), FileExistsError, "%s doesn't exist" % arg.output_dir)
        my_assert(os.path.isdir(_config_file), FileExistsError, "%s doesn't exist" % arg._config_file)
        # initial class instance, todo
        # MyCase = MYCASE(_inputs, config=_config['config'], test=_config['test'])
        # call function to return case names
        # check that these cases exit

        pass

    elif _commend == 'update_doc':
        # todo_future
        pass

    elif _commend == 'update':
        # update a case
        _project_dir = arg.output_dir
        _project_dict = UpdateProjectJson(_project_dir)  # update project json file
        UpdateProjectMd(_project_dict, _project_dir)  # update auto.md file for every case
        ProjectPlot(_project_dict, _project_dir, 'png', update=False)  # plot figures for every case
        UpdateProjectDoc(_project_dict, _project_dir, images=['Statistics' ,'DepthAverage', 'PvMesh', 'visit'])

    elif _commend == 'plot':
        # todo_future
        # plot something
        pass


# run script
if __name__ == '__main__':
    main()