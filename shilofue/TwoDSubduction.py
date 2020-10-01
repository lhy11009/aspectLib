import os
import sys
import json
import argparse
import re
import numpy as np
import shilofue.Parse as Parse
import shilofue.Doc as Doc
import shilofue.Plot as Plot
import shilofue.Rheology as Rheology
from numpy import linalg as LA
from shilofue.Utilities import my_assert, ggr2cart2, cart2sph2, Make2dArray


# global varibles
# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']

# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')


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


class BASH_OPTIONS(Parse.BASH_OPTIONS):
    """
    parse .prm file to a option file that bash can easily read
    """
    def Interpret(self):
        """
        Interpret the inputs, to be reloaded in children
        """
        # call function from parent
        Parse.BASH_OPTIONS.Interpret(self)

        # own implementations
        # initial adaptive refinement
        self.odict['INITIAL_ADAPTIVE_REFINEMENT'] = self.idict['Mesh refinement'].get('Initial adaptive refinement', '6')


class VISIT_XYZ(Parse.VISIT_XYZ):
    """
    todo
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
        self.output_header = {
            'Maximum depth': {'col': 0, 'unit': self.header['x']['unit']},
            'Trench position': {'col': 1},
            'Slab length': {'col': 2, 'unit': self.header['x']['unit']}
        }
        total_cols = 2
        for i in range(len(depth_ranges)):
            key = 'Dip angle %d_%d' % (int(depth_ranges[i][0]), int(depth_ranges[i][1]))
            self.output_header[key] = {'col': i+total_cols}

        # manage output
        output_data_temp = [max_depth, trench_position, slab_length]
        for i in range(len(depth_ranges)):
            output_data_temp.append(dips_in_ranges[i])
        self.output_data = Make2dArray(output_data_temp)


def SlabDip(r0, ph0, r1, ph1):
    """
    compute the dip angle between 2 adjacent point
    """
    alpha = np.arctan2((r0 - r1), (r1 * (ph1 - ph0)))
    return alpha

    
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
            _inputs = Parse.ParseFromDealiiInput(fin)
        if not os.path.isdir(arg.output_dir):
            os.mkdir(arg.output_dir)
        # create a directory under the name of the group
        _group_name = _config.get('name', 'foo')
        _odir = os.path.join(arg.output_dir, _group_name)
        my_assert(not os.path.isdir(_odir), ValueError, "The script doesn't support updating a pr-exiting group")
        os.mkdir(_odir)
        # initialte a class instance
        MyGroup = Parse.GROUP_CASE(MYCASE, _inputs, _config)
        # call __call__ function to generate
        _extra = _config.get('extra', {})
        # add an entry for parse_operations
        parse_operations = MY_PARSE_OPERATIONS()
        _case_names = MyGroup(parse_operations, _odir, extra=_extra, basename=_base_name)
        # generate auto.md
        Parse.AutoMarkdownCase(_group_name, _config, dirname=_odir)
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
            _inputs = Parse.ParseFromDealiiInput(fin)
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
    
    elif _commend == 'query':
        # for now, only out put the cases in this group
        print('Now we query into a group')
        _config_file = os.path.join(arg.output_dir, 'config.json')
        # check this group exist
        my_assert(os.path.isdir(arg.output_dir), FileExistsError, "%s doesn't exist" % arg.output_dir)
        my_assert(os.path.isdir(_config_file), FileExistsError, "%s doesn't exist" % arg._config_file)
        # initial class instance, future
        # MyCase = MYCASE(_inputs, config=_config['config'], test=_config['test'])
        # call function to return case names
        # check that these cases exit

        pass

    elif _commend == 'update_doc':
        # future
        pass

    elif _commend == 'update':
        # update a case
        # example usage:
        #   python -m shilofue.TwoDSubduction update -o /home/lochy/ASPECT_PROJECT/TwoDSubduction
        _project_dir = arg.output_dir
        _project_dict = Parse.UpdateProjectJson(_project_dir)  # update project json file
        Parse.UpdateProjectMd(_project_dict, _project_dir)  # update auto.md file for every case
        Plot.ProjectPlot(_project_dict, _project_dir, 'png', update=False)  # plot figures for every case
        Doc.UpdateProjectDoc(_project_dict, _project_dir, images=['Statistics' ,'DepthAverage', 'NewtonSolver', 'PvMesh', 'visit'])

    elif _commend == 'plot_newton_solver_step':
        # Plot one step from Newton solver
        # use -i option as input and -o option as output dir
        # example usage:
        #   python -m shilofue.TwoDSubduction plot_newton_solver_step -i tests/integration/fixtures/test-plot/newton_solver -o .test -s 1  
        filein = arg.input_dir
        output_dir = arg.output_dir
        step = arg.step
        ofile_route = os.path.join(output_dir, 'NewtonSolverStep.pdf')
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

    elif _commend == 'plot':
        # future
        # plot something
        pass

    elif _commend == 'bash_options':
        # output bash options to a file that could be
        # read by bash script
        # initiate class object
        case_dir = arg.input_dir
        Base_Options = BASH_OPTIONS(case_dir)

        # call function
        ofile = os.path.join(ASPECT_LAB_DIR, 'visit_keys_values')
        Base_Options(ofile)
        pass
    
    else:
        raise ValueError('Commend %s is not available.' % _commend)


# run script
if __name__ == '__main__':
    main()