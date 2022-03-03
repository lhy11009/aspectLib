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
from shilofue.PlotDepthAverage import DEPTH_AVERAGE_PLOT
import shilofue.Rheology as Rheology
from numpy import linalg as LA
from matplotlib import pyplot as plt

from shilofue.Schedule import ScheduleTask


# global varibles
# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
PROJECT_DIR = os.environ['TwoDSubduction_DIR']
# import utilities in subdirectiory
sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
from Utilities import my_assert, ggr2cart2, cart2sph2, Make2dArray, UNITCONVERT
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')

project = "TwoDSubduction"

def Usage():
    print("\
(One liner description\n\
\n\
Examples of usage: \n\
\n\
  - generate scheduler files for a single case: \n\
\n\
        Lib_TwoDSubduction schedule_case_pp -i EBA_CDPT1/eba_cdpt_SA80.0_OA40.0 \
\n\
  - generate scheduler files for a group of cases: \n\
\n\
        Lib_TwoDSubduction schedule_group_pp -i EBA_CDPT1 \
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
    DepthAverage = DEPTH_AVERAGE_PLOT('DepthAverage', unit_convert=UnitConvert)
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


def SchedulePPCase(case):
    '''
    Schedule post-procession tasks for a single case
    Inputs:
        case(str): case path relative to $TwoDSuduction_DIR
    '''
    case_dir = os.path.join(PROJECT_DIR, case)
    assert(os.path.isdir(case_dir))
    local_tasks = []
    local_tasks.append("Lib_rsync_case peloton TwoDSubduction %s" % case)  # sync data with server
    local_tasks.append("Lib_TwoDSubduction0_PlotCase plot_case -i %s" % case_dir)  # plot case results
    local_tasks.append("Lib_TwoDSubduction0_PlotCase morph_case -i %s" % case_dir)  # plot slab morphology
    ScheduleTask(case_dir, local_tasks) # add tasks to scheduler and generate scripts


def SchedulePPGroup(group):
    '''
    Schedule post-procession tasks for a group of cases
    Inputs:
        group(str): group path relative to $TwoDSuduction_DIR
    '''
    group_dir = os.path.join(PROJECT_DIR, group)
    assert(os.path.isdir(group_dir))
    local_tasks = []
    local_tasks.append("Lib_rsync_case peloton TwoDSubduction %s" % group)  # sync data with server
    local_tasks.append("Lib_TwoDSubduction0_PlotCase plot_case_in_dir -i %s" % group_dir)  # plot case results
    local_tasks.append("Lib_TwoDSubduction0_PlotCase morph_case_in_dir -i %s" % group_dir)  # plot slab morphology
    ScheduleTask(group_dir, local_tasks) # add tasks to scheduler and generate scripts


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
    parser.add_argument('-j', '--json_file', type=str,
                        default='./config_case.json',
                        help='Filename for json file')
    parser.add_argument('-o', '--output_dir', type=str,
                        default='../TwoDSubduction/',
                        help='Directory for output')
    parser.add_argument('-i', '--inputs', type=str,
                        default='foo',
                        help='inputs')
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
    if (_commend in ['-h', '--help']):
        # example:
        Usage()
    elif _commend == 'foo':
        pass
    elif _commend == 'plot_test_results':
        # plot the result of tests
        # example:
        # python -m shilofue.TwoDSubduction plot_test_results -i 
        #  /home/lochy/softwares/aspect/build_TwoDSubduction/tests/ -o $TwoDSubduction_DIR/test_results
        source_dir = arg.inputs
        # todo
        PlotTestResults(source_dir, output_dir=arg.output_dir)
    elif _commend == 'schedule_case_pp':
        SchedulePPCase(arg.inputs)
    elif _commend == 'schedule_group_pp':
        SchedulePPGroup(arg.inputs)
    else:
        raise ValueError('Commend %s is not available.' % _commend)


# run script
if __name__ == '__main__':
    main()