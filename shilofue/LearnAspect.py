import os
import sys
import json
import argparse
from shilofue.Parse import COMPOSITION
from shilofue.Parse import ParseFromDealiiInput
from shilofue.Parse import ParseToDealiiInput
from shilofue.Parse import CASE, GROUP_CASE, UpdateProjectMd, UpdateProjectJson, AutoMarkdownCase, AutoMarkdownGroup
from shilofue.Rheology import GetLowerMantleRheology
from shilofue.Utilities import my_assert
from shilofue.Doc import UpdateProjectDoc
from shilofue.Plot import ProjectPlot

_ALL_AVAILABLE_OPERATIONS = ['LowerMantle', 'MeshRefinement', 'Gravity', 'query', 'SinkingBlob']  # all the possible operations


def LowerMantle(Inputs, jump, T, P, V1):
    """
    LowerMantle(Inputs)

    calculate flow law parameters
    """
    # parse from input
    visco_plastic = Inputs["Material model"]['Visco Plastic']
    prefactors_for_diffusion_creep = COMPOSITION(visco_plastic["Prefactors for diffusion creep"])
    grain_size = COMPOSITION(visco_plastic["Grain size"])
    grain_size_exponents_for_diffusion_creep  = COMPOSITION(visco_plastic["Grain size exponents for diffusion creep"])
    activation_energies_for_diffusion_creep = COMPOSITION(visco_plastic["Activation energies for diffusion creep"])
    activation_volumes_for_diffusion_creep  = COMPOSITION(visco_plastic["Activation volumes for diffusion creep"])
    # call GetLowerMantleRheology to derive parameters for lower mantle flow law
    backgroud_upper_mantle_diffusion = {}
    backgroud_upper_mantle_diffusion['A'] = prefactors_for_diffusion_creep.data['background'][0]
    backgroud_upper_mantle_diffusion['d'] = grain_size.data['background'][0]
    backgroud_upper_mantle_diffusion['n'] = 1.0
    backgroud_upper_mantle_diffusion['m'] = grain_size_exponents_for_diffusion_creep.data['background'][0]
    backgroud_upper_mantle_diffusion['E'] = activation_energies_for_diffusion_creep.data['background'][0]
    backgroud_upper_mantle_diffusion['V'] = activation_volumes_for_diffusion_creep.data['background'][0]
    backgroud_lower_mantle_diffusion = GetLowerMantleRheology(backgroud_upper_mantle_diffusion, jump, T, P, V1=V1, strategy='d')
    # todo: add in choice of phases
    prefactors_for_diffusion_creep.data['background'] = [backgroud_upper_mantle_diffusion['A'], backgroud_lower_mantle_diffusion['A']]
    grain_size.data['background'] = [backgroud_upper_mantle_diffusion['d'], backgroud_lower_mantle_diffusion['d']]
    grain_size_exponents_for_diffusion_creep.data['background'] = [backgroud_upper_mantle_diffusion['m'], backgroud_lower_mantle_diffusion['m']]
    activation_energies_for_diffusion_creep.data['background'] = [backgroud_upper_mantle_diffusion['E'], backgroud_lower_mantle_diffusion['E']]
    activation_volumes_for_diffusion_creep.data['background'] = [backgroud_upper_mantle_diffusion['V'], backgroud_lower_mantle_diffusion['V']]
    # parse back
    visco_plastic["Prefactors for diffusion creep"] = prefactors_for_diffusion_creep.parse_back()
    visco_plastic["Grain size"] = grain_size.parse_back()
    visco_plastic["Grain size exponents for diffusion creep"] = grain_size_exponents_for_diffusion_creep.parse_back()
    visco_plastic["Activation energies for diffusion creep"] = activation_energies_for_diffusion_creep.parse_back()
    visco_plastic["Activation volumes for diffusion creep"] = activation_volumes_for_diffusion_creep.parse_back()
    return Inputs


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
        # longitude repetitions for chunk geometry
        _longitude_repetitions = int(_config['longitude_repetitions'])
    except KeyError:
        pass
    else:
        Inputs['Geometry model']['Chunk']['Longitude repetitions'] = str(_longitude_repetitions)

def Gravity(Inputs, _config):
    '''
    change the gravity properties
    '''
    try:
        # Initial global refinement
        _gravity_magnitude = int(_config['gravity_magnitude'])
    except KeyError:
        pass
    else:
        # todo
        Inputs['Gravity model']['Vertical']['Magnitude'] = str(_gravity_magnitude)

def SinkingBlob(Inputs, _config):
    '''
    change model width of sinking blob test case
    '''
    _cell_width = 0.1  # width of each cell in the model
    try:
        _sinking_blob_model_width = float(_config['sinking_blob_model_width'])
    except KeyError:
        pass
    else:
        Inputs['Geometry model']['Box']['X extent'] = str(_sinking_blob_model_width)
        Inputs['Geometry model']['Box']['X repetitions'] = str(int(_sinking_blob_model_width/_cell_width))
        Inputs['Mesh refinement']['Minimum refinement function']['Function constants'] = 'Xc=%s, Yc=0.8, Rc=0.1' % (_sinking_blob_model_width / 2.0)
        Inputs['Initial composition model']['Function']['Function constants'] = 'Xc=%s, Yc=0.8, Rc=0.1' % (_sinking_blob_model_width / 2.0)
                                                                        


def Parse(ifile, ofile):
    """
    todo
    """
    assert(os.access(ifile, os.R_OK))
    # todo
    with open(ifile, 'r') as fin:
        inputs = ParseFromDealiiInput(fin)
    # todo
    LowerMantle(inputs, 30.0, 1663.0, 21e9, 1.5e-6)
    # todo
    with open(ofile, 'w') as fout:
        ParseToDealiiInput(fout, inputs)


class MYCASE(CASE):
    '''
    Inherit from class CASE in Parse.py
    '''
    def Intepret(self, **kwargs):
        '''
        Intepret configuration for my TwoDSubduction Cases
        kwargs:
            extra(dict) - extra configuration
        '''
        _extra = kwargs.get('extra', {})
        _operations = kwargs.get('operations', [])
        _config = { **self.config, **self.test, **_extra }  # append _extra to self.config
        if 'LowerMantle' in _operations:
            # change lower mantle viscosity
            try:
                LowerMantle(self.idict, _config['upper_lower_viscosity'], _config['T660'], _config['P660'], _config['LowerV'])
            except KeyError:
                pass
        if 'MeshRefinement' in _operations:
            # change initial mesh refinement
            # check the key exists in the dict: fixed in the function
            MeshRefinement(self.idict, _config)
        if 'Gravity' in _operations:
            Gravity(self.idict, _config)
        if 'SinkingBlob' in _operations:
            # change model width for sinking blob model
            SinkingBlob(self.idict, _config)



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
    parser = argparse.ArgumentParser(description='LearnAspect Project')
    parser.add_argument('-b', '--base_file', type=str,
                        default='./files/LearnAspect/base.prm',
                        help='Filename for base file')
    parser.add_argument('-U', '--use_basename_as_base_file', type=int,
                        default=1,
                        help='Whether we use basename as base file')
    parser.add_argument('-j', '--json_file', type=str,
                        default='./config_case.json',
                        help='Filename for json file')
    parser.add_argument('-o', '--output_dir', type=str,
                        default='../LearnAspect/',
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
        # read base prm file
        if arg.use_basename_as_base_file == 1:
            _filename = './files/LearnAspect/%s.prm' % _base_name
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
        if arg.operations_file is None:
            # take all availale operations
            _operations = _ALL_AVAILABLE_OPERATIONS
        _case_names = MyGroup(_odir, operations=_operations, extra=_extra, basename=_base_name)
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
        # read base prm file
        _base_name = _config.get('basename', '')
        if arg.use_basename_as_base_file == 1:
            _filename = './files/LearnAspect/%s.prm' % _base_name
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
        if arg.operations_file is None:
            # take all availale operations
            _operations = _ALL_AVAILABLE_OPERATIONS
        # also get extra files
        _extra_files = _config.get('extra_file', {})
        _case_name = MyCase(dirname=arg.output_dir, extra=_config['extra'], operations=_operations, basename=_base_name, extra_file=_extra_files)
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
        UpdateProjectDoc(_project_dict, _project_dir, images=['Statistics' ,'DepthAverage', 'visit'])

    elif _commend == 'plot':
        # todo_future
        # plot something
        pass


# run script
if __name__ == '__main__':
    main()
