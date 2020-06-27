import os
import sys
import json
import argparse
from shilofue.Parse import COMPOSITION
from shilofue.Parse import ParseFromDealiiInput
from shilofue.Parse import ParseToDealiiInput
from shilofue.Parse import CASE, GROUP_CASE
from shilofue.Rheology import GetLowerMantleRheology
from shilofue.Utilities import my_assert


_ALL_AVAILABLE_OPERATIONS = ['LowerMantle', 'MeshRefinement', 'query']  # all the possible operations


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
        _initial_adaptive_refinement = int(_config['initial_adaptive_refinement'])
        Inputs['Mesh refinement']['Initial adaptive refinement'] = str(_initial_adaptive_refinement)
    except KeyError:
        pass


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
    _available_commends = ['create', 'create_group', 'plot']  # only these commends are available now
    _commend = sys.argv[1]
    if _commend not in _available_commends:
        raise ValueError('Commend %s is not available.' % _commend)
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
        with open(arg.base_file, 'r') as fin:
            _inputs = ParseFromDealiiInput(fin)
        with open(arg.json_file, 'r') as fin:
            _config = json.load(fin)
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
        _case_names = MyGroup(_odir, operations=_operations, extra=_extra)
        print(_group_name)
        for _case_name in _case_names:
            # ouptut to screen
            print(_case_name)

    elif _commend == 'create':
        print('Now we create a single case:')  # screen output
        # create a case
        # read files
        with open(arg.base_file, 'r') as fin:
            _inputs = ParseFromDealiiInput(fin)
        with open(arg.json_file, 'r') as fin:
            _config = json.load(fin)
        if not os.path.isdir(arg.output_dir):
            os.mkdir(arg.output_dir)
        # Initial a case
        MyCase = MYCASE(_inputs, config=_config['config'], test=_config['test'])
        # call __call__ function to generate
        _extra = _config.get('extra', {})
        if arg.operations_file is None:
            # take all availale operations
            _operations = _ALL_AVAILABLE_OPERATIONS
        _case_name = MyCase(dirname=arg.output_dir, extra=_config['extra'], operations=_operations)
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

    elif _commend == 'update':
        # update a case
        # todo
        pass
    elif _commend == 'plot':
        # todo
        # plot something
        pass


# run script
if __name__ == '__main__':
    main()