import os
import sys
import json
from shilofue.Parse import COMPOSITION
from shilofue.Parse import ParseFromDealiiInput
from shilofue.Parse import ParseToDealiiInput
from shilofue.Parse import GROUP_CASE
from shilofue.Rheology import GetLowerMantleRheology



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
    prefactors_for_diffusion_creep.data['background'].append(backgroud_lower_mantle_diffusion['A'])
    grain_size.data['background'].append(backgroud_lower_mantle_diffusion['d'])
    grain_size_exponents_for_diffusion_creep.data['background'].append(backgroud_lower_mantle_diffusion['m'])
    activation_energies_for_diffusion_creep.data['background'].append(backgroud_lower_mantle_diffusion['E']) 
    activation_volumes_for_diffusion_creep.data['background'].append(backgroud_lower_mantle_diffusion['V'])
    # parse back
    visco_plastic["Prefactors for diffusion creep"] = prefactors_for_diffusion_creep.parse_back()
    visco_plastic["Grain size"] = grain_size.parse_back()
    visco_plastic["Grain size exponents for diffusion creep"] = grain_size_exponents_for_diffusion_creep.parse_back()
    visco_plastic["Activation energies for diffusion creep"] = activation_energies_for_diffusion_creep.parse_back()
    visco_plastic["Activation volumes for diffusion creep"] = activation_volumes_for_diffusion_creep.parse_back()
    return Inputs

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


def GroupParse(_json_file, _ifile, _target_dir):
    '''
    Parse a group of input file with the parameters defined in the jsonfile
    Inputs:
        _json_flle(str): 
            filename for a json file
        _ifile(str):
            filename for a input .prm, as base for the group of cases
        _target_dir(str):
            name of the target directory
        kwarg:
            'test'(True or False):
                if this is a group of test case
    '''
    assert(os.access(_json_file, os.R_OK))
    _idict = json.load(_json_file)
    GroupCase = GROUP_CASE(_idict)
    GroupCase.Parse(_target_dir)
    with open(_ifile, 'r') as fin:
        _inputs_base = ParseFromDealiiInput(fin)

