# -*- coding: utf-8 -*-
r"""(one line description)

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
import shilofue.Cases as CasesP
import shilofue.ParsePrm as ParsePrm
from shilofue.PhaseTransition import ParsePTfromJson

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')
# import utilities in subdirectiory
sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities

year = 365 * 24 * 3600.0  # yr to s

class CASE_OPT(CasesP.CASE_OPT):
    '''
    Define a class to work with CASE
    List of keys:
    '''
    def __init__(self):
        '''
        Initiation, first perform parental class's initiation,
        then perform daughter class's initiation.
        see document (run with -h option) for detail
        '''
        CasesP.CASE_OPT.__init__(self)
        self.start = self.number_of_keys()
        # self.add_key("foo", int, ['phase transition model'], 0, nick='if_wb')
        self.add_key("Model to use for mantle phase transitions", str,\
         ["phase transition model"], 'default', nick="phase_model")
        self.add_key("Json file to use for mantle phase transitions", str,\
         ["phase transition", "json file"], 'foo.json', nick="phase_json_path")
        self.add_key("Material model to use", str,\
         ["material model"], 'visco plastic', nick="material_model")
        self.add_key("Vertical velocity. This is assigned to the top boundary.\
if the variable fix_velocity is true, then this is prescribed",\
        float, ["vertical velocity"], -6.7601e-4, nick="vy")
        self.add_key("Resolution", int,\
         ["resolution"], 100, nick="resolution")
        self.add_key("Fix vertical velocity", int,\
         ["fix vertical velocity"], 0, nick="fix_velocity")
        

    def check(self):
        phase_model = self.values[self.start]
        phase_json_path = Utilities.var_subs(self.values[self.start+1])
        material_model = self.values[self.start + 2]
        assert(phase_model in ["default", "CDPT"])

    def to_configure_prm(self):
        '''
        Interface to configure_prm
        '''
        phase_model = self.values[self.start]
        phase_json_path = Utilities.var_subs(self.values[self.start+1])
        material_model = self.values[self.start + 2]
        vy = self.values[self.start+3]
        resolution = self.values[self.start+4]
        potential_T = self.values[4]
        _type = self.values[9] 
        fix_velocity = self.values[self.start + 5]
        return _type, phase_model, phase_json_path, material_model, vy, resolution, potential_T, fix_velocity

    def to_configure_wb(self):
        '''
        Interface to configure_wb
        '''
        pass


class CASE(CasesP.CASE):
    '''
    class for a case
    More Attributes:
    '''
    def configure_prm(self, _type, phase_model, phase_json_path, material_model, vy, resolution, potential_T, fix_velocity):
        '''
        Configure prm file
        Inputs:
            phase_model (str): model to use for phase transition
            phase_json_path (str): path of file for CDPT model.
            potential_T (float): potential temperature
        '''
        o_dict = self.idict.copy()
        # resolution
        o_dict['Geometry model']['Box']['X repetitions'] = str(resolution)
        o_dict['Geometry model']['Box']['Y repetitions'] = str(resolution)
        # vertical velocity
        if _type == 'layered':
            o_dict['Boundary velocity model']['Function']['Function expression'] = "0; %.4e" % vy
        if fix_velocity:
            # fix velocity field
            # a. fix solver scheme
            # b. prescribe the velocity
            # c. get rid of boundary velocity conditions
            o_dict['Nonlinear solver scheme'] = 'single Advection, no Stokes'
            o_dict['Prescribed Stokes solution'] = {
                "Model name": "function",
                "Velocity function": {
                    "Function expression": "0; %.4e" % vy,
                    "Variable names": "x,y"
                }
            }
            o_dict.pop('Boundary velocity model')
        # Adiabatic surface temperature
        o_dict["Adiabatic surface temperature"] = str(potential_T)
        o_dict["Boundary temperature model"]["Box"]["Top temperature"] = str(potential_T)
        # modify material model
        # then merge it into the material model of the prm file.
        outputs = {}
        o_dict['Material model']['Model name'] = material_model  # material model to use
        # subsection setup
        if material_model == "visco plastic twod":
            material_model_subsection = "Visco Plastic TwoD"
        else:
            material_model_subsection = material_model.title()  # convert the first letter to capital
        o_dict['Material model'][material_model_subsection] = {**o_dict['Material model'][material_model_subsection], **outputs}  # prepare entries
        self.idict = o_dict
        pass

    def configure_wb(self):
        '''
        Configure wb file
        '''
        pass


def Usage():
    print("\
(One liner description\n\
\n\
Examples of usage: \n\
\n\
  - create case with json file: \n\
\n\
        Lib_LatentHeatBK0_Cases create_with_json -j foo.json \n"
    )


def ShowJsonOption():
    Case_Opt = CASE_OPT()
    print("\
  - options defined in the json file:\n\
        %s\n\
        " % Case_Opt.document_str())


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
    parser.add_argument('-j', '--json', type=str,
                        default='',
                        help='path to a json file')
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    # commands
    if (_commend in ['-h', '--help']):
        # example:
        Usage()
    elif (_commend in ['--json_option', '-jo']):
        # json options
        ShowJsonOption()
    elif _commend == 'create_with_json':
        # example:
        CasesP.create_case_with_json(arg.json, CASE, CASE_OPT)
    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()