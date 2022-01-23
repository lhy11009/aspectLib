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
from shilofue.PhaseTransition import ParsePhaseTransitionFile

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
         ["phase transition model"], 'CDPT', nick="phase_model")
        self.add_key("Json file to use for mantle phase transitions", str,\
         ["phase transition", "json file"], 'foo.json', nick="phase_json_path")
    
    def check(self):
        phase_model = self.values[self.start]
        phase_json_path = Utilities.var_subs(self.values[self.start+1])
        if phase_model == "CDPT":
            assert(os.path.isfile(phase_json_path))
        else:
            raise ValueError("check: value for phase model is not valid.")
        pass

    def to_configure_prm(self):
        '''
        Interface to configure_prm
        '''
        phase_model = self.values[self.start]
        phase_json_path = Utilities.var_subs(self.values[self.start+1])
        return phase_model, phase_json_path        

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
    def configure_prm(self, phase_model, phase_json_path):
        '''
        Configure prm file
        Inputs:
            phase_model (str): model to use for phase transition
            phase_json_path (str): path of file for CDPT model.
        '''
        o_dict = self.idict.copy()
        if phase_model == "CDPT":
            outputs = ParsePhaseTransitionFile(phase_json_path) 
        o_dict['Material model']['Visco Plastic TwoD'] = {**o_dict['Material model']['Visco Plastic TwoD'], **outputs}
        self.idict = o_dict
        pass

    def configure_wb(self):
        '''
        Configure wb file
        '''
        pass


def Usage():
    Case_Opt = CASE_OPT()
    print("\
(One liner description\n\
\n\
Examples of usage: \n\
\n\
  - create case with json file: \n\
\n\
        Lib_LatentHeatBK0_Cases create_with_json -j foo.json \n\
\n\
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
    elif _commend == 'create_with_json':
        # example:
        CasesP.create_case_with_json(arg.json, CASE, CASE_OPT)
    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()