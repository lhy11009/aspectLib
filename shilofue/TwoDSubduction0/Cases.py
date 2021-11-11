# -*- coding: utf-8 -*-
r"""class with interfaces for manipulating an aspect case in TwoDSubduction Project

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
import json, re
# import pathlib
# import subprocess
import numpy as np
# from matplotlib import cm
from matplotlib import pyplot as plt
import shilofue.Cases as CasesP

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')
# import utilities in subdirectiory
sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities


class CASE_OPT(CasesP.CASE_OPT):
    '''
    Define a class to work with CASE
    List of keys:
    '''
    def __init__(self):
        '''
        Initiation, first perform parental class's initiation,
        then perform daughter class's initiation.
        List of keys:
            0: if_wb (int) - If use world builder
            1: sp_age_trench (float) - Age of the subducting plate at trench (Ma)
            2: sp_rate - Spreading rate of the subducting plate (m/ yr)
            3: ov_age (float) - Age of the overiding plate (Ma)
        '''
        CasesP.CASE_OPT.__init__(self)
        self.add_key("If use world builder", int, ['use world builder'], 0)
        self.add_key("Age of the subducting plate at trench", float,\
            ['world builder', 'sp age trench'], 80e3)
        self.add_key("Spreading rate of the subducting plate", float,\
            ['world builder', 'sp rate'], 0.05)
        self.add_key("Age of the overiding plate", float,\
            ['world builder', 'ov age'], 40e3)
        pass
    
    def check(self):
        '''
        check to see if these values make sense
        '''
        CasesP.CASE_OPT.check(self)
        pass

    def to_configure_wb(self):
        '''
        Interface to configure_wb
        '''
        if_wb = self.values[0]
        sp_age_trench = self.values[1]
        sp_rate = self.values[2]
        ov_age = self.values[3]
        return if_wb, sp_age_trench, sp_rate, ov_age


class CASE(CasesP.CASE):
    '''
    class for a case
    More Attributes:
    '''
    def configure_wb(self, if_wb, sp_age_trench, sp_rate, ov_ag):
        '''
        Configure world builder file
        Inputs:
            see description of CASE_OPT
        '''
        Ro = self.idict['Geometry model']['Chunk']['Chunk outer radius']
        if not if_wb:
            # check first if we use wb file for this one
            return
        self.wb_dict = wb_configure_plates_sph(self.wb_dict, sp_age_trench, sp_rate, ov_ag, Ro) # plates
        pass


def wb_configure_plates_sph(wb_dict, sp_age_trench, sp_rate, ov_ag, Ro):
    '''
    configure plate in world builder
    '''
    return wb_dict

        


def Usage():
    print("\
(One liner description\n\
\n\
Examples of usage: \n\
\n\
  - default usage: \n\
\n\
        python -m \
        ")

def SomeFunction(foo):
    '''
    descriptions
    Inputs:
        -
    Returns:
        -
    '''
    pass


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
    elif _commend == 'foo':
        # example:
        SomeFunction('foo')
    else:
        # no such option, give an error message
        raise ValueError('No commend called %s, please run -h for help messages' % _commend)

# run script
if __name__ == '__main__':
    main()