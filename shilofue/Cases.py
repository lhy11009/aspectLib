# -*- coding: utf-8 -*-
r"""class with basic interfaces for manipulating an aspect case

Logistics:
    I mean to define every feature in subdfolders, under the folder of a project (also named Cases.py).
    I expect this processes to be different from model to model (e.g. how to change the size of the domain).
    As we don't want to repeat the process of changing prms in the prm file, but to summariz possible settings here.

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
from shutil import copy2, rmtree
from copy import deepcopy
# from matplotlib import cm
from matplotlib import pyplot as plt
import shilofue.ParsePrm as ParsePrm

# directory to the aspect Lab
ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']
RESULT_DIR = os.path.join(ASPECT_LAB_DIR, 'results')
# directory to shilofue
shilofue_DIR = os.path.join(ASPECT_LAB_DIR, 'shilofue')

sys.path.append(os.path.join(ASPECT_LAB_DIR, 'utilities', "python_scripts"))
import Utilities


class CASE_OPT(Utilities.JSON_OPT):
    '''
    Define a class to work with CASE
    List of keys:
    '''
    def __init__(self):
        '''
        Initiation, first perform parental class's initiation,
        then perform daughter class's initiation.
        '''
        Utilities.JSON_OPT.__init__(self)
        # todo
        self.add_key("Name of the case", str, ["name"], "foo", nick='name')
        self.add_key("Base directory (inputs)", str, ["base directory"], ".", nick='base_dir')
        self.add_key("Output directory", str, ["output directory"], ".", nick='o_dir')
        self.add_key("Geometry", str, ["geometry"], "chunk", nick='geometry')
        self.add_key("Potential temperature of the mantle", float,\
            ["Potential temperature"], 1673.0, nick='potential_T')
        pass
    
    def check(self):
        '''
        check to see if these values make sense
        '''
        # output and input dirs
        os.path.isdir(self.values[1])
        os.path.isdir(self.values[2])
        pass

    # todo 
    def to_init(self):
        '''
        Interface to init
        '''
        inputs = os.path.join(self.values[1], 'case.prm')
        return self.values[0], inputs

    def wb_inputs_path(self):
        '''
        Interface to wb_inputs
        '''
        wb_inputs = os.path.join(self.values[1], 'case.wb')
        return wb_inputs
    
    def o_dir(self):
        '''
        Interface to output dir
        '''
        return self.values[2]


class CASE():
    '''
    class for a case
    Attributes:
        name(str):
            list for name of variables to change
        idict(dict):
            dictionary of parameters
        wb_dict(idit):
            dictionary of world builder options
        extra_files(array):
            an array of extra files of this case
    '''
    # future: add interface for extra
    def __init__(self, case_name, inputs, **kwargs):
        '''
        initiate from a dictionary
        Inputs:
            idict(dict):
                dictionary import from a base file
            kwargs:
                wb_inputs(dict or str):
                    inputs from a world builder file
        '''
        self.case_name = case_name
        self.extra_files = []
        if type(inputs)==dict:
            # direct read if dict is given
            print("    Read inputs from a dictionary")
            self.idict = deepcopy(inputs)
        elif type(inputs)==str:
            # read from file if file path is given. This has the virtual that the new dict is indepent of the previous one.
            print("    Read inputs from %s" % inputs) 
            with open(inputs, 'r') as fin:
                self.idict = ParsePrm.ParseFromDealiiInput(fin)
            pass
        else:
            raise TypeError("Inputs must be a dictionary or a string")
        # read world builder
        wb_inputs = kwargs.get('wb_inputs', {})
        if type(wb_inputs) == dict:
            # direct read if dict is given
            print("    Read world builder options from a dictionary")
            self.wb_dict = deepcopy(wb_inputs)
        elif type(wb_inputs)==str:
            # read from file if file path is given. This has the virtual that the new dict is indepent of the previous one.
            print("    Read world builder options from %s" % wb_inputs)
            with open(wb_inputs, 'r') as fin:
                self.wb_dict = json.load(fin)
            pass
        else:
            raise TypeError("CASE:%s: wb_inputs must be a dictionary or a string" % Utilities.func_name())
    
    def create(self, _root, **kwargs):
        '''
        create a new case
        Inputs:
            _root(str): a directory to put the new case
            **kwargs:
                "fast_first_step": generate another file for fast running the 0th step
        '''
        # folder
        case_dir = os.path.join(_root, self.case_name)
        if os.path.isdir(case_dir):
           # remove old ones 
           rmtree(case_dir)
        os.mkdir(case_dir)
        # file output
        prm_out_path = os.path.join(case_dir, "case.prm")  # prm for running the case
        wb_out_path = os.path.join(case_dir, "case.wb")  # world builder file
        ParsePrm.WritePrmFile(prm_out_path, self.idict)
        # fast first step
        fast_first_step = kwargs.get('fast_first_step', 0) 
        if fast_first_step == 1:
            outputs = deepcopy(self.idict)
            prm_fast_out_path = os.path.join(case_dir, "case_f.prm")
            ParsePrm.FastZeroStep(outputs)  # generate another file for fast running the 0th step
            ParsePrm.WritePrmFile(prm_fast_out_path, outputs)
        for path in self.extra_files:
            base_name = os.path.basename(path)
            path_out = os.path.join(case_dir, base_name)
            copy2(path, path_out)
        # world builder
        if self.wb_dict != {}:
            with open(wb_out_path, 'w') as fout:
                json.dump(self.wb_dict, fout, indent=2)
        print("New case created: %s" % case_dir)

    def configure(self, func, config, **kwargs):
        '''
        applies configuration for this case
        Inputs:
            func(a function), the form of it is:
                outputs = func(inputs, config)
        '''
        rename = kwargs.get('rename', None)
        if rename != None:
            # rename the case with the configuration
            self.idict, appendix = func(self.idict, config)
            self.case_name += appendix
        else:
            # just apply the configuration
            self.idict = func(self.idict, config)
    
    def configure_wb(self):
        '''
        Configure world builder file
        '''
        pass

    def add_extra_file(self, path):
        '''
        add an extra file to list
        Inputs:
            path(str): an extra file
        '''
        self.extra_files.append(path)


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
    if _commend == 'foo':
        # example:
        SomeFunction('foo')

# run script
if __name__ == '__main__':
    main()