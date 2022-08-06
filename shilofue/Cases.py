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
from shilofue.Rheology import RHEOLOGY_OPR

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
        self.add_key("Name of the case", str, ["name"], "foo", nick='name')
        self.add_key("Base directory (inputs)", str, ["base directory"], ".", nick='base_dir')
        self.add_key("Output directory", str, ["output directory"], ".", nick='o_dir')
        self.add_key("Geometry", str, ["geometry"], "chunk", nick='geometry')
        self.add_key("potential temperature of the mantle", float,\
            ["potential temperature"], 1673.0, nick='potential_T')
        self.add_key("include fast first step", int,\
            ["include fast first step"], 0, nick='if_fast_first_step')
        self.add_key("Additional files to include", list,\
            ["additional files"], [], nick='additional_files')
        self.add_key("Root level from the project root", int,\
         ["root level"], 1, nick="root_level")
        self.add_key("If use world builder", int, ['use world builder'], 0, nick='if_wb')
        self.add_key("Type of the case", str, ["type"], '', nick='_type')
        self.add_key("Material model to use", str,\
         ["material model"], 'visco plastic', nick="material_model")
        self.add_key("Linear solver toleracne", float,\
         ["stokes solver", "linear solver tolerance"], 0.1, nick="stokes_linear_tolerance")
        self.add_key("End time", float, ["end time"], 60e6, nick="end_time")
        # todo_bc
        self.add_key("Type of velocity boundary condition\n\
            available options in [all fs, bt fs side ns]", str,\
            ["boundary condition", "velocity", "type"], "all fs", nick='type_bd_v')
        self.add_key("Dimension", int, ['dimension'], 2, nick='dimension')
        # todo_affinity
        self.add_key("Refinement level, note this is a summarized parameter of the refinement scheme assigned,\
it only takes effect if the input is positiveh",\
            int, ["refinement level"], -1, nick="refinement_level")
        self.add_key("Case Output directory", str, ["case output directory"], "output", nick='case_o_dir')
        self.add_key("mantle rheology", str, ['mantle rheology', 'scheme'], "HK03_wet_mod", nick='mantle_rheology_scheme')
        pass
    
    def check(self):
        '''
        check to see if these values make sense
        '''
        # output and input dirs
        base_dir = Utilities.var_subs(self.values[1])
        o_dir = Utilities.var_subs(self.values[2])
        Utilities.my_assert(os.path.isdir(base_dir), FileNotFoundError, "No such directory: %s" % base_dir)
        # in case this is "", we'll fix that later.
        Utilities.my_assert(o_dir=="" or os.path.isdir(o_dir), FileNotFoundError, "No such directory: %s" % o_dir)
        pass

    def to_init(self):
        '''
        Interface to init
        '''
        _type = self.values[9]
        base_dir = self.values[1]
        if _type == '':
            base_name = 'case.prm'
        else:
            base_name = 'case_%s.prm' % _type
        inputs = os.path.join(base_dir, base_name)
        if_wb = self.values[8]
        return self.values[0], inputs, if_wb

    def wb_inputs_path(self):
        '''
        Interface to wb_inputs
        '''
        _type = self.values[9]
        if _type == '':
            base_name = 'case.wb'
        else:
            base_name = 'case_%s.wb' % _type
        wb_inputs = os.path.join(self.values[1], base_name)
        return wb_inputs
    
    def o_dir(self):
        '''
        Interface to output dir
        '''
        return Utilities.var_subs(self.values[2])
    
    def case_name(self):
        '''
        Return name of the case
        Return:
            case name (str)
        '''
        return self.values[0]
    
    def get_additional_files(self):
        '''
        Interface to add_files
        '''
        files = []
        for additional_file in self.values[6]:
            _path = Utilities.var_subs(os.path.join(self.values[1], additional_file))
            Utilities.my_assert(os.access(_path, os.R_OK), FileNotFoundError,\
            "Additional file %s is not found" % _path)
            files.append(_path)
        return files
    
    def if_fast_first_step(self):
        '''
        If we generate a case with fast-first-step computation
        '''
        return self.values[5]
        pass

    def if_use_world_builder(self):
        '''
        if we use world builder
        '''
        if_wb = self.values[8]
        return  (if_wb==1)
    
    def fix_case_name(self, case_name):
        '''
        fix base dir with a new value
        '''
        self.values[0] = case_name

    def fix_base_dir(self, base_dir):
        '''
        fix base dir with a new value
        '''
        assert(os.path.isdir(base_dir))
        self.values[1] = base_dir
    
    def fix_output_dir(self, o_dir):
        '''
        fix directory to output
        '''
        self.values[2] = o_dir
        
    def reset_refinement(self, reset_refinement_level):
        '''
        reset refinement level
        '''
        # todo_affinity
        self.values[15] = reset_refinement_level
        pass
    
    def fix_case_output_dir(self, case_o_dir):
        '''
        reset refinement level
        '''
        # todo_affinity
        self.values[16] = case_o_dir
    

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
        model_stages(int):]
            stages in a model: if > 0, then create multiple prm files
    '''
    # future: add interface for extra
    def __init__(self, case_name, inputs, if_wb, **kwargs):
        '''
        initiate from a dictionary
        Inputs:
            idict(dict):
                dictionary import from a base file
            if_wb(True or False):
                if use world builder
            kwargs:
                wb_inputs(dict or str):
                    inputs from a world builder file
        '''
        self.case_name = case_name
        self.extra_files = []
        self.wb_dict = {}
        self.model_stages = 1
        self.additional_idicts = []
        self.output_files = [] # for saving the path of files and images output from this class
        self.output_imgs = []
        if type(inputs)==dict:
            # direct read if dict is given
            print("    Read inputs from a dictionary")
            self.idict = deepcopy(inputs)
        elif type(inputs)==str:
            # read from file if file path is given. This has the virtual that the new dict is indepent of the previous one.
            print("    Read inputs from %s" % Utilities.var_subs(inputs)) 
            with open(Utilities.var_subs(inputs), 'r') as fin:
                self.idict = ParsePrm.ParseFromDealiiInput(fin)
            pass
        else:
            raise TypeError("Inputs must be a dictionary or a string")
        # read world builder
        if if_wb:
            wb_inputs = kwargs.get('wb_inputs', {})
            if type(wb_inputs) == dict:
                # direct read if dict is given
                print("    Read world builder options from a dictionary")
                self.wb_dict = deepcopy(wb_inputs)
            elif type(wb_inputs)==str:
                # read from file if file path is given. This has the virtual that the new dict is indepent of the previous one.
                print("    Read world builder options from %s" % Utilities.var_subs(wb_inputs))
                with open(Utilities.var_subs(wb_inputs), 'r') as fin:
                    self.wb_dict = json.load(fin)
                pass
            else:
                raise TypeError("CASE:%s: wb_inputs must be a dictionary or a string" % Utilities.func_name())
        # operator of rheology
        self.Rheology_Opr = RHEOLOGY_OPR()
    
    def create(self, _root, **kwargs):
        '''
        create a new case
        Inputs:
            _root(str): a directory to put the new case
            **kwargs:
                "fast_first_step": generate another file for fast running the 0th step
        Return:
            case_dir(str): path to created case.
        '''
        # folder
        case_dir = os.path.join(_root, self.case_name)
        if os.path.isdir(case_dir):
           # remove old ones 
           rmtree(case_dir)
        os.mkdir(case_dir)
        output_files_dir = os.path.join(case_dir, 'configurations')
        os.mkdir(output_files_dir)
        img_dir = os.path.join(case_dir, 'img')
        os.mkdir(img_dir)
        output_img_dir = os.path.join(img_dir, 'initial_condition')
        os.mkdir(output_img_dir)
        # file output
        prm_out_path = os.path.join(case_dir, "case.prm")  # prm for running the case
        wb_out_path = os.path.join(case_dir, "case.wb")  # world builder file
        ParsePrm.WritePrmFile(prm_out_path, self.idict)
        if self.model_stages > 1:
            assert(len(self.additional_idicts) == self.model_stages-1)
            for i in range(self.model_stages-1):
                prm_out_path = os.path.join(case_dir, "case_%d.prm" % (i+1))  # prm for running the case
                ParsePrm.WritePrmFile(prm_out_path, self.additional_idicts[i])
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
        # copy paste files and figures generated
        for path in self.output_files:
            base_name = os.path.basename(path)
            path_out = os.path.join(output_files_dir, base_name)
            copy2(path, path_out)
        for path in self.output_imgs:
            base_name = os.path.basename(path)
            path_out = os.path.join(output_img_dir, base_name)
            copy2(path, path_out)
        return case_dir

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
    


def create_case_with_json(json_opt, CASE, CASE_OPT, **kwargs):
    '''
    A wrapper for the CASES class
    Inputs:
        json_opt(str, dict): path or dict a json file
        kwargs (dict):
            update (bool): update existing cases?
    Returns:
        case_dir: return case directory
    '''
    print("%s: Creating case" % Utilities.func_name())
    is_update = kwargs.get('update', True)
    fix_case_name = kwargs.get('fix_case_name', None)
    fix_base_dir = kwargs.get('fix_base_dir', None)
    fix_output_dir = kwargs.get('fix_output_dir', None)
    reset_refinement_level = kwargs.get('reset_refinement_level', None)
    fix_case_output_dir = kwargs.get('fix_case_output_dir', None)
    Case_Opt = CASE_OPT()
    if type(json_opt) == str:
        if not os.access(json_opt, os.R_OK):
            raise FileNotFoundError("%s doesn't exist" % json_opt)
        Case_Opt.read_json(json_opt)
    elif type(json_opt) == dict:
        Case_Opt.import_options(json_opt)
    else:
        raise TypeError("Type of json_opt must by str or dict")
    if fix_case_name != None:
        Case_Opt.fix_case_name(fix_case_name)  # fix base dir, useful when creating a group of case from a folder
    if fix_base_dir != None:
        Case_Opt.fix_base_dir(fix_base_dir)  # fix base dir, useful when creating a group of case from a folder
    if fix_output_dir != None:
        Case_Opt.fix_output_dir(fix_output_dir)  # fix output dir, useful when doing tests
    # todo_affinity
    if reset_refinement_level != None:
        Case_Opt.reset_refinement(reset_refinement_level)
    if fix_case_output_dir != None:
        Case_Opt.fix_case_output_dir(fix_case_output_dir)
    Case_Opt.check()
    # check if the case already exists. If so, only update if it is explicitly 
    # required
    case_dir_to_check = os.path.join(Case_Opt.o_dir(), Case_Opt.case_name())
    if os.path.isdir(case_dir_to_check):
        if is_update:
            print("Case %s already exists, updating" % case_dir_to_check)
        else:
            print("Case %s already exists, aborting" % case_dir_to_check)
            return case_dir_to_check
    # Manage case files
    Case = CASE(*Case_Opt.to_init(), wb_inputs=Case_Opt.wb_inputs_path())
    Case.configure_prm(*Case_Opt.to_configure_prm())
    if Case_Opt.if_use_world_builder():
        Case.configure_wb(*Case_Opt.to_configure_wb())
    for _path in Case_Opt.get_additional_files():
        Case.add_extra_file(_path)
    # create new case
    case_dir = Case.create(Case_Opt.o_dir(), fast_first_step=Case_Opt.if_fast_first_step())
    return case_dir

def SetBcVelocity(bc_dict, dimension, type_bc_v):
    '''
    set the boundary condition for velocity
    Inputs:
        bc_dict: an input with entries for Boundary velocity model in a prm file
        type_bc_v: type of velocity bc
    '''
    # todo_bc
    fs_indicator=[]
    ns_indicator = []
    if dimension == 2 and type_bc_v == 'all fs':
        fs_indicator = [0, 1, 2, 3]
    elif dimension == 3 and type_bc_v == 'all fs':
        fs_indicator = [0, 1, 2, 3, 4, 5]
    elif dimension == 2 and type_bc_v == 'bt fs side ns':
        fs_indicator = [1, 3]
        ns_indicator = [0, 2]
    elif dimension == 3 and type_bc_v == 'bt fs side ns':
        fs_indicator = [4, 5]
        ns_indicator = [0, 1, 2, 3]
    else:
        raise NotImplementedError("This combination of dimension (%d) and \
velocity boundary (%s) is not implemented yet." % (dimension, type_bc_v))
    if len(fs_indicator) > 0:
        bc_dict["Tangential velocity boundary indicators"] = str(fs_indicator).strip(']').lstrip('[')  # free slip
    else:
        if "Tangential velocity boundary indicators" in bc_dict:
            _ = bc_dict.pop("Tangential velocity boundary indicators")
    if len(ns_indicator) > 0:
        bc_dict["Zero velocity boundary indicators"] = str(ns_indicator).strip(']').lstrip('[')  # no slip
    else:
        if "Zero velocity boundary indicators" in bc_dict:
            _ = bc_dict.pop("Zero velocity boundary indicators")
    return bc_dict
    pass


def SetNewtonSolver(o_dict):
    '''
    Settings for the newton solver, starting from a combination that works
    '''
    o_dict["Nonlinear solver scheme"] = "single Advection, iterated Newton Stokes"
    o_dict["Max nonlinear iterations"] = "100"
    o_dict["Max nonlinear iterations in pre-refinement"] = "0"
    o_dict["Nonlinear solver tolerance"] = "1e-6"
    o_dict["Solver parameters"] = {}
    o_dict["Solver parameters"]["Newton solver parameters"] = {
        "Max pre-Newton nonlinear iterations" :"20",\
        "Nonlinear Newton solver switch tolerance": "1e-3",\
        "Max Newton line search iterations": "0",\
        "Maximum linear Stokes solver tolerance": "0.9",\
        "Use Newton residual scaling method": "true",\
        "Use Newton failsafe": "true",\
        "Stabilization preconditioner": "SPD",\
        "Stabilization velocity block": "SPD",\
        "Use Eisenstat Walker method for Picard iterations": "true"
    }
    o_dict["Solver parameters"]["Stokes solver parameters"] = {
        "Maximum number of expensive Stokes solver steps": "5000",\
        "Number of cheap Stokes solver steps": "500",\
        "Linear solver tolerance": "1e-1",\
        "GMRES solver restart length": "100"
    }
    return o_dict
    


def GROUP():
    '''
    A group of cases
    Attributes
        iCASE (class) - class of the case to use
    '''
    def __init__(self, iCASE):
        '''
        Input:
            iCASE(class of case)
        '''
        self.iCASE = iCASE


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