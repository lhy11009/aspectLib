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
        self.add_key("Compositon to model with", str, ['composition'], 'pyrolite', nick='composition')
        self.add_key("Width of the Box", float, ["box width"], 6.783e6, nick='box_width')
    
    def check(self):
        _type = self.values[9] 
        composition = self.values[self.start]
        if _type not in ["CDPT", "bd_lsolver", "bd_lsolver_test"]:
            raise NotImplementedError('Type %s is not implemented.' % _type)
        if _type == "CDPT":
            assert(composition in ['pyrolite', 'basalt', 'harzburgite'])
        pass

    def to_configure_prm(self):
        '''
        Interface to configure_prm
        '''
        if_wb = self.values[8]
        geometry = self.values[3]
        material_model = self.values[10]
        _type = self.values[9] 
        type_bc_v = self.values[13] # todo_bc
        dimension = self.values[14]
        composition = self.values[self.start]
        box_width = self.values[self.start + 1]

        return if_wb, geometry, material_model, _type, composition, type_bc_v, dimension, box_width
        

    def to_configure_wb(self):
        '''
        Interface to configure_wb
        '''
        if_wb = self.values[8]
        geometry = self.values[3]
        material_model = self.values[10]
        _type = self.values[9] 
        return if_wb, geometry, material_model, _type


class CASE(CasesP.CASE):
    '''
    class for a case
    More Attributes:
    '''
    def configure_prm(self, if_wb, geometry, material_model, _type, composition, type_bc_v, dimension, box_width):
        '''
        Configure prm file
        '''
        Ro = 6371e3 # parameters
        o_dict = self.idict.copy()
        o_dict["Dimension"] = str(dimension)
        # geometry options
        if _type in ["bd_lsolver", "bd_lsolver_test"]:
            if geometry == 'chunk' and dimension == 2:
                max_phi = box_width / Ro * 180.0 / np.pi  # extent in term of phi
                if abs((max_phi - 61.0)/max_phi):
                    str_max_phi = "61.0"
                else:
                    str_max_phi = "%.4e" % max_phi
                o_dict["Geometry model"] = {
                    "Model name": "chunk",
                    "Chunk": {
                        "Chunk inner radius": "3.481e6",\
                        "Chunk outer radius": "6.371e6",\
                        "Chunk maximum longitude": str_max_phi,\
                        "Chunk minimum longitude": "0.0",\
                        "Longitude repetitions": "2"
                    }
                }
            elif geometry == 'box' and dimension == 2:
                o_dict["Geometry model"] = {
                    "Model name": "box",
                    "Box": {
                        "X extent": "%.4e" % box_width,
                        "Y extent": "2.8900e6",
                        "X repetitions": "2"
                    }
                }
            elif geometry == 'box' and dimension == 3:
                o_dict["Geometry model"] = {
                    "Model name": "box",
                    "Box": {
                        "X extent": "%.4e" % box_width,
                        "Y extent": "2.8900e6",
                        "Z extent": "1.5000e6",
                        "X repetitions": "4",
                        "Y repetitions": "2",
                        "Z repetitions": "1"
                    }
                }
            elif geometry == 'chunk' and dimension == 3:
                max_phi = box_width / Ro * 180.0 / np.pi  # extent in term of phi
                if abs((max_phi - 61.0)/max_phi):
                    str_max_phi = "61.0"
                else:
                    str_max_phi = "%.4e" % max_phi
                o_dict["Geometry model"] = {
                    "Model name": "chunk",
                    "Chunk": {
                        "Chunk inner radius": "3.481e6",\
                        "Chunk outer radius": "6.371e6",\
                        "Chunk maximum longitude": str_max_phi,\
                        "Chunk minimum longitude": "0.0",\
                        "Chunk maximum latitude": str_max_phi,\
                        "Chunk minimum latitude": "0.0",\
                        "Longitude repetitions": "2",\
                        "Latitude repetitions": "2"
                    }
                }
            else:
                raise NotImplementedError()
        
        # refinement
        # don't forget to change the repitition when you change model domain
        # o_dict["Mesh refinement"] =
        # repitition
        # o_dict['Geometry model']['Box']['X repetitions'] = str(resolution)
        # o_dict['Geometry model']['Box']['Y repetitions'] = str(resolution)

        # velocity boundary
        if _type in ["bd_lsolver", "bd_lsolver_test"]:
            o_dict["Boundary velocity model"] = CasesP.SetBcVelocity(o_dict["Boundary velocity model"], dimension, type_bc_v)
        # boundary temperature model
        # o_dict['Boundary temperature model'] =
        
        # Adiabatic surface temperature
        # o_dict["Adiabatic surface temperature"] = str(potential_T)
        # o_dict["Boundary temperature model"]["Box"]["Top temperature"] = str(potential_T)

        # materical model
        # a. modify material model
        # b. then merge it into the material model of the prm file.
        outputs = {}
        o_dict['Material model']['Model name'] = material_model  # material model to use
        # subsection setup
        if material_model == "visco plastic twod":
            material_model_subsection = "Visco Plastic TwoD"
        else:
            material_model_subsection = material_model.title()  # convert the first letter to capital
        o_dict['Material model'][material_model_subsection] = {**o_dict['Material model'][material_model_subsection], **outputs}  # prepare entries
        
        # initial conditions 
        if _type == "CDPT":
            if composition == 'pyrolite':
                o_dict['Initial composition model']['Function']['Function expression'] = "0.0 ; 0.0; 0.0 ; 0.0"
            if composition == 'basalt':
                o_dict['Initial composition model']['Function']['Function expression'] = "1.0 ; 0.0; 0.0 ; 0.0"
            if composition == 'harzburgite':
                o_dict['Initial composition model']['Function']['Function expression'] = "0.0 ; 1.0; 0.0 ; 0.0"
        self.idict = o_dict


    def configure_wb(self, if_wb, geometry, material_model, _type):
        '''
        Configure wb file
        '''
        if not if_wb:
            # check first if we use wb file for this one
            return
        # geometry options

        pass


def Usage():
    print("\
(One liner description\n\
\n\
Examples of usage: \n\
\n\
  - create case with json file: \n\
\n\
        Lib_PhaseDiagram0_Cases create_with_json -j \
        /home/lochy/ASPECT_PROJECT/TwoDSubduction/wb_create_test/configure_1.json"
        )


def ShowJsonOption():
    Case_Opt = CASE_OPT()
    print("\
  - options defined in the json file:\n\
        %s\n\
        " % Case_Opt.document_str()
        )


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