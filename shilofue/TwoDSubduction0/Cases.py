# -*- coding: utf-8 -*-
r"""class with interfaces for manipulating an aspect case in TwoDSubduction Project

Logistics:
    I mean to define every feature here, under the folder of a project.
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
        List of keys:
            0: if_wb (int) - If use world builder
            1: sp_age_trench (float) - Age of the subducting plate at trench (yr)
            2: sp_rate - Spreading rate of the subducting plate (m/ yr)
            3: ov_age (float) - Age of the overiding plate (yr)
            4: ov_trans_age (float) - Age of the transit overiding plate (yr)
                this is assigned to -1.0 by default, which means we don't use
                transit plate
            5: ov_trans_length (float) - Length of the transit overiding plate (m)
            6: type_of_bd (str) - "Type of boundary condition"
            7: potential_T (float) - Potential temperature of the mantle
            8: box_width (float) - "Width of the Box"
            9: geometry (str) - "Geometry"
            10: base_dir - Base directory (inputs)
            11: o_dir - Output directory
            12 name - Name of the case
        '''
        CasesP.CASE_OPT.__init__(self)
        self.add_key("If use world builder", int, ['use world builder'], 0)
        self.add_key("Age of the subducting plate at trench", float,\
            ['world builder', 'subducting plate','age trench'], 80e6)
        self.add_key("Spreading rate of the subducting plate", float,\
            ['world builder', 'subducting plate', 'sp rate'], 0.05)
        self.add_key("Age of the overiding plate", float,\
            ['world builder', "overiding plate", 'age'], 40e6)
        self.add_key("Age of the transit overiding plate", float,\
            ['world builder', "overiding plate", "transit", 'age'], -1.0)
        self.add_key("Length of the transit overiding plate", float,\
            ['world builder', "overiding plate", "transit", 'length'], 300e3)
        self.add_key("Type of boundary condition\n\
            available options in [all free slip, ]", str,\
            ["boundary condition", "model"], "all free slip")
        self.add_key("Potential temperature of the mantle", float,\
            ["Potential temperature"], 1673.0)
        self.add_key("Width of the Box", float,\
            ["Box Width"], 6.783e6)
        self.add_key("Geometry, \"chunk\" or \"box\"", str, ["geometry"], "chunk")
        self.add_key("Base directory (inputs)", str, ["base directory"], ".")
        self.add_key("Output directory", str, ["output directory"], ".")
        self.add_key("Name of the case", str, ["name"], ".")
        pass
    
    def check(self):
        '''
        check to see if these values make sense
        '''
        CasesP.CASE_OPT.check(self)
        # geometry options
        assert(self.values[9] in ['chunk', 'box'])
        if self.values[9] == 'box':
            assert(self.values[0] == 1)  # use box geometry, wb is mandatory
        pass
        # output and input dirs
        os.path.isdir(self.values[10])
        os.path.isdir(self.values[11])

    def to_configure_prm(self):
        if_wb = self.values[0]
        type_of_bd = self.values[6]
        sp_rate = self.values[2]
        ov_age = self.values[3]
        potential_T = self.values[7]
        box_width = self.values[8]
        geometry = self.values[9]
        return if_wb, geometry, box_width, type_of_bd, potential_T, sp_rate, ov_age

    def to_configure_wb(self):
        '''
        Interface to configure_wb
        '''
        if_wb = self.values[0]
        geometry = self.values[9]
        potential_T = self.values[7]
        sp_age_trench = self.values[1]
        sp_rate = self.values[2]
        ov_age = self.values[3]
        ov_trans_age = self.values[4]
        ov_trans_length = self.values[5]
        if self.values[4] < 0.0:
            if_ov_trans = False
        else:
            if_ov_trans = True
        return if_wb, geometry, potential_T, sp_age_trench, sp_rate, ov_age,\
            if_ov_trans, ov_trans_age, ov_trans_length
    
    # todo
    def to_init():
        pass

    def wb_inputs_path():
        pass

    def prm_file_path():
        pass

    def wb_file_path():
        pass






class CASE(CasesP.CASE):
    '''
    class for a case
    More Attributes:
    '''
    def configure_prm(self, if_wb, geometry, box_width, type_of_bd, potential_T, sp_rate, ov_age):
        Ro = 6371e3
        if type_of_bd == "all free slip":  # boundary conditions
            if_fs_sides = True  # use free slip on both sides
        else:
            if_fs_sides = False
        o_dict = self.idict.copy()
        # Adiabatic surface temperature
        o_dict["Adiabatic surface temperature"] = str(potential_T)
        # geometry model
        if geometry == 'chunk':
            max_phi = box_width / Ro * 180.0 / np.pi  # extent in term of phi
            o_dict["Geometry model"] = prm_geometry_sph(max_phi)
        elif geometry == 'box':
            o_dict["Geometry model"] = prm_geometry_cart(box_width)
        # refinement
        if geometry == 'chunk':
            o_dict["Mesh refinement"]['Minimum refinement function'] = prm_minimum_refinement_sph()
        elif geometry == 'box':
            o_dict["Mesh refinement"]['Minimum refinement function'] = prm_minimum_refinement_cart()
        # boundary temperature model
        if geometry == 'chunk':
            o_dict['Boundary temperature model'] = prm_boundary_temperature_sph()
        elif geometry == 'box':
            o_dict['Boundary temperature model'] = prm_boundary_temperature_cart()
        # set up subsection reset viscosity function
        visco_plastic_twoD = self.idict['Material model']['Visco Plastic TwoD']
        if geometry == 'chunk':
            o_dict['Material model']['Visco Plastic TwoD'] =\
              prm_visco_plastic_TwoD_sph(visco_plastic_twoD, max_phi, if_fs_sides=if_fs_sides)
        elif geometry == 'box':
            o_dict['Material model']['Visco Plastic TwoD'] =\
              prm_visco_plastic_TwoD_cart(visco_plastic_twoD, box_width, if_fs_sides=if_fs_sides)
        # set up subsection Prescribed temperatures
        if type_of_bd == "all free slip":
            o_dict["Prescribe internal temperatures"] = "true"
            if geometry == 'chunk':
                o_dict['Prescribed temperatures'] =\
                    prm_prescribed_temperature_sph(max_phi, potential_T, sp_rate, ov_age)
            elif geometry == 'box':
                o_dict['Prescribed temperatures'] =\
                    prm_prescribed_temperature_cart(box_width, potential_T, sp_rate, ov_age)
        else:
            # remove this feature if otherwise
            pass
        self.idict = o_dict

    def configure_wb(self, if_wb, geometry, potential_T, sp_age_trench, sp_rate, ov_ag,\
        if_ov_trans, ov_trans_age, ov_trans_length):
        '''
        Configure world builder file
        Inputs:
            see description of CASE_OPT
        '''
        if not if_wb:
            # check first if we use wb file for this one
            return
        # potential T
        self.wb_dict['potential mantle temperature'] = potential_T
        # geometry
        if geometry == 'chunk':
            self.wb_dict["coordinate system"] = {"model": "spherical", "depth method": "begin segment"}
            self.wb_dict["cross section"] = [[0, 0], [180.0, 0.0]]
        elif geometry == 'box':
            self.wb_dict["coordinate system"] = {"model": "cartesian"}
            self.wb_dict["cross section"] = [[0, 0], [1e7, 0.0]]
        else:
            raise ValueError('%s: geometry must by one of \"chunk\" or \"box\"' % Utilities.func_name())
        # plates
        if geometry == 'chunk':
            Ro = float(self.idict['Geometry model']['Chunk']['Chunk outer radius'])
            self.wb_dict = wb_configure_plates(self.wb_dict, sp_age_trench,\
            sp_rate, ov_ag, Ro=Ro, if_ov_trans=if_ov_trans, ov_trans_age=ov_trans_age,\
            ov_trans_length=ov_trans_length, geometry=geometry) # plates
        elif geometry == 'box':
            Xmax = 7e6  # lateral extent of the box
            self.wb_dict = wb_configure_plates(self.wb_dict, sp_age_trench,\
            sp_rate, ov_ag, Xmax=Xmax, if_ov_trans=if_ov_trans, ov_trans_age=ov_trans_age,\
            ov_trans_length=ov_trans_length, geometry=geometry) # plates
        else:
            raise ValueError('%s: geometry must by one of \"chunk\" or \"box\"' % Utilities.func_name())


def wb_configure_plates(wb_dict, sp_age_trench, sp_rate, ov_age, **kwargs):
    '''
    configure plate in world builder
    '''
    Ro = kwargs.get('Ro', 6371e3)
    Xmax = kwargs.get('Xmax', 7e6)
    geometry = kwargs.get('geometry', 'chunk')
    o_dict = wb_dict.copy()
    trench_sph = (sp_age_trench * sp_rate / Ro) * 180.0 / np.pi
    trench_cart = sp_age_trench * sp_rate
    max_sph = 180.0  # maximum angle assigned to limit the extent of features.
    max_cart = 2 * Xmax
    side_angle = 5.0  # side angle to creat features in the 3rd dimension
    side_dist = 1e3
    if geometry == 'chunk':
        _side = side_angle
        _max = max_sph
        trench = trench_sph
    elif geometry == 'box':
        _side = side_dist
        _max = max_cart
        trench = trench_cart
        pass
    sp_ridge_coords = [[0, -_side], [0, _side]]
    # Overiding plate
    if_ov_trans = kwargs.get('if_ov_trans', False)  # transit to another age
    if if_ov_trans:
        i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate 1')
        ov_trans_feature, ov =\
            wb_configure_transit_ov_plates(wb_dict['features'][i0], trench,\
                ov_age, kwargs['ov_trans_age'], kwargs['ov_trans_length'],\
                Ro=Ro, geometry=geometry)
        o_dict['features'][i0] = ov_trans_feature
    else:
        # if no using transit plate, remove the feature
        try:
            i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate 1')
        except KeyError:
            pass
        o_dict = ParsePrm.RemoveWBFeatures(o_dict, i0)
        ov = trench
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate')
    op_dict = o_dict['features'][i0]
    op_dict["coordinates"] = [[ov, -_side], [ov, _side],\
        [_max, -_side], [_max, _side]] # trench position
    op_dict["temperature models"][0]["plate age"] = ov_age  # age of overiding plate
    o_dict['features'][i0] = op_dict
    # Subducting plate
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Subducting plate')
    sp_dict = o_dict['features'][i0]
    sp_dict["coordinates"] = [[0.0, -_side], [0.0, _side],\
        [trench, -_side], [trench, _side]] # trench position
    sp_dict["temperature models"][0]["spreading velocity"] = sp_rate
    sp_dict["temperature models"][0]["ridge coordinates"] = sp_ridge_coords
    o_dict['features'][i0] = sp_dict
    # Slab
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Slab')
    s_dict = o_dict['features'][i0]
    s_dict["coordinates"] = [[trench, -_side], [trench, _side]] 
    s_dict["dip point"] = [_max, 0.0]
    s_dict["temperature models"][0]["ridge coordinates"] = sp_ridge_coords
    s_dict["temperature models"][0]["plate velocity"] = sp_rate
    o_dict['features'][i0] = s_dict
    # mantle for substracting adiabat
    i0 = ParsePrm.FindWBFeatures(o_dict, 'mantle to substract')
    m_dict = o_dict['features'][i0]
    m_dict["coordinates"] =[[0.0, -_side], [0.0, _side],\
        [_max, -_side], [_max, _side]]
    o_dict['features'][i0] = m_dict
    return o_dict

def wb_configure_transit_ov_plates(i_feature, trench_sph, ov_age,\
    ov_trans_age, ov_trans_length, **kwargs):
    '''
    Transit overiding plate to a younger age at the trench
    See descriptions of the interface to_configure_wb
    '''
    geometry = kwargs.get('geometry', 'chunk')
    side_angle = 5.0  # side angle to creat features in the 3rd dimension
    Ro = kwargs.get("Ro", 6371e3)
    o_feature = i_feature.copy()
    trans_angle = ov_trans_length / Ro / np.pi * 180.0
    ov_sph = trench_sph  + trans_angle  # new ending point of the default overiding plage
    v = ov_trans_length / (ov_age - ov_trans_age)
    ridge_sph = trench_sph -\
        ov_trans_length * ov_trans_age / (ov_age - ov_trans_age) / Ro * 180.0 / np.pi
    o_feature["coordinates"] = [[trench_sph, side_angle], [trench_sph, -side_angle],\
        [ov_sph, side_angle], [ov_sph, -side_angle]]
    o_feature["temperature models"][0]["spreading velocity"] = v
    o_feature["temperature models"][0]["ridge coordinates"] =\
        [[ridge_sph, -side_angle], [ridge_sph, side_angle]]
    return o_feature, ov_sph


def prm_geometry_sph(max_phi):
    '''
    reset geometry for chunk geometry
    '''
    o_dict = {
        "Model name": "chunk",
        "Chunk": {
            "Chunk inner radius": "3.481e6",\
            "Chunk outer radius": "6.371e6",\
            "Chunk maximum longitude": "%.4e" % max_phi,\
            "Chunk minimum longitude": "0.0",\
            "Longitude repetitions": "2"
        }
    }
    return o_dict


def prm_minimum_refinement_sph(**kwargs):
    """
    minimum refinement function for spherical geometry
    """
    Ro = kwargs.get('Ro', 6371e3)
    o_dict = {
      "Coordinate system": "spherical",
      "Variable names": "r,phi,t",
      "Function constants": "Ro=%.4e, UM=670e3, DD=100e3" % Ro,
      "Function expression": "((Ro-r<UM)? \\\n                                   ((Ro-r<DD)? 8: 6): 0.0)"
    }
    return o_dict


def prm_minimum_refinement_cart(**kwargs):
    """
    minimum refinement function for cartesian geometry
    """
    Do = kwargs.get('Do', 2890e3)
    o_dict = {
      "Coordinate system": "cartesian",
      "Variable names": "x, y, t",
      "Function constants": "Do=%.4e, UM=670e3, DD=100e3" % Do,
      "Function expression": "((Do-y<UM)? \\\n                                   ((Do-y<DD)? 8: 6): 0.0)"
    }
    return o_dict


def prm_boundary_temperature_sph():
    '''
    boundary temperature model in spherical geometry
    '''
    o_dict = {
        "Fixed temperature boundary indicators": "bottom, top",
        "List of model names": "spherical constant",
        "Spherical constant": {
            "Inner temperature": "3500", 
            "Outer temperature": "273"
        }
    }
    return o_dict


def prm_boundary_temperature_cart():
    '''
    boundary temperature model in cartesian geometry
    '''
    o_dict = {
        "Fixed temperature boundary indicators": "bottom, top",
        "List of model names": "box",
        "Box": {
            "Bottom temperature": "3500",
            "Top temperature": "273"
            }
    }
    return o_dict


def prm_geometry_cart(box_width):
    '''
    reset geometry for box geometry
    '''
    o_dict = {
        "Model name": "box",
        "Box": {
            "X extent": "%.4e" % box_width,
            "Y extent": "2.8900e6",
            "X repetitions": "2"
        }
    }
    return o_dict


def prm_visco_plastic_TwoD_sph(visco_plastic_twoD, max_phi, **kwargs):
    '''
    reset subsection Visco Plastic TwoD
    Inputs:
        visco_plastic_twoD (dict): inputs for the "subsection Visco Plastic TwoD"
        part in a prm file
        kwargs(dict):
    '''
    o_dict = visco_plastic_twoD.copy()
    if_fs_sides = kwargs.get('if_fs_sides', True)
    if if_fs_sides:
        # use free slip on both sides, set ridges on both sides
        o_dict['Reset viscosity'] = 'true'
        o_dict['Reset viscosity function'] =\
            prm_reset_viscosity_function_sph(max_phi)
        o_dict["Reaction mor"] = 'true'
        o_dict["Reaction mor function"] =\
            prm_reaction_mor_function_sph(max_phi)
    else:
        # remove the related options
        pass
    return o_dict


def prm_visco_plastic_TwoD_cart(visco_plastic_twoD, box_width, **kwargs):
    '''
    reset subsection Visco Plastic TwoD
    Inputs:
        visco_plastic_twoD (dict): inputs for the "subsection Visco Plastic TwoD"
        part in a prm file
        kwargs(dict):
    '''
    o_dict = visco_plastic_twoD.copy()
    if_fs_sides = kwargs.get('if_fs_sides', True)
    if if_fs_sides:
        # use free slip on both sides, set ridges on both sides
        o_dict['Reset viscosity'] = 'true'
        o_dict['Reset viscosity function'] =\
            prm_reset_viscosity_function_cart(box_width)
        o_dict["Reaction mor"] = 'true'
        o_dict["Reaction mor function"] =\
            prm_reaction_mor_function_cart(box_width)
    else:
        # remove the related options
        pass
    return o_dict


def prm_reset_viscosity_function_sph(max_phi):
    '''
    Default setting for Reset viscosity function in spherical geometry
    '''
    max_phi_in_rad = max_phi * np.pi / 180.0
    odict = {
        "Coordinate system": "spherical",
        "Variable names": "r, phi",
        "Function constants": "Depth=1.45e5, Width=2.75e5, Ro=6.371e6, PHIM=%.4e, CV=1e20" % max_phi_in_rad,
        "Function expression": "(((r > Ro - Depth) && ((Ro*phi < Width) || (Ro*(PHIM-phi) < Width)))? CV: -1.0)"
      }
    return odict


def prm_reset_viscosity_function_cart(box_width):
    '''
    Default setting for Reset viscosity function in cartesian geometry
    '''
    odict = {
        "Coordinate system": "cartesian",
        "Variable names": "x, y",
        "Function constants": "Depth=1.45e5, Width=2.75e5, Do=2.890e6, xm=%.4e, CV=1e20" % box_width,
        "Function expression": "(((y > Do - Depth) && ((x < Width) || (xm-x < Width)))? CV: -1.0)"
    }
    return odict


def prm_reaction_mor_function_sph(max_phi):
    '''
    Default setting for Reaction mor function in cartesian geometry
    '''
    max_phi_in_rad = max_phi * np.pi / 180.0
    odict = {
        "Coordinate system": "spherical",\
        "Variable names": "r, phi",\
        "Function constants": "Width=2.75e5, Ro=6.371e6, PHIM=%.4e, DCS=7.500e+03, DHS=3.520e+04" % max_phi_in_rad,\
        "Function expression": "((r > Ro - DCS) && (Ro*phi < Width)) ? 0:\\\n                                        ((r < Ro - DCS) && (r > Ro - DHS) && (Ro*phi < Width)) ? 1:\\\n                                        ((r > Ro - DCS) && (Ro*(PHIM - phi) < Width)) ? 2:\\\n                                        ((r < Ro - DCS) && (r > Ro - DHS) && (Ro*(PHIM - phi) < Width)) ? 3: -1"
      }
    return odict


def prm_reaction_mor_function_cart(box_width):
    '''
    Default setting for Reaction mor function in cartesian geometry
    '''
    odict = {
        "Coordinate system": "cartesian",
        "Variable names": "x, y",
        "Function constants": "Width=2.75e5, Do=2.890e6, xm=%.4e, DCS=7.500e+03, DHS=3.520e+04" % box_width,
        "Function expression": "((y > Do - DCS) && (x < Width)) ? 0:\\\n                                        ((y < Do - DCS) && (y > Do - DHS) && (x < Width)) ? 1:\\\n                                        ((y > Do - DCS) && (xm - x < Width)) ? 2:\\\n                                        ((y < Do - DCS) && (y > Do - DHS) && (xm - x < Width)) ? 3: -1"
    }
    return odict


def prm_prescribed_temperature_sph(max_phi, potential_T, sp_rate, ov_age):
    '''
    Default setting for Prescribed temperatures in spherical geometry
    '''
    max_phi_in_rad = max_phi * np.pi / 180.0
    odict = {
        "Indicator function": {
          "Coordinate system": "spherical",\
          "Variable names": "r, phi",\
          "Function constants": "Depth=1.45e5, Width=2.75e5, Ro=6.371e6, PHIM=%.4e" % max_phi_in_rad,\
          "Function expression": "(((r>Ro-Depth)&&((r*phi<Width)||(r*(PHIM-phi)<Width))) ? 1:0)"\
        },\
        "Temperature function": {
          "Coordinate system": "spherical",\
          "Variable names": "r, phi",\
          "Function constants": "Depth=1.45e5, Width=2.75e5, Ro=6.371e6, PHIM=%.4e,\\\n                             AGEOP=%.4e, TS=2.730e+02, TM=%.4e, K=1.000e-06, VSUB=%.4e, PHILIM=1e-6" %\
              (max_phi_in_rad, ov_age * year, potential_T, sp_rate / year),\
          "Function expression": "((r*(PHIM-phi)<Width) ? TS+(TM-TS)*(1-erfc(abs(Ro-r)/(2*sqrt(K*AGEOP)))):\\\n\t(phi > PHILIM)? (TS+(TM-TS)*(1-erfc(abs(Ro-r)/(2*sqrt((K*Ro*phi)/VSUB))))): TM)"
        }
    }
    return odict


def prm_prescribed_temperature_cart(box_width, potential_T, sp_rate, ov_age):
    '''
    Default setting for Prescribed temperatures in cartesian geometry
    '''
    odict = {
        "Indicator function": {
          "Coordinate system": "cartesian",
          "Variable names": "x, y",
          "Function constants": "Depth=1.45e5, Width=2.75e5, Do=2.890e6, xm=%.4e" % box_width,
          "Function expression": "(((y>Do-Depth)&&((x<Width)||(xm-x<Width))) ? 1:0)"
        },
        "Temperature function": {
          "Coordinate system": "cartesian",
          "Variable names": "x, y",
          "Function constants": "Depth=1.45e5, Width=2.75e5, Do=2.890e6, xm=%.4e,\\\n                             AGEOP=%.4e, TS=2.730e+02, TM=%.4e, K=1.000e-06, VSUB=%.4e, XLIM=6" %\
                (box_width, ov_age * year, potential_T, sp_rate / year),\
          "Function expression": "(xm-x<Width) ? TS+(TM-TS)*(1-erfc(abs(Do-y)/(2*sqrt(K*AGEOP)))):\\\n\t((x > XLIM)? (TS+(TM-TS)*(1-erfc(abs(Do-y)/(2*sqrt((K*x)/VSUB))))): TM)"
        }
    }
    return odict


def create_case(json_file):
    '''
    A wrapper for the CASES class
    Inputs:
        json_file(str): path of a json file
        todo
    '''
    assert(os.access(json_file, os.R_OK))
    Case_Opt = CASE_OPT()
    Case_Opt.read_json(json_file)
    Case_Opt.check()
    # Case = CASE('wb_setup', prm_file, wb_inputs=wb_file)
    Case = CASE(Case_Opt.to_init(), wb_inputs=Case_Opt.wb_inputs_path())
    Case.configure_prm(*Case_Opt.to_configure_prm())
    Case.configure_wb(*Case_Opt.to_configure_wb())
    # create new case
    Case.create(Case_Opt.o_dir())
    assert(os.path.isfile(Case_Opt.prm_file_path())) # assert files generated
    assert(os.path.isfile(Case_Opt.wb_file_path()))




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