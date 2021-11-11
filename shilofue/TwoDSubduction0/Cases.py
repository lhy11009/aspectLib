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
import shilofue.ParsePrm as ParsePrm

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
            1: sp_age_trench (float) - Age of the subducting plate at trench (yr)
            2: sp_rate - Spreading rate of the subducting plate (m/ yr)
            3: ov_age (float) - Age of the overiding plate (yr)
            4: ov_trans_age (float) - Age of the transit overiding plate (yr)
                this is assigned to -1.0 by default, which means we don't use
                transit plate
            5: ov_trans_length (float) - Length of the transit overiding plate (m)
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
        ov_trans_age = self.values[4]
        ov_trans_length = self.values[5]
        if self.values[4] < 0.0:
            if_ov_trans = False
        else:
            if_ov_trans = True

        return if_wb, sp_age_trench, sp_rate, ov_age,\
            if_ov_trans, ov_trans_age, ov_trans_length


class CASE(CasesP.CASE):
    '''
    class for a case
    More Attributes:
    '''
    def configure_wb(self, if_wb, sp_age_trench, sp_rate, ov_ag,\
        if_ov_trans, ov_trans_age, ov_trans_length):
        '''
        Configure world builder file
        Inputs:
            see description of CASE_OPT
        '''
        Ro = float(self.idict['Geometry model']['Chunk']['Chunk outer radius'])
        if not if_wb:
            # check first if we use wb file for this one
            return
        self.wb_dict = wb_configure_plates_sph(self.wb_dict, sp_age_trench,\
            sp_rate, ov_ag, Ro=Ro, if_ov_trans=if_ov_trans, ov_trans_age=ov_trans_age,\
            ov_trans_length=ov_trans_length) # plates
        pass


def wb_configure_plates_sph(wb_dict, sp_age_trench, sp_rate, ov_age, **kwargs):
    '''
    configure plate in world builder
    '''
    Ro = kwargs.get('Ro', 6371e3)
    o_dict = wb_dict.copy()
    trench_sph = (sp_age_trench * sp_rate / Ro) * 180.0 / np.pi
    sp_ridge_coords = [[0, -5], [0, -5]]
    max_sph = 180.0  # maximum angle assigned to limit the extent of features.
    side_angle = 5.0  # side angle to creat features in the 3rd dimension
    # Overiding plate
    if_ov_trans = kwargs.get('if_ov_trans', False)  # transit to another age
    if if_ov_trans:
        i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate 1')
        ov_trans_feature, ov_sph =\
            wb_configure_transit_ov_plates_sph(wb_dict['features'][i0], trench_sph,\
                ov_age, kwargs['ov_trans_age'], kwargs['ov_trans_length'], Ro=Ro)
        o_dict['features'][i0] = ov_trans_feature
    else:
        # if no using transit plate, remove the feature
        try:
            i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate 1')
        except KeyError:
            pass
        o_dict = ParsePrm.RemoveWBFeatures(o_dict, i0)
        ov_sph = trench_sph
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Overiding plate')
    op_dict = o_dict['features'][i0]
    op_dict["coordinates"] = [[ov_sph, -side_angle], [ov_sph, side_angle],\
        [max_sph, -side_angle], [max_sph, side_angle]] # trench position
    op_dict["temperature models"][0]["plate age"] = ov_age  # age of overiding plate
    o_dict['features'][i0] = op_dict
    # Subducting plate
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Subducting plate')
    sp_dict = o_dict['features'][i0]
    sp_dict["coordinates"] = [[0.0, -side_angle], [0.0, side_angle],\
        [trench_sph, -side_angle], [trench_sph, side_angle]] # trench position
    sp_dict["temperature models"][0]["spreading velocity"] = sp_rate
    sp_dict["temperature models"][0]["ridge coordinates"] = sp_ridge_coords
    o_dict['features'][i0] = sp_dict
    # Slab
    i0 = ParsePrm.FindWBFeatures(o_dict, 'Slab')
    s_dict = o_dict['features'][i0]
    s_dict["coordinates"] = [[trench_sph, -side_angle], [trench_sph, side_angle]] 
    s_dict["dip point"] = [max_sph, 0.0]
    s_dict["temperature models"][0]["ridge coordinates"] = sp_ridge_coords
    s_dict["temperature models"][0]["plate velocity"] = sp_rate
    o_dict['features'][i0] = s_dict
    # mantle for substracting adiabat
    i0 = ParsePrm.FindWBFeatures(o_dict, 'mantle to substract')
    m_dict = o_dict['features'][i0]
    m_dict["coordinates"] =[[0.0, -side_angle], [0.0, side_angle],\
        [max_sph, -side_angle], [max_sph, side_angle]]
    o_dict['features'][i0] = m_dict
    return o_dict

def wb_configure_transit_ov_plates_sph(i_feature, trench_sph, ov_age,\
    ov_trans_age, ov_trans_length, **kwargs):
    '''
    Transit overiding plate to a younger age at the trench
    See descriptions of the interface to_configure_wb
    '''
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