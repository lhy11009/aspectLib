# -*- coding: utf-8 -*-
r"""(one line description)

This exports: 

  -

This depends on:

  -  

Examples of usage:

  - convert parameters of some type of rheology to aspect input:

        python -m shilofue.Rheology convert_to_ASPECT -r HK03

        options for "-r":
            HK03 - Hirth & Kohlstedt 2003
            AB17 - Arredondo & Billen 2017

descriptions
""" 

import json
import os
import sys
import math
import argparse
import numpy as np
from matplotlib import pyplot as plt

R = 8.314

ASPECT_LAB_DIR = os.environ['ASPECT_LAB_DIR']

class CheckValueError(Exception):
    pass


class RHEOLOGY_PRM():
    """
    class for rheologies
    """
    def __init__(self):
        '''
        Initiation, initiate rheology parameters
        '''
        # dislocation creep in Hirth & Kohlstedt 2003
        self.HK03_disl = \
            {
                "A": 90,
                "p": 0.0,
                "r": 1.2,
                "n": 3.5,
                "E": 480e3,
                "V": 11e-6,
                "d" : 1e4,
                "Coh" : 1000.0
            }

        # diffusion creep in Hirth & Kohlstedt 2003
        self.HK03_diff = \
            {
                "A" : 1.0e6,
                "p" : 3.0,
                "r" : 1.0,
                "n" : 1.0,
                "E" : 335e3,
                "V" : 4e-6,
                "d" : 1e4,
                "Coh" : 1000.0
            }
        
        # dislocation creep in Arredondo & Billen 2017
        self.AB17_disl = \
            {
                "A": 2.57e-20,
                "p": 0.0,
                "r": 1.2,
                "n": 3.5,
                "E": 496e3,
                "V": 11e-6,
                "d" : 1e4,
                "Coh" : 1000.0
            }
        
        # diffusion creep in Arredondo & Billen 2017
        self.AB17_diff = \
            {
                "A" : 2.85e6,
                "p" : 3.0,
                "r" : 1.0,
                "n" : 1.0,
                "E" : 317e3,
                "V" : 4e-6,
                "d" : 1e4,
                "Coh" : 1000.0
            }


def Config(_kwargs, _name, _default):
    """
    def Config(_kwargs, _name, _default)

    Read variable value and assign default if not found
    """
    try:
        value = _kwargs[_name]
    except KeyError:
        value = _default
    return value

def CreepStress(creep_type, strain_rate, P, T, d, Coh):
    """
    def DislocationCreep(strain_rate, P, T, d, Coh)

    Calculate stress by flow law in form of (strain_rate / B)^(1.0 / n) * exp((E + P * V) / (n * R * T))
    Units:
     - P: Pa
     - T: K
     - d: mu m
     - Coh: H / 10^6 Si
     - Return value: Mpa
    Pay attention to pass in the right value, this custom is inherited
    """
    A = creep_type['A']
    p = creep_type['p']
    r = creep_type['r']
    n = creep_type['n']
    E = creep_type['E']
    V = creep_type['V']
    # calculate B
    B = A * d**(-p) * Coh**r
    return (strain_rate / B)**(1.0 / n) * np.exp((E + P * V) / (n * R * T))


def CreepRheology(creep_type, strain_rate, P, T, d=None, Coh=None):
    """
    def CreepRheology(creep_type, strain_rate, P, T, d, Coh):

    Calculate viscosity by flow law in form of (strain_rate)**(1.0 / n - 1) * (B)**(-1.0 / n) * np.exp((E + P * V) / (n * R * T))
    Units:
     - P: Pa
     - T: K
     - d: mu m
     - Coh: H / 10^6 Si
     - Return value: Pa*s
    Pay attention to pass in the right value, this custom is inherited
    """
    A = creep_type['A']
    p = creep_type['p']
    r = creep_type['r']
    n = creep_type['n']
    E = creep_type['E']
    V = creep_type['V']
    if d is None:
        d = creep_type['d']
    if Coh is None:
        Coh = creep_type['Coh']
    # calculate B
    B = A * d**(-p) * Coh**r
    return (strain_rate)**(1.0 / n - 1) * (B)**(-1.0 / n) * np.exp((E + P * V) / (n * R * T)) * 1e6


def CreepRheologyInAspectViscoPlastic(creep_type, strain_rate, P, T):
    """
    def CreepRheologyInAspectVisoPlastic(creep_type, strain_rate, P, T)

    Calculate viscosity by way of Visco Plastic module in aspect
    flow law in form of 0.5 * A**(-1.0 / n) * d**(m / n) * (strain_rate)**(1.0 / n - 1) * np.exp((E + P * V) / (n * R * T))
    Units:
     - P: Pa
     - T: K
     - d: m
     - Return value: Pa*s
    """
    A = creep_type['A']
    m = creep_type['m']
    n = creep_type['n']
    E = creep_type['E']
    V = creep_type['V']
    d = creep_type['d']
    # calculate B
    return 0.5 * A**(-1.0 / n) * d**(m / n) * (strain_rate)**(1.0 / n - 1) * np.exp((E + P * V) / (n * R * T))


def Convert2AspectInput(creep_type):
    """
    Viscosity is calculated by flow law in form of (strain_rate)**(1.0 / n - 1) * (B)**(-1.0 / n) * np.exp((E + P * V) / (n * R * T)) * 1e6
    while in aspect, flow law in form of 0.5 * A**(-1.0 / n) * d**(m / n) * (strain_rate)**(1.0 / n - 1) * np.exp((E + P * V) / (n * R * T))
    Original Units:
     - P: Pa
     - T: K
     - d: mm
     - Coh: H / 10^6 Si
    Original Units:
     - P: Pa
     - T: K
     - d: m
    """
    # read in initial value
    A = creep_type['A']
    p = creep_type['p']
    r = creep_type['r']
    n = creep_type['n']
    E = creep_type['E']
    V = creep_type['V']
    d = creep_type['d']
    Coh = creep_type['Coh']
    # prepare values for aspect
    aspect_creep_type = {}
    aspect_creep_type['A'] = 1e6**(-p) * (2e6)**(-n) * Coh**r * A
    aspect_creep_type['d'] = d / 1e6
    aspect_creep_type['n'] = n
    aspect_creep_type['m'] = p
    aspect_creep_type['E'] = E
    aspect_creep_type['V'] = V
    return aspect_creep_type


def GetLowerMantleRheology(upper_mantle_creep_method, jump, T, P, **kwargs):
    """
    get flow law parameters in lower mantle based on upper mantle viscosity and jump in viscosity
    variables:
     - jump: viscosity jump at 660km
     - T: temperature at 660km
     - P: pressure at 660km
    """
    # extra inputs
    strategy = Config(kwargs, 'strategy', 'A')
    V1 = Config(kwargs, 'V1', upper_mantle_creep_method['V'])
    # read upper mantle values
    A = upper_mantle_creep_method['A']
    m = upper_mantle_creep_method['m']
    n = upper_mantle_creep_method['n']
    E = upper_mantle_creep_method['E']
    V = upper_mantle_creep_method['V']
    d = upper_mantle_creep_method['d']

    lower_mantle_creep_method = dict(upper_mantle_creep_method)
    lower_mantle_creep_method['V'] = V1
    if strategy is 'A':
        lower_mantle_creep_method['A'] = jump**(-n) * A * math.exp(P * (V1 - V) / (R * T))
    else:
        lower_mantle_creep_method['d'] = jump**(n / m) * d * math.exp(P * (V-V1) / (m * R * T))
    return lower_mantle_creep_method


def CheckConvert2AspectInput():
    """
    def CheckConvert2AspectInput(filename='Hirth_Kohlstedt.json')

    Check with parameters an example from Hirth and Kohlstedt, 2013. 
    Code converts to aspect input form, then calculate viscosity by the two forms and check to make sure everything is right
    Tolerence set to be 1%
    """
    
    # check file exist
    filename = os.path.join(ASPECT_LAB_DIR, "files", "Hirth_Kohlstedt.json")
    # read in standard flow law parameters
    check_result = [0.0, 0.0]
    creep_types = json.load(open(filename, 'r'))
    diffusion_creep = creep_types['diffusion_creep']
    dislocation_creep = creep_types['dislocation_creep']
    # calculate viscosity by standard form
    check_result[0] = CreepRheology(diffusion_creep, 1e-15, 10e9, 1300 + 273.15)
    check_result[1] = CreepRheology(dislocation_creep, 1e-15, 10e9, 1300 + 273.15)
    # convert to aspect inputs
    diffusion_creep_aspect = Convert2AspectInput(diffusion_creep)
    dislocation_creep_aspect = Convert2AspectInput(dislocation_creep)
    # check for viscosity
    check0 = CreepRheologyInAspectViscoPlastic(diffusion_creep_aspect, 1e-15, 10e9, 1300 + 273.15)
    if abs((check0 - check_result[0]) / check_result[0]) > 0.01:
        raise CheckValueError("%.4e vs %.4e check value for converting to aspect parametrs and yielding the same viscosity failed" % (check0, check_result[0]))
    check1 = CreepRheologyInAspectViscoPlastic(dislocation_creep_aspect, 1e-15, 10e9, 1300 + 273.15)
    if abs((check1 - check_result[1]) / check_result[1]) > 0.01:
        raise CheckValueError("%.4e vs %.4e check value for converting to aspect parametrs and yielding the same viscosity failed" % (check1, check_result[1]))
    print("check value for for converting to aspect parametrs and yielding the same viscosity: succeed")


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
    parser.add_argument('-r', '--rheology', type=str,
                        default='HK03',
                        help='Type of rheology to use')
    _options = []
    try:
        _options = sys.argv[2: ]
    except IndexError:
        pass
    arg = parser.parse_args(_options)

    # commands
    if _commend == 'convert_to_ASPECT':
        rheology = arg.rheology
        # read in standard flow law parameters
        RheologyPrm = RHEOLOGY_PRM()
        diffusion_creep = getattr(RheologyPrm, rheology + "_diff")
        dislocation_creep = getattr(RheologyPrm, rheology + "_disl")
        # convert 2 aspect
        diffusion_creep_aspect = Convert2AspectInput(diffusion_creep)
        dislocation_creep_aspect = Convert2AspectInput(dislocation_creep)
        # screen output
        print("ASPECT diffusion creep: ")
        print(diffusion_creep_aspect)
        print("ASPECT dislocation creep: ")
        print(dislocation_creep_aspect)


# run script
if __name__ == '__main__':
    main()
