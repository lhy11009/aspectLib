# -*- coding: utf-8 -*-
r"""Test for Rheology.py

This outputs:


This depends on:


Examples of usage:

  - default usage:

        python -m pytest test_rheology.py

descriptions:
    every function is a separate test, combined usage with pytest module
""" 


import os
# import pytest
# import filecmp  # for compare file contents
# import numpy as np
# import shilofue.Foo as Foo  # import test module
# from shilofue.Utilities import 
# from matplotlib import pyplot as plt
# from shutil import rmtree  # for remove directories
from shilofue.Rheology import *

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'parse')

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)


def test_HirthKohlstedt():
    """
    test the implementation of equations from Hirth & Kohlstedt, 2003(filename='Hirth_Kohlstedt.json')
    Tolerence set to be 1%
    Asserts:
        values computed from the rheology
    """
    # read parameters
    check_result = [.2991, .3006, 3.8348e19, 1.2024e17]
    RheologyPrm = RHEOLOGY_PRM()
    diffusion_creep = RheologyPrm.HK03_diff
    dislocation_creep = RheologyPrm.HK03_disl

    # check for stress
    tolerance = 0.01
    check0 = CreepStress(diffusion_creep, 7.8e-15, 1e9, 1400 + 273.15, 1e4, 1000.0)
    assert(abs(check0 - check_result[0]) / check_result[0] < tolerance)
    
    check1 = CreepStress(dislocation_creep, 2.5e-12, 1e9, 1400 + 273.15, 1e4, 1000.0)
    assert(abs((check1 - check_result[1]) / check_result[1]) < tolerance)
    
    # check for viscosity
    check2 = CreepRheology(diffusion_creep, 7.8e-15, 1e9, 1400 + 273.15, 1e4, 1000.0)
    assert(abs((check2 - check_result[2]) / check_result[2]) < tolerance)
    
    check3 = CreepRheology(dislocation_creep, 2.5e-12, 1e9, 1400 + 273.15, 1e4, 1000.0)
    assert(abs((check3 - check_result[3]) / check_result[3]) < tolerance)


def test_AspectExample():
    """
    check the implementation of aspcet rheology

    Tolerence set to be 1%
    """
    tolerance = 0.01
    check_result = [9.1049e+20]
    # these are two examples from aspect
    # I cannot track where this is from, as it's from a previous
    #version of codes, and I didn't document them properly
    diffusion_creep = \
        {
            "A": 8.571e-16,
            "d": 8.20e-3,
            "m": 3.0,
            "n": 1.0,
            "E": 335.0e3,
            "V": 4.0e-6
        }
    dislocation_creep = \
       {
           "A": 6.859e-15,
           "d": 8.20e-3,
           "m": 1.0,
           "n": 3.5,
           "E": 480.0e3,
           "V": 11.0e-6
       }
    check0 = CreepRheologyInAspectViscoPlastic(diffusion_creep, 1e-15, 10e9, 1300 + 273.15)
    assert(abs(check0 - check_result[0]) / check_result[0] < tolerance)


def test_Convert2AspectInput():
    """
    check the implementation of Convert2AspectInput(filename='Hirth_Kohlstedt.json')
    with parameters and example from Hirth and Kohlstedt, 2013. 
    Tolerence set to be 1%
    """
    tolerance = 0.01
    
    check_result = [0.0, 0.0]
    # read in standard flow law parameters
    RheologyPrm = RHEOLOGY_PRM()
    diffusion_creep = RheologyPrm.HK03_diff
    dislocation_creep = RheologyPrm.HK03_disl
    # calculate viscosity by standard form
    check_result[0] = CreepRheology(diffusion_creep, 1e-15, 10e9, 1300 + 273.15)
    check_result[1] = CreepRheology(dislocation_creep, 1e-15, 10e9, 1300 + 273.15)
    # convert to aspect inputs
    diffusion_creep_aspect = Convert2AspectInput(diffusion_creep)
    dislocation_creep_aspect = Convert2AspectInput(dislocation_creep)
    # check for viscosity
    check0 = CreepRheologyInAspectViscoPlastic(diffusion_creep_aspect, 1e-15, 10e9, 1300 + 273.15)
    assert (abs((check0 - check_result[0]) / check_result[0]) < tolerance)

    check1 = CreepRheologyInAspectViscoPlastic(dislocation_creep_aspect, 1e-15, 10e9, 1300 + 273.15)
    assert(abs((check1 - check_result[1]) / check_result[1]) < tolerance)
    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

