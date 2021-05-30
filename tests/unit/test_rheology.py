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
from shilofue.flow_law_functions import visc_diff_HK

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
    RheologyPrm = RHEOLOGY_PRM()
    diffusion_creep = RheologyPrm.HK03_diff
    dislocation_creep = RheologyPrm.HK03_disl
    
    # check for stress
    tolerance = 0.01
    check0 = CreepStress(diffusion_creep, 7.8e-15, 1e9, 1400 + 273.15, 1e4, 1000.0)
    assert(abs(check0 - 0.2991) / 0.2991 < tolerance)
    
    check1 = CreepStress(dislocation_creep, 2.5e-12, 1e9, 1400 + 273.15, 1e4, 1000.0)
    assert(abs((check1 - 0.3006) / 0.3006) < tolerance)
    
    # check for viscosity
    check2 = CreepRheology(diffusion_creep, 7.8e-15, 1e9, 1400 + 273.15, 1e4, 1000.0)
    assert(abs((check2 - 1.917382525013947e+19) / 1.917382525013947e+19) < tolerance)
    
    check2_mg = visc_diff_HK(1400 + 273.15,1e9,1e4,1000.0,'con','orig','mid','mid')
    
    check3 = CreepRheology(dislocation_creep, 2.5e-12, 1e9, 1400 + 273.15, 1e4, 1000.0)
    assert(abs((check3 - 6.0119447368476904e+16) / check3) < tolerance)

    # compute A
    check4 = CreepComputeA(diffusion_creep, 7.8e-15, 1e9, 1400 + 273.15, 3.8348e19, 1e4, 1000.0)
    assert(abs((check4 - 999990.8861030287) / check4) < tolerance)
    
    check5 = CreepComputeA(dislocation_creep, 2.5e-12, 1e9, 1400 + 273.15, 1.2024e17, 1e4, 1000.0)
    assert(abs((check5 - 89.99710450882286) / check5) < tolerance)


def test_HirthKohlstedt_wet_modified():
    """
    test the implementation of equations from Hirth & Kohlstedt, 2003(filename='Hirth_Kohlstedt.json')
    with a modified rheology, see magali's 
    Tolerence set to be 1%
    Asserts:
        My implementation yield same result with Magali's
    """
    # get rheology
    rheology = 'HK03_wet_mod'
    diffusion_creep, dislocation_creep = GetRheology(rheology)

    tolerance = 0.01

    check0 = CreepRheology(diffusion_creep, 7.8e-15, 1e9, 1400 + 273.15, 1e4, 1000.0)
    check0_std = visc_diff_HK(1400 + 273.15,1e9,1e4,1000.0,'wet','new','mid','mid')
    assert(abs(check0 - check0_std) / check0 < tolerance)
    pass


def test_creep_compute_V():
    # test CreepComputeV
    tolerance = 1e-6
    rheology = 'HK03_wet_mod'
    diffusion_creep, dislocation_creep = GetRheology(rheology)
    # test 1: diffusion creep
    eta = CreepRheology(diffusion_creep, 7.8e-15, 1e9, 1400 + 273.15, 1e4, 1000.0)
    V = CreepComputeV(diffusion_creep, 7.8e-15, 1e9, 1400 + 273.15, eta)
    assert(abs((V - diffusion_creep['V'])/V) < tolerance)

    # test 2: dislocation creep
    eta = CreepRheology(dislocation_creep, 7.8e-15, 1e9, 1400 + 273.15, 1e4, 1000.0, use_effective_strain_rate=1)
    V = CreepComputeV(dislocation_creep, 7.8e-15, 1e9, 1400 + 273.15, eta, use_effective_strain_rate=1)
    assert(abs((V - dislocation_creep['V'])/V) < tolerance)



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
    # test 1: convert without computing the strain tensor invarant
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
    

def test_Convert2AspectInput1():
    """
    check the implementation of Convert2AspectInput(filename='Hirth_Kohlstedt.json')
    with parameters and example from Hirth and Kohlstedt, 2013. 
    """
    tolerance = 0.01
    
    check_result = [0.0, 0.0]
    # read in standard flow law parameters
    RheologyPrm = RHEOLOGY_PRM()
    diffusion_creep, dislocation_creep = GetRheology("HK03_wet_mod1")
    # test 2: convert with computing the strain tensor invarant
    # calculate viscosity by standard form
    check_result[0] = CreepRheology(diffusion_creep, 1e-15, 10e9, 1300 + 273.15)
    check_result[1] = CreepRheology(dislocation_creep, 1e-15, 10e9, 1300 + 273.15, use_effective_strain_rate=True)
    # convert to aspect inputs
    diffusion_creep_aspect = Convert2AspectInput(diffusion_creep)
    dislocation_creep_aspect = Convert2AspectInput(dislocation_creep, use_effective_strain_rate=True)
    # check for viscosity
    check0 = CreepRheologyInAspectViscoPlastic(diffusion_creep_aspect, 1e-15, 10e9, 1300 + 273.15)
    assert (abs((check0 - check_result[0]) / check_result[0]) < tolerance)

    check1 = CreepRheologyInAspectViscoPlastic(dislocation_creep_aspect, 1e-15, 10e9, 1300 + 273.15)
    assert(abs((check1 - check_result[1]) / check_result[1]) < tolerance)


def test_Convert2AspectInputLowerMantle():
    """
    check the implementation of GetLowerMantleRheology
    """
    tolerance = 1e-6
    diffusion_creep = {
        "A": 8.17868945868946e-17,
        "d": 0.01,
        "n": 1.0,
        "m": 3.0,
        "E": 300000.0,
        "V": 6.899999999999998e-06}
    P660 = 20e9
    T660 = 2000
    V1 = 4e-6
    jump = 30.0
    diff_lm = GetLowerMantleRheology(diffusion_creep, jump, T660, P660, V1=V1, strategy='composite', eta=1e21)
    eta_lower = CreepRheologyInAspectViscoPlastic(diff_lm, 1e-15, P660, T660)

    assert(abs(eta_lower - jump * 1e21)/1e21 < tolerance)

    
def test_Convert2AspectInput():
    """
    compare different rheologies
    """
    pass


# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

