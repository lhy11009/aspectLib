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
import pytest
# import filecmp  # for compare file contents
# import numpy as np
# import shilofue.Foo as Foo  # import test module
# from shilofue.Utilities import 
# from matplotlib import pyplot as plt
# from shutil import rmtree  # for remove directories
from shilofue.Rheology import *
from shilofue.FlowLaws import visc_diff_HK, visc_disl_HK, peierls_visc_from_stress

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'parse')

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)


def test_RefitRheology():
    '''
    test the RefitRheology function
    '''
    # get the original rheology
    rheology = "WarrenHansen23"
    rheology_prm_dict = RHEOLOGY_PRM()
    diffusion_creep_ori = getattr(rheology_prm_dict, rheology + "_diff")
    dislocation_creep_ori = getattr(rheology_prm_dict, rheology + "_disl")
    rheology_dict = {'diffusion': diffusion_creep_ori, 'dislocation': dislocation_creep_ori}

    # prescribe the correction
    diff_correction = {'A': 1.0, 'p': 0.0, 'r': 0.0, 'n': 0.0, 'E': 0.0, 'V': -2.1e-6}
    disl_correction = {'A': 1.0, 'p': 0.0, 'r': 0.0, 'n': 0.0, 'E': 0.0, 'V': 3e-6}

    # prescribe the reference state
    ref_state = {}
    ref_state["Coh"] = 500.0 # H / 10^6 Si
    ref_state["stress"] = 50.0 # MPa
    ref_state["P"] = 100.0e6 # Pa
    ref_state["T"] = 1250.0 + 273.15 # K
    ref_state["d"] = 15.0 # mu m

    # refit rheology
    rheology_dict_refit = RefitRheology(rheology_dict, diff_correction, disl_correction, ref_state)

    # check
    diffusion_creep = rheology_dict_refit['diffusion']
    dislocation_creep = rheology_dict_refit['dislocation']
    diff_A_std = 285230.546793
    assert(abs(diffusion_creep['A'] - diff_A_std)/diff_A_std < 1e-6)
    diff_V_std = 1.9e-06
    assert(abs(diffusion_creep['V'] - diff_V_std)/diff_V_std < 1e-6)
    disl_A_std = 20.479460
    assert(abs(dislocation_creep['A'] - disl_A_std)/disl_A_std < 1e-6)
    disl_V_std = 1.4e-05
    assert(abs(dislocation_creep['V'] - disl_V_std)/disl_V_std < 1e-6)


def test_debug_Idrissy_convergence_large_P():
    '''
    test used to debug the Idrissy flow law with a very big pressure
    '''
    strain_rate_std = 1.0e-05
    stress = 2.516e5 # Pa
    T = 3500.0 # K
    P = 106608753761.386 # not dependent on P (V = 0)
    eta, strain_rate = peierls_visc_from_stress("Idrissi16", P, T, stress)
    eta_1, strain_rate1  = peierls_visc_from_stress("Idrissi16", P, T, 1.01*stress)
    strain_rate_1 = stress / 2.0 / eta_1
    strain_rate_dev = (strain_rate_1 - strain_rate) / (0.01*stress)
    print("T = ", T)
    print("strain_rate = ", strain_rate)  # debug
    print("strain_rate_dev = ", strain_rate_dev)
    T = 1600.0 # K
    eta, strain_rate = peierls_visc_from_stress("Idrissi16", P, T, stress)
    eta_1, strain_rate1  = peierls_visc_from_stress("Idrissi16", P, T, 1.01*stress)
    strain_rate_1 = stress / 2.0 / eta_1
    strain_rate_dev = (strain_rate_1 - strain_rate) / (0.01*stress)
    print("T = ", T)
    print("strain_rate = ", strain_rate)  # debug
    print("strain_rate_dev = ", strain_rate_dev)


def test_CreepStress_CreepRheology_CreepStrainRate_ComputeComposite():
    """
    test the implementation of equations from Hirth & Kohlstedt, 2003(filename='Hirth_Kohlstedt.json')
    Tolerence set to be 1%
    Asserts:
        1. values of stress
        2. values computed from the rheology
        3. values of strain rate
        4. composite rheology
    """
    # read parameters
    RheologyPrm = RHEOLOGY_PRM()
    diffusion_creep = RheologyPrm.HK03_diff
    dislocation_creep = RheologyPrm.HK03_disl
    
    # check for stress
    tolerance = 0.01
    check0 = CreepStress(diffusion_creep, 7.8e-15, 1e9, 1400 + 273.15, 1e4, 1000.0)
    assert(abs(check0 - 0.2991) / 0.2991 < tolerance)
    # using the second invariant
    # here the stress is converted to the 2nd invariant
    check0_ii_std = CreepStress(diffusion_creep, 7.8e-15, 1e9, 1400 + 273.15, 1e4, 1000.0) / 3**0.5 
    # here the strain rate in converted to 2nd invariant as input
    check0_ii = CreepStress(diffusion_creep, 7.8e-15/(2.0/3**0.5), 1e9,\
                            1400 + 273.15, 1e4, 1000.0, use_effective_strain_rate=True)
    assert(abs(check0_ii - check0_ii_std)/check0_ii_std < 1e-6)
    
    check1 = CreepStress(dislocation_creep, 2.5e-12, 1e9, 1400 + 273.15, 1e4, 1000.0)
    assert(abs((check1 - 0.3006) / 0.3006) < tolerance)
    # using the second invariant
    # here the stress is converted to the 2nd invariant
    check1_ii_std = CreepStress(dislocation_creep, 2.5e-12, 1e9, 1400 + 273.15, 1e4, 1000.0) / 3**0.5
    check1_ii = CreepStress(dislocation_creep, 2.5e-12/(2.0/3**0.5),\
                            1e9, 1400 + 273.15, 1e4, 1000.0, use_effective_strain_rate=True)
    assert(abs(check1_ii - check1_ii_std)/check1_ii_std < 1e-6)
    
    # check for viscosity
    check2 = CreepRheology(diffusion_creep, 7.8e-15, 1e9, 1400 + 273.15, 1e4, 1000.0)
    assert(abs((check2 - 1.917382525013947e+19) / 1.917382525013947e+19) < tolerance)
    
    check2_mg = visc_diff_HK(1400 + 273.15,1e9,1e4,1000.0,'con','orig','mid','mid')
    
    check3 = CreepRheology(dislocation_creep, 2.5e-12, 1e9, 1400 + 273.15, 1e4, 1000.0)
    assert(abs((check3 - 6.0119447368476904e+16) / check3) < tolerance)

    # use second invariant strain rate, then compute viscosity, second invariant stress and differential stress
    # accordingly, this checks the "use_effective_strain_rate" option is correctly implemented
    strain_rate_second_invariant = 2.5e-12 / (2.0/3**0.5)
    check3_1_std = 2 * 6.0119447368476904e+16 * 2.5e-12
    check3_1 = 2.0 * strain_rate_second_invariant *\
            CreepRheology(dislocation_creep, strain_rate_second_invariant, 1e9, 1400 + 273.15, 1e4, 1000.0, use_effective_strain_rate=True)\
            * 3**0.5
    assert(abs((check3_1 - check3_1_std)/check3_1_std)<1e-4)

    # check for strain rate
    check4 = CreepStrainRate(diffusion_creep, 0.2991, 1e9, 1400 + 273.15, 1e4, 1000.0)
    assert(abs(check4 - 7.8e-15) / 7.8e-15 < tolerance)
    # use the second order invariant
    check4_ii_std = CreepStrainRate(diffusion_creep, 0.2991, 1e9, 1400 + 273.15, 1e4, 1000.0) * 3**0.5 / 2.0
    check4_ii = CreepStrainRate(diffusion_creep, 0.2991/3**0.5, 1e9,\
                                1400 + 273.15, 1e4, 1000.0, use_effective_strain_rate=True)
    assert(abs(check4_ii - check4_ii_std)/check4_ii_std < 1e-6)
    
    check5 = CreepStrainRate(dislocation_creep, 0.3006, 1e9, 1400 + 273.15, 1e4, 1000.0)
    assert(abs(check5 - 2.5e-12) / 2.5e-12 < tolerance)
    # use the second order invariant
    check5_ii_std = CreepStrainRate(dislocation_creep, 0.3006, 1e9, 1400 + 273.15, 1e4, 1000.0) * 3**0.5/2.0
    check5_ii = CreepStrainRate(dislocation_creep, 0.3006/3**0.5, 1e9,\
                               1400 + 273.15, 1e4, 1000.0, use_effective_strain_rate=True)

    # check for the composite rheology
    check6_i = ComputeComposite(check2, None, None) # None inputs, should return the 1st value
    assert(abs(check6_i - check2)/check2 < 1e-6) 
    check6_ii_std = 5.99315e16
    check6_ii = ComputeComposite(check2, check3)  # returns the composite value
    assert(abs(check6_ii - check6_ii_std)/check6_ii_std < 1e-6)

    # fit new rheology: correction from variable differences
    stress_ref = 50.0 # MPa
    P_ref = 100.0e6 # Pa
    T_ref = 1250.0 + 273.15 # K
    Coh_ref = 1000.0 # H / 10^6 Si
    d_ref = 15.0 # mu m
    # test 1: trivial case: no correction 
    creep_correction = {'A': 1.0, 'p': 0.0, 'r': 0.0, 'n': 0.0, 'E': 0.0, 'V': 0.0}
    strain_rate_correction = CreepStrainRate(creep_correction, stress_ref, P_ref, T_ref, d_ref, Coh_ref)
    assert(abs(strain_rate_correction - 1.0)/1.0 < 1e-6)


def test_CreepComputeA():
    '''
    test the implementation of function CreepComputeA
    values of the rheology is taken from Hirth & Kohlstedt, 2003
    assert:

    '''
    # read parameters
    RheologyPrm = RHEOLOGY_PRM()
    diffusion_creep = RheologyPrm.HK03_diff
    dislocation_creep = RheologyPrm.HK03_disl
    check_diff_std = CreepRheology(diffusion_creep, 7.8e-15, 1e9, 1400 + 273.15, 1e4, 1000.0)
    check_disl_std = CreepRheology(dislocation_creep, 2.5e-12, 1e9, 1400 + 273.15, 1e4, 1000.0)
    # convert to aspect inputs
    diffusion_creep_A = CreepComputeA(diffusion_creep, 7.8e-15, 1e9, 1400 + 273.15, check_diff_std, d=1e4, Coh=1000.0)
    dislocation_creep_A =  CreepComputeA(dislocation_creep, 2.5e-12, 1e9, 1400 + 273.15, check_disl_std, d=1e4, Coh=1000.0)
    assert(abs(diffusion_creep_A-diffusion_creep['A'])/diffusion_creep['A']<1e-6)
    assert(abs(dislocation_creep_A-dislocation_creep['A'])/dislocation_creep['A']<1e-6)


def test_GetRheology():
    '''
    test the function of GetRheology
    assert:
        1. values of differences are correctly applied to the rheology 
        2.
    '''
    # get the original rheology
    rheology = 'HK03_wet_mod'
    diffusion_creep_0, dislocation_creep_0 = GetRheology(rheology)
    # get the rheology with modification on the variables 
    diffusion_creep, dislocation_creep = GetRheology(rheology,\
                                                     dAdiff_ratio=1.5, dEdiff=100e3, dVdiff=1e-6,\
                                                     dAdisl_ratio=2.0, dEdisl=50e3, dVdisl=2e-6)
    # values of differences are correctly applied to the rheology 
    tolerance = 1e-6
    assert(abs(diffusion_creep['A']/diffusion_creep_0['A'] - 1.5)/1.5 < tolerance) 
    assert(abs(diffusion_creep['E'] - diffusion_creep_0['E'] - 100e3)/100e3 < tolerance) 
    assert(abs(diffusion_creep['V'] - diffusion_creep_0['V'] - 1e-6)/1e-6 < tolerance) 
    assert(abs(dislocation_creep['A']/dislocation_creep_0['A'] - 2.0)/2.0 < tolerance) 
    assert(abs(dislocation_creep['E'] - dislocation_creep_0['E'] - 50e3)/50e3 < tolerance) 
    assert(abs(dislocation_creep['V'] - dislocation_creep_0['V'] - 2e-6)/2e-6 < tolerance) 
    # get the original rheology
    rheology = 'HK03_wet_mod'
    diffusion_creep, dislocation_creep = GetRheology(rheology, use_coh=False)  # use the experission of fh2o instead.
    check0 = CreepRheology(diffusion_creep, 7.8e-15, 1e9, 1400 + 273.15, 1e4, 888.38144735) # the last one is the fh2o instead of coh
    check0_std = visc_diff_HK(1400 + 273.15,1e9,1e4,1000.0,'wet','new','mid','mid')  # the tolerance has to be set bigger
    assert(abs(check0 - check0_std) / check0_std < 1e-3)


def test_HirthKohlstedt_wet_modified():
    """
    test the implementation of equations from Hirth & Kohlstedt, 2003(filename='Hirth_Kohlstedt.json')
    with a modified rheology, see magali's implementation (function visc_diff_HK)
    Tolerence set to be 1%
    Asserts:
        My implementation yield same result with Magali's
    """
    # get rheology
    rheology = 'HK03_wet_mod'
    diffusion_creep, dislocation_creep = GetRheology(rheology)
    tolerance = 0.01
    # check the value from the diffusion creep
    check0 = CreepRheology(diffusion_creep, 7.8e-15, 1e9, 1400 + 273.15, 1e4, 1000.0)
    check0_std = visc_diff_HK(1400 + 273.15,1e9,1e4,1000.0,'wet','new','mid','mid')
    assert(abs(check0 - check0_std) / check0_std < tolerance)
    # check the value from the diffusion creep
    check1 = CreepRheology(dislocation_creep, 1e-15, 1e9, 1400 + 273.15, 1e4, 1000.0)
    check1_std = visc_disl_HK(1400 + 273.15,1e9, 1e-15, 1000.0,'wet','new','mid','mid')
    assert(abs(check1 - check1_std) / check1_std < tolerance)
    pass


def test_AB17():
    '''
    test the implementation of flow low in Arredondo & Billen 2017
    '''
    rheology = 'AB17'
    diffusion_creep, dislocation_creep = GetRheology(rheology)
    # Assert 1: temperature around 300 km depth
    diff_eta_1 = CreepRheology(diffusion_creep, 1e-15, 10e9, 1400 + 273.15, 1e4, 1000.0, use_effective_strain_rate=True)
    diff_eta_1_std = 1.6358409496967286e+20
    assert(abs(diff_eta_1 - diff_eta_1_std)/diff_eta_1_std < 1e-6)
    disl_eta_1_std = 7.340152940603979e+19
    disl_eta_1 = CreepRheology(dislocation_creep, 1e-15, 10e9, 1400 + 273.15, 1e4, 1000.0, use_effective_strain_rate=True)
    eta_1 = ComputeComposite(diff_eta_1, disl_eta_1)
    eta_1_std = 7.736997237791452e+19
    assert(abs(eta_1 - eta_1_std)/eta_1_std < 1e-6)
    # Assert 2: temperature around 660 km depth
    diff_eta_1_std = 4.840046597750784e+21
    diff_eta_1 = CreepRheology(diffusion_creep, 1e-15, 6.6*3.3e9, 1400 + 273.15, 1e4, 1000.0, use_effective_strain_rate=True)
    assert(abs(diff_eta_1 - diff_eta_1_std)/diff_eta_1_std < 1e-6)
    disl_eta_1_std = 2.1018709278365573e+21
    disl_eta_1 = CreepRheology(dislocation_creep, 1e-15, 6.6*3.3e9, 1400 + 273.15, 1e4, 1000.0, use_effective_strain_rate=True)
    assert(abs(disl_eta_1 - disl_eta_1_std)/disl_eta_1_std < 1e-6)
    eta_1 = ComputeComposite(diff_eta_1, disl_eta_1)
    eta_1_std = 1.4654673144256182e+21
    assert(abs(eta_1 - eta_1_std)/eta_1_std < 1e-6)


def test_CreepComputeV():
    '''
    test CreepComputeV
    assert:
        the converted V value from flow-law derived viscosity are consistent with the original flow law
    '''
    tolerance = 1e-6
    rheology = 'HK03_wet_mod'
    diffusion_creep, dislocation_creep = GetRheology(rheology)
    # test 1: diffusion creep
    eta = CreepRheology(diffusion_creep, 7.8e-15, 1e9, 1400 + 273.15, 1e4, 1000.0, use_effective_strain_rate=1)
    V = CreepComputeV(diffusion_creep, 7.8e-15, 1e9, 1400 + 273.15, eta, use_effective_strain_rate=1)
    assert(abs((V - diffusion_creep['V'])/V) < tolerance)
    # test 2: dislocation creep
    eta = CreepRheology(dislocation_creep, 7.8e-15, 1e9, 1400 + 273.15, 1e4, 1000.0, use_effective_strain_rate=1)
    V = CreepComputeV(dislocation_creep, 7.8e-15, 1e9, 1400 + 273.15, eta, use_effective_strain_rate=1)
    assert(abs((V - dislocation_creep['V'])/V) < tolerance)


def test_CreepComputeAfromSS():
    '''
    test function CreepComputeAfromSS
    assert:
        the A value is correct
    '''
    rheology = "HK03_f" 
    diffusion_creep, dislocation_creep = GetRheology(rheology, use_coh=False)
    P = 300e6 # Pa
    T = 1250.0 + 273.15 # K
    fh2o = 300 # MPa
    d = 15.1 # mu m
    stress = 50.0 # MPa
    strain_rate = CreepStrainRate(diffusion_creep, stress, P, T, d, fh2o)
    A = CreepComputeAfromSS(diffusion_creep, strain_rate, stress, P, T, d, fh2o)
    A_std = diffusion_creep['A']
    assert(abs((A-A_std)/A_std) < 1e-6)
    strain_rate = CreepStrainRate(dislocation_creep, stress, P, T, d, fh2o)
    A = CreepComputeAfromSS(dislocation_creep, strain_rate, stress, P, T, d, fh2o)
    A_std = dislocation_creep['A']
    assert(abs((A-A_std)/A_std) < 1e-6)



def test_RheologyUpdateEV():
    '''
    test function RheologyUpdateEV
    assert:
        the updated creep law is correct
    '''
    stress = 50.0 # MPa
    P = 100.0e6 # Pa
    T = 1250.0 + 273.15 # K
    fh2o = 100.0 # MPa
    d = 15.0 # mu m
    # get the original rheology
    rheology = "HK03_f" 
    diffusion_creep, dislocation_creep = GetRheology(rheology, use_coh=False)
    # test 1: preserve the original rheology
    diffusion_creep_new = RheologyUpdateEV(diffusion_creep, stress, P, T, d, fh2o)
    A_std = diffusion_creep['A']
    A = diffusion_creep_new['A']
    assert(abs((A-A_std)/A_std) < 1e-6)
    # test 2: change the E and V values and update the rheology
    diffusion_creep_new = RheologyUpdateEV(diffusion_creep, stress, P, T, d, fh2o,\
                                            E=diffusion_creep['E'] - 40e3,\
                                            V=diffusion_creep['V'] - 5.5e-6)
    print("diffusion_creep_new:", diffusion_creep_new) # debug
    A_std = 1016892.7892374996
    A = diffusion_creep_new['A']
    assert(abs((A-A_std)/A_std) < 1e-6)


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
    # test 2: convert with the use_effective_strain_rate=True option
    check_result_std = CreepRheology(dislocation_creep, 1e-15, 10e9, 1300 + 273.15, use_effective_strain_rate=True)
    dislocation_creep_aspect = Convert2AspectInput(dislocation_creep, use_effective_strain_rate=True)
    check2 = CreepRheologyInAspectViscoPlastic(dislocation_creep_aspect, 1e-15, 10e9, 1300 + 273.15)
    assert(abs((check2 - check_result_std) / check_result_std) < tolerance)
    # test 3: the function of ConvertFromAspectInput
    diffusion_creep_1 = ConvertFromAspectInput(diffusion_creep_aspect)
    check3_diff = CreepRheology(diffusion_creep_1, 1e-15, 10e9, 1300 + 273.15)
    assert(abs((check3_diff - check_result[0]) / check_result[0]) < tolerance)
    diffusion_creep_2 = ConvertFromAspectInput(diffusion_creep_aspect, use_effective_strain_rate=True)
    check3_diff_2 = CreepRheology(diffusion_creep_2, 1e-15, 10e9, 1300 + 273.15, use_effective_strain_rate=True)
    assert(abs((check3_diff_2 - check_result[0]) / check_result[0]) < tolerance)
    dislocation_creep_1 = ConvertFromAspectInput(dislocation_creep_aspect, use_effective_strain_rate=True)
    check3_disl_1 = CreepRheology(dislocation_creep_1, 1e-15, 10e9, 1300 + 273.15)
    assert(abs((check3_disl_1 - check_result[1]) / check_result[1]) < tolerance)



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


def test_StrengthProfile():
    '''
    Test the implementation of STRENGTH_PROFILE class
    '''
    year = 365.0 * 24.0 * 3600.0
    rheology = 'HK03_wet_mod'
    Operator = STRENGTH_PROFILE(max_depth=40e3, T_type='ARCAY17')
    # Operator.SetRheologyByName(disl='ARCAY17', brittle='ARCAY17')
    Operator.SetRheologyByName(disl=rheology, brittle='ARCAY17')
    Operator.Execute(creep_type='disl')
    # assert 1: check the value of the viscosity from the viscous part
    # the P and T are taken from the variables in the STRENGTH_PROFILE functions
    _, dislocation_creep = GetRheology(rheology)
    eta_std = CreepRheology(dislocation_creep, 1e-14, 1335363814.8378572, 637.6361339382709, 1e4, 1000.0)
    assert(abs(eta_std - Operator.etas_viscous[-1])/eta_std < 1e-6)
    # assert 2: check the value of the viscosity using the 2nd invariant of stress and the strain rate
    # the P and T are taken from the variables in the STRENGTH_PROFILE functions
    eta_std = CreepRheology(dislocation_creep, 1e-14, 1335363814.8378572, 637.6361339382709, 1e4, 1000.0, use_effective_strain_rate=True)
    Operator.Execute(creep_type='disl', compute_second_invariant=True)
    assert(abs(eta_std - Operator.etas_viscous[-1])/eta_std < 1e-6)


def test_ST1981_basalt():
    '''
    test the basalt rheology from Shelton and Tullis 1981 
    and Hacker and Christie 1990
    Assert:
        The values of viscosity agrees with the values shown in fig 4 
        in Agard_etal_2016, note there are 3 lines below the mark "basalt"
        in fig4, use the one on the bottom.
    '''
    # get rheology
    rheology = 'ST1981_basalt'
    diffusion_creep, dislocation_creep = GetRheology(rheology)
    tolerance = 0.01
    # check the value from the dislocation creep
    # the last two variables (d and Coh) are not effective for this rheology
    check0 = CreepRheology(dislocation_creep, 1e-13, 0.0, 500.0 + 273.15, 1e4, 1000.0)  
    check0_std = 8.984374489053491e+20
    assert(abs(check0 - check0_std) / check0_std < tolerance)
    check1 = CreepRheology(dislocation_creep, 1e-13, 0.0, 400.0 + 273.15, 1e4, 1000.0)  
    check1_std = 4.681745069507827e+21
    assert(abs(check1 - check1_std) / check1_std < tolerance)
    

def test_RM1987_quartz():
    '''
    test the quartz rheology from Ranalli and Murphy 1987 
    Assert:
        The values of viscosity agrees with the values shown in fig 4 
        in Agard_etal_2016
    '''
    # get rheology
    rheology = 'RM1987_quartz'
    diffusion_creep, dislocation_creep = GetRheology(rheology)
    tolerance = 0.01
    # check the value from the dislocation creep
    # the last two variables (d and Coh) are not effective for this rheology
    check0 = CreepRheology(dislocation_creep, 1e-13, 0.0, 500.0 + 273.15, 1e4, 1000.0)  
    check0_std = 3.994134531813513e+19
    assert(abs(check0 - check0_std) / check0_std < tolerance)
    check1 = CreepRheology(dislocation_creep, 1e-13, 0.0, 400.0 + 273.15, 1e4, 1000.0)  
    check1_std = 10**(20.12)
    assert(abs(check1 - check1_std) / check1_std < tolerance)


def test_Rybachi_06_anorthite():
    '''
    test the anorthite flow law from the Rybachi 06
    Assert:
        1. For the dry rheologys the values of differential stresses agree with their Figure 8.
        2. For the wet rheologys the values of differential stresses agree with their Figure 8.
        3. Same as 2, but use fh2o
    Print:
        converted dry rheology in ASPECT (commenting out the relative lines)
    '''
    # assert 1: dry rheology
    rheology = 'Rybachi_06_anorthite_dry'
    diffusion_creep, dislocation_creep = GetRheology(rheology)
    strain_rate = 1e-12
    stress_diff_dry_0 = CreepStress(diffusion_creep, strain_rate, 1080e6, 841 + 273.15, 20, 5)  # 20um, 10^-12, 40km, 841 C
    assert(abs((stress_diff_dry_0 - 384.50841946490283) / 384.50841946490283) < 1e-6)  # stress = 384.5 MPa
    strain_rate = 1e-14
    stress_disl_dry_0 = CreepStress(dislocation_creep, strain_rate, 1080e6, 841 + 273.15, 20, 5)  # 20um, 10^-14, 40km, 841 C
    assert(abs((stress_disl_dry_0 - 33.32906329872165) / 33.32906329872165) < 1e-6)  # stress = 33.3 MPa
    # assert 2: wet rheology
    rheology = 'Rybachi_06_anorthite_wet'
    diffusion_creep, dislocation_creep = GetRheology(rheology)
    strain_rate = 1e-12
    # stress_diff_0 = CreepStress(diffusion_creep, strain_rate, 660e6, 482 + 273.15, 20, 1000)  # 20um, 10^-12, 20km, 482 C
    stress_diff_0 = CreepStress(diffusion_creep, strain_rate, 660e6, 482 + 273.15, 20, 5)  # 20um, 10^-12, 20km, 482 C
    assert(abs((stress_diff_0 - 433.71933526741276) / 433.71933526741276) < 1e-6)  # stress = 33.3 MPa
    stress_diff_1 = CreepStress(diffusion_creep, strain_rate, 990e6, 671 + 273.15, 20, 100)  # 20um, 10^-12, 30km, 671 C
    assert(abs((stress_diff_1 - 1.1931194545922006) / 433.71933526741276) < 1e-6)  # stress = 33.3 MPa
    strain_rate = 1e-14
    # stress_disl_0 = CreepStress(dislocation_creep, strain_rate, 660e6, 482 + 273.15, 20, 1000)  # 20um, 10^-14, 20km, 482 C
    stress_disl_0 = CreepStress(dislocation_creep, strain_rate, 660e6, 482 + 273.15, 20, 5)  # 20um, 10^-14, 20km, 482 C
    stress_disl_1 = CreepStress(dislocation_creep, strain_rate, 990e6, 671 + 273.15, 20, 100)  # 20um, 10^-14, 30km, 671 C
    assert(abs((stress_disl_0 - 792.7074001143128) / 792.7074001143128) < 1e-6)  # stress = 33.3 MPa
    assert(abs((stress_disl_1 - 15.384323814155401) / 15.384323814155401) < 1e-6)  # stress = 33.3 MPa
    tolerance = 0.1
    # assert 3: wet rheology with fh2o
    rheology = 'Rybachi_06_anorthite_wet'
    diffusion_creep, dislocation_creep = GetRheology(rheology, use_coh=False)
    strain_rate = 1e-12
    stress_diff_0 = CreepStress(diffusion_creep, strain_rate, 660e6, 482 + 273.15, 20, 500.0)  # 20um, 10^-12, 20km, 482 C
    assert(abs((stress_diff_0 - 433.0982142183824) / 433.0982142183824) < 1e-6)  # stress = 33.3 MPa


def test_Dimanov_Dresen():
    '''
    Test the plagioclase rheology from Dimanov Dresen 2005
    Asserts:
        1. Assert the 35 mu m rheology 
        2. Assert the 45 mu m rheology
        3.
        4. Assert the 35 mum wet rheology
    '''
    # assert 1: 35 mu m 
    rheology = 'Dimanov_Dresen_An50Di35D_dry'
    diffusion_creep, dislocation_creep = GetRheology(rheology)
    P = 300e6  # pa
    T = 1423.0
    d = 35.0  # in the fit, the d^(-p) term is set to 1
    stress = 200.0 # 200 Mpa
    disl_strain_rate =  CreepStrainRate(dislocation_creep, stress, P, T, d, 1.0)
    diff_strain_rate =  CreepStrainRate(diffusion_creep, stress, P, T, d, 1.0)
    strain_rate_200_35 = disl_strain_rate + diff_strain_rate
    assert(abs((strain_rate_200_35 - 8.156720213764019e-05)/8.156720213764019e-05) < 1e-6)
    stress = 50.0 # 50 Mpa
    disl_strain_rate =  CreepStrainRate(dislocation_creep, stress, P, T, d, 1.0)
    diff_strain_rate =  CreepStrainRate(diffusion_creep, stress, P, T, d, 1.0)
    strain_rate_50_35 = disl_strain_rate + diff_strain_rate
    assert(abs((strain_rate_50_35 - 6.182208873681617e-06)/6.182208873681617e-06) < 1e-6)
    # assert 2: 45 mu m
    rheology = 'Dimanov_Dresen_An50Di45D_dry'
    diffusion_creep, dislocation_creep = GetRheology(rheology)
    P = 300e6  # pa
    T = 1423.0
    d = 45.0  # in the fit, the d^(-p) term is set to 1
    stress = 200.0 # 200 Mpa
    disl_strain_rate =  CreepStrainRate(dislocation_creep, stress, P, T, d, 1.0)
    diff_strain_rate =  CreepStrainRate(diffusion_creep, stress, P, T, d, 1.0)
    strain_rate_200_45 = disl_strain_rate + diff_strain_rate
    assert(abs((strain_rate_200_45 - 6.606507553256844e-05)/6.606507553256844e-05) < 1e-6)
    stress = 50.0 # 50 Mpa
    disl_strain_rate =  CreepStrainRate(dislocation_creep, stress, P, T, d, 1.0)
    diff_strain_rate =  CreepStrainRate(diffusion_creep, stress, P, T, d, 1.0)
    strain_rate_50_45 = disl_strain_rate + diff_strain_rate
    assert(abs((strain_rate_50_45 - 2.3066772224136795e-06)/2.3066772224136795e-06) < 1e-6)
    # assert 3: 35 mu m, at a smaller differential stress and strain rate
    rheology = 'Dimanov_Dresen_An50Di35D_dry'
    diffusion_creep, dislocation_creep = GetRheology(rheology)
    P = 3000 * 10 * 3000.0  # pa, 3000 km depth, 3000 kg/m^3
    T = 925.0 + 273.0
    d = 35.0  # in the fit, the d^(-p) term is set to 1
    stress = 60.0 # 60 Mpa
    diff_strain_rate =  CreepStrainRate(diffusion_creep, stress, P, T, d, 1.0)
    disl_strain_rate =  CreepStrainRate(dislocation_creep, stress, P, T, d, 1.0)
    # print("diff_strain_rate: ", diff_strain_rate)  # debug
    # print("disl_strain_rate: ", disl_strain_rate)  # debug
    # assert 4:
    rheology = 'Dimanov_Dresen_An50Di35D_wet'
    diffusion_creep, dislocation_creep = GetRheology(rheology)
    P = 0.0  # independent
    T = 1423.0
    d = 35.0  # in the fit, the d^(-p) term is set to 1
    stress = 37.77 # Mpa
    diff_strain_rate =  CreepStrainRate(diffusion_creep, stress, P, T, d, 1.0)
    disl_strain_rate =  CreepStrainRate(dislocation_creep, stress, P, T, d, 1.0)
    strain_rate = diff_strain_rate + disl_strain_rate
    strain_rate_std = 1.22172311877202e-05
    assert(abs(strain_rate - strain_rate_std)/strain_rate_std < 1e-6)
    stress = 296.1602 # Mpa
    diff_strain_rate =  CreepStrainRate(diffusion_creep, stress, P, T, d, 1.0)
    disl_strain_rate =  CreepStrainRate(dislocation_creep, stress, P, T, d, 1.0)
    strain_rate = diff_strain_rate + disl_strain_rate
    strain_rate_std = 3.519310165608152e-4
    assert(abs(strain_rate - strain_rate_std)/strain_rate_std < 1e-6)



def test_Rybachi_2000_An100_dry():
    '''
    test the rheology of the Rybachi_2000_An100_dry rheology
    assert:
        1. diffusion creep
        2. dislocation creep
    '''
    rheology = "Rybachi_2000_An100_dry"
    diffusion_creep, dislocation_creep = GetRheology(rheology)
    # 1. diffusion creep
    P = 1.0  # pa, not dependent on the value of P
    T1 = 1428.57 # K
    d = 2.65 # mu m, a choice that I picked up myself
    stress = 10.0 # MPa
    diff_strain_rate =  CreepStrainRate(diffusion_creep, stress, P, T1, d, 0.0)
    assert(abs(diff_strain_rate - 5.6770749020979755e-06)/5.6770749020979755e-06 < 1e-6)
    T2 = 1388.88 # K, test another T
    diff_strain_rate =  CreepStrainRate(diffusion_creep, stress, P, T2, d, 0.0)
    assert(abs(diff_strain_rate - 1.8456108380428797e-06)/1.8456108380428797e-06 < 1e-6)
    # 2. dislocation creep
    T1 = 1428.57 # K
    disl_strain_rate =  CreepStrainRate(dislocation_creep, stress, P, T1, d, 0.0)
    assert(abs(disl_strain_rate - 1.012715951918425e-08)/1.012715951918425e-08 < 1e-6)
    T2 = 1388.88 # K, test another T
    disl_strain_rate =  CreepStrainRate(dislocation_creep, stress, P, T2, d, 0.0)
    assert(abs(disl_strain_rate - 2.129952815104616e-09)/2.129952815104616e-09 < 1e-6)

def test_MehlHirth08GabbroMylonite():
    '''
    test the piezometer of MehlHirth08GabbroMylonite
    fix this later when I need this relation
    Assert:
        1. a grain size is converted to the right stress
        2. the stress could be converted back to the original grain size
        3 & 4, with a different grain size of 100 mu m
        5 & 6, test the values beyond the range - a error will be raised
        7 & 8, invert 5 & 6, this time the values of d will be given with a stress beyond the range
    '''
    Piezometer = PIEZOMETER()
    d = 35 # mu m
    sigma = Piezometer.MehlHirth08GabbroMylonite(d)
    assert(abs((sigma - 43.01259895022442)/43.01259895022442) < 1e-6)
    sigma = 43.01259895022442 # invert the relation
    d = Piezometer.MehlHirth08GabbroMyloniteInvert(sigma)
    # assert(abs((d - 35.0)/35.0) < 1e-6)
    # 3 & 4
    d = 100 # mu m
    sigma = Piezometer.MehlHirth08GabbroMylonite(d)
    # assert(abs((sigma - 25.986533248593485)/25.986533248593485) < 1e-6)
    sigma = 25.986533248593485 # invert the relation
    d = Piezometer.MehlHirth08GabbroMyloniteInvert(sigma)
    assert(abs((d - 100.0)/100.0) < 1e-6)
    # 5 & 6 
    d = 4000 # mu m, value beyond the maximum
    with pytest.raises(ValueError) as excinfo:
        sigma = Piezometer.MehlHirth08GabbroMylonite(d)
        assert('maximum limit' in str(excinfo.value))
    d = 5 # mu m, value smaller than the minimum
    with pytest.raises(ValueError) as excinfo:
        sigma = Piezometer.MehlHirth08GabbroMylonite(d)
        assert('minimum limit' in str(excinfo.value))
    # 7 & 8
    sigma = 200 # invert the relation
    d = Piezometer.MehlHirth08GabbroMyloniteInvert(sigma)
    # assert(abs((d - 10.0)/10.0) < 1e-6)
    sigma = 5 # invert the relation
    d = Piezometer.MehlHirth08GabbroMyloniteInvert(sigma)
    assert(abs((d - 3000.0)/3000.0) < 1e-6)

def test_MK10_peierls():
    '''
    test the MK10 rheology
    assert:
        1. strain rates match with the values from the figure 5 in the original figure
        2. stress value could match vice versa
        3. the viscosity
    '''
    creep = GetPeierlsRheology("MK10")
    # assert 1.1, values at a strain rate = 3e-5
    strain_rate_std = 2.7778604629458894e-05
    stress = 3.82e3 # MPa
    T = 671.0 # K
    P = 0 # not dependent on P (V = 0)
    strain_rate = PeierlsCreepStrainRate(creep, stress, P, T)
    assert(abs(strain_rate_std - strain_rate)/strain_rate_std < 1e-6)
    # assert 1.2
    strain_rate_std = 2.6215973389278528e-05
    stress = 3.25e3 # MPa
    T = 907.0 # K
    P = 0 # not dependent on P (V = 0)
    strain_rate = PeierlsCreepStrainRate(creep, stress, P, T)
    assert(abs(strain_rate_std - strain_rate)/strain_rate_std < 1e-6)
    # assert 1.3: use a variation with activation volume
    strain_rate_std = 2.9979043137974457e-05
    T = 873.0
    P = 4.5e9
    stress = 3.35e3 # MPa
    strain_rate = PeierlsCreepStrainRate(creep, stress, P, T) 
    assert(abs(strain_rate_std - strain_rate)/strain_rate_std < 1e-6)
    strain_rate_std = 2.9979043137974457e-05
    T = 873.0
    P = 4.5e9
    dV = 30e-6
    stress = 3.35e3 # MPa
    strain_rate = PeierlsCreepStrainRate(creep, stress, P, T, dV=dV) 
    assert(abs(strain_rate_std - strain_rate)/strain_rate_std < 1e-6)
    strain_rate_std = 1.0823657173643094e-05
    T = 873.0
    P = 5.5e9
    dV = 30e-6
    stress = 3.35e3 # MPa
    strain_rate = PeierlsCreepStrainRate(creep, stress, P, T, dV=dV) 
    assert(abs(strain_rate_std - strain_rate)/strain_rate_std < 1e-6)
    # assert 2.1
    strain_rate = 2.7778604629458894e-05
    T = 671.0 # K
    P = 0 # not dependent on P (V = 0)
    stress_std = 3.82e3 # MPa
    stress = PeierlsCreepStress(creep, strain_rate, P, T)
    assert(abs(np.log(stress_std/stress)) < 0.05)  # a bigger tolerance, check the log value
    # assert 2.2
    strain_rate = 2.6215973389278528e-05
    T = 907.0 # K
    P = 0 # not dependent on P (V = 0)
    stress_std = 3.25e3 # MPa
    stress = PeierlsCreepStress(creep, strain_rate, P, T)
    assert(abs(np.log(stress_std/stress)) < 0.05)  # a bigger tolerance, check the log value
    # assert 3.1: viscosity
    strain_rate = 2.7778604629458894e-05
    T = 671.0 # K
    P = 0 # not dependent on P (V = 0)
    eta_std = 6.8757953e13
    eta = PeierlsCreepRheology(creep, strain_rate, P, T)
    assert(abs(np.log(eta_std/eta)) < 0.05)  # a bigger tolerance, check the log value
    # assert 3.2: viscosity, a realistic scenario
    # 1e-13
    strain_rate = 1e-13
    T = 800.0 + 273.15 # K
    P = 0 # not dependent on P (V = 0)
    eta_std = 2.125142090914006e+21
    eta = PeierlsCreepRheology(creep, strain_rate, P, T)
    assert(abs(np.log(eta_std/eta)) < 0.05)  # a bigger tolerance, check the log value
    # 1e-15
    strain_rate = 1e-15
    T = 800.0 + 273.15 # K
    P = 0 # not dependent on P (V = 0)
    eta_std = 9.893654111645399e+22
    eta = PeierlsCreepRheology(creep, strain_rate, P, T)
    assert(abs(np.log(eta_std/eta)) < 0.05)  # a bigger tolerance, check the log value
    

def test_Idrissi16_peierls():
    '''
    test the Idrissi16 rheology
    assert:
        1. strain rates match with the values from the figure 1 in the original figure
        2. stress value could match vice versa
        3. the viscosity
    '''
    creep = GetPeierlsRheology("Idrissi16")
    # assert 1.1, values at a strain rate = 3e-5
    strain_rate_std = 1.0e-05
    stress = 1.61e3 # MPa
    T = 318 # K
    P = 0 # not dependent on P (V = 0)
    strain_rate = PeierlsCreepStrainRate(creep, stress, P, T)
    assert(abs(strain_rate_std - strain_rate)/strain_rate_std < 1e-6)
    # assert 1.2, values at a strain rate = 
    strain_rate_std = 1.0e-05
    stress = 0.61e3 # MPa
    T = 909 # K
    P = 0 # not dependent on P (V = 0)
    strain_rate = PeierlsCreepStrainRate(creep, stress, P, T)
    assert(abs(strain_rate_std - strain_rate)/strain_rate_std < 1e-6)


def test_CoulumbYielding():
    '''
    test the function CoulumbYielding
    '''
    # test 1: dry profile
    P = 500e6
    cohesion=2e6
    friction=0.6
    tau_std = 302e6
    tau = CoulumbYielding(P, cohesion, friction)
    assert(abs(tau - tau_std)/tau_std < 1e-6)
    # test 2: Otherwise same as 1, hydrate profile with _lambda = 0.1
    P = 500e6
    cohesion=2e6
    friction=0.6
    _lambda=0.1
    tau_std = 271.8e6
    tau = CoulumbYielding(P, cohesion, friction, _lambda)
    assert(abs(tau - tau_std)/tau_std < 1e-6)


def test_Byerlee():
    '''
    test the Byerlee function
    '''
    # test 1: dry condition
    P = np.linspace(0.0, 100e6, 10)
    tau = Byerlee(P)
    assert(abs(tau[-1] - 85e6)/85e6 < 1e-6)
    # test 2: wet condition
    P = 330e6
    _lambda = 0.99
    tau = Byerlee(P, _lambda=_lambda)
    tau_std = 2.58e6 # 2.58 MPa
    assert(abs(tau - tau_std)/tau_std < 1e-6)


def test_PeierlsStressLimiter():
    '''
    test the stress limiter approach
    '''
    # test 1: strain rate = reference strain rate
    tau_y = 500e6 # 500 Mpa
    strain_rate_y = 1e-15
    ny = 10
    strain_rate = 1e-15
    tau = PeierlsStressLimiter(tau_y, strain_rate_y, ny, strain_rate)
    tau_std = 500e6
    assert(abs(tau-tau_std)/tau_std < 1e-6)
    # test 2: strain rate is 2 orders higher
    strain_rate = 1e-13
    tau = PeierlsStressLimiter(tau_y, strain_rate_y, ny, strain_rate)
    tau_std = 792.4465962305567e6
    assert(abs(tau-tau_std)/tau_std < 1e-6)



# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

