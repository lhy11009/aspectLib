# -*- coding: utf-8 -*-
r"""Test for foo.py

This outputs:

  - test results to value of variable "test_dir"

This depends on:

  - source files written from the value of "source_dir"

Examples of usage:

  - default usage:

        python -m pytest test_foo.py

descriptions:
    every function is a separate test, combined usage with pytest module
""" 


import os
# import pytest
# import filecmp  # for compare file contents
import numpy as np
from shilofue.ThermalModel import *  # import test module
# from shilofue.Utilities import 
# from matplotlib import pyplot as plt
# from shutil import rmtree  # for remove directories

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'parse')

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)

def test_plate_model():
    '''
    Test the implementation of the plate model
    Asserts:
        a. the value of PM_A is the same with my hand calculation
        a. the value of PM_B is the same with my hand calculation
    '''
    tolerance = 1e-6
    year = 365 * 24 * 3600.0  # s in year
    Pmodel = PLATE_MODEL(150e3, 1e-6, 273.0, 1673.0, 0.05/year) # initiate the plate model
    # assert the PM_A factor is consistent with my hand calculation
    assert(abs(Pmodel.PM_A(1, 1e6*year) - 0.6278754063131151)/0.6278754063131151 < 1e-6)
    assert(abs(Pmodel.PM_A(3, 1e6*year) - 0.1874020031998469)/0.1874020031998469 < 1e-6)
    # assert the PM_B factor is consistent with my hand calculation
    assert(abs(Pmodel.PM_B(0, 1e6*year) - 59957.684736338706)/59957.684736338706 < 1e-6)
    assert(abs(Pmodel.PM_B(1, 1e6*year) - 5965.19103094)/5965.19103094 < 1e-6)

    year = 365.25 * 24 * 3600.0  # s in year
    Pmodel1 = PLATE_MODEL(150e3, 0.804e-6, 273.0, 1673.0, 0.05/year) # initiate the plate model
    print(Pmodel1.PM_B(0, 1.38895e+15))


def test_mantle_adiabat():
    '''
    test the implementation of the class MANTLE_ADIABAT
    assert:
    1. the temperature match with the cmb temperature
    2. the approximated temperature from the constant gradient match with the theratical value
    '''
    # earth values
    cp = 1250.0
    alpha = 3e-5
    g = 9.8
    Ts = 1573.0
    zcmb = 2890e3
    MantleAdiabat = MANTLE_ADIABAT(cp, alpha, g, Ts, approx="constant variables")
    Tcmb = MantleAdiabat.Temperature(zcmb)
    Tcmb_std = 3104.0652
    assert(abs(Tcmb - Tcmb_std)/Tcmb_std < 1e-6)
    MantleAdiabatConstGradient = MANTLE_ADIABAT(cp, alpha, g, Ts, approx="constant variables and gradient")
    Tcmb_approx = MantleAdiabatConstGradient.Temperature(zcmb)
    Tcmb_approx_std = 2642.21
    assert(abs(Tcmb_approx - Tcmb_approx_std)/Tcmb_approx_std < 1e-6)


    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

