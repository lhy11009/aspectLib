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
from shilofue.Analytics import *  # import test module
# from shilofue.Utilities import 
# from matplotlib import pyplot as plt
# from shutil import rmtree  # for remove directories

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'parse')

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)

def test_HAGER_CONRAD1999():
    '''
    test the HAGER_CONRAD1999 class
    Asserts:
        the computed convergence rate (refer to figure 9 in Behr etal 2022)
    '''
    # test 1 & 2: a strong slab and a weak slab
    HC1999_strong = HAGER_CONRAD1999(1e23)
    HC1999_weak = HAGER_CONRAD1999(3e21)

    # length, thickness and the bending curvature of the slab
    L_l = 660e3 # m
    h_l = 80e3 # m
    R_l = 250e3 # m

    # aspect ratio of the shear zone
    zeta_f = 20.0

    # shear zone viscosity
    eta_sz = 2.5e20 # Pa s

    # compute the convergence rate 
    Vc_strong = HC1999_strong.ComputeConvergence(L_l, h_l, R_l, eta_sz, zeta_f)
    Vc_weak = HC1999_weak.ComputeConvergence(L_l, h_l, R_l, eta_sz, zeta_f)

    # assert something 
    Vc_strong_std = 1.2762879844550031e-09 # m / yr
    Vc_weak_std = 2.3304345002807513e-09
    assert(abs(Vc_strong - Vc_strong_std) / Vc_strong_std < 1e-6)
    assert(abs(Vc_weak - Vc_weak_std) / Vc_weak_std < 1e-6)

    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

