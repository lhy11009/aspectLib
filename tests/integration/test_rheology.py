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
from shilofue.Rheology import *  # import test module
# from shilofue.Utilities import 
# from matplotlib import pyplot as plt
# from shutil import rmtree  # for remove directories

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_rheology')

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)


def test_HK03_mod_whole_mantle_apsect_prm():
    """
    test converting the HK03 modified rheology to a whole-mantle rheology in aspect
    and converts to a prm file
    """
    da_file = os.path.join(source_dir, "depth_average.txt")
    assert(os.path.isfile(da_file))
    Operator = RHEOLOGY_OPR()
    # read profile
    Operator.ReadProfile(da_file)
    rheology_aspect = Operator.MantleRheology_v0(rheology="HK03_wet_mod", dEdiff=-40e3, dEdisl=20e3, dVdiff=-5.5e-6, dVdisl=0.0)
    diffusion_creep = rheology_aspect['diffusion_creep']
    diffusion_A = diffusion_creep['A']
    diffusion_E = diffusion_creep['E']
    diffusion_V = diffusion_creep['V']
    assert(abs(diffusion_A - 8.17868945868946e-17)/8.17868945868946e-17 < 1e-6)
    assert(abs(diffusion_E - 285000.0) / 285000.0 < 1e-6)
    assert(abs(diffusion_V - 6.9e-6)/6.9e-6 < 1e-6)
    dislocation_creep = rheology_aspect['dislocation_creep']
    dislocation_A = dislocation_creep['A']
    dislocation_E = dislocation_creep['E']
    dislocation_V = dislocation_creep['V']
    assert(abs(dislocation_A - 5.907603165099757e-16)/5.907603165099757e-16 < 1e-6)
    assert(abs(dislocation_E - 490000.0)/490000.0 < 1e-6)
    assert(abs(dislocation_V - 1.34e-5)/1.34e-5 < 1e-6)
    diffusion_creep_lm = rheology_aspect['diffusion_lm']
    diffusion_lm_A = diffusion_creep_lm['A']
    diffusion_lm_V = diffusion_creep_lm['V']
    assert(abs(diffusion_lm_A - 2.5003579819454594e-19)/2.5003579819454594e-19 < 1e-6)
    assert(abs(diffusion_lm_V - 3e-6)/3e-6 < 1e-6)
    pass

    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))
