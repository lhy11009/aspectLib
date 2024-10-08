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
import filecmp  # for compare file contents
# import numpy as np
# from shilofue.Utilities import 
# from matplotlib import pyplot as plt
# from shutil import rmtree  # for remove directories
from shilofue.PlotHeatFlow import *  # import test module

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_plot_heat_flux')

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)

def test_boundary_outputs():
    '''
    Test function for the BOUNDARYOUTPUTs class to ensure that boundary data is read and processed correctly.
    This test asserts that the boundary indicators for specific locations are correct and that the total count of 
    bottom and top boundaries matches the expected value of 2325.
    '''
    filein = os.path.join(source_dir, "heat_flux.00800")

    BdOutputs = BOUNDARYOUTPUTs(2, geometry="chunk")
    BdOutputs.ReadFile(filein)
    BdOutputs.ProcessBds()

    # Verify the location of boundary indicators:
    # 0 - left; 1 - right; 2 - bottom; 3 - top
    '''
    Assert that the first boundary indicator is "bottom" (2) and that the 
    230th boundary indicator is "top" (3). Additionally, assert that the total
    number of "bottom" and "top" boundaries combined equals 2325.
    '''
    assert(BdOutputs.bd_indicators[0] == 2)  
    assert(BdOutputs.bd_indicators[229] == 3)  
    assert(np.where(BdOutputs.bd_indicators == 2)[0].size + np.where(BdOutputs.bd_indicators == 3)[0].size == 2325)
        
    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

