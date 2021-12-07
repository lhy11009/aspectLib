# -*- coding: utf-8 -*-
r"""Test for PlotCombine.py

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
from shilofue.PlotCombine import *  # import test module
# from shilofue.Utilities import 
# from matplotlib import pyplot as plt
# from shutil import rmtree  # for remove directories

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_plot_combine')
case_dir1 = os.path.join(os.path.dirname(__file__), 'fixtures',\
'cases', 'test_multiple_cases', 'wb_sph_cdd50_substract_T')  # case directory for testing
assert(os.path.isdir(case_dir1))

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)

def test_plot_combine():
    '''
    test the PLOT_COMBINE class
    Asserts:
    '''
    # test the interface
    json_path = os.path.join(source_dir, 'test.json')
    assert(os.access(json_path, os.R_OK))
    Pc_opt = PC_OPT()
    Pc_opt.read_json(json_path)  # read options
    # assert something 
    result = Pc_opt.to_PC_init()
    result_std = os.path.join(os.environ['ASPECT_LAB_DIR'],\
    'tests/integration/fixtures/cases/test_multiple_cases/wb_sph_cdd50_substract_T')
    assert(len(result) == 2 and result[0] == result_std)
    # test the class
    figure_path = os.path.join(ASPECT_LAB_DIR, ".test", "combined", "foo.png")
    if os.path.isfile(figure_path):
        os.remove(figure_path)
    Plot_Combine = PLOT_COMBINE(Pc_opt.to_PC_init())
    Plot_Combine.set_plots(Pc_opt.to_PC_set_plots())
    Plot_Combine(*Pc_opt.to_PC_call())
    assert(os.path.isfile(figure_path))

    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

