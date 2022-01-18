# -*- coding: utf-8 -*-
r"""Test for PlotSlab.py in the TwoDSubduction0 subfolder

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
from shilofue.TwoDSubduction0.PlotSlab import *  # import test module
from shilofue.PlotVisit import RunVTKScripts, PrepareVTKOptions
from shilofue.TwoDSubduction0.PlotVisit import VISIT_OPTIONS
# from shilofue.Utilities import 
# from matplotlib import pyplot as plt
# from shutil import rmtree  # for remove directories

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'TwoDSubduction', 'test_plot_slab')
case_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'cases', 'test_vtk')

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)

def test_vtk_TwoDSubduction_SlabAnalysis():
    '''
    test the TwoDSubduction_SlabAnalysis module of vtk
    assert that we get the right value for trench position and slab depth
    '''
    option_path = os.path.join(test_dir, 'TwoDSubduction_SlabAnalysis.input')
    wedge_T_file_out = os.path.join(test_dir, 'wedge_T100_00001.txt')
    wedge_T_file_out_std = os.path.join(case_dir, 'wedge_T100_00001_std.txt')
    if os.path.isfile(wedge_T_file_out):
        os.remove(wedge_T_file_out)  # remove older file
    vtk_option_path, _, _ = PrepareVTKOptions(VISIT_OPTIONS, case_dir, 'TwoDSubduction_SlabAnalysis', vtk_step=0, output=option_path)
    _stdout = RunVTKScripts('TwoDSubduction_SlabAnalysis', vtk_option_path)
    outputs = slab_morph(_stdout)
    # compare trench & slab depth output
    trench_theta_std = 0.63215
    slab_depth_std = 180851
    assert(abs(outputs['trench_theta']-trench_theta_std)/trench_theta_std < 1e-6)
    assert(abs(outputs['slab_depth']-slab_depth_std)/slab_depth_std < 1e-6)
    # compare wedge temperature output
    assert(os.path.isfile(wedge_T_file_out))
    assert(filecmp.cmp(wedge_T_file_out, wedge_T_file_out_std))

    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

