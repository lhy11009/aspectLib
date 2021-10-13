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
    '''
    option_path = os.path.join(test_dir, 'TwoDSubduction_SlabAnalysis.input')
    vtk_option_path = PrepareVTKOptions(case_dir, 'TwoDSubduction_SlabAnalysis', vtk_step=0, output=option_path)
    output_file = os.path.join(test_dir, 'contour_slab00002.txt')
    output_std_file = os.path.join(case_dir, 'vtk_outputs', 'contour_slab_std.txt')
    assert(os.access(output_std_file, os.R_OK))
    if os.path.isfile(output_file):
        os.remove(output_file)
    RunVTKScripts('TwoDSubduction_SlabAnalysis', vtk_option_path)
    assert(os.path.isfile(output_file)) # assert file exists
    assert(filecmp.cmp(output_file, output_std_file))  # assert file contents


def test_slab_morph():
    '''
    test the SlabMorph function
    Asserts:
    '''
    tolerance = 1e-6
    vtk_file_path = os.path.join(case_dir, 'vtk_outputs', 'contour_slab_std.txt')
    assert(os.access(vtk_file_path, os.R_OK))  # assert file exists
    slab_morph_outputs = slab_morph(vtk_file_path)
    std_value = 0.6280194126721352
    assert(abs(slab_morph_outputs['trench']['theta']-std_value)/std_value < tolerance)  # assert trench position
    # assert something 
    assert(True)

    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

