# -*- coding: utf-8 -*-
r"""Test for ThDSubduction0/VtkPp.py

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
import numpy as np
# from shilofue.Utilities import 
# from matplotlib import pyplot as plt
# from shutil import rmtree  # for remove directories
from shutil import rmtree  # for remove directories
import vtk
from shilofue.ThDSubduction0.VtkPp import *  # import test module

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'ThDSubduction', 'test_vtk_pp_slab')
source_case_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'ThDSubduction', 'test_cases')

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)

# todo2
def test_prepare_slab():
    '''
    Test utilities from VtkPp.py
    0: read in pvtu file
    test 1: output slab grid
    '''
    case_dir = os.path.join(source_case_dir, 'test_prepare_slab')
    output_path = os.path.join(test_dir, "ThDSubduction_vtk_test_prepare_slab")
    if os.path.isdir(output_path):
        rmtree(output_path)  # remove old results
    os.mkdir(output_path)
    filein = os.path.join(case_dir, "output", "solution", "solution-00000.pvtu")
    assert(os.path.isfile(filein))
    # prepare visit options
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    geometry = Visit_Options.options['GEOMETRY']
    Ro =  Visit_Options.options['OUTER_RADIUS']
    VtkP = VTKP(geometry=geometry, Ro=Ro)
    VtkP.ReadFile(filein)
    field_names = ['T', 'density', 'sp_upper', 'sp_lower']
    VtkP.ConstructPolyData(field_names, include_cell_center=True, fix_cell_value=False)
    VtkP.PrepareSlabByPoints(['sp_upper', 'sp_lower'])
    assert(abs(VtkP.slab_depth - 125e3)/125e3 < 1e-6)
    print("slab depth: ", VtkP.slab_depth)  # debug

    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

