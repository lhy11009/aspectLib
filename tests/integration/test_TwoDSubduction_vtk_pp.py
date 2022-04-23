# -*- coding: utf-8 -*-
r"""Test for VtkPp.py

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
import shilofue.VtkPp as VtkPp
from shilofue.TwoDSubduction0.VtkPp import *  # import test module
# from shilofue.Utilities import 
# from matplotlib import pyplot as plt
from shutil import rmtree  # for remove directories
import vtk

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'parse')


if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)


def test_prepare_slab():
    '''
    Test utilities from VtkPp.py
    0: read in pvtu file
    test 1: output slab grid
    '''
    # assert something 
    case_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_vtk_pp')
    output_path = os.path.join(test_dir, "vtkp_readfile")
    if os.path.isdir(output_path):
        rmtree(output_path)  # remove old results
    os.mkdir(output_path)
    filein = os.path.join(case_dir, "output", "solution", "solution-00002.pvtu")
    assert(os.path.isfile(filein))
    VtkP = VTKP()
    VtkP.ReadFile(filein)
    field_names = ['T', 'density', 'spcrust', 'spharz']
    VtkP.ConstructPolyData(field_names, include_cell_center=True)
    VtkP.PrepareSlab(['spcrust', 'spharz'])
    # test 1 output slab grid & envelop
    fileout = os.path.join(output_path, 'slab.vtu')
    fileout_std = os.path.join(case_dir, 'slab_std.vtu')
    slab_grid = VtkP.ExportSlabInternal()
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputData(slab_grid)
    writer.SetFileName(fileout)
    writer.Update()
    writer.Write()
    assert(os.path.isfile(fileout))  # assert file existence
    assert(filecmp.cmp(fileout_std, fileout))  # compare file extent
    slab_envelop0, slab_envelop1 = VtkP.ExportSlabEnvelopCoord() # envelop
    assert(abs(slab_envelop0[0, 0] - 5086593.875)/5086593.875 < 1e-8)
    assert(abs(slab_envelop0[6, 1] - 3723199.125)/3723199.125< 1e-8)
    # test 2, slab buoyancy
    r0_range = [6371e3 - 2890e3, 6371e3]
    x1 = 0.01 
    n = 100
    v_profile = VtkP.VerticalProfile2D(r0_range, x1, n)
    total_buoyancy, b_profile = VtkP.SlabBuoyancy(v_profile, 5e3)
    assert(abs(total_buoyancy + 5394703810473.24)/5394703810473.24 < 1e-8)
    assert((b_profile[11, 1]-1.09996017e+12)/1.09996017e+12 < 1e-8)


def test_export_slab_info():
    '''
    Test slab properties from VtkPp.py
    assert:
        1. value of trench position, slab_depth, dip angle
    '''
    case_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'test_vtk_pp')
    output_path = os.path.join(test_dir, "vtkp_readfile")
    if os.path.isdir(output_path):
        rmtree(output_path)  # remove old results
    os.mkdir(output_path)
    filein = os.path.join(case_dir, "output", "solution", "solution-00002.pvtu")
    assert(os.path.isfile(filein))
    VtkP = VTKP()
    VtkP.ReadFile(filein)
    field_names = ['T', 'density', 'spcrust', 'spharz']
    VtkP.ConstructPolyData(field_names, include_cell_center=True)
    VtkP.PrepareSlab(['spcrust', 'spharz'])
    trench, slab_depth, dip_100 = VtkP.ExportSlabInfo()
    assert(abs(trench - 0.6342158389165757)/0.6342158389165757 < 1e-6)
    assert(abs(slab_depth - 191927.42159304488)/191927.42159304488 < 1e-6)
    assert(abs(dip_100 - 0.22073102845130024)/0.22073102845130024 < 1e-6)

# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

