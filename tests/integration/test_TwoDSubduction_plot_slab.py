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
from shutil import rmtree  # for remove directories

test_dir = ".test"


if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)

def test_vtk_TwoDSubduction_SlabAnalysis():
    '''
    test the TwoDSubduction_SlabAnalysis module of vtk
    assert that we get the right value for trench position and slab depth
    '''
    case_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'cases', 'test_TwoDSubduction_vtk_sph')
    local_test_dir = os.path.join(test_dir, 'vtk_TwoD0')
    if os.path.isdir(local_test_dir):
        rmtree(local_test_dir)
    os.mkdir(local_test_dir)
    option_path = os.path.join(local_test_dir, 'TwoDSubduction_SlabAnalysis.input')
    wedge_T_file_out = os.path.join(local_test_dir, 'wedge_T100_00001.txt')
    wedge_T_file_out_std = os.path.join(case_dir, 'wedge_T100_00001_std.txt')
    if os.path.isfile(wedge_T_file_out):
        os.remove(wedge_T_file_out)  # remove older file
    vtk_option_path, _, _ = PrepareVTKOptions(VISIT_OPTIONS, case_dir,\
        'TwoDSubduction_SlabAnalysis', vtk_step=0, output=option_path, operation='morphology')
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


def test_vtk_TwoDSubduction_SlabAnalysis_Buoyancy():
    '''
    test the TwoDSubduction_SlabAnalysis module of vtk
    assert that we get the right value for trench position and slab depth
    '''
    case_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'cases', 'test_TwoDSubduction_vtk_sph')
    local_test_dir = os.path.join(test_dir, 'vtk_TwoD1')
    if os.path.isdir(local_test_dir):
        rmtree(local_test_dir)
    os.mkdir(local_test_dir)
    option_path = os.path.join(local_test_dir, 'TwoDSubduction_SlabAnalysis_Buoyancy.input')
    wedge_T_file_out = os.path.join(local_test_dir, 'wedge_T100_00001.txt')
    wedge_T_file_out_std = os.path.join(case_dir, 'wedge_T100_00001_std.txt')
    slab_surface_file_out = os.path.join(local_test_dir, 'slab_surface_00001.txt')
    slab_surface_file_out_std = os.path.join(case_dir, 'slab_surface_00001_std.txt')
    slab_internal_file_out = os.path.join(local_test_dir, 'slab_internal_00001.txt')
    slab_internal_file_out_std = os.path.join(case_dir, 'slab_internal_00001_std.txt')
    slab_forces_file_out = os.path.join(local_test_dir, 'slab_forces_00001.txt')
    slab_forces_file_out_std = os.path.join(case_dir, 'slab_forces_00001_std.txt')
    vtk_option_path, _, _ = PrepareVTKOptions(VISIT_OPTIONS, case_dir, 'TwoDSubduction_SlabAnalysis', vtk_step=0, output=option_path)
    _stdout = RunVTKScripts('TwoDSubduction_SlabAnalysis', vtk_option_path)
    outputs = slab_morph(_stdout)
    outputs1 = slab_buoyancy(_stdout)
    # compare trench & slab depth output
    trench_theta_std = 0.63215
    slab_depth_std = 180851
    total_buoyancy_std = -1.32339e+12
    assert(abs(outputs['trench_theta']-trench_theta_std)/trench_theta_std < 1e-6)
    assert(abs(outputs['slab_depth']-slab_depth_std)/slab_depth_std < 1e-6)
    assert(abs(outputs1['total_buoyancy']-total_buoyancy_std)/total_buoyancy_std < 1e-6)
    # compare wedge temperature output
    assert(os.path.isfile(wedge_T_file_out))
    assert(filecmp.cmp(wedge_T_file_out, wedge_T_file_out_std))
    # compare slab surface output
    assert(os.path.isfile(slab_surface_file_out))
    assert(filecmp.cmp(slab_surface_file_out, slab_surface_file_out_std))
    # compare slab internal output
    assert(os.path.isfile(slab_internal_file_out))
    assert(filecmp.cmp(slab_internal_file_out, slab_internal_file_out_std))
    # compare slab internal output
    assert(os.path.isfile(slab_forces_file_out))
    assert(filecmp.cmp(slab_forces_file_out, slab_forces_file_out_std))


def test_vtk_TwoDSubduction_SlabAnalysis_Cart():
    '''
    test the TwoDSubduction_SlabAnalysis module of vtk in cartesian geometry
    assert that we get the right value for trench position and slab depth
    '''
    local_test_dir = os.path.join(test_dir, 'vtk_TwoD_cart')
    if os.path.isdir(local_test_dir):
        rmtree(local_test_dir)
    os.mkdir(local_test_dir)
    case_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'cases', 'test_TwoDSubduction_vtk_cart')
    option_path = os.path.join(local_test_dir, 'TwoDSubduction_SlabAnalysis.input')
    wedge_T_file_out = os.path.join(local_test_dir, 'wedge_T100_00001.txt')
    wedge_T_file_out_std = os.path.join(case_dir, 'wedge_T100_00001_std.txt')
    vtk_option_path, _, _ = PrepareVTKOptions(VISIT_OPTIONS, case_dir,\
        'TwoDSubduction_SlabAnalysis', vtk_step=0, output=option_path, operation='morphology')
    _stdout = RunVTKScripts('TwoDSubduction_SlabAnalysis', vtk_option_path)
    outputs = slab_morph(_stdout)
    # compare trench & slab depth output
    trench_theta_std = 4.00983e+06
    slab_depth_std = 203429
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

