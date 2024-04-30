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

test_dir = os.path.join(".test", "test_ThDSubduction_VtkPp")
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'ThDSubduction', 'test_vtk_pp_slab')
big_source_dir = os.path.join(ASPECT_LAB_DIR, "tests", "integration", "big_fixtures", "ThDSubduction")
source_case_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'ThDSubduction', 'test_cases')

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)

# todo_hv

def test_extract_slab_cross_section_at_depth():
    '''
    Test extract horizontal flow field
    Here, I test using the SlabMorphology function to get the cross section of the slab at a given depth
    '''
    case_dir = os.path.join(big_source_dir, "eba3d_width61_c23_AR4")
    assert(os.path.isdir(case_dir))
    std_vtk_outputs_dir = os.path.join(case_dir, "vtk_outputs_std") # directory contains the standard outputs
    assert(os.path.isdir(std_vtk_outputs_dir))
    odir = os.path.join(test_dir, "test_extract_slab_cross_section_at_depth")
    if os.path.isdir(odir):
        rmtree(odir)  # remove old results and make a new directory
    os.mkdir(odir)

    filein = os.path.join(case_dir, "output", "solution", "solution-00144.pvtu")
    assert(os.path.isfile(filein))

    o_vtk_dir = os.path.join(odir, 'vtk_outputs')
    if os.path.isdir(o_vtk_dir):
        rmtree(o_vtk_dir) # remove old results

    slab_envelop_interval_y = 20e3  # Interval along x axis to sort out the trench locations
    slab_envelop_interval_z = 20e3  # Interval along z axis to sort out the trench locations
    slab_shallow_cutoff = 40e3  # Minimum depth along z axis to sort out the trench locations
    crust_only = 1  # If we only use the crustal composition to sort out the trench locations

    SlabMorphology(case_dir, 144, slab_envelop_interval_y=slab_envelop_interval_y, slab_envelop_interval_z=slab_envelop_interval_z,\
        slab_shallow_cutoff=slab_shallow_cutoff, crust_only=crust_only, output=o_vtk_dir, horizontal_velocity_depths=[150e3])
    
    slab_surface_file = os.path.join(o_vtk_dir, "slab_surface_00144_d150.00km.txt")
    assert(os.path.isfile(slab_surface_file))


def test_slab_morphology():
    '''
    Test extract horizontal flow field
    Assert:
        check the file that contains the envelops of the slab surface and the
        file contains the trench locations
    '''
    case_dir = os.path.join(big_source_dir, "eba3d_width61_c23_AR4")
    assert(os.path.isdir(case_dir))
    std_vtk_outputs_dir = os.path.join(case_dir, "vtk_outputs_std") # directory contains the standard outputs
    assert(os.path.isdir(std_vtk_outputs_dir))
    odir = os.path.join(test_dir, "test_slab_morphology")
    if os.path.isdir(odir):
        rmtree(odir)  # remove old results and make a new directory
    os.mkdir(odir)

    filein = os.path.join(case_dir, "output", "solution", "solution-00144.pvtu")
    assert(os.path.isfile(filein))

    o_vtk_dir = os.path.join(odir, 'vtk_outputs')
    if os.path.isdir(o_vtk_dir):
        rmtree(o_vtk_dir) # remove old results

    slab_envelop_interval_y = 20e3  # Interval along x axis to sort out the trench locations
    slab_envelop_interval_z = 20e3  # Interval along z axis to sort out the trench locations
    slab_shallow_cutoff = 40e3  # Minimum depth along z axis to sort out the trench locations
    crust_only = 1  # If we only use the crustal composition to sort out the trench locations

    SlabMorphology(case_dir, 144, slab_envelop_interval_y=slab_envelop_interval_y, slab_envelop_interval_z=slab_envelop_interval_z,\
        slab_shallow_cutoff=slab_shallow_cutoff, crust_only=crust_only, output=o_vtk_dir)
    
    assert(os.path.isdir(o_vtk_dir))  # assert the vtk_outputs is generated
    o_env0 = os.path.join(o_vtk_dir, "slab_env0_00144.vtu")  # compare the contents of the envelops
    o_std_env0 = os.path.join(std_vtk_outputs_dir, "test_slab_morphology_std_env0.vtu")
    assert(filecmp.cmp(o_env0, o_std_env0))
    o_env1 = os.path.join(o_vtk_dir, "slab_env1_00144.vtu")  # compare the contents of the envelops
    o_std_env1 = os.path.join(std_vtk_outputs_dir, "test_slab_morphology_std_env1.vtu")
    assert(filecmp.cmp(o_env1, o_std_env1))
    o_trench = os.path.join(o_vtk_dir, "trench_00144.txt")  # compare the contents of the trench.txt file
    o_std_trench = os.path.join(std_vtk_outputs_dir, "test_slab_morphology_std_trench.txt")
    assert(filecmp.cmp(o_trench, o_std_trench))


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


def test_trench_positions():
    '''
    Test utilities from VtkPp.py
    0: read in pvtu file
    test 1: output slab grid
    '''
    case_dir = os.path.join(source_case_dir, 'test_prepare_slab')
    fileout_std = os.path.join(case_dir, "trench_std.txt")
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
    Xmax = Visit_Options.options['XMAX'] * np.pi / 180.0
    VtkP = VTKP(geometry=geometry, Ro=Ro, Xmax=Xmax, slab_shallow_cutoff=70e3,\
    slab_envelop_interval_y=100e3, slab_envelop_interval_z=100e3)
    VtkP.ReadFile(filein)
    field_names = ['T', 'density', 'sp_upper', 'sp_lower']
    VtkP.ConstructPolyData(field_names, include_cell_center=False)
    VtkP.PrepareSlabByPoints(['sp_upper', 'sp_lower'])
    # prepare outputs
    trench_coords_x, trench_coords_y, trench_coords_z = VtkP.ExportSlabInfo()
    outputs = "# trench coordinates: x, y and z\n"
    for i in range(len(trench_coords_x)):
        if i > 0:
            outputs += "\n"
        outputs += str(trench_coords_x[i]) + " " + str(trench_coords_y[i]) + " " + str(trench_coords_z[i])
    fileout = os.path.join(output_path, "trench.txt" )
    with open(fileout, 'w') as fout:
        fout.write(outputs)
    assert(filecmp.cmp(fileout, fileout_std))
    
# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

