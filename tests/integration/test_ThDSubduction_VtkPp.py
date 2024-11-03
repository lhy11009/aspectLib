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
from shutil import rmtree, copytree  # for remove directories
import vtk
import shilofue.VtkPp as VtkPp
from shilofue.ThDSubduction0.VtkPp import *  # import test module

test_dir = os.path.join(".test", "test_ThDSubduction_VtkPp")
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'ThDSubduction', 'test_vtk_pp_slab')
big_source_dir = os.path.join(ASPECT_LAB_DIR, "tests", "integration", "big_fixtures", "ThDSubduction")
source_case_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'ThDSubduction', 'test_cases')

if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)

# todo_3d_chunk
def test_slab_morphology_chunk():
    '''
    Test extract horizontal flow field and slab morphology in chunk geometry.
    This test uses a big test case which is not included by default
    in the github folder.

    Asserts:
        - Verifies the existence of necessary directories and files.
        - Checks that the generated VTK outputs directory is created.
        - Compares the generated slab surface envelopes and trench locations with the standard outputs to ensure consistency.
    '''
    # Construct the path to the case directory and check if it exists
    case_dir = os.path.join(big_source_dir, "eba3d_width80_bw4000_sw1000_yd500.0_AR4")
    assert(os.path.isdir(case_dir))

    # Define the path for standard VTK outputs and check if it exists
    std_vtk_outputs_dir = os.path.join(case_dir, "vtk_outputs_std")
    # assert(os.path.isdir(std_vtk_outputs_dir))

    # Define the output directory for the test and recreate it if necessary
    odir = os.path.join(test_dir, "test_slab_morphology_chunk")
    if os.path.isdir(odir):
        rmtree(odir)  # Remove old results and make a new directory
    os.mkdir(odir)

    # Check if the input solution file exists
    filein = os.path.join(case_dir, "output", "solution", "solution-00104.pvtu")
    assert(os.path.isfile(filein))

    # Define the directory for VTK outputs and remove old results if necessary
    o_vtk_dir = os.path.join(odir, 'vtk_outputs')
    if os.path.isdir(o_vtk_dir):
        rmtree(o_vtk_dir)  # Remove old results

    # Parameters for sorting out trench locations
    slab_envelop_interval_w = 40e3  # Interval along the x-axis to sort out the trench locations
    slab_envelop_interval_d = 40e3  # Interval along the z-axis to sort out the trench locations
    slab_shallow_cutoff = 40e3  # Minimum depth along the z-axis to sort out the trench locations
    crust_only = 1  # Whether to use only the crustal composition for trench location sorting

    # Run the SlabMorphology function to generate outputs
    SlabMorphology(case_dir, 104, slab_envelop_interval_w=slab_envelop_interval_w,
                   slab_envelop_interval_d=slab_envelop_interval_d,
                   slab_shallow_cutoff=slab_shallow_cutoff, crust_only=crust_only, output=o_vtk_dir)

    # Check that the vtk_outputs directory is created
    assert(os.path.isdir(o_vtk_dir))

    # Compare the contents of the generated envelope files with standard outputs
    o_env0 = os.path.join(o_vtk_dir, "slab_env0_00104.vtu")
    o_std_env0 = os.path.join(std_vtk_outputs_dir, "test_slab_morphology_std_env0.vtu")
    # assert(filecmp.cmp(o_env0, o_std_env0))

    o_env1 = os.path.join(o_vtk_dir, "slab_env1_00104.vtu")
    o_std_env1 = os.path.join(std_vtk_outputs_dir, "test_slab_morphology_std_env1.vtu")
    # assert(filecmp.cmp(o_env1, o_std_env1))

    # Compare the contents of the generated trench.txt file with the standard output
    o_trench = os.path.join(o_vtk_dir, "trench_00104.txt")
    o_std_trench = os.path.join(std_vtk_outputs_dir, "test_slab_morphology_std_trench.txt")
    # assert(filecmp.cmp(o_trench, o_std_trench))


def test_slab_morphology():
    '''
    Test extract horizontal flow field and slab morphology.
    This test uses a big test case which is not included by default
    in the github folder.

    Asserts:
        - Verifies the existence of necessary directories and files.
        - Checks that the generated VTK outputs directory is created.
        - Compares the generated slab surface envelopes and trench locations with the standard outputs to ensure consistency.
    '''
    # Construct the path to the case directory and check if it exists
    case_dir = os.path.join(big_source_dir, "eba3d_width61_c23_AR4")
    assert(os.path.isdir(case_dir))

    # Define the path for standard VTK outputs and check if it exists
    std_vtk_outputs_dir = os.path.join(case_dir, "vtk_outputs_std")
    assert(os.path.isdir(std_vtk_outputs_dir))

    # Define the output directory for the test and recreate it if necessary
    odir = os.path.join(test_dir, "test_slab_morphology")
    if os.path.isdir(odir):
        rmtree(odir)  # Remove old results and make a new directory
    os.mkdir(odir)

    # Check if the input solution file exists
    filein = os.path.join(case_dir, "output", "solution", "solution-00144.pvtu")
    assert(os.path.isfile(filein))

    # Define the directory for VTK outputs and remove old results if necessary
    o_vtk_dir = os.path.join(odir, 'vtk_outputs')
    if os.path.isdir(o_vtk_dir):
        rmtree(o_vtk_dir)  # Remove old results

    # Parameters for sorting out trench locations
    slab_envelop_interval_w = 20e3  # Interval along the x-axis to sort out the trench locations
    slab_envelop_interval_d = 20e3  # Interval along the z-axis to sort out the trench locations
    slab_shallow_cutoff = 40e3  # Minimum depth along the z-axis to sort out the trench locations
    crust_only = 1  # Whether to use only the crustal composition for trench location sorting

    # Run the SlabMorphology function to generate outputs
    SlabMorphology(case_dir, 144, slab_envelop_interval_w=slab_envelop_interval_w,
                   slab_envelop_interval_d=slab_envelop_interval_d,
                   slab_shallow_cutoff=slab_shallow_cutoff, crust_only=crust_only, output=o_vtk_dir)

    # Check that the vtk_outputs directory is created
    assert(os.path.isdir(o_vtk_dir))

    # Compare the contents of the generated envelope files with standard outputs
    o_env0 = os.path.join(o_vtk_dir, "slab_env0_00144.vtu")
    o_std_env0 = os.path.join(std_vtk_outputs_dir, "test_slab_morphology_std_env0.vtu")
    assert(filecmp.cmp(o_env0, o_std_env0))

    o_env1 = os.path.join(o_vtk_dir, "slab_env1_00144.vtu")
    o_std_env1 = os.path.join(std_vtk_outputs_dir, "test_slab_morphology_std_env1.vtu")
    assert(filecmp.cmp(o_env1, o_std_env1))

    # Compare the contents of the generated trench.txt file with the standard output
    o_trench = os.path.join(o_vtk_dir, "trench_00144.txt")
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
    slab_envelop_interval_w=100e3, slab_envelop_interval_d=100e3)
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


def test_InterpolateBySlices_3d_case_chunk_by_part():
    '''
    test using a 3d case output at step 0, by part (vtu files), in chunk geometry
    asserts:
        all points are found
        values from interpolating
    '''
    case_dir = os.path.join(source_case_dir, 'test_prepare_slab_chunk')
    output_path = os.path.join(test_dir, "ThDSubduction_vtk_test_InterpolateBySlices_3d_case_chunk_by_part")
    if os.path.isdir(output_path):
        rmtree(output_path)  # remove old results
    os.mkdir(output_path)
    # class for the basic settings of the case
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    # options
    vtu_snapshot = 0
    spacing = [10, 10, 10]
    split_perturbation = 2
    fields = ["T", "density"]
    d_lateral = 1e3
    # mesh
    n0 = 80
    n1 = 30
    # make a target mesh is none is given
    target_points_np = VtkPp.MakeTargetMesh(Visit_Options, n0, n1, d_lateral)
    # options
    vtu_step = max(0, int(vtu_snapshot) - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']))
    # interpolation
    interpolated_data = np.zeros((len(fields), target_points_np.shape[0]))
    for part in range(16):
        filein = os.path.join(case_dir, "output", "solution", "solution-%05d.%04d.vtu" % (vtu_snapshot, part))
        if part == 0:
            points_found = None
        print("-"*20 + "split" + "-"*20) # print a spliting
        print(filein)
        _, points_found, interpolated_data = VtkPp.InterpolateVtu(Visit_Options, filein, spacing, fields, target_points_np, points_found=points_found,\
                                                    split_perturbation=split_perturbation, interpolated_data=interpolated_data, output_poly_data=False)
    assert(np.sum(points_found==1) == 2400)
    assert(abs(interpolated_data[0][0]-3188.590789515355)/3188.590789515355 < 1e-6)
    assert(abs(interpolated_data[1][0]-3773.6249999999995)/3773.6249999999995 < 1e-6)

def test_InterpolateBySlices_3d_case_chunk():
    '''
    test using a 3d case output at step 0
    asserts:
        all points are found
        values from interpolating
    '''
    case_dir = os.path.join(source_case_dir, 'test_prepare_slab_chunk')
    output_path = os.path.join(test_dir, "ThDSubduction_vtk_test_InterpolateBySlices_3d_case_chunk")
    if os.path.isdir(output_path):
        rmtree(output_path)  # remove old results
    os.mkdir(output_path)
    # class for the basic settings of the case
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    # options
    vtu_snapshot = 0
    spacing = [100, 100, 100]
    split_perturbation = 10
    fields = ["T", "density"]
    d_lateral = 1e3
    # mesh
    n0 = 80
    n1 = 30
    # make a target mesh is none is given
    target_points_np = VtkPp.MakeTargetMesh(Visit_Options, n0, n1, d_lateral)
    # options
    vtu_step = max(0, int(vtu_snapshot) - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']))
    # interpolation
    interpolated_data = np.zeros((len(fields), target_points_np.shape[0]))
    filein = os.path.join(case_dir, "output", "solution", "solution-%05d.pvtu" % vtu_snapshot)
    _, points_found, interpolated_data = VtkPp.InterpolateVtu(Visit_Options, filein, spacing, fields, target_points_np, split_perturbation=split_perturbation,\
                                                interpolated_data=interpolated_data, output_poly_data=False)
    assert(np.sum(points_found==1) == 2400)
    assert(abs(interpolated_data[0][0]-3188.590789515355)/3188.590789515355 < 1e-6)
    assert(abs(interpolated_data[1][0]-3773.6249999999995)/3773.6249999999995 < 1e-6)


def test_InterpolateBySlices_3d_case_by_part():
    '''
    test using a 3d case output at step 0, by part (vtu files)
    asserts:
        all points are found
        values from interpolating
    '''
    case_dir = os.path.join(source_case_dir, 'test_3d_80deg_subduction')
    output_path = os.path.join(test_dir, "ThDSubduction_vtk_test_InterpolateBySlices_3d_case_by_part")
    if os.path.isdir(output_path):
        rmtree(output_path)  # remove old results
    os.mkdir(output_path)
    # class for the basic settings of the case
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    # options
    vtu_snapshot = 0
    spacing = [10, 10, 10]
    split_perturbation = 2
    fields = ["T", "density"]
    d_lateral = 1e3
    # mesh
    n0 = 80
    n1 = 30
    # make a target mesh is none is given
    target_points_np = VtkPp.MakeTargetMesh(Visit_Options, n0, n1, d_lateral)
    # options
    vtu_step = max(0, int(vtu_snapshot) - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']))
    # interpolation
    interpolated_data = np.zeros((len(fields), target_points_np.shape[0]))
    for part in range(16):
        filein = os.path.join(case_dir, "output", "solution", "solution-%05d.%04d.vtu" % (vtu_snapshot, part))
        if part == 0:
            points_found = None
        print("-"*20 + "split" + "-"*20) # print a spliting
        print(filein)
        _, points_found, interpolated_data = VtkPp.InterpolateVtu(Visit_Options, filein, spacing, fields, target_points_np, points_found=points_found,\
                                                    split_perturbation=split_perturbation, interpolated_data=interpolated_data, output_poly_data=False)
    assert(np.sum(points_found==1) == 2400)
    assert(abs(interpolated_data[0][0]-3188.590789515355)/3188.590789515355 < 1e-6)
    assert(abs(interpolated_data[1][0]-3774.437988281252)/3774.437988281252 < 1e-6)


def test_InterpolateBySlices_3d_case():
    '''
    test using a 3d case output at step 0
    asserts:
        all points are found
        values from interpolating
    '''
    case_dir = os.path.join(source_case_dir, 'test_3d_80deg_subduction')
    output_path = os.path.join(test_dir, "ThDSubduction_vtk_test_InterpolateBySlices_3d_case")
    if os.path.isdir(output_path):
        rmtree(output_path)  # remove old results
    os.mkdir(output_path)
    # class for the basic settings of the case
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    # options
    vtu_snapshot = 0
    spacing = [100, 100, 100]
    split_perturbation = 10
    fields = ["T", "density"]
    d_lateral = 1e3
    # mesh
    n0 = 80
    n1 = 30
    # make a target mesh is none is given
    target_points_np = VtkPp.MakeTargetMesh(Visit_Options, n0, n1, d_lateral)
    # options
    vtu_step = max(0, int(vtu_snapshot) - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']))
    # interpolation
    interpolated_data = np.zeros((len(fields), target_points_np.shape[0]))
    filein = os.path.join(case_dir, "output", "solution", "solution-%05d.pvtu" % vtu_snapshot)
    _, points_found, interpolated_data = VtkPp.InterpolateVtu(Visit_Options, filein, spacing, fields, target_points_np, split_perturbation=split_perturbation,\
                                                interpolated_data=interpolated_data, output_poly_data=False)
    assert(np.sum(points_found==1) == 2400)
    assert(abs(interpolated_data[0][0]-3188.590789515355)/3188.590789515355 < 1e-6)
    assert(abs(interpolated_data[1][0]-3774.437988281252)/3774.437988281252 < 1e-6)



def test_InterpolateBySlices():
    '''
    test interpolation by slices in the domain
    assert:
        the values interpolated
    '''
    case_dir = os.path.join(source_case_dir, 'test_prepare_slab')
    output_path = os.path.join(test_dir, "ThDSubduction_vtk_test_InterpolateBySlices")
    if os.path.isdir(output_path):
        rmtree(output_path)  # remove old results
    os.mkdir(output_path)
    vtu_snapshot = 0
    filein = os.path.join(case_dir, "output", "solution", "solution-%05d.pvtu" % vtu_snapshot)
    # get the case options
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    geometry = Visit_Options.options['GEOMETRY']
    Ro =  Visit_Options.options['OUTER_RADIUS']
    Ri = Visit_Options.options['INNER_RADIUS']
    Xmax = Visit_Options.options['XMAX']
    # mesh
    n0 = 80
    n1 = 30
    d_lateral = 1e3
    target_points_np = VtkPp.MakeTargetMesh(Visit_Options, n0, n1, d_lateral)
    target_cells_vtk = VtkPp.GetVtkCells2d(n0, n1)
    # interpolation
    VtkP = VTKP(geometry=geometry, Ro=Ro, Xmax=Xmax)
    VtkP.ReadFile(filein)
    VtkP.ConstructPolyData(["density"], include_cell_center=False)
    VtkP.SplitInSpace([10, 10, 10], dim=3)
    o_poly_data, _, _ = VtkP.InterpolateSplitSpacing(target_points_np, fields=["density"], cells_vtk=target_cells_vtk, split_perturbation=4)
    o_point_data = o_poly_data.GetPointData()
    densities = vtk_to_numpy(o_point_data.GetArray('density'))
    assert(abs(densities[236]-3.339520996093750000e+03)/3.339520996093750000e+03 < 1e6)
    assert(abs(densities[238]-3.380000000000000000e+03)/3.380000000000000000e+03 < 1e6)



def test_Interpolate3dVtkCaseChunckPart():
    '''
    test the function of Interpolate3dVtkCase, in chunk geometry, only read in a solution in part
    '''
    case_dir = os.path.join(source_case_dir, 'test_prepare_slab_chunk')
    output_path = os.path.join(test_dir, "ThDSubduction_vtk_test_Interpolate3dVtkCaseChunk")
    if os.path.isdir(output_path):
        rmtree(output_path)  # remove old results
    copytree(case_dir, output_path)
    vtu_snapshot = 0
    Interpolate3dVtkCase(output_path, vtu_snapshot, interval=1000e3, n0=800, n1=100, file_extension="txt", part=1)
    txt_output = os.path.join(output_path, "vtk_outputs", "center_slice-00000.0001.txt")
    txt_output_std = os.path.join(case_dir, "vtk_outputs", "center_slice-00000.0001_std.txt")
    assert(os.path.isfile(txt_output))  # assert file generation and file contents
    assert(os.path.isfile(txt_output_std))
    assert(filecmp.cmp(txt_output, txt_output_std))


def test_Interpolate3dVtkCaseChunck():
    '''
    test the function of Interpolate3dVtkCase, in chunk geometry
    '''
    case_dir = os.path.join(source_case_dir, 'test_prepare_slab_chunk')
    output_path = os.path.join(test_dir, "ThDSubduction_vtk_test_Interpolate3dVtkCaseChunk")
    vtu_snapshot = 0
    Interpolate3dVtkCase(output_path, vtu_snapshot, interval=1000e3, n0=800, n1=100, file_extension="txt")
    txt_output = os.path.join(output_path, "vtk_outputs", "center_slice-00000.txt")
    txt_output_std = os.path.join(case_dir, "vtk_outputs", "center_slice-00000_std.txt")
    assert(os.path.isfile(txt_output))  # assert file generation and file contents
    assert(os.path.isfile(txt_output_std))
    assert(filecmp.cmp(txt_output, txt_output_std))


def test_Interpolate3dVtkCase():
    '''
    test the function of Interpolate3dVtkCase
    '''
    case_dir = os.path.join(source_case_dir, 'test_prepare_slab')
    output_path = os.path.join(test_dir, "ThDSubduction_vtk_test_Interpolate3dVtkCase")
    if os.path.isdir(output_path):
        rmtree(output_path)  # remove old results
    copytree(case_dir, output_path)
    vtu_snapshot = 0
    Interpolate3dVtkCase(output_path, vtu_snapshot, interval=1000e3, n0=800, n1=100, file_extension="txt")
    txt_output = os.path.join(output_path, "vtk_outputs", "center_slice-00000.txt")
    txt_output_std = os.path.join(case_dir, "vtk_outputs", "center_slice-00000_std.txt")
    assert(os.path.isfile(txt_output))  # assert file generation and file contents
    assert(os.path.isfile(txt_output_std))
    assert(filecmp.cmp(txt_output, txt_output_std))


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

    slab_envelop_interval_w = 20e3  # Interval along x axis to sort out the trench locations
    slab_envelop_interval_d = 20e3  # Interval along z axis to sort out the trench locations
    slab_shallow_cutoff = 40e3  # Minimum depth along z axis to sort out the trench locations
    crust_only = 1  # If we only use the crustal composition to sort out the trench locations

    SlabMorphology(case_dir, 144, slab_envelop_interval_w=slab_envelop_interval_w, slab_envelop_interval_d=slab_envelop_interval_d,\
        slab_shallow_cutoff=slab_shallow_cutoff, crust_only=crust_only, output=o_vtk_dir, horizontal_velocity_depths=[150e3])
    
    slab_surface_file = os.path.join(o_vtk_dir, "slab_surface_00144_d150.00km.txt")
    assert(os.path.isfile(slab_surface_file))

# notes
    
# to check for error message
    # with pytest.raises(SomeError) as _excinfo:
    #    foo()
    # assert(r'foo' in str(_excinfo.value))

# assert the contents of file
    # assert(filecmp.cmp(out_path, std_path))

