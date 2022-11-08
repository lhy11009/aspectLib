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
from shilofue.ParsePrm import ParseFromDealiiInput
# from shilofue.Utilities import 
# from matplotlib import pyplot as plt
from shutil import rmtree  # for remove directories
import vtk

test_dir = ".test"
source_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'TwoDSubduction', 'test_vtk_pp_slab')
TwoDSubduction_DIR = os.environ['TwoDSubduction_DIR']
has_project_root = (os.path.isdir(TwoDSubduction_DIR))
ASPECTLIB_PERFORM_TEST_ON_LOCAL_DATA = True


if not os.path.isdir(test_dir):
    # check we have the directory to store test result
    os.mkdir(test_dir)


def test_find_mdd():
    # todo_mdd
    case_dir = os.path.join(ASPECT_LAB_DIR, "tests/integration/fixtures/big_files/test_TwoD_vtk_pp_full")
    assert(os.path.isdir(case_dir))
    vtu_snapshot = 97
    filein = os.path.join(case_dir, "output", "solution", "solution-%05d.pvtu" % vtu_snapshot)
    if not os.path.isfile(filein):
        raise FileExistsError("input file (pvtu) doesn't exist: %s" % filein)
    else:
        print("SlabMorphology: processing %s" % filein)
    Visit_Options = VISIT_OPTIONS(case_dir)
    Visit_Options.Interpret()
    # vtk_option_path, _time, step = PrepareVTKOptions(VISIT_OPTIONS, case_dir, 'TwoDSubduction_SlabAnalysis',\
    # vtu_step=vtu_step, include_step_in_filename=True, generate_horiz=True)
    vtu_step = max(0, int(vtu_snapshot) - int(Visit_Options.options['INITIAL_ADAPTIVE_REFINEMENT']))
    _time, step = Visit_Options.get_time_and_step(vtu_step)
    geometry = Visit_Options.options['GEOMETRY']
    Ro =  Visit_Options.options['OUTER_RADIUS']
    Xmax = Visit_Options.options['XMAX'] * np.pi / 180.0
    VtkP = VTKP(geometry=geometry, Ro=Ro, Xmax=Xmax)
    VtkP.ReadFile(filein)
    field_names = ['T', 'density', 'spcrust', 'spharz', 'velocity']
    VtkP.ConstructPolyData(field_names, include_cell_center=True)
    VtkP.PrepareSlab(['spcrust', 'spharz'])
    mdd = VtkP.FindMDD()
    assert(abs(mdd - 7.9729e+04) / 7.9729e+04 < 1e-3) # check the mdd value
    pass


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


def test_export_slab_info_sph():
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
    assert(abs(trench - 0.6466922210397674)/0.6466922210397674 < 1e-6)
    assert(abs(slab_depth - 191927.42159304488)/191927.42159304488 < 1e-6)
    assert(abs(dip_100 - 1.061791071552635)/1.061791071552635< 1e-6)


def test_export_slab_info_cart():
    '''
    Test slab properties from VtkPp.py
    assert:
        1. value of trench position, slab_depth, dip angle
    Note:
        here the error in dip100 has a lot to do with resolution. If
        we take both the 5th adaptive refinement, both of these angles (chunk and box)
        are around 0.63
    '''
    case_dir = os.path.join(source_dir, 'cartesian')
    prm_file = os.path.join(case_dir, 'case.prm')
    assert(os.path.isfile(prm_file))
    with open(prm_file, 'r') as fin:
        idict = ParseFromDealiiInput(fin)
    geometry = idict['Geometry model']['Model name']
    if geometry == 'chunk':
        Ro = float(idict['Geometry model']['Chunk']['Chunk outer radius'])
    elif geometry == 'box':
        Ro = float(idict['Geometry model']['Box']['Y extent'])
    output_path = os.path.join(test_dir, "TwoDSubduction_vtk_pp_slab")
    if os.path.isdir(output_path):
        rmtree(output_path)  # remove old results
    os.mkdir(output_path)
    filein = os.path.join(case_dir, "output", "solution", "solution-00002.pvtu")
    assert(os.path.isfile(filein))
    VtkP = VTKP(geometry=geometry, Ro=Ro)
    VtkP.ReadFile(filein)
    field_names = ['T', 'density', 'spcrust', 'spharz']
    VtkP.ConstructPolyData(field_names, include_cell_center=True)
    VtkP.PrepareSlab(['spcrust', 'spharz'])
    trench, slab_depth, dip_100 = VtkP.ExportSlabInfo()
    assert(abs(trench - 4210541.75)/4210541.75 < 1e-6)
    assert(abs(slab_depth - 220136.75)/220136.75 < 1e-6)
    assert(abs(dip_100- 0.6702772823940486)/0.6702772823940486 < 1e-6)


def test_export_velocity():
    '''
    test function ExportVelocity
    assert:
        1. velocity of the overiding plate and the subducting plate
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
    field_names = ['T', 'density', 'spcrust', 'spharz', 'velocity']
    VtkP.ConstructPolyData(field_names, include_cell_center=True)
    VtkP.PrepareSlab(['spcrust', 'spharz'])
    vsp, vov = VtkP.ExportVelocity()
    assert(abs(vsp[0]) < 1e-6 and abs(vsp[1]) < 1e-6 and abs(vsp[2]) < 1e-6)  # thest two values are 0.0
    assert(abs(vov[0]) < 1e-6 and abs(vov[1]) < 1e-6 and abs(vov[2]) < 1e-6)
    # assert


def test_slab_analysis():
    ''' 
    test the SlabAnalysis class, only works if the project files are presented
    '''
    source_dir1 = os.path.join(source_dir, 'slab_analysis') 
    # test 1: compute dynamic pressure by post-processing
    if has_project_root and ASPECTLIB_PERFORM_TEST_ON_LOCAL_DATA:
        o_file_std = os.path.join(source_dir1, 'slab_forces_std')
        assert(os.path.isfile(o_file_std))
        o_dir = os.path.join(test_dir, "TwoDSubduction_vtk_pp")
        if not os.path.isdir(o_dir):
            os.mkdir(o_dir)
        o_file = os.path.join(o_dir, "slab_forces")
        if os.path.isfile(o_file):
            os.remove(o_file)
        vtu_snapshot = 105 # 10 Ma
        case_dir = os.path.join(TwoDSubduction_DIR, 'EBA_CDPT3', 'eba_cdpt_SA80.0_OA40.0')
        assert(os.path.isdir(case_dir))
        SlabAnalysis(case_dir, vtu_snapshot, o_file, output_slab=True)
        assert(os.path.isfile(o_file))
        assert(filecmp.cmp(o_file, o_file_std))  # compare file contents
        fig_ofile = os.path.join(o_dir, "slab_forces.png")
        PlotSlabForces(o_file, fig_ofile)
        assert(os.path.isfile(fig_ofile))
    # test 2: with dynamic pressure outputed in vtu file.
    if has_project_root and ASPECTLIB_PERFORM_TEST_ON_LOCAL_DATA:
        o_file_std = os.path.join(source_dir1, 'slab_forces_dp_std')
        assert(os.path.isfile(o_file_std))
        o_dir = os.path.join(test_dir, "TwoDSubduction_vtk_pp")
        if not os.path.isdir(o_dir):
            os.mkdir(o_dir)
        o_file = os.path.join(o_dir, "slab_forces_dp")
        if os.path.isfile(o_file):
            os.remove(o_file)
        vtu_snapshot = 105 # 10 Ma
        case_dir = os.path.join(TwoDSubduction_DIR, 'EBA_CDPT3_dp', 'eba_cdpt_SA80.0_OA40.0')
        assert(os.path.isdir(case_dir))
        SlabAnalysis(case_dir, vtu_snapshot, o_file, output_slab=True)
        assert(os.path.isfile(o_file))
        assert(filecmp.cmp(o_file, o_file_std))  # compare file contents
        fig_ofile = os.path.join(o_dir, "slab_forces_dp.png")
        PlotSlabForces(o_file, fig_ofile)
        assert(os.path.isfile(fig_ofile))
    # test 3: with dynamic pressure outputed in vtu file, use temperature differences
    # as the criteria for slab, dT = 100 K.
    # Note there are files containing more information generated, in "eba_cdpt_SA80.0_OA40.0/vtk_outputs"
    # e.g. the slab_env1_00101.vtp and slab_env0_00101.vtp contains a profile of the slab surface
    if has_project_root and ASPECTLIB_PERFORM_TEST_ON_LOCAL_DATA:
        o_file_std = os.path.join(source_dir1, 'slab_forces_dp_T100_std')
        o_env_std = os.path.join(source_dir1, 'slab_env0_00101_std.vtp')
        assert(os.path.isfile(o_file_std))
        o_dir = os.path.join(test_dir, "TwoDSubduction_vtk_pp")
        if not os.path.isdir(o_dir):
            os.mkdir(o_dir)
        o_file = os.path.join(o_dir, "slab_forces_dp_T100")
        if os.path.isfile(o_file):
            os.remove(o_file)
        vtu_snapshot = 106 # 10.1 Ma
        _time = 10.1e6
        case_dir = os.path.join(TwoDSubduction_DIR, 'EBA_CDPT3_dp', 'eba_cdpt_SA80.0_OA40.0')
        assert(os.path.isdir(case_dir))
        o_env = os.path.join(case_dir, "vtk_outputs", "slab_env0_00101.vtp")
        if os.path.isfile(o_env):
            os.remove(o_env)
        # get the interpolation function for a horizontally averaged T profile
        ha_file = os.path.join(case_dir, "output", "depth_average.txt")
        DepthAverage = DEPTH_AVERAGE_PLOT('DepthAverage')
        DepthAverage.Import(ha_file)
        itp_func = DepthAverage.GetInterpolateFunc(_time, "temperature")
        SlabAnalysis(case_dir, vtu_snapshot, o_file, use_dT=True, output_slab=True, output_poly_data=True, slab_envelop_interval=20e3)  # use the vertical profile from field data
        assert(os.path.isfile(o_env))
        assert(filecmp.cmp(o_env, o_env_std))  # compare file contents
        assert(os.path.isfile(o_file))
        assert(filecmp.cmp(o_file, o_file_std))  # compare file contents
        fig_ofile = os.path.join(o_dir, "slab_forces_dp_T100.png")
        PlotSlabForces(o_file, fig_ofile)
        assert(os.path.isfile(fig_ofile))

